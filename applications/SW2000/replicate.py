import html
import os
import re
import sys

import numpy as np
import scipy.optimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

DATA_PATH    = config.PATHS['data']
RESULTS_PATH = config.output('replicate.txt')

DEFLATOR  = 'ppi_annual_average_1967'
DEFLATORS = ('ppi_annual_average_1967', 'ppi_january_1967',
             'ppi_january_1982', 'consumption_deflator_1972')
LAGS      = 1

SERIES = {'year': 10, 'stock_price_january': 11, 'dividends': 12,
          'ppi_january_1982': 15, 'ppi_january_1967': 16, 'ppi_annual_average_1967': 17,
          'consumption_deflator_1972': 18, 'real_consumption_per_capita': 19}

PUBLISHED = {'delta': 0.903, 'delta_se': 0.022,
             'gamma': 0.132, 'gamma_se': 1.037, 'J': 0.69}

START  = np.array([0.95, 1.0])
GRID   = (np.linspace(0.6, 1.1, 51), np.linspace(-6.0, 60.0, 133))
SEARCH = {'xatol': 1e-12, 'fatol': 1e-14, 'maxiter': 20000, 'maxfev': 20000}


def Main():
    os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
    outfile = open(RESULTS_PATH, 'w')
    sys.stdout = _Tee(sys.__stdout__, outfile)

    panel = LoadPanel(DATA_PATH, DEFLATOR)
    estimate = ComputeEstimates(panel)
    deflators = ComputeDeflatorSensitivity()
    PrintSample(panel)
    PrintEstimates(estimate)
    PrintDeflators(deflators)
    tolerance = 2.0
    CheckAssertions(estimate, tolerance)


class _Tee:
    def __init__(self, *files):
        self._files = files

    def write(self, text):
        for handle in self._files:
            handle.write(text)

    def flush(self):
        for handle in self._files:
            handle.flush()


def LoadPanel(path, deflator=DEFLATOR):
    table = _ReadTable(path)
    year, price, dividend = table['year'], table['stock_price_january'], table['dividends']
    level, prices = table['real_consumption_per_capita'], table[deflator]

    nominal = np.full(len(year), np.nan)
    nominal[:-1] = (price[1:] + dividend[:-1]) / price[:-1]
    real = nominal * prices / np.roll(prices, -1)

    growth = level / np.roll(level, 1)
    rows = [(int(year[t]), growth[t], real[t - 1], growth[t - 1], real[t - 2],
             level[t], level[t - 1])
            for t in range(2, len(year))]
    rows = [row for row in rows if np.all(np.isfinite(row[1:]))]
    columns = np.array(rows, dtype=float)
    return {'year': columns[:, 0].astype(int), 'growth': columns[:, 1],
            'ret': columns[:, 2], 'growth_lag': columns[:, 3],
            'ret_lag': columns[:, 4], 'level': columns[:, 5],
            'level_lag': columns[:, 6], 'T': len(rows),
            'instruments': np.column_stack([np.ones(len(rows)), columns[:, 3],
                                            columns[:, 4]])}


def _ReadTable(path):
    """Each series is one table cell holding a <br>-separated list, aligned on 1871."""
    with open(path, encoding='utf-8', errors='replace') as handle:
        cells = re.findall(r'<td\b[^>]*>(.*?)</td>', handle.read(), flags=re.S | re.I)
    series = {name: _ReadCell(cells[index]) for name, index in SERIES.items()}
    length = len(series['year'])
    return {name: np.array([float(item) if item else np.nan
                            for item in (values + [''] * length)[:length]])
            for name, values in series.items()}


def _ReadCell(cell):
    text = re.sub(r'<br\s*/?>', '|', cell, flags=re.I)
    text = html.unescape(re.sub(r'<[^>]+>', '', text)).replace('\xa0', ' ')
    return [item.strip().replace(',', '') for item in text.split('|')]


def MomentRows(panel, theta):
    discount, curvature = float(theta[0]), float(theta[1])
    pricing_error = discount * panel['growth'] ** (-curvature) * panel['ret'] - 1.0
    return pricing_error[:, None] * panel['instruments']


def Moments(panel, theta):
    return MomentRows(panel, theta).sum(axis=0) / np.sqrt(panel['T'])


def HacCovariance(rows, lags=LAGS):
    """Newey-West, uncentered, Bartlett weights."""
    sample_size = rows.shape[0]
    covariance = rows.T @ rows / sample_size
    for lag in range(1, lags + 1):
        autocovariance = rows[lag:].T @ rows[:-lag] / sample_size
        weight = 1.0 - lag / (lags + 1.0)
        covariance = covariance + weight * (autocovariance + autocovariance.T)
    return covariance


def TwoStepWeight(panel, theta):
    return np.linalg.inv(HacCovariance(MomentRows(panel, theta)))


def FixedWeightObjective(panel, theta, weight):
    value = Moments(panel, theta)
    return float(value @ weight @ value)


def CueObjective(panel, theta):
    value = Moments(panel, theta)
    try:
        return float(value @ np.linalg.solve(HacCovariance(MomentRows(panel, theta)),
                                             value))
    except np.linalg.LinAlgError:
        return np.inf


def FitFirstStep(panel):
    instruments = panel['instruments']
    weight = np.linalg.inv(instruments.T @ instruments / panel['T'])
    return FitFixedWeight(panel, weight, START)


def FitFixedWeight(panel, weight, start):
    result = scipy.optimize.minimize(
        lambda theta: FixedWeightObjective(panel, theta, weight),
        start, method='Nelder-Mead', options=SEARCH)
    if not result.success:
        raise RuntimeError(f'fixed-weight GMM did not converge: {result.message}')
    return result.x


def FitCue(panel, start, bounds=None):
    result = scipy.optimize.minimize(lambda theta: CueObjective(panel, theta),
                                     start, method='Nelder-Mead', bounds=bounds, options=SEARCH)
    if not result.success:
        raise RuntimeError(f'CUE did not converge: {result.message}')
    return result.x


def Jacobian(panel, theta, step=1e-6):
    columns = []
    for index in range(len(theta)):
        forward, backward = np.array(theta, float), np.array(theta, float)
        spacing = step * max(1.0, abs(theta[index]))
        forward[index] += spacing
        backward[index] -= spacing
        columns.append((Moments(panel, forward) - Moments(panel, backward))
                       / (2.0 * spacing))
    return np.column_stack(columns)


def StandardErrors(panel, theta, weight):
    slope = Jacobian(panel, theta)
    bread = np.linalg.inv(slope.T @ weight @ slope)
    meat = slope.T @ weight @ HacCovariance(MomentRows(panel, theta)) @ weight @ slope
    return np.sqrt(np.diag(bread @ meat @ bread))


def CueGridMinimum(panel):
    best = (np.nan, np.nan, np.inf)
    for discount in GRID[0]:
        for curvature in GRID[1]:
            value = CueObjective(panel, (discount, curvature))
            if value < best[2]:
                best = (discount, curvature, value)
    return best


def ComputeEstimates(panel):
    first = FitFirstStep(panel)
    weight = TwoStepWeight(panel, first)
    two_step = FitFixedWeight(panel, weight, first)
    cue = FitCue(panel, two_step)
    coarse = CueGridMinimum(panel)
    refined = FitCue(panel, np.array(coarse[:2]))
    return {'first': first,
            'two_step': two_step,
            'two_step_se': StandardErrors(panel, two_step, weight),
            'J': FixedWeightObjective(panel, two_step, weight),
            'budget': CueObjective(panel, two_step),
            'cue': cue,
            'cue_objective': CueObjective(panel, cue),
            'coarse': coarse,
            'refined': refined,
            'refined_objective': CueObjective(panel, refined)}


def ComputeDeflatorSensitivity():
    sensitivity = {}
    for name in DEFLATORS:
        panel = LoadPanel(DATA_PATH, name)
        first = FitFirstStep(panel)
        weight = TwoStepWeight(panel, first)
        two_step = FitFixedWeight(panel, weight, first)
        sensitivity[name] = (two_step, FixedWeightObjective(panel, two_step, weight),
                             CueObjective(panel, two_step))
    return sensitivity


def PrintSample(panel):
    print(f'sample             {panel["year"][0]}-{panel["year"][-1]}, T = {panel["T"]}')
    print(f'deflator           {DEFLATOR}')
    print(f'consumption growth [{panel["growth"].min():.6f}, {panel["growth"].max():.6f}]')
    print(f'gross real return  [{panel["ret"].min():.6f}, {panel["ret"].max():.6f}]')


def PrintEstimates(estimate):
    first, two_step, cue = estimate['first'], estimate['two_step'], estimate['cue']
    print(f'\nfirst step         delta={first[0]:.6f}  gamma={first[1]:+.6f}')
    print(f'two-step GMM       delta={two_step[0]:.6f}  gamma={two_step[1]:+.6f}  '
          f'J={estimate["J"]:.6f}  (targets delta={PUBLISHED["delta"]}  '
          f'gamma={PUBLISHED["gamma"]}  J={PUBLISHED["J"]})')
    print(f'  standard errors        {estimate["two_step_se"][0]:.6f}        '
          f'{estimate["two_step_se"][1]:.6f}  '
          f'(published {PUBLISHED["delta_se"]}, {PUBLISHED["gamma_se"]})')
    print(f'  deviation, in published standard errors  '
          f'delta {(two_step[0] - PUBLISHED["delta"]) / PUBLISHED["delta_se"]:+.2f}  '
          f'gamma {(two_step[1] - PUBLISHED["gamma"]) / PUBLISHED["gamma_se"]:+.2f}')
    print(f'\nCUE at the two-step point            S={estimate["budget"]:.6f}')
    print(f'CUE minimum        delta={cue[0]:.6f}  gamma={cue[1]:+.6f}  '
          f'S={estimate["cue_objective"]:.6f}')
    print('\nCUE sweep of delta in [0.6, 1.1], gamma in [-6, 60]')
    print('  best grid point  delta={:.6f}  gamma={:+.6f}  S={:.6f}'.format(*estimate['coarse']))
    print(f'  refined          delta={estimate["refined"][0]:.6f}  '
          f'gamma={estimate["refined"][1]:+.6f}  S={estimate["refined_objective"]:.6f}')


def PrintDeflators(sensitivity):
    print('\nsensitivity to the deflator used to make returns real (two-step)')
    for name, (two_step, j_statistic, cue_value) in sensitivity.items():
        print(f'  {name:<26} delta={two_step[0]:.4f}  gamma={two_step[1]:+.4f}  '
              f'J={j_statistic:.4f}  S={cue_value:.4f}')


def CheckAssertions(estimate, tolerance):
    for name, index in (('delta', 0), ('gamma', 1)):
        computed = estimate['two_step'][index]
        deviation = abs(computed - PUBLISHED[name]) / PUBLISHED[f'{name}_se']
        assert deviation < tolerance, (
            f'FAIL {name}: computed {computed:.5g}, target {PUBLISHED[name]:.5g}, '
            f'{deviation:.2f} published standard errors apart')

    j_tolerance = 0.15
    assert np.isclose(estimate['J'], PUBLISHED['J'], rtol=j_tolerance), (
        f'FAIL J: computed {estimate["J"]:.5g}, target {PUBLISHED["J"]:.5g} '
        f'(rtol={j_tolerance})')


if __name__ == '__main__':
    Main()
