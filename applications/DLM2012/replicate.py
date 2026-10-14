import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import specification as spec

RESULTS_PATH = config.output('replicate.txt')

CONDITIONAL_SD = {'lar': 10.545, 'ecu': 7.994}      # Table IV, column 1
HAPPY_SHARE    = {'lar': 0.518, 'ecu': 0.528}       # Table VI, Panel B
PUBLISHED      = {'criterion': 86.618, 'lar': -1.102, 'ecu': -0.442}


def Main():
    os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
    outfile = open(RESULTS_PATH, 'w')
    sys.stdout = _Tee(sys.__stdout__, outfile)

    shape = {charity: spec.RecoverShape(**spec.ALTRUISM[charity]) for charity in spec.CHARITIES}
    implied = {charity: Implied(*shape[charity]) for charity in spec.CHARITIES}
    happy = {charity: HappyShare(*shape[charity], spec.SCALARS[f's_{charity}'])
             for charity in spec.CHARITIES}
    fit = ComputeFit(spec)
    PrintScalars()
    PrintShape(spec, shape, implied, happy)
    PrintFit(fit)
    tolerance = 0.005  # the published values are rounded to three decimals
    CheckAssertions(implied, happy, fit, tolerance)


class _Tee:
    def __init__(self, *files):
        self._files = files

    def write(self, text):
        for f in self._files:
            f.write(text)

    def flush(self):
        for f in self._files:
            f.flush()


def Implied(mu, sigma):
    alpha = -mu / sigma
    ratio = norm.pdf(alpha) / (1.0 - norm.cdf(alpha))
    return {'share_at_zero': float(norm.cdf(alpha)),
            'conditional_mean': mu + sigma * ratio,
            'conditional_sd': sigma * np.sqrt(1.0 - ratio * (ratio - alpha))}


def HappyShare(mu, sigma, s):
    a_minmin = (1.0 - s) * spec.SCALARS['gi']
    return float((1.0 - norm.cdf((OptOutThreshold(s) - mu) / sigma))
                 / (1.0 - norm.cdf((a_minmin - mu) / sigma)))


def OptOutThreshold(s):
    gi = spec.SCALARS['gi']
    a_minmin, a_min, a_plus = (1 - s) * gi, (1 - s) * (gi + spec.GS), gi + spec.GS
    a_0 = spec.GS / np.log1p(spec.GS / gi)
    if a_min <= a_0 <= a_plus:
        return a_0

    def net(a):
        return a * np.log(a / ((1 - s) * gi)) - a + (1 - s) * gi - s * spec.GS

    return brentq(net, a_minmin - 1e-5, a_min + 1e-5, xtol=1e-14, rtol=1e-15)


def ComputeFit(spec):
    working = spec.PUBLISHED_MODEL
    gap = spec.EMPIRICAL - spec.Moments(working)
    blocks = {name: float(gap[k] @ spec.WEIGHT[np.ix_(k, k)] @ gap[k])
              for name, k in spec.BLOCK.items()}
    return {'criterion': spec.Criterion(working), 'blocks': blocks,
            'lar': spec.TargetLaRabida(working), 'ecu': spec.TargetECU(working)}


def PrintScalars():
    print('Table IV, column 1: the eleven scalars')
    for name, value in spec.SCALARS.items():
        print(f'  {name:<8} {value:g}')


def PrintShape(spec, shape, implied, happy):
    print('\nAltruism distribution, recovered from the Table IV implied quantities')
    for charity in spec.CHARITIES:
        mu, sigma = shape[charity]
        published = spec.ALTRUISM[charity]
        print(f'  {charity}  mu={mu:.5f}  sigma={sigma:.5f}  a_max={spec.A_MAX[charity]:.3f}')
        print(f'       share at zero={implied[charity]["share_at_zero"]:.4f} '
              f'(target {published["share_at_zero"]})  '
              f'E[a|a>0]={implied[charity]["conditional_mean"]:.4f} '
              f'(target {published["conditional_mean"]})')
        print(f'       sd[a|a>0]={implied[charity]["conditional_sd"]:.4f} '
              f'(target {CONDITIONAL_SD[charity]})  '
              f'happy givers={happy[charity]:.4f} (target {HAPPY_SHARE[charity]})')
        print(f'       thresholds a_minmin, a_min, a_plus, a_0 = '
              + ', '.join(f'{x:.4f}' for x in spec.Thresholds(spec.BENCHMARK, charity)))


def PrintFit(fit):
    print('\nFit and welfare at the benchmark estimates')
    print(f'  criterion    {fit["criterion"]:9.4f} (target {PUBLISHED["criterion"]})')
    for name, value in fit['blocks'].items():
        print(f'    {name:<7}   {value:9.4f}')
    for charity in ('lar', 'ecu'):
        print(f'  welfare {charity}  {fit[charity]:+9.4f} (target {PUBLISHED[charity]:+.3f})')


def CheckAssertions(implied, happy, fit, tolerance):
    def check(name, computed, target):
        assert np.isclose(computed, target, rtol=tolerance), (
            f'FAIL {name}: computed {computed:.5g}, target {target:.5g} (rtol={tolerance})'
        )

    for charity in ('lar', 'ecu'):
        for key, target in spec.ALTRUISM[charity].items():
            check(f'{key} ({charity})', implied[charity][key], target)
        check(f'conditional_sd ({charity})', implied[charity]['conditional_sd'],
              CONDITIONAL_SD[charity])
        check(f'happy givers ({charity})', happy[charity], HAPPY_SHARE[charity])
        check(f'welfare ({charity})', fit[charity], PUBLISHED[charity])

    check('criterion', fit['criterion'], PUBLISHED['criterion'])


if __name__ == '__main__':
    Main()
