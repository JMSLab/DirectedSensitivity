import io
import os
import re
import sys
import zipfile

import numpy as np
import pandas as pd
from scipy.stats import norm, qmc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

MODEL      = config.CONFIG['model']
SIMULATION = config.CONFIG['simulation']
TOLERANCES = config.CONFIG['tolerances']
ITERATIONS = config.CONFIG['iterations']

COALITIONS    = ('coal', 'oil_gas', 'nuclear', 'renewable')
TARGET        = COALITIONS.index('coal')
SUPPORT_TERMS = ('constant', 'more_regulation', 'less_regulation', 'more_spending',
                 'salience', 'opinion')

PROFILES  = np.array([[(code >> (3 - player)) & 1 for player in range(4)]
                      for code in range(16)], dtype=bool)
DEVIATION = np.array([[code ^ (1 << (3 - player)) for player in range(4)]
                      for code in range(16)])
PLAYERS   = np.arange(4)


class RatioContest:

    def __init__(self, effectiveness_for, effectiveness_against, spending_exponent):
        self.effectiveness_for     = float(effectiveness_for)
        self.effectiveness_against = float(effectiveness_against)
        self.spending_exponent     = float(spending_exponent)

    def Value(self, initial, effort_for, effort_against):
        return ((initial + self.effectiveness_for * effort_for)
                / self.Denominator(effort_for, effort_against))

    def MarginalFor(self, initial, effort_for, effort_against):
        return (self.effectiveness_for
                * (1.0 - initial + self.effectiveness_against * effort_against)
                / self.Denominator(effort_for, effort_against) ** 2)

    def MarginalAgainst(self, initial, effort_for, effort_against):
        return (self.effectiveness_against
                * (initial + self.effectiveness_for * effort_for)
                / self.Denominator(effort_for, effort_against) ** 2)

    def Denominator(self, effort_for, effort_against):
        return (1.0 + self.effectiveness_for * effort_for
                + self.effectiveness_against * effort_against)


def LoadData(path=config.PATHS['data'], directory=config.CONFIG['members']['data']):
    with zipfile.ZipFile(path) as archive:
        def Read(name):
            return np.loadtxt(io.BytesIO(archive.read(f'{directory}/{name}')))
        outcome   = Read('data_outcome.txt')
        content   = Read('data_explanatory.txt')
        coalition = Read('data_explanatory2.txt')
        opinion   = Read('data_polling.txt')
    return {'enactment':     outcome[:, 0].astype(bool),
            'participation': outcome[:, 1:].astype(bool),
            'support':       np.column_stack([np.ones(len(content)), content, opinion]),
            'position':      np.where(coalition[:, :4] == 1, 1.0, -1.0),
            'relevance':     coalition[:, 4:]}


def LoadEstimates(path=config.PATHS['estimates']):
    estimate = pd.read_csv(path).set_index('parameter')['estimate']
    return {'effectiveness_for':     float(estimate['effectiveness_for']),
            'effectiveness_against': float(estimate['effectiveness_against']),
            'spending_exponent':     float(estimate['spending_exponent']),
            'support_coefficients':  estimate[[f'initial_support_{name}'
                                               for name in SUPPORT_TERMS]].values,
            'support_scale':         float(estimate['initial_support_scale']),
            'value_intercept':       estimate[[f'value_intercept_{name}'
                                               for name in COALITIONS]].values,
            'value_relevance':       estimate[[f'value_relevance_{name}'
                                               for name in COALITIONS]].values,
            'value_scale':           estimate[[f'value_scale_{name}'
                                               for name in COALITIONS]].values}


def LoadSpending(path=config.PATHS['data'], member=config.CONFIG['members']['spending']):
    """Total spending of each coalition in dollars; dataread.f90 states it in millions."""
    with zipfile.ZipFile(path) as archive:
        source = archive.read(member).decode()
    millions = dict(re.findall(r'rst_spd\((\d)\)\s*=\s*([\d.]+)d0\*unitm', source))
    return np.array([float(millions[str(index)]) * 1e6
                     for index in range(1, len(COALITIONS) + 1)])


def Baseline():
    return PUBLISHED_MODEL


def Draws(block, count=None, burn_in=None):
    """Kang's shocks: block b starts burn_in + count * b points into the Halton sequence and is
    rounded to the digits of her shock files."""
    count    = SIMULATION['draws'] if count is None else count
    burn_in  = SIMULATION['burn_in'] if burn_in is None else burn_in
    sequence = qmc.Halton(d=5, scramble=False)
    sequence.fast_forward(burn_in + count * block)
    points   = sequence.random(count)
    return {'support': np.round(2.0 * points[:, 0] - 1.0, SIMULATION['digits']),
            'value':   np.round(norm.ppf(points[:, 1:]), SIMULATION['digits'])}


def SimulationBlocks(data):
    """Policy k uses the draws of block k mod blocks."""
    policies = np.arange(len(data['enactment']))
    blocks = []
    for block in range(SIMULATION['blocks']):
        rows = policies[policies % SIMULATION['blocks'] == block]
        subset = {key: value[rows] for key, value in data.items()}
        blocks.append({'rows': rows, 'patterns': CovariatePatterns(subset),
                       'shocks': Draws(block)})
    return blocks


def CovariatePatterns(data):
    covariates = np.column_stack([data['support'], data['position'], data['relevance']])
    unique, index = np.unique(covariates, axis=0, return_inverse=True)
    width = data['support'].shape[1]
    return {'support':   unique[:, :width],
            'position':  unique[:, width:width + 4],
            'relevance': unique[:, width + 4:],
            'index':     index,
            'count':     np.bincount(index, minlength=len(unique)).astype(float)}


DATA           = LoadData()
ESTIMATES      = LoadEstimates()
TOTAL_SPENDING = LoadSpending()
BLOCKS         = SimulationBlocks(DATA)

POLICIES       = len(DATA['enactment'])
LOBBIED        = int(DATA['participation'].any(axis=1).sum())
OBSERVED_CODE  = DATA['participation'] @ (1 << (3 - PLAYERS))

PUBLISHED_MODEL = RatioContest(ESTIMATES['effectiveness_for'],
                               ESTIMATES['effectiveness_against'],
                               ESTIMATES['spending_exponent'])


def Main():
    summary = Simulate(Baseline())
    for name in ('criterion', 'target', 'log_likelihood', 'penalty', 'enactment_rate',
                 'match_rate', 'equilibrium_rate', 'equilibrium_residual'):
        print(f'{name:22s} {summary[name]!r}')
    print(f"{'participation':22s} {summary['participation']!r}")
    print(f"{'spending ($m)':22s} {summary['spending'] / 1e6!r}")


def Criterion(primitive, summary=None):
    return (summary or Simulate(primitive))['criterion']


def Target(primitive, summary=None):
    return (summary or Simulate(primitive))['target']


def Simulate(primitive, estimates=None):
    estimates = ESTIMATES if estimates is None else estimates
    equilibrium = PolicyEquilibrium(primitive, estimates)

    matched = equilibrium['solved'] & (equilibrium['code'] == OBSERVED_CODE[:, None])
    probability = equilibrium['probability']
    initial     = equilibrium['initial']
    outcome     = np.where(DATA['enactment'][:, None], probability, 1.0 - probability)

    likelihood = np.where(matched, outcome, 0.0).mean(axis=1)
    effect = (np.where(matched, np.abs(probability - initial), 0.0).sum(axis=1)
              / np.maximum(matched.sum(axis=1), 1))
    lobbying = matched[:, :, None] & DATA['participation'][:, None, :]
    returns = 100.0 * PolicyAverage(equilibrium['returns'], lobbying)
    value = PolicyAverage(equilibrium['value'],
                           np.repeat(matched[:, :, None], len(COALITIONS), axis=2))
    log_likelihood = float(np.mean(np.log(np.maximum(likelihood, MODEL['likelihood_floor']))))

    spending = equilibrium['spending'].sum(axis=0)
    shortfall = 1.0 - spending / TOTAL_SPENDING
    penalty = MODEL['penalty_weight'] / POLICIES * float(shortfall @ shortfall)

    solved = equilibrium['solved'].mean(axis=1)
    return {'criterion':            -log_likelihood + penalty,
            'target':               float(returns[TARGET]),
            'returns':              returns,
            'value':                value,
            'returns_coverage':     (lobbying.any(axis=1).sum(axis=0)
                                     / np.maximum(DATA['participation'].sum(axis=0), 1)),
            'effect':               100.0 * float(effect.sum()) / LOBBIED,
            'log_likelihood':       log_likelihood,
            'penalty':              penalty,
            'spending':             spending,
            'enactment_rate':       float((equilibrium['enactment']
                                           / np.maximum(solved, 1e-12)).sum() / POLICIES),
            'participation':        equilibrium['participation'].sum(axis=0) / POLICIES,
            'match_rate':           float(matched.mean()),
            'equilibrium_rate':     float(solved.mean()),
            'equilibrium_residual': equilibrium['residual']}


def PolicyEquilibrium(primitive, estimates):
    solved = [(block, SimulateEquilibrium(primitive, block['patterns'], block['shocks'],
                                          estimates))
              for block in BLOCKS]
    order = np.argsort(np.concatenate([block['rows'] for block, _ in solved]))
    expanded = {key: np.concatenate([equilibrium[key][block['patterns']['index']]
                                     for block, equilibrium in solved])[order]
                for key in solved[0][1] if key != 'residual'}
    expanded['residual'] = max(equilibrium['residual'] for _, equilibrium in solved)
    return expanded


def PolicyAverage(values, lobbying):
    """Average over the draws flagged for each coalition, then over the policies with at least
    one such draw."""
    draws = lobbying.sum(axis=1)
    per_policy = np.where(lobbying, values, 0.0).sum(axis=1) / np.maximum(draws, 1)
    covered = draws > 0
    return (per_policy * covered).sum(axis=0) / np.maximum(covered.sum(axis=0), 1)


def SimulateEquilibrium(primitive, patterns, shocks, estimates):
    blocks = max(1, len(patterns['count']) // SIMULATION['pattern_chunk'])
    chunks = [EquilibriumChunk(primitive, patterns, shocks, estimates, rows)
              for rows in np.array_split(np.arange(len(patterns['count'])), blocks)]
    return {'initial':       np.concatenate([chunk['initial'] for chunk in chunks]),
            'probability':   np.concatenate([chunk['probability'] for chunk in chunks]),
            'code':          np.concatenate([chunk['code'] for chunk in chunks]),
            'solved':        np.concatenate([chunk['solved'] for chunk in chunks]),
            'enactment':     np.concatenate([chunk['enactment'] for chunk in chunks]),
            'spending':      np.concatenate([chunk['spending'] for chunk in chunks]),
            'participation': np.concatenate([chunk['participation'] for chunk in chunks]),
            'returns':       np.concatenate([chunk['returns'] for chunk in chunks]),
            'value':         np.concatenate([chunk['value'] for chunk in chunks]),
            'residual':      max(chunk['residual'] for chunk in chunks)}


def EquilibriumChunk(primitive, patterns, shocks, estimates, rows):
    draws = len(shocks['support'])
    initial = norm.cdf((patterns['support'][rows] @ estimates['support_coefficients'])[:, None]
                       + estimates['support_scale'] * shocks['support'][None, :])
    magnitude = (np.exp(estimates['value_intercept']
                        + estimates['value_relevance'] * patterns['relevance'][rows])[:, None, :]
                 * np.exp(estimates['value_scale'] * shocks['value'])[None, :, :])
    position = patterns['position'][rows][:, None, :]

    exponent = primitive.spending_exponent
    power = exponent / (1.0 - exponent)
    weight = (exponent * magnitude) ** power
    capacity_for     = np.einsum('dsj,cj->dsc', weight * (position > 0), PROFILES.astype(float))
    capacity_against = np.einsum('dsj,cj->dsc', weight * (position < 0), PROFILES.astype(float))

    nested = initial[:, :, None]
    effort_for, effort_against = SolveEquilibrium(primitive, nested,
                                                  capacity_for, capacity_against)
    residual = EquilibriumResidual(primitive, nested, capacity_for, capacity_against,
                                   effort_for, effort_against)

    probability = primitive.Value(nested, effort_for, effort_against)
    marginal = np.where(
        position[:, :, None, :] > 0,
        primitive.MarginalFor(nested, effort_for, effort_against)[..., None],
        primitive.MarginalAgainst(nested, effort_for, effort_against)[..., None])
    spend = np.where(PROFILES,
                     (exponent * magnitude[:, :, None, :] * marginal) ** (1.0 / (1.0 - exponent))
                     + MODEL['entry_cost'], 0.0)
    payoff = probability[..., None] * (position * magnitude)[:, :, None, :] - spend

    is_equilibrium = ~(payoff[:, :, DEVIATION, PLAYERS[None, :]] > payoff).any(axis=3)
    code = np.argmax(np.where(is_equilibrium, payoff.sum(axis=3), -np.inf), axis=2)
    solved = is_equilibrium.any(axis=2)

    selected = code[:, :, None]
    chosen_probability = np.take_along_axis(probability, selected, axis=2)[:, :, 0]
    chosen_spend = np.take_along_axis(spend, selected[..., None], axis=2)[:, :, 0, :]
    chosen_payoff = np.take_along_axis(payoff, selected[..., None], axis=2)[:, :, 0, :]
    staying_out = np.take_along_axis(payoff, DEVIATION[code][:, :, None, :], axis=2)[:, :, 0, :]
    participating = PROFILES[code]
    returns = np.where(participating,
                       (chosen_payoff - staying_out) / np.where(participating, chosen_spend, 1.0),
                       0.0)
    return {'initial':       initial,
            'probability':   chosen_probability,
            'code':          code,
            'solved':        solved,
            'enactment':     np.where(solved, chosen_probability, 0.0).sum(axis=1) / draws,
            'spending':      np.where(solved[..., None], chosen_spend, 0.0).sum(axis=1) / draws,
            'participation': np.where(solved[..., None], PROFILES[code], 0.0).sum(axis=1) / draws,
            'returns':       returns,
            'value':         magnitude,
            'residual':      residual}


def SolveEquilibrium(primitive, initial, capacity_for, capacity_against):
    power = primitive.spending_exponent / (1.0 - primitive.spending_exponent)
    effort_for = np.zeros_like(capacity_for)
    effort_against = np.zeros_like(capacity_against)
    for _ in range(ITERATIONS['equilibrium']):
        effort_for = capacity_for * primitive.MarginalFor(
            initial, effort_for, effort_against) ** power
        effort_against = capacity_against * primitive.MarginalAgainst(
            initial, effort_for, effort_against) ** power
    return effort_for, effort_against


def EquilibriumResidual(primitive, initial, capacity_for, capacity_against,
                        effort_for, effort_against):
    power = primitive.spending_exponent / (1.0 - primitive.spending_exponent)
    gap_for = effort_for - capacity_for * primitive.MarginalFor(
        initial, effort_for, effort_against) ** power
    gap_against = effort_against - capacity_against * primitive.MarginalAgainst(
        initial, effort_for, effort_against) ** power
    scale = np.maximum(effort_for + effort_against, 1.0)
    return float(max(np.max(np.abs(gap_for) / scale), np.max(np.abs(gap_against) / scale)))


if __name__ == '__main__':
    Main()
