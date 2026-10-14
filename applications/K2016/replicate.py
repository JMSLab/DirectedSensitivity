import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
import specification as spec

RESULTS_PATH = config.output('replicate.txt')

SPENDING_STEPS = (1000.0, 66000.0, 3000000.0)
RATE_TOLERANCE   = 0.50
EFFECT_TOLERANCE = 0.01
RETURNS_TOLERANCE = 5.0
SHAPE_TOLERANCE  = 1e-3


def Main():
    os.makedirs(os.path.dirname(RESULTS_PATH), exist_ok=True)
    with open(RESULTS_PATH, 'w') as stream:
        try:
            sys.stdout = _Tee(sys.__stdout__, stream)
            targets = LoadTargets()
            summary = spec.Simulate(spec.PUBLISHED_MODEL)
            PrintFit(summary, targets)
            PrintEffect(summary, targets)
            PrintReturns(summary, targets)
            PrintMarginalChange(spec.PUBLISHED_MODEL, targets)
            PrintCriterion(summary)
        finally:
            sys.stdout = sys.__stdout__
    CheckAssertions(summary, targets)


class _Tee:
    def __init__(self, *streams):
        self._streams = streams

    def write(self, text):
        for stream in self._streams:
            stream.write(text)

    def flush(self):
        for stream in self._streams:
            stream.flush()


def LoadTargets(path=config.PATHS['targets']):
    return pd.read_csv(path).set_index('moment')['value'].to_dict()


def PrintFit(summary, targets):
    print(f"Model fit ({spec.SIMULATION['draws']} draws)")
    print(f"  Enactment (%)     simulated {100 * summary['enactment_rate']:.2f}"
          f"   published {targets['enactment_rate_predicted']}"
          f"   observed {targets['enactment_rate_observed']}")
    for position, name in enumerate(spec.COALITIONS):
        print(f"  Participation {name:<10} simulated {100 * summary['participation'][position]:.2f}"
              f"   published {targets[f'participation_predicted_{name}']}"
              f"   observed {targets[f'participation_observed_{name}']}")
    for position, name in enumerate(spec.COALITIONS):
        print(f"  Spending {name:<10} ($m) simulated {summary['spending'][position] / 1e6:.2f}"
              f"   published {targets[f'spending_predicted_{name}']}"
              f"   observed {spec.TOTAL_SPENDING[position] / 1e6:.2f}")
    print('  The spending moment is driven by the upper tail of the value distribution;'
          ' it is not asserted below.')


def PrintEffect(summary, targets):
    print('\nEffect of lobbying on the enactment probability (percentage points)')
    print(f"  All lobbied policies  simulated {summary['effect']:.4f}"
          f"   published {targets['lobbying_effect_all']}")
    print(f"  Draws matching the observed participation profile: {summary['match_rate']:.4f}")
    print(f"  Draws with a pure-strategy equilibrium: {summary['equilibrium_rate']:.4f}")


def PrintReturns(summary, targets):
    print('\nAverage value of a policy ($m) and average returns to lobbying (%)')
    for position, name in enumerate(spec.COALITIONS):
        print(f"  {name:<10} value simulated {summary['value'][position] / 1e6:8.2f}"
              f"   published {targets[f'value_{name}']:8.2f}"
              f"      returns simulated {summary['returns'][position]:7.2f}"
              f"   published {targets[f'returns_{name}']:7.2f}")
    print('  The value is driven by the upper tail of the value distribution, like the spending'
          ' moment; it is not asserted below.')


def PrintMarginalChange(primitive, targets):
    print('\nChange in the enactment probability from a sole lobbyist (percentage points)')
    for spending in SPENDING_STEPS:
        support, opposition = MarginalChange(primitive, spending)
        print(f"  ${spending:>12,.0f}  in support {support:.3f}"
              f"   published {targets[f'change_support_{spending:.0f}']}"
              f"      in opposition {opposition:.3f}"
              f"   published {targets[f'change_opposition_{spending:.0f}']}")


def MarginalChange(primitive, spending):
    effort = spending ** primitive.spending_exponent
    support = primitive.Value(0.0, effort, 0.0) - primitive.Value(0.0, 0.0, 0.0)
    opposition = primitive.Value(1.0, 0.0, effort) - primitive.Value(1.0, 0.0, 0.0)
    return 100.0 * support, 100.0 * opposition


def PrintCriterion(summary):
    print('\nCriterion at the published estimates')
    print(f"  criterion {summary['criterion']:.6f}"
          f"   mean log likelihood {summary['log_likelihood']:.6f}"
          f"   spending penalty {summary['penalty']:.6f}")
    print(f"  equilibrium residual {summary['equilibrium_residual']:.2e}")


def CheckAssertions(summary, targets):
    def check(name, simulated, published, tolerance):
        assert abs(simulated - published) < tolerance, (
            f'FAIL {name}: simulated {simulated:.5g}, published {published:.5g} '
            f'(tolerance {tolerance:g})')

    check('enactment rate', 100 * summary['enactment_rate'],
          targets['enactment_rate_predicted'], RATE_TOLERANCE)
    for position, name in enumerate(spec.COALITIONS):
        check(f'participation {name}', 100 * summary['participation'][position],
              targets[f'participation_predicted_{name}'], RATE_TOLERANCE)
    check('lobbying effect', summary['effect'], targets['lobbying_effect_all'], EFFECT_TOLERANCE)
    for position, name in enumerate(spec.COALITIONS):
        check(f'returns {name}', summary['returns'][position], targets[f'returns_{name}'],
              RETURNS_TOLERANCE)
    for spending in SPENDING_STEPS:
        support, opposition = MarginalChange(spec.PUBLISHED_MODEL, spending)
        check(f'change in support at {spending:.0f}', support,
              targets[f'change_support_{spending:.0f}'], SHAPE_TOLERANCE)
        check(f'change in opposition at {spending:.0f}', opposition,
              targets[f'change_opposition_{spending:.0f}'], SHAPE_TOLERANCE)
    assert summary['equilibrium_residual'] < spec.TOLERANCES['equilibrium'], (
        'FAIL equilibrium residual')


if __name__ == '__main__':
    Main()
