import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config
from replicate import (DEFLATOR, FitCue, FitFirstStep, FitFixedWeight, HacCovariance,
                       LoadPanel, TwoStepWeight)

DATA_PATH   = config.PATHS['data']
GRID_POINTS = 20_001

RISK_AVERSION_CEILING = 10.0
DISCOUNT_CEILING      = 1.0
CURVATURE_FLOOR       = 1e-12
CURVATURE_BOUND       = 1e-9

PANEL = LoadPanel(DATA_PATH, DEFLATOR)
GROWTH, RETURN = PANEL['growth'], PANEL['ret']
LEVEL, LEVEL_LAG = PANEL['level'], PANEL['level_lag']
INSTRUMENTS, SAMPLE_SIZE = PANEL['instruments'], PANEL['T']

FIRST_STEP = FitFirstStep(PANEL)
TWO_STEP = FitFixedWeight(PANEL, TwoStepWeight(PANEL, FIRST_STEP), FIRST_STEP)
ESTIMATE = FitCue(PANEL, TWO_STEP, bounds=[(None, DISCOUNT_CEILING),
                                          (CURVATURE_BOUND, RISK_AVERSION_CEILING)])
DISCOUNT, CURVATURE = float(ESTIMATE[0]), float(ESTIMATE[1])

LEVEL_INTERVAL = (float(min(LEVEL.min(), LEVEL_LAG.min())),
                  float(max(LEVEL.max(), LEVEL_LAG.max())))
REFERENCE_LEVEL = float(np.sqrt(LEVEL_INTERVAL[0] * LEVEL_INTERVAL[1]))
LEVEL_GRID = np.union1d(np.linspace(*LEVEL_INTERVAL, GRID_POINTS),
                        np.union1d(np.union1d(LEVEL, LEVEL_LAG), [REFERENCE_LEVEL]))


class PowerMarginalUtility:

    def __init__(self, discount, curvature):
        self.discount, self.curvature = float(discount), float(curvature)

    def Value(self, level):
        return (np.asarray(level, dtype=float) / REFERENCE_LEVEL) ** (-self.curvature)

    def Derivative(self, level):
        level = np.asarray(level, dtype=float)
        return -self.curvature / level * (level / REFERENCE_LEVEL) ** (-self.curvature)


WORKING_MODEL = PowerMarginalUtility(DISCOUNT, CURVATURE)


def Main():
    print(f'sample                {PANEL["year"][0]}-{PANEL["year"][-1]}, '
          f'T = {SAMPLE_SIZE}')
    print(f'discount, curvature   {DISCOUNT!r}, {CURVATURE!r}')
    print(f'level interval        {LEVEL_INTERVAL!r}')
    print(f'reference level       {REFERENCE_LEVEL!r}')
    print(f'criterion budget      {Criterion(WORKING_MODEL)!r}')
    print(f'target                {Target(WORKING_MODEL)!r}')
    for name, constraint in ConstraintSet().items():
        print(f'{name:<38}{constraint(WORKING_MODEL)!r}')


def MarginalUtilityRatio(primitive):
    return (np.asarray(primitive.Value(LEVEL), dtype=float)
            / np.asarray(primitive.Value(LEVEL_LAG), dtype=float))


def MomentRows(primitive):
    error = primitive.discount * MarginalUtilityRatio(primitive) * RETURN - 1.0
    return np.asarray(error, dtype=float)[:, None] * INSTRUMENTS


def Moments(primitive):
    return MomentRows(primitive).sum(axis=0) / np.sqrt(SAMPLE_SIZE)


def Covariance(primitive):
    return HacCovariance(MomentRows(primitive))


def Criterion(primitive):
    """Continuously updated: the covariance is re-estimated at every candidate."""
    value = Moments(primitive)
    try:
        return float(value @ np.linalg.solve(Covariance(primitive), value))
    except np.linalg.LinAlgError:
        return np.inf


def Target(primitive):
    return float(primitive.discount)


def RelativeRiskAversion(primitive, level):
    level = np.asarray(level, dtype=float)
    return (-level * np.asarray(primitive.Derivative(level), dtype=float)
            / np.asarray(primitive.Value(level), dtype=float))


def ConstraintSet(grid=LEVEL_GRID):
    def discount_at_most_one(primitive):
        return float(DISCOUNT_CEILING - primitive.discount)

    def positive_marginal_utility(primitive):
        return _Least(np.asarray(primitive.Value(grid), dtype=float))

    def strictly_concave_utility(primitive):
        return _Least(RelativeRiskAversion(primitive, grid) - CURVATURE_FLOOR)

    def relative_risk_aversion_at_most_ten(primitive):
        return _Least(RISK_AVERSION_CEILING - RelativeRiskAversion(primitive, grid))

    def defined_on_level_interval(primitive):
        value = np.asarray(primitive.Value(grid), dtype=float)
        slope = np.asarray(primitive.Derivative(grid), dtype=float)
        return 0.0 if np.all(np.isfinite(value)) and np.all(np.isfinite(slope)) else -np.inf

    def unit_marginal_utility_at_reference(primitive):
        return -float(abs(np.asarray(primitive.Value(REFERENCE_LEVEL), dtype=float) - 1.0))

    return {'discount_at_most_one': discount_at_most_one,
            'positive_marginal_utility': positive_marginal_utility,
            'strictly_concave_utility': strictly_concave_utility,
            'relative_risk_aversion_at_most_ten': relative_risk_aversion_at_most_ten,
            'defined_on_level_interval': defined_on_level_interval,
            'unit_marginal_utility_at_reference': unit_marginal_utility_at_reference}


def _Least(margin):
    """A margin that is not a number is a violation, not a pass."""
    margin = np.asarray(margin, dtype=float)
    return -np.inf if np.any(np.isnan(margin)) else float(np.min(margin))


if __name__ == '__main__':
    Main()
