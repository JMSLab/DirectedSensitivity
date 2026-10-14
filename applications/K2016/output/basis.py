"""Basis model for K2016: Box-Cox impact functions, and a cached parallel evaluator."""
import json
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
APPLICATION = HERE.parent
sys.path.insert(0, str(APPLICATION))

import specification as spec  # noqa: E402

PUBLISHED = spec.Baseline()

LOG = HERE / 'evaluations.jsonl'
WORKERS = 2
NAMES = ('spending_exponent', 'log_shift_for', 'log_shift_against',
         'curvature_for', 'curvature_against', 'log_scale_for', 'log_scale_against',
         'bump_for', 'bump_against', 'log_bump_scale_for', 'log_bump_scale_against')
DEFAULT_TAIL = (0.0, 0.0, np.log(20.0), np.log(20.0))


class BoxCoxContest:
    """f(T) = beta [kappa B_rho(T / kappa) + phi lambda (1 - exp(-T / lambda))], with
    B_rho(x) = ((1 + x)^rho - 1) / rho, B_0 = log(1 + x), and beta the published slope times
    exp(shift)."""

    def __init__(self, point):
        point = Pad(point)
        (self.spending_exponent, shift_for, shift_against, self.curvature_for,
         self.curvature_against, log_scale_for, log_scale_against, self.bump_for,
         self.bump_against, log_bump_for, log_bump_against) = (float(v) for v in point)
        self.bump_scale_for = np.exp(log_bump_for)
        self.bump_scale_against = np.exp(log_bump_against)
        self.effectiveness_for = PUBLISHED.effectiveness_for * np.exp(shift_for)
        self.effectiveness_against = PUBLISHED.effectiveness_against * np.exp(shift_against)
        self.scale_for = np.exp(log_scale_for)
        self.scale_against = np.exp(log_scale_against)

    def Impact(self, effort, side):
        beta, rho, kappa, phi, scale = self._Side(side)
        main = effort if rho == 1.0 else kappa * BoxCox(effort / kappa, rho)
        if phi == 0.0:
            return beta * main
        return beta * (main - phi * scale * np.expm1(-effort / scale))

    def Slope(self, effort, side):
        beta, rho, kappa, phi, scale = self._Side(side)
        slope = (1.0 + effort / kappa) ** (rho - 1.0)
        if phi == 0.0:
            return beta * slope
        return beta * (slope + phi * np.exp(-effort / scale))

    def Value(self, initial, effort_for, effort_against):
        impact_for = self.Impact(effort_for, 'for')
        return (initial + impact_for) / (1.0 + impact_for + self.Impact(effort_against, 'against'))

    def MarginalFor(self, initial, effort_for, effort_against):
        return (self.Slope(effort_for, 'for')
                * (1.0 - initial + self.Impact(effort_against, 'against'))
                / self.Denominator(effort_for, effort_against) ** 2)

    def MarginalAgainst(self, initial, effort_for, effort_against):
        return (self.Slope(effort_against, 'against')
                * (initial + self.Impact(effort_for, 'for'))
                / self.Denominator(effort_for, effort_against) ** 2)

    def Denominator(self, effort_for, effort_against):
        return 1.0 + self.Impact(effort_for, 'for') + self.Impact(effort_against, 'against')

    def _Side(self, side):
        if side == 'for':
            return (self.effectiveness_for, self.curvature_for, self.scale_for, self.bump_for,
                    self.bump_scale_for)
        return (self.effectiveness_against, self.curvature_against, self.scale_against,
                self.bump_against, self.bump_scale_against)


def BoxCox(x, rho):
    if abs(rho) < 1e-12:
        return np.log1p(x)
    return np.expm1(rho * np.log1p(x)) / rho


def BaselinePoint(log_scale=np.log(50.0)):
    return Pad([PUBLISHED.spending_exponent, 0.0, 0.0, 1.0, 1.0, log_scale, log_scale])


def Pad(point):
    point = np.asarray(point, dtype=float)
    if len(point) == 7:
        point = np.concatenate([point, DEFAULT_TAIL])
    return point


def Constraints():
    """Each nonnegative exactly where the Box-Cox pair lies in the nesting model."""
    return {'exponent_above_zero': lambda x: float(x[0]),
            'exponent_below_one': lambda x: float(1.0 - x[0]),
            'concave_for': lambda x: float(1.0 - x[3]),
            'concave_against': lambda x: float(1.0 - x[4]),
            'bump_for': lambda x: float(Pad(x)[7]),
            'bump_against': lambda x: float(Pad(x)[8])}


def Key(point):
    return tuple(float(f'{v:.12g}') for v in Pad(point))


def _Work(point):
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    summary = spec.Simulate(BoxCoxContest(point))
    return {'point': list(Key(point)), 'criterion': summary['criterion'],
            'target': summary['target'], 'log_likelihood': summary['log_likelihood'],
            'penalty': summary['penalty'], 'match_rate': summary['match_rate'],
            'equilibrium_rate': summary['equilibrium_rate'],
            'spending': (summary['spending'] / 1e6).tolist(),
            'returns': summary['returns'].tolist()}


class Evaluator:
    """Simulations, run in parallel and logged, so each point is simulated once across runs."""

    def __init__(self, path=LOG):
        self.path = path
        self.records = {}
        if path.exists():
            for line in path.read_text().splitlines():
                record = json.loads(line)
                self.records[Key(record['point'])] = record

    def Batch(self, points):
        missing = list(dict.fromkeys(Key(p) for p in points if Key(p) not in self.records))
        if missing:
            with ProcessPoolExecutor(WORKERS, mp_context=multiprocessing.get_context("fork")) as pool:
                for record in pool.map(_Work, missing):
                    self.records[Key(record['point'])] = record
                    with open(self.path, 'a') as stream:
                        stream.write(json.dumps(record) + '\n')
        return [self.records[Key(p)] for p in points]

    def __call__(self, point):
        return self.Batch([point])[0]

    def __len__(self):
        return len(self.records)
