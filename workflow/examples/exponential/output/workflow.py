"""Directed sensitivity for the exponential example: how far the standard deviation moves."""
import json
import math
import sys
from collections import namedtuple
from pathlib import Path

import numpy as np
from scipy.special import gammaincc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[3]))

from DirectedSensitivity import (DeformationPath, EvaluationTrace,
                                 FeasibleSequentialQuadraticProgramming,
                                 GeneralizedReducedGradient, NestsBaseline,
                                 SensitivityProblem, WriteAnimation)

MEAN, DEVIATION, WEIGHT = 1.5, 2.0, 0.5
ESTIMATE = 1.75
LIMITS = {'criterion': 5000, 'target': 5000}
FEASIBILITY_TOLERANCE = 1e-10
TIE_TOLERANCE = 1e-9
POSITIVITY_FLOOR = 1e-9
ITERATIONS = 25
STOPS = 25
GRID = np.linspace(0.0, 10.0, 201)
SOLVERS = (GeneralizedReducedGradient, FeasibleSequentialQuadraticProgramming)

Primitive = namedtuple('Primitive', ('mix', 'alpha', 'mean', 'sd'))

CONSTRAINTS = {'mix_lower': lambda primitive: primitive.mix,
               'mix_upper': lambda primitive: 1.0 - primitive.mix,
               'mean_positive': lambda primitive: primitive.mean - POSITIVITY_FLOOR,
               'sd_positive': lambda primitive: primitive.sd - POSITIVITY_FLOOR}


def Criterion(primitive):
    return ((1.0 - WEIGHT) * (MEAN - primitive.mean) ** 2
            + WEIGHT * (DEVIATION - primitive.sd) ** 2)


def Target(primitive):
    return primitive.sd


def Anchor(sign):
    """The extremal moments of Claim 3, as a point of the nesting model."""
    radius = abs(DEVIATION - MEAN) * math.sqrt(1.0 - WEIGHT)
    return np.array([MEAN, DEVIATION + sign * radius])


def OnDeformation(anchor):
    """Basis-model chart in (mix, alpha): the deformation from A towards ``anchor``."""
    def Map(point):
        mix, alpha = float(point[0]), float(point[1])
        return Primitive(mix, alpha,
                         mix * anchor[0] + (1.0 - mix) * alpha,
                         mix * anchor[1] + (1.0 - mix) * alpha)
    return Map


def OnMoments(point):
    """Basis-model chart in (mean, sd) along alpha = ESTIMATE, for the path."""
    mean, sd = float(point[0]), float(point[1])
    return Primitive((ESTIMATE - mean) / (ESTIMATE - MEAN), ESTIMATE, mean, sd)


def Problem(to_primitive, estimate, trace):
    return SensitivityProblem(Criterion, Target, np.asarray(estimate, dtype=float),
                              constraints=CONSTRAINTS, to_primitive=to_primitive,
                              feasibility_tolerance=FEASIBILITY_TOLERANCE, trace=trace)


def GammaParameters(primitive):
    return primitive.mean ** 2 / primitive.sd ** 2, primitive.sd ** 2 / primitive.mean


def Survival(point):
    """P(X > x) on ``GRID`` for the gamma distribution with the moments of ``point``."""
    shape, scale = GammaParameters(OnMoments(point))
    return gammaincc(shape, GRID / scale)


def Record(problem, point, source):
    primitive = problem.to_primitive(point)
    shape, scale = GammaParameters(primitive)
    return {'target': float(problem.target(point)),
            'criterion': float(problem.criterion(point)),
            'slack': float(problem.Slack(point)),
            'deviation': float(problem.SignedDeviation(point)),
            'feasible': bool(problem.IsFeasible(point)),
            'primitive': {'family': 'gamma', 'mix': primitive.mix,
                          'alpha': primitive.alpha, 'mean': primitive.mean,
                          'sd': primitive.sd, 'shape': shape, 'scale': scale},
            'source': source}


def WorkingModelFit(problem, grid):
    """The fit of A, read off the mix = 0 edge of the basis model."""
    edge = [problem.to_primitive(np.array([0.0, alpha])) for alpha in grid]
    values = [float(problem.criterion(np.array([0.0, alpha]))) for alpha in grid]
    nesting_gap = max(max(abs(primitive.mean - alpha), abs(primitive.sd - alpha))
                      for primitive, alpha in zip(edge, grid))
    return {'minimizer': float(grid[int(np.argmin(values))]),
            'minimum': float(min(values)), 'nesting_gap': float(nesting_gap)}


def Certify(problem, sign):
    """The anchor and each solver's best, where ``sign`` is the direction of interest."""
    candidates = [Record(problem, np.array([1.0, ESTIMATE]), 'anchor')]
    solvers = {}
    for Search in SOLVERS:
        solver = Search(problem).Solve(ITERATIONS)
        reached = solver.positive_primitive if sign > 0 else solver.negative_primitive
        candidates.append(Record(problem, np.asarray(reached, dtype=float),
                                 Search.__name__))
        solvers[Search.__name__] = {'iterations': solver.iteration,
                                    'feasible_candidates': solver.feasible_candidates,
                                    'infeasible_candidates': solver.infeasible_candidates,
                                    'converged': solver.converged}
    feasible = [candidate for candidate in candidates if candidate['feasible']]
    best = max(feasible, key=lambda candidate: sign * candidate['target'])
    anchor = candidates[0]
    if sign * (best['target'] - anchor['target']) <= TIE_TOLERANCE:
        best = anchor
    return best, candidates, solvers


def Path(lower_anchor, upper_anchor):
    """Lower certificate to baseline to upper certificate, along the two deformations."""
    baseline = np.array([ESTIMATE, ESTIMATE])
    def Leg(anchor, mixes):
        return [mix * anchor + (1.0 - mix) * baseline for mix in mixes]
    return (Leg(lower_anchor, np.linspace(1.0, 0.0, STOPS))
            + Leg(upper_anchor, np.linspace(0.0, 1.0, STOPS))[1:])


def Main():
    trace = EvaluationTrace(limits=LIMITS)
    lower_anchor, upper_anchor = Anchor(-1.0), Anchor(1.0)
    problems = {'lower': Problem(OnDeformation(lower_anchor), [0.0, ESTIMATE], trace),
                'upper': Problem(OnDeformation(upper_anchor), [0.0, ESTIMATE], trace)}
    path_problem = Problem(OnMoments, [ESTIMATE, ESTIMATE], trace)

    nesting = {name: NestsBaseline(problem, np.array([0.0, ESTIMATE]))
               for name, problem in problems.items()}
    nesting['path'] = NestsBaseline(path_problem, np.array([ESTIMATE, ESTIMATE]))
    working = WorkingModelFit(problems['upper'], np.linspace(0.75, 2.75, 41))

    lower, lower_candidates, lower_solvers = Certify(problems['lower'], -1.0)
    upper, upper_candidates, upper_solvers = Certify(problems['upper'], 1.0)
    baseline = Record(path_problem, np.array([ESTIMATE, ESTIMATE]), 'estimate')

    deformation = DeformationPath(path_problem, Path(lower_anchor, upper_anchor),
                                  curve=Survival)
    WriteAnimation(HERE / 'animation.gif', deformation, grid=GRID,
                   labels={'x': 'x', 'y': 'P(X > x)', 'object': 'SD',
                           'lower': 'lowest SD', 'upper': 'highest SD'},
                   figure=HERE / 'certificates.png',
                   annotate=lambda axes: axes.set_title(
                       'Certificates: the survival curves at the extremes of the '
                       'standard deviation', fontsize=10, loc='left'))

    certificates = {
        'baseline': {'target': baseline['target'], 'criterion': baseline['criterion'],
                     'primitive': baseline['primitive']},
        'lower': lower, 'upper': upper,
        'analytic': {'lower': float(lower_anchor[1]), 'upper': float(upper_anchor[1])},
        'candidates': {'lower': lower_candidates, 'upper': upper_candidates},
        'nesting': nesting, 'working_model': working,
        'solvers': {'lower': lower_solvers, 'upper': upper_solvers},
        'spend': trace.Snapshot(), 'limits': LIMITS}
    (HERE / 'certificates.json').write_text(json.dumps(certificates, indent=2) + '\n')
    print(json.dumps(certificates, indent=2))


if __name__ == '__main__':
    Main()
