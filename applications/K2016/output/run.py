"""K2016 directed sensitivity: assemble certificates, verify them, and draw the deformation path.

Reads the simulations logged in evaluations.jsonl by explore.py and search.py; a point missing
from the log would be simulated and logged.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from basis import NAMES, BoxCoxContest, Constraints, Evaluator, Key, LOG, Pad  # noqa: E402
from search import BUDGET, ESTIMATE  # noqa: E402

from DirectedSensitivity import (DeformationPath, EvaluationTrace, NestsBaseline,  # noqa: E402
                                 SensitivityProblem, WriteAnimation)

LIMITS = {'criterion': 2000, 'target': 2000}
PUBLISHED = Pad([0.2885, 0.0, 0.0, 1.0, 1.0, np.log(50.0), np.log(50.0)])
GRID = np.linspace(0.0, 40.0, 401)


def Describe(problem, evaluator, point, source=None):
    record = evaluator(point)
    item = {'target': float(problem.target(point)), 'criterion': float(problem.criterion(point)),
            'slack': float(problem.Slack(point)), 'deviation': float(problem.SignedDeviation(point)),
            'feasible': bool(problem.IsFeasible(point)),
            'primitive': dict(zip(NAMES, (float(v) for v in Pad(point)))),
            'log_likelihood': record['log_likelihood'], 'penalty': record['penalty'],
            'spending_millions': record['spending'],
            'constraints': {k: float(v) for k, v in problem.ConstraintReport(point).items()}}
    if source is not None:
        item['source'] = source
    return item


def Duplicate(point):
    """The estimate logged under another inert Box-Cox scale."""
    return np.array_equal(point[:5], ESTIMATE[:5]) and not np.any(point[7:9]) \
        and not np.allclose(point, ESTIMATE, rtol=0.0, atol=1e-9)


def OpponentImpact(point):
    contest = BoxCoxContest(point)
    return contest.Impact(GRID, 'against') / BoxCoxContest(ESTIMATE).effectiveness_against


def Main():
    evaluator = Evaluator()
    simulations = sum(1 for _ in open(LOG))
    trace = EvaluationTrace(limits=LIMITS).Restore({'criterion': simulations,
                                                    'target': simulations})
    problem = SensitivityProblem(criterion=lambda x: evaluator(x)['criterion'],
                                 target=lambda x: evaluator(x)['target'],
                                 estimate=ESTIMATE, constraints=Constraints())
    assert problem.criterion_at_estimate == BUDGET
    nests = NestsBaseline(problem, ESTIMATE)

    points = [Pad(r['point']) for r in evaluator.records.values()
              if problem.IsFeasible(Pad(r['point'])) and not Duplicate(Pad(r['point']))]
    points = sorted(points, key=lambda p: problem.target(p))
    lower, upper = points[0], points[-1]

    deformation = DeformationPath(problem, points, curve=OpponentImpact)
    WriteAnimation(HERE / 'animation.gif', deformation, GRID,
                   labels={'x': "opponents' effort $T_a$",
                           'y': r"opponents' impact $f_a(T_a) / \hat\beta_a$",
                           'object': 'coal return (%)'},
                   figure=HERE / 'certificates.png')

    certificates = {
        'baseline': {'target': problem.target_at_estimate,
                     'criterion': problem.criterion_at_estimate,
                     'primitive': dict(zip(NAMES, (float(v) for v in ESTIMATE)))},
        'lower': Describe(problem, evaluator, lower, 'solver'),
        'upper': Describe(problem, evaluator, upper, 'solver'),
        'analytic': {'lower': 0.0, 'upper': None},
        'spend': trace.Snapshot(),
        'limits': LIMITS,
        'published_estimate': {'target': evaluator(PUBLISHED)['target'],
                               'criterion': evaluator(PUBLISHED)['criterion'],
                               'primitive': dict(zip(NAMES, (float(v) for v in PUBLISHED)))},
        'nests_baseline': nests,
        'path_points': len(points),
        'feasible_points_found': len(points),
    }
    with open(HERE / 'certificates.json', 'w') as stream:
        json.dump(certificates, stream, indent=2)
        stream.write('\n')
    for side in ('lower', 'upper'):
        item = certificates[side]
        print(side, item['target'], item['criterion'], item['slack'], item['feasible'],
              item['primitive'])
    print('baseline', certificates['baseline']['target'], certificates['baseline']['criterion'])
    print('nests', nests, 'spend', trace.Snapshot(), 'path', len(points))


if __name__ == '__main__':
    Main()
