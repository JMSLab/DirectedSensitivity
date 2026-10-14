"""Directed sensitivity of the present bias to the logistic net-benefit assumption."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import basis  # noqa: E402
import specification as spec  # noqa: E402

from DirectedSensitivity import (DeformationPath, EvaluationTrace,  # noqa: E402
                                 GeneralizedReducedGradient, NestsBaseline,
                                 SensitivityProblem, WriteAnimation)

PRESENT_BIAS_FLOOR = 1e-6
STAGE_PULL = 0.5
SOLVER_ITERATIONS = 150
LEG_POINTS = 21
CDF_GRID = np.linspace(-8.0, 4.0, 481)


def LoadEstimate():
    with open(HERE / 'estimate.json') as stream:
        estimate = json.load(stream)
    return float(estimate['present_bias']), -float(estimate['location'])


def Fit(primitive):
    try:
        return spec.Criterion(primitive)
    except ValueError:
        return np.inf


def CheckWorkingModel(trace, present_bias, kappa):
    """A second local search of the working model from the estimate, counted in the trace."""
    fit = trace.Watch('criterion', lambda x: Fit(spec.Primitive(x[0], spec.Logistic(-x[1]))))
    result = minimize(fit, np.array([present_bias, kappa]), method='L-BFGS-B',
                      bounds=[(1e-3, 1.0), (0.0, 20.0)],
                      options={'ftol': 1e-15, 'gtol': 1e-10, 'maxiter': 200})
    return {'present_bias': float(result.x[0]), 'kappa': float(result.x[1]),
            'criterion': float(result.fun), 'at_estimate': float(fit(np.array([present_bias, kappa])))}


def Safe(function):
    def wrapped(point):
        try:
            value = float(function(point))
        except (ValueError, FloatingPointError, ZeroDivisionError, np.linalg.LinAlgError):
            return -np.inf
        return value if np.isfinite(value) else -np.inf
    return wrapped


def BuildBasisProblem(model, estimate, trace):
    def Criterion(x):
        if not (np.all(np.isfinite(x)) and x[0] > 0.0):
            return np.inf
        try:
            return Fit(model.Primitive(x))
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            return np.inf

    def Target(x):
        return float(x[0])

    constraints = {
        'present_bias_floor': lambda x: float(x[0] - PRESENT_BIAS_FLOOR),
        'kappa_nonnegative': lambda x: float(x[1]),
        'weight_nonnegative': lambda x: float(x[2]),
        'weight_at_most_one': lambda x: float(1.0 - x[2]),
        'pull_nonnegative': lambda x: float(x[3]),
        'pull_at_most_one': lambda x: float(1.0 - x[3]),
        'complement_masses': Safe(lambda x: model.MassFloor(x) if np.all(np.isfinite(x)) and x[0] > 0 else -np.inf),
    }
    for name, constraint in spec.ConstraintSet().items():
        constraints[name] = Safe(lambda x, constraint=constraint:
                                 constraint(model.Primitive(x))
                                 if np.all(np.isfinite(x)) and x[0] > 0 else -np.inf)
    return SensitivityProblem(criterion=Criterion, target=Target, estimate=estimate,
                              constraints=constraints, trace=trace)


def DescribeBasis(problem, model, x, source):
    x = np.asarray(x, dtype=float)
    primitive, values, masses = model.Solve(x)
    distribution = primitive.distribution
    item = {
        'target': float(problem.target(x)),
        'criterion': float(problem.criterion(x)),
        'slack': float(problem.Slack(x)),
        'deviation': float(problem.SignedDeviation(x)),
        'feasible': bool(problem.IsFeasible(x)),
        'primitive': {'present_bias': float(x[0]), 'kappa': float(x[1]),
                      'weight': float(x[2]), 'pull': float(x[3])},
        'constraints': {k: float(v) for k, v in problem.ConstraintReport(x).items()},
        'source': source,
    }
    if x[2] > 0:
        item['primitive']['complement'] = {'atoms': distribution.first.atoms.tolist(),
                                           'masses': distribution.first.masses.tolist()}
        item['primitive']['continuation_values'] = values.tolist()
    return item


def DescribeTwin(problem, primitive, values, source):
    return {
        'target': float(problem.target(primitive)),
        'criterion': float(problem.criterion(primitive)),
        'slack': float(problem.Slack(primitive)),
        'deviation': float(problem.SignedDeviation(primitive)),
        'feasible': bool(problem.IsFeasible(primitive)),
        'primitive': {'present_bias': primitive.present_bias,
                      'atoms': primitive.distribution.atoms.tolist(),
                      'masses': primitive.distribution.masses.tolist(),
                      'continuation_values': values.tolist()},
        'constraints': {k: float(v) for k, v in problem.ConstraintReport(primitive).items()},
        'source': source,
    }


def Leg(start, end, points=LEG_POINTS):
    start, end = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
    return [start + t * (end - start) for t in np.linspace(0.0, 1.0, points)[1:]]


def Main():
    trace = EvaluationTrace(limits=None)
    present_bias, kappa = LoadEstimate()
    working = spec.Primitive(present_bias, spec.Logistic(-kappa))
    budget = spec.Criterion(working)
    recheck = CheckWorkingModel(trace, present_bias, kappa)
    print('working model recheck', recheck, 'budget', budget)

    observed = spec.FILED / spec.AT_RISK
    isotonic = basis.Isotonic(observed, spec.AT_RISK)
    baseline_hazards = spec.FilingProbabilities(working)
    analytic = {
        'isotonic_criterion': basis.NegativeLogLikelihood(isotonic),
        'baseline_hazards_nondecreasing': bool(np.all(np.diff(baseline_hazards) >= 0)),
        'baseline_hazards_range': [float(baseline_hazards.min()), float(baseline_hazards.max())],
        'isotonic_range': [float(isotonic.min()), float(isotonic.max())],
        'lower': 0.0, 'upper': 1.0,
        'lean_hypotheses': {
            'isotonic_nondecreasing': bool(np.all(np.diff(isotonic) >= 0)),
            'last_day_below_one': bool(isotonic[-1] < 1.0),
            'first_day_positive': bool(isotonic[0] > 0.0),
            'mean_bound': float((1.0 - isotonic[-1]) * basis.LOW_ATOM
                                + isotonic[-1] * basis.HIGH_ATOM),
            'isotonic_within_budget': bool(basis.NegativeLogLikelihood(isotonic) <= budget)},
    }
    print('analytic', analytic)

    direct = SensitivityProblem(criterion=Fit, target=spec.Target, estimate=working,
                                constraints=spec.ConstraintSet(), trace=trace)
    twins = {}
    for side, value in (('upper', 1.0), ('lower', PRESENT_BIAS_FLOOR)):
        primitive, values = basis.Twin(value, isotonic)
        twins[side] = DescribeTwin(direct, primitive, values, 'analytic twin')
        twins[side]['hazard_error'] = float(np.abs(spec.FilingProbabilities(primitive)
                                                   - isotonic).max())
        print('twin', side, twins[side]['target'], twins[side]['criterion'],
              twins[side]['slack'], twins[side]['feasible'], twins[side]['hazard_error'])

    model = basis.BasisModel(isotonic)
    estimate = np.array([present_bias, kappa, 0.0, STAGE_PULL])
    problem = BuildBasisProblem(model, estimate, trace)
    nests = NestsBaseline(problem, estimate)
    print('nests', nests)

    anchors = {'upper': np.array([1.0, kappa, 1.0, 1.0]),
               'lower': np.array([PRESENT_BIAS_FLOOR, kappa, 1.0, 1.0])}

    start = np.array([present_bias, kappa, 1.0, STAGE_PULL])
    solved, runs = {}, {}
    for side, sign in (('upper', 1.0), ('lower', -1.0)):
        solver = GeneralizedReducedGradient(problem)
        solver.StartPhase(sign)
        solver._point = start.copy()
        solver.Solve(max_iterations=SOLVER_ITERATIONS)
        found = solver.positive_primitive if sign > 0 else solver.negative_primitive
        solved[side] = np.asarray(found, dtype=float)
        runs[side] = {'start': start.tolist(), 'iterations': solver.iteration,
                      'converged': bool(solver.converged)}
        print('solver', side, runs[side], solved[side], trace.Snapshot())

    candidates = {side: [DescribeBasis(problem, model, anchors[side], 'anchor'),
                         DescribeBasis(problem, model, solved[side], 'solver')]
                  for side in ('upper', 'lower')}
    for side, items in candidates.items():
        for item in items:
            print(side, item['source'], item['primitive']['present_bias'], item['criterion'],
                  item['slack'], item['feasible'])

    def Best(side, sign):
        feasible = [item for item in candidates[side] if item['feasible']]
        return max(feasible, key=lambda item: (sign * item['target'],
                                               item['source'] == 'anchor'))

    upper, lower = Best('upper', 1), Best('lower', -1)

    def Point(item):
        p = item['primitive']
        return np.array([p['present_bias'], p['kappa'], p['weight'], p['pull']])

    improve = Leg(estimate, [present_bias, kappa, 1.0, STAGE_PULL])
    pull = Leg(improve[-1], [present_bias, kappa, 1.0, 1.0])

    def Branch(item, anchor):
        legs = improve + pull + Leg(pull[-1], anchor)
        certificate = Point(item)
        return legs if np.array_equal(certificate, anchor) else legs + [certificate]

    points = ([estimate] + Branch(lower, anchors['lower']))[::-1] + Branch(upper, anchors['upper'])
    path_feasible = all(problem.IsFeasible(x) for x in points)
    print('path feasible', path_feasible, len(points))

    deformation = DeformationPath(problem, points,
                                  curve=lambda x: model.Primitive(x).distribution.Cdf(CDF_GRID))
    WriteAnimation(HERE / 'animation.gif', deformation, CDF_GRID,
                   labels={'x': r'net benefit of filing $x$ (dollars)',
                           'y': r'distribution function $G(x)$',
                           'object': r'present bias $\beta$'},
                   title='Net-benefit distribution along the deformation path',
                   figure=HERE / 'certificates.png')

    snapshot = trace.Snapshot()
    certificates = {
        'baseline': {'target': float(problem.target_at_estimate), 'criterion': budget,
                     'feasible': bool(problem.IsFeasible(estimate)),
                     'primitive': {'present_bias': present_bias, 'kappa': kappa,
                                   'weight': 0.0, 'pull': STAGE_PULL}},
        'lower': lower,
        'upper': upper,
        'analytic': {'lower': analytic['lower'], 'upper': analytic['upper']},
        'spend': snapshot,
        'limits': {'criterion': None, 'target': None},
        'analytic_details': analytic,
        'analytic_certificates': twins,
        'working_model_recheck': recheck,
        'nests_baseline': nests,
        'solver': runs,
        'path_feasible': path_feasible,
        'candidates': candidates,
    }
    with open(HERE / 'certificates.json', 'w') as stream:
        json.dump(certificates, stream, indent=2,
                  default=lambda v: v.tolist() if hasattr(v, 'tolist') else str(v))
        stream.write('\n')
    print(json.dumps(snapshot))


if __name__ == '__main__':
    Main()
