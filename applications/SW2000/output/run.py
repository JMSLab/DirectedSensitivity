"""Directed sensitivity of the discount factor to the CRRA felicity, Stock and Wright (2000)."""
import json
import os
import sys
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import numpy as np
import scipy.optimize

HERE = Path(__file__).resolve().parent
APPLICATION = HERE.parent
sys.path.insert(0, str(APPLICATION))

import specification as spec  # noqa: E402
from replicate import HacCovariance  # noqa: E402

from DirectedSensitivity import (DeformationPath, EvaluationTrace,  # noqa: E402
                                 NestsBaseline, SensitivityProblem, WriteAnimation)

SEED = 20000
RAMP = 1e-3
RRA_FLOOR = spec.CURVATURE
RRA_CEILING = spec.RISK_AVERSION_CEILING - 1e-9
BUDGET_MARGIN = 1e-12
LOWER_STARTS = 40
FLIP_ROUNDS = 3
UPPER_STARTS = 8
PROFILE_POINTS = 101
PATH_POINTS = 25
PLOT_POINTS = 400
LIMITS = {'criterion': 100_000, 'target': 100_000}
BARTLETT = 0.5

RETURN, INSTRUMENTS, T = spec.RETURN, spec.INSTRUMENTS, spec.SAMPLE_SIZE
LEVELS = np.union1d(spec.LEVEL, spec.LEVEL_LAG)
KNOTS = np.log(LEVELS)
GAPS = np.diff(KNOTS)
K = len(GAPS)
LOG_REFERENCE = np.log(spec.REFERENCE_LEVEL)
DIM = 1 + 1 + K + 2


def Incidence():
    """Row t holds the signed gaps between C_{t-1} and C_t."""
    now, lag = np.searchsorted(LEVELS, spec.LEVEL), np.searchsorted(LEVELS, spec.LEVEL_LAG)
    incidence = np.zeros((T, K))
    for t in range(T):
        if now[t] > lag[t]:
            incidence[t, lag[t]:now[t]] = 1.0
        elif lag[t] > now[t]:
            incidence[t, now[t]:lag[t]] = -1.0
    return incidence


INCIDENCE = Incidence()


def KnotMatrix():
    """Knot values of rho from the plateaus: ends copy their plateau, interior knots average."""
    matrix = np.zeros((K + 1, K))
    matrix[0, 0], matrix[K, K - 1] = 1.0, 1.0
    for j in range(1, K):
        matrix[j, j - 1], matrix[j, j] = 0.5, 0.5
    return matrix


KNOT_MATRIX = KnotMatrix()
GAP_MATRIX = GAPS[:, None] * ((1 - RAMP) * np.eye(K)
                              + 0.5 * RAMP * (KNOT_MATRIX[:-1] + KNOT_MATRIX[1:]))


class Basis:
    """x = (w, delta_theta, p_1..p_K, delta_alpha, gamma): w * nu(theta) + (1 - w) * alpha."""

    def __init__(self, x):
        x = np.asarray(x, dtype=float)
        self.weight, self.anchor_discount = float(x[0]), float(x[1])
        self.plateau = x[2:2 + K]
        self.working_discount, self.curvature = float(x[2 + K]), float(x[3 + K])
        self.discount = (self.weight * self.anchor_discount
                         + (1 - self.weight) * self.working_discount)
        knots = KNOT_MATRIX @ self.plateau
        nodes = np.column_stack([KNOTS[:-1], KNOTS[:-1] + RAMP * GAPS,
                                 KNOTS[1:] - RAMP * GAPS]).ravel()
        values = np.column_stack([knots[:-1], self.plateau, self.plateau]).ravel()
        self.nodes = np.append(nodes, KNOTS[-1])
        self.rho = np.append(values, knots[-1])
        area = np.diff(self.nodes) * (self.rho[1:] + self.rho[:-1]) / 2
        self.cumulative = np.concatenate([[0.0], np.cumsum(area)])
        self.offset = self._Integral(np.array([LOG_REFERENCE]))[0]

    def _Rho(self, log_level):
        return np.interp(log_level, self.nodes, self.rho)

    def _Integral(self, log_level):
        index = np.clip(np.searchsorted(self.nodes, log_level, side='right') - 1,
                        0, len(self.nodes) - 2)
        start = self.nodes[index]
        return (self.cumulative[index]
                + (log_level - start) * (self.rho[index] + self._Rho(log_level)) / 2)

    def AnchorValue(self, level):
        log_level = np.log(np.asarray(level, dtype=float))
        return np.exp(-(self._Integral(np.atleast_1d(log_level)) - self.offset)).reshape(
            np.shape(level))

    def AnchorRho(self, level):
        return self._Rho(np.log(np.asarray(level, dtype=float)))

    def WorkingValue(self, level):
        return (np.asarray(level, dtype=float) / spec.REFERENCE_LEVEL) ** (-self.curvature)

    def Value(self, level):
        return (self.weight * self.AnchorValue(level)
                + (1 - self.weight) * self.WorkingValue(level))

    def Derivative(self, level):
        level = np.asarray(level, dtype=float)
        return -(self.weight * self.AnchorValue(level) * self.AnchorRho(level)
                 + (1 - self.weight) * self.WorkingValue(level) * self.curvature) / level


def Point(weight, anchor_discount, plateau, working_discount, curvature):
    return np.concatenate([[weight, anchor_discount], plateau,
                           [working_discount, curvature]])


def ReducedParts(z):
    """Criterion and gradient on the w = 1 slice, in z = (delta, p), read off the gaps."""
    discount, plateau = z[0], z[1:]
    ratio = np.exp(-INCIDENCE @ (GAP_MATRIX @ plateau))
    scaled = discount * ratio * RETURN
    error = scaled - 1.0
    rows = error[:, None] * INSTRUMENTS
    moments = rows.sum(axis=0) / np.sqrt(T)
    try:
        direction = np.linalg.solve(HacCovariance(rows), moments)
    except np.linalg.LinAlgError:
        return np.inf, np.zeros_like(z)
    value = float(moments @ direction)
    loading = INSTRUMENTS @ direction
    score = error * loading
    neighbours = score.copy()
    neighbours[1:] += BARTLETT * score[:-1]
    neighbours[:-1] += BARTLETT * score[1:]
    slope = 2 * loading / np.sqrt(T) - (2 / T) * loading * neighbours
    gradient = np.concatenate([[slope @ (ratio * RETURN)],
                               -GAP_MATRIX.T @ (INCIDENCE.T @ (slope * scaled))])
    return value, gradient


def BasisConstraints():
    def weight_in_unit_interval(x):
        return float(min(x[0], 1 - x[0]))

    def plateaus_in_rra_range(x):
        plateau = x[2:2 + K]
        return float(min(np.min(plateau - RRA_FLOOR), np.min(RRA_CEILING - plateau)))

    def discounts_in_unit_interval(x):
        return float(min(x[1], 1 - x[1], x[2 + K], 1 - x[2 + K]))

    def curvature_in_range(x):
        return float(min(x[3 + K], RRA_CEILING - x[3 + K]))

    constraints = {'weight_in_unit_interval': weight_in_unit_interval,
                   'plateaus_in_rra_range': plateaus_in_rra_range,
                   'discounts_in_unit_interval': discounts_in_unit_interval,
                   'curvature_in_range': curvature_in_range}
    for name, check in spec.ConstraintSet().items():
        constraints[name] = (lambda x, check=check: check(Basis(x)))
    return constraints


def WorkingProfile(problem):
    """The criterion over the working model, profiled in delta on a grid of gamma."""
    rows = []
    for curvature in np.concatenate([[spec.CURVATURE], np.linspace(0.1, 10, PROFILE_POINTS - 1)]):
        def Fit(discount):
            return problem.criterion(Point(0.0, spec.DISCOUNT, np.full(K, spec.CURVATURE),
                                           discount, curvature))
        result = scipy.optimize.minimize_scalar(Fit, bounds=(0.5, 1.0), method='bounded',
                                                options={'xatol': 1e-9})
        rows.append((float(curvature), float(result.x), float(result.fun)))
    return rows


def Watched(trace):
    value = trace.Watch('criterion', lambda z: ReducedParts(z)[0])
    gradient = trace.Watch('criterion', lambda z: ReducedParts(z)[1])
    return value, gradient


def Extremize(sign, start, budget, value, gradient):
    bounds = [(1e-3, 1.0)] + [(RRA_FLOOR, RRA_CEILING)] * K
    fit = {'type': 'ineq', 'fun': lambda z: budget - BUDGET_MARGIN - value(z), 'jac': lambda z: -gradient(z)}
    result = scipy.optimize.minimize(lambda z: -sign * z[0], start,
                                     jac=lambda z: np.r_[-sign, np.zeros(K)], bounds=bounds,
                                     constraints=[fit], method='SLSQP',
                                     options={'maxiter': 400, 'ftol': 1e-13})
    return np.clip(result.x, [b[0] for b in bounds], [b[1] for b in bounds])


def BestFitAtUnitDiscount(start, value, gradient):
    result = scipy.optimize.minimize(lambda p: value(np.r_[1.0, p]), start,
                                     jac=lambda p: gradient(np.r_[1.0, p])[1:],
                                     bounds=[(RRA_FLOOR, RRA_CEILING)] * K, method='L-BFGS-B',
                                     options={'maxiter': 2000, 'ftol': 1e-15, 'gtol': 1e-10})
    return np.r_[1.0, result.x]


def Starts(rng, count, low, high):
    starts = []
    for index in range(count):
        discount = rng.uniform(low, high)
        if index % 2:
            plateau = rng.choice([RRA_FLOOR, RRA_CEILING], size=K)
        else:
            plateau = rng.uniform(RRA_FLOOR, RRA_CEILING, size=K)
        starts.append(np.r_[discount, plateau])
    return starts


def Describe(problem, x, source):
    basis = Basis(x)
    rho = basis.AnchorRho(LEVELS)
    return {'target': float(problem.target(x)),
            'criterion': float(problem.criterion(x)),
            'slack': float(problem.Slack(x)),
            'deviation': float(problem.SignedDeviation(x)),
            'feasible': bool(problem.IsFeasible(x)),
            'primitive': {'weight': basis.weight, 'anchor_discount': basis.anchor_discount,
                          'plateau_rra': basis.plateau.tolist(),
                          'working_discount': basis.working_discount,
                          'working_curvature': basis.curvature,
                          'discount': basis.discount,
                          'levels': LEVELS.tolist(),
                          'rra_at_levels': rho.tolist(),
                          'marginal_utility_at_levels': basis.Value(LEVELS).tolist()},
            'constraints': {name: float(v) for name, v in problem.ConstraintReport(x).items()},
            'source': source}


def Leg(problem, certificate):
    """Mix the anchor into the working model, re-fitting the working discount at each weight."""
    plateau, anchor_discount = certificate[2:2 + K], certificate[1]
    points = []
    for weight in np.linspace(0, 1, PATH_POINTS)[1:-1]:
        def Fit(discount):
            return problem.criterion(Point(weight, anchor_discount, plateau, discount,
                                           spec.CURVATURE))
        result = scipy.optimize.minimize_scalar(Fit, bounds=(0.5, 1.0), method='bounded',
                                                options={'xatol': 1e-10})
        points.append(Point(weight, anchor_discount, plateau, result.x, spec.CURVATURE))
    return points + [certificate]


def Main():
    rng = np.random.default_rng(SEED)
    trace = EvaluationTrace(limits=LIMITS)
    estimate = Point(0.0, spec.DISCOUNT, np.full(K, spec.CURVATURE), spec.DISCOUNT,
                     spec.CURVATURE)
    problem = SensitivityProblem(criterion=lambda x: spec.Criterion(Basis(x)),
                                 target=lambda x: spec.Target(Basis(x)), estimate=estimate,
                                 constraints=BasisConstraints(), trace=trace)
    budget = problem.criterion_at_estimate
    nests = NestsBaseline(problem, estimate)
    print('budget', budget, 'nests', nests)

    profile = WorkingProfile(problem)
    best_working = min(profile, key=lambda row: row[2])
    print('working profile minimum (gamma, delta, q)', best_working)

    value, gradient = Watched(trace)
    z_check = np.r_[0.95, rng.uniform(RRA_FLOOR, RRA_CEILING, K)]
    agreement = abs(value(z_check) - problem.criterion(Point(1.0, z_check[0], z_check[1:],
                                                             spec.DISCOUNT, spec.CURVATURE)))
    print('reduced criterion agreement', agreement)

    lower_runs = []
    floor_start = np.r_[spec.DISCOUNT, np.full(K, RRA_FLOOR)]
    for start in [floor_start] + Starts(rng, LOWER_STARTS - 1, 0.88, 0.96):
        z = Extremize(-1, start, budget, value, gradient)
        lower_runs.append((float(z[0]), float(value(z)), z))
    z_lower = min((run for run in lower_runs if run[1] <= budget), key=lambda run: run[0])[2]
    flips = 0
    for _ in range(FLIP_ROUNDS):
        improved = False
        for gap in range(K):
            start = z_lower.copy()
            start[1 + gap] = RRA_CEILING if start[1 + gap] < RRA_CEILING / 2 else RRA_FLOOR
            z = Extremize(-1, start, budget, value, gradient)
            flips += 1
            lower_runs.append((float(z[0]), float(value(z)), z))
            if value(z) <= budget and z[0] < z_lower[0] - 1e-10:
                z_lower, improved = z, True
        if not improved:
            break
    feasible_lower = [run for run in lower_runs if run[1] <= budget]
    lower_values = sorted(run[0] for run in feasible_lower)
    print('lower runs', len(feasible_lower), 'of', len(lower_runs), 'flips', flips,
          lower_values[:10], trace.Snapshot())

    upper_runs = []
    for start in [np.r_[spec.DISCOUNT, np.full(K, spec.CURVATURE)]] + Starts(
            rng, UPPER_STARTS - 1, 0.9, 1.0):
        z = Extremize(1, start, budget, value, gradient)
        upper_runs.append((float(z[0]), float(value(z)), z))
    feasible_upper = [run for run in upper_runs if run[1] <= budget]
    z_solver_upper = max(feasible_upper, key=lambda run: (run[0], -run[1]))[2]
    z_upper = BestFitAtUnitDiscount(z_solver_upper[1:], value, gradient)
    print('upper solver', z_solver_upper[0], value(z_solver_upper),
          'unit-discount best fit', value(z_upper))

    def Full(z):
        return Point(1.0, z[0], z[1:], spec.DISCOUNT, spec.CURVATURE)

    lower = Describe(problem, Full(z_lower), 'solver')
    upper_candidates = [Describe(problem, Full(z_solver_upper), 'solver'),
                        Describe(problem, Full(z_upper), 'analytic')]
    upper = max((item for item in upper_candidates if item['feasible']),
                key=lambda item: (item['target'], item['slack']))
    for item in [lower] + upper_candidates:
        print(item['source'], item['target'], item['criterion'], item['slack'], item['feasible'])

    lower_x, upper_x = Full(z_lower), Full(z_upper if upper is upper_candidates[1]
                                          else z_solver_upper)
    points = Leg(problem, lower_x)[::-1] + [estimate] + Leg(problem, upper_x)
    path_feasible = [bool(problem.IsFeasible(x)) for x in points]
    print('path feasible', all(path_feasible), [i for i, ok in enumerate(path_feasible) if not ok],
          [float(problem.Slack(x)) for x in points])

    plot_grid = np.linspace(*spec.LEVEL_INTERVAL, PLOT_POINTS)
    deformation = DeformationPath(problem, points,
                                  curve=lambda x: np.log(Basis(x).Value(plot_grid)))
    WriteAnimation(HERE / 'animation.gif', deformation, plot_grid,
                   labels={'x': r'consumption level $C$', 'y': r"$\log u'(C)$",
                           'object': r'$\delta$'},
                   title='Log marginal utility along the deformation path',
                   figure=HERE / 'certificates.png')

    snapshot = trace.Snapshot()
    certificates = {
        'baseline': {'target': float(problem.target(estimate)), 'criterion': budget,
                     'feasible': bool(problem.IsFeasible(estimate)),
                     'primitive': {'discount': spec.DISCOUNT, 'curvature': spec.CURVATURE}},
        'lower': lower,
        'upper': upper,
        'analytic': {'lower': None, 'upper': 1.0},
        'spend': snapshot,
        'limits': LIMITS,
        'nests_baseline': nests,
        'working_profile_minimum': {'curvature': best_working[0], 'discount': best_working[1],
                                    'criterion': best_working[2]},
        'reduced_criterion_agreement': agreement,
        'lower_local_optima': lower_values,
        'upper_candidates': upper_candidates,
        'path_feasible': all(path_feasible),
        'lower_plateau_counts': {'floor': int(np.sum(z_lower[1:] < 1e-6)),
                                 'ceiling': int(np.sum(z_lower[1:] > RRA_CEILING - 1e-6)),
                                 'interior': int(np.sum((z_lower[1:] >= 1e-6)
                                                        & (z_lower[1:] <= RRA_CEILING - 1e-6)))},
        'settings': {'seed': SEED, 'ramp': RAMP, 'rra_floor': RRA_FLOOR,
                     'lower_starts': LOWER_STARTS, 'upper_starts': UPPER_STARTS},
    }
    with open(HERE / 'certificates.json', 'w') as stream:
        json.dump(certificates, stream, indent=2,
                  default=lambda v: v.tolist() if hasattr(v, 'tolist') else str(v))
        stream.write('\n')
    print(json.dumps(snapshot))


if __name__ == '__main__':
    Main()
