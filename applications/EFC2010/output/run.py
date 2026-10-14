"""Directed sensitivity of the EFC2010 deadweight loss to the linear demand assumption."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
APPLICATION = HERE.parent
sys.path.insert(0, str(APPLICATION))

import specification as spec  # noqa: E402

from DirectedSensitivity import (DeformationPath, EvaluationTrace,  # noqa: E402
                                 GeneralizedReducedGradient, NestsBaseline,
                                 SensitivityProblem, WriteAnimation)

GAMMA, DELTA = float(spec.AVG_COST_INTERCEPT), float(spec.AVG_COST_SLOPE)
P_C = float(spec.COMPETITIVE_PRICE)
KAPPA = DELTA / (1.0 - DELTA)
PRICES, SHARES, WEIGHTS = spec.PRICE, spec.SHARE, spec.WEIGHTS
LOW, HIGH = spec.PRICE_INTERVAL
LIMITS = {'criterion': 100_000, 'target': 100_000}

SHRINK = 1e-5
TILT = 1e-9
FLAT_TILT = 1e-7
SHARPNESS = 10000.0
CHECK_GRID = np.linspace(LOW, HIGH, 3201)
PLOT_GRID = np.linspace(LOW, HIGH, 801)
CROSSING_TOLERANCE = 1e-9
PATH_POINTS = 25


def AverageCost(price):
    return GAMMA + DELTA * np.asarray(price, dtype=float)


def Phi(start):
    """Deadweight loss per unit of total cost when marginal cost is zero on [start, p_C]."""
    cost = AverageCost(start)
    return float((np.log(P_C / cost) - GAMMA * (1.0 / cost - 1.0 / P_C)) / DELTA)


def DecreasingFit(values, weights):
    """Weighted least-squares nonincreasing fit, by pooling adjacent violators."""
    blocks = []
    for value, weight in zip(values, weights):
        blocks.append([value, weight, 1])
        while len(blocks) > 1 and blocks[-2][0] < blocks[-1][0]:
            later, earlier = blocks.pop(), blocks.pop()
            total = earlier[1] + later[1]
            blocks.append([(earlier[0] * earlier[1] + later[0] * later[1]) / total, total,
                           earlier[2] + later[2]])
    return np.concatenate([[value] * count for value, _, count in blocks])


COSTS = AverageCost(PRICES)
TOTAL_SHARES = SHARES * COSTS
TOTAL_WEIGHTS = WEIGHTS / COSTS ** 2


def TailFit(cap=np.inf):
    """Best total costs H_3..H_6 below ``cap``, and their criterion contribution."""
    fitted = np.minimum(DecreasingFit(TOTAL_SHARES[2:], TOTAL_WEIGHTS[2:]), cap)
    return fitted, float(TOTAL_WEIGHTS[2:] @ (TOTAL_SHARES[2:] - fitted) ** 2)


def UpperReduced(budget):
    fitted = DecreasingFit(TOTAL_SHARES[1:], TOTAL_WEIGHTS[1:])
    rest = float(TOTAL_WEIGHTS[1:] @ (TOTAL_SHARES[1:] - fitted) ** 2)
    first = SHARES[0] + np.sqrt((budget - rest) / WEIGHTS[0])
    values = np.concatenate([[first], fitted / COSTS[1:]])
    return values, float(first * COSTS[0] * Phi(LOW))


def Junction(first, level):
    """Where the marginal-cost-equals-price curve through (p_1, d_1) meets level / AC."""
    flat = first / (P_C - PRICES[0]) ** KAPPA
    return brentq(lambda price: level / (AverageCost(price) * (P_C - price) ** KAPPA) - flat,
                  PRICES[0], P_C - 1e-12, xtol=1e-13)


def LowerValue(first, level):
    if level <= 0.0:
        return 0.0
    return level * Phi(Junction(first, level))


def LowerProfile(level, budget):
    tail, tail_loss = TailFit(level)
    loss = TOTAL_WEIGHTS[1] * (TOTAL_SHARES[1] - level) ** 2 + tail_loss
    if loss > budget:
        return np.inf, None
    first = SHARES[0] + np.sqrt((budget - loss) / WEIGHTS[0])
    if level >= first * COSTS[0]:
        return np.inf, None
    values = np.concatenate([[first, level / COSTS[1]], tail / COSTS[2:]])
    return LowerValue(first, level), values


def LowerReduced(budget):
    levels = np.linspace(0.8 * TOTAL_SHARES[1], TOTAL_SHARES[1], 401)
    values = [LowerProfile(level, budget)[0] for level in levels]
    best = int(np.argmin(values))
    bracket = (levels[max(best - 1, 0)], levels[min(best + 1, len(levels) - 1)])
    result = minimize_scalar(lambda level: LowerProfile(level, budget)[0], bounds=bracket,
                             method='bounded', options={'xatol': 1e-10})
    value, fitted = LowerProfile(result.x, budget)
    return fitted, value, float(result.x)


def FullLowerCheck(budget, starts=20, seed=0):
    """The lower reduced problem over all six values, as a cross-check of the profile."""
    rng = np.random.default_rng(seed)

    def Objective(values):
        try:
            return LowerValue(values[0], values[1] * COSTS[1])
        except ValueError:
            return 1e3

    constraints = [{'type': 'ineq', 'fun': lambda d: budget - WEIGHTS @ (SHARES - d) ** 2}]
    constraints += [{'type': 'ineq',
                     'fun': lambda d, j=j: d[j] * COSTS[j] - d[j + 1] * COSTS[j + 1]}
                    for j in range(5)]
    best = np.inf
    for _ in range(starts):
        start = SHARES + rng.normal(0, 0.01, 6)
        start = DecreasingFit(start * COSTS, np.ones(6)) / COSTS
        result = minimize(Objective, start, method='SLSQP', constraints=constraints,
                          options={'ftol': 1e-12, 'maxiter': 500})
        if all(c['fun'](result.x) >= -1e-9 for c in constraints):
            best = min(best, float(result.fun))
    return best


class Curve:
    """A demand curve from its value and log-derivative; subclasses set ``_Parts``."""

    def Value(self, price):
        return self._Parts(price)[0]

    def Derivative(self, price):
        value, slope = self._Parts(price)
        return value * slope


def Tilted(knots, levels, tilt, origin):
    """Total cost: a monotone cubic through the knots, times a strictly decreasing factor."""
    cubic = PchipInterpolator(knots, levels)
    slope = cubic.derivative()

    def Parts(price):
        factor = np.exp(-tilt * (price - origin))
        level = cubic(price) * factor
        return level, (slope(price) / cubic(price) - tilt)
    return Parts


class UpperAnchor(Curve):
    """Marginal cost (nearly) zero below p_C, total cost through the tail values above."""

    def __init__(self, values):
        top = values[0] * COSTS[0]
        knots = np.concatenate([[LOW, P_C], PRICES[1:], [HIGH]])
        levels = np.concatenate([[top, top], values[1:] * COSTS[1:], [0.5 * values[-1] * COSTS[-1]]])
        self.total = Tilted(knots, levels, TILT, PRICES[0])

    def _Parts(self, price):
        price = np.asarray(price, dtype=float)
        level, log_slope = self.total(price)
        return level / AverageCost(price), log_slope - DELTA / AverageCost(price)


class LowerAnchor(Curve):
    """The smooth maximum of a marginal-cost-at-price curve and a zero-marginal-cost curve."""

    def __init__(self, values):
        self.flat = values[0] / (P_C - PRICES[0]) ** KAPPA
        knots = np.concatenate([[LOW], PRICES[1:], [HIGH]])
        levels = np.concatenate([[values[1] * COSTS[1]], values[1:] * COSTS[1:],
                                 [0.5 * values[-1] * COSTS[-1]]])
        self.total = Tilted(knots, levels, TILT, PRICES[1])

    def _Parts(self, price):
        shape = np.shape(price)
        price = np.atleast_1d(np.asarray(price, dtype=float))
        below = price < P_C
        gap = np.where(below, P_C - price, 1.0)
        log_first = np.where(below, np.log(self.flat) + KAPPA * np.log(gap)
                             - FLAT_TILT * (price - PRICES[0]), -np.inf)
        slope_first = np.where(below, -KAPPA / gap - FLAT_TILT, 0.0)
        level, log_slope = self.total(price)
        log_second = np.log(level / AverageCost(price))
        slope_second = log_slope - DELTA / AverageCost(price)
        stacked = SHARPNESS * np.vstack([log_first, log_second])
        log_value = logsumexp(stacked, axis=0) / SHARPNESS
        shares = np.exp(stacked - SHARPNESS * log_value)
        value = np.exp(log_value)
        slope = shares[0] * slope_first + shares[1] * slope_second
        return value.reshape(shape), slope.reshape(shape)


class Mixture:
    def __init__(self, weight, anchor, working):
        self.weight, self.anchor, self.working = weight, anchor, working

    def Value(self, price):
        return self.weight * self.anchor.Value(price) + (1 - self.weight) * self.working.Value(price)

    def Derivative(self, price):
        return (self.weight * self.anchor.Derivative(price)
                + (1 - self.weight) * self.working.Derivative(price))


def MarginalCost(demand, price):
    return AverageCost(price) + DELTA * demand.Value(price) / demand.Derivative(price)


def CrossingMargin(demand):
    """Nonnegative when p - MC(p) on the grid never turns negative after it is nonnegative."""
    gap = CHECK_GRID - MarginalCost(demand, CHECK_GRID)
    started = np.maximum.accumulate(gap >= -CROSSING_TOLERANCE)
    after = gap[started]
    return float(min(after.min() + CROSSING_TOLERANCE, 1.0)) if after.size else -1.0


def NestingConstraints(demand_of):
    return {
        'decreasing': lambda x: float(-np.max(demand_of(x).Derivative(CHECK_GRID))),
        'ceiling': lambda x: float(1.0 - demand_of(x).Value(LOW)),
        'floor': lambda x: float(demand_of(x).Value(HIGH)),
        'marginal_cost': lambda x: float(np.min(MarginalCost(demand_of(x), CHECK_GRID))),
        'single_crossing': lambda x: CrossingMargin(demand_of(x)),
    }


def BuildBasis(anchors, side=None):
    def Demand(x):
        weight, intercept, slope = (float(v) for v in x)
        if side is None:
            side_, mix = (1 if weight >= 0 else -1), abs(weight)
        else:
            side_, mix = side, weight
        return Mixture(mix, anchors[side_], spec.LinearDemand(intercept, slope))

    constraints = {'weight': lambda x: 1.0 - abs(x[0])}
    if side is not None:
        constraints['branch'] = lambda x: x[0]
    constraints.update(NestingConstraints(Demand))
    return (Demand, lambda x: spec.Target(Demand(x)),
            lambda x: spec.Criterion(Demand(x)), constraints)


def Describe(problem, demand_of, x, source, extra=None):
    weight, intercept, slope = (float(v) for v in x)
    demand = demand_of(x)
    return {'target': float(problem.target(x)),
            'criterion': float(problem.criterion(x)),
            'slack': float(problem.Slack(x)),
            'deviation': float(problem.SignedDeviation(x)),
            'feasible': bool(problem.IsFeasible(x)),
            'primitive': dict({'weight': weight, 'intercept': intercept, 'slope': slope,
                               'values_at_prices': demand.Value(PRICES).tolist(),
                               'efficient_price': float(spec.EfficientPrice(demand))},
                              **(extra or {})),
            'constraints': {k: float(v) for k, v in problem.ConstraintReport(x).items()},
            'source': source}


def Main():
    baseline_demand = spec.Baseline()
    intercept, slope = baseline_demand.intercept, baseline_demand.slope
    budget = spec.Criterion(baseline_demand)

    upper_values, analytic_upper = UpperReduced(budget)
    lower_values, analytic_lower, lower_level = LowerReduced(budget)
    lower_check = FullLowerCheck(budget)
    print('analytic', analytic_lower, analytic_upper, 'full check', lower_check)
    print('upper values', upper_values, 'lower values', lower_values,
          'junction', Junction(lower_values[0], lower_level))

    anchor_upper_values, _ = UpperReduced((1 - SHRINK) * budget)
    anchor_lower_values, _, _ = LowerReduced((1 - SHRINK) * budget)
    anchors = {1: UpperAnchor(anchor_upper_values), -1: LowerAnchor(anchor_lower_values)}

    trace = EvaluationTrace(limits=LIMITS)
    Demand, Target, Fit, constraints = BuildBasis(anchors)
    estimate = np.array([0.0, intercept, slope])
    problem = SensitivityProblem(criterion=Fit, target=Target, estimate=estimate,
                                 constraints=constraints, trace=trace)
    nests = NestsBaseline(problem, estimate)
    print('nests', nests)

    solved = {}
    for sign in (1, -1):
        _, branch_target, branch_fit, branch_constraints = BuildBasis(anchors, sign)
        branch = SensitivityProblem(criterion=branch_fit, target=branch_target,
                                    estimate=estimate, constraints=branch_constraints,
                                    trace=trace)
        solver = GeneralizedReducedGradient(branch)
        solver.Solve(max_iterations=200)
        found = np.asarray(solver.positive_primitive if sign > 0
                           else solver.negative_primitive, dtype=float)
        solved[sign] = np.array([sign * found[0], found[1], found[2]])
        print('solver', sign, solved[sign], solver.iteration, solver.converged)

    anchor_upper = np.array([1.0, intercept, slope])
    anchor_lower = np.array([-1.0, intercept, slope])
    extras = {1: {'anchor_values_at_prices': anchor_upper_values.tolist()},
              -1: {'anchor_values_at_prices': anchor_lower_values.tolist()}}
    candidates = {'upper': [(anchor_upper, 'anchor', 1), (solved[1], 'solver', 1)],
                  'lower': [(anchor_lower, 'anchor', -1), (solved[-1], 'solver', -1)]}
    described = {side: [Describe(problem, Demand, x, source, extras[sign])
                        for x, source, sign in pairs]
                 for side, pairs in candidates.items()}
    for side, items in described.items():
        for item in items:
            print(side, item['source'], item['primitive']['weight'], item['target'],
                  item['criterion'], item['slack'], item['feasible'])

    def Best(side, sign):
        feasible = [item for item in described[side] if item['feasible']]
        return max(feasible, key=lambda item: sign * item['target'])

    upper, lower = Best('upper', 1), Best('lower', -1)

    def Point(item):
        p = item['primitive']
        return np.array([p['weight'], p['intercept'], p['slope']])

    def Leg(certificate, anchor):
        leg = [estimate + t * (anchor - estimate) for t in np.linspace(0, 1, PATH_POINTS)[1:]]
        return leg if np.array_equal(certificate, anchor) else leg + [certificate]

    points = ([estimate] + Leg(Point(lower), anchor_lower))[::-1] + Leg(Point(upper), anchor_upper)
    path_feasible = all(problem.IsFeasible(x) for x in points)
    print('path feasible', path_feasible)
    deformation = DeformationPath(problem, points, curve=lambda x: Demand(x).Value(PLOT_GRID))
    WriteAnimation(HERE / 'animation.gif', deformation, PLOT_GRID,
                   labels={'x': 'relative price $p$', 'y': 'demand $D(p)$',
                           'object': 'deadweight loss'},
                   markers={'x': PRICES, 'y': SHARES, 'label': 'observed shares'},
                   title='Demand along the deformation path',
                   figure=HERE / 'certificates.png')

    snapshot = trace.Snapshot()
    certificates = {
        'baseline': {'target': float(problem.target_at_estimate), 'criterion': budget,
                     'feasible': bool(problem.IsFeasible(estimate)),
                     'primitive': {'weight': 0.0, 'intercept': intercept, 'slope': slope}},
        'lower': lower,
        'upper': upper,
        'analytic': {'lower': analytic_lower, 'upper': analytic_upper},
        'spend': snapshot,
        'limits': LIMITS,
        'reduced': {'upper_values_at_prices': upper_values.tolist(),
                    'lower_values_at_prices': lower_values.tolist(),
                    'lower_total_cost_at_p2': lower_level,
                    'lower_junction': Junction(lower_values[0], lower_level),
                    'lower_full_check': lower_check,
                    'competitive_price': P_C, 'avg_cost_intercept': GAMMA,
                    'avg_cost_slope': DELTA},
        'nests_baseline': nests,
        'path_feasible': path_feasible,
        'candidates': described,
    }
    with open(HERE / 'certificates.json', 'w') as stream:
        json.dump(certificates, stream, indent=2,
                  default=lambda v: v.tolist() if hasattr(v, 'tolist') else str(v))
        stream.write('\n')
    print(json.dumps(snapshot))


if __name__ == '__main__':
    Main()
