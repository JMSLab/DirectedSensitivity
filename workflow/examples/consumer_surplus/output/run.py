"""Directed sensitivity of the consumer surplus loss to the isoelastic demand assumption."""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parent
sys.path.insert(0, str(EXAMPLE))

from consumer_surplus import ConsumerSurplusLoss, PRICE_RANGE  # noqa: E402
from estimator import Criterion, Isoelastic, LoadDemand, SLOPE_FLOOR  # noqa: E402

from DirectedSensitivity import (DeformationPath, EvaluationTrace,  # noqa: E402
                                 GeneralizedReducedGradient, NestsBaseline,
                                 SensitivityProblem, WriteAnimation)

LOW, HIGH = PRICE_RANGE
SHRINK = 1e-6
FLAT_SLOPE = 1e-8
DROP_WIDTH = 1e-4
GRID = np.linspace(LOW, HIGH, 2001)
PATH_POINTS = 41


def LoadEstimate():
    data = np.genfromtxt(EXAMPLE / 'estimates.csv', delimiter=',', names=True)
    return float(data['intercept']), float(data['slope'])


def AnalyticBounds(prices, shares, budget):
    widths = np.diff(prices)
    radius = np.sqrt(len(shares) * budget)
    upper = widths @ shares[:-1] + np.linalg.norm(widths) * radius
    lower = widths @ shares[1:] - np.linalg.norm(widths) * radius
    return float(lower), float(upper)


def AnchorKnots(prices, shares, budget, sign):
    widths = np.diff(prices)
    radius = (1 - SHRINK) * np.sqrt(len(shares) * budget)
    values = shares.astype(float).copy()
    if sign > 0:
        values[:-1] += radius * widths / np.linalg.norm(widths)
    else:
        values[1:] -= radius * widths / np.linalg.norm(widths)
    nodes, levels = [], []
    for i in range(len(prices) - 1):
        if sign > 0:
            nodes += [prices[i], prices[i + 1] - DROP_WIDTH]
            levels += [values[i], values[i] - FLAT_SLOPE * (widths[i] - DROP_WIDTH)]
        else:
            nodes += [prices[i], prices[i] + DROP_WIDTH]
            levels += [values[i], values[i + 1] + FLAT_SLOPE * (widths[i] - DROP_WIDTH)]
    nodes.append(prices[-1])
    levels.append(values[-1])
    return np.array(nodes), np.array(levels), values


def BuildBasis(prices, shares, budget, side=None):
    """The union basis model in a signed weight, or one branch of it when ``side`` is given."""
    anchors = {sign: AnchorKnots(prices, shares, budget, sign) for sign in (1, -1)}

    def Split(weight):
        if side is None:
            return (1 if weight >= 0 else -1), abs(weight)
        return side, weight

    def Demand(x):
        weight, intercept, slope = x
        sign, mix = Split(weight)
        nodes, levels, _ = anchors[sign]
        working = Isoelastic(intercept, slope)
        return lambda price: mix * np.interp(price, nodes, levels) + (1 - mix) * working(price)

    def Target(x):
        demand = Demand(x)
        nodes = anchors[Split(x[0])[0]][0]
        return sum(ConsumerSurplusLoss(demand, lo, hi) for lo, hi in zip(nodes[:-1], nodes[1:]))

    def Fit(x):
        return Criterion(Demand(x), prices, shares)

    constraints = {
        'weight': lambda x: 1 - abs(x[0]),
        'slope': lambda x: x[2] - SLOPE_FLOOR,
        'intercept': lambda x: x[2] * np.log(LOW) - x[1],
        'decreasing': lambda x: float(np.min(-np.diff(Demand(x)(GRID)))),
        'ceiling': lambda x: float(1 - Demand(x)(LOW)),
        'floor': lambda x: float(Demand(x)(HIGH)),
    }
    if side is not None:
        constraints['branch'] = lambda x: x[0]
    return Demand, Target, Fit, constraints, anchors


def Describe(problem, x, anchors, source):
    weight, intercept, slope = (float(v) for v in x)
    nodes, levels, values = anchors[1 if weight >= 0 else -1]
    return {
        'target': float(problem.target(x)),
        'criterion': float(problem.criterion(x)),
        'slack': float(problem.Slack(x)),
        'deviation': float(problem.SignedDeviation(x)),
        'feasible': bool(problem.IsFeasible(x)),
        'primitive': {'weight': weight, 'intercept': intercept, 'slope': slope,
                      'anchor_values_at_prices': values.tolist(),
                      'anchor_nodes': nodes.tolist(), 'anchor_levels': levels.tolist()},
        'constraints': {k: float(v) for k, v in problem.ConstraintReport(x).items()},
        'source': source,
    }


def Main():
    prices, shares = LoadDemand()
    intercept, slope = LoadEstimate()
    baseline_demand = Isoelastic(intercept, slope)
    budget = Criterion(baseline_demand, prices, shares)

    trace = EvaluationTrace(limits=None)
    Demand, Target, Fit, constraints, anchors = BuildBasis(prices, shares, budget)
    estimate = np.array([0.0, intercept, slope])
    problem = SensitivityProblem(criterion=Fit, target=Target, estimate=estimate,
                                 constraints=constraints, trace=trace)
    nests = NestsBaseline(problem, estimate)
    print('nests', nests)

    analytic_lower, analytic_upper = AnalyticBounds(prices, shares, budget)
    print('analytic', analytic_lower, analytic_upper)

    anchor_upper = np.array([1.0, intercept, slope])
    anchor_lower = np.array([-1.0, intercept, slope])

    solved = {}
    for sign in (1, -1):
        _, branch_target, branch_fit, branch_constraints, _ = BuildBasis(prices, shares, budget, sign)
        branch = SensitivityProblem(criterion=branch_fit, target=branch_target, estimate=estimate,
                                    constraints=branch_constraints, trace=trace)
        solver = GeneralizedReducedGradient(branch)
        solver.Solve()
        found = solver.positive_primitive if sign > 0 else solver.negative_primitive
        found = np.asarray(found, dtype=float)
        solved[sign] = np.array([sign * found[0], found[1], found[2]])
    solved_upper, solved_lower = solved[1], solved[-1]

    candidates = {
        'upper': [(anchor_upper, 'anchor'), (solved_upper, 'solver')],
        'lower': [(anchor_lower, 'anchor'), (solved_lower, 'solver')],
    }
    described = {side: [Describe(problem, x, anchors, source) for x, source in pairs]
                 for side, pairs in candidates.items()}
    for side, items in described.items():
        for item in items:
            print(side, item['source'], item['primitive']['weight'], item['target'],
                  item['criterion'], item['slack'], item['feasible'])

    def Best(side, sign):
        feasible = [item for item in described[side] if item['feasible']]
        return max(feasible, key=lambda item: sign * item['target'])

    upper, lower = Best('upper', 1), Best('lower', -1)
    baseline = {'target': float(problem.target(estimate)), 'criterion': budget,
                'feasible': bool(problem.IsFeasible(estimate)),
                'primitive': {'weight': 0.0, 'intercept': intercept, 'slope': slope}}

    def Point(item):
        p = item['primitive']
        return np.array([p['weight'], p['intercept'], p['slope']])

    def Leg(start, end):
        return [start + t * (end - start) for t in np.linspace(0, 1, PATH_POINTS)[1:]]

    def Branch(certificate, anchor):
        if np.array_equal(certificate, anchor):
            return Leg(estimate, anchor)
        return Leg(estimate, anchor) + [certificate]

    points = ([estimate] + Branch(Point(lower), anchor_lower))[::-1] + Branch(Point(upper), anchor_upper)
    path_feasible = all(problem.IsFeasible(x) for x in points)
    print('path feasible', path_feasible)

    deformation = DeformationPath(problem, points, curve=lambda x: Demand(x)(GRID))
    WriteAnimation(HERE / 'animation.gif', deformation, GRID,
                   labels={'x': 'price $p$', 'y': 'demand $D(p)$'},
                   markers={'x': prices, 'y': shares},
                   title='Demand along the deformation path',
                   figure=HERE / 'certificates.png')

    snapshot = trace.Snapshot()
    certificates = {
        'baseline': baseline,
        'lower': lower,
        'upper': upper,
        'analytic': {'lower': analytic_lower, 'upper': analytic_upper},
        'spend': snapshot,
        'limits': {'criterion': None, 'target': None},
        'nests_baseline': nests,
        'path_feasible': path_feasible,
        'candidates': described,
    }
    with open(HERE / 'certificates.json', 'w') as stream:
        json.dump(certificates, stream, indent=2, default=lambda v: v.tolist() if hasattr(v, 'tolist') else str(v))
        stream.write('\n')
    print(json.dumps(snapshot))


if __name__ == '__main__':
    Main()
