"""Surrogate-guided search for the K2016 certificates in the bump basis; all simulations logged."""
import sys
import numpy as np
from scipy.optimize import minimize
from scipy.stats import qmc

sys.path.insert(0, '.')
from basis import Evaluator, Pad

ESTIMATE = Pad([0.288509110693, -2.98192892931e-05, 9.50846758198e-05, 1.0, 1.0,
                np.log(500.0), np.log(500.0)])
ACTIVE = (0, 1, 2, 8, 10)
CENTER = np.array([ESTIMATE[0], 0.0, 0.0, 0.8, np.log(1.5)])
SCALE = np.array([0.002, 0.04, 0.08, 0.8, 0.6])


def Embed(z):
    point = ESTIMATE.copy()
    point[list(ACTIVE)] = CENTER + np.asarray(z) * SCALE
    return point


def Coordinates(point):
    return (np.asarray(Pad(point))[list(ACTIVE)] - CENTER) / SCALE


def Features(z):
    z = np.atleast_2d(z)
    d = z.shape[1]
    return np.column_stack([np.ones(len(z))] + [z[:, i] for i in range(d)]
                           + [z[:, i] * z[:, j] for i in range(d) for j in range(i, d)])


def Records(evaluator, radius=2.5):
    rows = []
    for record in evaluator.records.values():
        point = Pad(record['point'])
        if point[3] != 1.0 or point[4] != 1.0 or point[7] != 0.0:
            continue
        if point[8] == 0.0:
            point = point.copy(); point[10] = CENTER[4]
        z = Coordinates(point)
        if np.max(np.abs(z)) <= radius:
            rows.append((z, record['criterion'], record['target']))
    return rows


def Fit(rows, centre=None, width=None):
    z = np.array([r[0] for r in rows])
    weights = np.ones(len(z)) if centre is None else np.exp(
        -0.5 * np.sum((z - centre) ** 2, axis=1) / width ** 2)
    features = Features(z) * np.sqrt(weights)[:, None]
    models = []
    for column in (1, 2):
        y = np.array([r[column] for r in rows]) * np.sqrt(weights)
        coefficients = np.linalg.lstsq(features, y, rcond=None)[0]
        models.append(coefficients)
        residual = (y - features @ coefficients) / np.sqrt(weights)
        print('fit', column, 'rmse', float(np.sqrt(np.average(residual ** 2, weights=weights))))
    return models


def Design(evaluator, count, seed):
    unit = qmc.LatinHypercube(d=len(ACTIVE), seed=seed).random(count)
    z = -1.0 + 2.0 * unit
    z[:, 3] = -1.0 + 3.0 * unit[:, 3] * (2.5 / 2.4)
    return evaluator.Batch([Embed(row) for row in z])


def SurrogateOptimum(models, sign, budget, margin, start, bounds):
    criterion, target = models
    result = minimize(lambda z: -sign * float(Features(z) @ target), start, method='SLSQP',
                      bounds=bounds,
                      constraints=[{'type': 'ineq', 'fun': lambda z: budget - margin
                                    - float(Features(z) @ criterion)}])
    return result.x


if __name__ == '__main__':
    evaluator = Evaluator()
    Design(evaluator, 40, 11)
    rows = Records(evaluator)
    print('rows', len(rows))
    Fit(rows)


BUDGET = 2.797250440267154
LOWER_BOUNDS = np.array([-np.inf, -np.inf, -np.inf, -1.0, -np.inf])


def Local(rows, centre, width):
    z = np.array([r[0] for r in rows]) - centre
    weights = np.exp(-0.5 * np.sum(z ** 2, axis=1) / width ** 2)
    features = np.column_stack([np.ones(len(z)), z]) * np.sqrt(weights)[:, None]
    ridge = 1e-3 * np.eye(features.shape[1]); ridge[0, 0] = 0.0
    gradients = []
    for column in (1, 2):
        y = np.array([r[column] for r in rows]) * np.sqrt(weights)
        coefficients = np.linalg.solve(features.T @ features + ridge, features.T @ y)
        gradients.append(coefficients[1:])
    return gradients


def Climb(evaluator, sign, rounds, seed, width=0.6, step=0.3, batch=4):
    rng = np.random.default_rng(seed)
    feasible = [(r['target'], Coordinates(r['point'])) for r in evaluator.records.values()
                if r['criterion'] <= BUDGET and Pad(r['point'])[3] == 1.0
                and Pad(r['point'])[7] == 0.0]
    value, current = max(feasible, key=lambda item: sign * item[0])
    for round_ in range(rounds):
        rows = Records(evaluator, radius=np.inf)
        criterion_gradient, target_gradient = Local(rows, current, width)
        unit_q = criterion_gradient / np.linalg.norm(criterion_gradient)
        ascent = sign * target_gradient
        tangent = ascent - (ascent @ unit_q) * unit_q
        tangent /= np.linalg.norm(tangent)
        candidates = []
        for k in range(batch):
            length = step * (0.5 + k / (batch - 1))
            pull = rng.uniform(0.0, 0.3) * length
            z = current + length * tangent - pull * unit_q + rng.normal(0, 0.03, 5)
            z = np.maximum(z, LOWER_BOUNDS)
            candidates.append(z)
        records = evaluator.Batch([Embed(z) for z in candidates])
        better = [(r['target'], z) for z, r in zip(candidates, records)
                  if r['criterion'] <= BUDGET and sign * (r['target'] - value) > 0]
        if better:
            value, current = max(better, key=lambda item: sign * item[0])
            step *= 1.3
        else:
            step *= 0.6
        print(f'{sign:+d} round {round_} c {value:.3f} step {step:.3f} n {len(evaluator)}',
              np.round(Embed(current)[list(ACTIVE)], 4), flush=True)
    return value, current


EFFICIENT = np.array([-735.76272189, 423.97886236, 267.65786572]) * np.array([0.002, 0.04, 0.08])
CURVATURE = 417.32


def Trade(sign, bump, scale, criterion_cost, jitter):
    """The estimate with an opponent bump, moved along the working model's efficient direction
    for the target by the step whose quadratic-model criterion cost is ``criterion_cost``."""
    point = ESTIMATE.copy()
    point[8], point[10] = bump, np.log(scale)
    point[:3] += sign * np.sqrt(criterion_cost / CURVATURE) * EFFICIENT + jitter
    return point


def Sample(evaluator, sign, count, bumps, scales, costs, seed):
    rng = np.random.default_rng(seed)
    points = [Trade(sign, rng.choice(bumps), rng.choice(scales), rng.uniform(*costs),
                    rng.normal(0.0, [5e-5, 5e-4, 1e-3])) for _ in range(count)]
    records = evaluator.Batch(points)
    feasible = [(r['target'], p) for p, r in zip(points, records) if r['criterion'] <= BUDGET]
    print(sign, 'feasible', len(feasible), 'of', count,
          'best', max((sign * f[0] for f in feasible), default=None), flush=True)
    return feasible
