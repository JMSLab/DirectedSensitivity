"""Directed sensitivity of the ECU welfare effect to the censored normal altruism distribution."""
import json
import os
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..')))

import reduction as rd  # noqa: E402
import specification as spec  # noqa: E402
from DirectedSensitivity import (DeformationPath, EvaluationTrace,  # noqa: E402
                                 GeneralizedReducedGradient, NestsBaseline,
                                 SensitivityProblem, WriteAnimation)

LIMITS = {'criterion': 20_000, 'target': 50_000}
WIDTH = 1e-4
SIGMA_FLOOR = 1e-6
SOLVER_ITERATIONS = 15
PATH_POINTS = 25
GRID = np.linspace(0.0, rd.A_MAX, 1401)
LAR = spec.Baseline().lar


class Mixture(spec.AltruismDistribution):
    """w times an anchor (an atom at zero plus narrow uniform densities) plus 1 - w times
    the censored normal F_{mu, sigma}."""

    def __init__(self, weight, anchor, mu, sigma):
        self.weight, self.anchor = float(weight), anchor
        self.base = spec.CensoredNormal(mu, sigma, rd.A_MAX)
        self.a_max = rd.A_MAX

    @property
    def Atom(self):
        atom = (1.0 - self.weight) * self.base.Atom
        return atom + (self.weight * self.anchor['atom'] if self.weight else 0.0)

    def Density(self, a):
        a = np.asarray(a, dtype=float)
        value = (1.0 - self.weight) * self.base.Density(a)
        if self.weight:
            for left, right, mass in self.anchor['pieces']:
                value = value + self.weight * mass / (right - left) * ((a >= left) & (a <= right))
        return value

    def _Continuous(self, psi, lower, upper, breakpoints=()):
        total = (1.0 - self.weight) * self.base._Continuous(psi, lower, upper, breakpoints) if self.weight < 1 else 0.0
        if not self.weight:
            return total
        for left, right, mass in self.anchor['pieces']:
            lo, hi = max(lower, left), min(upper, right)
            if hi <= lo:
                continue
            edges = [lo] + sorted(b for b in breakpoints if lo < b < hi) + [hi]
            piece = sum(quad(lambda a: float(psi(a)), x, y, epsabs=1e-14, epsrel=1e-12, limit=200)[0]
                        for x, y in zip(edges[:-1], edges[1:]))
            total += self.weight * mass / (right - left) * piece
        return total

    def CdfCurve(self, grid):
        base = np.where(grid >= self.a_max, 1.0, norm.cdf(grid, self.base.mu, self.base.sigma) + self.base._tail)
        curve = (1.0 - self.weight) * base
        if self.weight:
            curve = curve + self.weight * self.anchor['atom']
            for left, right, mass in self.anchor['pieces']:
                curve = curve + self.weight * mass * np.clip((grid - left) / (right - left), 0.0, 1.0)
        return curve


def Anchor(z, sign, width=WIDTH):
    """A primitive in N whose region statistics approach z as the width shrinks: narrow
    uniform densities at the points the extremal face selects."""
    p = np.clip(z[:rd.REGIONS], 0.0, None)
    p = p / p.sum()
    ratio = np.where(p > 0, z[rd.REGIONS:] / np.where(p > 0, p, 1.0), 0.0)
    pieces = []
    for j in range(1, rd.REGIONS):
        if p[j] < 1e-12:
            continue
        left, right = rd.EDGES[j], rd.EDGES[j + 1]
        if j == rd.TOP:
            spots = [(right - width if sign > 0 else left, p[j])]
        elif j == rd.CLIP and sign > 0:
            share = np.clip((1.0 - ratio[j]) / (1.0 - rd.H_LOW[j]), 0.0, 1.0)
            spots = [(left, p[j] * share), (right - width, p[j] * (1.0 - share))]
        else:
            spots = [(Locate(ratio[j], left, right, width), p[j])]
        pieces += [(start, start + width, mass) for start, mass in spots if mass > 1e-14]
    return {'atom': float(p[0]), 'pieces': pieces, 'width': width, 'sign': sign}


def Locate(ratio, left, right, width):
    """Start of a width-wide interval in (left, right] on which h* averages ``ratio``."""
    def mean(start):
        return quad(lambda a: float(rd.Home(a)), start, start + width, epsabs=1e-15)[0] / width - ratio
    low, high = left, right - width
    if mean(low) >= 0:
        return low
    if mean(high) <= 0:
        return high
    return brentq(mean, low, high, xtol=1e-14)


def Main():
    os.makedirs(HERE, exist_ok=True)
    trace = EvaluationTrace(limits=LIMITS)

    def WorkingPrimitive(x):
        return spec.Primitive(spec.BENCHMARK, lar=LAR, ecu=spec.CensoredNormal(x[0], x[1], rd.A_MAX))

    working_criterion = trace.Watch('criterion', lambda x: spec.Criterion(WorkingPrimitive(x)))
    published = spec.Baseline().ecu
    start = np.array([published.mu, published.sigma])
    published_criterion = working_criterion(start)
    fit = minimize(working_criterion, start, method='Nelder-Mead',
                   options={'xatol': 1e-8, 'fatol': 1e-11, 'maxfev': 600})
    mu, sigma = (fit.x if fit.fun < published_criterion else start)
    rebaselined = bool(fit.fun < published_criterion)
    print('published', published_criterion, 'refit', fit.fun, fit.x, rebaselined)

    baseline_primitive = WorkingPrimitive([mu, sigma])
    fixed = rd.FixedCriterion(baseline_primitive)
    p, h, u = rd.Statistics(baseline_primitive.ecu)
    z_baseline = np.r_[p, h]
    budget_reduced = rd.ReducedCriterion(z_baseline, fixed)
    _, minimum_reduced = rd.MinimumCriterion(fixed)
    extremes = {sign: rd.Extremum(sign, fixed, budget_reduced, z_baseline) for sign in (1, -1)}
    for sign, item in extremes.items():
        print('reduced', sign, item['primal'], item['dual'], item['gap'], item['criterion'] - budget_reduced)
    anchors = {sign: Anchor(extremes[sign]['z'], sign) for sign in (1, -1)}

    def ToPrimitive(x, side=None):
        weight, mu_x, sigma_x = (float(v) for v in x)
        sign = side if side is not None else (1 if weight >= 0 else -1)
        share = weight if side is not None else abs(weight)
        return spec.Primitive(spec.BENCHMARK, lar=LAR,
                              ecu=Mixture(share, anchors[sign], mu_x, sigma_x))

    constraints = dict(spec.ConstraintSet('ecu'))
    point_constraints = {'weight': lambda x: 1.0 - abs(x[0]), 'sigma': lambda x: x[2] - SIGMA_FLOOR}

    def Build(side=None):
        wrapped = {name: (lambda x, f=f: f(ToPrimitive(x, side))) for name, f in constraints.items()}
        wrapped.update(point_constraints)
        if side is not None:
            wrapped['branch'] = lambda x: x[0]
        return SensitivityProblem(criterion=lambda x: spec.Criterion(ToPrimitive(x, side)),
                                  target=lambda x: spec.TargetECU(ToPrimitive(x, side)),
                                  estimate=np.array([0.0, mu, sigma]),
                                  constraints=wrapped, trace=trace)

    problem = Build()
    estimate = problem.estimate
    budget = problem.criterion_at_estimate
    nests = NestsBaseline(problem, estimate)
    working_gap = abs(budget - working_criterion([mu, sigma]))
    print('budget', budget, 'reduced', budget_reduced, 'nests', nests, 'working gap', working_gap)

    def SegmentWeight(sign):
        """Largest w with q(w * anchor + (1 - w) * baseline) <= budget, from the quadratic in w."""
        end = problem.criterion(np.array([sign * 1.0, mu, sigma]))
        middle = problem.criterion(np.array([sign * 0.5, mu, sigma]))
        curvature = 2.0 * (end + budget - 2.0 * middle)
        drift = end - budget
        if drift <= 0:
            return 1.0
        weight = 1.0 - drift / (curvature / 2.0)
        for _ in range(60):
            if problem.IsFeasible(np.array([sign * weight, mu, sigma])):
                return weight
            weight -= 1e-9 + 1e-3 * (1.0 - weight)
        return 0.0

    segment = {sign: SegmentWeight(sign) for sign in (1, -1)}
    print('segment weights', segment)

    solved = {}
    for sign in (1, -1):
        branch = Build(sign)
        solver = GeneralizedReducedGradient(branch, step_size=0.05)
        solver.StartPhase(float(sign))
        solver.Solve(max_iterations=SOLVER_ITERATIONS)
        found = np.asarray(solver.positive_primitive if sign > 0 else solver.negative_primitive, dtype=float)
        solved[sign] = {'point': np.array([sign * found[0], found[1], found[2]]),
                        'iterations': solver.iteration, 'converged': bool(solver.converged)}
        print('solver', sign, solved[sign], branch.target(found))

    def Describe(x, source):
        x = np.asarray(x, dtype=float)
        sign = 1 if x[0] >= 0 else -1
        anchor = anchors[sign]
        return {'target': float(problem.target(x)), 'criterion': float(problem.criterion(x)),
                'slack': float(problem.Slack(x)), 'deviation': float(problem.SignedDeviation(x)),
                'feasible': bool(problem.IsFeasible(x)),
                'primitive': {'weight': float(x[0]), 'mu': float(x[1]), 'sigma': float(x[2]),
                              'anchor_atom_at_zero': anchor['atom'],
                              'anchor_pieces': [list(map(float, piece)) for piece in anchor['pieces']]},
                'constraints': {k: float(v) for k, v in problem.ConstraintReport(x).items()},
                'source': source}

    candidates = {}
    for name, sign in (('upper', 1), ('lower', -1)):
        candidates[name] = [Describe([sign * segment[sign], mu, sigma], 'anchor'),
                            Describe(solved[sign]['point'], 'solver')]

    def Best(name, sign):
        feasible = [item for item in candidates[name] if item['feasible']]
        return max(feasible, key=lambda item: sign * item['target'])

    upper, lower = Best('upper', 1), Best('lower', -1)
    for name in ('lower', 'upper'):
        for item in candidates[name]:
            print(name, item['source'], item['primitive']['weight'], item['target'], item['slack'], item['feasible'])

    def Point(item):
        primitive = item['primitive']
        return np.array([primitive['weight'], primitive['mu'], primitive['sigma']])

    def Leg(end):
        return [estimate + t * (end - estimate) for t in np.linspace(0.0, 1.0, PATH_POINTS)[1:]]

    points = Leg(Point(lower))[::-1] + [estimate] + Leg(Point(upper))
    path_feasible = all(problem.IsFeasible(x) for x in points)
    print('path feasible', path_feasible)
    deformation = DeformationPath(problem, points,
                                  curve=lambda x: ToPrimitive(x).ecu.CdfCurve(GRID))
    WriteAnimation(os.path.join(HERE, 'animation.gif'), deformation, GRID,
                   labels={'x': 'altruism $a$', 'y': 'distribution function $F(a)$',
                           'object': 'c', 'lower': 'lower certificate',
                           'upper': 'upper certificate'},
                   title='ECU altruism distribution along the deformation path',
                   figure=os.path.join(HERE, 'certificates.png'))

    snapshot = trace.Snapshot()
    certificates = {
        'baseline': {'target': float(problem.target_at_estimate), 'criterion': float(budget),
                     'feasible': bool(problem.IsFeasible(estimate)),
                     'primitive': {'weight': 0.0, 'mu': float(mu), 'sigma': float(sigma)}},
        'lower': lower,
        'upper': upper,
        'analytic': {'lower': extremes[-1]['primal'], 'upper': extremes[1]['primal']},
        'spend': snapshot,
        'limits': LIMITS,
        'published_baseline': {'mu': float(published.mu), 'sigma': float(published.sigma),
                               'criterion': float(published_criterion)},
        'rebaselined': rebaselined,
        'nests_baseline': nests,
        'path_feasible': path_feasible,
        'reduced': {
            'edges': rd.EDGES.tolist(), 'h_low': rd.H_LOW.tolist(), 'h_high': rd.H_HIGH.tolist(),
            'u_edges': rd.U_EDGE.tolist(), 'budget': budget_reduced, 'fixed': fixed,
            'minimum_criterion': minimum_reduced, 'baseline_z': z_baseline.tolist(),
            **{name: {'z': extremes[sign]['z'].tolist(), 'primal': extremes[sign]['primal'],
                      'dual': extremes[sign]['dual'], 'gap': extremes[sign]['gap'],
                      'multiplier': extremes[sign]['multiplier'],
                      'objective': rd.ObjectiveVector(sign).tolist()}
               for name, sign in (('upper', 1), ('lower', -1))}},
        'segment_weights': {'upper': segment[1], 'lower': segment[-1]},
        'solver': {name: {'iterations': solved[sign]['iterations'], 'converged': solved[sign]['converged']}
                   for name, sign in (('upper', 1), ('lower', -1))},
        'candidates': candidates,
    }
    with open(os.path.join(HERE, 'certificates.json'), 'w') as stream:
        json.dump(certificates, stream, indent=2, default=lambda v: v.tolist() if hasattr(v, 'tolist') else str(v))
        stream.write('\n')
    print(json.dumps(snapshot))


if __name__ == '__main__':
    Main()
