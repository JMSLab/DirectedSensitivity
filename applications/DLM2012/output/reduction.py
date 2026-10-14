"""Finite-dimensional reduction of the ECU problem: region masses and home-presence integrals."""
import os
import sys

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import brentq, minimize, nnls

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import specification as spec  # noqa: E402

CHARITY = 'ecu'
P = spec.BENCHMARK
H0, R, ETA, GI = P.h0, P.r, P.eta, P.gi
S = P.s_ecu
A_MAX = float(spec.A_MAX[CHARITY])
A_MINMIN, A_MIN, A_PLUS, A_ZERO = spec.Thresholds(P, CHARITY)
EDGES = np.array([0.0, A_MINMIN, A_ZERO, A_MIN, A_PLUS, GI + 20.0, GI + 50.0, A_MAX])
REGIONS = len(EDGES) - 1
CLIP = 5
TOP = 6
U_STAR = (1.0 - H0) / ETA
ECU = spec.BLOCK[CHARITY]
NAMES = [name for block, name in spec.MOMENT_ORDER if block == CHARITY]


def Utility(a):
    return spec.VisitUtility(a, P, CHARITY)


def Home(a):
    return spec.HomeStar(a, P, CHARITY)


U_EDGE = np.array([float(Utility(a)) for a in EDGES])
U_EDGE[0] = -S * spec.GS
H_EDGE = np.array([float(Home(a)) for a in EDGES])
A_STAR = brentq(lambda a: float(Utility(a)) - U_STAR, EDGES[CLIP], EDGES[CLIP + 1], xtol=1e-14)


def Statistics(distribution):
    """Region masses p, home-presence integrals H, and utility integrals U."""
    kinks = tuple(EDGES[1:-1]) + (A_STAR,)
    p, h, u = np.zeros(REGIONS), np.zeros(REGIONS), np.zeros(REGIONS)
    for j in range(REGIONS):
        lower = -np.inf if j == 0 else EDGES[j]
        upper = EDGES[j + 1]
        p[j] = distribution.Integrate(lambda a: np.ones_like(np.asarray(a, dtype=float)), lower, upper, kinks)
        h[j] = distribution.Integrate(Home, lower, upper, kinks)
        u[j] = distribution.Integrate(Utility, lower, upper, kinks)
    return p, h, u


def MomentMap():
    """ECU moments as b + A z, with z = (p_0..p_6, H_0..H_6)."""
    rows = {}

    def row(p=(), h=(), constant=0.0):
        vector = np.zeros(2 * REGIONS)
        for j, weight in p:
            vector[j] += weight
        for j, weight in h:
            vector[REGIONS + j] += weight
        return constant, vector

    survives = [(j, 1.0) for j in range(1, REGIONS)]
    rows['PH_NF'] = row(constant=H0)
    rows['PH_F'] = row(h=[(j, R) for j in range(REGIONS)], constant=(1 - R) * H0)
    rows['PH_OO'] = row(h=[(j, R) for j in range(2, REGIONS)], constant=(1 - R) * H0)
    rows['PG_NF'] = row(p=[(j, H0) for j, _ in survives])
    rows['PG_F'] = row(p=[(j, (1 - R) * H0) for j, _ in survives], h=[(j, R) for j in range(1, REGIONS)])
    rows['PG_OO'] = row(p=[(j, (1 - R) * H0) for j, _ in survives], h=[(j, R) for j in range(2, REGIONS)])
    rows['POO_OO'] = row(p=[(0, R * H0), (1, R * H0)])
    bins = {'0gs': ((1, 2), (2,)), 'gs': ((3,), (3,)), 'int': ((4,), (4,)),
            'int1': ((5,), (5,)), 'int2': ((6,), (6,))}
    for name, (inside, opted) in bins.items():
        rows[f'PG_NF_{name}'] = row(p=[(j, H0) for j in inside])
        rows[f'PG_F_{name}'] = row(p=[(j, (1 - R) * H0) for j in inside], h=[(j, R) for j in inside])
        rows[f'PG_OO_{name}'] = row(p=[(j, (1 - R) * H0) for j in inside], h=[(j, R) for j in opted])
    b = np.array([rows[name][0] for name in NAMES])
    A = np.array([rows[name][1] for name in NAMES])
    return b, A


B_MAP, A_MAP = MomentMap()
WEIGHTS = np.diag(spec.WEIGHT)[ECU]
TARGET_MOMENTS = spec.EMPIRICAL[ECU]


def FixedCriterion(primitive):
    """The part of the criterion that does not depend on the ECU distribution."""
    gap = spec.EMPIRICAL - spec.Moments(primitive)
    weights = np.diag(spec.WEIGHT)
    return float(np.sum((weights * gap ** 2)[~ECU]))


def ReducedCriterion(z, fixed):
    gap = TARGET_MOMENTS - B_MAP - A_MAP @ z
    return fixed + float(np.sum(WEIGHTS * gap ** 2))


def Bounds():
    """Per-region bounds on H_j / p_j, from the range of h* on the closure of the region."""
    lower, upper = H_EDGE[:-1].copy(), H_EDGE[1:].copy()
    lower[0] = upper[0] = H0 - ETA * S * spec.GS
    upper[CLIP] = 1.0
    lower[TOP] = upper[TOP] = 1.0
    return lower, upper


H_LOW, H_HIGH = Bounds()
U_TOP = (U_EDGE[TOP], U_EDGE[TOP + 1])
CHORD_SLOPE = (U_EDGE[CLIP + 1] - U_EDGE[CLIP]) / (1.0 - H_EDGE[CLIP])


def ObjectiveVector(sign):
    """c = l'z on the extremal face for direction ``sign``: U_j = (H_j - h0 p_j)/eta below the
    clip, the chord of the clipped region for the maximum, and the top region at an end."""
    l = np.zeros(2 * REGIONS)
    for j in range(REGIONS):
        l[j] = -H0 / ETA
        l[REGIONS + j] = 1.0 / ETA
    if sign > 0:
        l[CLIP] = U_EDGE[CLIP] - H_EDGE[CLIP] * CHORD_SLOPE
        l[REGIONS + CLIP] = CHORD_SLOPE
        l[TOP], l[REGIONS + TOP] = U_TOP[1], 0.0
    else:
        l[TOP], l[REGIONS + TOP] = U_TOP[0], 0.0
    return H0 * l


def Constraints():
    """Inequalities G z >= 0 and equalities E z = e describing the closure of the reduced set."""
    rows = []
    for j in range(REGIONS):
        if H_LOW[j] == H_HIGH[j]:
            continue
        row = np.zeros(2 * REGIONS)
        row[REGIONS + j], row[j] = 1.0, -H_LOW[j]
        rows.append(row)
        row = np.zeros(2 * REGIONS)
        row[REGIONS + j], row[j] = -1.0, H_HIGH[j]
        rows.append(row)
    for j in range(REGIONS):
        row = np.zeros(2 * REGIONS)
        row[j] = 1.0
        rows.append(row)
    equalities = [np.r_[np.ones(REGIONS), np.zeros(REGIONS)]]
    for j in range(REGIONS):
        if H_LOW[j] == H_HIGH[j]:
            row = np.zeros(2 * REGIONS)
            row[REGIONS + j], row[j] = 1.0, -H_LOW[j]
            equalities.append(row)
    rhs = np.zeros(len(equalities))
    rhs[0] = 1.0
    return np.array(rows), np.array(equalities), rhs


G_INEQ, E_EQ, E_RHS = Constraints()
Z_PARTICULAR = np.linalg.lstsq(E_EQ, E_RHS, rcond=None)[0]
NULL = null_space(E_EQ)


def Quadratic(fixed):
    """q(z0 + N y) = y'Q y - 2 g'y + k."""
    M = A_MAP @ NULL
    resid = TARGET_MOMENTS - B_MAP - A_MAP @ Z_PARTICULAR
    Q = M.T @ (WEIGHTS[:, None] * M)
    g = M.T @ (WEIGHTS * resid)
    k = fixed + float(np.sum(WEIGHTS * resid ** 2))
    return Q, g, k


def MinimumCriterion(fixed):
    """Minimum of the reduced criterion over the closed reduced set."""
    Q, g, k = Quadratic(fixed)
    Gy, Gc = G_INEQ @ NULL, G_INEQ @ Z_PARTICULAR
    result = minimize(lambda y: y @ Q @ y - 2 * g @ y + k, np.zeros(NULL.shape[1]),
                      jac=lambda y: 2 * Q @ y - 2 * g, method='SLSQP',
                      constraints=[{'type': 'ineq', 'fun': lambda y: Gy @ y + Gc, 'jac': lambda y: Gy}],
                      options={'ftol': 1e-15, 'maxiter': 1000})
    return Z_PARTICULAR + NULL @ result.x, float(result.fun)


def Extremum(sign, fixed, budget, start):
    """Primal maximiser of sign * c over the reduced feasible set, and a dual upper bound."""
    Q, g, k = Quadratic(fixed)
    l = sign * ObjectiveVector(sign)
    ly, lc = NULL.T @ l, l @ Z_PARTICULAR
    Gy, Gc = G_INEQ @ NULL, G_INEQ @ Z_PARTICULAR

    def q(y):
        return y @ Q @ y - 2 * g @ y + k

    y0 = np.linalg.lstsq(NULL, start - Z_PARTICULAR, rcond=None)[0]
    result = minimize(lambda y: -(ly @ y + lc), y0, jac=lambda y: -ly, method='SLSQP',
                      constraints=[{'type': 'ineq', 'fun': lambda y: Gy @ y + Gc, 'jac': lambda y: Gy},
                                   {'type': 'ineq', 'fun': lambda y: budget - q(y),
                                    'jac': lambda y: -(2 * Q @ y - 2 * g)}],
                      options={'ftol': 1e-15, 'maxiter': 2000})
    best = None
    for tolerance in (1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        y = Polish(result.x, Q, g, k, ly, Gy, Gc, budget, tolerance)
        if np.min(Gy @ y + Gc) < -1e-13 or q(y) > budget + 1e-10:
            continue
        primal = float(ly @ y + lc)
        multiplier, prices = Multipliers(y, Q, g, ly, Gy, Gc, max(tolerance, 1e-9))
        dual = DualBound(multiplier, prices, Q, g, k, ly, lc, Gy, Gc, budget)
        if best is None or dual - primal < best[0]:
            best = (dual - primal, y, primal, multiplier, prices, dual)
    _, y, primal, multiplier, prices, dual = best
    z = Z_PARTICULAR + NULL @ y
    return {'z': z, 'primal': sign * primal, 'dual': sign * dual, 'gap': dual - primal,
            'criterion': float(q(y)), 'multiplier': multiplier, 'prices': prices,
            'active': np.flatnonzero(np.abs(Gy @ y + Gc) < 1e-10).tolist(),
            'margin': float(np.min(Gy @ y + Gc)), 'status': result.message}


def Polish(y, Q, g, k, ly, Gy, Gc, budget, tolerance=1e-9):
    """Solve the KKT system on the active set found by SLSQP, with q = budget."""
    active = np.flatnonzero(np.abs(Gy @ y + Gc) < tolerance)
    Ga, ca = Gy[active], Gc[active]
    n, m = Q.shape[0], len(active)

    def Stationary(multiplier):
        system = np.block([[2 * multiplier * Q, -Ga.T], [Ga, np.zeros((m, m))]])
        rhs = np.r_[ly + 2 * multiplier * g, -ca]
        return np.linalg.lstsq(system, rhs, rcond=None)[0][:n]

    def Excess(log_multiplier):
        z = Stationary(np.exp(log_multiplier))
        return z @ Q @ z - 2 * g @ z + k - budget

    try:
        log_multiplier = brentq(Excess, -30.0, 30.0, xtol=1e-15)
    except ValueError:
        return y
    polished = Stationary(np.exp(log_multiplier))
    if np.min(Gy @ polished + Gc) < -1e-12:
        return y
    return polished


def Multipliers(y, Q, g, ly, Gy, Gc, tolerance=1e-8):
    """Nonnegative (lambda, nu) with ly = 2 lambda (Q y - g) - Gy' nu on the active set."""
    active = np.flatnonzero(np.abs(Gy @ y + Gc) < tolerance)
    columns = np.column_stack([2 * (Q @ y - g)] + [-Gy[i] for i in active])
    weights, _ = nnls(columns, ly)
    prices = np.zeros(Gy.shape[0])
    prices[active] = weights[1:]
    return float(weights[0]), prices


def DualBound(multiplier, prices, Q, g, k, ly, lc, Gy, Gc, budget):
    """Weak duality: sup_y [l'y + nu'(Gy y + Gc) - lambda (q(y) - budget)], valid for any
    lambda > 0 and nu >= 0, in closed form."""
    linear = ly + Gy.T @ prices + 2 * multiplier * g
    best = np.linalg.solve(2 * multiplier * Q, linear)
    value = (ly @ best + lc + prices @ (Gy @ best + Gc)
             - multiplier * (best @ Q @ best - 2 * g @ best + k - budget))
    return float(value)
