"""Strictly convex quadratic programs with linear inequality constraints."""
import numpy as np


def QuadraticProgram(metric, linear, rows, bounds, margin=0.0):
    """Minimize x'Mx/2 - c'x subject to rows @ x >= bounds; None if infeasible.

    A positive ``margin`` asks each constraint to hold by that much more, measured as a
    distance in the metric.
    """
    metric = np.asarray(metric, dtype=float)
    linear = np.asarray(linear, dtype=float)
    rows = np.atleast_2d(np.asarray(rows, dtype=float))
    bounds = np.asarray(bounds, dtype=float).ravel()

    factor = np.linalg.cholesky(metric)
    shift = np.linalg.solve(factor, linear)
    unconstrained = np.linalg.solve(factor.T, shift)
    if rows.size == 0:
        return unconstrained

    transformed = np.linalg.solve(factor, rows.T).T
    offset = bounds - rows @ unconstrained
    scale = np.linalg.norm(transformed, axis=1)
    keep = scale > 0.0
    if np.any(offset[~keep] > 0.0):
        return None
    distance = LeastDistance(transformed[keep] / scale[keep, None],
                             offset[keep] / scale[keep] + margin)
    if distance is None:
        return None
    return unconstrained + np.linalg.solve(factor.T, distance)


def LeastDistance(rows, bounds):
    """Minimize |z| subject to rows @ z >= bounds, for rows of unit norm; None if infeasible.

    Rows whose bound is far below zero cannot bind near the solution, so they are set aside
    while solving and checked afterwards.
    """
    rows = np.atleast_2d(np.asarray(rows, dtype=float))
    bounds = np.asarray(bounds, dtype=float).ravel()
    if bounds.size == 0 or np.all(bounds <= 0.0):
        return np.zeros(rows.shape[1])

    reach = 1e3 * (1.0 + bounds.max())
    while True:
        near = bounds >= -reach
        solution = _LeastDistanceOf(rows[near], bounds[near])
        if solution is None or np.all(rows[~near] @ solution >= bounds[~near]):
            return solution
        reach *= 1e3


def _LeastDistanceOf(rows, bounds):
    system = np.vstack([rows.T, bounds])
    goal = np.zeros(system.shape[0])
    goal[-1] = 1.0
    weights = NonNegativeLeastSquares(system, goal)
    residual = system @ weights - goal
    if abs(residual[-1]) <= 1e-12:
        return None
    solution = -residual[:-1] / residual[-1]
    slack = 1e-6 * (1.0 + np.linalg.norm(solution))
    return solution if np.all(rows @ solution >= bounds - slack) else None


def NonNegativeLeastSquares(matrix, vector, tolerance=None, max_iterations=None):
    """Lawson-Hanson active-set solution of min |Ax - b| subject to x >= 0."""
    matrix = np.asarray(matrix, dtype=float)
    vector = np.asarray(vector, dtype=float)
    columns = matrix.shape[1]
    if tolerance is None:
        tolerance = 10.0 * np.finfo(float).eps * np.linalg.norm(matrix, 1) * max(matrix.shape)
    if max_iterations is None:
        max_iterations = 3 * columns

    solution = np.zeros(columns)
    passive = np.zeros(columns, dtype=bool)
    gradient = matrix.T @ (vector - matrix @ solution)
    for _ in range(max_iterations):
        candidates = ~passive & (gradient > tolerance)
        if not candidates.any():
            break
        passive[np.argmax(np.where(candidates, gradient, -np.inf))] = True
        while True:
            trial = np.zeros(columns)
            trial[passive] = np.linalg.lstsq(matrix[:, passive], vector, rcond=None)[0]
            if np.all(trial[passive] > 0.0):
                solution = trial
                break
            blocking = passive & (trial <= 0.0)
            ratio = solution[blocking] / (solution[blocking] - trial[blocking])
            solution = solution + ratio.min() * (trial - solution)
            passive &= solution > tolerance
            solution[~passive] = 0.0
        gradient = matrix.T @ (vector - matrix @ solution)
    return solution
