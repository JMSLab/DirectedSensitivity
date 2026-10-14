"""Generalized reduced gradient solver."""
import numpy as np

from .GradientSolver import GradientSolver
from .Numerics import CentralDifference


class GeneralizedReducedGradient(GradientSolver):
    """Reduced-gradient walk along the surface of the fit and parameter constraints.

    Each constraint carries a slack and solves for one basic variable; the walk ascends
    in the others, never pushing the slack of a binding constraint below zero. After each
    step, Newton steps in the basic coordinates restore any violated constraint.
    """

    _PHASE_STATE = ('_step', '_restore_index')

    def __init__(self, problem, step_size=1.0, growth_factor=2.0,
                 backtracking_factor=0.5, max_backtracks=50,
                 max_restorations=25, feasibility_tolerance=1e-9,
                 convergence_tolerance=1e-8, finite_difference_step=1e-6):
        super().__init__(problem, finite_difference_step)
        self._initial_step = float(step_size)
        self._growth_factor = float(growth_factor)
        self._backtracking_factor = float(backtracking_factor)
        self._max_backtracks = int(max_backtracks)
        self._max_restorations = int(max_restorations)
        self._feasibility_tolerance = float(feasibility_tolerance)
        self._convergence_tolerance = float(convergence_tolerance)
        self._Begin()

    def _ResetState(self):
        self._step = self._initial_step

    def _Propose(self):
        direction = self._AscentDirection(self._point, self._sign)
        if direction is None:
            return None
        accepted = self._SearchAlong(self._point, direction, self._sign, self._Accept)
        if accepted is None:
            return None
        self._point = accepted[0]
        return self._point

    def _AscentDirection(self, point, sign):
        target_gradient = self._TargetGradient(point)
        values, gradients = self._ConstraintSystem(point, self._CriterionGradient(point))
        rows, dimension = gradients.shape

        system = np.hstack([gradients, np.eye(rows)])
        objective_gradient = np.concatenate([sign * target_gradient, np.zeros(rows)])
        at_bound = np.concatenate([np.zeros(dimension, dtype=bool),
                                   np.maximum(-values, 0.0) <= self._feasibility_tolerance])
        basic = self._Basis(system, at_bound)
        if basic is None:
            return None

        basis = system[:, basic]
        prices = _Solve(basis.T, objective_gradient[basic])
        reduced = objective_gradient - system.T @ prices
        reduced[basic] = 0.0
        reduced[at_bound] = np.maximum(reduced[at_bound], 0.0)
        if np.linalg.norm(reduced) < self._convergence_tolerance:
            return None

        direction = reduced.copy()
        direction[basic] = -_Solve(basis, system @ reduced)
        self._restore_index = _RestoreIndex(gradients, basic)
        return direction[:dimension]

    def _Basis(self, system, at_bound):
        """One basic column per row, by largest pivot, avoiding slacks held at their bound.

        A binding parameter constraint left without a pivot keeps its own slack basic; a
        fit without one has nowhere to walk.
        """
        remaining = system.copy()
        dimension = system.shape[1] - system.shape[0]
        eligible = ~at_bound
        basic = []
        for row in range(len(remaining)):
            pivots = np.where(eligible, np.abs(remaining[row]), -1.0)
            column = int(np.argmax(pivots))
            if pivots[column] < self._convergence_tolerance:
                if row == 0:
                    return None
                column = dimension + row
            basic.append(column)
            eligible[column] = False
            remaining[row + 1:] -= np.outer(remaining[row + 1:, column] / remaining[row, column],
                                            remaining[row])
        return basic

    def _Accept(self, trial, sign, current_objective):
        restored = self._Restore(trial)
        if restored is None or not self._IsFeasible(restored):
            return None
        objective = sign * self._Target(restored)
        return (restored, objective) if objective > current_objective else None

    def _Restore(self, point):
        """Newton steps on the restoring coordinates of the violated constraints."""
        point = point.copy()
        for _ in range(self._max_restorations):
            excess = self._Excess(point)
            if not np.all(np.isfinite(excess)):
                return None
            violated = np.flatnonzero(excess > 0.0)
            if violated.size == 0:
                return point
            columns = self._restore_index[violated]
            slopes = np.array([[self._Partial(point, row, column) for column in columns]
                               for row in violated])
            if np.linalg.svd(slopes, compute_uv=False).min() < self._convergence_tolerance:
                return None
            point[columns] -= _Solve(slopes, excess[violated])
        return point if np.all(self._Excess(point) <= 0.0) else None

    def _Excess(self, point):
        return np.concatenate([[self._Criterion(point) - self._budget],
                               -self._ConstraintValues(point)])

    def _Partial(self, point, row, column):
        """The partial of constraint ``row`` of the constraint system in ``column``."""
        if row == 0:
            if self.problem.criterion_gradient is not None:
                return self._CriterionGradient(point)[column]
            return CentralDifference(self._Criterion, point, column,
                                     self._finite_difference_step)
        if self.problem.constraint_jacobian is not None:
            return -self._ConstraintJacobian(point)[row - 1, column]
        constraint = self.problem.constraints[row - 1]
        return -CentralDifference(lambda vector: float(constraint(self._Wrap(vector))),
                                  point, column, self._finite_difference_step)


def _Solve(matrix, vector):
    """``np.linalg.solve``, dividing exactly when there is a single row."""
    if matrix.shape == (1, 1):
        return vector / matrix[0, 0]
    return np.linalg.solve(matrix, vector)


def _RestoreIndex(gradients, basic):
    """For each row, the coordinate that restores it: its basic one, else its steepest free one."""
    dimension = gradients.shape[1]
    taken = {column for column in basic if column < dimension}
    indices = []
    for row, column in enumerate(basic):
        if column >= dimension:
            free = [i for i in range(dimension) if i not in taken] or list(range(dimension))
            column = max(free, key=lambda i: abs(gradients[row, i]))
            taken.add(column)
        indices.append(column)
    return np.array(indices)
