"""Predictor-corrector solver that ascends within {criterion <= budget}."""
import numpy as np

from .GradientSolver import GradientSolver
from .QuadraticProgram import QuadraticProgram


class BoundaryWalk(GradientSolver):
    """Walk along {criterion <= budget} in the (M + lambda D) metric."""

    _PHASE_STATE = ('_step', '_settled', 'stop_reason')
    CONVERGED = ('stationary', 'value')

    def __init__(self, problem, damping=1e-2, scaling_floor=1e-6, step_size=1.0,
                 growth_factor=2.0, backtracking_factor=0.5, max_backtracks=50,
                 max_corrections=25, convergence_tolerance=1e-8, value_tolerance=None,
                 correction_margin=1e-8, finite_difference_step=1e-6):
        super().__init__(problem, finite_difference_step)
        self._damping = float(damping)
        self._scaling_floor = float(scaling_floor)
        self._initial_step = float(step_size)
        self._growth_factor = float(growth_factor)
        self._backtracking_factor = float(backtracking_factor)
        self._max_backtracks = int(max_backtracks)
        self._max_corrections = int(max_corrections)
        self._convergence_tolerance = float(convergence_tolerance)
        self._value_tolerance = None if value_tolerance is None else float(value_tolerance)
        self._correction_margin = float(correction_margin)
        self._Begin()

    def _ResetState(self):
        self._step = self._initial_step
        self._settled = False
        self.stop_reason = None

    def _Propose(self):
        if self._settled:
            return self._Stop('value')
        predicted = self._Predictor(self._point, self._sign)
        if predicted is None:
            return self._Stop('no direction')
        direction, ascent = predicted
        if ascent < self._convergence_tolerance:
            return self._Stop('stationary')
        accepted = self._SearchAlong(self._point, direction, self._sign, self._Accept)
        if accepted is None:
            return self._Stop('no improving step')
        self._point, improvement = accepted
        self._settled = (self._value_tolerance is not None
                         and improvement < self._value_tolerance)
        return self._point

    def _Stop(self, reason):
        self.stop_reason = reason
        return None

    def _Predictor(self, point, sign):
        """Best ascent step in the metric keeping every linearized constraint, and its ascent."""
        grad_target = sign * self._TargetGradient(point)
        system = self._LinearizedConstraints(point)
        if system is None or not np.all(np.isfinite(grad_target)):
            return None
        direction = QuadraticProgram(system[0], grad_target, *system[1:])
        return None if direction is None else (direction, float(direction @ grad_target))

    def _Accept(self, trial, sign, current_objective):
        """The corrected trial and its objective if feasible and better; None otherwise or if
        the problem cannot evaluate it."""
        try:
            corrected = self._Correct(trial)
            if corrected is None or not self._IsFeasible(corrected):
                return None
            objective = sign * self._Target(corrected)
        except (ArithmeticError, ValueError, np.linalg.LinAlgError):
            return None
        return (corrected, objective) if objective > current_objective else None

    def _Correct(self, point):
        """Shortest moves in the metric onto the linearized constraints, until feasible."""
        point = point.copy()
        for _ in range(self._max_corrections):
            if self._IsFeasible(point):
                return point
            system = self._LinearizedConstraints(point)
            if system is None:
                return None
            move = QuadraticProgram(system[0], np.zeros(len(point)), *system[1:],
                                    margin=self._correction_margin)
            if move is None:
                return None
            point = point + move
        return point if self._IsFeasible(point) else None

    def _LinearizedConstraints(self, point):
        """The metric and rows @ move >= bounds for the fit and every parameter constraint."""
        gauss_newton, grad_criterion = self.problem.LocalStructure(self._Wrap(point))
        metric = gauss_newton + self._damping * np.diag(self._Scaling(gauss_newton))
        rows = [-grad_criterion[None, :]]
        bounds = [[self._Criterion(point) - self._budget]]
        values = self._ConstraintValues(point)
        if values.size:
            rows.append(self._ConstraintJacobian(point))
            bounds.append(-values)
        system = (metric, np.vstack(rows), np.concatenate(bounds))
        if not all(np.all(np.isfinite(part)) for part in system):
            return None
        return system

    def _Scaling(self, gauss_newton):
        diagonal = np.diag(gauss_newton).astype(float).copy()
        largest = diagonal.max()
        floor = self._scaling_floor * largest if largest > 0.0 else self._scaling_floor
        return np.maximum(diagonal, floor)
