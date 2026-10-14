"""Shared base for two-sided, gradient-driven solvers over vector primitives."""
import abc

import numpy as np

from .Numerics import CentralDifference, CopyValue
from .Solver import Solver


class GradientSolver(Solver):
    """Gradient plumbing and the two-sided sign loop, walked alternately."""

    _PHASE_STATE = ()

    def __init__(self, problem, finite_difference_step=1e-6):
        super().__init__(problem)
        self._budget = problem.criterion_at_estimate
        self._estimate_vector = self._ToVector(problem.estimate)
        self._dimension = len(self._estimate_vector)
        self._finite_difference_step = float(finite_difference_step)

    def _Begin(self):
        self._signs = [1.0, -1.0]
        self._saved, self._exhausted = {}, {}
        for sign in self._signs:
            self._EnterPhase(sign)
            self._saved[sign] = self._CaptureState()
            self._exhausted[sign] = False
        self._turn = 0
        self._RestoreState(self._saved[self._signs[0]])

    def _EnterPhase(self, sign):
        self._sign = sign
        self._point = self._estimate_vector.copy()
        before = set(vars(self))
        self._ResetState()
        undeclared = (set(vars(self)) - before) - set(self._PHASE_STATE)
        if undeclared:
            raise TypeError(f"{type(self).__name__} must declare "
                            f"{sorted(undeclared)} in _PHASE_STATE")

    def StartPhase(self, sign):
        """Walk ``sign`` alone, leaving the caller free to move ``_point`` after."""
        self._EnterPhase(sign)
        self._signs = [sign]
        self._saved, self._exhausted = {}, {sign: False}
        self._turn = 0

    @abc.abstractmethod
    def _ResetState(self):
        """Clear per-sign state at the start of a phase."""

    @abc.abstractmethod
    def _Propose(self):
        """Advance one step; return the new point, or None when the sign is spent."""

    def Step(self):
        if len(self._signs) == 1:
            return self._Advance()
        for _ in range(len(self._signs)):
            sign = self._signs[self._turn % len(self._signs)]
            self._turn += 1
            if self._exhausted[sign]:
                continue
            self._RestoreState(self._saved[sign])
            candidate = self._Advance()
            self._saved[sign] = self._CaptureState()
            if candidate is not None:
                return candidate
            self._exhausted[sign] = True
        return None

    def _Advance(self):
        candidate = self._Propose()
        return None if candidate is None else self._Wrap(candidate)

    def _CaptureState(self):
        names = ('_point', '_sign') + tuple(self._PHASE_STATE)
        return {name: CopyValue(getattr(self, name, None)) for name in names}

    def _RestoreState(self, state):
        for name, value in state.items():
            setattr(self, name, value)

    def _IsFeasible(self, point):
        return self.problem.IsFeasible(self._Wrap(point))

    def _Criterion(self, point):
        return float(self.problem.criterion(self._Wrap(point)))

    def _Target(self, point):
        return float(self.problem.target(self._Wrap(point)))

    def _TargetGradient(self, point):
        return self._Gradient(self.problem.target_gradient, self._Target, point)

    def _CriterionGradient(self, point):
        return self._Gradient(self.problem.criterion_gradient, self._Criterion, point)

    def _Gradient(self, derivative, function, point):
        if derivative is not None:
            return self._ToVector(derivative(self._Wrap(point)))
        return self._NumericalGradient(function, point)

    def _ConstraintValues(self, point):
        return self.problem.ConstraintValues(self._Wrap(point))

    def _ConstraintJacobian(self, point):
        constraints = self.problem.constraints
        if not constraints:
            return np.empty((0, len(point)))
        if self.problem.constraint_jacobian is not None:
            jacobian = self.problem.constraint_jacobian(self._Wrap(point))
            return np.atleast_2d(np.asarray(jacobian, dtype=float))
        return np.array([
            self._NumericalGradient(
                lambda vector, constraint=constraint: float(constraint(self._Wrap(vector))),
                point)
            for constraint in constraints])

    def _ConstraintSystem(self, point, criterion_gradient):
        """The fit and every parameter constraint, each feasible iff its value is <= 0."""
        values = [self._Criterion(point) - self._budget]
        gradients = [criterion_gradient]
        constraint_values = self._ConstraintValues(point)
        if constraint_values.size:
            jacobian = self._ConstraintJacobian(point)
            for value, gradient in zip(constraint_values, jacobian):
                values.append(-float(value))
                gradients.append(-gradient)
        return np.array(values), np.array(gradients)

    def _SearchAlong(self, point, direction, sign, accept):
        """The first trial ``accept`` takes along ``direction``, and its improvement.

        Starts from the step that last succeeded, grows it after a first-try success, and
        shrinks it after each refusal. ``accept(trial, sign, current_objective)`` returns
        the accepted point and its objective, or None.
        """
        norm = np.linalg.norm(direction)
        if norm < self._convergence_tolerance:
            return None
        unit = direction / norm
        current_objective = sign * self._Target(point)

        step = self._step
        for backtrack in range(self._max_backtracks):
            accepted = accept(point + step * unit, sign, current_objective)
            if accepted is not None:
                self._step = step * self._growth_factor if backtrack == 0 else step
                return accepted[0], accepted[1] - current_objective
            step *= self._backtracking_factor

        self._step = step
        return None

    def _NumericalGradient(self, function, point):
        return np.array([CentralDifference(function, point, i, self._finite_difference_step)
                         for i in range(len(point))])

    def _ToVector(self, primitive):
        return np.atleast_1d(np.asarray(primitive, dtype=float)).ravel()

    def _Wrap(self, vector):
        estimate = self.problem.estimate
        if isinstance(estimate, np.ndarray):
            return np.asarray(vector, dtype=float).reshape(estimate.shape)
        if isinstance(estimate, (list, tuple)):
            return type(estimate)(float(component) for component in vector)
        if np.ndim(estimate) == 0:
            return float(vector[0])
        return np.asarray(vector, dtype=float)
