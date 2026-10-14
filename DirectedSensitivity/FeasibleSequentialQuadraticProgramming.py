"""Feasible sequential quadratic programming solver."""
import numpy as np

from .GradientSolver import GradientSolver


class FeasibleSequentialQuadraticProgramming(GradientSolver):
    """Quasi-Newton QP steps tilted to stay inside the feasible set.

    The feasible set is the fit constraint together with any parameter
    constraints carried by the problem. Each step projects the quasi-Newton
    direction onto the constraints predicted to bind and tilts inward from them.
    """

    _PHASE_STATE = ('_inverse_hessian',)

    def __init__(self, problem, hessian_scale=1.0, tilt=0.1, armijo=1e-4,
                 backtracking_factor=0.5, max_backtracks=60,
                 convergence_tolerance=1e-8, curvature_tolerance=1e-8,
                 finite_difference_step=1e-6):
        super().__init__(problem, finite_difference_step)
        self._hessian_scale = float(hessian_scale)
        self._tilt = float(tilt)
        self._armijo = float(armijo)
        self._backtracking_factor = float(backtracking_factor)
        self._max_backtracks = int(max_backtracks)
        self._convergence_tolerance = float(convergence_tolerance)
        self._curvature_tolerance = float(curvature_tolerance)
        self._Begin()

    def _ResetState(self):
        self._inverse_hessian = self._hessian_scale * np.eye(self._dimension)

    def _Propose(self):
        target_gradient, criterion_gradient = self._Gradients(self._point)
        values, gradients = self._ConstraintSystem(self._point, criterion_gradient)
        step = self._DescentStep(self._point, target_gradient, values, gradients,
                                 self._sign)
        if step is None:
            return None
        new_point, multipliers = step
        self._UpdateInverseHessian(self._point, target_gradient, gradients,
                                   new_point, multipliers, self._sign)
        self._point = new_point
        return new_point

    def _DescentStep(self, point, target_gradient, values, gradients, sign):
        objective_gradient = -sign * target_gradient
        direction, multipliers = self._QuadraticProgramDirection(
            objective_gradient, values, gradients)

        slope = sign * float(target_gradient @ direction)
        if slope < self._convergence_tolerance:
            return None

        new_point = self._LineSearch(point, direction, slope, sign)
        if new_point is None:
            return None
        return new_point, multipliers

    def _QuadraticProgramDirection(self, objective_gradient, values, gradients):
        multipliers = np.zeros(len(values))
        unconstrained = -(self._inverse_hessian @ objective_gradient)

        linearized = values + gradients @ unconstrained
        active = np.where(linearized > 0.0)[0]
        if active.size == 0:
            return unconstrained, multipliers

        active_gradients = gradients[active]
        projected = self._inverse_hessian @ active_gradients.T
        gram = active_gradients @ projected
        if np.linalg.eigvalsh(gram).min() < self._convergence_tolerance:
            return unconstrained, multipliers

        active_multipliers = np.linalg.solve(gram, linearized[active])
        direction = unconstrained - projected @ active_multipliers
        multipliers[active] = active_multipliers

        inward = projected @ np.linalg.solve(gram, np.ones(active.size))
        return self._TiltInward(direction, objective_gradient, inward), multipliers

    def _TiltInward(self, direction, objective_gradient, inward):
        inward = inward / max(1.0, float(np.linalg.norm(inward)))
        magnitude = self._tilt * np.linalg.norm(direction)

        descent = float(objective_gradient @ direction)
        erosion = float(objective_gradient @ inward)
        if erosion > 0.0:
            magnitude = min(magnitude, 0.5 * max(-descent, 0.0) / erosion)
        return direction - magnitude * inward

    def _LineSearch(self, point, direction, slope, sign):
        current_objective = sign * self._Target(point)
        step = 1.0
        for _ in range(self._max_backtracks):
            trial = point + step * direction
            if self._IsFeasible(trial):
                improvement = sign * self._Target(trial) - current_objective
                if improvement >= self._armijo * step * slope:
                    return trial
            step *= self._backtracking_factor
        return None

    def _UpdateInverseHessian(self, point, target_gradient, gradients, new_point,
                              multipliers, sign):
        new_target_gradient, new_criterion_gradient = self._Gradients(new_point)
        _, new_gradients = self._ConstraintSystem(new_point, new_criterion_gradient)
        old_lagrangian = -sign * target_gradient + multipliers @ gradients
        new_lagrangian = -sign * new_target_gradient + multipliers @ new_gradients

        step = new_point - point
        gradient_change = new_lagrangian - old_lagrangian
        curvature = float(step @ gradient_change)
        if curvature <= self._curvature_tolerance:
            return

        rho = 1.0 / curvature
        identity = np.eye(self._dimension)
        left = identity - rho * np.outer(step, gradient_change)
        right = identity - rho * np.outer(gradient_change, step)
        self._inverse_hessian = (left @ self._inverse_hessian @ right
                                 + rho * np.outer(step, step))

    def _Gradients(self, point):
        return self._TargetGradient(point), self._CriterionGradient(point)
