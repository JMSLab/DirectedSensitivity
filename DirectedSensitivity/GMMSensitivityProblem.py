"""GMM sensitivity problem exposing Gauss-Newton structure."""
import numpy as np

from .EvaluationTrace import EvaluationTrace
from .Numerics import CentralDifference
from .SensitivityProblem import SensitivityProblem


class GMMSensitivityProblem(SensitivityProblem):
    """Criterion g'Wg with its moments, weight, Jacobian, and Gauss-Newton matrix."""

    def __init__(self, moments, weight, target, estimate, jacobian=None,
                 finite_difference_step=1e-6, constraints=None,
                 constraint_jacobian=None, feasibility_tolerance=0.0, trace=None,
                 target_gradient=None):
        trace = trace if trace is not None else EvaluationTrace()
        self.moments = trace.Watch('criterion', moments)
        self.weight = np.asarray(weight, dtype=float)
        self._jacobian = trace.Watch('criterion', jacobian)
        self._finite_difference_step = float(finite_difference_step)
        estimate = np.atleast_1d(np.asarray(estimate, dtype=float)).ravel()
        super().__init__(self._GMMCriterion, target, estimate,
                         constraints=constraints,
                         constraint_jacobian=constraint_jacobian,
                         feasibility_tolerance=feasibility_tolerance, trace=trace,
                         target_gradient=target_gradient)
        if jacobian is not None:
            self.criterion_gradient = self.CriterionGradient

    def _WatchCriterion(self, criterion):
        """Counted through the moments: each evaluation of them is one of the criterion."""
        return criterion

    def _GMMCriterion(self, theta):
        moments = self.Moments(theta)
        return float(moments @ self.weight @ moments)

    def Moments(self, theta):
        return np.atleast_1d(np.asarray(self.moments(theta), dtype=float)).ravel()

    def Jacobian(self, theta):
        if self._jacobian is not None:
            return np.asarray(self._jacobian(theta), dtype=float)
        return self._NumericalJacobian(theta)

    def _NumericalJacobian(self, theta):
        point = np.asarray(theta, dtype=float).ravel()
        return np.column_stack([
            CentralDifference(self.Moments, point, i, self._finite_difference_step)
            for i in range(len(point))])

    def CriterionGradient(self, theta):
        return self.LocalStructure(theta)[1]

    def LocalStructure(self, theta):
        moments = self.Moments(theta)
        jacobian = self.Jacobian(theta)
        weighted = self.weight @ jacobian
        gauss_newton = jacobian.T @ weighted
        gradient = 2.0 * jacobian.T @ (self.weight @ moments)
        return gauss_newton, gradient
