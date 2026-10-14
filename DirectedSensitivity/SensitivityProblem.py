"""Definition of a sensitivity-bound problem."""
import numpy as np

from .EvaluationTrace import EvaluationTrace


class SensitivityProblem:
    """A criterion, a target, an estimate, and optional named parameter constraints.

    ``constraints`` is a dict of named callables, or a sequence of them, each nonnegative
    where the constraint holds. ``to_primitive`` maps the solver's coordinates to the
    point the criterion, target, and constraints read. ``feasibility_tolerance`` is how far
    the fit and each constraint may be violated and still count as feasible.

    ``criterion_gradient``/``target_gradient`` are optional analytic or automatic
    derivatives in the solver's coordinates; each new point they are run at counts as
    one evaluation of the criterion/target.
    """

    def __init__(self, criterion, target, estimate, constraints=None,
                 constraint_jacobian=None, to_primitive=None,
                 feasibility_tolerance=0.0, trace=None, criterion_gradient=None,
                 target_gradient=None):
        if not callable(criterion):
            raise TypeError("criterion must be callable")
        if not callable(target):
            raise TypeError("target must be callable")

        self.trace = trace if trace is not None else EvaluationTrace()
        self.to_primitive = to_primitive
        self.feasibility_tolerance = float(feasibility_tolerance)
        self.criterion = self._WatchCriterion(_OnPrimitive(criterion, to_primitive))
        self.target = self.trace.Watch('target', _OnPrimitive(target, to_primitive))
        self.criterion_gradient = self.trace.Watch('criterion', criterion_gradient)
        self.target_gradient = self.trace.Watch('target', target_gradient)
        self.estimate = estimate
        self.constraint_names, functions = _NameConstraints(constraints)
        self.constraints = [_OnPrimitive(function, to_primitive) for function in functions]
        self.constraint_jacobian = constraint_jacobian
        self.criterion_at_estimate = float(self.criterion(estimate))
        self.target_at_estimate = float(self.target(estimate))

        if self.ConstraintMargin(estimate) < -self.feasibility_tolerance:
            raise ValueError("estimate violates parameter constraints: "
                             + ", ".join(self.ViolatedConstraints(estimate)))

    def SignedDeviation(self, primitives):
        return float(self.target(primitives)) - self.target_at_estimate

    def Objective(self, primitives):
        return self.SignedDeviation(primitives) ** 2

    def Slack(self, primitives):
        """The criterion at the estimate less the criterion at ``primitives``."""
        return self.criterion_at_estimate - float(self.criterion(primitives))

    def ConstraintValues(self, primitives):
        if not self.constraints:
            return np.empty(0)
        return np.array([float(constraint(primitives))
                         for constraint in self.constraints])

    def ConstraintMargin(self, primitives):
        values = self.ConstraintValues(primitives)
        return float(values.min()) if values.size else float("inf")

    def ConstraintReport(self, primitives):
        return dict(zip(self.constraint_names, self.ConstraintValues(primitives)))

    def ViolatedConstraints(self, primitives):
        return [name for name, margin in self.ConstraintReport(primitives).items()
                if margin < -self.feasibility_tolerance]

    def IsFeasible(self, primitives):
        """Whether ``primitives`` fit at least as well as the estimate and meet every
        constraint, up to ``feasibility_tolerance``."""
        allowance = -self.feasibility_tolerance
        return (self.Slack(primitives) >= allowance
                and self.ConstraintMargin(primitives) >= allowance)

    def _WatchCriterion(self, criterion):
        return self.trace.Watch('criterion', criterion)


def _OnPrimitive(function, to_primitive):
    if function is None or to_primitive is None:
        return function
    return lambda point: function(to_primitive(point))


def _NameConstraints(constraints):
    if constraints is None:
        return (), ()
    if isinstance(constraints, dict):
        return tuple(constraints), tuple(constraints.values())
    constraints = tuple(constraints)
    names = tuple(f'constraint[{index}]' for index in range(len(constraints)))
    return names, constraints
