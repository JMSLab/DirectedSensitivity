"""Abstract base class for sensitivity-bound solvers."""
import abc


class Solver(abc.ABC):
    """Iterative search that keeps the best feasible certificate in each direction.

    ``positive_primitive``/``negative_primitive`` are the feasible points that push
    the target farthest above/below its estimate value; ``best_primitive`` is the more
    extreme of the two. The matching ``*_deviation`` attributes hold the signed target
    deviation each attains (``best_objective`` is its square).

    ``feasible_candidates``/``infeasible_candidates`` count the proposals the search
    kept and discarded; ``trace`` exposes the problem's evaluation counts.
    """

    def __init__(self, problem):
        self.problem = problem
        self.best_primitive = problem.estimate
        self.best_objective = problem.Objective(problem.estimate)
        self.positive_primitive = problem.estimate
        self.negative_primitive = problem.estimate
        self.positive_deviation = 0.0
        self.negative_deviation = 0.0
        self.iteration = 0
        self.feasible_candidates = 0
        self.infeasible_candidates = 0
        self.converged = False

    @property
    def trace(self):
        return self.problem.trace

    @abc.abstractmethod
    def Step(self):
        """Return the next candidate primitive, or None to stop."""

    def Solve(self, max_iterations=10_000, target_objective=None):
        """Step until ``Step`` returns None, ``max_iterations`` is reached, or
        ``best_objective`` reaches ``target_objective``; ``converged`` is set only in the
        first case, and only if some candidate was feasible."""
        while not self._ShouldStop(max_iterations, target_objective):
            self.iteration += 1
            candidate = self.Step()
            if candidate is None:
                self.converged = self.feasible_candidates > 0
                break
            if self.problem.IsFeasible(candidate):
                self.feasible_candidates += 1
                self._Record(candidate)
            else:
                self.infeasible_candidates += 1
        return self

    def _Record(self, candidate):
        deviation = self.problem.SignedDeviation(candidate)
        if deviation > self.positive_deviation:
            self.positive_deviation = deviation
            self.positive_primitive = candidate
        if deviation < self.negative_deviation:
            self.negative_deviation = deviation
            self.negative_primitive = candidate
        objective = deviation ** 2
        if objective > self.best_objective:
            self.best_objective = objective
            self.best_primitive = candidate

    def _ShouldStop(self, max_iterations, target_objective):
        reached_iterations = (max_iterations is not None
                              and self.iteration >= max_iterations)
        reached_objective = (target_objective is not None
                             and self.best_objective >= target_objective)
        return reached_iterations or reached_objective
