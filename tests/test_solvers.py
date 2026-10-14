"""Each solver recovers known bounds while every certificate it reports stays feasible."""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
from DirectedSensitivity import (BoundaryWalk, FeasibleSequentialQuadraticProgramming,
                                 GeneralizedReducedGradient, GMMSensitivityProblem,
                                 SensitivityProblem, Solver)

SOLVERS = [GeneralizedReducedGradient, FeasibleSequentialQuadraticProgramming,
           BoundaryWalk]
CONSTRAINED_SOLVERS = SOLVERS

CENTER = np.array([0.5, -0.2, 0.1])
TARGET = np.array([1.0, 2.0, -0.5])
ESTIMATE = np.array([1.0, 0.0, 0.0])
CORRELATED = np.array([[2.0, 0.3, 0.1], [0.3, 1.5, 0.2], [0.1, 0.2, 1.0]])

SAMPLE_MEAN, SAMPLE_SD = 1.5, 2.0
WEIGHTS = [0.25, 0.5, 0.75]


def ellipsoid(weight=CORRELATED, constraints=None):
    return GMMSensitivityProblem(
        moments=lambda x: np.asarray(x, float) - CENTER, weight=weight,
        target=lambda x: float(TARGET @ np.asarray(x, float)), estimate=ESTIMATE,
        constraints=constraints)


def ellipsoid_bounds(weight=CORRELATED):
    """Closed-form upper and lower deviations of the target, and the farther argmax."""
    inverse = np.linalg.inv(weight)
    budget = float((ESTIMATE - CENTER) @ weight @ (ESTIMATE - CENTER))
    scaled = float(TARGET @ inverse @ TARGET)
    center = float(TARGET @ (CENTER - ESTIMATE))
    reach = math.sqrt(budget * scaled)
    argmax = CENTER + np.sign(center) * math.sqrt(budget / scaled) * (inverse @ TARGET)
    return {'upper': center + reach, 'lower': center - reach,
            'objective': (abs(center) + reach) ** 2, 'argmax': argmax}


def exponential(weight_on_sd):
    """Gamma family nesting the exponential; the target is the standard deviation."""
    def moments(theta):
        scale, excess = theta
        return np.array([SAMPLE_MEAN - scale * (1.0 + excess),
                         SAMPLE_SD - scale * np.sqrt(max(1.0 + excess, 0.0))])

    estimate = (1.0 - weight_on_sd) * SAMPLE_MEAN + weight_on_sd * SAMPLE_SD
    return GMMSensitivityProblem(
        moments=moments, weight=np.diag([1.0 - weight_on_sd, weight_on_sd]),
        target=lambda theta: float(theta[0] * np.sqrt(1.0 + theta[1])),
        estimate=np.array([estimate, 0.0]))


def exponential_bounds(weight_on_sd):
    """Closed-form deviations and argmax for 0 < w < 1."""
    root = math.sqrt(1.0 - weight_on_sd)
    estimate = (1.0 - weight_on_sd) * SAMPLE_MEAN + weight_on_sd * SAMPLE_SD
    upper = (SAMPLE_SD - SAMPLE_MEAN) * root * (1.0 + root)
    lower = (SAMPLE_SD - SAMPLE_MEAN) * root * (root - 1.0)
    farthest = estimate + upper
    return {'upper': upper, 'lower': lower,
            'argmax': np.array([farthest ** 2 / SAMPLE_MEAN,
                                (SAMPLE_MEAN / farthest) ** 2 - 1.0])}


def numerical_gradient(function, point, spacing=1e-6):
    point = np.asarray(point, float)
    gradient = np.empty_like(point)
    for i in range(len(point)):
        step = spacing * max(1.0, abs(point[i]))
        forward, backward = point.copy(), point.copy()
        forward[i] += step
        backward[i] -= step
        gradient[i] = (function(forward) - function(backward)) / (2.0 * step)
    return gradient


def recorded_proposals(solver, **solve_options):
    proposed = []
    step = solver.Step

    def recording_step():
        candidate = step()
        if candidate is not None:
            proposed.append(candidate)
        return candidate

    solver.Step = recording_step
    solver.Solve(**solve_options)
    return proposed


# --- known answers ------------------------------------------------------------

@pytest.mark.parametrize('solver_class', SOLVERS)
def test_matches_the_closed_form_ellipsoid_bound(solver_class):
    problem, bounds = ellipsoid(), ellipsoid_bounds()
    solver = solver_class(problem).Solve()
    assert solver.converged is True
    assert solver.best_objective == pytest.approx(bounds['objective'], rel=1e-6)
    assert np.allclose(solver.best_primitive, bounds['argmax'], atol=1e-4)
    assert solver.positive_deviation == pytest.approx(bounds['upper'], rel=1e-2)
    assert solver.negative_deviation == pytest.approx(bounds['lower'], rel=1e-2)
    assert problem.IsFeasible(solver.positive_primitive)
    assert problem.IsFeasible(solver.negative_primitive)


@pytest.mark.parametrize('solver_class', [GeneralizedReducedGradient,
                                          FeasibleSequentialQuadraticProgramming])
def test_analytic_gradients_reach_the_ellipsoid_bound(solver_class):
    def criterion(x):
        gap = np.asarray(x, float) - CENTER
        return float(gap @ CORRELATED @ gap)

    problem = SensitivityProblem(
        criterion, lambda x: float(TARGET @ np.asarray(x, float)), ESTIMATE,
        criterion_gradient=lambda x: 2.0 * CORRELATED @ (np.asarray(x, float) - CENTER),
        target_gradient=lambda x: TARGET)
    solver = solver_class(problem).Solve()
    assert solver.best_objective == pytest.approx(ellipsoid_bounds()['objective'], rel=1e-6)


@pytest.mark.parametrize('solver_class', SOLVERS)
@pytest.mark.parametrize('w', WEIGHTS)
def test_matches_the_closed_form_exponential_bound(solver_class, w):
    problem, bounds = exponential(w), exponential_bounds(w)
    solver = solver_class(problem).Solve()
    assert solver.positive_deviation == pytest.approx(bounds['upper'], rel=1e-3)
    assert solver.negative_deviation == pytest.approx(bounds['lower'], rel=1e-2)
    np.testing.assert_allclose(solver.best_primitive, bounds['argmax'], atol=1e-3)
    assert problem.IsFeasible(solver.positive_primitive)
    assert problem.IsFeasible(solver.negative_primitive)


@pytest.mark.parametrize('solver_class', [GeneralizedReducedGradient,
                                          FeasibleSequentialQuadraticProgramming])
def test_returned_point_satisfies_the_kkt_conditions(solver_class):
    problem = SensitivityProblem(lambda p: p[0] ** 2 + p[1] ** 2,
                                 lambda p: p[0] ** 2 + p[1], estimate=[1.0, 0.0])
    point = np.asarray(solver_class(problem).Solve().best_primitive, float)
    assert problem.criterion(list(point)) == pytest.approx(
        problem.criterion_at_estimate, abs=1e-6)
    criterion_gradient = numerical_gradient(lambda x: problem.criterion(list(x)), point)
    target_gradient = numerical_gradient(lambda x: problem.target(list(x)), point)
    tangential = (target_gradient
                  - (target_gradient @ criterion_gradient)
                  / (criterion_gradient @ criterion_gradient) * criterion_gradient)
    assert np.linalg.norm(tangential) <= 1e-6 * np.linalg.norm(target_gradient)


def test_boundary_walk_satisfies_the_kkt_conditions():
    problem = GMMSensitivityProblem(
        moments=lambda t: np.array([t[0] - 1.0, t[0] * t[1]]),
        weight=np.eye(2), target=lambda t: float(t[1]),
        estimate=np.array([1.0, 0.5]))
    point = np.asarray(BoundaryWalk(problem).Solve().best_primitive, float)
    assert problem.criterion(point) == pytest.approx(
        problem.criterion_at_estimate, abs=1e-6)
    criterion_gradient = problem.CriterionGradient(point)
    target_gradient = numerical_gradient(problem.target, point)
    tangential = (target_gradient
                  - (target_gradient @ criterion_gradient)
                  / (criterion_gradient @ criterion_gradient) * criterion_gradient)
    assert np.linalg.norm(tangential) <= 1e-3 * np.linalg.norm(target_gradient)


def test_boundary_walk_matches_the_bound_under_weak_identification():
    angle = 0.7
    rotation = np.array([[math.cos(angle), -math.sin(angle), 0.0],
                         [math.sin(angle), math.cos(angle), 0.0],
                         [0.0, 0.0, 1.0]])
    weight = rotation @ np.diag([3.0, 1.0, 0.01]) @ rotation.T
    solver = BoundaryWalk(ellipsoid(weight)).Solve()
    assert solver.best_objective == pytest.approx(
        ellipsoid_bounds(weight)['objective'], rel=1e-6)
    assert solver.best_objective > 10.0


# --- feasibility along the way ------------------------------------------------

@pytest.mark.parametrize('solver_class', SOLVERS)
def test_every_proposed_point_on_the_path_is_feasible(solver_class):
    problem = ellipsoid()
    proposed = recorded_proposals(solver_class(problem))
    assert len(proposed) > 5
    assert all(problem.IsFeasible(point) for point in proposed)


def test_boundary_walk_follows_an_unidentified_direction_without_leaving_the_fit():
    problem = GMMSensitivityProblem(moments=lambda t: np.array([t[0] + t[1]]),
                                    weight=np.array([[1.0]]),
                                    target=lambda t: float(t[0] - t[1]),
                                    estimate=np.array([0.5, 0.5]))
    solver = BoundaryWalk(problem)
    proposed = recorded_proposals(solver, max_iterations=100)
    assert all(problem.IsFeasible(point) for point in proposed)
    assert solver.best_objective > 1e3


def test_infeasible_candidates_never_become_certificates():
    class SequenceSolver(Solver):
        def __init__(self, problem, candidates):
            super().__init__(problem)
            self._candidates = iter(candidates)

        def Step(self):
            return next(self._candidates, None)

    problem = SensitivityProblem(lambda x: (x - 1.0) ** 2, lambda x: x, estimate=2.0)
    solver = SequenceSolver(problem, candidates=[-1.0, 3.0, 0.5]).Solve()
    assert solver.negative_primitive == pytest.approx(0.5)
    assert solver.positive_primitive == pytest.approx(2.0)
    assert solver.infeasible_candidates == 2


# --- search behavior ----------------------------------------------------------

@pytest.mark.parametrize('solver_class', [GeneralizedReducedGradient,
                                          FeasibleSequentialQuadraticProgramming])
def test_recovers_when_the_first_direction_is_a_dead_end(solver_class):
    problem = SensitivityProblem(lambda x: (x - 1.0) ** 2, lambda x: x, estimate=2.0)
    solver = solver_class(problem).Solve()
    assert solver.best_objective == pytest.approx(4.0, rel=1e-6)
    assert solver.best_primitive == pytest.approx(0.0, abs=1e-5)


@pytest.mark.parametrize('solver_class', CONSTRAINED_SOLVERS)
def test_a_binding_constraint_caps_one_side_and_leaves_the_other(solver_class):
    bounds = ellipsoid_bounds()
    cap = 1.0 + 0.5 * bounds['upper']
    problem = ellipsoid(constraints={'cap': lambda x: cap - TARGET @ x})
    solver = solver_class(problem).Solve()
    assert problem.IsFeasible(solver.positive_primitive)
    assert problem.IsFeasible(solver.negative_primitive)
    assert solver.positive_deviation == pytest.approx(cap - 1.0, rel=1e-2)
    assert solver.negative_deviation == pytest.approx(bounds['lower'], rel=1e-2)


@pytest.mark.parametrize('solver_class', [GeneralizedReducedGradient,
                                          FeasibleSequentialQuadraticProgramming])
def test_slides_along_a_binding_constraint_to_the_closed_form_bound(solver_class):
    problem = SensitivityProblem(lambda x: x[0] ** 2 + x[1] ** 2, lambda x: x[0] + x[1],
                                 estimate=np.array([0.0, -1.0]),
                                 constraints={'cap': lambda x: 0.5 - x[0]})
    solver = solver_class(problem).Solve()
    assert solver.positive_deviation == pytest.approx(1.5 + math.sqrt(0.75), rel=1e-6)
    assert problem.IsFeasible(solver.positive_primitive)


@pytest.mark.parametrize('solver_class', CONSTRAINED_SOLVERS)
def test_a_short_run_carries_a_certificate_in_both_directions(solver_class):
    solver = solver_class(ellipsoid()).Solve(max_iterations=6)
    assert solver.positive_deviation > 0.0
    assert solver.negative_deviation < 0.0
