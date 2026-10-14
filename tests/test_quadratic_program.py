"""The quadratic programs behind BoundaryWalk, against closed forms."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
from DirectedSensitivity.QuadraticProgram import (LeastDistance, NonNegativeLeastSquares,
                                                  QuadraticProgram)

METRIC = np.array([[2.0, 0.5], [0.5, 1.0]])
ROW = np.array([1.0, 2.0])


def halfspace_solution(bound, margin=0.0):
    """argmin x'Mx/2 subject to ROW @ x >= bound > 0, with the bound raised by the margin."""
    direction = np.linalg.solve(METRIC, ROW)
    scale = np.sqrt(ROW @ direction)
    return (bound / scale + margin) / scale * direction


def test_an_inactive_constraint_leaves_the_unconstrained_optimum():
    linear = np.array([1.0, -1.0])
    solution = QuadraticProgram(METRIC, linear, [ROW], [-100.0])
    np.testing.assert_allclose(solution, np.linalg.solve(METRIC, linear))


def test_one_binding_halfspace_matches_its_closed_form():
    solution = QuadraticProgram(METRIC, np.zeros(2), [ROW], [3.0])
    np.testing.assert_allclose(solution, halfspace_solution(3.0), atol=1e-12)
    assert ROW @ solution == pytest.approx(3.0)


def test_a_margin_moves_the_solution_that_far_inside_in_the_metric():
    solution = QuadraticProgram(METRIC, np.zeros(2), [ROW], [3.0], margin=0.5)
    np.testing.assert_allclose(solution, halfspace_solution(3.0, margin=0.5), atol=1e-12)


def test_an_infeasible_pair_returns_none():
    rows = np.array([[1.0, 0.0], [-1.0, 0.0]])
    assert QuadraticProgram(METRIC, np.zeros(2), rows, [1.0, 0.0]) is None
    assert LeastDistance(rows, [1.0, 0.0]) is None


def test_a_zero_row_with_a_positive_bound_is_infeasible():
    assert QuadraticProgram(METRIC, np.zeros(2), [[0.0, 0.0]], [1.0]) is None


def test_least_distance_meets_two_binding_rows_at_their_intersection():
    rows = np.array([[1.0, 0.0], [0.0, 1.0]])
    np.testing.assert_allclose(LeastDistance(rows, [1.0, 2.0]), [1.0, 2.0])


TILT = 2e-4
NEAR = np.array([[1.0, 0.0, 0.0], np.array([-1.0, TILT, 0.0]) / np.hypot(1.0, TILT)])


def test_a_far_row_set_aside_is_brought_back_when_violated():
    assert LeastDistance(NEAR, [1.0, 0.0])[2] == 0.0
    far = np.array([0.0, -1.0, 1.0]) / np.sqrt(2.0)
    solution = LeastDistance(np.vstack([NEAR, far]), [1.0, 0.0, -3000.0 / np.sqrt(2.0)])
    np.testing.assert_allclose(solution, [1.0, 1.0 / TILT, 1.0 / TILT - 3000.0], rtol=1e-5)


def test_a_far_row_that_makes_the_problem_infeasible_returns_none():
    far = np.array([0.0, -1.0, 0.0])
    assert LeastDistance(np.vstack([NEAR, far]), [1.0, 0.0, -3000.0]) is None


def test_non_negative_least_squares_matches_a_hand_solved_case():
    matrix = np.array([[1.0, 0.0], [1.0, 1.0]])
    vector = np.array([2.0, -1.0])
    solution = NonNegativeLeastSquares(matrix, vector)
    np.testing.assert_allclose(solution, [0.5, 0.0])
    gradient = matrix.T @ (vector - matrix @ solution)
    assert np.all(gradient <= 1e-12)
    assert np.all(np.abs(gradient * solution) <= 1e-12)


def test_non_negative_least_squares_returns_an_interior_least_squares_solution():
    matrix = np.array([[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    vector = np.array([1.0, 3.0, 2.0])
    np.testing.assert_allclose(NonNegativeLeastSquares(matrix, vector),
                               np.linalg.lstsq(matrix, vector, rcond=None)[0])
