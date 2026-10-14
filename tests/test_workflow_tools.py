"""The checks and outputs the workflow relies on: limits, nesting, and the animation."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
from DirectedSensitivity import (DeformationPath, EvaluationBudgetExceeded,
                                 EvaluationTrace, GMMSensitivityProblem, NestsBaseline,
                                 SensitivityProblem, WriteAnimation)

ESTIMATE = np.array([2.0, 0.0])
LOWER = np.array([1.0, -2.0])
UPPER = np.array([1.0, 3.0])
GRID = np.linspace(0.0, 1.0, 7)


def criterion(x):
    return float((np.asarray(x)[0] - 1.0) ** 2)


def target(x):
    return float(np.asarray(x)[1])


def curve(x):
    point = np.asarray(x, dtype=float)
    return point[1] + point[0] * GRID


def make_problem(constraints=None, **kwargs):
    return SensitivityProblem(criterion, target, estimate=ESTIMATE,
                              constraints=constraints, **kwargs)


def route(steps=3):
    """A feasible path: the lower certificate, the baseline, then the upper one."""
    down = [LOWER + weight * (ESTIMATE - LOWER) for weight in np.linspace(0.0, 1.0, steps)]
    up = [ESTIMATE + weight * (UPPER - ESTIMATE)
          for weight in np.linspace(0.0, 1.0, steps)[1:]]
    return down + up


def test_an_evaluation_limit_stops_the_criterion_once_it_is_passed():
    problem = make_problem(trace=EvaluationTrace(limits={'criterion': 3}))
    problem.criterion(np.array([1.5, 0.0]))
    problem.criterion(np.array([1.6, 0.0]))
    with pytest.raises(EvaluationBudgetExceeded, match='criterion'):
        problem.criterion(np.array([1.7, 0.0]))


def make_gmm_problem(trace):
    return GMMSensitivityProblem(moments=lambda x: np.asarray(x, dtype=float) - 1.0,
                                 weight=np.eye(2), target=target, estimate=ESTIMATE,
                                 trace=trace)


def test_each_new_point_of_the_moments_counts_once_as_a_criterion_evaluation():
    trace = EvaluationTrace()
    problem = make_gmm_problem(trace)
    problem.LocalStructure(ESTIMATE)
    assert trace.Count('criterion') == 1 + 2 * ESTIMATE.size


def test_a_point_already_evaluated_is_not_evaluated_or_counted_again():
    calls = []
    problem = SensitivityProblem(lambda x: calls.append(1) or criterion(x), target,
                                 estimate=ESTIMATE)
    point = np.array([1.5, 0.0])
    problem.criterion(point)
    problem.IsFeasible(point)
    problem.IsFeasible(ESTIMATE)
    assert len(calls) == problem.trace.Count('criterion') == 2


def test_an_analytic_gradient_counts_once_per_new_point_as_its_function():
    problem = make_problem(criterion_gradient=lambda x: np.array([2.0 * (x[0] - 1.0), 0.0]),
                           target_gradient=lambda x: np.array([0.0, 1.0]))
    point = np.array([1.5, 0.0])
    for _ in range(2):
        problem.criterion_gradient(point)
        problem.target_gradient(point)
    assert problem.trace.Snapshot() == {'criterion': 2, 'target': 2}


def test_an_analytic_jacobian_counts_once_per_new_point_as_the_criterion():
    trace = EvaluationTrace()
    problem = GMMSensitivityProblem(moments=lambda x: np.asarray(x, dtype=float) - 1.0,
                                    weight=np.eye(2), target=target, estimate=ESTIMATE,
                                    jacobian=lambda x: np.eye(2), trace=trace)
    problem.LocalStructure(ESTIMATE)
    problem.Jacobian(ESTIMATE)
    assert trace.Count('criterion') == 2


def test_a_criterion_limit_caps_the_moments_behind_a_numerical_jacobian():
    problem = make_gmm_problem(EvaluationTrace(limits={'criterion': 3}))
    with pytest.raises(EvaluationBudgetExceeded, match='criterion'):
        problem.Jacobian(ESTIMATE)


def test_nesting_accepts_the_estimate_and_rejects_a_baseline_that_is_not_it():
    assert NestsBaseline(make_problem(), ESTIMATE)['nests'] is True
    displaced = NestsBaseline(make_problem(), ESTIMATE + np.array([0.5, 0.0]))
    assert displaced['nests'] is False
    assert displaced['criterion_gap'] > 0.0
    capped = make_problem(constraints={'cap': lambda x: 2.5 - x[0]})
    breaking = NestsBaseline(capped, ESTIMATE + np.array([1.0, 0.0]))
    assert breaking['nests'] is False
    assert breaking['violated'] == ['cap']


def test_the_path_must_pass_through_the_baseline_and_stay_feasible():
    with pytest.raises(ValueError, match='must pass through the baseline'):
        DeformationPath(make_problem(), [LOWER, UPPER])
    with pytest.raises(ValueError, match='not feasible'):
        DeformationPath(make_problem(), [LOWER, ESTIMATE, np.array([4.0, 3.0])])


def test_the_animation_and_figure_show_the_baseline_and_both_certificates(tmp_path):
    from PIL import Image

    path = DeformationPath(make_problem(), route(), curve=curve)
    seen = {}
    destination = WriteAnimation(
        tmp_path / 'animation.gif', path, grid=GRID, figure=tmp_path / 'certificates.png',
        annotate=lambda axes: seen.update(
            labels=[line.get_label() for line in axes.get_lines()]))
    with Image.open(destination) as loop:
        assert loop.n_frames == 2 * len(path['frames']) - 2
    assert (tmp_path / 'certificates.png').stat().st_size > 0
    assert [label.split(':')[0] for label in seen['labels']] == [
        'baseline', 'lower certificate', 'upper certificate']
