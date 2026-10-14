"""Checks the workflow run on the exponential example against its closed form."""
import json
import math
from pathlib import Path

import pytest

RUN = Path(__file__).parents[1] / 'workflow' / 'examples' / 'exponential' / 'output'
SAMPLE_MEAN, SAMPLE_SD, WEIGHT = 1.5, 2.0, 0.5
ESTIMATE = 1.75
LIMITS = {'criterion': 5000, 'target': 5000}
TOLERANCE = 1e-10

REACH = math.sqrt(1.0 - WEIGHT) * abs(SAMPLE_MEAN - SAMPLE_SD)
EXACT = {'budget': WEIGHT * (1.0 - WEIGHT) * (SAMPLE_MEAN - SAMPLE_SD) ** 2,
         'baseline': (1.0 - WEIGHT) * SAMPLE_MEAN + WEIGHT * SAMPLE_SD,
         'lower': SAMPLE_SD - REACH,
         'upper': SAMPLE_SD + REACH}

pytestmark = pytest.mark.skipif(
    not RUN.is_dir(), reason='the workflow has not been run on the exponential example')


@pytest.fixture(scope='module')
def certificates():
    return json.loads((RUN / 'certificates.json').read_text())


def Object(mean, standard_deviation):
    """Section 3, at a primitive with these two moments."""
    return standard_deviation


def Criterion(mean, standard_deviation):
    """Section 4, at a primitive with these two moments."""
    return ((1.0 - WEIGHT) * (SAMPLE_MEAN - mean) ** 2
            + WEIGHT * (SAMPLE_SD - standard_deviation) ** 2)
def test_the_run_stayed_within_the_declared_limits(certificates):
    assert certificates['limits'] == LIMITS
    for name, limit in LIMITS.items():
        assert certificates['spend'][name] <= limit


def test_the_baseline_is_the_working_model_at_the_estimate(certificates):
    baseline = certificates['baseline']
    assert baseline['primitive']['mean'] == pytest.approx(ESTIMATE, rel=1e-9)
    assert baseline['primitive']['sd'] == pytest.approx(ESTIMATE, rel=1e-9)
    assert baseline['criterion'] == pytest.approx(EXACT['budget'], rel=1e-9)
    assert baseline['target'] == pytest.approx(EXACT['baseline'], rel=1e-9)


def test_each_certificate_is_feasible_and_reports_its_own_numbers(certificates):
    baseline = certificates['baseline']
    for direction in ('lower', 'upper'):
        certificate = certificates[direction]
        moments = (certificate['primitive']['mean'], certificate['primitive']['sd'])
        assert Object(*moments) == pytest.approx(certificate['target'], abs=1e-9)
        assert Criterion(*moments) == pytest.approx(certificate['criterion'], abs=1e-9)
        assert baseline['criterion'] - certificate['criterion'] == pytest.approx(
            certificate['slack'], abs=1e-9)
        assert certificate['slack'] >= -TOLERANCE
        assert certificate['target'] - baseline['target'] == pytest.approx(
            certificate['deviation'], abs=1e-9)


def test_the_run_found_the_exact_extremes(certificates):
    for direction in ('lower', 'upper'):
        assert certificates['analytic'][direction] == pytest.approx(EXACT[direction],
                                                                    rel=1e-6)
        assert certificates[direction]['target'] == pytest.approx(EXACT[direction],
                                                                  rel=1e-3)
