import io
import os
import sys
import zipfile

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.special import expit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

DATA_PATH   = config.PATHS['data']
DATA_MEMBER = config.CONFIG['members']['data']
X_INTERVAL  = (-25.0, 25.0)
LIMITS      = {}

DELTA                = 1.0
PRESENT_BIAS_CEILING = 1.0
SEASON_DAYS          = 65

GRID_POINTS  = 8001      # spacing 0.00625
PROBE_POINTS = 81

TAIL_TOLERANCE      = 1e-8
TAIL_MEAN_TOLERANCE = 1e-6


class OutsideTheClass(ValueError):
    pass


class NetBenefitDistribution:

    def Cdf(self, x):
        raise NotImplementedError

    def Density(self, x):
        raise NotImplementedError

    def PartialMean(self, lower):
        """E[X 1{X > lower}]."""
        lower = max(float(lower), X_INTERVAL[0])
        value, _ = quad(lambda x: x * float(self.Density(x)), lower, X_INTERVAL[1],
                        epsabs=1e-12, limit=200)
        return value

    def Mean(self):
        return self.PartialMean(-np.inf)


class Logistic(NetBenefitDistribution):
    """The working model. Location is minus the cost of filing."""

    def __init__(self, location, scale=1.0):
        self.location, self.scale = float(location), float(scale)
        if not self.scale > 0.0:
            raise OutsideTheClass('the logistic scale must be positive')

    def Cdf(self, x):
        return expit((np.asarray(x, dtype=float) - self.location) / self.scale)

    def Density(self, x):
        z = expit((np.asarray(x, dtype=float) - self.location) / self.scale)
        return z * (1.0 - z) / self.scale

    def PartialMean(self, lower):
        if lower == -np.inf:
            return self.location
        if lower == np.inf:
            return 0.0
        z = (float(lower) - self.location) / self.scale
        tail = z * expit(-z) + np.logaddexp(0.0, -z)
        return float(self.location * expit(-z) + self.scale * tail)


class Primitive:

    def __init__(self, present_bias, distribution):
        self.present_bias = float(present_bias)
        self.distribution = distribution

    def Check(self):
        if not self.present_bias > 0.0:
            raise OutsideTheClass('present bias must be positive: the naive plan divides '
                                  'the continuation value by it')
        return self


def ToPrimitive(vector):
    present_bias, location = np.asarray(vector, dtype=float).ravel()
    return Primitive(present_bias, Logistic(location))


def LoadData(path=DATA_PATH, member=DATA_MEMBER):
    with zipfile.ZipFile(path) as archive:
        frame = pd.read_excel(io.BytesIO(archive.read(member))).iloc[:SEASON_DAYS]
    return {'filed': frame['filed'].values.astype(float),
            'at_risk': frame['frequency'].values.astype(float)}


DATA = LoadData()
FILED, AT_RISK = DATA['filed'], DATA['at_risk']


def Main():
    print(f'working model      beta={WORKING_MODEL.present_bias!r}, '
          f'location={WORKING_MODEL.distribution.location!r}')
    print(f'criterion budget   {Criterion(WORKING_MODEL)!r}')
    print(f'target             {Target(WORKING_MODEL)!r}')
    threshold = Thresholds(WORKING_MODEL)
    print(f'thresholds         [{threshold.min()!r}, {threshold.max()!r}]')
    probability = FilingProbabilities(WORKING_MODEL)
    print(f'hazard day 0, 64   {probability[0]!r}, {probability[-1]!r}')
    for name, constraint in ConstraintSet().items():
        print(f'  {name:<28} {constraint(WORKING_MODEL)!r}')


def Thresholds(primitive):
    """The daily filing thresholds, day 0 first; the day's draw must exceed one to file."""
    primitive = primitive.Check()
    beta, distribution = primitive.present_bias, primitive.distribution
    threshold = np.empty(SEASON_DAYS)
    threshold[-1] = beta * DELTA * distribution.PartialMean(-np.inf)
    for day in range(SEASON_DAYS - 2, -1, -1):
        planned = threshold[day + 1] / beta
        threshold[day] = (beta * DELTA * distribution.PartialMean(planned)
                          + float(distribution.Cdf(planned)) * DELTA * threshold[day + 1])
    return threshold


def FilingProbabilities(primitive):
    return 1.0 - np.asarray(primitive.distribution.Cdf(Thresholds(primitive)), dtype=float)


def Criterion(primitive):
    probability = FilingProbabilities(primitive)
    if not np.all(np.isfinite(probability)) or np.any(probability <= 0.0) \
            or np.any(probability >= 1.0):
        raise OutsideTheClass('the implied filing probabilities leave (0, 1)')
    return -float(FILED @ np.log(probability)
                  + (AT_RISK - FILED) @ np.log1p(-probability))


def Target(primitive):
    return primitive.Check().present_bias


def FitMaximumLikelihood(start=(0.8, -3.0)):
    result = minimize(
        lambda theta: Criterion(ToPrimitive(theta)), np.asarray(start, dtype=float),
        method='Nelder-Mead',
        options={'xatol': 1e-10, 'fatol': 1e-12, 'maxiter': 20000, 'maxfev': 20000})
    if not result.success:
        raise RuntimeError(f'maximum likelihood did not converge: {result.message}')
    return result.x


ESTIMATE = FitMaximumLikelihood()
WORKING_MODEL = ToPrimitive(ESTIMATE)


def ConstraintSet(grid_points=GRID_POINTS):
    grid = np.linspace(*X_INTERVAL, grid_points)
    probe = np.linspace(*X_INTERVAL, PROBE_POINTS)

    def present_bias_at_most_one(primitive):
        return float(PRESENT_BIAS_CEILING - primitive.present_bias)

    def non_decreasing(primitive):
        return float(np.diff(primitive.distribution.Cdf(grid)).min())

    def bounded_in_unit_interval(primitive):
        value = np.asarray(primitive.distribution.Cdf(grid), dtype=float)
        return float(min(value.min(), (1.0 - value).min()))

    def tails_reach_zero_and_one(primitive):
        distribution = primitive.distribution
        low = float(np.asarray(distribution.Cdf(X_INTERVAL[0])))
        high = float(np.asarray(distribution.Cdf(X_INTERVAL[1])))
        return TAIL_TOLERANCE - max(low, 1.0 - high)

    def mean_exists(primitive):
        """Grid proxy: integrability itself cannot be certified on a grid."""
        distribution = primitive.distribution
        mean = distribution.Mean()
        upper = distribution.PartialMean(X_INTERVAL[1])
        lower = distribution.PartialMean(X_INTERVAL[0])
        if not all(np.isfinite(value) for value in (mean, upper, lower)):
            return -np.inf
        return TAIL_MEAN_TOLERANCE - max(abs(upper), abs(mean - lower))

    def costly_on_average(primitive):
        return -float(primitive.distribution.Mean())

    def one_distribution_per_season(primitive):
        first = np.asarray(primitive.distribution.Cdf(probe), dtype=float)
        for _ in range(SEASON_DAYS - 1):
            if not np.array_equal(np.asarray(primitive.distribution.Cdf(probe),
                                             dtype=float), first):
                return -1.0
        return 0.0

    return {'present_bias_at_most_one': present_bias_at_most_one,
            'non_decreasing': non_decreasing,
            'bounded_in_unit_interval': bounded_in_unit_interval,
            'tails_reach_zero_and_one': tails_reach_zero_and_one,
            'mean_exists': mean_exists,
            'costly_on_average': costly_on_average,
            'one_distribution_per_season': one_distribution_per_season}


if __name__ == '__main__':
    Main()
