"""The analytic twin and the basis model for MMS2023."""
import os
import sys

import numpy as np
from scipy.optimize import fsolve
from scipy.special import expit

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import specification as spec  # noqa: E402

LOW_ATOM, HIGH_ATOM, GAP_SHARE = -24.0, 25.0, 0.5
NARROW_CELL = 1e-4
FIXED_POINT_TOLERANCE = 1e-13
FIXED_POINT_ITERATIONS = 400
ROUNDING_MASS = 1e-12


class Discrete(spec.NetBenefitDistribution):

    def __init__(self, atoms, masses):
        self.atoms = np.asarray(atoms, dtype=float)
        self.masses = np.asarray(masses, dtype=float)

    def Cdf(self, x):
        x = np.asarray(x, dtype=float)
        order = np.argsort(self.atoms, kind='stable')
        cumulative = np.concatenate([[0.0], np.cumsum(self.masses[order])])
        return np.minimum(cumulative[np.searchsorted(self.atoms[order], x, side='right')],
                          1.0)

    def PartialMean(self, lower):
        above = self.atoms > lower
        return float(self.masses[above] @ self.atoms[above])

    def Mean(self):
        return float(self.masses @ self.atoms)


class Mixture(spec.NetBenefitDistribution):

    def __init__(self, weight, first, second):
        self.weight, self.first, self.second = float(weight), first, second

    def Cdf(self, x):
        return np.minimum(self.weight * self.first.Cdf(x)
                          + (1.0 - self.weight) * self.second.Cdf(x), 1.0)

    def PartialMean(self, lower):
        return (self.weight * self.first.PartialMean(lower)
                + (1.0 - self.weight) * self.second.PartialMean(lower))

    def Mean(self):
        return self.weight * self.first.Mean() + (1.0 - self.weight) * self.second.Mean()


def Isotonic(hazards, weights):
    """Weighted pool-adjacent-violators: the nondecreasing hazards of best Bernoulli fit."""
    blocks = []
    for value, weight in zip(hazards, weights):
        blocks.append([float(value), float(weight), 1])
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            upper, lower = blocks.pop(), blocks.pop()
            total = upper[1] + lower[1]
            blocks.append([(upper[0] * upper[1] + lower[0] * lower[1]) / total, total,
                           upper[2] + lower[2]])
    return np.concatenate([[value] * size for value, _, size in blocks])


def NegativeLogLikelihood(hazards):
    hazards = np.asarray(hazards, dtype=float)
    return -float(spec.FILED @ np.log(hazards) + (spec.AT_RISK - spec.FILED) @ np.log1p(-hazards))


def ExpectedMaximum(atoms, masses, value):
    return float(masses @ np.maximum(atoms, value))


def TwinAtoms(present_bias, values, masses):
    middle = present_bias * ((1.0 - GAP_SHARE) * values[1:] + GAP_SHARE * values[:-1])
    return np.concatenate([[LOW_ATOM], middle, [HIGH_ATOM]])


def Twin(present_bias, hazards):
    """The analytic anchor: a discrete law whose filing hazards under ``present_bias`` are
    ``hazards``, found as the least fixed point of a monotone map by iteration from below.

    ``values[k]`` is the continuation value k days before the deadline, so the threshold
    on day 64 - k is ``present_bias * values[k]``.
    """
    remaining = np.asarray(hazards, dtype=float)[::-1]
    masses = np.concatenate([[1.0 - remaining[0]], remaining[:-1] - remaining[1:],
                             [remaining[-1]]])
    values = np.full(len(remaining), LOW_ATOM)
    for _ in range(100000):
        atoms = TwinAtoms(present_bias, values, masses)
        updated = np.empty_like(values)
        updated[0] = float(masses @ atoms)
        for k in range(1, len(values)):
            updated[k] = ExpectedMaximum(atoms, masses, updated[k - 1])
        if np.array_equal(updated, values):
            break
        values = updated
    else:
        raise RuntimeError('the twin iteration did not settle')
    return spec.Primitive(present_bias, Discrete(TwinAtoms(present_bias, values, masses),
                                                 masses)), values


def _LogisticPartialMean(location, x):
    x = np.asarray(x, dtype=float)
    out = np.where(np.isneginf(x), location, 0.0)
    finite = np.isfinite(x)
    z = x[finite] - location
    out[finite] = location * expit(-z) + z * expit(-z) + np.logaddexp(0.0, -z)
    return out


def _LogisticCdf(location, x):
    return expit(np.asarray(x, dtype=float) - location)


class BasisModel:
    """The basis model, in coordinates x = (present_bias, kappa, weight, pull).

    The law is (1 - weight) G_kappa + weight H, where H spreads G_kappa over the cells cut
    by the thresholds and the continuation values, one atom per cell at its conditional
    mean, and rescales the threshold cells so that the hazards equal
    p(present_bias, kappa) + weight * pull * (isotonic - p(present_bias, kappa)).
    """

    def __init__(self, isotonic):
        self.isotonic = np.asarray(isotonic, dtype=float)
        self._cache = {}

    def WorkingHazards(self, present_bias, kappa):
        return spec.FilingProbabilities(spec.Primitive(present_bias, spec.Logistic(-kappa)))

    def TargetHazards(self, x):
        present_bias, kappa, weight, pull = (float(v) for v in x)
        share = weight * pull
        return (1.0 - share) * self.WorkingHazards(present_bias, kappa) + share * self.isotonic

    def Complement(self, x, values, remaining):
        present_bias, kappa, weight, _ = (float(v) for v in x)
        location = -kappa
        thresholds = present_bias * values
        cuts = np.unique(np.concatenate([thresholds, values]))
        lower = np.concatenate([[-np.inf], cuts])
        upper = np.concatenate([cuts, [np.inf]])
        cell = np.searchsorted(thresholds, upper, side='left')
        mass = _LogisticCdf(location, upper) - _LogisticCdf(location, lower)
        mass = np.where(np.isinf(upper), 1.0 - _LogisticCdf(location, lower), mass)
        first = _LogisticPartialMean(location, lower) - _LogisticPartialMean(location, upper)
        with np.errstate(invalid='ignore', divide='ignore'):
            atoms = first / mass
        narrow = np.isfinite(lower) & np.isfinite(upper) & (upper - lower < NARROW_CELL)
        atoms = np.where(narrow, 0.5 * (lower + upper), atoms)
        atoms = np.clip(atoms, lower, upper)
        target = np.concatenate([[1.0 - remaining[0]], remaining[:-1] - remaining[1:],
                                 [remaining[-1]]])
        cell_mass = np.bincount(cell, weights=mass, minlength=len(target))
        scale = (target - (1.0 - weight) * cell_mass) / weight / cell_mass
        return atoms, mass * scale[cell]

    def _Update(self, x, values, remaining):
        _, kappa, weight, _ = (float(v) for v in x)
        logistic = spec.Logistic(-kappa)
        atoms, masses = self.Complement(x, values, remaining)
        updated = np.empty_like(values)
        updated[0] = (1.0 - weight) * (-kappa) + weight * float(masses @ atoms)
        for k in range(1, len(values)):
            v = updated[k - 1]
            working = v * float(logistic.Cdf(v)) + logistic.PartialMean(v)
            updated[k] = (1.0 - weight) * working + weight * ExpectedMaximum(atoms, masses, v)
        return updated

    def Solve(self, x):
        """The primitive at x, its continuation values, and its complement masses."""
        key = tuple(float(v) for v in x)
        if key in self._cache:
            return self._cache[key]
        present_bias, kappa, weight, _ = key
        logistic = spec.Logistic(-kappa)
        working = spec.Primitive(present_bias, logistic)
        values = spec.Thresholds(working)[::-1] / present_bias
        if weight <= 0.0:
            result = (working, values, np.zeros(1))
            self._cache[key] = result
            return result
        remaining = self.TargetHazards(key)[::-1]
        for _ in range(FIXED_POINT_ITERATIONS):
            updated = self._Update(key, values, remaining)
            if np.max(np.abs(updated - values)) <= FIXED_POINT_TOLERANCE:
                values = updated
                break
            values = updated
        else:
            values = fsolve(lambda v: self._Update(key, v, remaining) - v, values,
                            xtol=1e-14)
        atoms, masses = self.Complement(key, values, remaining)
        primitive = spec.Primitive(present_bias, Mixture(weight, Discrete(atoms, _Clean(masses)),
                                                         logistic))
        result = (primitive, values, masses)
        self._cache[key] = result
        return result

    def Primitive(self, x):
        return self.Solve(x)[0]

    def MassFloor(self, x):
        """The least complement mass, less a rounding allowance of ROUNDING_MASS."""
        masses = self.Solve(x)[2]
        if not np.all(np.isfinite(masses)):
            return -np.inf
        return float(masses.min()) + ROUNDING_MASS


def _Clean(masses):
    """Masses with rounding-level negatives set to zero and the total rescaled to one."""
    masses = np.where(masses < 0.0, 0.0, masses)
    return masses / masses.sum()
