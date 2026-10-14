import io
import json
import os
import sys
import zipfile

import numpy as np
from scipy.integrate import quad
from scipy.io import loadmat
from scipy.optimize import brentq
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

GRIDS      = config.CONFIG['grids']
TOLERANCES = config.CONFIG['tolerances']

GS = 10.0

# Table IV, column 1. The table reports social pressure as s * GS.
SCALARS = {'h0': 0.414, 'r': 0.341, 'eta': 0.040, 'gi': 10.606,
           's_lar': 3.751 / GS, 's_ecu': 1.438 / GS, 's_svy': 6.197,
           'mu_s': -26.863, 'sigma_s': 29.591, 'timeval': 80.656, 'h_2009': 0.414}
ALTRUISM = {'lar': {'share_at_zero': 0.753, 'conditional_mean': 12.786},
            'ecu': {'share_at_zero': 0.763, 'conditional_mean': 9.659}}

CHARITIES    = ('lar', 'ecu')
BINS         = ((10.0, 20.0), (20.0, 50.0), (50.0, np.inf))
TREATMENTS_2008 = (('NF0d10m', 0.0, 10, 'NF'), ('F0d5m', 0.0, 5, 'F'),
                   ('F0d10m', 0.0, 10, 'F'), ('F10d10m', 10.0, 10, 'F'))
TREATMENTS_2009 = (('NF_0d5m_09', 0.0, 5, 'NF'), ('NF_5d5m_09', 5.0, 5, 'NF'),
                   ('SOO_0d5m_09', 0.0, 5, 'SOO'), ('SOO_5d5m_09', 5.0, 5, 'SOO'),
                   ('F_0d10m_09', 0.0, 10, 'F'), ('F_0d5m_09', 0.0, 5, 'F'),
                   ('F_10d5m_09', 10.0, 5, 'F'), ('F_5d5m_09', 5.0, 5, 'F'))

PARAMETER_NAMES = ('h0', 'r', 'eta', 'gi', 's_lar', 's_ecu',
                   's_svy', 'mu_s', 'sigma_s', 'timeval', 'h_2009')

MOMENT_ORDER = (
    [('lar', 'PH_NF'), ('ecu', 'PH_NF'), ('lar', 'PH_F'), ('ecu', 'PH_F'),
     ('lar', 'PH_OO'), ('ecu', 'PH_OO')]
    + [('survey', k) for k in ('PH_NF0d10m', 'PH_F0d5m', 'PH_F0d10m', 'PH_F10d10m')]
    + [('lar', 'PG_NF'), ('ecu', 'PG_NF'), ('lar', 'PG_F'), ('ecu', 'PG_F'),
       ('lar', 'PG_OO'), ('ecu', 'PG_OO')]
    + [('survey', k) for k in ('PSV_NF0d10m', 'PSV_F0d5m', 'PSV_F0d10m', 'PSV_F10d10m')]
    + [('lar', 'POO_OO'), ('ecu', 'POO_OO')]
    + [('survey', k) for k in (
        'PH_NF_0d5m_09', 'PH_NF_5d5m_09', 'PH_SOO_0d5m_09', 'PH_SOO_5d5m_09',
        'PH_F_0d10m_09', 'PH_F_0d5m_09', 'PH_F_10d5m_09', 'PH_F_5d5m_09',
        'PSV_NF_0d5m_09', 'PSV_NF_5d5m_09', 'PSV_SOO_0d5m_09', 'PSV_SOO_5d5m_09',
        'PSV_F_0d10m_09', 'PSV_F_0d5m_09', 'PSV_F_10d5m_09', 'PSV_F_5d5m_09',
        'POO_SOO_0d5m_09', 'POO_SOO_5d5m_09')]
    + [(charity, f'{treatment}_{bin_name}')
       for bin_name in ('0gs', 'gs', 'int', 'int1', 'int2')
       for treatment in ('PG_NF', 'PG_OO', 'PG_F')
       for charity in ('ecu', 'lar')])

BLOCK = {name: np.array([block == name for block, _ in MOMENT_ORDER])
         for name in ('lar', 'ecu', 'survey')}


def LoadMoments(path=config.PATHS['data'], member=config.CONFIG['members']['data']):
    with zipfile.ZipFile(path) as archive:
        return loadmat(io.BytesIO(archive.read(member)))


MOMENTS   = LoadMoments()
EMPIRICAL = MOMENTS['Moments'].ravel()
WEIGHT    = np.linalg.inv(np.diag(np.diag(MOMENTS['VCcontrol'])))


class OutsideTheClass(ValueError):
    pass


class Parameters:

    def __init__(self, h0, r, eta, gi, s_lar, s_ecu, s_svy, mu_s, sigma_s,
                 timeval, h_2009):
        self.h0, self.r, self.eta, self.gi = float(h0), float(r), float(eta), float(gi)
        self.s_lar, self.s_ecu, self.s_svy = float(s_lar), float(s_ecu), float(s_svy)
        self.mu_s, self.sigma_s = float(mu_s), float(sigma_s)
        self.timeval, self.h_2009 = float(timeval), float(h_2009)

    def Pressure(self, charity):
        return self.s_lar if charity == 'lar' else self.s_ecu

    def AsVector(self):
        return np.array([getattr(self, name) for name in PARAMETER_NAMES])

    @classmethod
    def FromVector(cls, vector):
        return cls(*np.asarray(vector, dtype=float).ravel())

    def Check(self):
        if not 0.0 < self.gi:
            raise OutsideTheClass('gi must be positive: log(gi) appears in every threshold')
        if not all(0.0 <= s < 1.0 for s in (self.s_lar, self.s_ecu)):
            raise OutsideTheClass('only the 0 <= S < 1 branch of the Matlab code is ported')
        if not (0.0 <= self.h0 <= 1.0 and 0.0 <= self.h_2009 <= 1.0
                and 0.0 <= self.r <= 1.0):
            raise OutsideTheClass('h0, h_2009 and r are probabilities')
        if not (self.eta > 0.0 and self.sigma_s > 0.0 and self.s_svy >= 0.0
                and self.timeval >= 0.0):
            raise OutsideTheClass('eta and sigma_s must be positive; s_svy, timeval '
                                  'non-negative')
        return self


class AltruismDistribution:

    a_max = None

    @property
    def Atom(self):
        raise NotImplementedError

    def Density(self, a):
        raise NotImplementedError

    def Cdf(self, a):
        if a <= 0.0:
            return 0.0 if a < 0.0 else self.Atom
        return self.Atom + self._Continuous(lambda x: 1.0, 0.0, min(a, self.a_max))

    def Mass(self, lower, upper):
        return self.Integrate(lambda a: np.ones_like(a), lower, upper)

    def Integrate(self, psi, lower, upper, breakpoints=()):
        upper = min(upper, self.a_max)
        if upper <= lower:
            return 0.0
        total = self.Atom * float(psi(0.0)) if lower <= 0.0 <= upper else 0.0
        return total + self._Continuous(psi, max(lower, 0.0), upper, breakpoints)

    def _Continuous(self, psi, lower, upper, breakpoints=()):
        edges = [lower] + sorted(b for b in breakpoints if lower < b < upper) + [upper]

        def integrand(a):
            return float(psi(a)) * float(self.Density(a))

        total = 0.0
        for left, right in zip(edges[:-1], edges[1:]):
            if right > left:
                total += quad(integrand, left, right, epsabs=1e-12, limit=200)[0]
        return total


class CensoredNormal(AltruismDistribution):

    def __init__(self, mu, sigma, a_max):
        self.mu, self.sigma, self.a_max = float(mu), float(sigma), float(a_max)
        self._atom = float(norm.cdf(-self.mu / self.sigma))
        self._tail = float(norm.sf((self.a_max - self.mu) / self.sigma))

    @property
    def Atom(self):
        """Mass at zero, plus the 2.9e-7 above a_max, so the distribution is proper."""
        return self._atom + self._tail

    def Density(self, a):
        return norm.pdf(np.asarray(a, dtype=float), self.mu, self.sigma)

    def Cdf(self, a):
        if a < 0.0:
            return 0.0
        if a >= self.a_max:
            return 1.0
        return float(norm.cdf(a, self.mu, self.sigma)) + self._tail


class GridDistribution(AltruismDistribution):

    def __init__(self, atom, nodes, values, a_max):
        self._atom, self.a_max = float(atom), float(a_max)
        self.nodes = np.asarray(nodes, dtype=float)
        self.values = np.asarray(values, dtype=float)

    @property
    def Atom(self):
        return self._atom

    def Density(self, a):
        return np.interp(np.asarray(a, dtype=float), self.nodes, self.values,
                         left=0.0, right=0.0)

    def _Continuous(self, psi, lower, upper, breakpoints=()):
        edges = [lower] + sorted(b for b in breakpoints if lower < b < upper) + [upper]
        abscissa, weights = np.polynomial.legendre.leggauss(GRIDS['quadrature_nodes'])
        total = 0.0
        for left, right in zip(edges[:-1], edges[1:]):
            if right <= left:
                continue
            middle, half = 0.5 * (left + right), 0.5 * (right - left)
            points = middle + half * abscissa
            total += half * float(weights @ (np.asarray(psi(points), dtype=float)
                                             * self.Density(points)))
        return total


class Primitive:

    def __init__(self, parameters, lar, ecu):
        self.parameters = parameters
        self.lar, self.ecu = lar, ecu

    def __getitem__(self, charity):
        return getattr(self, charity)


def RecoverShape(share_at_zero, conditional_mean):
    """(mu, sigma) from the share at a = 0 and the mean of a conditional on a > 0."""
    alpha = norm.ppf(share_at_zero)
    ratio = norm.pdf(alpha) / (1.0 - share_at_zero)
    sigma = conditional_mean / (ratio - alpha)
    return -alpha * sigma, sigma


BENCHMARK       = Parameters(**SCALARS)
BENCHMARK_SHAPE = {charity: RecoverShape(**ALTRUISM[charity])
                   for charity in CHARITIES}
A_MAX = {charity: mu + 5.0 * sigma for charity, (mu, sigma) in BENCHMARK_SHAPE.items()}

PUBLISHED_MODEL = Primitive(BENCHMARK,
                            lar=CensoredNormal(*BENCHMARK_SHAPE['lar'], A_MAX['lar']),
                            ecu=CensoredNormal(*BENCHMARK_SHAPE['ecu'], A_MAX['ecu']))


def Main():
    print(f'scalars            {BENCHMARK.AsVector()!r}')
    for charity in CHARITIES:
        print(f'{charity} mu, sigma, a_max  {BENCHMARK_SHAPE[charity]!r}, {A_MAX[charity]!r}')
        print(f'{charity} thresholds        {Thresholds(BENCHMARK, charity)!r}')
    baseline = Baseline()
    print(f'criterion budget   {Criterion(baseline)!r}')
    print(f'target lar         {TargetLaRabida(baseline)!r}')
    print(f'target ecu         {TargetECU(baseline)!r}')


def Baseline(path=None):
    path = config.output('estimate.json') if path is None else path
    with open(path) as stream:
        estimate = json.load(stream)
    return Primitive(BENCHMARK, **{
        charity: CensoredNormal(estimate[charity]['mu'], estimate[charity]['sigma'],
                                A_MAX[charity])
        for charity in CHARITIES})


def Criterion(primitive):
    gap = EMPIRICAL - Moments(primitive)
    return float(gap @ WEIGHT @ gap)


def Moments(primitive):
    parameters = primitive.parameters.Check()
    blocks = {'survey': SurveyMoments(parameters)}
    for charity in CHARITIES:
        blocks[charity] = CharityMoments(primitive[charity], parameters, charity)
    return np.array([blocks[block][name] for block, name in MOMENT_ORDER])


def SurveyMoments(parameters):
    return SurveyBlock(parameters.h0, parameters.r, parameters.eta, parameters.mu_s,
                       parameters.sigma_s, parameters.s_svy, parameters.timeval,
                       parameters.h_2009)


def SurveyBlock(h0, r, eta, mu_s, sigma_s, s_svy, timeval, h_2009):
    block = _Block(h0, r, eta, mu_s, sigma_s, s_svy, timeval, TREATMENTS_2008)
    block.update(_Block(h_2009, r, eta, mu_s, sigma_s, s_svy, timeval, TREATMENTS_2009))
    return block


def _Block(h0, r, eta, mu_s, sigma_s, s_svy, timeval, treatments):
    out = {}
    for name, pay, minutes, kind in treatments:
        time_cost = 0.0 if minutes == 10 else -timeval * minutes / 60.0
        s_yes = time_cost - pay - s_svy
        s_zero = time_cost - pay
        s_h0 = time_cost - pay - h0 / eta
        s_h1 = time_cost - pay + (1.0 - h0) / eta

        def cdf(x):
            return float(norm.cdf((x - mu_s) / sigma_s))

        def pdf(x):
            return float(norm.pdf((x - mu_s) / sigma_s))

        def truncated(lower):
            return (mu_s * (cdf(s_h1) - cdf(lower))
                    + sigma_s * (pdf(lower) - pdf(s_h1)))

        completes = 1.0 - cdf(s_yes)
        if kind == 'NF':
            out[f'PH_{name}'] = h0
            out[f'PSV_{name}'] = h0 * completes
            continue

        threshold = s_zero if kind == 'SOO' else s_yes
        floor = h0 - eta * s_svy
        if kind == 'F' and floor <= 0.0:
            threshold = s_h0
        below = (r * floor * cdf(s_yes)) if kind == 'F' and floor > 0.0 else 0.0

        sorted_in = (r * (h0 + eta * (pay - time_cost)) * (cdf(s_h1) - cdf(threshold))
                     + r * eta * truncated(threshold) + r * (1.0 - cdf(s_h1)))
        out[f'PH_{name}'] = (1 - r) * h0 + below + sorted_in
        out[f'PSV_{name}'] = (1 - r) * h0 * completes + sorted_in
        if kind == 'SOO':
            out[f'POO_{name}'] = r * h0 * cdf(s_zero)
    return out


def CharityMoments(distribution, parameters, charity):
    h0, r, a_max = parameters.h0, parameters.r, distribution.a_max
    a_minmin, a_min, a_plus, a_0 = Thresholds(parameters, charity)
    kinks = (a_minmin, a_min, a_plus, a_0)

    def home(lower, upper, opt_out=False):
        return distribution.Integrate(
            lambda a: HomeStar(a, parameters, charity, opt_out), lower, upper, kinks)

    survives = 1.0 - distribution.Cdf(a_minmin)
    moments = {
        'PH_NF': h0,
        'PH_F': (1 - r) * h0 + r * home(-np.inf, a_max),
        'PH_OO': (1 - r) * h0 + r * home(-np.inf, a_max, opt_out=True),
        'PG_NF': h0 * survives,
        'PG_F': (1 - r) * h0 * survives + r * home(a_minmin, a_max),
        'POO_OO': r * h0 * (distribution.Cdf(a_0) if np.isfinite(a_0) else 0.0),
    }
    moments['PG_OO'] = moments['PH_OO'] - (1 - r) * h0 + h0 * (1 - r) * survives

    intervals = [('0gs', (a_minmin, a_min)), ('gs', (a_min, a_plus))]
    intervals += [(name, _BinEdges(low, high, parameters, charity))
                  for name, (low, high) in zip(('int', 'int1', 'int2'), BINS)]
    for name, (lower, upper) in intervals:
        mass = distribution.Cdf(min(upper, a_max)) - distribution.Cdf(max(lower, 0.0))
        opt_out_lower = min(max(lower, a_0), upper) if np.isfinite(a_0) else lower
        moments[f'PG_NF_{name}'] = h0 * mass
        moments[f'PG_F_{name}'] = (1 - r) * h0 * mass + r * home(lower, upper)
        moments[f'PG_OO_{name}'] = (1 - r) * h0 * mass + r * home(opt_out_lower, upper)
    return moments


def Thresholds(parameters, charity):
    """(a_minmin, a_min, a_plus, a_0): refuse below, give GS at, give more above, seek above."""
    s, gi = parameters.Pressure(charity), parameters.gi
    a_minmin = (1.0 - s) * gi
    a_min = (1.0 - s) * (gi + GS)
    a_plus = gi + GS
    a_0 = GS / np.log1p(GS / gi)
    if s == 0.0:
        return a_minmin, a_min, a_plus, -np.inf
    if not a_min <= a_0 <= a_plus:
        def net(a):
            return a * np.log(a / ((1.0 - s) * gi)) - a + (1.0 - s) * gi - s * GS
        a_0 = brentq(net, a_minmin - 1e-5, a_min + 1e-5, xtol=1e-14, rtol=1e-15)
    return a_minmin, a_min, a_plus, a_0


def GivingStar(a, parameters, charity):
    a = np.asarray(a, dtype=float)
    s, gi = parameters.Pressure(charity), parameters.gi
    a_minmin, a_min, a_plus, _ = Thresholds(parameters, charity)
    return np.select([a <= a_minmin, a <= a_min, a <= a_plus],
                     [0.0, a / (1.0 - s) - gi, GS], default=a - gi)


def HomeStar(a, parameters, charity, opt_out=False):
    a = np.asarray(a, dtype=float)
    s, gi, eta = parameters.Pressure(charity), parameters.gi, parameters.eta
    a_minmin, a_min, a_plus, a_0 = Thresholds(parameters, charity)
    safe = np.where(a > 0, a, 1.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        interior = (a * np.log(safe / (1.0 - s)) - a * (1.0 + np.log(gi))
                    - s * GS + gi * (1.0 - s))
        bunching = a * np.log1p(GS / gi) - GS
        above = gi - a + a * np.log(safe / gi)
    h = parameters.h0 + eta * np.select(
        [a <= a_minmin, a <= a_min, a <= a_plus], [-s * GS, interior, bunching],
        default=above)
    if opt_out:
        h = np.where(a <= a_0, 0.0, h)
    return np.clip(h, 0.0, 1.0)


def _BinEdges(low, high, parameters, charity):
    s, gi = parameters.Pressure(charity), parameters.gi
    a_minmin, a_min, a_plus, _ = Thresholds(parameters, charity)

    def edge(g, upper):
        if g <= 0.0:
            return a_minmin if upper else -np.inf
        if g < GS:
            return (g + gi) * (1.0 - s)
        if g == GS:
            return a_plus
        return g + gi if np.isfinite(g) else np.inf

    return edge(low, False), edge(high, True)


def TargetLaRabida(primitive):
    return Welfare(primitive, 'lar')


def TargetECU(primitive):
    return Welfare(primitive, 'ecu')


def Welfare(primitive, charity):
    parameters = primitive.parameters.Check()
    distribution = primitive[charity]
    return parameters.h0 * distribution.Integrate(
        lambda a: VisitUtility(a, parameters, charity), -np.inf, distribution.a_max,
        Thresholds(parameters, charity))


def VisitUtility(a, parameters, charity):
    a = np.asarray(a, dtype=float)
    s, gi = parameters.Pressure(charity), parameters.gi
    g = GivingStar(a, parameters, charity)
    return -g + a * np.log(g + gi) - s * (GS - g) * (g < GS) - a * np.log(gi)


def ConstraintSet(charity, density_points=GRIDS['density_points']):
    grid = np.linspace(0.0, A_MAX[charity], density_points)[1:]

    def atom_non_negative(primitive):
        return float(primitive[charity].Atom)

    def density_non_negative(primitive):
        distribution = primitive[charity]
        points = grid
        nodes = getattr(distribution, 'nodes', None)
        if nodes is not None:
            points = np.union1d(points, np.asarray(nodes, dtype=float))
            points = points[(points > 0.0) & (points <= distribution.a_max)]
        return float(np.min(np.asarray(distribution.Density(points), dtype=float)))

    def total_mass_one(primitive):
        distribution = primitive[charity]
        return float(TOLERANCES['mass']
                     - abs(distribution.Mass(0.0, distribution.a_max) - 1.0))

    def support_bounded(primitive):
        return float(A_MAX[charity] - primitive[charity].a_max)

    def probabilities_in_unit_interval(primitive):
        scalars = primitive.parameters
        return float(min(min(x, 1.0 - x)
                         for x in (scalars.h0, scalars.h_2009, scalars.r)))

    def pressure_in_unit_interval(primitive):
        scalars = primitive.parameters
        return float(min(min(s, 1.0 - TOLERANCES['strict_floor'] - s)
                         for s in (scalars.s_lar, scalars.s_ecu)))

    def warm_glow_positive(primitive):
        return float(primitive.parameters.gi - TOLERANCES['strict_floor'])

    def presence_elasticity_positive(primitive):
        return float(primitive.parameters.eta - TOLERANCES['strict_floor'])

    def survey_scale_non_negative(primitive):
        scalars = primitive.parameters
        return float(min(scalars.sigma_s - TOLERANCES['strict_floor'], scalars.s_svy, scalars.timeval))

    return {'atom_non_negative': atom_non_negative,
            'density_non_negative': density_non_negative,
            'total_mass_one': total_mass_one,
            'support_bounded': support_bounded,
            'probabilities_in_unit_interval': probabilities_in_unit_interval,
            'pressure_in_unit_interval': pressure_in_unit_interval,
            'warm_glow_positive': warm_glow_positive,
            'presence_elasticity_positive': presence_elasticity_positive,
            'survey_scale_non_negative': survey_scale_non_negative}


if __name__ == '__main__':
    Main()
