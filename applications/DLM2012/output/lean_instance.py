"""Writes lean/DLM2012.lean: the reduced ECU problem with its data as exact decimals, the
dual bounds, and feasible witnesses within a stated width of them."""
import json
import os
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reduction as rd  # noqa: E402

N = 2 * rd.REGIONS
SHRINK = Fraction(1, 10 ** 8)
CENTER_DIGITS = 12
MARGIN = Fraction(1, 10 ** 9)
WIDTH = Fraction(1, 10 ** 8)


def Exact(value):
    return Fraction(repr(float(value)))


def Decimal(fraction, digits=None):
    if digits is not None:
        fraction = Fraction(round(fraction * 10 ** digits), 10 ** digits)
    if fraction.denominator == 1:
        return f'({fraction.numerator} : ℝ)'
    return f'(({fraction.numerator} : ℝ) / {fraction.denominator})'


def Lit(fraction):
    """A decimal literal when the denominator has no prime factors but 2 and 5, else a quotient."""
    d, k = fraction.denominator, 0
    while 10 ** k % d:
        k += 1
        if k > 400:
            return Decimal(fraction)
    n = abs(fraction.numerator) * (10 ** k // d)
    sign = '-' if fraction < 0 else ''
    if k == 0:
        return f'({sign}{n} : ℝ)'
    digits = str(n).rjust(k + 1, '0')
    return f'({sign}{digits[:-k]}.{digits[-k:]} : ℝ)'


def Linear(coefficients, constant=Fraction(0), variable='z'):
    terms = [f'{Lit(c)} * {variable} {i}' for i, c in enumerate(coefficients) if c != 0]
    if constant != 0 or not terms:
        terms.append(Lit(constant))
    return ' + '.join(terms)


class Data:

    def __init__(self, certificates):
        reduced = certificates['reduced']
        self.fixed = Exact(reduced['fixed'])
        self.budget = Exact(reduced['budget'])
        self.weights = [Exact(w) for w in rd.WEIGHTS]
        self.targets = [Exact(m) for m in rd.TARGET_MOMENTS]
        self.b = [Exact(v) for v in rd.B_MAP]
        self.A = [[Exact(v) for v in row] for row in rd.A_MAP]
        self.low = [Exact(v) for v in rd.H_LOW]
        self.high = [Exact(v) for v in rd.H_HIGH]
        self.rows = [k for k in range(len(self.b)) if any(self.A[k])]
        self.inequalities, self.equalities = self.Constraints()

    def Constraints(self):
        ineq, eq = [], []
        R = rd.REGIONS
        for j in range(R):
            row = [Fraction(0)] * N
            row[j] = Fraction(1)
            ineq.append((row, Fraction(0)))
        for j in range(R):
            if self.low[j] == self.high[j]:
                row = [Fraction(0)] * N
                row[R + j], row[j] = Fraction(1), -self.low[j]
                eq.append((row, Fraction(0)))
                continue
            row = [Fraction(0)] * N
            row[R + j], row[j] = Fraction(1), -self.low[j]
            ineq.append((row, Fraction(0)))
            row = [Fraction(0)] * N
            row[R + j], row[j] = Fraction(-1), self.high[j]
            ineq.append((row, Fraction(0)))
        eq.append(([Fraction(1)] * R + [Fraction(0)] * R, Fraction(-1)))
        return ineq, eq

    def Criterion(self, z):
        total = self.fixed
        for k in range(len(self.b)):
            gap = self.targets[k] - self.b[k] - sum(a * x for a, x in zip(self.A[k], z))
            total += self.weights[k] * gap * gap
        return total

    def Feasible(self, z):
        return (all(sum(a * x for a, x in zip(row, z)) + c >= 0 for row, c in self.inequalities)
                and all(sum(a * x for a, x in zip(row, z)) + c == 0 for row, c in self.equalities)
                and self.Criterion(z) <= self.budget)

    def CriterionLean(self):
        terms = [Lit(self.fixed)]
        for k in range(len(self.b)):
            inner = Linear([-a for a in self.A[k]], self.targets[k] - self.b[k])
            terms.append(f'{Lit(self.weights[k])} * ({inner}) ^ 2')
        return ' +\n    '.join(terms)


def Witness(data, z_opt, z_min):
    """A decimal point of the reduced feasible set near the extremal point."""
    z = [(1 - SHRINK) * Exact(a) + SHRINK * Exact(b) for a, b in zip(z_opt, z_min)]
    z = [Fraction(round(x * 10 ** 16), 10 ** 16) for x in z]
    R = rd.REGIONS
    z = [max(x, Fraction(0)) if i < R else x for i, x in enumerate(z)]
    z[0] = 1 - sum(z[1:R])
    for j in range(R):
        if data.low[j] == data.high[j]:
            z[R + j] = data.low[j] * z[j]
        else:
            z[R + j] = min(max(z[R + j], data.low[j] * z[j]), data.high[j] * z[j])
    assert data.Feasible(z), 'witness infeasible'
    return z


def Objective(sign):
    return [Exact(v) for v in rd.ObjectiveVector(sign)]


def Value(objective, z):
    return sum(a * x for a, x in zip(objective, z))


def Bound(sign, data, center, objective, primal):
    """The dual bound: the primal value, rounded outward by MARGIN."""
    value = Exact(primal) + sign * MARGIN
    return Fraction(round(value * 10 ** 10 + sign * Fraction(1, 2)), 10 ** 10)


def Render(data, sides):
    lines = ['import DirectedSensitivity.Framework', '',
             '/-!',
             'DLM2012, ECU: the reduced problem.',
             '',
             'A primitive F is summarised by z = (p_0, ..., p_6, H_0, ..., H_6): the mass of F in',
             'the seven altruism regions and the integral of h* over each.',
             'The 22 ECU moments are affine in z and the criterion is the quadratic `crit`.',
             '`K` is the closure of the set of z that primitives in N produce, intersected with',
             'the budget.',
             'The extremal objectives `cU` and `cL` are the welfare effect on the face of the',
             'utility integrals that the maximum and the minimum select.',
             'All data are the double-precision values the Python code uses, written as exact',
             'decimals.',
             '-/', '',
             'namespace DLM2012', '', 'open DirectedSensitivity', '',
             '/-- The criterion as a function of the region statistics. -/',
             'def crit (z : Fin 14 → ℝ) : ℝ :=',
             '    ' + data.CriterionLean(), '',
             f'/-- The budget, the criterion of the refitted censored normal. -/',
             f'def qbar : ℝ := {Lit(data.budget)}', '',
             '/-- Closure of the reduced feasible set. -/',
             'def K : Set (Fin 14 → ℝ) :=',
             '  {z | crit z ≤ qbar ∧']
    conditions = [f'0 ≤ {Linear(row, c)}' for row, c in data.inequalities]
    conditions += [f'{Linear(row, c)} = 0' for row, c in data.equalities]
    lines.append('    ' + ' ∧\n    '.join(conditions) + '}')
    lines.append('')
    for name, sign, objective, center, bound, witness, low in sides:
        tag = 'U' if sign > 0 else 'L'
        lines += [f'/-- The welfare effect on the {"upper" if sign > 0 else "lower"} face. -/',
                  f'def c{tag} (z : Fin 14 → ℝ) : ℝ := {Linear(objective)}', '',
                  f'/-- The {"witness" } for the {name} side. -/',
                  f'def w{tag} : Fin 14 → ℝ := ![' + ', '.join(Lit(x) for x in witness) + ']', '']
        hints = ', '.join(f'sq_nonneg ({Linear(data.A[k], -sum(a * x for a, x in zip(data.A[k], center)))})'
                          for k in data.rows)
        relation = '≤' if sign > 0 else '≥'
        reached = '≥' if sign > 0 else '≤'
        lines += [f'theorem {name}_bound : ∀ z ∈ K, c{tag} z {relation} {Lit(bound)} := by',
                  '  intro z hz',
                  '  simp only [K, crit, Set.mem_ofPred_eq] at hz',
                  '  obtain ⟨hq, ' + ', '.join(f'h{i}' for i in range(len(data.inequalities) + len(data.equalities))) + '⟩ := hz',
                  '  simp only [qbar] at hq',
                  f'  simp only [c{tag}]',
                  f'  linarith [{hints}]', '',
                  f'theorem {name}_witness : w{tag} ∈ K ∧ c{tag} w{tag} {reached} {Lit(low)} := by',
                  f'  refine ⟨?_, ?_⟩',
                  f'  · simp only [K, crit, qbar, w{tag}, Set.mem_ofPred_eq, Matrix.cons_val, Fin.isValue]',
                  '    norm_num',
                  f'  · simp only [c{tag}, w{tag}, Matrix.cons_val, Fin.isValue]',
                  '    norm_num', '']
    upper, lower = sides[0], sides[1]
    lines += ['/-- The supremum of the welfare effect over the reduced feasible set lies in an',
              'interval of width 1e-8 around the reported value. -/',
              'theorem upper_enclosure (v : ℝ) (hv : IsUpperSensitivity K cU v) :',
              f'    {Lit(upper[6])} ≤ v ∧ v ≤ {Lit(upper[4])} := by',
              '  constructor',
              '  · exact le_trans upper_witness.2 (hv.1 ⟨wU, upper_witness.1, rfl⟩)',
              '  · exact hv.2 (by rintro _ ⟨z, hz, rfl⟩; exact upper_bound z hz)', '',
              'theorem lower_enclosure (v : ℝ) (hv : IsLowerSensitivity K cL v) :',
              f'    {Lit(lower[4])} ≤ v ∧ v ≤ {Lit(lower[6])} := by',
              '  constructor',
              '  · exact hv.2 (by rintro _ ⟨z, hz, rfl⟩; exact lower_bound z hz)',
              '  · exact le_trans (hv.1 ⟨wL, lower_witness.1, rfl⟩) lower_witness.2', '',
              'end DLM2012', '']
    return '\n'.join(lines)


def Main():
    with open(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'certificates.json')) as stream:
        certificates = json.load(stream)
    data = Data(certificates)
    reduced = certificates['reduced']
    fixed = reduced['fixed']
    z_min, _ = rd.MinimumCriterion(fixed)
    sides = []
    for name, sign in (('upper', 1), ('lower', -1)):
        item = reduced[name]
        objective = Objective(sign)
        center = [Fraction(round(Exact(x) * 10 ** CENTER_DIGITS), 10 ** CENTER_DIGITS) for x in item['z']]
        bound = Bound(sign, data, center, objective, item['primal'])
        witness = Witness(data, item['z'], z_min)
        value = Value(objective, witness)
        low = Fraction(round(value * 10 ** 10 - sign * Fraction(1, 2)), 10 ** 10)
        assert sign * (value - low) >= 0 and sign * (bound - low) <= WIDTH
        print(name, float(bound), float(value), float(low), float(data.budget - data.Criterion(witness)))
        sides.append((name, sign, objective, center, bound, witness, low))
    os.makedirs(os.path.join(HERE, 'lean'), exist_ok=True)
    with open(sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'lean', 'DLM2012.lean'), 'w') as stream:
        stream.write(Render(data, sides))


if __name__ == '__main__':
    Main()
