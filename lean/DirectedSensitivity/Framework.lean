import Mathlib

/-!
Shared vocabulary for formalizing directed-sensitivity claims.

A problem has a set of feasible primitives in the nesting model and an object of
interest `c : N → ℝ`. The upper (lower) sensitivity bound is the supremum
(infimum) of `c` over the feasible set. The workflow's central claims are that a
given value is that supremum or infimum; this file states them as `IsLUB`/`IsGLB`
and gives the two natural ways to discharge them — an attained extremum, or a
bound that is approached.
-/

namespace DirectedSensitivity

open Set

variable {N ι : Type*}

/-- The feasible set: primitives in the nesting model whose criterion `q` does not
exceed the `budget` and that satisfy every constraint (each a nonnegative quantity
on the nesting model). Mirrors `SensitivityProblem.IsFeasible` at zero tolerance. -/
def feasibleSet (nesting : Set N) (q : N → ℝ) (budget : ℝ)
    (constraints : ι → N → ℝ) : Set N :=
  {ν | ν ∈ nesting ∧ q ν ≤ budget ∧ ∀ i, 0 ≤ constraints i ν}

/-- The objective values attained over a feasible set. -/
def objectiveValues (feasible : Set N) (c : N → ℝ) : Set ℝ := c '' feasible

/-- `v` is the upper sensitivity bound: the least upper bound of the object of
interest `c` over the feasible set. -/
def IsUpperSensitivity (feasible : Set N) (c : N → ℝ) (v : ℝ) : Prop :=
  IsLUB (objectiveValues feasible c) v

/-- `v` is the lower sensitivity bound: the greatest lower bound of `c` over the
feasible set. -/
def IsLowerSensitivity (feasible : Set N) (c : N → ℝ) (v : ℝ) : Prop :=
  IsGLB (objectiveValues feasible c) v

/-- Certify the upper bound by an attained maximum: a feasible primitive reaches
`v` and no feasible primitive exceeds it. -/
theorem isUpperSensitivity_of_isGreatest {feasible : Set N} {c : N → ℝ} {v : ℝ}
    (h : IsGreatest (objectiveValues feasible c) v) :
    IsUpperSensitivity feasible c v :=
  h.isLUB

/-- Certify the lower bound by an attained minimum. -/
theorem isLowerSensitivity_of_isLeast {feasible : Set N} {c : N → ℝ} {v : ℝ}
    (h : IsLeast (objectiveValues feasible c) v) :
    IsLowerSensitivity feasible c v :=
  h.isGLB

/-- Certify the upper bound when it is not attained but approached: every feasible
primitive satisfies `c ν ≤ v`, and for every `ε > 0` some feasible primitive gets
within `ε` of `v` from below. -/
theorem isUpperSensitivity_of_le_of_approx {feasible : Set N} {c : N → ℝ} {v : ℝ}
    (hub : ∀ ν ∈ feasible, c ν ≤ v)
    (happrox : ∀ ε > 0, ∃ ν ∈ feasible, v - ε < c ν) :
    IsUpperSensitivity feasible c v := by
  constructor
  · rintro _ ⟨ν, hν, rfl⟩
    exact hub ν hν
  · intro w hw
    by_contra h
    push Not at h
    obtain ⟨ν, hν, hlt⟩ := happrox (v - w) (by linarith)
    have hle : c ν ≤ w := hw ⟨ν, hν, rfl⟩
    linarith

/-- Certify the lower bound when it is approached: every feasible primitive
satisfies `v ≤ c ν`, and for every `ε > 0` some feasible primitive gets within `ε`
of `v` from above. -/
theorem isLowerSensitivity_of_ge_of_approx {feasible : Set N} {c : N → ℝ} {v : ℝ}
    (hlb : ∀ ν ∈ feasible, v ≤ c ν)
    (happrox : ∀ ε > 0, ∃ ν ∈ feasible, c ν < v + ε) :
    IsLowerSensitivity feasible c v := by
  constructor
  · rintro _ ⟨ν, hν, rfl⟩
    exact hlb ν hν
  · intro w hw
    by_contra h
    push Not at h
    obtain ⟨ν, hν, hlt⟩ := happrox (w - v) (by linarith)
    have hge : w ≤ c ν := hw ⟨ν, hν, rfl⟩
    linarith

end DirectedSensitivity
