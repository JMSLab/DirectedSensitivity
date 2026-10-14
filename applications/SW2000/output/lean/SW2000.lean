import DirectedSensitivity.Framework

/-!
Stock and Wright (2000): directed sensitivity of the discount factor.
Claim 1 (the per-gap increment bound), Claim 2 (the upper value is one), and
Claim 4 (the basis model mixes within the nesting model).
The numerical fact that the reported upper certificate fits the data, `q ν₀ ≤ budget`,
enters Claim 2 as the hypothesis `hfeas`.
-/

namespace SW2000

open DirectedSensitivity Set Real

/-- Claim 1: with `Φ = -log u'` and `Φ'(s) = ρ(s)/s`, `0 < ρ ≤ 10`, the increment of `Φ`
over `[a, b]` lies in `(0, 10 log(b/a)]`. -/
theorem increment_bounds {Φ ρ : ℝ → ℝ} {a b : ℝ} (ha : 0 < a) (hab : a < b)
    (hderiv : ∀ s ∈ Icc a b, HasDerivAt Φ (ρ s / s) s)
    (hpos : ∀ s ∈ Icc a b, 0 < ρ s) (hcap : ∀ s ∈ Icc a b, ρ s ≤ 10) :
    0 < Φ b - Φ a ∧ Φ b - Φ a ≤ 10 * (log b - log a) := by
  have hcont : ContinuousOn Φ (Icc a b) := fun s hs => (hderiv s hs).continuousAt.continuousWithinAt
  have hs_pos : ∀ s ∈ interior (Icc a b), 0 < s := by
    intro s hs
    rw [interior_Icc] at hs
    exact lt_trans ha hs.1
  have hin : ∀ s ∈ interior (Icc a b), s ∈ Icc a b := fun s hs => interior_subset hs
  constructor
  · have hmono : StrictMonoOn Φ (Icc a b) := by
      apply strictMonoOn_of_deriv_pos (convex_Icc a b) hcont
      intro s hs
      rw [(hderiv s (hin s hs)).deriv]
      exact div_pos (hpos s (hin s hs)) (hs_pos s hs)
    have := hmono (left_mem_Icc.2 hab.le) (right_mem_Icc.2 hab.le) hab
    linarith
  · set G : ℝ → ℝ := fun s => 10 * log s - Φ s with hG
    have hGderiv : ∀ s ∈ Icc a b, HasDerivAt G (10 * s⁻¹ - ρ s / s) s := by
      intro s hs
      have hs0 : s ≠ 0 := ne_of_gt (lt_of_lt_of_le ha hs.1)
      exact ((hasDerivAt_log hs0).const_mul 10).sub (hderiv s hs)
    have hGcont : ContinuousOn G (Icc a b) :=
      fun s hs => (hGderiv s hs).continuousAt.continuousWithinAt
    have hmono : MonotoneOn G (Icc a b) := by
      apply monotoneOn_of_deriv_nonneg (convex_Icc a b) hGcont
      · intro s hs
        exact (hGderiv s (hin s hs)).differentiableAt.differentiableWithinAt
      · intro s hs
        rw [(hGderiv s (hin s hs)).deriv]
        have h0 := hs_pos s hs
        have h1 := hcap s (hin s hs)
        rw [div_eq_mul_inv]
        have hinv : 0 < s⁻¹ := inv_pos.2 h0
        nlinarith
    have := hmono (left_mem_Icc.2 hab.le) (right_mem_Icc.2 hab.le) hab.le
    simp only [hG] at this
    linarith

/-- Claim 2: every primitive in the nesting model has `δ ≤ 1`, and a feasible primitive
has `δ = 1`, so one is the upper sensitivity bound, and it is attained. -/
theorem upper_is_one {N ι : Type*} (nesting : Set N) (q : N → ℝ) (budget : ℝ)
    (constraints : ι → N → ℝ) (δ : N → ℝ)
    (hcap : ∀ ν ∈ nesting, δ ν ≤ 1)
    (ν₀ : N) (hfeas : ν₀ ∈ feasibleSet nesting q budget constraints) (hone : δ ν₀ = 1) :
    IsUpperSensitivity (feasibleSet nesting q budget constraints) δ 1 := by
  apply isUpperSensitivity_of_isGreatest
  refine ⟨⟨ν₀, hfeas, hone⟩, ?_⟩
  rintro _ ⟨ν, hν, rfl⟩
  exact hcap ν hν.1

/-- Claim 4: at a level where two marginal utilities `a, b > 0` have relative risk aversion
`ρ₁, ρ₂ ∈ (0, 10]`, their mixture with weight `w ∈ [0, 1]` has relative risk aversion
`(w a ρ₁ + (1 - w) b ρ₂) / (w a + (1 - w) b) ∈ (0, 10]`. -/
theorem mixture_rra_bounds {a b ρ₁ ρ₂ w : ℝ} (ha : 0 < a) (hb : 0 < b)
    (hw0 : 0 ≤ w) (hw1 : w ≤ 1) (h₁ : 0 < ρ₁) (h₁' : ρ₁ ≤ 10) (h₂ : 0 < ρ₂) (h₂' : ρ₂ ≤ 10) :
    0 < (w * a * ρ₁ + (1 - w) * b * ρ₂) / (w * a + (1 - w) * b) ∧
    (w * a * ρ₁ + (1 - w) * b * ρ₂) / (w * a + (1 - w) * b) ≤ 10 := by
  have hv : 0 ≤ 1 - w := by linarith
  have hden : 0 < w * a + (1 - w) * b := by
    have e1 : 0 ≤ w * (a - min a b) := mul_nonneg hw0 (by linarith [min_le_left a b])
    have e2 : 0 ≤ (1 - w) * (b - min a b) := mul_nonneg hv (by linarith [min_le_right a b])
    have hm : 0 < min a b := lt_min ha hb
    nlinarith
  have hnum : 0 < w * a * ρ₁ + (1 - w) * b * ρ₂ := by
    set m := min (a * ρ₁) (b * ρ₂)
    have hm : 0 < m := lt_min (mul_pos ha h₁) (mul_pos hb h₂)
    have e1 : 0 ≤ w * (a * ρ₁ - m) := mul_nonneg hw0 (by linarith [min_le_left (a * ρ₁) (b * ρ₂)])
    have e2 : 0 ≤ (1 - w) * (b * ρ₂ - m) :=
      mul_nonneg hv (by linarith [min_le_right (a * ρ₁) (b * ρ₂)])
    nlinarith
  refine ⟨div_pos hnum hden, ?_⟩
  rw [div_le_iff₀ hden]
  have e1 : w * a * ρ₁ ≤ w * a * 10 := mul_le_mul_of_nonneg_left h₁' (mul_nonneg hw0 ha.le)
  have e2 : (1 - w) * b * ρ₂ ≤ (1 - w) * b * 10 :=
    mul_le_mul_of_nonneg_left h₂' (mul_nonneg hv hb.le)
  nlinarith

end SW2000
