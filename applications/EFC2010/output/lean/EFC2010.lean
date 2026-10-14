import DirectedSensitivity.Framework
import Mathlib

/-!
EFC2010: the relaxed problems behind the exact bounds of the report.

`welfare γ δ D e` is the deadweight loss `∫_e^{p_C} D - (1 - δ)(p_C - e) D(e)` with efficient
price `e`, and `D * ac` is total cost, nonincreasing exactly when marginal cost is nonnegative.
`upper_isGreatest` is Claim 2 after the reduction of Claim 1; `lower_isLeast` is the
bound in Claim 3 for fixed data values.
-/

open Set intervalIntegral

namespace EFC2010

noncomputable def ac (γ δ p : ℝ) : ℝ := γ + δ * p

noncomputable def pC (γ δ : ℝ) : ℝ := γ / (1 - δ)

noncomputable def welfare (γ δ : ℝ) (D : ℝ → ℝ) (e : ℝ) : ℝ :=
  (∫ x in e..pC γ δ, D x) - (1 - δ) * (pC γ δ - e) * D e

noncomputable def unitLoss (γ δ e : ℝ) : ℝ :=
  (∫ x in e..pC γ δ, 1 / ac γ δ x) - (1 - δ) * (pC γ δ - e) / ac γ δ e

def upperSet (γ δ p₁ K : ℝ) : Set ((ℝ → ℝ) × ℝ) :=
  {ν | 0 ≤ ν.2 ∧ ν.2 ≤ pC γ δ ∧ ContinuousOn ν.1 (Icc 0 (pC γ δ)) ∧
    AntitoneOn (fun p => ν.1 p * ac γ δ p) (Icc 0 (pC γ δ)) ∧
    0 ≤ ν.1 p₁ ∧ ν.1 p₁ * ac γ δ p₁ ≤ K}

section

variable {γ δ : ℝ} (hγ : 0 < γ) (hδ : 0 < δ) (hδ1 : δ < 1)
include hγ hδ hδ1

lemma ac_pos {p : ℝ} (hp : 0 ≤ p) : 0 < ac γ δ p := by
  unfold ac; nlinarith

lemma pC_mul : (1 - δ) * pC γ δ = γ := by
  have h : (1 - δ) ≠ 0 := by linarith
  unfold pC; field_simp

lemma gamma_lt_pC : γ < pC γ δ := by
  have h := pC_mul hγ hδ hδ1
  nlinarith

lemma ac_pC : ac γ δ (pC γ δ) = pC γ δ := by
  have h := pC_mul hγ hδ hδ1
  unfold ac; linarith

lemma ac_ge {p : ℝ} (hp : 0 ≤ p) : γ ≤ ac γ δ p := by
  unfold ac; nlinarith

lemma ac_mono {a b : ℝ} (hab : a ≤ b) : ac γ δ a ≤ ac γ δ b := by
  unfold ac; nlinarith

lemma inv_cont {a b : ℝ} (ha : 0 ≤ a) : ContinuousOn (fun x => 1 / ac γ δ x) (Icc a b) := by
  apply ContinuousOn.div continuousOn_const
  · unfold ac; fun_prop
  · intro x hx; exact ne_of_gt (ac_pos hγ hδ hδ1 (ha.trans hx.1))

lemma inv_integrable {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) :
    IntervalIntegrable (fun x => 1 / ac γ δ x) MeasureTheory.volume a b :=
  (inv_cont hγ hδ hδ1 (b := b) ha).intervalIntegrable_of_Icc hab

lemma integral_inv_le {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) :
    (∫ x in a..b, 1 / ac γ δ x) ≤ (b - a) * (1 / ac γ δ a) := by
  have h := integral_mono_on hab (inv_integrable hγ hδ hδ1 ha hab)
    (intervalIntegrable_const (c := 1 / ac γ δ a)) (fun x hx => by
      apply one_div_le_one_div_of_le (ac_pos hγ hδ hδ1 ha) (ac_mono hγ hδ hδ1 hx.1))
  simpa [intervalIntegral.integral_const, smul_eq_mul] using h

lemma integral_inv_ge {a b : ℝ} (ha : 0 ≤ a) (hab : a ≤ b) :
    (b - a) * (1 / ac γ δ b) ≤ ∫ x in a..b, 1 / ac γ δ x := by
  have h := integral_mono_on hab (intervalIntegrable_const (c := 1 / ac γ δ b))
    (inv_integrable hγ hδ hδ1 ha hab) (fun x hx => by
      apply one_div_le_one_div_of_le (ac_pos hγ hδ hδ1 (ha.trans hx.1)) (ac_mono hγ hδ hδ1 hx.2))
  simpa [intervalIntegral.integral_const, smul_eq_mul] using h

/-- Demand below `M / ac` integrates below `M ∫ 1 / ac`. -/
lemma integral_le_of_total {D : ℝ → ℝ} {a b M : ℝ} (ha : 0 ≤ a) (hab : a ≤ b)
    (hD : ContinuousOn D (Icc a b)) (hM : ∀ x ∈ Icc a b, D x * ac γ δ x ≤ M) :
    (∫ x in a..b, D x) ≤ M * ∫ x in a..b, 1 / ac γ δ x := by
  rw [← intervalIntegral.integral_const_mul]
  apply integral_mono_on hab (hD.intervalIntegrable_of_Icc hab)
  · exact (continuousOn_const.mul (inv_cont hγ hδ hδ1 ha)).intervalIntegrable_of_Icc hab
  · intro x hx
    have hpos := ac_pos hγ hδ hδ1 (ha.trans hx.1)
    rw [mul_one_div, le_div_iff₀ hpos]
    exact hM x hx

lemma unitLoss_zero_nonneg : 0 ≤ unitLoss γ δ 0 := by
  have hC := gamma_lt_pC hγ hδ hδ1
  have h := integral_inv_ge hγ hδ hδ1 (le_refl 0) (by linarith : (0:ℝ) ≤ pC γ δ)
  rw [ac_pC hγ hδ hδ1] at h
  have h0 : ac γ δ 0 = γ := by unfold ac; ring
  have hm := pC_mul hγ hδ hδ1
  unfold unitLoss
  rw [h0, sub_zero, hm, div_self (ne_of_gt hγ)]
  have hne : pC γ δ ≠ 0 := by linarith
  have : (pC γ δ - 0) * (1 / pC γ δ) = 1 := by rw [sub_zero, mul_one_div, div_self hne]
  linarith

lemma unitLoss_nonneg {e : ℝ} (he : 0 ≤ e) (heC : e ≤ pC γ δ) : 0 ≤ unitLoss γ δ e := by
  have hC := gamma_lt_pC hγ hδ hδ1
  have h := integral_inv_ge hγ hδ hδ1 he heC
  rw [ac_pC hγ hδ hδ1] at h
  have hm := pC_mul hγ hδ hδ1
  have hpos := ac_pos hγ hδ hδ1 he
  have hge := ac_ge hγ hδ hδ1 he
  have hpC : 0 < pC γ δ := by linarith
  have h1 : (1 - δ) * (pC γ δ - e) / ac γ δ e ≤ (1 - δ) * (pC γ δ - e) / γ :=
    div_le_div_of_nonneg_left (by nlinarith) hγ hge
  have h3 : (1 - δ) / γ = 1 / pC γ δ := by
    rw [div_eq_div_iff (ne_of_gt hγ) (ne_of_gt hpC)]; linarith
  have key : (1 - δ) * (pC γ δ - e) / ac γ δ e ≤ (pC γ δ - e) * (1 / pC γ δ) := by
    rw [show (1 - δ) * (pC γ δ - e) / γ = (pC γ δ - e) * ((1 - δ) / γ) by ring, h3] at h1
    exact h1
  unfold unitLoss; linarith

lemma unitLoss_le_zero {e : ℝ} (he : 0 ≤ e) (heC : e ≤ pC γ δ) :
    unitLoss γ δ e ≤ unitLoss γ δ 0 := by
  have hpos := ac_pos hγ hδ hδ1 he
  have hm := pC_mul hγ hδ hδ1
  have hsplit := integral_add_adjacent_intervals
    (inv_integrable hγ hδ hδ1 (le_refl 0) he) (inv_integrable hγ hδ hδ1 he heC)
  have hlow := integral_inv_ge hγ hδ hδ1 (le_refl (0:ℝ)) he
  have h0 : ac γ δ 0 = γ := by unfold ac; ring
  unfold unitLoss
  rw [← hsplit, h0, sub_zero, hm, div_self (ne_of_gt hγ)]
  have : (e - 0) * (1 / ac γ δ e) - 1 + (1 - δ) * (pC γ δ - e) / ac γ δ e = 0 := by
    field_simp
    unfold ac
    linarith
  linarith

theorem welfare_le {p₁ K : ℝ} (hp₁ : 0 ≤ p₁) (hp₁γ : p₁ < γ) {ν : (ℝ → ℝ) × ℝ}
    (hν : ν ∈ upperSet γ δ p₁ K) :
    welfare γ δ ν.1 ν.2 ≤ K * unitLoss γ δ 0 := by
  obtain ⟨he, heC, hcont, hanti, hD1, hK⟩ := hν
  set D := ν.1
  set e := ν.2
  set H := fun p => D p * ac γ δ p with hH
  have hC := gamma_lt_pC hγ hδ hδ1
  have hp₁C : p₁ ≤ pC γ δ := by linarith
  have hpe := ac_pos hγ hδ hδ1 he
  have hHp₁ : 0 ≤ H p₁ := mul_nonneg hD1 (ac_pos hγ hδ hδ1 hp₁).le
  have hDe : (1 - δ) * (pC γ δ - e) * D e
      = H e * ((1 - δ) * (pC γ δ - e) / ac γ δ e) := by
    simp only [hH]; field_simp
  have hcontOn {a b : ℝ} (ha : 0 ≤ a) (hb : b ≤ pC γ δ) : ContinuousOn D (Icc a b) :=
    hcont.mono (Icc_subset_Icc ha hb)
  have hmid : H p₁ * unitLoss γ δ e ≤ K * unitLoss γ δ 0 :=
    calc H p₁ * unitLoss γ δ e ≤ H p₁ * unitLoss γ δ 0 :=
          mul_le_mul_of_nonneg_left (unitLoss_le_zero hγ hδ hδ1 he heC) hHp₁
      _ ≤ K * unitLoss γ δ 0 :=
          mul_le_mul_of_nonneg_right hK (unitLoss_zero_nonneg hγ hδ hδ1)
  rcases le_or_gt e p₁ with hep | hpe₁
  · have hsplit := integral_add_adjacent_intervals
      ((hcontOn he hp₁C).intervalIntegrable_of_Icc (μ := MeasureTheory.volume) hep)
      ((hcontOn hp₁ le_rfl).intervalIntegrable_of_Icc (μ := MeasureTheory.volume) hp₁C)
    have hsplit' := integral_add_adjacent_intervals
      (inv_integrable hγ hδ hδ1 he hep) (inv_integrable hγ hδ hδ1 hp₁ hp₁C)
    have hA := integral_le_of_total hγ hδ hδ1 (M := H e) he hep (hcontOn he hp₁C)
      (fun x hx => hanti ⟨he, heC⟩ ⟨he.trans hx.1, hx.2.trans hp₁C⟩ hx.1)
    have hB := integral_le_of_total hγ hδ hδ1 (M := H p₁) hp₁ hp₁C (hcontOn hp₁ le_rfl)
      (fun x hx => hanti ⟨hp₁, hp₁C⟩ ⟨hp₁.trans hx.1, hx.2⟩ hx.1)
    have hAle := integral_inv_le hγ hδ hδ1 he hep
    have hHe : H p₁ ≤ H e := hanti ⟨he, heC⟩ ⟨hp₁, hp₁C⟩ hep
    have hneg : (∫ x in e..p₁, 1 / ac γ δ x) - (1 - δ) * (pC γ δ - e) / ac γ δ e ≤ 0 := by
      have hm := pC_mul hγ hδ hδ1
      have : (p₁ - e) * (1 / ac γ δ e) ≤ (1 - δ) * (pC γ δ - e) / ac γ δ e := by
        rw [mul_one_div]
        apply div_le_div_of_nonneg_right _ hpe.le
        nlinarith
      linarith
    have hprod := mul_nonpos_of_nonneg_of_nonpos (sub_nonneg.mpr hHe) hneg
    have : welfare γ δ D e ≤ H p₁ * unitLoss γ δ e := by
      unfold welfare unitLoss
      rw [← hsplit, ← hsplit', hDe]
      nlinarith
    linarith
  · have hI := integral_le_of_total hγ hδ hδ1 (M := H e) he heC (hcontOn he le_rfl)
      (fun x hx => hanti ⟨he, heC⟩ ⟨he.trans hx.1, hx.2⟩ hx.1)
    have hHe : H e ≤ H p₁ := hanti ⟨hp₁, hp₁C⟩ ⟨he, heC⟩ hpe₁.le
    have hU := unitLoss_nonneg hγ hδ hδ1 he heC
    have : welfare γ δ D e ≤ H e * unitLoss γ δ e := by
      unfold welfare unitLoss
      rw [hDe]
      nlinarith
    nlinarith

theorem upper_isGreatest {p₁ K : ℝ} (hp₁ : 0 ≤ p₁) (hp₁γ : p₁ < γ) (hK : 0 ≤ K) :
    IsGreatest ((fun ν : (ℝ → ℝ) × ℝ => welfare γ δ ν.1 ν.2) '' upperSet γ δ p₁ K)
      (K * unitLoss γ δ 0) := by
  have hC := gamma_lt_pC hγ hδ hδ1
  refine ⟨⟨(fun p => K / ac γ δ p, 0), ?_, ?_⟩, ?_⟩
  · refine ⟨le_rfl, by linarith, ?_, ?_, ?_, ?_⟩
    · apply ContinuousOn.div continuousOn_const
      · unfold ac; fun_prop
      · intro x hx; exact ne_of_gt (ac_pos hγ hδ hδ1 hx.1)
    · intro a ha b hb _
      simp only
      rw [div_mul_cancel₀ _ (ne_of_gt (ac_pos hγ hδ hδ1 ha.1)),
        div_mul_cancel₀ _ (ne_of_gt (ac_pos hγ hδ hδ1 hb.1))]
    · exact div_nonneg hK (ac_pos hγ hδ hδ1 hp₁).le
    · simp only
      rw [div_mul_cancel₀ _ (ne_of_gt (ac_pos hγ hδ hδ1 hp₁))]
  · simp only [welfare, unitLoss]
    have : (fun x => K / ac γ δ x) = fun x => K * (1 / ac γ δ x) := by
      funext x; ring
    rw [this, intervalIntegral.integral_const_mul]
    ring
  · rintro _ ⟨ν, hν, rfl⟩
    exact welfare_le hγ hδ hδ1 hp₁ hp₁γ hν

theorem upper_sensitivity {p₁ K : ℝ} (hp₁ : 0 ≤ p₁) (hp₁γ : p₁ < γ) (hK : 0 ≤ K) :
    DirectedSensitivity.IsUpperSensitivity (upperSet γ δ p₁ K)
      (fun ν => welfare γ δ ν.1 ν.2) (K * unitLoss γ δ 0) :=
  DirectedSensitivity.isUpperSensitivity_of_isGreatest
    (upper_isGreatest hγ hδ hδ1 hp₁ hp₁γ hK)

end

noncomputable def envelope (γ δ κ p₁ d₁ h p : ℝ) : ℝ :=
  (d₁ * ((pC γ δ - p) / (pC γ δ - p₁)) ^ κ) ⊔ (h / ac γ δ p)

def lowerSet (γ δ κ p₁ d₁ h : ℝ) : Set (ℝ → ℝ) :=
  {D | ContinuousOn D (Icc p₁ (pC γ δ)) ∧ D p₁ = d₁ ∧
    ∀ p ∈ Icc p₁ (pC γ δ), envelope γ δ κ p₁ d₁ h p ≤ D p}

noncomputable def tailWelfare (γ δ p₁ : ℝ) (D : ℝ → ℝ) : ℝ :=
  (∫ x in p₁..pC γ δ, D x) - (1 - δ) * (pC γ δ - p₁) * D p₁

theorem lower_isLeast {γ δ κ p₁ d₁ h : ℝ} (hγ : 0 < γ) (hδ : 0 < δ) (hδ1 : δ < 1)
    (hp₁ : 0 ≤ p₁) (hp₁C : p₁ < pC γ δ) (hκ : 0 < κ) (hh : h / ac γ δ p₁ ≤ d₁) :
    IsLeast (tailWelfare γ δ p₁ '' lowerSet γ δ κ p₁ d₁ h)
      (tailWelfare γ δ p₁ (envelope γ δ κ p₁ d₁ h)) := by
  have hcont : ContinuousOn (envelope γ δ κ p₁ d₁ h) (Icc p₁ (pC γ δ)) := by
    apply ContinuousOn.sup
    · apply continuousOn_const.mul
      apply ContinuousOn.rpow_const
      · fun_prop
      · intro x _; exact Or.inr hκ.le
    · apply ContinuousOn.div continuousOn_const
      · unfold ac; fun_prop
      · intro x hx; exact ne_of_gt (ac_pos hγ hδ hδ1 (hp₁.trans hx.1))
  have hat : envelope γ δ κ p₁ d₁ h p₁ = d₁ := by
    unfold envelope
    rw [div_self (by linarith : pC γ δ - p₁ ≠ 0), Real.one_rpow, mul_one]
    exact sup_of_le_left hh
  refine ⟨⟨envelope γ δ κ p₁ d₁ h, ⟨hcont, hat, fun p _ => le_rfl⟩, rfl⟩, ?_⟩
  rintro _ ⟨D, ⟨hD, hD1, hge⟩, rfl⟩
  unfold tailWelfare
  rw [hat, hD1]
  have := integral_mono_on hp₁C.le
    (hcont.intervalIntegrable_of_Icc (μ := MeasureTheory.volume) hp₁C.le)
    (hD.intervalIntegrable_of_Icc (μ := MeasureTheory.volume) hp₁C.le) hge
  linarith

/-- The unit loss at zero in closed form. -/
theorem unitLoss_zero_eq {γ δ : ℝ} (hγ : 0 < γ) (hδ : 0 < δ) (hδ1 : δ < 1) :
    unitLoss γ δ 0 = -Real.log (1 - δ) / δ - 1 := by
  have hm := pC_mul hγ hδ hδ1
  have hC := gamma_lt_pC hγ hδ hδ1
  have h0 : ac γ δ 0 = γ := by unfold ac; ring
  have hint : (∫ x in (0:ℝ)..pC γ δ, 1 / ac γ δ x) = δ⁻¹ * Real.log (pC γ δ / γ) := by
    have : (fun x => 1 / ac γ δ x) = fun x => (fun y : ℝ => y⁻¹) (δ * x + γ) := by
      funext x; unfold ac; simp [add_comm]
    have e1 : δ * pC γ δ + γ = pC γ δ := by
      have := ac_pC hγ hδ hδ1; unfold ac at this; linarith
    have e2 : δ * (0:ℝ) + γ = γ := by ring
    rw [this, intervalIntegral.integral_comp_mul_add (fun y : ℝ => y⁻¹) (ne_of_gt hδ), e1, e2,
      integral_inv_of_pos hγ (by linarith), smul_eq_mul]
  unfold unitLoss
  rw [hint, h0, sub_zero, hm, div_self (ne_of_gt hγ)]
  have : pC γ δ / γ = (1 - δ)⁻¹ := by
    unfold pC; field_simp
  rw [this, Real.log_inv]
  field_simp

end EFC2010
