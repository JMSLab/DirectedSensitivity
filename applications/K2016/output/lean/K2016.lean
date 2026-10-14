import DirectedSensitivity.Framework

namespace K2016

open Finset

/-- Claim 1, one draw: at a pure-strategy entry equilibrium a participant's payoff is at least
its stay-out payoff, so its return per dollar of cost is nonnegative. -/
theorem return_nonneg {payoff stayOut cost : ℝ} (hNash : stayOut ≤ payoff) (hcost : 0 < cost) :
    0 ≤ (payoff - stayOut) / cost :=
  div_nonneg (sub_nonneg.mpr hNash) hcost.le

/-- Claim 1, aggregate: the object, a scaled average over policies of averages over matched
draws of nonnegative returns, is nonnegative (and is `0` when no policy qualifies). -/
theorem target_nonneg {ι κ : Type*} (K : Finset ι) (D : ι → Finset κ) (r : ι → κ → ℝ)
    (hr : ∀ k ∈ K, ∀ d ∈ D k, 0 ≤ r k d) :
    0 ≤ 100 / (K.card : ℝ) * ∑ k ∈ K, (∑ d ∈ D k, r k d) / ((D k).card : ℝ) := by
  apply mul_nonneg
  · positivity
  · apply sum_nonneg
    intro k hk
    exact div_nonneg (sum_nonneg fun d hd => hr k hk d hd) (by positivity)

/-- Claim 2, impact: a concave impact function with `f 0 = 0` has average product at least its
marginal product, `T f'(T) ≤ f(T)`. -/
theorem mul_deriv_le_of_concave {f : ℝ → ℝ} {f' T : ℝ} (hf : ConcaveOn ℝ (Set.Ici 0) f)
    (h0 : f 0 = 0) (hT : 0 < T) (hd : HasDerivAt f f' T) : T * f' ≤ f T := by
  have h := hf.le_slope_of_hasDerivAt (Set.mem_Ici.mpr le_rfl) (Set.mem_Ici.mpr hT.le) hT hd
  rw [slope_def_field, h0, sub_zero, sub_zero, le_div_iff₀ hT] at h
  linarith [mul_comm T f']

/-- Claim 2, contest: with `F, A ≥ 0` the enactment gain over the status quo is at least
`F` times the marginal effect `∂P/∂F`. -/
theorem gain_ge_marginal {π A F : ℝ} (hπ : π ≤ 1) (hA : 0 ≤ A) (hF : 0 ≤ F) :
    F * ((1 - π + A) / (1 + F + A) ^ 2) ≤ (π + F) / (1 + F + A) - π / (1 + A) := by
  have hA1 : 0 < 1 + A := by linarith
  have hD : 0 < 1 + F + A := by linarith
  have key : (π + F) / (1 + F + A) - π / (1 + A)
      = F * (1 - π + A) / ((1 + A) * (1 + F + A)) := by
    field_simp
    ring
  rw [key, mul_div_assoc]
  apply mul_le_mul_of_nonneg_left _ hF
  apply div_le_div_of_nonneg_left (by linarith) (by positivity)
  nlinarith

/-- Claim 2, return: with the first-order condition `s = γ v P'(F) T f'(T)`, a gain of at
least `P'(F) F`, and `T f'(T) ≤ F`, the gross return is at least `(1/γ) s/(s+c)`. -/
theorem sole_return_lower {v γ s c gain marginal F slopeEffort : ℝ} (hv : 0 ≤ v)
    (hγ : 0 < γ) (hc : 0 < c) (hs : 0 ≤ s) (hm : 0 ≤ marginal)
    (hgain : marginal * F ≤ gain) (hF : slopeEffort ≤ F)
    (hfoc : s = γ * v * marginal * slopeEffort) :
    s / γ / (s + c) ≤ gain * v / (s + c) := by
  have h1 : s / γ = v * marginal * slopeEffort := by
    rw [hfoc]
    field_simp
  have h2 : v * marginal * slopeEffort ≤ v * marginal * F :=
    mul_le_mul_of_nonneg_left hF (mul_nonneg hv hm)
  have h3 : v * (marginal * F) ≤ v * gain := mul_le_mul_of_nonneg_left hgain hv
  apply div_le_div_of_nonneg_right _ (by linarith)
  nlinarith

end K2016
