import DirectedSensitivity.Framework

/-!
MMS2023: the filing hazards do not bound the present bias away from 0 or 1.

A primitive is a present bias `β` and a finite discrete law of the daily net benefit.
Days are counted back from the deadline: `cont k` is the value, `k` days before the
deadline, of a filer with no present bias, and the naive filer files that day when the
draw exceeds `β * cont k`. Every hazard sequence that is nondecreasing toward the deadline,
lies in `(0, 1)`, and meets a mean condition is reproduced exactly, under every
`β ∈ (0, 1]`, by a law with `n + 1` atoms built as a least fixed point.
-/

namespace MMS2023

open Finset DirectedSensitivity

noncomputable section

def mean (m : ℕ) (pos mass : ℕ → ℝ) : ℝ := ∑ i ∈ range m, mass i * pos i

def expMax (m : ℕ) (pos mass : ℕ → ℝ) (v : ℝ) : ℝ :=
  ∑ i ∈ range m, mass i * max (pos i) v

def cdf (m : ℕ) (pos mass : ℕ → ℝ) (t : ℝ) : ℝ :=
  ∑ i ∈ range m, if pos i ≤ t then mass i else 0

def cont (m : ℕ) (pos mass : ℕ → ℝ) : ℕ → ℝ
  | 0 => mean m pos mass
  | k + 1 => expMax m pos mass (cont m pos mass k)

/-- The filing hazard `k` days before the deadline. -/
def hazard (β : ℝ) (m : ℕ) (pos mass : ℕ → ℝ) (k : ℕ) : ℝ :=
  1 - cdf m pos mass (β * cont m pos mass k)

section Law

variable {m : ℕ} {pos pos' mass : ℕ → ℝ}

theorem wavg_le (hmass : ∀ i ∈ range m, 0 ≤ mass i) (hsum : ∑ i ∈ range m, mass i = 1)
    {f : ℕ → ℝ} {c : ℝ} (hf : ∀ i ∈ range m, f i ≤ c) :
    ∑ i ∈ range m, mass i * f i ≤ c := by
  calc ∑ i ∈ range m, mass i * f i ≤ ∑ i ∈ range m, mass i * c :=
        sum_le_sum fun i hi => mul_le_mul_of_nonneg_left (hf i hi) (hmass i hi)
    _ = c := by rw [← sum_mul, hsum, one_mul]

theorem le_wavg (hmass : ∀ i ∈ range m, 0 ≤ mass i) (hsum : ∑ i ∈ range m, mass i = 1)
    {f : ℕ → ℝ} {c : ℝ} (hf : ∀ i ∈ range m, c ≤ f i) :
    c ≤ ∑ i ∈ range m, mass i * f i := by
  calc c = ∑ i ∈ range m, mass i * c := by rw [← sum_mul, hsum, one_mul]
    _ ≤ ∑ i ∈ range m, mass i * f i :=
        sum_le_sum fun i hi => mul_le_mul_of_nonneg_left (hf i hi) (hmass i hi)

/-- One atom's share of the gap between a bound and a weighted average. -/
theorem gap_below (hmass : ∀ i ∈ range m, 0 ≤ mass i) (hsum : ∑ i ∈ range m, mass i = 1)
    {f : ℕ → ℝ} {c : ℝ} (hf : ∀ i ∈ range m, f i ≤ c) {j : ℕ} (hj : j ∈ range m) :
    mass j * (c - f j) ≤ c - ∑ i ∈ range m, mass i * f i := by
  have h : c - ∑ i ∈ range m, mass i * f i = ∑ i ∈ range m, mass i * (c - f i) := by
    simp only [mul_sub, sum_sub_distrib, ← sum_mul, hsum, one_mul]
  rw [h]
  exact single_le_sum (f := fun i => mass i * (c - f i))
    (fun i hi => mul_nonneg (hmass i hi) (sub_nonneg.2 (hf i hi))) hj

theorem gap_above (hmass : ∀ i ∈ range m, 0 ≤ mass i) (hsum : ∑ i ∈ range m, mass i = 1)
    {f : ℕ → ℝ} {c : ℝ} (hf : ∀ i ∈ range m, c ≤ f i) {j : ℕ} (hj : j ∈ range m) :
    mass j * (f j - c) ≤ ∑ i ∈ range m, mass i * f i - c := by
  have h : ∑ i ∈ range m, mass i * f i - c = ∑ i ∈ range m, mass i * (f i - c) := by
    simp only [mul_sub, sum_sub_distrib, ← sum_mul, hsum, one_mul]
  rw [h]
  exact single_le_sum (f := fun i => mass i * (f i - c))
    (fun i hi => mul_nonneg (hmass i hi) (sub_nonneg.2 (hf i hi))) hj

theorem cont_mono (hmass : ∀ i ∈ range m, 0 ≤ mass i)
    (h : ∀ i ∈ range m, pos i ≤ pos' i) : ∀ k, cont m pos mass k ≤ cont m pos' mass k
  | 0 => sum_le_sum fun i hi => mul_le_mul_of_nonneg_left (h i hi) (hmass i hi)
  | k + 1 => sum_le_sum fun i hi =>
      mul_le_mul_of_nonneg_left (max_le_max (h i hi) (cont_mono hmass h k)) (hmass i hi)

theorem cont_mem (hmass : ∀ i ∈ range m, 0 ≤ mass i) (hsum : ∑ i ∈ range m, mass i = 1)
    {a b : ℝ} (hpos : ∀ i ∈ range m, a ≤ pos i ∧ pos i ≤ b) :
    ∀ k, a ≤ cont m pos mass k ∧ cont m pos mass k ≤ b
  | 0 => ⟨le_wavg hmass hsum fun i hi => (hpos i hi).1,
          wavg_le hmass hsum fun i hi => (hpos i hi).2⟩
  | k + 1 =>
      ⟨le_wavg hmass hsum fun i hi => le_max_of_le_left (hpos i hi).1,
       wavg_le hmass hsum fun i hi => max_le (hpos i hi).2 (cont_mem hmass hsum hpos k).2⟩

/-- With mass on an atom below `b` and none above it, every continuation value is below `b`. -/
theorem cont_lt (hmass : ∀ i ∈ range m, 0 ≤ mass i) (hsum : ∑ i ∈ range m, mass i = 1)
    {b : ℝ} (hpos : ∀ i ∈ range m, pos i ≤ b) {j : ℕ} (hj : j ∈ range m)
    (hmj : 0 < mass j) (hlow : pos j < b) : ∀ k, cont m pos mass k < b
  | 0 => by
      have := gap_below hmass hsum hpos hj
      have : 0 < mass j * (b - pos j) := mul_pos hmj (by linarith)
      show ∑ i ∈ range m, mass i * pos i < b
      linarith
  | k + 1 => by
      have hk := cont_lt hmass hsum hpos hj hmj hlow k
      have := gap_below hmass hsum (f := fun i => max (pos i) (cont m pos mass k))
        (fun i hi => max_le (hpos i hi) hk.le) hj
      have : 0 < mass j * (b - max (pos j) (cont m pos mass k)) :=
        mul_pos hmj (by have := max_lt hlow hk; linarith)
      show ∑ i ∈ range m, mass i * max (pos i) (cont m pos mass k) < b
      linarith

/-- With mass on an atom above every continuation value, the values strictly increase. -/
theorem cont_lt_succ (hmass : ∀ i ∈ range m, 0 ≤ mass i)
    (hsum : ∑ i ∈ range m, mass i = 1) {j : ℕ} (hj : j ∈ range m) (hmj : 0 < mass j)
    (hhigh : ∀ k, cont m pos mass k < pos j) (k : ℕ) :
    cont m pos mass k < cont m pos mass (k + 1) := by
  have := gap_above hmass hsum (f := fun i => max (pos i) (cont m pos mass k))
    (c := cont m pos mass k) (fun i _ => le_max_right _ _) hj
  have : 0 < mass j * (max (pos j) (cont m pos mass k) - cont m pos mass k) :=
    mul_pos hmj (by rw [max_eq_left (hhigh k).le]; linarith [hhigh k])
  show cont m pos mass k < ∑ i ∈ range m, mass i * max (pos i) (cont m pos mass k)
  linarith

end Law

/-! ### The twin law -/

/-- Masses: `1 - P 0` on the low atom, `P (i - 1) - P i` on middle atom `i`, `P (n - 1)` on
the high atom. -/
def twinMass (n : ℕ) (P : ℕ → ℝ) (i : ℕ) : ℝ :=
  if i = 0 then 1 - P 0 else if i < n then P (i - 1) - P i else P (n - 1)

/-- Atoms: `a`, then middle atom `i` between the thresholds `β V (i - 1)` and `β V i`,
then `b`. -/
def twinPos (n : ℕ) (β η a b : ℝ) (V : ℕ → ℝ) (i : ℕ) : ℝ :=
  if i = 0 then a else if i < n then β * ((1 - η) * V i + η * V (i - 1)) else b

structure Hyp (n : ℕ) (β η a b : ℝ) (P : ℕ → ℝ) : Prop where
  hn : 1 ≤ n
  hβ0 : 0 < β
  hβ1 : β ≤ 1
  hη0 : 0 ≤ η
  hη1 : η < 1
  ha : a ≤ 0
  hb : 0 ≤ b
  hab : a < b
  hP0 : P 0 < 1
  hPn : 0 < P (n - 1)
  hanti : ∀ k, k + 1 < n → P (k + 1) ≤ P k

variable {n : ℕ} {β η a b : ℝ} {P : ℕ → ℝ}

theorem partial_sum (H : Hyp n β η a b P) :
    ∀ k, k < n → ∑ i ∈ range (k + 1), twinMass n P i = 1 - P k
  | 0, _ => by simp [twinMass]
  | k + 1, hk => by
      rw [sum_range_succ, partial_sum H k (by omega)]
      have h1 : k + 1 ≠ 0 := by omega
      simp only [twinMass, h1, if_false, hk, if_true, Nat.add_sub_cancel]
      ring

theorem mass_nonneg (H : Hyp n β η a b P) : ∀ i ∈ range (n + 1), 0 ≤ twinMass n P i := by
  intro i _
  unfold twinMass
  split_ifs with h0 hn
  · linarith [H.hP0]
  · have := H.hanti (i - 1) (by omega)
    rw [Nat.sub_add_cancel (by omega)] at this
    linarith
  · exact H.hPn.le

theorem mass_sum (H : Hyp n β η a b P) : ∑ i ∈ range (n + 1), twinMass n P i = 1 := by
  have hn := H.hn
  rw [sum_range_succ]
  have h := partial_sum H (n - 1) (by omega)
  rw [Nat.sub_add_cancel hn] at h
  rw [h]
  have h1 : n ≠ 0 := by omega
  simp only [twinMass, h1, if_false, lt_irrefl, sub_add_cancel]

theorem scale_mem (H : Hyp n β η a b P) {y : ℝ} (h1 : a ≤ y) (h2 : y ≤ b) :
    a ≤ β * y ∧ β * y ≤ b := by
  have := H.hβ0; have := H.hβ1; have := H.ha; have := H.hb
  rcases le_or_gt 0 y with hy | hy
  · constructor <;> nlinarith
  · constructor <;> nlinarith

theorem pos_mem (H : Hyp n β η a b P) {V : ℕ → ℝ} (hV : ∀ i, a ≤ V i ∧ V i ≤ b) :
    ∀ i, a ≤ twinPos n β η a b V i ∧ twinPos n β η a b V i ≤ b := by
  intro i
  unfold twinPos
  have := H.hη0; have := H.hη1
  split_ifs
  · exact ⟨le_rfl, H.hab.le⟩
  · obtain ⟨h1, h2⟩ := hV i
    obtain ⟨h3, h4⟩ := hV (i - 1)
    exact scale_mem H (by nlinarith) (by nlinarith)
  · exact ⟨H.hab.le, le_rfl⟩

theorem pos_mono (H : Hyp n β η a b P) {V W : ℕ → ℝ} (h : V ≤ W) :
    ∀ i, twinPos n β η a b V i ≤ twinPos n β η a b W i := by
  intro i
  unfold twinPos
  split_ifs
  · exact le_rfl
  · apply mul_le_mul_of_nonneg_left _ H.hβ0.le
    exact add_le_add (mul_le_mul_of_nonneg_left (h i) (by linarith [H.hη1]))
      (mul_le_mul_of_nonneg_left (h (i - 1)) H.hη0)
  · exact le_rfl

/-- The continuation values of the twin law built on the values `V`. -/
def update (n : ℕ) (β η a b : ℝ) (P V : ℕ → ℝ) : ℕ → ℝ :=
  cont (n + 1) (twinPos n β η a b V) (twinMass n P)

theorem update_mem (H : Hyp n β η a b P) {V : ℕ → ℝ} (hV : ∀ i, a ≤ V i ∧ V i ≤ b) :
    ∀ k, a ≤ update n β η a b P V k ∧ update n β η a b P V k ≤ b :=
  cont_mem (mass_nonneg H) (mass_sum H) fun i _ => pos_mem H hV i

theorem update_mono (H : Hyp n β η a b P) {V W : ℕ → ℝ} (h : V ≤ W) :
    update n β η a b P V ≤ update n β η a b P W :=
  cont_mono (mass_nonneg H) fun i _ => pos_mono H h i

/-- Knaster–Tarski on the box `[a, b]^ℕ`. -/
theorem exists_fixed (H : Hyp n β η a b P) :
    ∃ V : ℕ → ℝ, (∀ k, a ≤ V k ∧ V k ≤ b) ∧ update n β η a b P V = V := by
  have : Fact ((fun _ : ℕ => a) ≤ fun _ => b) := ⟨fun _ => H.hab.le⟩
  let F : Set.Icc (fun _ : ℕ => a) (fun _ : ℕ => b) →o
      Set.Icc (fun _ : ℕ => a) (fun _ : ℕ => b) :=
    { toFun := fun V => ⟨update n β η a b P V.1,
        fun k => (update_mem H (fun i => ⟨V.2.1 i, V.2.2 i⟩) k).1,
        fun k => (update_mem H (fun i => ⟨V.2.1 i, V.2.2 i⟩) k).2⟩
      monotone' := fun _ _ h => update_mono H h }
  exact ⟨F.lfp.1, fun k => ⟨F.lfp.2.1 k, F.lfp.2.2 k⟩, congrArg Subtype.val F.map_lfp⟩

/-- The twin law reproduces the hazards `P` under the present bias `β`. -/
theorem realize (H : Hyp n β η a b P) :
    ∃ pos mass : ℕ → ℝ, (∀ i ∈ range (n + 1), 0 ≤ mass i) ∧
      ∑ i ∈ range (n + 1), mass i = 1 ∧
      (∀ i ∈ range (n + 1), a ≤ pos i ∧ pos i ≤ b) ∧
      mean (n + 1) pos mass ≤ (1 - P 0) * a + P 0 * b ∧
      ∀ k < n, hazard β (n + 1) pos mass k = P k := by
  obtain ⟨V, hVmem, hfix⟩ := exists_fixed H
  set pos := twinPos n β η a b V with hposdef
  set mass := twinMass n P with hmassdef
  have hn := H.hn
  have hmass := mass_nonneg H
  have hsum := mass_sum H
  have hpos : ∀ i, a ≤ pos i ∧ pos i ≤ b := pos_mem H hVmem
  have hcont : ∀ k, cont (n + 1) pos mass k = V k := fun k => congrFun hfix k
  have h0 : (0 : ℕ) ∈ range (n + 1) := mem_range.2 (by omega)
  have hnr : n ∈ range (n + 1) := mem_range.2 (by omega)
  have hm0 : 0 < mass 0 := by simp only [hmassdef, twinMass, if_true]; linarith [H.hP0]
  have hmn : 0 < mass n := by
    have h1 : n ≠ 0 := by omega
    simp only [hmassdef, twinMass, h1, if_false, lt_irrefl]; exact H.hPn
  have hp0 : pos 0 = a := by simp [hposdef, twinPos]
  have hpn : pos n = b := by
    have h1 : n ≠ 0 := by omega
    simp only [hposdef, twinPos, h1, if_false, lt_irrefl]
  have hltb : ∀ k, V k < b := fun k => by
    rw [← hcont k]
    exact cont_lt hmass hsum (fun i _ => (hpos i).2) h0 hm0 (by rw [hp0]; exact H.hab) k
  have hstep : ∀ k, V k < V (k + 1) := fun k => by
    rw [← hcont k, ← hcont (k + 1)]
    exact cont_lt_succ hmass hsum hnr hmn (fun k => by rw [hcont k, hpn]; exact hltb k) k
  have hmono : StrictMono V := strictMono_nat_of_lt_succ hstep
  refine ⟨pos, mass, hmass, hsum, fun i _ => hpos i, ?_, ?_⟩
  · have := gap_below hmass hsum (fun i _ => (hpos i).2) h0
    rw [hp0] at this
    have hm : twinMass n P 0 = 1 - P 0 := by simp [twinMass]
    rw [hm] at this
    show ∑ i ∈ range (n + 1), twinMass n P i * pos i ≤ (1 - P 0) * a + P 0 * b
    nlinarith
  · intro k hk
    have hβ0 := H.hβ0
    have hβ1 := H.hβ1
    have key : ∀ i ∈ range (n + 1),
        (if pos i ≤ β * V k then mass i else 0) = if i ≤ k then mass i else 0 := by
      intro i hi
      have hi' := mem_range.1 hi
      by_cases hik : i ≤ k
      · rw [if_pos hik, if_pos]
        by_cases hi0 : i = 0
        · rw [hi0, hp0]
          have := (hVmem k).1
          nlinarith [H.ha]
        · have hin : i < n := by omega
          simp only [hposdef, twinPos, hi0, if_false, hin, if_true]
          apply mul_le_mul_of_nonneg_left _ hβ0.le
          have h1 : V (i - 1) ≤ V i := hmono.monotone (by omega)
          have h2 : V i ≤ V k := hmono.monotone hik
          nlinarith [H.hη0]
      · rw [if_neg hik, if_neg]
        by_cases hin : i < n
        · have hi0 : i ≠ 0 := by omega
          simp only [hposdef, twinPos, hi0, if_false, hin, if_true, not_le]
          have h1 : V (i - 1) < V i := hmono (by omega)
          have h2 : V k ≤ V (i - 1) := hmono.monotone (by omega)
          have h3 : V (i - 1) < (1 - η) * V i + η * V (i - 1) := by nlinarith [H.hη1]
          exact mul_lt_mul_of_pos_left (by linarith) hβ0
        · have hin' : i = n := by omega
          rw [hin', hpn, not_le]
          have := hltb k
          rcases le_or_gt 0 (V k) with hv | hv
          · nlinarith
          · nlinarith [H.hb]
    unfold hazard cdf
    rw [hcont k, sum_congr rfl key, ← sum_filter]
    have hfilter : (range (n + 1)).filter (· ≤ k) = range (k + 1) := by
      ext i; simp only [mem_filter, mem_range]; omega
    rw [hfilter, hmassdef, partial_sum H k hk]
    ring

/-! ### Sensitivity of the present bias -/

structure Primitive where
  β : ℝ
  m : ℕ
  pos : ℕ → ℝ
  mass : ℕ → ℝ

/-- The discrete laws of the nesting model: support in `(-25, 25]` and a nonpositive mean. -/
def nesting : Set Primitive :=
  {ν | 0 < ν.β ∧ ν.β ≤ 1 ∧ (∀ i ∈ range ν.m, 0 ≤ ν.mass i) ∧
    ∑ i ∈ range ν.m, ν.mass i = 1 ∧ (∀ i ∈ range ν.m, -25 < ν.pos i ∧ ν.pos i ≤ 25) ∧
    mean ν.m ν.pos ν.mass ≤ 0}

def hazards (n : ℕ) (ν : Primitive) : Fin n → ℝ := fun k => hazard ν.β ν.m ν.pos ν.mass k

def presentBias (ν : Primitive) : ℝ := ν.β

/-- Feasible primitives for a criterion `Q` that reads the hazards only. -/
def feasible (n : ℕ) (Q : (Fin n → ℝ) → ℝ) (budget : ℝ) : Set Primitive :=
  feasibleSet nesting (fun ν => Q (hazards n ν)) budget (Empty.elim : Empty → Primitive → ℝ)

theorem exists_feasible (n : ℕ) (hn : 1 ≤ n) (Q : (Fin n → ℝ) → ℝ) (budget : ℝ)
    (P : ℕ → ℝ) (hP0 : P 0 < 1) (hPn : 0 < P (n - 1))
    (hanti : ∀ k, k + 1 < n → P (k + 1) ≤ P k) (a b : ℝ) (hlo : -25 < a) (ha : a ≤ 0)
    (hb : 0 ≤ b) (hhi : b ≤ 25) (hab : a < b) (hmean : (1 - P 0) * a + P 0 * b ≤ 0)
    (hQ : Q (fun k => P k) ≤ budget) (β : ℝ) (hβ0 : 0 < β) (hβ1 : β ≤ 1) :
    ∃ ν ∈ feasible n Q budget, presentBias ν = β := by
  have H : Hyp n β (1 / 2) a b P :=
    ⟨hn, hβ0, hβ1, by norm_num, by norm_num, ha, hb, hab, hP0, hPn, hanti⟩
  obtain ⟨pos, mass, hmass, hsum, hpos, hmeanle, hhaz⟩ := realize H
  refine ⟨⟨β, n + 1, pos, mass⟩, ⟨⟨hβ0, hβ1, hmass, hsum,
    fun i hi => ⟨by linarith [(hpos i hi).1], by linarith [(hpos i hi).2]⟩,
    le_trans hmeanle hmean⟩, ?_, fun e => e.elim⟩, rfl⟩
  have : hazards n ⟨β, n + 1, pos, mass⟩ = fun k : Fin n => P k :=
    funext fun k : Fin n => hhaz k k.isLt
  simp only [this]
  exact hQ

/-- The present bias ranges over all of `(0, 1]`: its supremum is 1, attained, and its
infimum is 0, approached but not attained. -/
theorem present_bias_sensitivity (n : ℕ) (hn : 1 ≤ n) (Q : (Fin n → ℝ) → ℝ) (budget : ℝ)
    (P : ℕ → ℝ) (hP0 : P 0 < 1) (hPn : 0 < P (n - 1))
    (hanti : ∀ k, k + 1 < n → P (k + 1) ≤ P k) (a b : ℝ) (hlo : -25 < a) (ha : a ≤ 0)
    (hb : 0 ≤ b) (hhi : b ≤ 25) (hab : a < b) (hmean : (1 - P 0) * a + P 0 * b ≤ 0)
    (hQ : Q (fun k => P k) ≤ budget) :
    IsUpperSensitivity (feasible n Q budget) presentBias 1 ∧
      IsLowerSensitivity (feasible n Q budget) presentBias 0 ∧
      0 ∉ presentBias '' feasible n Q budget := by
  have exists_at := exists_feasible n hn Q budget P hP0 hPn hanti a b hlo ha hb hhi hab
    hmean hQ
  refine ⟨isUpperSensitivity_of_isGreatest ⟨exists_at 1 one_pos le_rfl, ?_⟩,
    isLowerSensitivity_of_ge_of_approx (fun ν hν => hν.1.1.le) ?_, ?_⟩
  · rintro _ ⟨ν, hν, rfl⟩
    exact hν.1.2.1
  · intro ε hε
    obtain ⟨ν, hν, hβ⟩ := exists_at (min (ε / 2) 1) (lt_min (by linarith) one_pos)
      (min_le_right _ _)
    refine ⟨ν, hν, ?_⟩
    rw [hβ]
    linarith [min_le_left (ε / 2) 1]
  · rintro ⟨ν, hν, h⟩
    exact absurd h (ne_of_gt hν.1.1)

end

end MMS2023
