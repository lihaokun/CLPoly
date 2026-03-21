/-
  CLPoly/Algorithm/SquarefreeZp.lean — L2 SQF 算法模型与正确性证明

  Phase 3 T3.2: Squarefree Decomposition in F_p[x]
  对应 C++: polynomial_factorize_zp.hh:108-180 __squarefree_Zp

  证明结构：
  1. yunLoop: Yun 内循环（迭代 gcd 分离不同重数因子）
  2. sqfZp: 顶层 SQF（Yun + p-th root 递归）
  3. 正确性证明（TODO）
-/
import CLPoly.Spec
import CLPoly.Math.FiniteFieldFact
import Mathlib.RingTheory.UniqueFactorizationDomain.Defs
import Mathlib.Algebra.Polynomial.Expand

set_option autoImplicit false
set_option maxHeartbeats 800000

open Polynomial

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

-- ============================================================
-- 0. 辅助引理
-- ============================================================

/-- normalize a ≠ 0 当且仅当 a ≠ 0 -/
private lemma normalize_ne_zero_iff {R : Type*} [CommMonoidWithZero R] [NormalizationMonoid R]
    {a : R} : normalize a ≠ 0 ↔ a ≠ 0 := by
  rw [ne_eq, normalize_eq_zero, ne_eq]

/-- normalize 保持多项式的 natDegree -/
private lemma natDegree_normalize_eq
    (f : Polynomial (ZMod p)) : (normalize f).natDegree = f.natDegree := by
  rcases eq_or_ne f 0 with rfl | hne
  · simp
  · -- normalize f = f * C(leadingCoeff f)⁻¹ (up to unit), same degree
    have hne' : normalize f ≠ 0 := normalize_ne_zero_iff.mpr hne
    have hassoc := normalize_associated f  -- Associated f (normalize f)
    have hdeg := degree_eq_degree_of_associated hassoc
    -- hdeg : degree f = degree (normalize f)
    rw [Polynomial.degree_eq_natDegree hne, Polynomial.degree_eq_natDegree hne'] at hdeg
    exact Nat.cast_injective hdeg

/-- Monic divByMonic 的精确除法给出非零商（当被除数非零时） -/
private lemma divByMonic_ne_zero_of_ne_zero
    (f g : Polynomial (ZMod p)) (hg : Monic g) (hdvd : g ∣ f) (hf : f ≠ 0) :
    f /ₘ g ≠ 0 := by
  intro h
  have := modByMonic_add_div f hg
  rw [(modByMonic_eq_zero_iff_dvd hg).mpr hdvd, h, mul_zero, zero_add] at this
  exact hf this.symm

-- ============================================================
-- 1. Yun 内循环
-- ============================================================

/-- Yun 内循环：迭代 gcd(w, c) 分离不同重数的因子。
    前提 c ≠ 0（SQF 算法保证）。
    返回 (提取的因子列表, 残余 c)。 -/
noncomputable def yunLoop
    (w c : Polynomial (ZMod p)) (i : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    (hc : c ≠ 0)
    : List (Polynomial (ZMod p) × ℕ) × Polynomial (ZMod p) :=
  if _hw : w.natDegree = 0 then
    (acc, c)
  else
    let y := normalize (EuclideanDomain.gcd w c)
    let z := normalize (w /ₘ y)
    let acc' := if 0 < z.natDegree then acc ++ [(z, i)] else acc
    let c' := normalize (c /ₘ y)
    have hw_ne : w ≠ 0 := by intro h; simp [h] at _hw
    have hgcd_ne : EuclideanDomain.gcd w c ≠ 0 := by
      intro h; exact hw_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_left w c))
    have hy_monic : Monic y := Polynomial.monic_normalize hgcd_ne
    have hy_dvd_c : y ∣ c := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right w c)
    have hc' : c' ≠ 0 := by
      apply normalize_ne_zero_iff.mpr
      exact divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
    yunLoop y c' (i + 1) acc' hc'
termination_by w.natDegree + c.natDegree
decreasing_by
  -- New measure: deg(y) + deg(c') < deg(w) + deg(c)
  -- Key: deg(c) = deg(y) + deg(c/y), and deg(w) ≥ 1
  -- So new measure = deg(y) + deg(c') ≤ deg(y) + deg(c/y) = deg(c) < deg(w) + deg(c)
  have hw_ne : w ≠ 0 := by intro h; simp [h] at _hw
  have hgcd_ne' : EuclideanDomain.gcd w c ≠ 0 := by
    intro h; exact hw_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_left w c))
  have hy_monic' : Monic (normalize (EuclideanDomain.gcd w c)) :=
    Polynomial.monic_normalize hgcd_ne'
  have hy_ne : normalize (EuclideanDomain.gcd w c) ≠ 0 := Monic.ne_zero hy_monic'
  have hy_dvd_c' : normalize (EuclideanDomain.gcd w c) ∣ c :=
    normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right w c)
  have hy_dvd_w : normalize (EuclideanDomain.gcd w c) ∣ w :=
    normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left w c)
  -- c = y * (c /ₘ y) → deg(c) = deg(y) + deg(c /ₘ y)
  have hc_eq : normalize (EuclideanDomain.gcd w c) *
      (c /ₘ normalize (EuclideanDomain.gcd w c)) = c := by
    have := modByMonic_add_div c hy_monic'
    rw [(modByMonic_eq_zero_iff_dvd hy_monic').mpr hy_dvd_c', zero_add] at this
    exact this
  have hq_ne : c /ₘ normalize (EuclideanDomain.gcd w c) ≠ 0 :=
    divByMonic_ne_zero_of_ne_zero c _ hy_monic' hy_dvd_c' hc
  have hdeg_c := Polynomial.natDegree_mul hy_ne hq_ne
  rw [hc_eq] at hdeg_c
  -- deg(normalize(c /ₘ y)) = deg(c /ₘ y)
  have hdeg_norm : (normalize (c /ₘ normalize (EuclideanDomain.gcd w c))).natDegree =
      (c /ₘ normalize (EuclideanDomain.gcd w c)).natDegree :=
    natDegree_normalize_eq _
  -- deg(y) ≤ deg(w) (y | w)
  have hy_le := Polynomial.natDegree_le_of_dvd hy_dvd_w hw_ne
  -- Combine: new measure = deg(y) + deg(c') = deg(y) + deg(c/y) = deg(c)
  -- old = deg(w) + deg(c) ≥ 1 + deg(c) > deg(c)
  omega

/-- yunLoop 不增加 c 的 natDegree -/
private lemma yunLoop_c_natDegree_le
    (w c : Polynomial (ZMod p)) (i : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    (hc : c ≠ 0) :
    (yunLoop w c i acc hc).2.natDegree ≤ c.natDegree := by
  -- Strong induction on the termination measure
  suffices ∀ n, ∀ w c : Polynomial (ZMod p), ∀ i acc (hc : c ≠ 0),
      n = w.natDegree + c.natDegree →
      (yunLoop w c i acc hc).2.natDegree ≤ c.natDegree from
    this _ w c i acc hc rfl
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih =>
    intro w c i acc hc hn
    rw [yunLoop]
    split
    · -- base: w.natDegree = 0 → return c
      exact Nat.le_refl _
    · -- recursive case
      rename_i hw
      dsimp only
      -- Apply IH: the recursive call has smaller measure
      -- y = normalize(gcd(w, c)), c' = normalize(c /ₘ y)
      -- deg(y) + deg(c') < deg(w) + deg(c) (from yunLoop termination)
      -- IH gives: (yunLoop y c' ...).2.natDegree ≤ c'.natDegree
      -- And c'.natDegree = (c /ₘ y).natDegree ≤ c.natDegree
      set y := normalize (EuclideanDomain.gcd w c)
      set c' := normalize (c /ₘ y) with hc'_def
      -- hc' proof needed for recursive call
      have hw_ne : w ≠ 0 := by intro h; simp [h] at hw
      have hgcd_ne : EuclideanDomain.gcd w c ≠ 0 := by
        intro h; exact hw_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_left w c))
      have hy_monic : Monic y := Polynomial.monic_normalize hgcd_ne
      have hy_dvd_c : y ∣ c := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right w c)
      have hy_dvd_w : y ∣ w := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left w c)
      have hc'_ne : c' ≠ 0 := by
        apply normalize_ne_zero_iff.mpr
        exact divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
      -- deg(c') ≤ deg(c)
      have hcq_dvd : c /ₘ y ∣ c := by
        have hmod := modByMonic_add_div c hy_monic
        rw [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_c, zero_add] at hmod
        exact ⟨y, by rw [mul_comm]; exact hmod.symm⟩
      have hc'_le : c'.natDegree ≤ c.natDegree := by
        rw [hc'_def, natDegree_normalize_eq]
        exact Polynomial.natDegree_le_of_dvd hcq_dvd hc
      -- New measure < n (same argument as yunLoop termination)
      have hy_ne : y ≠ 0 := Monic.ne_zero hy_monic
      have hq_ne : c /ₘ y ≠ 0 := divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
      have hdeg_c := Polynomial.natDegree_mul hy_ne hq_ne
      have hc_eq : y * (c /ₘ y) = c := by
        have := modByMonic_add_div c hy_monic
        rwa [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_c, zero_add] at this
      rw [hc_eq] at hdeg_c
      have hdeg_norm_c := natDegree_normalize_eq (c /ₘ y)
      have hy_le_w := Polynomial.natDegree_le_of_dvd hy_dvd_w hw_ne
      have hmeas : y.natDegree + c'.natDegree < n := by
        rw [hn, hc'_def, hdeg_norm_c, ← hdeg_c]; omega
      -- Apply IH
      have ih_result := ih _ hmeas y c' (i + 1)
          (if 0 < (normalize (w /ₘ y)).natDegree then acc ++ [(normalize (w /ₘ y), i)] else acc)
          hc'_ne rfl
      exact Nat.le_trans ih_result hc'_le

/-- yunLoop 只追加 acc，不删除：输入 acc 的元素都在输出 .1 中 -/
private lemma yunLoop_acc_subset
    (w c : Polynomial (ZMod p)) (i : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    (hc : c ≠ 0) :
    ∀ pr ∈ acc, pr ∈ (yunLoop w c i acc hc).1 := by
  suffices ∀ n, ∀ w c : Polynomial (ZMod p), ∀ i acc (hc : c ≠ 0),
      n = w.natDegree + c.natDegree →
      ∀ pr ∈ acc, pr ∈ (yunLoop w c i acc hc).1 from
    this _ w c i acc hc rfl
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih =>
    intro w c i acc hc hn pr hpr
    rw [yunLoop]; split
    · exact hpr  -- base: return acc unchanged
    · rename_i hw; dsimp only
      -- Recursive: acc' = acc or acc ++ [(z, i)], both contain acc
      -- Result = yunLoop(y, c', i+1, acc', hc').1 ⊇ acc' ⊇ acc
      set y := normalize (EuclideanDomain.gcd w c)
      set z := normalize (w /ₘ y) with hz_def
      set c' := normalize (c /ₘ y) with hc'_def
      have hw_ne : w ≠ 0 := by intro h; simp [h] at hw
      have hgcd_ne : EuclideanDomain.gcd w c ≠ 0 := by
        intro h; exact hw_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_left w c))
      have hy_monic : Monic y := Polynomial.monic_normalize hgcd_ne
      have hy_dvd_c : y ∣ c := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right w c)
      have hy_dvd_w : y ∣ w := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left w c)
      have hc'_ne : c' ≠ 0 := by
        rw [hc'_def]; apply normalize_ne_zero_iff.mpr
        exact divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
      have hmeas : y.natDegree + c'.natDegree < n := by
        have hy_ne := Monic.ne_zero hy_monic
        have hq_ne := divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
        have hc_eq : y * (c /ₘ y) = c := by
          have := modByMonic_add_div c hy_monic
          rwa [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_c, zero_add] at this
        have hdeg_c := Polynomial.natDegree_mul hy_ne hq_ne
        rw [hc_eq] at hdeg_c
        have hdeg_norm_c := natDegree_normalize_eq (c /ₘ y)
        rw [hn, hc'_def, hdeg_norm_c, ← hdeg_c]
        have := Polynomial.natDegree_le_of_dvd hy_dvd_w hw_ne; omega
      -- pr ∈ acc → pr ∈ acc' (acc' ⊇ acc)
      have hpr_acc' : pr ∈ (if 0 < z.natDegree then acc ++ [(z, i)] else acc) := by
        split_ifs <;> simp [hpr]
      -- pr ∈ acc' → pr ∈ result (by IH: acc' ⊆ result)
      exact ih _ hmeas y c' (i + 1) _ hc'_ne rfl pr hpr_acc'

/-- yunLoop 保持 c ≠ 0 -/
private lemma yunLoop_c_ne_zero
    (w c : Polynomial (ZMod p)) (i : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    (hc : c ≠ 0) :
    (yunLoop w c i acc hc).2 ≠ 0 := by
  -- By strong induction, same as yunLoop_c_natDegree_le
  suffices ∀ n, ∀ w c : Polynomial (ZMod p), ∀ i acc (hc : c ≠ 0),
      n = w.natDegree + c.natDegree →
      (yunLoop w c i acc hc).2 ≠ 0 from
    this _ w c i acc hc rfl
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih =>
    intro w c i acc hc hn
    rw [yunLoop]; split
    · exact hc  -- base: return c itself
    · rename_i hw; dsimp only
      -- Recursive: same measure decrease as yunLoop_c_natDegree_le
      set y := normalize (EuclideanDomain.gcd w c)
      set c' := normalize (c /ₘ y) with hc'_def
      have hw_ne : w ≠ 0 := by intro h; simp [h] at hw
      have hgcd_ne : EuclideanDomain.gcd w c ≠ 0 := by
        intro h; exact hw_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_left w c))
      have hy_monic : Monic y := Polynomial.monic_normalize hgcd_ne
      have hy_dvd_c : y ∣ c := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right w c)
      have hy_dvd_w : y ∣ w := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left w c)
      have hc'_ne : c' ≠ 0 := by
        rw [hc'_def]; apply normalize_ne_zero_iff.mpr
        exact divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
      -- Measure decrease: same as yunLoop termination proof
      have hy_ne := Monic.ne_zero hy_monic
      have hq_ne := divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
      have hc_eq : y * (c /ₘ y) = c := by
        have := modByMonic_add_div c hy_monic
        rwa [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_c, zero_add] at this
      have hdeg_c := Polynomial.natDegree_mul hy_ne hq_ne
      rw [hc_eq] at hdeg_c
      have hmeas : y.natDegree + c'.natDegree < n := by
        -- Same pattern as yunLoop_c_natDegree_le measure decrease
        have hdeg_norm_c := natDegree_normalize_eq (c /ₘ y)
        rw [hn, hc'_def, hdeg_norm_c, ← hdeg_c]
        have := Polynomial.natDegree_le_of_dvd hy_dvd_w hw_ne; omega
      exact ih _ hmeas y c' (i + 1) _ hc'_ne rfl

/-- yunLoop 提取 w 的所有因子：若 q | w 且 q 不可约，则 q 出现在某个 acc 条目中
    (nl-proof §3.2 Step 3) -/
private lemma yunLoop_extracts_factor
    (q w c : Polynomial (ZMod p)) (i : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    (hc : c ≠ 0) (hq : Irreducible q) (hq_dvd : q ∣ w) (hY2 : Squarefree w) :
    ∃ pr ∈ (yunLoop w c i acc hc).1, q ∣ pr.1 := by
  -- Strong induction on termination measure
  suffices ∀ n, ∀ w c : Polynomial (ZMod p), ∀ i acc (hc : c ≠ 0),
      n = w.natDegree + c.natDegree →
      Squarefree w → q ∣ w →
      ∃ pr ∈ (yunLoop w c i acc hc).1, q ∣ pr.1 from
    this _ w c i acc hc rfl hY2 hq_dvd
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih_n =>
    intro w c i acc hc hn hY2 hq_dvd
    rw [yunLoop]; split
    · -- Base: w.natDegree = 0, q | w but w is unit → q is unit, contradicts Irreducible
      rename_i hw
      exfalso; exact hq.1 (isUnit_of_dvd_unit hq_dvd (by
        rw [Polynomial.eq_C_of_natDegree_eq_zero hw]
        exact Polynomial.isUnit_C.mpr (IsUnit.mk0 _ (by
          intro h; exact (Squarefree.ne_zero hY2) (by
            rw [Polynomial.eq_C_of_natDegree_eq_zero hw, h, map_zero])))))
    · -- Recursive: q | w, either q | z (extracted) or q | y (continue in recursive call)
      rename_i hw; dsimp only
      set y := normalize (EuclideanDomain.gcd w c)
      set z := normalize (w /ₘ y) with hz_def
      set c' := normalize (c /ₘ y) with hc'_def
      have hw_ne : w ≠ 0 := by intro h; simp [h] at hw
      have hgcd_ne : EuclideanDomain.gcd w c ≠ 0 := by
        intro h; exact hw_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_left w c))
      have hy_monic : Monic y := Polynomial.monic_normalize hgcd_ne
      have hy_dvd_w : y ∣ w := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left w c)
      have hy_dvd_c : y ∣ c := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right w c)
      have hc'_ne : c' ≠ 0 := by
        rw [hc'_def]; apply normalize_ne_zero_iff.mpr
        exact divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
      -- w = y * (w/y), Squarefree w → IsRelPrime y (w/y)
      -- q | w, q prime → q | y or q | (w/y)
      have hw_eq : w = y * (w /ₘ y) := by
        have := modByMonic_add_div w hy_monic
        rw [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_w, zero_add] at this; exact this.symm
      have hq_prime := hq.prime
      have hq_dvd_prod : q ∣ y * (w /ₘ y) := hw_eq ▸ hq_dvd
      rcases hq_prime.dvd_or_dvd hq_dvd_prod with hq_y | hq_wq
      · -- q | y: q persists in next w (= y), recurse
        -- Measure decrease (same as before)
        have hy_ne := Monic.ne_zero hy_monic
        have hq_ne := divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
        have hc_eq : y * (c /ₘ y) = c := by
          have := modByMonic_add_div c hy_monic
          rwa [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_c, zero_add] at this
        have hdeg_c := Polynomial.natDegree_mul hy_ne hq_ne
        rw [hc_eq] at hdeg_c
        have hdeg_norm_c := natDegree_normalize_eq (c /ₘ y)
        have hmeas : y.natDegree + c'.natDegree < n := by
          rw [hn, hc'_def, hdeg_norm_c, ← hdeg_c]
          have := Polynomial.natDegree_le_of_dvd hy_dvd_w hw_ne; omega
        have hy_sqf := Squarefree.squarefree_of_dvd hy_dvd_w hY2
        obtain ⟨pr, hpr, hqdvd⟩ := ih_n _ hmeas y c' (i + 1) _ hc'_ne rfl hy_sqf hq_y
        exact ⟨pr, hpr, hqdvd⟩
      · -- q | (w/y): q | z (since z = normalize(w/y)), z gets extracted
        have hz_dvd : q ∣ z := dvd_trans hq_wq (normalize_associated _).symm.dvd
        have hz_ne : z ≠ 0 := by
          rw [hz_def]; exact normalize_ne_zero_iff.mpr
            (divByMonic_ne_zero_of_ne_zero w y hy_monic hy_dvd_w hw_ne)
        have hz_pos : 0 < z.natDegree := by
          rcases Nat.eq_zero_or_pos z.natDegree with h | h
          · exfalso; exact hq.1 (isUnit_of_dvd_unit hz_dvd (by
              rw [Polynomial.eq_C_of_natDegree_eq_zero h]
              exact Polynomial.isUnit_C.mpr (IsUnit.mk0 _ (by
                intro habs; exact hz_ne (by
                  rw [Polynomial.eq_C_of_natDegree_eq_zero h, habs, map_zero])))))
          · exact h
        -- (z, i) ∈ acc' (since hz_pos → acc' = acc ++ [(z,i)])
        -- yunLoop_acc_subset: acc' ⊆ result.1
        have hpr_in : (z, i) ∈ (if 0 < z.natDegree then acc ++ [(z, i)] else acc) := by
          simp [hz_pos]
        exact ⟨(z, i),
          yunLoop_acc_subset y c' (i + 1) _ hc'_ne _ hpr_in,
          hz_dvd⟩

/-- yunLoop 保持不可约因子的幂次：若 q 不可约且 q ∤ w，则 q^k | c → q^k | c_rem。
    （反方向 q^k | c_rem → q^k | c 由 c_rem | c 的传递性即得。） -/
private lemma yunLoop_preserves_pow_dvd
    (q : Polynomial (ZMod p)) (hq : Irreducible q)
    (w c : Polynomial (ZMod p)) (i : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    (hc : c ≠ 0) (hq_w : ¬(q ∣ w)) (k : ℕ) (hqk : q ^ k ∣ c) :
    q ^ k ∣ (yunLoop w c i acc hc).2 := by
  suffices ∀ n, ∀ w c : Polynomial (ZMod p), ∀ i acc (hc : c ≠ 0),
      n = w.natDegree + c.natDegree → ¬(q ∣ w) → q ^ k ∣ c →
      q ^ k ∣ (yunLoop w c i acc hc).2 from
    this _ w c i acc hc rfl hq_w hqk
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih =>
    intro w c i acc hc hn hq_w hqk_c
    rw [yunLoop]; split
    · exact hqk_c  -- base: return c
    · rename_i hw; dsimp only
      set y := normalize (EuclideanDomain.gcd w c)
      set c' := normalize (c /ₘ y) with hc'_def
      -- q ∤ w → q ∤ gcd(w,c) → q ∤ y
      have hw_ne : w ≠ 0 := by intro h; simp [h] at hw
      have hgcd_ne : EuclideanDomain.gcd w c ≠ 0 := by
        intro h; exact hw_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_left w c))
      have hy_monic : Monic y := Polynomial.monic_normalize hgcd_ne
      have hy_dvd_w : y ∣ w := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left w c)
      have hy_dvd_c : y ∣ c := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right w c)
      have hq_y : ¬(q ∣ y) := fun h => hq_w (dvd_trans h hy_dvd_w)
      -- q^k | c and IsCoprime q^k y → q^k | c/y
      -- q irreducible, q ∤ y → IsCoprime q y (in PID)
      have hcop_q_y : IsCoprime q y := by
        rw [← isRelPrime_iff_isCoprime]
        exact fun d hd_q hd_y => by
          obtain ⟨e, he⟩ := hd_q
          exact (hq.isUnit_or_isUnit he).elim id fun hu =>
            -- IsUnit e, q = d * e → Associated d q → q | d → q | y → contradiction
            absurd (dvd_trans ⟨↑hu.unit⁻¹, by
              rw [he, mul_assoc, IsUnit.mul_val_inv hu, mul_one]⟩ hd_y) hq_y
      have hcop : IsCoprime (q ^ k) y := hcop_q_y.pow_left
      have hc_eq : y * (c /ₘ y) = c := by
        have := modByMonic_add_div c hy_monic
        rwa [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_c, zero_add] at this
      have hqk_yc : q ^ k ∣ y * (c /ₘ y) := ⟨hqk_c.choose, by rw [hc_eq]; exact hqk_c.choose_spec⟩
      have hqk_cq : q ^ k ∣ (c /ₘ y) := hcop.dvd_of_dvd_mul_right (by rwa [mul_comm] at hqk_yc)
      -- q^k | normalize(c/y) = c'
      have hqk_c' : q ^ k ∣ c' :=
        dvd_trans hqk_cq (normalize_associated _).symm.dvd
      -- Recursive: q ∤ y (new w)
      have hc'_ne : c' ≠ 0 := by
        rw [hc'_def]; apply normalize_ne_zero_iff.mpr
        exact divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
      -- Measure decrease (same pattern)
      have hmeas : y.natDegree + c'.natDegree < n := by
        have hy_ne := Monic.ne_zero hy_monic
        have hq_ne' := divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
        have hdeg_c := Polynomial.natDegree_mul hy_ne hq_ne'
        rw [hc_eq] at hdeg_c
        have hdeg_norm := natDegree_normalize_eq (c /ₘ y)
        rw [hn, hc'_def, hdeg_norm, ← hdeg_c]
        have := Polynomial.natDegree_le_of_dvd hy_dvd_w hw_ne; omega
      exact ih _ hmeas y c' (i + 1) _ hc'_ne rfl hq_y hqk_c'

/-- contract p 的 natDegree 不超过原多项式 natDegree 除以 p -/
private lemma natDegree_contract_le
    (f : Polynomial (ZMod p)) :
    (Polynomial.contract p f).natDegree ≤ f.natDegree / p := by
  have hp_ne : p ≠ 0 := Nat.Prime.ne_zero hp.out
  apply Polynomial.natDegree_le_iff_coeff_eq_zero.mpr
  intro n hn
  rw [Polynomial.coeff_contract hp_ne]
  apply Polynomial.coeff_eq_zero_of_natDegree_lt
  exact (Nat.div_lt_iff_lt_mul (Nat.Prime.pos hp.out)).mp hn

-- ============================================================
-- 2. 顶层 SQF
-- ============================================================

/-- SQF 顶层：Yun 循环 + char p 的 p-th root 递归 -/
noncomputable def sqfZp
    (f : Polynomial (ZMod p))
    : List (Polynomial (ZMod p) × ℕ) :=
  if _hf : f.natDegree = 0 then []
  else if _hderiv : derivative f = 0 then
    let g := Polynomial.contract p f
    (sqfZp g).map (fun pr => (pr.1, pr.2 * p))
  else
    let c := normalize (EuclideanDomain.gcd f (derivative f))
    let w := normalize (f /ₘ c)
    have hf_ne : f ≠ 0 := by intro h; simp [h] at _hf
    have hc_ne : c ≠ 0 := by
      apply normalize_ne_zero_iff.mpr
      intro h; exact hf_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_left f (derivative f)))
    let yun_output := yunLoop w c 1 [] hc_ne
    let yun_result := yun_output.1
    let c_rem := yun_output.2
    have hcrem_le : c_rem.natDegree ≤ c.natDegree :=
      yunLoop_c_natDegree_le w c 1 [] hc_ne
    if 0 < c_rem.natDegree then
      let g := Polynomial.contract p c_rem
      yun_result ++ (sqfZp g).map (fun pr => (pr.1, pr.2 * p))
    else
      yun_result
termination_by f.natDegree
decreasing_by
  · -- Branch 1: f' = 0 → deg(contract p f) < deg(f)
    have hprime := hp.out
    -- f = expand p (contract p f) when f' = 0
    have hexp := @expand_contract _ _ p _ _ _ _hderiv hprime.ne_zero
    -- natDegree f = natDegree (contract p f) * p
    have hdeg : f.natDegree = (Polynomial.contract p f).natDegree * p := by
      conv_lhs => rw [← hexp]; rw [natDegree_expand]
    -- f.natDegree ≥ 1, p ≥ 2, so contract's degree < f's degree
    have hf_pos : 0 < f.natDegree := Nat.pos_of_ne_zero _hf
    have hg_pos : 0 < (Polynomial.contract p f).natDegree := by
      by_contra h; push_neg at h; rw [Nat.le_zero.mp h, zero_mul] at hdeg
      exact _hf hdeg
    -- a < a * p when a ≥ 1, p ≥ 2
    have hp2 := hprime.two_le
    have key : (Polynomial.contract p f).natDegree * 1 <
        (Polynomial.contract p f).natDegree * p :=
      Nat.mul_lt_mul_of_pos_left (by omega) hg_pos
    simp only [Nat.mul_one] at key; omega
  · -- Branch 2: deg(contract p c_rem) < deg(f)
    -- Need to relate c_rem (from let destructuring) to yunLoop output
    -- and chain: contract(c_rem)/p ≤ c_rem/p ≤ c₀/p ≤ f/p < f
    have hf_pos := Nat.pos_of_ne_zero _hf
    have hf_ne : f ≠ 0 := by intro h; simp [h] at _hf
    -- hcrem_le from function body: c_rem.natDegree ≤ c.natDegree
    -- c.natDegree ≤ f.natDegree (gcd divides f)
    have hc_le_f : c.natDegree ≤ f.natDegree :=
      Polynomial.natDegree_le_of_dvd
        (normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left f _)) hf_ne
    have hcrem_le_f : c_rem.natDegree ≤ f.natDegree := by
      exact Nat.le_trans hcrem_le hc_le_f
    -- contract(c_rem).natDegree ≤ c_rem.natDegree / p ≤ f.natDegree / p < f.natDegree
    calc (Polynomial.contract p c_rem).natDegree
        ≤ c_rem.natDegree / p := natDegree_contract_le _
      _ ≤ f.natDegree / p := Nat.div_le_div_right hcrem_le_f
      _ < f.natDegree := Nat.div_lt_self hf_pos hp.out.one_lt

-- ============================================================
-- 3. 正确性辅助引理
-- ============================================================

/-- 核心引理：d^k ∣ f (k ≥ 1) → d^{k-1} ∣ derivative f。
    证明：f = d^k * e → f' = d^{k-1} * (C(k)*d'*e + d*e') -/
private lemma pow_dvd_derivative_of_pow_succ_dvd
    (d f : Polynomial (ZMod p)) (n : ℕ) (h : d ^ (n + 1) ∣ f) :
    d ^ n ∣ derivative f := by
  obtain ⟨e, he⟩ := h
  rw [he, derivative_mul, Polynomial.derivative_pow_succ]
  -- f' = C(n+1) * d^n * d' * e + d^{n+1} * e'
  -- = d^n * (C(n+1) * d' * e + d * e')
  -- Goal after obtain: d^n ∣ C(n+1)*d^n*d'*e + d^{n+1}*e' (already expanded)
  apply dvd_add
  · -- d^n | C(n+1) * d^n * d' * e
    exact dvd_mul_of_dvd_left (dvd_mul_of_dvd_left
      (dvd_mul_of_dvd_right (dvd_refl _) _) _) _
  · -- d^n | d^{n+1} * derivative(e)
    exact dvd_mul_of_dvd_left (pow_dvd_pow d (Nat.le_succ n)) _

/-- f / gcd(f, f') 是 squarefree（Yun 算法的数学基石）。
    证明：若 d² | w，则 d^n | f 对所有 n（通过反复 "d^k|f → d^{k-1}|f' → d^{k-1}|c → d^{k+1}|f"），
    与 deg(f) 有限矛盾。 -/
private lemma squarefree_div_gcd_derivative
    (f : Polynomial (ZMod p)) (hf : f ≠ 0) :
    let c := normalize (EuclideanDomain.gcd f (derivative f))
    let w := normalize (f /ₘ c)
    Squarefree w := by
  intro c w
  -- c | f (gcd divides first arg)
  have hc_dvd_f : c ∣ f := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left f _)
  have hc_dvd_f' : c ∣ derivative f := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right f _)
  have hc_ne : c ≠ 0 := by
    intro h; exact hf (zero_dvd_iff.mp (h ▸ hc_dvd_f))
  have hc_monic : Monic c := Polynomial.monic_normalize (by
    intro h; exact hf (zero_dvd_iff.mp ((normalize_eq_zero.mpr h) ▸
      normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left f _))))
  -- f = c * (f /ₘ c) exactly
  have hf_eq : c * (f /ₘ c) = f := by
    have := modByMonic_add_div f hc_monic
    rwa [(modByMonic_eq_zero_iff_dvd hc_monic).mpr hc_dvd_f, zero_add] at this
  -- w = normalize(f /ₘ c), so w | f (via Associated + f/c | f)
  have hw_dvd_f : w ∣ f := dvd_trans (normalize_associated _).dvd
    ⟨c, by rw [mul_comm]; exact hf_eq.symm⟩
  -- Squarefree w: suppose d * d | w, show IsUnit d
  intro d hdd_w
  -- If d = 0: d² = 0 | w → w = 0 → f = 0, contradiction
  by_cases hd : d = 0
  · exfalso; exact hf (by
      have : w = 0 := eq_zero_of_zero_dvd (hd ▸ (dvd_trans (dvd_mul_right d d) hdd_w))
      exact eq_zero_of_zero_dvd (this ▸ hw_dvd_f))
  -- d ≠ 0. Show IsUnit by contradiction: if ¬IsUnit d, then deg(d) ≥ 1 → degree argument
  by_contra hd_nu
  -- deg(d) ≥ 1 (non-zero non-unit polynomial over a field)
  have hd_deg : 1 ≤ d.natDegree := by
    rcases Nat.eq_zero_or_pos d.natDegree with h | h
    · exfalso; apply hd_nu
      rw [Polynomial.eq_C_of_natDegree_eq_zero h]
      exact Polynomial.isUnit_C.mpr (IsUnit.mk0 _ (by
        intro habs; exact hd (by rw [Polynomial.eq_C_of_natDegree_eq_zero h, habs, map_zero])))
    · exact h
  -- Key: d^n | f for all n ≥ 2 (by induction)
  have hd_dvd_f : d * d ∣ f := dvd_trans hdd_w hw_dvd_f
  -- d | gcd(f, f') because d² | f → d | f' and d | f
  have hd_dvd_c : d ∣ c := by
    have hd_dvd_f' := pow_dvd_derivative_of_pow_succ_dvd d f 1 (by rwa [pow_succ, pow_one])
    rw [pow_one] at hd_dvd_f'
    exact dvd_trans
      (EuclideanDomain.dvd_gcd (dvd_trans (dvd_mul_left d d) hd_dvd_f) hd_dvd_f')
      (normalize_associated _).symm.dvd
  -- ∀ n, d^{n+2} | f
  have hall : ∀ n, d ^ (n + 2) ∣ f := by
    intro n; induction n with
    | zero => rwa [pow_succ, pow_one]
    | succ n ih =>
      -- ih : d^{n+2} | f
      -- d^{n+1} | f' (from pow_dvd_derivative)
      have h1 := pow_dvd_derivative_of_pow_succ_dvd d f (n + 1) ih
      -- d^{n+1} | c (from d^{n+1} | f and d^{n+1} | f')
      have h2 : d ^ (n + 1) ∣ c := dvd_trans
        (EuclideanDomain.dvd_gcd (dvd_trans (pow_dvd_pow d (by omega : n + 1 ≤ n + 2)) ih) h1)
        (normalize_associated _).symm.dvd
      -- d² | (f /ₘ c) (from d² | w ~ (f/c)) and d^{n+1} | c
      have hdd_fc : d ^ 2 ∣ (f /ₘ c) := by
        rw [sq]; exact dvd_trans hdd_w (normalize_associated _).dvd
      -- d^{n+1} * d^2 | c * (f /ₘ c) = f
      rw [show n + 3 = (n + 1) + 2 by omega, pow_add, ← hf_eq]
      exact mul_dvd_mul h2 hdd_fc
  -- Now: d^n | f for all n ≥ 2 → natDegree contradiction
  -- Take n = natDegree f: d^{natDegree f + 2} | f
  -- But natDegree(d^{natDegree f + 2}) = (natDegree f + 2) * natDegree d > natDegree f
  have := hall f.natDegree
  have hdeg := Polynomial.natDegree_le_of_dvd this hf
  rw [Polynomial.natDegree_pow] at hdeg
  -- hdeg : (f.natDegree + 2) * d.natDegree ≤ f.natDegree
  -- hd_deg : 1 ≤ d.natDegree
  -- (f.natDegree + 2) * 1 ≤ (f.natDegree + 2) * d.natDegree ≤ f.natDegree
  -- But (f.natDegree + 2) * 1 = f.natDegree + 2 > f.natDegree. Contradiction.
  have : f.natDegree + 2 ≤ f.natDegree :=
    le_trans (Nat.le_mul_of_pos_right _ hd_deg) hdeg
  omega


/-- 精确幂次 + 可分 + p∤v → q^v 不整除导数。
    f = q^v * h（q ∤ h），q' ≠ 0，p ∤ v → f' = q^{v-1}*(C(v)*q'*h + q*h')，q ∤ 内因子 → q^v ∤ f' -/
private lemma not_pow_dvd_derivative_of_separable
    (q f : Polynomial (ZMod p)) (v : ℕ) (hv : 1 ≤ v)
    (hq : Irreducible q) (hq_sep : derivative q ≠ 0)
    (hp_ndvd : ¬((p : ℕ) ∣ v))
    (h : Polynomial (ZMod p)) (hf_eq : f = q ^ v * h) (hq_h : ¬(q ∣ h)) :
    ¬(q ^ v ∣ derivative f) := by
  intro hdvd_f'
  -- f' = derivative(q^v) * h + q^v * h' = C(v)*q^{v-1}*q' * h + q^v * h'
  rw [hf_eq, derivative_mul, Polynomial.derivative_pow] at hdvd_f'
  -- q^{v-1} | both terms, factor it out
  -- Second term: q^v * h' = q^{v-1} * (q * h'). q^v | this via q^{v-1} * q.
  -- If q^v | (sum), and q^v | second term, then q^v | first term.
  -- First term: C(v) * q^{v-1} * q' * h. q^v | this means q | C(v) * q' * h.
  -- q prime: q | C(v) or q | q' or q | h.
  -- C(v) is a constant C((v:ZMod p)). q has deg ≥ 1, so q ∤ C(v) (unless C(v) = 0).
  -- p ∤ v → (v : ZMod p) ≠ 0 → C(v) ≠ 0 → C(v) is unit → q ∤ C(v).
  -- q ∤ q': deg(q') < deg(q), q irreducible → q ∤ q' (if q' ≠ 0; deg argument).
  -- q ∤ h: given.
  -- All three → q ∤ C(v)*q'*h → q^v ∤ first term → q^v ∤ sum. Contradiction.
  have hv_cast_ne : (v : ZMod p) ≠ 0 := by
    intro h; exact hp_ndvd (CharP.cast_eq_zero_iff (ZMod p) p v |>.mp h)
  have hCv_ne : Polynomial.C (v : ZMod p) ≠ 0 := by
    intro h; exact hv_cast_ne (Polynomial.C_eq_zero.mp h)
  -- q ∤ q' : deg q' < deg q and q irreducible with deg ≥ 1
  have hq_ne_zero : q ≠ 0 := hq.ne_zero
  have hq_deg_pos : 0 < q.natDegree := Irreducible.natDegree_pos hq
  have hq_ndvd_q' : ¬(q ∣ derivative q) := by
    intro hdvd
    have h1 := Polynomial.natDegree_le_of_dvd hdvd hq_sep
    have h2 := Polynomial.natDegree_derivative_lt (by omega : q.natDegree ≠ 0)
    omega
  -- q ∤ C(v) * q' * h
  have hq_ndvd_prod : ¬(q ∣ Polynomial.C ((v : ZMod p)) * derivative q * h) := by
    intro hdvd
    -- q prime, q | C(v) * q' * h → q | C(v)*q' or q | h
    rcases hq.prime.dvd_or_dvd hdvd with h1 | h2
    · -- q | C(v) * q' → q | C(v) or q | q'
      rcases hq.prime.dvd_or_dvd h1 with h3 | h4
      · -- q | C(v): impossible, C(v) is unit (non-zero constant, deg(q) ≥ 1)
        exact absurd h3 (by
          intro h; have := Polynomial.natDegree_le_of_dvd h hCv_ne
          simp at this; omega)
      · exact hq_ndvd_q' h4
    · exact hq_h h2
  -- hdvd_f' : q^v | (C(v)*q^{v-1}*q' * h + q^v * h')  (after rw)
  -- q^v | second term. So q^v | first term C(v)*q^{v-1}*q'*h.
  -- q^{v-1}*q | q^{v-1} * (C(v)*q'*h) → q | C(v)*q'*h. Contradicts hq_ndvd_prod.
  have hq_ne : q ≠ 0 := hq.ne_zero
  -- Extract: q^v | C(v)*q^{v-1}*q' * h  (by subtracting q^v * h' from both sides)
  have hdvd_first : q ^ v ∣ Polynomial.C ((v : ZMod p)) * q ^ (v - 1) * derivative q * h := by
    have hdvd_second : q ^ v ∣ q ^ v * derivative h := dvd_mul_right _ _
    have : q ^ v ∣ Polynomial.C ((v : ZMod p)) * q ^ (v - 1) * derivative q * h +
        q ^ v * derivative h - q ^ v * derivative h :=
      dvd_sub hdvd_f' hdvd_second
    rwa [add_sub_cancel_right] at this
  -- q^v = q^{v-1} * q, and C(v)*q^{v-1}*q'*h = q^{v-1} * (C(v)*q'*h)
  -- So q^{v-1} * q | q^{v-1} * (C(v)*q'*h) → q | C(v)*q'*h (cancel q^{v-1} in domain)
  have hqv_eq : q ^ v = q ^ (v - 1) * q := by rw [← pow_succ, Nat.sub_add_cancel hv]
  rw [hqv_eq] at hdvd_first
  have : q ^ (v - 1) * q ∣ q ^ (v - 1) * (Polynomial.C ((v : ZMod p)) * derivative q * h) := by
    convert hdvd_first using 1; ring
  have hq_dvd := (mul_dvd_mul_iff_left (by
    exact pow_ne_zero _ hq_ne : q ^ (v - 1) ≠ 0)).mp this
  exact hq_ndvd_prod hq_dvd

/-- Yun 残余的导数为零（独立引理，nl-proof §3.2.1）。
    所有参数是普通多项式，不涉及 yunLoop 的 set-bound 变量。 -/
private lemma derivative_of_yun_remainder_eq_zero
    (f w₀ c₀ crem P : Polynomial (ZMod p))
    (hf_ne : f ≠ 0) (hcrem_ne : crem ≠ 0)
    (hw₀_sqf : Squarefree w₀)
    (hf_decomp : Associated f (w₀ * c₀))
    (hc₀_dvd_f' : c₀ ∣ derivative f)
    (hY1 : Associated f (P * crem))
    (hcop : IsCoprime P crem)
    (hextract : ∀ q : Polynomial (ZMod p), Irreducible q → q ∣ w₀ → q ∣ P)
    (hpreserve : ∀ q : Polynomial (ZMod p), Irreducible q → ¬(q ∣ w₀) →
        ∀ k, q ^ k ∣ c₀ → q ^ k ∣ crem)
    : derivative crem = 0 := by
  by_contra hderiv_ne
  -- Step (a): sqf part of crem has irred factor q
  have hsqf_w' := squarefree_div_gcd_derivative crem hcrem_ne
  have hcg_ne : EuclideanDomain.gcd crem (derivative crem) ≠ 0 := by
    intro h; exact hcrem_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_left _ _))
  have hcg_monic : Monic (normalize (EuclideanDomain.gcd crem (derivative crem))) :=
    Polynomial.monic_normalize hcg_ne
  have hcg_dvd : normalize (EuclideanDomain.gcd crem (derivative crem)) ∣ crem :=
    normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left _ _)
  have hw'_ne := normalize_ne_zero_iff.mpr
    (divByMonic_ne_zero_of_ne_zero crem _ hcg_monic hcg_dvd hcrem_ne)
  have hw'_nu : ¬IsUnit (normalize (crem /ₘ normalize
      (EuclideanDomain.gcd crem (derivative crem)))) := by
    -- w' unit → natDegree(crem/gcd) = 0 → natDegree(gcd) = natDegree(crem)
    -- → gcd ~ crem → crem | derivative(crem) → natDegree(derivative) ≥ natDegree(crem)
    -- But natDegree(derivative) ≤ natDegree(crem) - 1. So derivative = 0. Contradiction.
    intro hu
    have hdeg0 := natDegree_eq_zero_of_isUnit hu
    rw [natDegree_normalize_eq] at hdeg0
    -- crem = gcd * (crem/gcd). deg(crem/gcd) = 0 → deg(gcd) = deg(crem).
    have hgcd_dvd_deriv := normalize_dvd_iff.mpr
      (EuclideanDomain.gcd_dvd_right crem (derivative crem))
    -- gcd | derivative(crem). deg(gcd) = deg(crem). deg(derivative) ≤ deg(crem) - 1.
    -- So gcd ∤ derivative when derivative ≠ 0 (degree too small). But gcd | derivative. So derivative = 0.
    have hcrem_eq' : normalize (EuclideanDomain.gcd crem (derivative crem)) *
        (crem /ₘ normalize (EuclideanDomain.gcd crem (derivative crem))) = crem := by
      have := modByMonic_add_div crem hcg_monic
      rwa [(modByMonic_eq_zero_iff_dvd hcg_monic).mpr hcg_dvd, zero_add] at this
    have hgcd_deg : (normalize (EuclideanDomain.gcd crem (derivative crem))).natDegree =
        crem.natDegree := by
      have hdeg_mul := Polynomial.natDegree_mul (Monic.ne_zero hcg_monic)
        (divByMonic_ne_zero_of_ne_zero crem _ hcg_monic hcg_dvd hcrem_ne)
      rw [hcrem_eq'] at hdeg_mul; omega
    have hcrem_deg_pos : 0 < crem.natDegree := by
      by_contra h; push_neg at h
      have hd0 : crem.natDegree = 0 := by omega
      have := Polynomial.eq_C_of_natDegree_eq_zero hd0
      rw [this, Polynomial.derivative_C] at hderiv_ne
      exact hderiv_ne rfl
    have h1 := Polynomial.natDegree_le_of_dvd hgcd_dvd_deriv hderiv_ne
    have h2 := Polynomial.natDegree_derivative_lt (by omega : crem.natDegree ≠ 0)
    omega
  -- Step (b): ∃ irred q | w'
  obtain ⟨q, hq_irr, hq_w'⟩ := WfDvdMonoid.exists_irreducible_factor hw'_nu hw'_ne
  -- Step (c): q | crem (q | w' | crem/gcd | crem)
  have hq_crem : q ∣ crem := by
    have hcrem_eq : normalize (EuclideanDomain.gcd crem (derivative crem)) *
        (crem /ₘ normalize (EuclideanDomain.gcd crem (derivative crem))) = crem := by
      have := modByMonic_add_div crem hcg_monic
      rwa [(modByMonic_eq_zero_iff_dvd hcg_monic).mpr hcg_dvd, zero_add] at this
    exact dvd_trans hq_w' (dvd_trans (normalize_associated _).dvd
      ⟨_, hcrem_eq.symm.trans (mul_comm _ _)⟩)
  -- Step (d): q ∤ P
  have hq_notP : ¬(q ∣ P) := fun h => hq_irr.1 (hcop.isUnit_of_dvd' h hq_crem)
  -- Step (e): q | f
  have hq_f : q ∣ f := dvd_trans hq_crem (dvd_trans (dvd_mul_left crem P) hY1.symm.dvd)
  -- Step (f): q | w₀ or q | c₀
  rcases hq_irr.prime.dvd_or_dvd (dvd_trans hq_f hf_decomp.dvd) with hq_w0 | hq_c0
  · -- Case 1: q | w₀ → hextract → q | P → contradiction hq_notP
    exact hq_notP (hextract q hq_irr hq_w0)
  · -- Case 2: q | c₀, q ∤ w₀. Infinite ascent → degree contradiction.
    have hq_nw0 : ¬(q ∣ w₀) := fun h => hq_notP (hextract q hq_irr h)
    -- Infinite ascent: q | (crem/gcd) + q^n | gcd → q^{n+1} | crem → pow_dvd_deriv → repeat
    have hq_cdg : q ∣ (crem /ₘ normalize (EuclideanDomain.gcd crem (derivative crem))) :=
      dvd_trans hq_w' (normalize_associated _).dvd
    have hcrem_eq' : normalize (EuclideanDomain.gcd crem (derivative crem)) *
        (crem /ₘ normalize (EuclideanDomain.gcd crem (derivative crem))) = crem := by
      have := modByMonic_add_div crem hcg_monic
      rwa [(modByMonic_eq_zero_iff_dvd hcg_monic).mpr hcg_dvd, zero_add] at this
    have hall : ∀ n, q ^ (n + 1) ∣ crem := by
      intro n; induction n with
      | zero => rwa [pow_one]
      | succ n ih =>
        have hd := pow_dvd_derivative_of_pow_succ_dvd q crem n ih
        have hg := EuclideanDomain.dvd_gcd (dvd_trans (pow_dvd_pow q (by omega)) ih) hd
        have hng := dvd_trans hg (normalize_associated _).symm.dvd
        -- q^n | normalize(gcd), q | (crem/gcd) → q^{n+1} | gcd * (crem/gcd) = crem
        rw [show n + 2 = (n + 1) + 1 from rfl, pow_succ, ← hcrem_eq']
        exact mul_dvd_mul hng hq_cdg
    exact absurd (Polynomial.natDegree_le_of_dvd (hall crem.natDegree) hcrem_ne) (by
      rw [Polynomial.natDegree_pow]; have := Irreducible.natDegree_pos hq_irr; omega)

/-- 列表积的幂 = 各元素幂的积 -/
private lemma list_prod_pow
    (l : List (Polynomial (ZMod p))) (n : ℕ) :
    l.prod ^ n = (l.map (· ^ n)).prod := by
  induction l with
  | nil => simp
  | cons a t ih => simp [List.prod_cons, mul_pow, ih]

/-- Frobenius: expand p f = f^p in F_p[X]。
    由 map_frobenius_expand + frobenius = id on ZMod p。 -/
private lemma expand_eq_pow
    (f : Polynomial (ZMod p)) : expand (ZMod p) p f = f ^ p := by
  have hprime := hp.out
  -- map_frobenius_expand: map (frobenius (ZMod p) p) (expand (ZMod p) p f) = f ^ p
  have h := map_frobenius_expand p f
  -- frobenius (ZMod p) p = RingHom.id (ZMod p) (since a^p = a in F_p)
  have hfrob : frobenius (ZMod p) p = RingHom.id (ZMod p) := by
    ext a; simp [frobenius_def, ZMod.pow_card a]
  rw [hfrob, Polynomial.map_id] at h
  exact h

/-- SquarefreeDecomp 在 p-th root 下的提升：
    若 SquarefreeDecomp g sub，则 SquarefreeDecomp (g^p) (sub.map (fun pr => (pr.1, pr.2 * p))) -/
private lemma sqfDecomp_pow_lift
    (g : Polynomial (ZMod p)) (sub : List (Polynomial (ZMod p) × ℕ))
    (hsub : SquarefreeDecomp g sub) :
    SquarefreeDecomp (g ^ p) (sub.map (fun pr => (pr.1, pr.2 * p))) := by
  obtain ⟨hassoc, hfactors, hmult, hcoprime⟩ := hsub
  refine ⟨?_, ?_, ?_, ?_⟩
  · -- Associated (g^p) (mapped result product)
    -- g ~ ∏ sⱼ^{eⱼ} → g^p ~ (∏ sⱼ^{eⱼ})^p = ∏ sⱼ^{eⱼ*p}
    -- Step 1: g ~ prod → g^p ~ prod^p
    have h1 := hassoc.pow_pow (n := p)
    -- Step 2: prod^p = ∏ (sⱼ^eⱼ)^p = ∏ sⱼ^{eⱼ*p}
    have h2 : (sub.map (fun pr => pr.1 ^ pr.2)).prod ^ p =
        (sub.map (fun pr => pr.1 ^ (pr.2 * p))).prod := by
      rw [list_prod_pow, List.map_map]; congr 1
      exact List.map_congr_left (fun pr _ => (pow_mul pr.1 pr.2 p).symm)
    -- h1 uses sub.map (fun pr => pr.1 ^ (pr.2 * p)), goal uses sub.map (fun pr => (pr.1, pr.2 * p))
    -- Need to align: pr.1 ^ (pr.2 * p) = (fun pr => pr.1 ^ pr.2) (pr.1, pr.2 * p)
    rw [h2] at h1
    convert h1 using 1
    congr 1; simp only [List.map_map]; rfl
  · -- Each factor squarefree and monic (same factors, just different multiplicities)
    intro pr hpr
    obtain ⟨pr₀, hpr₀_mem, rfl⟩ := List.mem_map.mp hpr
    exact hfactors pr₀ hpr₀_mem
  · -- Multiplicity ≥ 1: eⱼ * p ≥ 1
    intro pr hpr
    obtain ⟨pr₀, hpr₀_mem, rfl⟩ := List.mem_map.mp hpr
    show pr₀.2 * p ≥ 1
    exact Nat.one_le_iff_ne_zero.mpr (Nat.mul_ne_zero
      (Nat.one_le_iff_ne_zero.mp (hmult pr₀ hpr₀_mem)) hp.out.ne_zero)
  · -- Pairwise coprime (same factors)
    intro pr₁ hpr₁ pr₂ hpr₂ hne
    obtain ⟨q₁, hq₁_mem, rfl⟩ := List.mem_map.mp hpr₁
    obtain ⟨q₂, hq₂_mem, rfl⟩ := List.mem_map.mp hpr₂
    have hne_orig : q₁ ≠ q₂ := fun h => hne (by rw [h])
    exact hcoprime q₁ hq₁_mem q₂ hq₂_mem hne_orig

/-- Yun 循环终止时的正确性：捕获全部不变量在 w=1 时的状态。
    f_orig = w₀ * c₀ 是原始多项式的分解。 -/
private theorem yunLoop_correct
    (w c : Polynomial (ZMod p)) (i : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    (hc : c ≠ 0)
    (f_orig : Polynomial (ZMod p))
    -- Invariants at entry:
    (hY1 : Associated f_orig ((acc.map (fun pr => pr.1 ^ pr.2)).prod * w ^ i * c))
    (hY2 : Squarefree w)
    (hY3 : i ≥ 1)
    (hY4 : ∀ pr ∈ acc, Squarefree pr.1 ∧ Monic pr.1 ∧ 0 < pr.1.natDegree ∧ pr.2 ≥ 1)
    (hY5 : ∀ pr₁ ∈ acc, ∀ pr₂ ∈ acc, pr₁ ≠ pr₂ → IsCoprime pr₁.1 pr₂.1)
    (hY6 : ∀ pr ∈ acc, IsCoprime pr.1 w)
    (hY10 : ∀ pr ∈ acc, IsCoprime pr.1 c)
    : let result := yunLoop w c i acc hc
      -- At termination, the result satisfies:
      Associated f_orig ((result.1.map (fun pr => pr.1 ^ pr.2)).prod * result.2)
      ∧ (∀ pr ∈ result.1, Squarefree pr.1 ∧ Monic pr.1 ∧ 0 < pr.1.natDegree ∧ pr.2 ≥ 1)
      ∧ (∀ pr₁ ∈ result.1, ∀ pr₂ ∈ result.1, pr₁ ≠ pr₂ → IsCoprime pr₁.1 pr₂.1)
      ∧ (∀ pr ∈ result.1, IsCoprime pr.1 result.2) := by
  -- Strong induction on termination measure (same as yunLoop_c_natDegree_le)
  suffices ∀ n, ∀ w c : Polynomial (ZMod p), ∀ i acc (hc : c ≠ 0),
      n = w.natDegree + c.natDegree →
      ∀ f_orig,
      Associated f_orig ((acc.map (fun pr => pr.1 ^ pr.2)).prod * w ^ i * c) →
      Squarefree w → i ≥ 1 →
      (∀ pr ∈ acc, Squarefree pr.1 ∧ Monic pr.1 ∧ 0 < pr.1.natDegree ∧ pr.2 ≥ 1) →
      (∀ pr₁ ∈ acc, ∀ pr₂ ∈ acc, pr₁ ≠ pr₂ → IsCoprime pr₁.1 pr₂.1) →
      (∀ pr ∈ acc, IsCoprime pr.1 w) →
      (∀ pr ∈ acc, IsCoprime pr.1 c) →
      let result := yunLoop w c i acc hc
      Associated f_orig ((result.1.map (fun pr => pr.1 ^ pr.2)).prod * result.2)
      ∧ (∀ pr ∈ result.1, Squarefree pr.1 ∧ Monic pr.1 ∧ 0 < pr.1.natDegree ∧ pr.2 ≥ 1)
      ∧ (∀ pr₁ ∈ result.1, ∀ pr₂ ∈ result.1, pr₁ ≠ pr₂ → IsCoprime pr₁.1 pr₂.1)
      ∧ (∀ pr ∈ result.1, IsCoprime pr.1 result.2) from
    this _ w c i acc hc rfl f_orig hY1 hY2 hY3 hY4 hY5 hY6 hY10
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih_n =>
    intro w c i acc hc hn f_orig hY1 hY2 hY3 hY4 hY5 hY6 hY10
    rw [yunLoop]
    split
    · -- Base: w.natDegree = 0 → return (acc, c)
      rename_i hw
      dsimp only
      refine ⟨?_, hY4, hY5, hY10⟩
      have hw_ne := Squarefree.ne_zero hY2
      have hw_unit : IsUnit w := by
        rw [Polynomial.eq_C_of_natDegree_eq_zero hw]
        exact Polynomial.isUnit_C.mpr (IsUnit.mk0 _ (by
          intro h; exact hw_ne (by rw [Polynomial.eq_C_of_natDegree_eq_zero hw, h, map_zero])))
      exact hY1.trans (by
        rw [show (acc.map (fun pr => pr.1 ^ pr.2)).prod * w ^ i * c =
          (acc.map (fun pr => pr.1 ^ pr.2)).prod * c * w ^ i by ring]
        exact ((associated_mul_isUnit_right_iff (IsUnit.pow i hw_unit)).mpr
          (Associated.refl _)).symm)
    · -- Recursive: one Yun step, then apply ih_n
      rename_i hw
      dsimp only
      -- Set names for intermediate values
      set y := normalize (EuclideanDomain.gcd w c) with hy_def
      set z := normalize (w /ₘ y) with hz_def
      set c' := normalize (c /ₘ y) with hc'_def
      set acc' := (if 0 < z.natDegree then acc ++ [(z, i)] else acc) with hacc'_def
      -- Key facts about y
      have hw_ne : w ≠ 0 := by intro h; simp [h] at hw
      have hgcd_ne : EuclideanDomain.gcd w c ≠ 0 := by
        intro h; exact hw_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_left w c))
      have hy_monic : Monic y := Polynomial.monic_normalize hgcd_ne
      have hy_dvd_w : y ∣ w := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left w c)
      have hy_dvd_c : y ∣ c := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right w c)
      have hy_ne : y ≠ 0 := Monic.ne_zero hy_monic
      -- c' ≠ 0
      have hc'_ne : c' ≠ 0 := by
        rw [hc'_def]; apply normalize_ne_zero_iff.mpr
        exact divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
      -- Measure decrease
      have hmeas : y.natDegree + c'.natDegree < n := by
        rw [hn]
        -- Same as yunLoop termination proof
        have hq_ne := divByMonic_ne_zero_of_ne_zero c y hy_monic hy_dvd_c hc
        have hc_eq : y * (c /ₘ y) = c := by
          have := modByMonic_add_div c hy_monic
          rwa [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_c, zero_add] at this
        have hdeg_c := Polynomial.natDegree_mul hy_ne hq_ne
        rw [hc_eq] at hdeg_c
        have hdeg_c' := natDegree_normalize_eq (c /ₘ y)
        rw [← hc'_def] at hdeg_c'
        rw [hdeg_c', ← hdeg_c]
        have := Polynomial.natDegree_le_of_dvd hy_dvd_w hw_ne
        omega
      -- Apply ih_n to the recursive call
      -- Need to verify all invariants for (y, c', i+1, acc')
      -- Useful facts for coprimality
      have hz_dvd_w : z ∣ w := by
        rw [hz_def, normalize_dvd_iff]
        have hmod := modByMonic_add_div w hy_monic
        rw [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_w, zero_add] at hmod
        exact ⟨y, by rw [mul_comm]; exact hmod.symm⟩
      have hcq_dvd_c : (c /ₘ y) ∣ c := by
        have hmod := modByMonic_add_div c hy_monic
        rw [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_c, zero_add] at hmod
        exact ⟨y, by rw [mul_comm]; exact hmod.symm⟩
      have hc'_dvd_c : c' ∣ c := by
        rw [hc'_def]; exact dvd_trans (normalize_associated (c /ₘ y)).dvd hcq_dvd_c
      apply ih_n _ hmeas y c' (i + 1) acc' hc'_ne rfl f_orig
      -- Y1: Associated f_orig (acc'_prod * y^(i+1) * c')
      -- From hY1: f_orig ~ P * w^i * c
      -- acc' is acc or acc ++ [(z, i)]. In either case:
      -- acc'_prod * y^(i+1) * c' ~ P * w^i * c (nl-proof §2 Y1)
      · -- Y1: f_orig ~ acc'_prod * y^(i+1) * c'
        -- From hY1: f_orig ~ P * w^i * c
        -- Key: w ~ z * y (normalize), c = y * c' (normalize + divByMonic)
        -- Split on whether z was added to acc
        -- In both cases: acc'_prod * y^(i+1) * c' ~ P * w^i * c
        -- (nl-proof §2 Y1 preservation: z^i * y^i * c = (z*y)^i * c ~ w^i * c)
        have hw_eq : y * (w /ₘ y) = w := by
          have := modByMonic_add_div w hy_monic
          rwa [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_w, zero_add] at this
        have hc_eq : y * (c /ₘ y) = c := by
          have := modByMonic_add_div c hy_monic
          rwa [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_c, zero_add] at this
        -- Associated (w /ₘ y) z and (c /ₘ y) c' (via normalize)
        have hz_assoc : Associated (w /ₘ y) z := (normalize_associated _).symm
        have hc'_assoc : Associated (c /ₘ y) c' := (normalize_associated _).symm
        -- acc'_prod * y^(i+1) * c'
        -- Case split on acc'
        rw [hacc'_def]; split_ifs with hz_pos
        · -- acc' = acc ++ [(z, i)]: acc'_prod = P * z^i
          simp only [List.map_append, List.prod_append, List.map_singleton, List.prod_singleton]
          -- Goal: P * z^i * y^(i+1) * c' ~ f_orig... need to show ~ P * w^i * c
          -- P * w^i * c ~ P * z^i * y^(i+1) * c' (split case)
          -- Strategy: rewrite w and c with exact equalities, then use Associated for normalize
          conv at hY1 => rhs; rw [show (acc.map (fun pr => pr.1 ^ pr.2)).prod * w ^ i * c =
            (acc.map (fun pr => pr.1 ^ pr.2)).prod * (y * (w /ₘ y)) ^ i * (y * (c /ₘ y)) by
            rw [hw_eq, hc_eq]]
          -- Now hY1 has exact w/y and c/y terms
          -- Need: Associated (P * (y*(w/y))^i * (y*(c/y))) (P * z^i * y^(i+1) * c')
          -- (y*(w/y))^i = y^i * (w/y)^i, Associated (w/y) z, Associated (c/y) c'
          exact hY1.trans (by
            have h1 := hz_assoc.pow_pow (n := i)  -- Associated ((w/y)^i) (z^i)
            have h2 := hc'_assoc                  -- Associated (c/y) c'
            -- P * (y*(w/y))^i * (y*(c/y))
            -- = P * y^i * (w/y)^i * y * (c/y)  (by mul_pow + ring)
            -- ~ P * y^i * z^i * y * c'          (by h1, h2)
            -- = P * z^i * y^(i+1) * c'          (by ring)
            -- Goal: Associated (P * (y*(w/y))^i * (y*(c/y))) (P * z^i * y^{i+1} * c')
            -- Rewrite: (y*(w/y))^i = y^i * (w/y)^i, then group
            have lhs_eq : (acc.map (fun pr => pr.1 ^ pr.2)).prod * (y * (w /ₘ y)) ^ i * (y * (c /ₘ y)) =
                (acc.map (fun pr => pr.1 ^ pr.2)).prod * ((w /ₘ y) ^ i * ((c /ₘ y) * y ^ (i + 1))) := by
              rw [mul_pow]; ring
            have rhs_eq : (acc.map (fun pr => pr.1 ^ pr.2)).prod * z ^ i * y ^ (i + 1) * c' =
                (acc.map (fun pr => pr.1 ^ pr.2)).prod * (z ^ i * (c' * y ^ (i + 1))) := by ring
            rw [lhs_eq, rhs_eq]
            exact (Associated.refl _).mul_mul
              (h1.mul_mul (hc'_assoc.mul_right (y ^ (i + 1)))))
        · -- acc' = acc: z has deg 0 (no split case)
          -- (w/y) is unit → w ~ y → w^i ~ y^i
          conv at hY1 => rhs; rw [show (acc.map (fun pr => pr.1 ^ pr.2)).prod * w ^ i * c =
            (acc.map (fun pr => pr.1 ^ pr.2)).prod * (y * (w /ₘ y)) ^ i * (y * (c /ₘ y)) by
            rw [hw_eq, hc_eq]]
          exact hY1.trans (by
            have lhs_eq : (acc.map (fun pr => pr.1 ^ pr.2)).prod * (y * (w /ₘ y)) ^ i * (y * (c /ₘ y)) =
                (acc.map (fun pr => pr.1 ^ pr.2)).prod * ((w /ₘ y) ^ i * ((c /ₘ y) * y ^ (i + 1))) := by
              rw [mul_pow]; ring
            have rhs_eq : (acc.map (fun pr => pr.1 ^ pr.2)).prod * y ^ (i + 1) * c' =
                (acc.map (fun pr => pr.1 ^ pr.2)).prod * (1 * (c' * y ^ (i + 1))) := by ring
            rw [lhs_eq, rhs_eq]
            -- (w/y) is unit (deg 0, non-zero) → (w/y)^i is unit → Associated (w/y)^i 1
            have hwq_ne := divByMonic_ne_zero_of_ne_zero w y hy_monic hy_dvd_w hw_ne
            have hwq_unit : IsUnit (w /ₘ y) := by
              rw [hz_def, natDegree_normalize_eq] at hz_pos
              have hwq_deg : (w /ₘ y).natDegree = 0 := by omega
              rw [Polynomial.eq_C_of_natDegree_eq_zero hwq_deg]
              exact Polynomial.isUnit_C.mpr (IsUnit.mk0 _ (by
                intro h; exact hwq_ne (by rw [Polynomial.eq_C_of_natDegree_eq_zero hwq_deg, h, map_zero])))
            have h1 : Associated ((w /ₘ y) ^ i) 1 :=
              associated_one_iff_isUnit.mpr (IsUnit.pow i hwq_unit)
            exact (Associated.refl _).mul_mul
              (h1.mul_mul (hc'_assoc.mul_right (y ^ (i + 1)))))
      -- Y2: Squarefree y
      · exact Squarefree.squarefree_of_dvd hy_dvd_w hY2
      -- Y3: i+1 ≥ 1
      · omega
      -- Y4: acc' entries correct (old entries + possibly new (z, i))
      · intro pr hpr
        rw [hacc'_def] at hpr
        split_ifs at hpr with hz_pos
        · -- z was added: pr ∈ acc ++ [(z, i)]
          simp at hpr; rcases hpr with h_old | h_new
          · exact hY4 pr h_old
          · rw [h_new]; refine ⟨?_, ?_, ?_, ?_⟩
            · -- Squarefree z: z | w + Squarefree w
              have hz_dvd : z ∣ w := by
                rw [hz_def, normalize_dvd_iff]
                have hmod := modByMonic_add_div w hy_monic
                rw [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_w, zero_add] at hmod
                exact ⟨y, by rw [mul_comm]; exact hmod.symm⟩
              exact Squarefree.squarefree_of_dvd hz_dvd hY2
            · -- Monic z: normalize of non-zero
              exact Polynomial.monic_normalize (by
                exact divByMonic_ne_zero_of_ne_zero w y hy_monic hy_dvd_w hw_ne)
            · exact hz_pos
            · exact hY3
        · exact hY4 pr hpr
      -- Y5: pairwise coprime in acc'
      · intro pr₁ hpr₁ pr₂ hpr₂ hne
        rw [hacc'_def] at hpr₁ hpr₂
        split_ifs at hpr₁ hpr₂ with hz_pos
        · -- Both in acc ++ [(z, i)]
          simp at hpr₁ hpr₂
          rcases hpr₁ with h₁_old | h₁_new <;> rcases hpr₂ with h₂_old | h₂_new
          · exact hY5 _ h₁_old _ h₂_old hne
          · rw [h₂_new]; exact (hY6 _ h₁_old).of_isCoprime_of_dvd_right hz_dvd_w
          · rw [h₁_new]; exact (hY6 _ h₂_old).of_isCoprime_of_dvd_right hz_dvd_w |>.symm
          · exact absurd (h₁_new.trans h₂_new.symm) hne
        · exact hY5 _ hpr₁ _ hpr₂ hne
      -- Y6: acc' coprime with w_new = y
      · intro pr hpr
        rw [hacc'_def] at hpr
        split_ifs at hpr with hz_pos
        · simp at hpr; rcases hpr with h_old | h_new
          · exact (hY6 _ h_old).of_isCoprime_of_dvd_right hy_dvd_w
          · -- z coprime with y: from Squarefree w, w = y * (w/y), z ~ w/y
            rw [h_new]
            -- IsCoprime z y from Squarefree w = Squarefree (y * (w/y))
            have hw_eq' : w = y * (w /ₘ y) := by
              have := modByMonic_add_div w hy_monic
              rw [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_w, zero_add] at this
              exact this.symm
            have hsqf_wy : Squarefree (y * (w /ₘ y)) := by rw [← hw_eq']; exact hY2
            have hrel := (squarefree_mul_iff.mp hsqf_wy).1  -- IsRelPrime y (w/y)
            -- IsRelPrime y (w/y) → IsRelPrime (w/y) y → IsRelPrime z y (via Associated)
            -- → IsCoprime z y
            exact (hrel.symm.of_dvd_left (by rw [hz_def]; exact (normalize_associated _).dvd)).isCoprime
        · exact (hY6 _ hpr).of_isCoprime_of_dvd_right hy_dvd_w
      -- Y10: coprime with c'
      · intro pr hpr
        rw [hacc'_def] at hpr
        split_ifs at hpr with hz_pos
        · simp at hpr; rcases hpr with h_old | h_new
          · exact (hY10 _ h_old).of_isCoprime_of_dvd_right hc'_dvd_c
          · -- z coprime with c': nl-proof §2 Y10 gcd argument
            rw [h_new]
            -- IsCoprime z c': d | z → d | w, d | c' → d | c → d | gcd(w,c) → d | y
            -- But IsCoprime z y → IsUnit d
            have hz_cop_y : IsCoprime z y := by
              have hw_eq' : w = y * (w /ₘ y) := by
                have := modByMonic_add_div w hy_monic
                rw [(modByMonic_eq_zero_iff_dvd hy_monic).mpr hy_dvd_w, zero_add] at this
                exact this.symm
              have hsqf_wy : Squarefree (y * (w /ₘ y)) := by rw [← hw_eq']; exact hY2
              have hrel := (squarefree_mul_iff.mp hsqf_wy).1
              exact (hrel.symm.of_dvd_left (by rw [hz_def]; exact (normalize_associated _).dvd)).isCoprime
            exact isRelPrime_iff_isCoprime.mp (fun d hd_z hd_c' =>
              hz_cop_y.isUnit_of_dvd' hd_z
                (dvd_trans (EuclideanDomain.dvd_gcd (dvd_trans hd_z hz_dvd_w)
                  (dvd_trans hd_c' hc'_dvd_c))
                  (normalize_associated _).symm.dvd))
        · exact (hY10 _ hpr).of_isCoprime_of_dvd_right hc'_dvd_c

-- ============================================================
-- 4. 主定理
-- ============================================================

/-- SQF 顶层正确性 -/
theorem sqf_correct
    (f : Polynomial (ZMod p)) (hf : f ≠ 0) :
    SquarefreeDecomp f (sqfZp f) := by
  -- Strong induction on f.natDegree (same measure as sqfZp)
  suffices ∀ n, ∀ f : Polynomial (ZMod p), f ≠ 0 → n = f.natDegree →
      SquarefreeDecomp f (sqfZp f) from
    this _ f hf rfl
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih =>
    intro f hf hn
    unfold sqfZp
    split
    · -- Case A: f.natDegree = 0 → result = []
      rename_i hf_deg
      refine ⟨?_, by simp, by simp, by simp⟩
      have hf_coeff : f.coeff 0 ≠ 0 := by
        intro h; exact hf (by rw [Polynomial.eq_C_of_natDegree_eq_zero hf_deg, h, map_zero])
      have hunit : IsUnit f := by
        rw [Polynomial.eq_C_of_natDegree_eq_zero hf_deg]
        exact Polynomial.isUnit_C.mpr (IsUnit.mk0 _ hf_coeff)
      simp; exact hunit
    · -- else: f.natDegree ≠ 0
      rename_i hf_deg
      split
      · -- Case B: derivative f = 0 → p-th root recursion
        rename_i hderiv
        rename_i hderiv
        have hprime := hp.out
        have hexp := @expand_contract (ZMod p) _ p _ _ f hderiv hprime.ne_zero
        -- hexp : expand p (contract p f) = f
        -- contract p f ≠ 0
        have hg_ne : Polynomial.contract p f ≠ 0 := by
          intro h; rw [h, map_zero] at hexp; exact hf hexp.symm
        -- contract p f has smaller natDegree
        have hdeg_eq : f.natDegree = (Polynomial.contract p f).natDegree * p := by
          conv_lhs => rw [← hexp]; exact natDegree_expand p (Polynomial.contract p f)
        have hg_pos : 0 < (Polynomial.contract p f).natDegree := by
          by_contra h; push_neg at h; rw [Nat.le_zero.mp h, zero_mul] at hdeg_eq
          exact hf_deg hdeg_eq
        have hg_lt : (Polynomial.contract p f).natDegree < n := by
          rw [hn, hdeg_eq]
          have key := Nat.mul_lt_mul_of_pos_left (show 1 < p from hprime.one_lt) hg_pos
          simp only [Nat.mul_one] at key; exact key
        -- IH gives SquarefreeDecomp for the contraction
        have ih_g := ih _ hg_lt _ hg_ne rfl
        -- f = expand p g = g^p (Frobenius)
        have hf_eq : f = (Polynomial.contract p f) ^ p := by
          rw [← expand_eq_pow, hexp]
        -- Decompose ih_g and build SquarefreeDecomp f result
        obtain ⟨hassoc_g, hfactors_g, hmult_g, hcoprime_g⟩ := ih_g
        refine ⟨?_, ?_, ?_, ?_⟩
        · -- Associated f (result.map (^).prod)
          -- f = g^p, use sqfDecomp_pow_lift then rewrite back
          have h := (sqfDecomp_pow_lift _ _ ⟨hassoc_g, hfactors_g, hmult_g, hcoprime_g⟩).1
          rw [← hf_eq] at h; exact h
        · -- Each factor squarefree and monic
          intro pr hpr
          obtain ⟨pr₀, hpr₀_mem, rfl⟩ := List.mem_map.mp hpr
          exact hfactors_g pr₀ hpr₀_mem
        · -- Multiplicity ≥ 1
          intro pr hpr
          obtain ⟨pr₀, hpr₀_mem, rfl⟩ := List.mem_map.mp hpr
          show pr₀.2 * p ≥ 1
          exact Nat.one_le_iff_ne_zero.mpr (Nat.mul_ne_zero
            (Nat.one_le_iff_ne_zero.mp (hmult_g pr₀ hpr₀_mem)) hprime.ne_zero)
        · -- Pairwise coprime
          intro pr₁ hpr₁ pr₂ hpr₂ hne
          obtain ⟨q₁, hq₁_mem, rfl⟩ := List.mem_map.mp hpr₁
          obtain ⟨q₂, hq₂_mem, rfl⟩ := List.mem_map.mp hpr₂
          exact hcoprime_g q₁ hq₁_mem q₂ hq₂_mem (fun h => hne (by rw [h]))
      · -- Case C: derivative f ≠ 0 → Yun + optional p-th root
        rename_i hderiv
        -- Case C: f' ≠ 0. All math lemmas proved; assemble the result.
        dsimp only
        -- Get yunLoop_correct result with initial invariants
        have hf_ne' : f ≠ 0 := by intro h; simp [h] at hf_deg
        have hsqf_w := squarefree_div_gcd_derivative f hf_ne'
        -- Y1 initial: Associated f (w * c) where w = normalize(f/c), c = normalize(gcd(f,f'))
        have hY1_init : Associated f
            ((List.map (fun (pr : Polynomial (ZMod p) × ℕ) => pr.1 ^ pr.2) []).prod *
              (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f)))) ^ 1 *
              normalize (EuclideanDomain.gcd f (derivative f))) := by
          simp only [List.map_nil, List.prod_nil, one_mul, pow_one]
          -- Associated f (w * c) from (f/c)*c = f and w ~ f/c
          have hc₀_monic : Monic (normalize (EuclideanDomain.gcd f (derivative f))) :=
            Polynomial.monic_normalize (by
              intro h; exact hf_ne' (zero_dvd_iff.mp ((normalize_eq_zero.mpr h) ▸
                normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left f _))))
          have hc₀_dvd : normalize (EuclideanDomain.gcd f (derivative f)) ∣ f :=
            normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left f _)
          have hfc_eq : normalize (EuclideanDomain.gcd f (derivative f)) *
              (f /ₘ normalize (EuclideanDomain.gcd f (derivative f))) = f := by
            have := modByMonic_add_div f hc₀_monic
            rwa [(modByMonic_eq_zero_iff_dvd hc₀_monic).mpr hc₀_dvd, zero_add] at this
          -- w * c ~ (f/c) * c = f, so Associated f (w * c)
          exact ((normalize_associated (f /ₘ normalize (EuclideanDomain.gcd f (derivative f)) :
              Polynomial (ZMod p))).mul_right
            (normalize (EuclideanDomain.gcd f (derivative f)))
            |>.trans ⟨1, by simp [mul_comm]; exact hfc_eq⟩).symm
        have hc_ne' : normalize (EuclideanDomain.gcd f (derivative f)) ≠ 0 := by
          apply normalize_ne_zero_iff.mpr
          intro h; exact hf_ne' (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_left f _))
        have hyun := yunLoop_correct
          (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f))))
          (normalize (EuclideanDomain.gcd f (derivative f)))
          1 [] hc_ne' f hY1_init hsqf_w (by omega)
          (by simp) (by simp) (by simp) (by simp)
        obtain ⟨hyun1, hyun4, hyun5, hyun10⟩ := hyun
        -- Split on c_rem.natDegree
        split
        · -- Case C.1: c_rem.natDegree > 0, result = yun_result ++ pth_root
          rename_i hcrem_pos
          -- derivative(c_rem) = 0 (nl-proof §3.2: Y10 + valuation argument)
          -- This is the one remaining mathematical claim not yet formalized
          -- Abbreviate c_rem for readability
          set crem := (yunLoop
              (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f))))
              (normalize (EuclideanDomain.gcd f (derivative f))) 1 [] hc_ne').2
            with hcrem_def
          have hcrem_ne := yunLoop_c_ne_zero
            (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f))))
            (normalize (EuclideanDomain.gcd f (derivative f))) 1 [] hc_ne'
          -- derivative(c_rem) = 0: by contradiction (nl-proof §3.2.2)
          -- If derivative ≠ 0, crem/gcd(crem,crem') has an irred factor q.
          -- q | crem, q ∤ P (IsCoprime). Either q | w₀ (→ q | P, contradiction)
          -- or q ∤ w₀ (→ emultiplicity argument → q | w₀, contradiction).
          -- Full proof requires emultiplicity API (~40 lines). See nl-proof §3.2.2.
          -- derivative(c_rem) = 0: nl-proof §3.2.1 反证法
          -- All sub-lemmas proved: squarefree_div_gcd_derivative, yunLoop_extracts_factor,
          -- yunLoop_preserves_pow_dvd, not_pow_dvd_derivative_of_separable.
          -- Proof: by_contra → get irred q from sqf part → q | w₀ (Case 1 direct,
          -- Case 2 via yunLoop_preserves_pow_dvd + not_pow_dvd_derivative contradiction)
          -- → yunLoop_extracts_factor → q | P → contradiction IsCoprime P c_rem.
          -- Blocked by: set-bound variable type matching in deeply nested context.
          -- TODO: extract as standalone lemma with clean signature.
          have hcrem_deriv : derivative crem = 0 := by sorry
          -- c_rem = (contract p c_rem)^p (Frobenius + expand_contract)
          have hcrem_ne := yunLoop_c_ne_zero
            (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f))))
            (normalize (EuclideanDomain.gcd f (derivative f))) 1 [] hc_ne'
          have hcrem_eq_pow := by
            have := @expand_contract (ZMod p) _ p _ _
              (yunLoop (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f)))) (normalize (EuclideanDomain.gcd f (derivative f))) 1 [] hc_ne').2 hcrem_deriv hp.out.ne_zero
            rw [expand_eq_pow] at this; exact this.symm
          -- Recursive IH
          have hcrem_le := yunLoop_c_natDegree_le
            (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f))))
            (normalize (EuclideanDomain.gcd f (derivative f))) 1 [] hc_ne'
          have hc_le_f : (normalize (EuclideanDomain.gcd f (derivative f))).natDegree ≤
              f.natDegree :=
            Polynomial.natDegree_le_of_dvd
              (normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_left f _)) hf_ne'
          have hcrem_le_f := Nat.le_trans hcrem_le hc_le_f
          have hg_ne : Polynomial.contract p (yunLoop (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f)))) (normalize (EuclideanDomain.gcd f (derivative f))) 1 [] hc_ne').2 ≠ 0 := by
            intro h; exact hcrem_ne (by rw [hcrem_eq_pow, h, zero_pow hp.out.ne_zero])
          have hg_lt : (Polynomial.contract p (yunLoop (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f)))) (normalize (EuclideanDomain.gcd f (derivative f))) 1 [] hc_ne').2).natDegree < n := by
            rw [hn]; calc (Polynomial.contract p _).natDegree
              ≤ (yunLoop (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f)))) (normalize (EuclideanDomain.gcd f (derivative f))) 1 [] hc_ne').2.natDegree / p := natDegree_contract_le _
              _ ≤ f.natDegree / p := Nat.div_le_div_right hcrem_le_f
              _ < f.natDegree := Nat.div_lt_self (Nat.pos_of_ne_zero hf_deg) hp.out.one_lt
          have ih_g := ih _ hg_lt _ hg_ne rfl
          -- Assemble SquarefreeDecomp for yun_result ++ pth_root
          have hpth := sqfDecomp_pow_lift _ _ ih_g
          -- c_rem ~ (contract p c_rem)^p ~ pth_prod (from hcrem_eq_pow + hpth)
          refine ⟨?_, ?_, ?_, ?_⟩
          · -- Associated f (yun_prod * pth_prod)
            simp only [List.map_append, List.prod_append]
            -- c_rem = g^p, hpth.1 : Associated (g^p) pth_prod
            -- hyun1 : Associated f (yun_prod * c_rem)
            -- Need: Associated f (yun_prod * pth_prod)
            -- Via: Associated c_rem (g^p) → Associated (g^p) pth_prod → Associated c_rem pth_prod
            -- c_rem = g^p ~ pth_prod. Chain: f ~ yun * c_rem ~ yun * pth_prod
            sorry
          · -- Squarefree + Monic
            intro pr hpr
            rcases List.mem_append.mp hpr with h | h
            · exact ⟨(hyun4 pr h).1, (hyun4 pr h).2.1⟩
            · exact hpth.2.1 pr h
          · -- Multiplicity ≥ 1
            intro pr hpr
            rcases List.mem_append.mp hpr with h | h
            · exact (hyun4 pr h).2.2.2
            · exact hpth.2.2.1 pr h
          · -- Pairwise coprime
            intro pr₁ hpr₁ pr₂ hpr₂ hne
            rcases List.mem_append.mp hpr₁ with h₁ | h₁ <;>
              rcases List.mem_append.mp hpr₂ with h₂ | h₂
            · exact hyun5 pr₁ h₁ pr₂ h₂ hne
            · -- Cross coprime: yun entry (coprime c_rem via Y10) vs pth entry (divides c_rem)
              -- pr₂ from pth: pr₂.1 | (contract p c_rem)^p = c_rem
              -- hyun10: IsCoprime pr₁.1 c_rem → IsCoprime pr₁.1 pr₂.1
              sorry
            · sorry
            · exact hpth.2.2.2 pr₁ h₁ pr₂ h₂ hne
        · -- c_rem.natDegree = 0: result = yun_result
          -- hyun1-hyun10 directly give SquarefreeDecomp (with c_rem unit from deg 0)
          exact ⟨hyun1.trans (((associated_mul_isUnit_right_iff (by
            rename_i hcrem_zero; push_neg at hcrem_zero
            have hcrem_ne := yunLoop_c_ne_zero
              (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f))))
              (normalize (EuclideanDomain.gcd f (derivative f))) 1 [] hc_ne'
            have hcrem_deg : (yunLoop
              (normalize (f /ₘ normalize (EuclideanDomain.gcd f (derivative f))))
              (normalize (EuclideanDomain.gcd f (derivative f))) 1 [] hc_ne').2.natDegree = 0 :=
              by omega
            rw [Polynomial.eq_C_of_natDegree_eq_zero hcrem_deg]
            exact Polynomial.isUnit_C.mpr (IsUnit.mk0 _ (by
              intro h; exact hcrem_ne (by
                rw [Polynomial.eq_C_of_natDegree_eq_zero hcrem_deg, h, map_zero]))))).mpr (Associated.refl _)).symm),
                  fun pr hpr => ⟨(hyun4 pr hpr).1, (hyun4 pr hpr).2.1⟩,
                  fun pr hpr => (hyun4 pr hpr).2.2.2, hyun5⟩
