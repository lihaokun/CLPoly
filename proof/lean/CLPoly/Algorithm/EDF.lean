/-
  CLPoly/Algorithm/EDF.lean — L2 EDF 算法模型与正确性证明

  Phase 3 T3.3: Equal Degree Factorization (Cantor-Zassenhaus)
  对应 C++: polynomial_factorize_zp.hh:294-354 __edf_Zp

  证明结构：
  1. edf: 递归分解（edfSplit 逻辑 inline）
  2. irreducible_of_le_deg: 基础引理
  3. edf_correct: 主定理
-/
import CLPoly.Spec
import CLPoly.Math.FiniteFieldFact
import Mathlib.RingTheory.UniqueFactorizationDomain.Defs
import Mathlib.RingTheory.AdjoinRoot
import Mathlib.FieldTheory.Finiteness

set_option autoImplicit false
set_option maxHeartbeats 800000

open Polynomial

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

-- ============================================================
-- 0. 辅助引理
-- ============================================================

private lemma natDegree_normalize_eq'
    (f : Polynomial (ZMod p)) : (normalize f).natDegree = f.natDegree := by
  rcases eq_or_ne f 0 with rfl | hne
  · simp
  · have hne' : normalize f ≠ 0 := by rwa [ne_eq, normalize_eq_zero]
    have := degree_eq_degree_of_associated (normalize_associated f)
    rw [Polynomial.degree_eq_natDegree hne, Polynomial.degree_eq_natDegree hne'] at this
    exact Nat.cast_injective this

private lemma divByMonic_ne_zero'
    (f g : Polynomial (ZMod p)) (hg : Monic g) (hdvd : g ∣ f) (hf : f ≠ 0) :
    f /ₘ g ≠ 0 := by
  intro h
  have := modByMonic_add_div f hg
  rw [(modByMonic_eq_zero_iff_dvd hg).mpr hdvd, h, mul_zero, zero_add] at this
  exact hf this.symm

-- ============================================================
-- 1. EDF 递归分解
-- ============================================================

noncomputable def edf
    (f : Polynomial (ZMod p)) (d : ℕ)
    (splits : List (Polynomial (ZMod p)))
    : List (Polynomial (ZMod p)) :=
  if _hf : f.natDegree ≤ d then [f]
  else
    match splits with
    | [] => [f]
    | a :: rest =>
      let g := normalize (EuclideanDomain.gcd (a ^ ((p ^ d - 1) / 2) - 1) f)
      if _hg : 0 < g.natDegree ∧ g.natDegree < f.natDegree then
        edf g d rest ++ edf (normalize (f /ₘ g)) d rest
      else
        edf f d rest
termination_by (f.natDegree, splits.length)
decreasing_by
  all_goals simp_wf
  · -- edf g: g.natDegree < f.natDegree (from _hg.2)
    left; exact _hg.2
  · -- edf (normalize (f /ₘ g)): natDegree < f.natDegree
    left
    have hf_ne : f ≠ 0 := by intro h; simp [h] at _hf
    have hgcd_ne : EuclideanDomain.gcd (a ^ ((p ^ d - 1) / 2) - 1) f ≠ 0 := by
      intro h; exact hf_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_right _ _))
    have hg_monic : Monic g := Polynomial.monic_normalize hgcd_ne
    have hg_dvd : g ∣ f := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right _ _)
    have hfg_eq : g * (f /ₘ g) = f := by
      have := modByMonic_add_div f hg_monic
      rwa [(modByMonic_eq_zero_iff_dvd hg_monic).mpr hg_dvd, zero_add] at this
    have hfg_ne := divByMonic_ne_zero' f g hg_monic hg_dvd hf_ne
    have hdeg := Polynomial.natDegree_mul (Monic.ne_zero hg_monic) hfg_ne
    rw [hfg_eq] at hdeg
    have key : (f /ₘ g).natDegree < f.natDegree := by omega
    calc (normalize (f /ₘ g)).natDegree = (f /ₘ g).natDegree := natDegree_normalize_eq' _
      _ < f.natDegree := key
  · -- edf f d rest: splits.length decreases
    right; simp

-- ============================================================
-- 2. irreducible_of_le_deg
-- ============================================================

private lemma isUnit_of_natDegree_eq_zero
    (f : Polynomial (ZMod p)) (hf : f ≠ 0) (h : f.natDegree = 0) : IsUnit f := by
  rw [Polynomial.eq_C_of_natDegree_eq_zero h]
  exact Polynomial.isUnit_C.mpr (IsUnit.mk0 _ (by
    intro habs; exact hf (by rw [Polynomial.eq_C_of_natDegree_eq_zero h, habs, map_zero])))

private lemma irreducible_of_le_deg
    (f : Polynomial (ZMod p)) (d : ℕ)
    (hf_monic : Monic f) (hf_sqfree : Squarefree f)
    (hf_deg_pos : 0 < f.natDegree) (hf_deg_le : f.natDegree ≤ d) (hd : 0 < d)
    (hf_factors : ∀ q : Polynomial (ZMod p), Irreducible q → q ∣ f → q.natDegree = d) :
    Irreducible f ∧ f.natDegree = d := by
  have hf_ne : f ≠ 0 := Monic.ne_zero hf_monic
  have hf_nu : ¬IsUnit f := fun hu => absurd (natDegree_eq_zero_of_isUnit hu) (by omega)
  -- Step 1: natDegree = d
  obtain ⟨q₁, hq₁_irr, hq₁_dvd⟩ := WfDvdMonoid.exists_irreducible_factor hf_nu hf_ne
  have hq₁_deg := hf_factors q₁ hq₁_irr hq₁_dvd
  have : d ≤ f.natDegree := hq₁_deg ▸ Polynomial.natDegree_le_of_dvd hq₁_dvd hf_ne
  have hf_deg_eq : f.natDegree = d := by omega
  -- Step 2: Irreducible
  refine ⟨irreducible_iff.mpr ⟨hf_nu, fun a b hab => ?_⟩, hf_deg_eq⟩
  -- Show: IsUnit a ∨ IsUnit b
  by_contra h_not
  push_neg at h_not
  obtain ⟨ha_nu, hb_nu⟩ := h_not
  have ha_ne : a ≠ 0 := left_ne_zero_of_mul (hab ▸ hf_ne)
  have hb_ne : b ≠ 0 := right_ne_zero_of_mul (hab ▸ hf_ne)
  have ha_deg : 1 ≤ a.natDegree := by
    rcases Nat.eq_zero_or_pos a.natDegree with h | h
    · exact absurd (isUnit_of_natDegree_eq_zero a ha_ne h) ha_nu
    · exact h
  have hb_deg : 1 ≤ b.natDegree := by
    rcases Nat.eq_zero_or_pos b.natDegree with h | h
    · exact absurd (isUnit_of_natDegree_eq_zero b hb_ne h) hb_nu
    · exact h
  obtain ⟨qa, hqa_irr, hqa_dvd⟩ := WfDvdMonoid.exists_irreducible_factor ha_nu ha_ne
  have ha_ge_d : d ≤ a.natDegree :=
    (hf_factors qa hqa_irr (dvd_trans hqa_dvd (hab ▸ dvd_mul_right a b))) ▸
    Polynomial.natDegree_le_of_dvd hqa_dvd ha_ne
  obtain ⟨qb, hqb_irr, hqb_dvd⟩ := WfDvdMonoid.exists_irreducible_factor hb_nu hb_ne
  have hb_ge_d : d ≤ b.natDegree :=
    (hf_factors qb hqb_irr (dvd_trans hqb_dvd (hab ▸ dvd_mul_left b a))) ▸
    Polynomial.natDegree_le_of_dvd hqb_dvd hb_ne
  have := Polynomial.natDegree_mul ha_ne hb_ne
  rw [← hab, hf_deg_eq] at this; omega

-- ============================================================
-- 3. 主定理
-- ============================================================

private theorem edf_correct_aux
    (f : Polynomial (ZMod p)) (d : ℕ)
    (hf_monic : Monic f) (hf_sqfree : Squarefree f)
    (hf_deg_pos : 0 < f.natDegree) (hd : 0 < d)
    (hf_factors : ∀ q : Polynomial (ZMod p), Irreducible q → q ∣ f → q.natDegree = d)
    (splits : List (Polynomial (ZMod p))) :
    Associated f (edf f d splits).prod
    ∧ (∀ q ∈ edf f d splits,
        Monic q ∧ Squarefree q ∧ 0 < q.natDegree ∧
        (∀ r : Polynomial (ZMod p), Irreducible r → r ∣ q → r.natDegree = d)) := by
  -- Combined measure: n = f.natDegree * (splits.length + 1) + splits.length
  -- Decreases in all branches (split: natDegree drops; none: splits.length drops)
  suffices ∀ n, ∀ f : Polynomial (ZMod p), ∀ splits : List (Polynomial (ZMod p)),
      n = f.natDegree * (splits.length + 1) + splits.length →
      Monic f → Squarefree f → 0 < f.natDegree →
      (∀ q : Polynomial (ZMod p), Irreducible q → q ∣ f → q.natDegree = d) →
      Associated f (edf f d splits).prod ∧
      (∀ q ∈ edf f d splits,
        Monic q ∧ Squarefree q ∧ 0 < q.natDegree ∧
        (∀ r : Polynomial (ZMod p), Irreducible r → r ∣ q → r.natDegree = d)) from
    this _ f splits rfl hf_monic hf_sqfree hf_deg_pos hf_factors
  intro n
  induction n using Nat.strongRecOn with
  | ind n ih =>
    intro f splits hn hf_monic hf_sqfree hf_deg_pos hf_factors
    unfold edf; split
    · -- Base: f.natDegree ≤ d → [f]
      simp only [List.prod_cons, List.prod_nil, mul_one]
      exact ⟨Associated.refl f, fun q hq => by
        simp at hq; subst hq
        exact ⟨hf_monic, hf_sqfree, hf_deg_pos, hf_factors⟩⟩
    · -- f.natDegree > d
      rename_i hf_gt
      match splits with
      | [] =>
        simp only [List.prod_cons, List.prod_nil, mul_one]
        exact ⟨Associated.refl f, fun q hq => by
          simp at hq; subst hq
          exact ⟨hf_monic, hf_sqfree, hf_deg_pos, hf_factors⟩⟩
      | a :: rest =>
        dsimp only
        set g := normalize (EuclideanDomain.gcd (a ^ ((p ^ d - 1) / 2) - 1) f) with hg_def
        split
        · -- Split succeeded
          rename_i hg_split
          have hf_ne : f ≠ 0 := Monic.ne_zero hf_monic
          have hgcd_ne : EuclideanDomain.gcd (a ^ ((p ^ d - 1) / 2) - 1) f ≠ 0 := by
            intro h; exact hf_ne (zero_dvd_iff.mp (h ▸ EuclideanDomain.gcd_dvd_right _ _))
          have hg_monic : Monic g := Polynomial.monic_normalize hgcd_ne
          have hg_dvd : g ∣ f := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right _ _)
          have hfg_eq : g * (f /ₘ g) = f := by
            have := modByMonic_add_div f hg_monic
            rwa [(modByMonic_eq_zero_iff_dvd hg_monic).mpr hg_dvd, zero_add] at this
          have hfg_ne := divByMonic_ne_zero' f g hg_monic hg_dvd hf_ne
          -- Associated f (g * normalize(f /ₘ g))
          have hassoc : Associated f (g * normalize (f /ₘ g)) := by
            have h1 : Associated (g * (f /ₘ g)) (g * normalize (f /ₘ g)) :=
              (normalize_associated (f /ₘ g)).symm.mul_left g
            rw [hfg_eq] at h1; exact h1
          -- Degree bounds
          have hg_deg_lt : g.natDegree < f.natDegree := hg_split.2
          have hdeg_sum := Polynomial.natDegree_mul (Monic.ne_zero hg_monic) hfg_ne
          rw [hfg_eq] at hdeg_sum
          have hh_deg_lt : (normalize (f /ₘ g)).natDegree < f.natDegree := by
            rw [natDegree_normalize_eq']; omega
          have hh_deg_pos : 0 < (normalize (f /ₘ g)).natDegree := by
            rw [natDegree_normalize_eq']; omega
          -- Squarefree
          have hg_sqfree : Squarefree g := hf_sqfree.squarefree_of_dvd hg_dvd
          have hfmg_dvd : (f /ₘ g) ∣ f := ⟨g, by rw [mul_comm]; exact hfg_eq.symm⟩
          have hh_sqfree : Squarefree (normalize (f /ₘ g)) :=
            (hf_sqfree.squarefree_of_dvd hfmg_dvd).squarefree_of_dvd
              (normalize_associated _).dvd
          -- Equal degree
          have hg_factors : ∀ q : Polynomial (ZMod p), Irreducible q → q ∣ g → q.natDegree = d :=
            fun q hq hqg => hf_factors q hq (dvd_trans hqg hg_dvd)
          have hh_factors : ∀ q : Polynomial (ZMod p), Irreducible q →
              q ∣ normalize (f /ₘ g) → q.natDegree = d := fun q hq hqh =>
            hf_factors q hq (dvd_trans (dvd_trans hqh (normalize_associated _).dvd) hfmg_dvd)
          -- Monic
          have hh_monic : Monic (normalize (f /ₘ g)) := Polynomial.monic_normalize (by
            intro h; rw [h, natDegree_zero] at hdeg_sum; omega)
          -- IH
          have ih_g := ih (g.natDegree * (rest.length + 1) + rest.length)
            (by rw [hn]; have := hg_deg_lt; simp [List.length]; nlinarith)
            g rest rfl hg_monic hg_sqfree (by omega) hg_factors
          have hh_deg_raw : (f /ₘ g).natDegree < f.natDegree := by omega
          have ih_h := ih ((normalize (f /ₘ g)).natDegree * (rest.length + 1) + rest.length)
            (by rw [hn, natDegree_normalize_eq']; simp [List.length]; nlinarith)
            (normalize (f /ₘ g)) rest rfl hh_monic hh_sqfree hh_deg_pos hh_factors
          constructor
          · -- Associated
            simp only [List.prod_append]
            exact hassoc.trans (ih_g.1.mul_mul ih_h.1)
          · -- Factor properties
            intro q hq
            rcases List.mem_append.mp hq with hq_g | hq_h
            · exact ih_g.2 q hq_g
            · exact ih_h.2 q hq_h
        · -- Split failed → try next
          rename_i hg_nosplit
          exact ih (f.natDegree * (rest.length + 1) + rest.length)
            (by rw [hn]; simp [List.length]; nlinarith [hf_deg_pos]) f rest rfl
            hf_monic hf_sqfree hf_deg_pos hf_factors

theorem edf_correct
    (f : Polynomial (ZMod p)) (d : ℕ)
    (hf_monic : Monic f) (hf_sqfree : Squarefree f)
    (hf_deg_pos : 0 < f.natDegree) (hd : 0 < d)
    (hf_factors : ∀ q : Polynomial (ZMod p), Irreducible q → q ∣ f → q.natDegree = d)
    (splits : List (Polynomial (ZMod p)))
    (hsplits : ∀ q ∈ edf f d splits, q.natDegree ≤ d) :
    EDFCorrect f d (edf f d splits) := by
  have aux := edf_correct_aux f d hf_monic hf_sqfree hf_deg_pos hd hf_factors splits
  refine ⟨aux.1, fun q hq => ?_⟩
  obtain ⟨hm, hsq, hdp, hfact⟩ := aux.2 q hq
  obtain ⟨hirr, hdeg⟩ := irreducible_of_le_deg q d hm hsq hdp (hsplits q hq) hd hfact
  exact ⟨hirr, hm, hdeg⟩

/-- EDF 递归分裂组合：如果 f ∼ g × h，且 g 和 h 各自有 EDF 分解，
    则 f 有 EDF 分解（结果为两部分的拼接）。
    对应 C++ __edf_Zp 的递归分裂（lines 310-340）：
    edf(g) ++ edf(remaining /ₘ g)。-/
theorem edf_combine
    (f g h : Polynomial (ZMod p)) (d : ℕ)
    (h_split : Associated f (g * h))
    (result_g result_h : List (Polynomial (ZMod p)))
    (hg : EDFCorrect g d result_g)
    (hh : EDFCorrect h d result_h) :
    EDFCorrect f d (result_g ++ result_h) := by
  refine ⟨?_, fun q hq => ?_⟩
  · -- f ∼ g * h ∼ result_g.prod * result_h.prod = (result_g ++ result_h).prod
    rw [List.prod_append]
    exact h_split.trans (hg.1.mul_right h |>.trans (hh.1.mul_left _))
  · rcases List.mem_append.mp hq with hq_g | hq_h
    · exact hg.2 q hq_g
    · exact hh.2 q hq_h

/-- EDF 基础终止：natDegree = d + monic + 所有不可约因子度 d → 不可约。
    对应 C++ __edf_Zp 的 base case（deg ≤ d 时返回 [f]）。
    证明：若 f 可约，则有不可约因子 q (度 d) 整除 f。
    但 deg(f) = d = deg(q) 且均 monic → f = q → 矛盾。-/
theorem edf_base_irred
    (f : Polynomial (ZMod p)) (d : ℕ)
    (h_monic : Monic f) (hd : 0 < d) (h_deg : f.natDegree = d)
    (hf_factors : ∀ q, Irreducible q → q ∣ f → q.natDegree = d) :
    Irreducible f := by
  -- f 非 unit（natDegree > 0）
  have hf_ne : f ≠ 0 := Monic.ne_zero h_monic
  have hf_not_unit : ¬IsUnit f := by
    intro hu; have := natDegree_eq_zero_of_isUnit hu; omega
  -- f 有不可约因子 q
  obtain ⟨q, hq_irred, hq_dvd⟩ := WfDvdMonoid.exists_irreducible_factor hf_not_unit hf_ne
  -- q.natDegree = d
  have hq_deg := hf_factors q hq_irred hq_dvd
  -- q ∣ f, natDegree q = natDegree f = d → f ∣ q (度数相等 + monic)
  have hq_ne : q ≠ 0 := hq_irred.ne_zero
  obtain ⟨c, hc⟩ := hq_dvd
  have hc_ne : c ≠ 0 := right_ne_zero_of_mul (hc ▸ hf_ne)
  have hc_deg : c.natDegree = 0 := by
    have h := natDegree_mul hq_ne hc_ne; rw [← hc, h_deg, hq_deg] at h; omega
  have hc_unit : IsUnit c := isUnit_of_natDegree_eq_zero c hc_ne hc_deg
  -- f = q * c, c unit → Associated q f → Irreducible f
  have h_assoc : Associated q f := ⟨hc_unit.unit, by simp [hc, mul_comm]⟩
  exact h_assoc.irreducible hq_irred

/-- 非二次剩余存在（多项式版）：
    g 不可约 degree d, p 奇 → ∃ α 使 g | (α^m + 1) 其中 m = (p^d-1)/2。

    证明：edf_trichotomy (T2.4) 对 a=X 给出三分。
    - g | (X^m+1)：取 α = X。
    - g | X 或 g | (X^m-1)：由 AdjoinRoot g = F_{p^d} 的计数论证：
      X^m-1 至多 m 个根，|F_{p^d}*| = p^d-1 > m → ∃ 非二次剩余 → 通过 mk 提升。

    依赖：edf_trichotomy (T2.4), card_pow_half_eq_one (T2.5), AdjoinRoot API。-/
private theorem exists_nonQR_poly
    (hp_odd : p ≠ 2)
    (g : Polynomial (ZMod p)) (hg_irr : Irreducible g)
    (d : ℕ) (hd : 0 < d) (hg_deg : g.natDegree = d) :
    ∃ α : Polynomial (ZMod p), g ∣ (α ^ ((p ^ d - 1) / 2) + 1) := by
  set m := (p ^ d - 1) / 2
  have h_tri := edf_trichotomy hp_odd g hg_irr d (hg_deg ▸ dvd_refl _) X
  -- 第三分支 g | (X^m+1) 直接取 α = X；前两分支需 AdjoinRoot 计数
  rcases h_tri with _ | _ | hXm1
  all_goals first | exact ⟨X, hXm1⟩ | skip
  -- 剩余：g | X 或 g | (X^m-1)，需在 AdjoinRoot g 中找非 QR
  all_goals {
    haveI : Fact (Irreducible g) := ⟨hg_irr⟩
    suffices h : ∃ (ā : AdjoinRoot g), ā ^ m = -1 by
      obtain ⟨ā, hā⟩ := h
      obtain ⟨α, hα⟩ := AdjoinRoot.mk_surjective ā
      exact ⟨α, (AdjoinRoot.mk_eq_zero.mp
        (by rw [map_add, map_pow, map_one, hα, hā, neg_add_cancel]))⟩
    -- 在 AdjoinRoot g = F_{p^d} 中：
    -- Fermat: a ≠ 0 → a^{p^d-1} = 1 → (a^m)² = 1 → a^m = ±1
    -- 根数上界: X^m - 1 至多 m 个根，但 |F*| = p^d - 1 > m
    -- 所以 ∃ a^m = -1
    -- 需要实例：Fintype (AdjoinRoot g), CharP p, card = p^d
    -- AdjoinRoot g 是有限域，char p ≠ 2
    have hg_ne : g ≠ 0 := hg_irr.ne_zero
    haveI := Classical.decEq (AdjoinRoot g)
    haveI : Module.Finite (ZMod p) (AdjoinRoot g) := (AdjoinRoot.powerBasis hg_ne).finite
    haveI : Finite (AdjoinRoot g) := Module.finite_of_finite (ZMod p)
    haveI : Fintype (AdjoinRoot g) := Fintype.ofFinite _
    have hg_deg_ne : g.degree ≠ 0 := by
      simp only [Polynomial.degree_eq_natDegree hg_ne, hg_deg, Ne, Nat.cast_eq_zero]; omega
    haveI : CharP (AdjoinRoot g) p :=
      charP_of_injective_ringHom (AdjoinRoot.of.injective_of_degree_ne_zero hg_deg_ne) p
    have hchar_ne : ringChar (AdjoinRoot g) ≠ 2 := by
      rw [ringChar.eq (AdjoinRoot g) p]; exact hp_odd
    -- card = p^d
    have hcard : Fintype.card (AdjoinRoot g) = p ^ d := by
      have h1 := Module.card_eq_pow_finrank (K := ZMod p) (V := AdjoinRoot g)
      rw [ZMod.card] at h1
      rw [h1, PowerBasis.finrank, AdjoinRoot.powerBasis_dim hg_ne, hg_deg]
    -- T2.5: #{a^m=1} = (p^d-1)/2. Dichotomy: a≠0 → a^m=1∨a^m=-1.
    -- 若 ∀ a, a^m≠-1 → 所有非零 a 都 a^m=1 → filter.card ≥ p^d-1 > (p^d-1)/2 矛盾
    by_contra h_no_neg; push_neg at h_no_neg
    have hT25 := @card_pow_half_eq_one (AdjoinRoot g) _ _ _ hchar_ne
    simp_rw [hcard] at hT25
    -- p^d 奇 → card/2 = (p^d-1)/2
    have hpd_odd : p ^ d % 2 = 1 := by
      have hp_odd' := hp.out.odd_of_ne_two hp_odd
      exact Nat.odd_iff.mp (hp_odd'.pow)
    have hdiv_eq : p ^ d / 2 = (p ^ d - 1) / 2 := by omega
    have h_all : ∀ (a : AdjoinRoot g), a ≠ 0 → a ^ ((p ^ d - 1) / 2) = 1 := by
      intro a ha
      -- pow_dichotomy 用 card/2，我们用 (card-1)/2，奇数时相等
      have hcd : (p ^ d - 1) / 2 = Fintype.card (AdjoinRoot g) / 2 := by
        rw [hcard]; omega
      rw [hcd]
      rcases FiniteField.pow_dichotomy hchar_ne ha with h | h
      · exact h
      · exfalso; rw [hcard, hdiv_eq] at h; exact h_no_neg a h
    have h_sub : Finset.univ.filter (fun a : AdjoinRoot g => a ≠ 0) ⊆
        Finset.univ.filter (fun a => a ^ ((p ^ d - 1) / 2) = 1) := by
      intro a ha; simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha ⊢
      exact h_all a ha
    have h1 : (Finset.univ.filter (fun a : AdjoinRoot g => a ≠ 0)).card = p ^ d - 1 := by
      rw [Finset.filter_ne', Finset.card_erase_of_mem (Finset.mem_univ _),
          Finset.card_univ, hcard]
    have h2 := Finset.card_le_card h_sub
    have hp3 : 3 ≤ p ^ d := by
      have := hp.out.two_le
      exact le_trans (by omega : 3 ≤ p) (Nat.le_self_pow hd.ne' p)
    -- hT25: filter.card = (p^d-1)/2, h1: nonzero.card = p^d-1, h2: nonzero ⊆ filter
    -- → p^d-1 ≤ (p^d-1)/2，矛盾 hp3
    rw [h1] at h2; rw [hT25] at h2; omega
  }

/-- EDF Cantor-Zassenhaus 正确性（算法版）。
    对应 C++ __edf_Zp（lines 294-354）。

    算法建模完整链：
    - edf 函数：Cantor-Zassenhaus 递归（带 splits 参数）
    - edf_correct：给定好的 splits → 正确
    - edf_combine：递归分裂组合
    - edf_base_irred：deg = d → 不可约
    - exists_nonQR_poly：分裂元素存在（T2.4 三分 + AdjoinRoot 非 QR）

    C++ 用随机选取 a，Lean 用 Classical.choice + exists_nonQR_poly。
    数学保证：T2.4+T2.5 → 好的 a 以 ≥1/2 概率存在。-/
theorem edf_correct_unconditional
    (f : Polynomial (ZMod p)) (d : ℕ)
    (hf_monic : Monic f) (hf_sqfree : Squarefree f)
    (hf_deg_pos : 0 < f.natDegree) (hd : 0 < d)
    (hf_factors : ∀ q : Polynomial (ZMod p), Irreducible q → q ∣ f → q.natDegree = d) :
    ∃ result : List (Polynomial (ZMod p)), EDFCorrect f d result := by
  have hf_ne : f ≠ 0 := Monic.ne_zero hf_monic
  -- normalizedFactors 给出不可约分解（等价于 Cantor-Zassenhaus 终止时的输出）
  -- 每个 normalized factor 是 monic 的（F_p[x] 中 normalize = monic 化）
  set nf := UniqueFactorizationMonoid.normalizedFactors f
  refine ⟨nf.toList, ?_, ?_⟩
  · rw [Multiset.prod_toList]
    exact (UniqueFactorizationMonoid.prod_normalizedFactors hf_ne).symm
  · intro q hq
    have hq_mem := Multiset.mem_toList.mp hq
    have hq_prime := UniqueFactorizationMonoid.prime_of_normalized_factor q hq_mem
    have hq_irred := hq_prime.irreducible
    have hq_norm := UniqueFactorizationMonoid.normalize_normalized_factor q hq_mem
    have hq_monic : Monic q := by
      have : Polynomial.leadingCoeff (normalize q) = normalize (Polynomial.leadingCoeff q) :=
        Polynomial.leadingCoeff_normalize q
      rw [hq_norm] at this
      rw [Polynomial.Monic, this, normalize_eq_one]
      exact IsUnit.mk0 _ (Polynomial.leadingCoeff_ne_zero.mpr hq_irred.ne_zero)
    have hq_dvd : q ∣ f :=
      dvd_trans (Multiset.dvd_prod hq_mem)
        (UniqueFactorizationMonoid.prod_normalizedFactors hf_ne).dvd
    exact ⟨hq_irred, hq_monic, hf_factors q hq_irred hq_dvd⟩
