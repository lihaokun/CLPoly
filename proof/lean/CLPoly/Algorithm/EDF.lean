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
