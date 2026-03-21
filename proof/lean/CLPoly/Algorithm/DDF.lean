/-
  CLPoly/Algorithm/DDF.lean — L2 DDF 算法模型与正确性证明

  Phase 3 T3.1: Distinct Degree Factorization
  对应 C++: polynomial_factorize_zp.hh:247-292 __ddf_Zp

  证明结构：
  1. ddfLoop: 递归函数定义 + 终止性证明
  2. ddf: 入口包装
  3. 辅助引理: h-congruence step, gd characterization, etc.
  4. ddfLoop_correct: 归纳证明（携带不变量 P0-P6）
  5. ddf_correct: 顶层定理
-/
import CLPoly.Spec
import CLPoly.Math.FiniteFieldFact
import Mathlib.RingTheory.UniqueFactorizationDomain.Defs

set_option autoImplicit false
set_option maxHeartbeats 800000

open Polynomial

variable {p : ℕ} [hp : Fact (Nat.Prime p)]

-- ============================================================
-- 0. decreasing_by 辅助引理（需在 ddfLoop 之前定义）
-- ============================================================

/-- Monic g 整除 f 且 deg(g) > 0 且 f ≠ 0 → deg(f /ₘ g) < deg(f) -/
private lemma natDegree_divByMonic_lt
    (f g : Polynomial (ZMod p)) (hg : Monic g)
    (hdvd : g ∣ f) (hg_pos : 0 < g.natDegree) (hf : f ≠ 0) :
    (f /ₘ g).natDegree < f.natDegree := by
  have hg_ne : g ≠ 0 := Monic.ne_zero hg
  have hmod : f %ₘ g = 0 := (modByMonic_eq_zero_iff_dvd hg).mpr hdvd
  have heq : g * (f /ₘ g) = f := by
    have := modByMonic_add_div f hg; rwa [hmod, zero_add] at this
  have hfn_ne : f /ₘ g ≠ 0 := by
    intro h; rw [h, mul_zero] at heq; exact hf heq.symm
  have hdeg := Polynomial.natDegree_mul hg_ne hfn_ne
  rw [heq] at hdeg
  omega

-- ============================================================
-- 1. 算法模型
-- ============================================================

/-- DDF 内循环：递归模型对应 C++ __ddf_Zp 的 for 循环。
    h 是 X^{p^{d-1}} mod f_star 的代表元，d 当前度数，acc 累积结果。-/
noncomputable def ddfLoop
    (h f_star : Polynomial (ZMod p)) (d : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    : List (Polynomial (ZMod p) × ℕ) :=
  if _hterm : f_star.natDegree < 2 * d then
    if 0 < f_star.natDegree then acc ++ [(f_star, f_star.natDegree)]
    else acc
  else
    let h' := (h ^ p) %ₘ f_star
    let gd := normalize (EuclideanDomain.gcd (h' - X) f_star)
    if _hsplit : 0 < gd.natDegree then
      let f_new := f_star /ₘ gd
      ddfLoop (h' %ₘ f_new) f_new (d + 1) (acc ++ [(gd, d)])
    else
      ddfLoop h' f_star (d + 1) acc
termination_by f_star.natDegree + 1 - 2 * d
decreasing_by
  · -- Split case: f_new = f_star /ₘ gd, gd.natDegree > 0
    set gd_val := normalize (EuclideanDomain.gcd ((h ^ p) %ₘ f_star - X) f_star) with hgd_def
    have hgd_pos : 0 < gd_val.natDegree := _hsplit
    have hgd_raw_ne : EuclideanDomain.gcd ((h ^ p) %ₘ f_star - X) f_star ≠ 0 := by
      intro heq; simp [gd_val, heq, normalize_zero] at hgd_pos
    have hgd_monic : Monic gd_val := by
      rw [hgd_def]; exact Polynomial.monic_normalize hgd_raw_ne
    have hgd_dvd : gd_val ∣ f_star := by
      rw [hgd_def]; exact normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right _ _)
    by_cases hfs : f_star = 0
    · simp [hfs] at _hterm ⊢; omega
    · have := natDegree_divByMonic_lt f_star gd_val hgd_monic hgd_dvd hgd_pos hfs
      omega
  · -- No-split case: same f_star, d → d+1
    omega

/-- DDF 入口：对首一无平方多项式 f 执行 DDF -/
noncomputable def ddf (f : Polynomial (ZMod p)) :
    List (Polynomial (ZMod p) × ℕ) :=
  ddfLoop X f 1 []

-- ============================================================
-- 2. 辅助引理
-- ============================================================

/-- normalize 非零多项式是 monic -/
private lemma monic_normalize_of_ne_zero
    (f : Polynomial (ZMod p)) (hf : f ≠ 0) :
    Monic (normalize f) :=
  Polynomial.monic_normalize hf

/-- Monic 的 divByMonic 给出精确乘积分解 -/
private lemma monic_divByMonic_mul_eq
    (f g : Polynomial (ZMod p)) (hg : Monic g) (hdvd : g ∣ f) :
    g * (f /ₘ g) = f := by
  have h0 : f %ₘ g = 0 := (modByMonic_eq_zero_iff_dvd hg).mpr hdvd
  have := modByMonic_add_div f hg
  rw [h0, zero_add] at this
  exact this

/-- Prime 元素整除列表积 → 整除某个元素 -/
private lemma prime_dvd_list_prod
    {q : Polynomial (ZMod p)} (hq : Prime q)
    {l : List (Polynomial (ZMod p))} (h : q ∣ l.prod) :
    ∃ a ∈ l, q ∣ a := by
  induction l with
  | nil => simp at h; exact absurd (isUnit_of_dvd_one h) hq.not_unit
  | cons x xs ih =>
    rw [List.prod_cons] at h
    rcases hq.dvd_or_dvd h with hx | hxs
    · exact ⟨x, List.Mem.head xs, hx⟩
    · obtain ⟨a, ha, hdvd⟩ := ih hxs
      exact ⟨a, List.Mem.tail x ha, hdvd⟩

-- ============================================================
-- 3. 核心桥接引理
-- ============================================================

/-- h-congruence 递推：若 f_star ∣ (h - X^{p^{d-1}})，
    则 f_star ∣ (h^p %ₘ f_star - X^{p^d})。
    证明：(a-b) | (a^p - b^p) + modByMonic 分解 -/
private lemma h_cong_step
    (h f_star : Polynomial (ZMod p)) (d : ℕ) (hd : d ≥ 1)
    (hmonic : Monic f_star)
    (hcong : f_star ∣ (h - X ^ (p ^ (d - 1)))) :
    f_star ∣ ((h ^ p) %ₘ f_star - X ^ (p ^ d)) := by
  -- Step A: f_star ∣ (h^p - X^{p^d})
  have hdvd_sub : (h - X ^ (p ^ (d - 1))) ∣
      (h ^ p - (X ^ (p ^ (d - 1))) ^ p) :=
    Commute.sub_dvd_pow_sub_pow (Commute.all _ _) p
  have hpow_eq : ((X : Polynomial (ZMod p)) ^ (p ^ (d - 1))) ^ p =
      (X : Polynomial (ZMod p)) ^ (p ^ d) := by
    rw [← pow_mul]
    congr 1
    conv_rhs => rw [show d = d - 1 + 1 from (Nat.sub_add_cancel hd).symm]
    exact (pow_succ p (d - 1)).symm
  rw [hpow_eq] at hdvd_sub
  have step_a : f_star ∣ (h ^ p - X ^ (p ^ d)) := dvd_trans hcong hdvd_sub
  -- Step B: f_star ∣ (h^p - h^p %ₘ f_star)
  have step_b : f_star ∣ (h ^ p - (h ^ p) %ₘ f_star) := by
    refine ⟨(h ^ p) /ₘ f_star, ?_⟩
    have := (modByMonic_add_div (h ^ p) hmonic).symm
    rw [add_comm] at this
    exact sub_eq_of_eq_add this
  -- Step C: combine
  have hdecomp : ((h ^ p) %ₘ f_star - X ^ (p ^ d)) =
      -(h ^ p - (h ^ p) %ₘ f_star) + (h ^ p - X ^ (p ^ d)) := by ring
  rw [hdecomp]
  exact dvd_add (dvd_neg.mpr step_b) step_a

/-- 核心桥接：gd 的不可约因子恰为 f_star 中度 = d 的不可约因子。
    q ∣ gd ↔ q ∣ f_star ∧ q.natDegree = d -/
private lemma gd_irred_characterization
    (f_star : Polynomial (ZMod p))
    (h' : Polynomial (ZMod p))
    (d : ℕ) (hd : d ≥ 1)
    (_hmonic : Monic f_star) (_hsqf : Squarefree f_star)
    (hcong : f_star ∣ (h' - X ^ (p ^ d)))
    (hfac : ∀ q : Polynomial (ZMod p), Irreducible q → q ∣ f_star → q.natDegree ≥ d)
    (q : Polynomial (ZMod p)) (hq : Irreducible q) :
    let gd := normalize (EuclideanDomain.gcd (h' - X) f_star)
    q ∣ gd ↔ q ∣ f_star ∧ q.natDegree = d := by
  intro gd
  -- Step 1: q ∣ gd ↔ q ∣ gcd_raw
  have h_norm : q ∣ gd ↔ q ∣ EuclideanDomain.gcd (h' - X) f_star :=
    dvd_normalize_iff
  -- Step 2: q ∣ gcd_raw ↔ q ∣ (h'-X) ∧ q ∣ f_star
  have h_gcd : q ∣ EuclideanDomain.gcd (h' - X) f_star ↔
      q ∣ (h' - X) ∧ q ∣ f_star := by
    constructor
    · intro hdvd
      exact ⟨dvd_trans hdvd (EuclideanDomain.gcd_dvd_left _ _),
             dvd_trans hdvd (EuclideanDomain.gcd_dvd_right _ _)⟩
    · intro ⟨h1, h2⟩; exact EuclideanDomain.dvd_gcd h1 h2
  -- Step 3: swap h'-X ↔ X^{p^d}-X using congruence
  have h_swap : q ∣ (h' - X) ∧ q ∣ f_star ↔
      q ∣ (X ^ (p ^ d) - X) ∧ q ∣ f_star := by
    constructor
    · intro ⟨hqhx, hqfs⟩
      have hq_cong : q ∣ (h' - X ^ (p ^ d)) := dvd_trans hqfs hcong
      have hkey : (X ^ (p ^ d) - X : Polynomial (ZMod p)) =
          -(h' - X ^ (p ^ d)) + (h' - X) := by ring
      exact ⟨by rw [hkey]; exact dvd_add (dvd_neg.mpr hq_cong) hqhx, hqfs⟩
    · intro ⟨hqxp, hqfs⟩
      have hq_cong : q ∣ (h' - X ^ (p ^ d)) := dvd_trans hqfs hcong
      have hkey : (h' - X : Polynomial (ZMod p)) =
          (h' - X ^ (p ^ d)) + (X ^ (p ^ d) - X) := by ring
      exact ⟨by rw [hkey]; exact dvd_add hq_cong hqxp, hqfs⟩
  -- Step 4: T2.2 + T2.2'
  have h_t22 : q ∣ (X ^ (p ^ d) - X : Polynomial (ZMod p)) ↔ q.natDegree ∣ d :=
    ⟨dvd_X_pow_sub_X_imp_natDegree_dvd q hq d,
     fun hdvd => irreducible_dvd_X_pow_sub_X q hq d hdvd⟩
  -- Step 5: deg ∣ d ∧ q ∣ f* ↔ deg = d ∧ q ∣ f* (using P5)
  have h_eq : q.natDegree ∣ d ∧ q ∣ f_star ↔ q.natDegree = d ∧ q ∣ f_star := by
    constructor
    · intro ⟨hdvd_d, hqfs⟩
      exact ⟨Nat.le_antisymm (Nat.le_of_dvd (by omega) hdvd_d) (hfac q hq hqfs), hqfs⟩
    · intro ⟨heq, hqfs⟩; exact ⟨heq ▸ dvd_refl _, hqfs⟩
  -- Combine
  rw [h_norm, h_gcd]
  constructor
  · intro ⟨hqhx, hqfs⟩
    have hqxp := (h_swap.mp ⟨hqhx, hqfs⟩).1
    exact ⟨hqfs, (h_eq.mp ⟨h_t22.mp hqxp, hqfs⟩).1⟩
  · intro ⟨hqfs, heq⟩
    have hqxp := h_t22.mpr (h_eq.mpr ⟨heq, hqfs⟩).1
    exact (h_swap.mpr ⟨hqxp, hqfs⟩)

-- ============================================================
-- 4. 主定理
-- ============================================================

/-- ddfLoop 正确性：携带不变量 P0-P6 的归纳证明 -/
theorem ddfLoop_correct
    (f h f_star : Polynomial (ZMod p)) (d : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    (hd : d ≥ 1)                                          -- P0
    (hprod : f = f_star * (acc.map Prod.fst).prod)         -- P1
    (hcong : f_star ∣ (h - X ^ (p ^ (d - 1))))            -- P2
    (hmonic : Monic f_star)                                -- P3
    (hsqf : Squarefree f_star)                             -- P4
    (hfac : ∀ q : Polynomial (ZMod p),
        Irreducible q → q ∣ f_star → q.natDegree ≥ d)    -- P5
    (hacc : ∀ pr ∈ acc, pr.1 ∣ f ∧ Monic pr.1 ∧
            (∀ q : Polynomial (ZMod p),
             Irreducible q → q ∣ pr.1 → q.natDegree = pr.2)) -- P6
    : DDFCorrect f (ddfLoop h f_star d acc) := by
  -- Revert invariants for functional induction motive
  revert hd hprod hcong hmonic hsqf hfac hacc
  induction h, f_star, d, acc using ddfLoop.induct
  -- Case 1: Termination, deg(f_star) > 0 (6 inaccessible names)
  · rename_i h f_star d acc hterm hdeg_pos
    intro hd hprod hcong hmonic hsqf hfac hacc
    unfold ddfLoop; rw [dif_pos hterm, if_pos hdeg_pos]
    -- §5.1: f_star is irreducible
    -- deg < 2d, squarefree, all irred factors ≥ d → can't have 2+ factors
    have hfs_ne : f_star ≠ 0 := Monic.ne_zero hmonic
    have hfs_nu : ¬ IsUnit f_star := by
      intro hu; exact absurd (natDegree_eq_zero_of_isUnit hu ▸ hdeg_pos) (lt_irrefl 0)
    have hirr : Irreducible f_star := by
      rw [irreducible_iff]; exact ⟨hfs_nu, fun a b hab => by
        by_contra h_both; push_neg at h_both
        obtain ⟨ha_nu, hb_nu⟩ := h_both
        have ha_ne : a ≠ 0 := by rintro rfl; exact hfs_ne (by rw [hab, zero_mul])
        have hb_ne : b ≠ 0 := by rintro rfl; exact hfs_ne (by rw [hab, mul_zero])
        obtain ⟨qa, hqa_irr, hqa_dvd⟩ := WfDvdMonoid.exists_irreducible_factor ha_nu ha_ne
        obtain ⟨qb, hqb_irr, hqb_dvd⟩ := WfDvdMonoid.exists_irreducible_factor hb_nu hb_ne
        have hqa_ge := hfac qa hqa_irr (dvd_trans hqa_dvd (hab ▸ dvd_mul_right a b))
        have hqb_ge := hfac qb hqb_irr (dvd_trans hqb_dvd (hab ▸ dvd_mul_left b a))
        have hqa_le := Polynomial.natDegree_le_of_dvd hqa_dvd ha_ne
        have hqb_le := Polynomial.natDegree_le_of_dvd hqb_dvd hb_ne
        have hdeg_sum := hab ▸ Polynomial.natDegree_mul ha_ne hb_ne
        omega⟩
    -- DDFCorrect for acc ++ [(f_star, f_star.natDegree)]
    refine ⟨?_, ?_, ?_, ?_, ?_⟩
    · -- Cond 1: pr.1 ∣ f
      intro pr hpr; simp at hpr
      rcases hpr with h_acc | h_new
      · exact (hacc pr h_acc).1
      · rw [h_new]; exact ⟨(acc.map Prod.fst).prod, hprod⟩
    · -- Cond 2: q irred, q ∣ pr.1 → deg q = pr.2
      intro pr hpr q hq hqdvd; simp at hpr
      rcases hpr with h_acc | h_new
      · exact (hacc pr h_acc).2.2 q hq hqdvd
      · -- pr = (f_star, natDegree), q irred ∣ f_star irred → same deg
        rw [h_new] at hqdvd ⊢; simp at hqdvd ⊢
        -- q ∣ f_star, f_star irreducible: f_star = q * r → IsUnit q ∨ IsUnit r
        obtain ⟨r, hr⟩ := hqdvd  -- hr : f_star = q * r
        rcases hirr.isUnit_or_isUnit hr with hu | hu
        · exact absurd hu hq.1  -- q is not a unit
        · -- IsUnit r → r.natDegree = 0 → q.natDegree = f_star.natDegree
          have hq_ne : q ≠ 0 := by rintro rfl; exact hfs_ne (by rw [hr, zero_mul])
          have hr_ne : r ≠ 0 := by rintro rfl; exact hfs_ne (by rw [hr, mul_zero])
          have hdeg := Polynomial.natDegree_mul hq_ne hr_ne
          rw [← hr, natDegree_eq_zero_of_isUnit hu] at hdeg
          omega
    · -- Cond 3: completeness
      intro q hq hqf; rw [hprod] at hqf
      rcases hq.prime.dvd_or_dvd hqf with hq_fs | hq_acc
      · exact ⟨(f_star, f_star.natDegree),
          List.mem_append.mpr (Or.inr (List.mem_singleton.mpr rfl)), hq_fs⟩
      · obtain ⟨a, ha_mem, hqa⟩ := prime_dvd_list_prod hq.prime hq_acc
        obtain ⟨pr, hpr_mem, rfl⟩ := List.mem_map.mp ha_mem
        exact ⟨pr, List.mem_append.mpr (Or.inl hpr_mem), hqa⟩
    · -- Cond 4: Associated f result.prod
      simp only [List.map_append, List.prod_append, List.map_singleton, List.prod_singleton]
      rw [hprod, mul_comm]
    · -- Cond 5: Monic
      intro pr hpr; simp at hpr
      rcases hpr with h_acc | h_new
      · exact (hacc pr h_acc).2.1
      · rw [h_new]; exact hmonic
  -- Case 2: Termination, deg(f_star) = 0 (6 inaccessible names)
  · rename_i h f_star d acc hterm hdeg_zero
    intro hd hprod hcong hmonic hsqf hfac hacc
    unfold ddfLoop; rw [dif_pos hterm, if_neg hdeg_zero]
    have hfs_one : f_star = 1 :=
      Polynomial.eq_one_of_monic_natDegree_zero hmonic (by omega)
    rw [hfs_one, one_mul] at hprod
    exact ⟨
      fun pr hpr => (hacc pr hpr).1,
      fun pr hpr => (hacc pr hpr).2.2,
      fun q hq hqf => by
        rw [hprod] at hqf
        obtain ⟨a, ha_mem, hqa⟩ := prime_dvd_list_prod hq.prime hqf
        obtain ⟨pr, hpr_mem, rfl⟩ := List.mem_map.mp ha_mem
        exact ⟨pr, hpr_mem, hqa⟩,
      by rw [hprod],
      fun pr hpr => (hacc pr hpr).2.1⟩
  -- Case 3: Split (gd.natDegree > 0) with IH
  · rename_i hv fv dv av hterm h'v gdv hsplit fnv ih
    intro hd hprod hcong hmonic hsqf hfac hacc
    -- Unfold ddfLoop to the split recursive call
    unfold ddfLoop; rw [dif_neg hterm]; dsimp only; rw [dif_pos hsplit]
    -- Establish F1-F4 facts about gd and f_new
    have hfv_ne : fv ≠ 0 := Monic.ne_zero hmonic
    have hgcd_ne : EuclideanDomain.gcd (h'v - X) fv ≠ 0 := by
      intro habs; exact hfv_ne (zero_dvd_iff.mp (habs ▸ EuclideanDomain.gcd_dvd_right _ _))
    -- F1: gd ∣ f_star
    have hgd_dvd : gdv ∣ fv := normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right _ _)
    -- F2: Monic gd
    have hgd_monic : Monic gdv := Polynomial.monic_normalize hgcd_ne
    -- F3: f_star = gd * f_new
    have hF3 : fv = gdv * (fv /ₘ gdv) :=
      (monic_divByMonic_mul_eq fv gdv hgd_monic hgd_dvd).symm
    -- F4: Monic f_new
    have hfn_monic : Monic (fv /ₘ gdv) := by
      have hfn_ne : fv /ₘ gdv ≠ 0 := by
        intro h; rw [h, mul_zero] at hF3; exact hfv_ne hF3
      exact Polynomial.Monic.of_mul_monic_left hgd_monic (hF3 ▸ hmonic)
    -- h-congruence for the new d
    have hcong_d := h_cong_step hv fv dv hd hmonic hcong
    -- Apply IH with updated invariants
    apply ih
    · -- P0: d+1 ≥ 1
      omega
    · -- P1: f = f_new * (acc_new.map Prod.fst).prod
      simp only [List.map_append, List.prod_append, List.map_singleton, List.prod_singleton]
      rw [hprod, hF3]; ring
    · -- P2: fnv ∣ (h'v %ₘ fnv - X^{p^(dv+1-1)})
      -- Simplify dv + 1 - 1 = dv
      show fnv ∣ (h'v %ₘ fnv - X ^ (p ^ (dv + 1 - 1)))
      rw [show dv + 1 - 1 = dv by omega]
      -- From hcong_d: fv ∣ (h'v - X^{p^dv})
      -- fnv ∣ fv → fnv ∣ (h'v - X^{p^dv})
      have hfn_dvd_fs : fnv ∣ fv := ⟨gdv, hF3.trans (mul_comm _ _)⟩
      have hfn_dvd_cong : fnv ∣ (h'v - X ^ (p ^ dv)) := dvd_trans hfn_dvd_fs hcong_d
      -- fnv ∣ (h'v - h'v %ₘ fnv) from modByMonic
      have hmod_dvd : fnv ∣ (h'v - h'v %ₘ fnv) := by
        refine ⟨h'v /ₘ fnv, ?_⟩
        have := (modByMonic_add_div h'v hfn_monic).symm
        rw [add_comm] at this; exact sub_eq_of_eq_add this
      have hdecomp : (h'v %ₘ fnv - X ^ (p ^ dv)) =
          -(h'v - h'v %ₘ fnv) + (h'v - X ^ (p ^ dv)) := by ring
      rw [hdecomp]; exact dvd_add (dvd_neg.mpr hmod_dvd) hfn_dvd_cong
    · -- P3: Monic f_new
      exact hfn_monic
    · -- P4: Squarefree f_new
      exact Squarefree.squarefree_of_dvd ⟨gdv, hF3.trans (mul_comm _ _)⟩ hsqf
    · -- P5: ∀ q irred, q ∣ f_new → q.natDegree ≥ d+1
      intro q hq hqfn
      have hqfs : q ∣ fv := dvd_trans hqfn ⟨gdv, hF3.trans (mul_comm _ _)⟩
      have hge := hfac q hq hqfs
      by_contra hlt; push_neg at hlt
      have heq : q.natDegree = dv := by omega
      -- q ∣ f_star ∧ deg = d → q ∣ gd (by characterization)
      have hq_gd := (gd_irred_characterization fv h'v dv hd hmonic hsqf hcong_d hfac q hq).mpr
        ⟨hqfs, heq⟩
      -- q ∣ gd and q ∣ f_new → q² ∣ f_star (squarefree contradiction)
      have hq_sq : q * q ∣ fv := by
        rw [hF3]; exact mul_dvd_mul hq_gd hqfn
      exact absurd (hsqf q hq_sq) hq.1
    · -- P6: acc_new = acc ++ [(gd, d)]
      intro pr hpr; simp at hpr
      rcases hpr with h_acc | h_new
      · exact hacc pr h_acc
      · rw [h_new]; simp
        refine ⟨dvd_trans hgd_dvd ⟨_, hprod⟩, hgd_monic, ?_⟩
        intro q hq hqdvd
        exact ((gd_irred_characterization fv h'v dv hd hmonic hsqf hcong_d hfac q hq).mp hqdvd).2
  -- Case 4: No-split with IH
  · rename_i hv fv dv av hterm h'v gdv hnsplit ihv
    intro hd hprod hcong hmonic hsqf hfac hacc
    -- Unfold and simplify to recursive call
    unfold ddfLoop; rw [dif_neg hterm]; dsimp only; rw [dif_neg hnsplit]
    -- Apply IH
    apply ihv
    · omega                                    -- P0
    · exact hprod                              -- P1
    · exact h_cong_step hv fv dv hd hmonic hcong  -- P2
    · exact hmonic                             -- P3
    · exact hsqf                               -- P4
    · -- P5: ∀ q irred, q ∣ f_star → q.natDegree ≥ d+1
      intro q hq hqfs
      by_contra hlt; push_neg at hlt
      have hge := hfac q hq hqfs
      have heq : q.natDegree = dv := by omega
      have hcong_d := h_cong_step hv fv dv hd hmonic hcong
      have hchar := (gd_irred_characterization fv ((hv ^ p) %ₘ fv)
        dv hd hmonic hsqf hcong_d hfac q hq).mpr ⟨hqfs, heq⟩
      -- hchar : q ∣ gdv, hnsplit : ¬(0 < gdv.natDegree)
      -- gdv ≠ 0 (f_star ≠ 0 → gcd_dvd_right → gcd ≠ 0 → normalize ≠ 0)
      have hfv_ne : fv ≠ 0 := Monic.ne_zero hmonic
      have hgcd_ne : EuclideanDomain.gcd (h'v - X) fv ≠ 0 := by
        intro habs
        have := habs ▸ EuclideanDomain.gcd_dvd_right (h'v - X) fv
        exact hfv_ne (zero_dvd_iff.mp this)
      have hgd_ne : gdv ≠ 0 := by
        intro habs; exact hgcd_ne (normalize_eq_zero.mp habs)
      -- q.natDegree ≤ gdv.natDegree = 0, but q.natDegree ≥ 1
      have hle := Polynomial.natDegree_le_of_dvd hchar hgd_ne
      have hgd_zero : gdv.natDegree = 0 := Nat.le_zero.mp (not_lt.mp hnsplit)
      rw [hgd_zero] at hle
      exact absurd (Irreducible.natDegree_pos hq) (not_lt.mpr hle)
    · exact hacc                               -- P6

/-- DDF 顶层正确性 -/
theorem ddf_correct
    (f : Polynomial (ZMod p)) (hm : Monic f) (hsq : Squarefree f) :
    DDFCorrect f (ddf f) := by
  unfold ddf
  apply ddfLoop_correct f X f 1 []
  · omega
  · simp
  · simp
  · exact hm
  · exact hsq
  · intro q hq _; exact Irreducible.natDegree_pos hq
  · intro _ hpr; simp at hpr
