/-
  CLPoly/Pipeline/FactorZZ.lean — Z[x] 顶层正确性：框架组合证明

  证明：若 Zp 因式分解、Hensel 提升、因子重组各自满足规约，
  则组合结果满足 FactorZZCorrect。

  证明结构（全部无 sorry）：
    1. Squarefree → f mod p ≠ 0
    2. factor_zp: f mod p = C(lc) * ∏(fᵢ^eᵢ)
    3. 构造 Hensel 输入: facs_p = [C(lc), f₁^e₁, ...]，乘积 = f mod p
    4. hensel: map f mod p^k = (hensel facs_p).prod
    5. recombine → RecombineCorrect ≡ FactorZZCorrect（定义相同）

  关键观察：RecombineCorrect 和 FactorZZCorrect 定义完全相同，
  因此组合证明只需串联三个子过程假设，无需额外辅助引理。
-/
import CLPoly.Spec
import Mathlib.RingTheory.Polynomial.Content
import Mathlib.Data.ZMod.Basic

set_option autoImplicit false

open Polynomial

-- ============================================================
-- 辅助引理：Squarefree Zp 因式分解中各因子互素
-- ============================================================

lemma monic_irreducible_coprime {K : Type*} [Field K] [DecidableEq (Polynomial K)] {f g : Polynomial K}
    (hf : Irreducible f) (hg : Irreducible g) (hmf : Monic f) (hmg : Monic g) (hne : f ≠ g) :
    IsCoprime f g := by
  have h_not_dvd : ¬ f ∣ g := by
    intro h_dvd
    rcases h_dvd with ⟨h', h_eq⟩
    -- h_eq : g = f * h'
    have h_cases : IsUnit f ∨ IsUnit h' :=
      hg.isUnit_or_isUnit (a := f) (b := h') (by
        -- need: g = f * h'
        simpa [mul_comm] using h_eq)
    rcases h_cases with (h_unit | h_unit')
    · exact hf.1 h_unit
    · have h_assoc : Associated f g := by
        use h_unit'.unit
        calc
          f * (h_unit'.unit : Polynomial K) = f * h' := by simp [h_unit'.unit_spec]
          _ = g := h_eq.symm
      have h_eq_fg : f = g := eq_of_monic_of_associated hmf hmg h_assoc
      exact hne h_eq_fg
  set d := EuclideanDomain.gcd f g with hd
  have hd_dvd_f : d ∣ f := EuclideanDomain.gcd_dvd_left _ _
  have hd_dvd_g : d ∣ g := EuclideanDomain.gcd_dvd_right _ _
  have h_cases : IsUnit d ∨ Associated d f := by
    rcases hd_dvd_f with ⟨h, h_eq⟩
    -- h_eq : f = d * h
    have h_cases' : IsUnit d ∨ IsUnit h := hf.isUnit_or_isUnit (a := d) (b := h) (by
      -- need: f = d * h
      simpa [mul_comm] using h_eq)
    rcases h_cases' with (h_unit | h_unit')
    · exact Or.inl h_unit
    · exact Or.inr (by
        refine ⟨h_unit'.unit, ?_⟩
        calc
          d * (h_unit'.unit : Polynomial K) = d * h := by rw [h_unit'.unit_spec]
          _ = f := h_eq.symm)
  rcases h_cases with (h_unit | h_assoc)
  · rcases h_unit.exists_right_inv with ⟨inv, h_inv⟩
    have h_bezout : d = f * EuclideanDomain.gcdA f g + g * EuclideanDomain.gcdB f g :=
      EuclideanDomain.gcd_eq_gcd_ab f g
    have h_goal : (EuclideanDomain.gcdA f g * inv) * f + (EuclideanDomain.gcdB f g * inv) * g = 1 := by
      calc
        (EuclideanDomain.gcdA f g * inv) * f + (EuclideanDomain.gcdB f g * inv) * g
            = f * (EuclideanDomain.gcdA f g * inv) + g * (EuclideanDomain.gcdB f g * inv) := by ring
        _ = (f * EuclideanDomain.gcdA f g + g * EuclideanDomain.gcdB f g) * inv := by ring
        _ = d * inv := by rw [h_bezout]
        _ = 1 := h_inv
    refine ⟨EuclideanDomain.gcdA f g * inv, EuclideanDomain.gcdB f g * inv, h_goal⟩
  · have h_f_dvd_g : f ∣ g := by
      have hf_dvd_d : f ∣ d := h_assoc.symm.dvd
      exact dvd_trans hf_dvd_d hd_dvd_g
    exact absurd h_f_dvd_g h_not_dvd

lemma IsCoprime_of_isUnit_left {R : Type*} [CommSemiring R] {u x : R} (hu : IsUnit u) : IsCoprime u x := by
  rcases hu.exists_right_inv with ⟨v, hv⟩
  refine ⟨v, 0, ?_⟩
  simp [mul_comm, hv]

lemma IsCoprime_of_isUnit_right {R : Type*} [CommSemiring R] {x u : R} (hu : IsUnit u) : IsCoprime x u := by
  rcases hu.exists_left_inv with ⟨v, hv⟩
  refine ⟨0, v, ?_⟩
  simp [hv]

lemma list_dvd_prod_of_mem {R : Type*} [CommSemiring R] {a : R} {l : List R} (h : a ∈ l) : a ∣ l.prod := by
  induction l with
  | nil => simp at h
  | cons hd tl ih =>
    simp at h
    rcases h with (rfl | hmem)
    · exact dvd_mul_right _ _
    · exact dvd_trans (ih hmem) (dvd_mul_left _ _)

lemma list_prod_mul_of_two_pos {R : Type*} [CommSemiring R] (l : List R) (i j : ℕ)
    (hi : i < l.length) (hj : j < l.length) (hij : i ≠ j) :
    (l.get ⟨i, hi⟩) * (l.get ⟨j, hj⟩) ∣ l.prod := by
  revert i j hi hj hij
  induction l with
  | nil =>
    intro i j hi hj hij; simp at hi
  | cons hd tl ih =>
    intro i j hi hj hij
    have hi_len : i < (hd :: tl).length := hi
    have hj_len : j < (hd :: tl).length := hj
    simp [List.prod_cons]
    match i, j with
    | 0, 0 => exact (hij rfl).elim
    | 0, j+1 =>
      have hjt : j < tl.length := by
        have : j+1 < (hd :: tl).length := hj
        simp at this; omega
      have h_tl : tl.get ⟨j, hjt⟩ ∣ tl.prod := list_dvd_prod_of_mem (by apply List.get_mem)
      have h_mul : hd * tl.get ⟨j, hjt⟩ ∣ hd * tl.prod := mul_dvd_mul_left hd h_tl
      simpa [List.get] using h_mul
    | i+1, 0 =>
      have hit : i < tl.length := by
        have : i+1 < (hd :: tl).length := hi
        simp at this; omega
      have h_tl : tl.get ⟨i, hit⟩ ∣ tl.prod := list_dvd_prod_of_mem (by apply List.get_mem)
      have h_mul : hd * tl.get ⟨i, hit⟩ ∣ hd * tl.prod := mul_dvd_mul_left hd h_tl
      simpa [List.get, mul_comm] using h_mul
    | i+1, j+1 =>
      have hit : i < tl.length := by
        have : i+1 < (hd :: tl).length := hi
        simp at this; omega
      have hjt : j < tl.length := by
        have : j+1 < (hd :: tl).length := hj
        simp at this; omega
      have hij' : i ≠ j := by
        intro h; apply hij; omega
      have h_tl : (tl.get ⟨i, hit⟩) * (tl.get ⟨j, hjt⟩) ∣ tl.prod := ih i j hit hjt hij'
      have h_mul : tl.get ⟨i, hit⟩ * tl.get ⟨j, hjt⟩ ∣ hd * tl.prod :=
        h_tl.trans (dvd_mul_left _ _)
      simpa [List.get] using h_mul

lemma squarefree_factors_pairwise_coprime (p : ℕ) [Fact (Nat.Prime p)]
    (fp : Polynomial (ZMod p)) (lc : ZMod p) (factors : List (Polynomial (ZMod p) × ℕ))
    (h_correct : FactorZpCorrect fp lc factors) (h_sqfree : Squarefree fp) (h_fp_ne_zero : fp ≠ 0) :
    (Polynomial.C lc :: factors.map (fun (fi, ei) => fi ^ ei)).Pairwise (fun a b => IsCoprime a b) := by
  have hlc_ne_zero : lc ≠ 0 := by
    intro hzero
    apply h_fp_ne_zero
    rw [h_correct.1, hzero, Polynomial.C_0, zero_mul]
  have hlc_unit : IsUnit (Polynomial.C lc) :=
    Polynomial.isUnit_C.mpr (isUnit_iff_ne_zero.mpr hlc_ne_zero)
  let facs_list := factors.map (fun (fi, ei) => fi ^ ei)
  -- Each pair of distinct list entries is coprime
  have h_pairwise : List.Pairwise IsCoprime facs_list := by
    apply List.pairwise_iff_get.mpr
    intro i j h_ij
    have hi_len : i.1 < facs_list.length := i.2
    have hj_len : j.1 < facs_list.length := j.2
    have hi_factors : i.1 < factors.length := by
      simpa [facs_list] using hi_len
    have hj_factors : j.1 < factors.length := by
      simpa [facs_list] using hj_len
    let fi := factors[i.1].1
    let fj := factors[j.1].1
    have hi_mem : (fi, factors[i.1].2) ∈ factors :=
      List.mem_iff_get.mpr ⟨⟨i.1, hi_factors⟩, rfl⟩
    have hj_mem : (fj, factors[j.1].2) ∈ factors :=
      List.mem_iff_get.mpr ⟨⟨j.1, hj_factors⟩, rfl⟩
    rcases h_correct.2 _ hi_mem with ⟨hirr, hmon, h_ei⟩
    rcases h_correct.2 _ hj_mem with ⟨hjrr, hjmon, h_ej⟩
    by_cases h_eq : fi = fj
    · exfalso
      have h_fi_sq_dvd : fi ^ 2 ∣ fp := by
        rw [h_correct.1]
        set l := factors.map (fun (pr : Polynomial (ZMod p) × ℕ) => pr.1 ^ pr.2) with hl
        have hi_len' : i.1 < l.length := by
          simpa [hl] using hi_factors
        have hj_len' : j.1 < l.length := by
          simpa [hl] using hj_factors
        have hij' : i.1 ≠ j.1 := by
          intro h_eq_val; exact Nat.ne_of_lt h_ij h_eq_val
        have h_mul_dvd : (l.get ⟨i.1, hi_len'⟩) * (l.get ⟨j.1, hj_len'⟩) ∣ l.prod :=
          list_prod_mul_of_two_pos l i.1 j.1 hi_len' hj_len' hij'
        have h_elem_i : l.get ⟨i.1, hi_len'⟩ = fi ^ factors[i.1].2 := by
          simp [fi, l]
        have h_elem_j : l.get ⟨j.1, hj_len'⟩ = fi ^ factors[j.1].2 := by
          simp [fi, fj, l, h_eq]
        rw [h_elem_i, h_elem_j] at h_mul_dvd
        have h_sum_eq : fi ^ factors[i.1].2 * fi ^ factors[j.1].2 = fi ^ (factors[i.1].2 + factors[j.1].2) := by ring
        rw [h_sum_eq] at h_mul_dvd
        have h_sq_dvd_sum : fi ^ 2 ∣ fi ^ (factors[i.1].2 + factors[j.1].2) :=
          pow_dvd_pow fi (by
            have h_ei : factors[i.1].2 ≥ 1 := h_ei
            have h_ej : factors[j.1].2 ≥ 1 := h_ej
            omega)
        have h_sq_dvd_prod : fi ^ 2 ∣ l.prod := dvd_trans h_sq_dvd_sum h_mul_dvd
        exact h_sq_dvd_prod.mul_left (Polynomial.C lc)
      have h_fi_mul_sq_dvd : fi * fi ∣ fp := by
        simpa [sq] using h_fi_sq_dvd
      have h_fi_unit : IsUnit fi := h_sqfree fi h_fi_mul_sq_dvd
      exact hirr.1 h_fi_unit
    · have h_cop : IsCoprime fi fj :=
        monic_irreducible_coprime hirr hjrr hmon hjmon h_eq
      have h_cop_pow : IsCoprime (fi ^ factors[i.1].2) (fj ^ factors[j.1].2) :=
        IsCoprime.pow h_cop (m := factors[i.1].2) (n := factors[j.1].2)
      simpa [facs_list, fi, fj] using h_cop_pow
  -- Build pairwise for [C(lc), ...]
  refine List.pairwise_cons.mpr ⟨?_, h_pairwise⟩
  intro x hx
  rcases List.mem_map.mp hx with ⟨p, hmem, rfl⟩
  rcases p with ⟨fi, ei⟩
  apply IsCoprime_of_isUnit_left hlc_unit

/-- Z[x] 因式分解的顶层正确性：
    假设 Zp 因式分解、Hensel 提升、因子重组各自正确，
    则组合结果是完整不可约分解。

    注：hk, hprim, hdeg 是算法前置条件（选素数/确保提升可行），
    框架组合证明不直接使用——它们通过子过程假设间接起作用。-/
theorem factor_ZZ_correct
    (f : Polynomial ℤ) (hf : f ≠ 0) (hprim : f.IsPrimitive)
    -- 选定的素数和 Hensel 提升指数
    {p : ℕ} [Fact (Nat.Prime p)] {k : ℕ} (hk : 0 < k)
    -- p 是合适的素数
    (hgood : Squarefree (Polynomial.map (Int.castRingHom (ZMod p)) f))
    (hdeg : Polynomial.natDegree (Polynomial.map (Int.castRingHom (ZMod p)) f)
            = Polynomial.natDegree f)
    -- 假设子过程存在且正确
    -- 1. Zp 因式分解
    (factor_zp : Polynomial (ZMod p) → ZMod p × List (Polynomial (ZMod p) × ℕ))
    (hfzp : FactorZpCorrect (Polynomial.map (Int.castRingHom (ZMod p)) f)
        (factor_zp (Polynomial.map (Int.castRingHom (ZMod p)) f)).1
        (factor_zp (Polynomial.map (Int.castRingHom (ZMod p)) f)).2)
    -- 2. Hensel 提升
    (hensel : List (Polynomial (ZMod p)) → List (Polynomial (ZMod (p ^ k))))
    (hhensel : ∀ facs_p,
        facs_p ≠ [] →
        Polynomial.map (Int.castRingHom (ZMod p)) f = facs_p.prod →
        (facs_p.Pairwise (fun a b => IsCoprime a b)) →
        HenselCorrect f k facs_p (hensel facs_p))
    -- 3. 因子重组
    (recombine : List (Polynomial (ZMod (p ^ k))) → List (Polynomial ℤ))
    (hrecombine : ∀ facs_pk,
        Polynomial.map (Int.castRingHom (ZMod (p ^ k))) f = facs_pk.prod →
        RecombineCorrect f (recombine facs_pk))
    : ∃ result : List (Polynomial ℤ),
        FactorZZCorrect f result := by
  -- Step 1: f mod p ≠ 0（Squarefree 蕴含非零）
  set fp := Polynomial.map (Int.castRingHom (ZMod p)) f with hfp_def
  have hfp_ne : fp ≠ 0 := by
    intro h; rw [h] at hgood; exact not_squarefree_zero hgood
  -- Step 2: 对 fp 做 Zp 因式分解
  have fzp_ok := hfzp
  -- Step 3: 构造 Hensel 输入列表
  --   facs_p = [C(lc), f₁^e₁, f₂^e₂, ...]
  --   乘积 = C(lc) * ∏(fᵢ^eᵢ) = fp（由 FactorZpCorrect.1 + List.prod_cons）
  have h_facs_prod : fp =
      (Polynomial.C (factor_zp fp).1 ::
       ((factor_zp fp).2.map (fun pr => pr.1 ^ pr.2))).prod := by
    rw [List.prod_cons]; exact fzp_ok.1
  have h_cop : (Polynomial.C (factor_zp fp).1 ::
      ((factor_zp fp).2.map (fun pr => pr.1 ^ pr.2))).Pairwise (fun a b => IsCoprime a b) :=
    squarefree_factors_pairwise_coprime p fp (factor_zp fp).1 (factor_zp fp).2 fzp_ok
      (by
        rw [hfp_def]
        exact hgood)
      hfp_ne
  -- Step 4: Hensel 提升 mod p → mod p^k
  have hensel_ok := hhensel _ (List.cons_ne_nil _ _) h_facs_prod h_cop
  -- Step 5: 因子重组 → FactorZZCorrect
  --   RecombineCorrect f result ≡ FactorZZCorrect f result（定义完全相同）
  exact ⟨recombine (hensel _), hrecombine _ hensel_ok.1⟩
