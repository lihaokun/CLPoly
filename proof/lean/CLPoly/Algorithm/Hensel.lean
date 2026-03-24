/-
  CLPoly/Algorithm/Hensel.lean — L2 Hensel 提升正确性证明

  Phase 3 T3.4: Hensel Lifting (ℤ[x] approach)
  对应 C++: polynomial_factorize_univar.hh:307-597

  证明策略：在 ℤ[x] 中构造 g' = g + C(m)*τ, h' = h + C(m)*σ，
  验证 g'*h' - f = C(m)²*(something)，投射到 ZMod(m²) 得 0。
-/
import CLPoly.Spec
import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.Polynomial.Div
import Mathlib.RingTheory.Nilpotent.Basic

set_option autoImplicit false
set_option maxHeartbeats 1600000

open Polynomial

-- ============================================================
-- 0. 辅助引理
-- ============================================================

/-- ℤ → ZMod m 映射整数 m 到 0 -/
lemma int_cast_m_eq_zero (m : ℕ) : ((m : ℤ) : ZMod m) = 0 := by
  exact_mod_cast ZMod.natCast_self m

/-- 若 map_m p = 0，则 ∃ q, p = C(m) * q（ℤ[x] 中整除提取） -/
lemma exists_C_mul_of_map_eq_zero (m : ℕ) (hm : 0 < m) (p : Polynomial ℤ)
    (hp : Polynomial.map (Int.castRingHom (ZMod m)) p = 0) :
    ∃ q : Polynomial ℤ, p = Polynomial.C ((m : ℤ)) * q := by
  -- Each coefficient of p is divisible by m
  have hdvd : ∀ i, (m : ℤ) ∣ p.coeff i := by
    intro i
    have : (p.coeff i : ZMod m) = 0 := by
      have h := congr_arg (fun q => q.coeff i) hp
      simpa [Polynomial.coeff_map] using h
    exact (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mp this
  -- Construct q with q.coeff i = p.coeff i / m
  refine ⟨Polynomial.ofFinsupp (p.toFinsupp.mapRange (· / (m : ℤ)) (by simp)), ?_⟩
  ext n
  simp only [coeff_C_mul, Polynomial.coeff_ofFinsupp, Finsupp.mapRange_apply]
  -- Goal: p.coeff n = (m : ℤ) * (p.toFinsupp n / (m : ℤ))
  -- p.coeff n = p.toFinsupp n (definitionally equal)
  show p.toFinsupp n = (m : ℤ) * (p.toFinsupp n / (m : ℤ))
  exact (Int.mul_ediv_cancel' (hdvd n)).symm

-- ============================================================
-- 1. hensel_step: 2-factor 单步 (ℤ[x] 方法)
-- ============================================================

theorem hensel_step
    (m : ℕ) (hm : 1 < m)
    (f g h : Polynomial ℤ)
    (hprod : Polynomial.map (Int.castRingHom (ZMod m)) f =
             Polynomial.map (Int.castRingHom (ZMod m)) g *
             Polynomial.map (Int.castRingHom (ZMod m)) h)
    (hcop : IsCoprime (Polynomial.map (Int.castRingHom (ZMod m)) g)
                       (Polynomial.map (Int.castRingHom (ZMod m)) h))
    : ∃ g' h' : Polynomial ℤ,
        Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) f =
        Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) g' *
        Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) h'
        ∧ Polynomial.map (Int.castRingHom (ZMod m)) g' =
          Polynomial.map (Int.castRingHom (ZMod m)) g
        ∧ Polynomial.map (Int.castRingHom (ZMod m)) h' =
          Polynomial.map (Int.castRingHom (ZMod m)) h := by
  -- Step A: f - g*h = C(m) * e_int
  have hfgh : Polynomial.map (Int.castRingHom (ZMod m)) (f - g * h) = 0 := by
    simp only [Polynomial.map_sub, Polynomial.map_mul, hprod, sub_self]
  obtain ⟨e_int, he_int⟩ := exists_C_mul_of_map_eq_zero m (by omega) (f - g * h) hfgh
  -- f = g*h + C(m) * e_int
  have hf_eq : f = g * h + C ((m : ℤ)) * e_int :=
    sub_eq_iff_eq_add'.mp he_int
  -- Step B: Bézout lift to ℤ[x]
  obtain ⟨s_bar, t_bar, hbez⟩ := hcop
  obtain ⟨s, hs⟩ := Polynomial.map_surjective (Int.castRingHom (ZMod m)) ZMod.intCast_surjective s_bar
  obtain ⟨t, ht⟩ := Polynomial.map_surjective (Int.castRingHom (ZMod m)) ZMod.intCast_surjective t_bar
  have hbez_m : Polynomial.map (Int.castRingHom (ZMod m)) (s * g + t * h - 1) = 0 := by
    simp only [Polynomial.map_sub, Polynomial.map_add, Polynomial.map_mul,
               Polynomial.map_one, hs, ht, hbez, sub_self]
  obtain ⟨w, hw⟩ := exists_C_mul_of_map_eq_zero m (by omega) _ hbez_m
  -- s*g + t*h = 1 + C(m)*w
  have hbez_eq : s * g + t * h = 1 + C ((m : ℤ)) * w :=
    sub_eq_iff_eq_add'.mp hw
  -- Step C: Define g', h'
  refine ⟨g + C ((m : ℤ)) * (t * e_int),
          h + C ((m : ℤ)) * (s * e_int), ?_, ?_, ?_⟩
  · -- (H1) map_{m²}(g'*h') = map_{m²}(f)
    -- Key identity: g'*h' - f = C(m)² * (w*e_int + (t*e_int)*(s*e_int))
    suffices h_key : (g + C (↑m : ℤ) * (t * e_int)) * (h + C (↑m : ℤ) * (s * e_int)) - f =
        C (↑m : ℤ) * C (↑m : ℤ) * (w * e_int + t * e_int * (s * e_int)) by
      -- map_{m²} of C(m)*C(m)*anything = 0 since (m : ZMod m²) * (m : ZMod m²) = m² = 0
      have h_zero : Polynomial.map (Int.castRingHom (ZMod (m ^ 2)))
          (C (↑m : ℤ) * C (↑m : ℤ) * (w * e_int + t * e_int * (s * e_int))) = 0 := by
        have hmm : (Int.castRingHom (ZMod (m ^ 2))) ((↑m : ℤ) * (↑m : ℤ)) = 0 := by
          rw [map_mul]
          show ((m : ℤ) : ZMod (m ^ 2)) * ((m : ℤ) : ZMod (m ^ 2)) = 0
          have : ((m * m : ℕ) : ZMod (m ^ 2)) = 0 := by
            rw [show m * m = m ^ 2 from by ring]; exact ZMod.natCast_self _
          exact_mod_cast this
        simp only [Polynomial.map_mul, Polynomial.map_C, Polynomial.map_add, ← C_mul, hmm,
                   map_zero, zero_mul]
      -- g'*h' - f = h_key, so map(g'*h' - f) = h_zero = 0, hence map(g'*h') = map(f)
      have : Polynomial.map (Int.castRingHom (ZMod (m ^ 2)))
          ((g + C (↑m : ℤ) * (t * e_int)) * (h + C (↑m : ℤ) * (s * e_int)) - f) = 0 := by
        rw [h_key]; exact h_zero
      rw [Polynomial.map_sub, Polynomial.map_mul] at this
      exact (sub_eq_zero.mp this).symm
    -- Prove the key identity by algebra
    rw [hf_eq]; linear_combination C ((m : ℤ)) * e_int * hbez_eq
  · -- (H2) map_m(g') = map_m(g)
    simp only [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_C]; simp
  · -- (H3) map_m(h') = map_m(h)
    simp only [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_C]; simp

-- ============================================================
-- 1b. canonical_lift: ZMod m[x] → ℤ[x] with degree control
-- ============================================================

/-- Canonical lift from ZMod m[x] to ℤ[x] via ZMod.val coefficients.
    Key property: natDegree(canonical_lift m p) ≤ natDegree(p). -/
noncomputable def canonical_lift (m : ℕ) (p : Polynomial (ZMod m)) : Polynomial ℤ :=
  (Finset.range (p.natDegree + 1)).sum (fun i => C ((ZMod.val (p.coeff i) : ℤ)) * X ^ i)

/-- The canonical lift maps back to the original polynomial under map_m. -/
lemma canonical_lift_map (m : ℕ) (_hm : 1 < m) (p : Polynomial (ZMod m)) :
    Polynomial.map (Int.castRingHom (ZMod m)) (canonical_lift m p) = p := by
  haveI : NeZero m := ⟨by omega⟩
  -- map distributes over Finset.sum
  simp only [canonical_lift, Polynomial.map_sum, Polynomial.map_mul, Polynomial.map_pow,
             Polynomial.map_C, Polynomial.map_X]
  ext i
  rw [Polynomial.finset_sum_coeff]
  by_cases hi : i ≤ p.natDegree
  · -- i in range: the sum picks out exactly the i-th term
    have hi_mem : i ∈ Finset.range (p.natDegree + 1) := Finset.mem_range.mpr (by omega)
    rw [Finset.sum_eq_single i]
    · simp only [coeff_C_mul, coeff_X_pow, ite_true, mul_one]
      -- Goal: (Int.castRingHom (ZMod m)) ↑(p.coeff i).val = p.coeff i
      change ((ZMod.val (p.coeff i) : ℤ) : ZMod m) = p.coeff i
      rw [Int.cast_natCast, ZMod.natCast_zmod_val]
    · intro j _ hji
      simp only [coeff_C_mul, coeff_X_pow, if_neg (Ne.symm hji), mul_zero]
    · intro hi_nmem; exact absurd hi_mem hi_nmem
  · -- i out of range: sum is zero, and p.coeff i = 0 by natDegree bound
    push_neg at hi
    have hcoeff_zero : p.coeff i = 0 :=
      Polynomial.coeff_eq_zero_of_natDegree_lt hi
    rw [hcoeff_zero]
    apply Finset.sum_eq_zero
    intro j hj
    have hji : j ≠ i := by
      intro heq; subst heq
      exact absurd (Finset.mem_range.mp hj) (by omega)
    simp only [coeff_C_mul, coeff_X_pow, if_neg (Ne.symm hji), mul_zero]

/-- The canonical lift has controlled natDegree: natDegree ≤ natDegree of the input. -/
lemma canonical_lift_natDegree_le (m : ℕ) (p : Polynomial (ZMod m)) :
    (canonical_lift m p).natDegree ≤ p.natDegree := by
  unfold canonical_lift
  calc ((Finset.range (p.natDegree + 1)).sum
          (fun i => C ((ZMod.val (p.coeff i) : ℤ)) * X ^ i)).natDegree
      ≤ (Finset.range (p.natDegree + 1)).sup
          (fun i => (C ((ZMod.val (p.coeff i) : ℤ)) * X ^ i).natDegree) :=
        natDegree_sum_le _ _
    _ ≤ p.natDegree := by
        apply Finset.sup_le
        intro i hi
        have hi_lt := Finset.mem_range.mp hi
        calc (C ((ZMod.val (p.coeff i) : ℤ)) * X ^ i).natDegree
            ≤ (C ((ZMod.val (p.coeff i) : ℤ))).natDegree + (X ^ i).natDegree := natDegree_mul_le
          _ ≤ 0 + i := by
              apply Nat.add_le_add
              · exact (natDegree_C _).le
              · exact natDegree_X_pow_le i
          _ = i := Nat.zero_add i
          _ ≤ p.natDegree := by omega

-- ============================================================
-- 1c. hensel_step_with_degree: 2-factor with degree preservation
-- ============================================================

theorem hensel_step_with_degree
    (m : ℕ) (hm : 1 < m)
    (f g h : Polynomial ℤ)
    (hprod : Polynomial.map (Int.castRingHom (ZMod m)) f =
             Polynomial.map (Int.castRingHom (ZMod m)) g *
             Polynomial.map (Int.castRingHom (ZMod m)) h)
    (hcop : IsCoprime (Polynomial.map (Int.castRingHom (ZMod m)) g)
                       (Polynomial.map (Int.castRingHom (ZMod m)) h))
    (hh_monic : Monic (Polynomial.map (Int.castRingHom (ZMod m)) h))
    (hh_deg : h.natDegree = (Polynomial.map (Int.castRingHom (ZMod m)) h).natDegree)
    : ∃ g' h' : Polynomial ℤ,
        -- H1: product mod m²
        Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) f =
        Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) g' *
        Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) h'
        -- H2: g' ≡ g (mod m)
        ∧ Polynomial.map (Int.castRingHom (ZMod m)) g' =
          Polynomial.map (Int.castRingHom (ZMod m)) g
        -- H3: h' ≡ h (mod m)
        ∧ Polynomial.map (Int.castRingHom (ZMod m)) h' =
          Polynomial.map (Int.castRingHom (ZMod m)) h
        -- H4: degree preservation (ZMod(m²) side)
        ∧ (Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) h').natDegree =
          (Polynomial.map (Int.castRingHom (ZMod m)) h).natDegree
        -- H5: degree preservation (Z[x] side, for iteration)
        ∧ h'.natDegree = h.natDegree
        -- H6: leading coefficient preservation
        ∧ h'.leadingCoeff = h.leadingCoeff := by
  -- Notation abbreviations
  set map_m := Polynomial.map (Int.castRingHom (ZMod m)) with hmap_m_def
  set map_m2 := Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) with hmap_m2_def
  -- Step A: extract error — f - g*h = C(m) * e_int
  have hfgh : map_m (f - g * h) = 0 := by
    simp only [map_m, Polynomial.map_sub, Polynomial.map_mul, hprod, sub_self]
  obtain ⟨e_int, he_int⟩ := exists_C_mul_of_map_eq_zero m (by omega) (f - g * h) hfgh
  have hf_eq : f = g * h + C ((m : ℤ)) * e_int :=
    sub_eq_iff_eq_add'.mp he_int
  -- Step B: Bézout coefficients
  obtain ⟨s_bar, t_bar, hbez⟩ := hcop
  obtain ⟨s, hs⟩ := Polynomial.map_surjective (Int.castRingHom (ZMod m)) ZMod.intCast_surjective s_bar
  obtain ⟨t, ht⟩ := Polynomial.map_surjective (Int.castRingHom (ZMod m)) ZMod.intCast_surjective t_bar
  have hbez_m : map_m (s * g + t * h - 1) = 0 := by
    simp only [map_m, Polynomial.map_sub, Polynomial.map_add, Polynomial.map_mul,
               Polynomial.map_one, hs, ht, hbez, sub_self]
  obtain ⟨w, hw⟩ := exists_C_mul_of_map_eq_zero m (by omega) _ hbez_m
  have hbez_eq : s * g + t * h = 1 + C ((m : ℤ)) * w :=
    sub_eq_iff_eq_add'.mp hw
  -- Step C: modByMonic in ZMod m[x]
  haveI : Fact (1 < m) := ⟨hm⟩
  haveI : NeZero m := ⟨by omega⟩
  haveI : Nontrivial (ZMod m) := ZMod.nontrivial m
  set se_bar := map_m s * map_m e_int with hse_bar_def
  set σ_bar := se_bar %ₘ map_m h with hσ_bar_def
  set q_bar := se_bar /ₘ map_m h with hq_bar_def
  -- Euclidean division: se_bar = (map_m h) * q_bar + σ_bar
  have hse_div : σ_bar + map_m h * q_bar = se_bar :=
    modByMonic_add_div se_bar hh_monic
  -- τ_bar = (map_m t) * (map_m e_int) + q_bar * (map_m g)
  set τ_bar := map_m t * map_m e_int + q_bar * map_m g with hτ_bar_def
  -- Key Bézout identity in ZMod m[x]:
  -- σ_bar * map_m g + τ_bar * map_m h = map_m e_int
  have hbez_bar : σ_bar * map_m g + τ_bar * map_m h = map_m e_int := by
    -- Rewrite τ_bar and use se_bar = σ_bar + (map_m h) * q_bar
    have hq_eq : map_m h * q_bar = se_bar - σ_bar := by
      have := hse_div  -- σ_bar + map_m h * q_bar = se_bar
      linear_combination this
    -- Expand step by step using calc
    calc σ_bar * map_m g + τ_bar * map_m h
        = σ_bar * map_m g + (map_m t * map_m e_int + q_bar * map_m g) * map_m h := by
          rfl
      _ = σ_bar * map_m g + map_m t * map_m e_int * map_m h +
          q_bar * map_m g * map_m h := by ring
      _ = σ_bar * map_m g + map_m t * map_m e_int * map_m h +
          (map_m h * q_bar) * map_m g := by ring
      _ = σ_bar * map_m g + map_m t * map_m e_int * map_m h +
          (se_bar - σ_bar) * map_m g := by rw [hq_eq]
      _ = map_m t * map_m e_int * map_m h + se_bar * map_m g := by ring
      _ = map_m t * map_m e_int * map_m h +
          (map_m s * map_m e_int) * map_m g := by rfl
      _ = (map_m s * map_m g + map_m t * map_m h) * map_m e_int := by ring
      _ = _ := by
          have hs' : map_m s = s_bar := hs
          have ht' : map_m t = t_bar := ht
          rw [hs', ht', hbez, one_mul]
  -- Step D: lift σ_bar and τ_bar back to ℤ[x] via canonical_lift
  set σ_int := canonical_lift m σ_bar with hσ_int_def
  -- For τ, we can use map_surjective (degree control not needed for g')
  obtain ⟨τ_int, hτ_int⟩ := Polynomial.map_surjective (Int.castRingHom (ZMod m)) ZMod.intCast_surjective τ_bar
  -- σ_int maps back to σ_bar
  have hσ_map : map_m σ_int = σ_bar := canonical_lift_map m hm σ_bar
  -- σ_int has controlled degree
  have hσ_deg_le : σ_int.natDegree ≤ σ_bar.natDegree :=
    canonical_lift_natDegree_le m σ_bar
  -- σ_bar has degree < natDegree(map_m h) by modByMonic
  -- Need: map_m h ≠ 1 (follows from monic + degree consideration... actually just need natDegree > 0)
  -- We'll handle the edge case h=1 separately if needed. For now:
  -- Step E: define g', h'
  set g' := g + C ((m : ℤ)) * τ_int with hg'_def
  set h' := h + C ((m : ℤ)) * σ_int with hh'_def
  refine ⟨g', h', ?_, ?_, ?_, ?_⟩
  · -- (H1) map_{m²}(g'*h') = map_{m²}(f)
    -- Algebra: g'·h' - f = C(m)·(σ_int·g + τ_int·h + C(m)·σ_int·τ_int - e_int)
    -- map_m(σ_int·g + τ_int·h - e_int) = σ_bar·map_m(g) + τ_bar·map_m(h) - map_m(e_int) = 0
    -- So σ_int·g + τ_int·h - e_int = C(m)·d for some d
    -- g'·h' - f = C(m²)·(d + σ_int·τ_int) → map_{m²} = 0
    have hmod_zero : map_m (σ_int * g + τ_int * h - e_int) = 0 := by
      simp only [map_m, Polynomial.map_sub, Polynomial.map_add, Polynomial.map_mul,
                 hσ_map, hτ_int, hbez_bar, sub_self]
    obtain ⟨d, hd⟩ := exists_C_mul_of_map_eq_zero m (by omega) _ hmod_zero
    have hd_eq : σ_int * g + τ_int * h = e_int + C ((m : ℤ)) * d :=
      sub_eq_iff_eq_add'.mp hd
    -- Key identity: g'*h' - f = C(m²) * (d + σ_int * τ_int)
    suffices h_key : g' * h' - f =
        C (↑m : ℤ) * C (↑m : ℤ) * (d + σ_int * τ_int) by
      have h_zero : map_m2 (C (↑m : ℤ) * C (↑m : ℤ) * (d + σ_int * τ_int)) = 0 := by
        have hmm : (Int.castRingHom (ZMod (m ^ 2))) ((↑m : ℤ) * (↑m : ℤ)) = 0 := by
          rw [map_mul]
          show ((m : ℤ) : ZMod (m ^ 2)) * ((m : ℤ) : ZMod (m ^ 2)) = 0
          have : ((m * m : ℕ) : ZMod (m ^ 2)) = 0 := by
            rw [show m * m = m ^ 2 from by ring]; exact ZMod.natCast_self _
          exact_mod_cast this
        simp only [map_m2, Polynomial.map_mul, Polynomial.map_add, Polynomial.map_C,
                   ← C_mul, hmm, map_zero, zero_mul]
      have : map_m2 (g' * h' - f) = 0 := by rw [h_key]; exact h_zero
      -- map_m2 is definitionally Polynomial.map (...)
      simp only [hmap_m2_def, Polynomial.map_sub, Polynomial.map_mul] at this
      exact (sub_eq_zero.mp this).symm
    -- Prove the key identity
    simp only [hg'_def, hh'_def]
    rw [hf_eq]
    -- g' * h' = (g + C m * τ)(h + C m * σ)
    --         = g*h + C(m)*(σ*g + τ*h) + C(m)²*σ*τ
    -- f = g*h + C(m)*e
    -- g'*h' - f = C(m)*(σ*g + τ*h - e) + C(m)²*σ*τ
    --           = C(m)*C(m)*d + C(m)²*σ*τ = C(m)²*(d + σ*τ)
    linear_combination C ((m : ℤ)) * hd_eq
  · -- (H2) map_m(g') = map_m(g)
    show map_m (g + C ((m : ℤ)) * τ_int) = map_m g
    change Polynomial.map (Int.castRingHom (ZMod m)) (g + C ((m : ℤ)) * τ_int) =
           Polynomial.map (Int.castRingHom (ZMod m)) g
    simp only [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_C]
    have : (Int.castRingHom (ZMod m)) ((m : ℤ)) = 0 := int_cast_m_eq_zero m
    rw [this]; simp
  · -- (H3) map_m(h') = map_m(h)
    show map_m (h + C ((m : ℤ)) * σ_int) = map_m h
    change Polynomial.map (Int.castRingHom (ZMod m)) (h + C ((m : ℤ)) * σ_int) =
           Polynomial.map (Int.castRingHom (ZMod m)) h
    simp only [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_C]
    have : (Int.castRingHom (ZMod m)) ((m : ℤ)) = 0 := int_cast_m_eq_zero m
    rw [this]; simp
  · -- (H4) natDegree(map_{m²} h') = natDegree(map_m h)
    -- and (H5) h'.natDegree = h.natDegree — proved together
    -- Case split: if map_m h = 1, then σ_bar = 0 so h' = h trivially
    -- If map_m h ≠ 1, use modByMonic degree bound
    -- First establish that σ_int = 0 or has small degree
    by_cases hh_one : map_m h = 1
    · -- Trivial case: map_m h = 1 means h.natDegree = 0
      -- σ_bar = se_bar %ₘ 1 = 0, so σ_int = canonical_lift m 0
      have hσ_bar_zero : σ_bar = 0 := by
        simp only [hσ_bar_def, hh_one, Polynomial.modByMonic_one]
      have hσ_int_zero : σ_int = 0 := by
        rw [hσ_int_def, hσ_bar_zero]
        simp only [canonical_lift, natDegree_zero, Nat.zero_add, Finset.range_one,
                   Finset.sum_singleton, coeff_zero, ZMod.val_zero, Nat.cast_zero,
                   map_zero, zero_mul]
      have hh'_eq_h : h' = h := by rw [hh'_def, hσ_int_zero, mul_zero, add_zero]
      constructor
      · -- H4
        rw [hh'_eq_h, hh_one]
        simp only [natDegree_one]
        -- Need: natDegree(map_{m²} h) = 0
        -- h.natDegree = 0 from hh_deg + hh_one
        have hh_deg0 : h.natDegree = 0 := by rw [hh_deg, hh_one, natDegree_one]
        -- leadingCoeff(h) maps to nonzero mod m² (since it's 1 mod m)
        have hlc_ne : (Int.castRingHom (ZMod (m ^ 2))) h.leadingCoeff ≠ 0 := by
          rw [Polynomial.leadingCoeff, hh_deg0]
          intro h_abs
          -- m² | h.coeff 0
          have h1 : (m ^ 2 : ℤ) ∣ h.coeff 0 :=
            (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mp h_abs
          -- m | m² | h.coeff 0 → m | h.coeff 0
          have h2 : (m : ℤ) ∣ h.coeff 0 :=
            dvd_trans (dvd_pow_self (m : ℤ) (by omega : 2 ≠ 0)) h1
          -- But h.coeff 0 ≡ 1 (mod m) from hh_monic
          have hmonic := hh_monic
          rw [Polynomial.Monic, Polynomial.leadingCoeff] at hmonic
          rw [hmap_m_def, Polynomial.coeff_map, ← hh_deg, hh_deg0] at hmonic
          -- hmonic : (Int.castRingHom (ZMod m)) (h.coeff 0) = 1
          have h3 : (Int.castRingHom (ZMod m)) (h.coeff 0) = 0 :=
            (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mpr h2
          exact one_ne_zero (hmonic ▸ h3)
        have := (Polynomial.natDegree_map_of_leadingCoeff_ne_zero
            (Int.castRingHom (ZMod (m ^ 2))) hlc_ne).symm
        rw [hh_deg0] at this
        rw [hmap_m2_def]; exact this.symm
      · -- H5 ∧ H6
        exact ⟨by rw [hh'_eq_h], by rw [hh'_eq_h]⟩
    · -- Main case: map_m h ≠ 1, use modByMonic
      have hσ_bar_deg : σ_bar.natDegree < (map_m h).natDegree :=
        natDegree_modByMonic_lt se_bar hh_monic hh_one
      have hσ_deg_lt_h : σ_int.natDegree < h.natDegree := by
        calc σ_int.natDegree ≤ σ_bar.natDegree := hσ_deg_le
          _ < (map_m h).natDegree := hσ_bar_deg
          _ = h.natDegree := hh_deg.symm
      have hCm_σ_deg : (C ((m : ℤ)) * σ_int).natDegree < h.natDegree :=
        calc (C ((m : ℤ)) * σ_int).natDegree
            ≤ σ_int.natDegree := natDegree_C_mul_le _ _
          _ < h.natDegree := hσ_deg_lt_h
      have hh'_deg : h'.natDegree = h.natDegree := by
        show (h + C ((m : ℤ)) * σ_int).natDegree = h.natDegree
        exact natDegree_add_eq_left_of_natDegree_lt hCm_σ_deg
      -- leadingCoeff(h') = leadingCoeff(h)
      have hlc_h' : h'.leadingCoeff = h.leadingCoeff := by
        simp only [Polynomial.leadingCoeff, hh'_deg]
        show (h + C ((m : ℤ)) * σ_int).coeff h.natDegree = h.coeff h.natDegree
        rw [coeff_add]
        have : (C ((m : ℤ)) * σ_int).coeff h.natDegree = 0 :=
          Polynomial.coeff_eq_zero_of_natDegree_lt hCm_σ_deg
        rw [this, add_zero]
      -- leadingCoeff(h) maps to 1 under map_m (by hh_monic)
      have hlc_mod_m : (Int.castRingHom (ZMod m)) h.leadingCoeff = 1 := by
        have hmonic := hh_monic
        rw [Polynomial.Monic, Polynomial.leadingCoeff] at hmonic
        rw [hmap_m_def, Polynomial.coeff_map, ← hh_deg] at hmonic
        rw [Polynomial.leadingCoeff]; exact hmonic
      -- So m ∤ leadingCoeff(h), hence m² ∤ leadingCoeff(h')
      have hlc_ne_zero_m2 : (Int.castRingHom (ZMod (m ^ 2))) h'.leadingCoeff ≠ 0 := by
        rw [hlc_h']
        intro h_abs
        have h1 : (m ^ 2 : ℤ) ∣ h.leadingCoeff :=
          (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mp h_abs
        have h2 : (m : ℤ) ∣ h.leadingCoeff :=
          dvd_trans (dvd_pow_self (m : ℤ) (by omega : 2 ≠ 0)) h1
        have h3 : (Int.castRingHom (ZMod m)) h.leadingCoeff = 0 :=
          (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mpr h2
        exact one_ne_zero (hlc_mod_m ▸ h3)
      constructor
      · -- H4: natDegree(map_{m²} h') = natDegree(map_m h)
        -- natDegree_map_of_leadingCoeff_ne_zero gives:
        -- natDegree(Polynomial.map f h') = natDegree(h')
        have h_natdeg_map := Polynomial.natDegree_map_of_leadingCoeff_ne_zero
            (Int.castRingHom (ZMod (m ^ 2))) hlc_ne_zero_m2
        -- h_natdeg_map : (map (Int.castRingHom (ZMod (m^2))) h').natDegree = h'.natDegree
        -- map_m2 h' is definitionally map (Int.castRingHom (ZMod (m^2))) h'
        -- map_m h is definitionally map (Int.castRingHom (ZMod m)) h
        show (map_m2 h').natDegree = (map_m h).natDegree
        change (Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) h').natDegree =
               (Polynomial.map (Int.castRingHom (ZMod m)) h).natDegree
        rw [h_natdeg_map, hh'_deg, hh_deg]
      · -- H5 ∧ H6
        exact ⟨hh'_deg, hlc_h'⟩

-- ============================================================
-- 2. ker(π)² = 0 in ZMod(m²) — ring level
-- ============================================================

/-- 在 ZMod(m²) 中，若 castHom(a) = 0 且 castHom(b) = 0，则 a*b = 0。
    证明：a = m*k, b = m*l → a*b = m²*k*l = 0。 -/
private lemma zmod_ker_mul_eq_zero (m : ℕ) (hm : 1 < m)
    (a b : ZMod (m ^ 2))
    (ha : (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)) a = 0)
    (hb : (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)) b = 0) :
    a * b = 0 := by
  haveI : NeZero (m ^ 2) := ⟨by positivity⟩
  -- a.val and b.val are divisible by m
  simp only [ZMod.castHom_apply] at ha hb
  rw [ZMod.cast_eq_val, ZMod.natCast_eq_zero_iff] at ha
  rw [ZMod.cast_eq_val, ZMod.natCast_eq_zero_iff] at hb
  obtain ⟨ka, hka⟩ := ha; obtain ⟨kb, hkb⟩ := hb
  -- a = (a.val : ZMod(m²)) = (m * ka : ZMod(m²))
  have ha' : a = ((m * ka : ℕ) : ZMod (m ^ 2)) := by
    rw [← hka, ZMod.natCast_zmod_val]
  have hb' : b = ((m * kb : ℕ) : ZMod (m ^ 2)) := by
    rw [← hkb, ZMod.natCast_zmod_val]
  rw [ha', hb']
  push_cast
  have : (m : ZMod (m ^ 2)) * (m : ZMod (m ^ 2)) = 0 := by
    have : ((m * m : ℕ) : ZMod (m ^ 2)) = 0 := by
      rw [show m * m = m ^ 2 from by ring]; exact ZMod.natCast_self _
    exact_mod_cast this
  have key : (m : ZMod (m ^ 2)) * (ka : ZMod (m ^ 2)) * ((m : ZMod (m ^ 2)) * (kb : ZMod (m ^ 2))) =
    ((m : ZMod (m ^ 2)) * (m : ZMod (m ^ 2))) * ((ka : ZMod (m ^ 2)) * (kb : ZMod (m ^ 2))) := by ring
  rw [key, this, zero_mul]

-- ============================================================
-- 3. IsCoprime 传播 + 迭代
-- ============================================================

/-- 若 map π (p * q) = 0 当 p, q 的系数都在 ker(π) 中 -/
private lemma poly_ker_mul_eq_zero (m : ℕ) (hm : 1 < m)
    (p q : Polynomial (ZMod (m ^ 2)))
    (hp : Polynomial.map (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)) p = 0)
    (hq : Polynomial.map (ZMod.castHom (dvd_pow_self m (by omega : 2 ≠ 0)) (ZMod m)) q = 0) :
    p * q = 0 := by
  ext n
  simp only [Polynomial.coeff_mul, Polynomial.coeff_zero]
  apply Finset.sum_eq_zero
  intro ⟨i, j⟩ _
  apply zmod_ker_mul_eq_zero m hm
  · have := congr_arg (fun r => r.coeff i) hp; simpa [Polynomial.coeff_map] using this
  · have := congr_arg (fun r => r.coeff j) hq; simpa [Polynomial.coeff_map] using this

/-- IsCoprime 从 mod m 传播到 mod m² -/
private lemma isCoprime_lift_sq (m : ℕ) (hm : 1 < m) (g h : Polynomial ℤ)
    (hcop : IsCoprime (Polynomial.map (Int.castRingHom (ZMod m)) g)
                       (Polynomial.map (Int.castRingHom (ZMod m)) h))
    : IsCoprime (Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) g)
                 (Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) h) := by
  obtain ⟨s_bar, t_bar, hbez⟩ := hcop
  obtain ⟨s, hs⟩ := Polynomial.map_surjective (Int.castRingHom (ZMod m)) ZMod.intCast_surjective s_bar
  obtain ⟨t, ht⟩ := Polynomial.map_surjective (Int.castRingHom (ZMod m)) ZMod.intCast_surjective t_bar
  let φ2 := Int.castRingHom (ZMod (m ^ 2))
  have hdvd : (m : ℕ) ∣ m ^ 2 := dvd_pow_self m (by omega : 2 ≠ 0)
  let π := ZMod.castHom hdvd (ZMod m)
  -- bez = s*g + t*h in ZMod(m²)[x], maps to 1 under π
  set bez := Polynomial.map φ2 s * Polynomial.map φ2 g +
             Polynomial.map φ2 t * Polynomial.map φ2 h
  have hbez_mod : Polynomial.map π bez = 1 := by
    simp only [bez, Polynomial.map_add, Polynomial.map_mul, Polynomial.map_map]
    have : π.comp φ2 = Int.castRingHom (ZMod m) := by
      ext x; simp [π, φ2, ZMod.castHom_apply, ZMod.cast_intCast hdvd]
    rw [this, hs, ht]; exact hbez
  -- (bez - 1)² = 0 (nilpotent)
  have hnil : (bez - 1) ^ 2 = 0 := by
    have h0 : Polynomial.map π (bez - 1) = 0 := by
      rw [Polynomial.map_sub, Polynomial.map_one, hbez_mod, sub_self]
    have hmul := poly_ker_mul_eq_zero m hm _ _ h0 h0
    rwa [← sq] at hmul
  -- bez = 1 + nilpotent → unit
  have hunit : IsUnit bez := by
    have : IsNilpotent (bez - 1) := ⟨2, hnil⟩
    rw [show bez = 1 + (bez - 1) from by ring]
    exact this.isUnit_one_add
  -- bez unit → divide to get exact Bézout
  obtain ⟨u, hu⟩ := hunit
  -- u⁻¹ * bez = 1, so (u⁻¹ * s̃) * ḡ + (u⁻¹ * t̃) * h̄ = 1
  have : ↑u⁻¹ * bez = 1 := by rw [← hu]; simp [Units.inv_mul]
  exact ⟨↑u⁻¹ * Polynomial.map φ2 s, ↑u⁻¹ * Polynomial.map φ2 t, by
    rw [show ↑u⁻¹ * Polynomial.map φ2 s * Polynomial.map φ2 g +
            ↑u⁻¹ * Polynomial.map φ2 t * Polynomial.map φ2 h =
            ↑u⁻¹ * (Polynomial.map φ2 s * Polynomial.map φ2 g +
                     Polynomial.map φ2 t * Polynomial.map φ2 h) from by ring]
    rw [show Polynomial.map φ2 s * Polynomial.map φ2 g +
            Polynomial.map φ2 t * Polynomial.map φ2 h = bez from rfl]
    exact this⟩

-- ============================================================
-- 4. hensel_two_factor (迭代) — 对 k 归纳
-- ============================================================

theorem hensel_two_factor
    (p : ℕ) (hp : Nat.Prime p) (k : ℕ) (hk : 0 < k)
    (f g h : Polynomial ℤ)
    (hprod : Polynomial.map (Int.castRingHom (ZMod p)) f =
             Polynomial.map (Int.castRingHom (ZMod p)) g *
             Polynomial.map (Int.castRingHom (ZMod p)) h)
    (hcop : IsCoprime (Polynomial.map (Int.castRingHom (ZMod p)) g)
                       (Polynomial.map (Int.castRingHom (ZMod p)) h))
    : ∃ g' h' : Polynomial ℤ,
        Polynomial.map (Int.castRingHom (ZMod (p ^ k))) f =
        Polynomial.map (Int.castRingHom (ZMod (p ^ k))) g' *
        Polynomial.map (Int.castRingHom (ZMod (p ^ k))) h'
        ∧ Polynomial.map (Int.castRingHom (ZMod p)) g' =
          Polynomial.map (Int.castRingHom (ZMod p)) g
        ∧ Polynomial.map (Int.castRingHom (ZMod p)) h' =
          Polynomial.map (Int.castRingHom (ZMod p)) h := by
  -- Induction: mod p → mod p² → mod p⁴ → ... → mod p^{2^n} → project to p^k
  -- For simplicity, use linear induction: mod p^j → mod p^{j+1} (via hensel_step with m = p^j)
  -- Actually, hensel_step goes m → m², so we double each time.
  -- But we can also go linearly using hensel_step with m = p^j, getting p^{2j} which may overshoot.
  -- Simplest: induction on k directly.
  induction k with
  | zero => omega
  | succ k ih =>
    by_cases hk' : k = 0
    · -- k = 0: p^1 = p. Take g' = g, h' = h.
      subst hk'
      have : p ^ (0 + 1) = p := by simp
      refine ⟨g, h, ?_, rfl, rfl⟩
      rw [show (0 + 1) = 1 from rfl, pow_one]
      exact hprod
    · -- k ≥ 1: use IH to get mod p^k, then hensel_step to get mod p^{2k} ⊇ mod p^{k+1}
      have hk_pos : 0 < k := Nat.pos_of_ne_zero hk'
      obtain ⟨g_k, h_k, hprod_k, hg_k, hh_k⟩ := ih hk_pos
      -- IsCoprime in ZMod(p^k): propagate from mod p
      -- First, IsCoprime in ZMod p is given. We need it in ZMod(p^k).
      -- For now, we prove hensel_step mod p^k → mod p^{2k} and project.
      -- hensel_step needs IsCoprime in ZMod(p^k).
      -- This requires iterative propagation which is complex.
      -- Alternative: directly prove the result using the ℤ[x] identities.
          -- Use hensel_step with m = p^k (requires IsCoprime in ZMod(p^k))
      -- Build IsCoprime in ZMod(p^k) by iterating isCoprime_lift_sq
      have hcop_k : IsCoprime (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) g_k)
                               (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) h_k) := by
        -- IsCoprime(map_p g_k, map_p h_k) → IsCoprime(map_{p^k} g_k, map_{p^k} h_k)
        -- Use: map_p g_k = map_p g (from hg_k), IsCoprime(map_p g, map_p h) (from hcop)
        -- Lift Bézout from ZMod p to ZMod(p^k) via the same 1+nilpotent argument
        -- (p is nilpotent in ZMod(p^k): p^k = 0)
        have hcop_pk : IsCoprime (Polynomial.map (Int.castRingHom (ZMod p)) g_k)
                                  (Polynomial.map (Int.castRingHom (ZMod p)) h_k) := by
          rw [hg_k, hh_k]; exact hcop
        -- Lift from mod p to mod p^k
        -- Use: the Bézout sum maps to 1 mod p, hence = 1 + ε where ε is nilpotent in ZMod(p^k)
        -- p^k = 0 in ZMod(p^k) → p nilpotent → (ε coefficients are p-multiples) → ε nilpotent
        -- For k ≤ 1 this is trivial. For k > 1 we use the general nilpotent argument.
        obtain ⟨s_bar, t_bar, hbez_p⟩ := hcop_pk
        obtain ⟨s0, hs0⟩ := Polynomial.map_surjective (Int.castRingHom (ZMod p)) ZMod.intCast_surjective s_bar
        obtain ⟨t0, ht0⟩ := Polynomial.map_surjective (Int.castRingHom (ZMod p)) ZMod.intCast_surjective t_bar
        let φk := Int.castRingHom (ZMod (p ^ k))
        set bez_k := Polynomial.map φk s0 * Polynomial.map φk g_k +
                     Polynomial.map φk t0 * Polynomial.map φk h_k
        -- bez_k maps to 1 mod p
        have hdvd_pk_p : p ∣ p ^ k := dvd_pow_self p (by omega : k ≠ 0)
        have hbez_k_mod : Polynomial.map (ZMod.castHom hdvd_pk_p (ZMod p)) bez_k = 1 := by
          simp only [bez_k, Polynomial.map_add, Polynomial.map_mul, Polynomial.map_map]
          have : (ZMod.castHom hdvd_pk_p (ZMod p)).comp φk = Int.castRingHom (ZMod p) := by
            ext x; simp [ZMod.castHom_apply, ZMod.cast_intCast hdvd_pk_p]
          rw [this, hs0, ht0]; exact hbez_p
        -- bez_k - 1 is nilpotent: (bez_k - 1)^k = 0
        -- Because each coeff of bez_k - 1 maps to 0 mod p, hence is p * something in ZMod(p^k).
        -- (p * x)^k = p^k * x^k = 0 in ZMod(p^k). So each coeff^k = 0.
        -- Polynomial nilpotent when all coeffs nilpotent (can be shown by degree bound).
        -- For simplicity, use IsNilpotent directly:
        have hbk0 : Polynomial.map (ZMod.castHom hdvd_pk_p (ZMod p)) (bez_k - 1) = 0 := by
          rw [Polynomial.map_sub, Polynomial.map_one, hbez_k_mod, sub_self]
        -- All coefficients of (bez_k - 1) are nilpotent in ZMod(p^k)
        -- (bez_k - 1) itself is nilpotent (polynomial over commutative ring, all coeffs nilpotent)
        haveI : NeZero (p ^ k) := ⟨pow_ne_zero k hp.ne_zero⟩
        have hunit : IsUnit bez_k := by
          rw [show bez_k = 1 + (bez_k - 1) from by ring]
          apply IsNilpotent.isUnit_one_add
          -- (bez_k - 1) is nilpotent: use Polynomial.isNilpotent_iff or degree argument
          rw [Polynomial.isNilpotent_iff]
          intro i
          have : (ZMod.castHom hdvd_pk_p (ZMod p)) ((bez_k - 1).coeff i) = 0 := by
            have h := congr_arg (fun r => r.coeff i) hbk0
            simp only [Polynomial.coeff_map, Polynomial.coeff_zero] at h; exact h
          -- (bez_k - 1).coeff i maps to 0 mod p → it's a multiple of p → nilpotent in ZMod(p^k)
          rw [ZMod.castHom_apply, ZMod.cast_eq_val,
              ZMod.natCast_eq_zero_iff] at this
          obtain ⟨c, hc⟩ := this
          refine ⟨k, ?_⟩
          have hval : (bez_k - 1).coeff i = ((p * c : ℕ) : ZMod (p ^ k)) := by
            rw [← hc, ZMod.natCast_zmod_val]
          rw [hval]; push_cast; rw [mul_pow]
          have : (p : ZMod (p ^ k)) ^ k = 0 := by exact_mod_cast ZMod.natCast_self (p ^ k)
          rw [this, zero_mul]
        obtain ⟨u, hu⟩ := hunit
        exact ⟨↑u⁻¹ * Polynomial.map φk s0, ↑u⁻¹ * Polynomial.map φk t0, by
          rw [show ↑u⁻¹ * Polynomial.map φk s0 * Polynomial.map φk g_k +
                  ↑u⁻¹ * Polynomial.map φk t0 * Polynomial.map φk h_k =
                  ↑u⁻¹ * bez_k from by ring, ← hu]; simp [Units.inv_mul]⟩
      have hpk_gt : 1 < p ^ k := by
        exact Nat.one_lt_pow hk' hp.one_lt
      obtain ⟨g', h', hprod', hg', hh'⟩ :=
        hensel_step (p ^ k) hpk_gt _ g_k h_k hprod_k hcop_k
      -- hprod' : mod (p^k)² = mod p^{2k}
      -- Need: mod p^{k+1}. Since k+1 ≤ 2k (when k ≥ 1), project via castHom.
      have h_le : k + 1 ≤ 2 * k := by omega
      have hdvd_pk : p ^ (k + 1) ∣ p ^ (2 * k) := by
        exact Nat.pow_dvd_pow p h_le
      have hdvd_pk2 : p ^ (2 * k) = (p ^ k) ^ 2 := by ring
      refine ⟨g', h', ?_, ?_, ?_⟩
      · -- Product mod p^{k+1}: project from mod p^{2k}
        have hcast : ∀ q : Polynomial ℤ,
            Polynomial.map (Int.castRingHom (ZMod (p ^ (k + 1)))) q =
            Polynomial.map (ZMod.castHom hdvd_pk (ZMod (p ^ (k + 1))))
              (Polynomial.map (Int.castRingHom (ZMod (p ^ (2 * k)))) q) := by
          intro q; rw [Polynomial.map_map]
          congr 1; ext x
          simp [ZMod.castHom_apply, ZMod.cast_intCast hdvd_pk]
        rw [hcast f, hcast g', hcast h', ← Polynomial.map_mul]
        congr 1; rw [hdvd_pk2]; exact hprod'
      · -- mod p preservation: g' ≡ g_k ≡ g (mod p)
        have : Polynomial.map (Int.castRingHom (ZMod p)) g' =
               Polynomial.map (Int.castRingHom (ZMod p)) g_k := by
          -- map_p g' = map_p g_k because g' ≡ g_k (mod p^k) and p | p^k
          have hcast : ∀ q : Polynomial ℤ,
              Polynomial.map (Int.castRingHom (ZMod p)) q =
              Polynomial.map (ZMod.castHom (dvd_pow_self p (by omega : k ≠ 0)) (ZMod p))
                (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) q) := by
            intro q; rw [Polynomial.map_map]; congr 1; ext x
            simp [ZMod.castHom_apply, ZMod.cast_intCast (dvd_pow_self p (by omega : k ≠ 0))]
          rw [hcast g', hcast g_k, hg']
        rw [this, hg_k]
      · -- mod p preservation: h'
        have : Polynomial.map (Int.castRingHom (ZMod p)) h' =
               Polynomial.map (Int.castRingHom (ZMod p)) h_k := by
          have hcast : ∀ q : Polynomial ℤ,
              Polynomial.map (Int.castRingHom (ZMod p)) q =
              Polynomial.map (ZMod.castHom (dvd_pow_self p (by omega : k ≠ 0)) (ZMod p))
                (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) q) := by
            intro q; rw [Polynomial.map_map]; congr 1; ext x
            simp [ZMod.castHom_apply, ZMod.cast_intCast (dvd_pow_self p (by omega : k ≠ 0))]
          rw [hcast h', hcast h_k, hh']
        rw [this, hh_k]

-- ============================================================
-- 4b. IsCoprime from mod p to mod p^k (extracted for reuse)
-- ============================================================

/-- IsCoprime lifts from mod p to mod p^k via nilpotent argument. -/
theorem isCoprime_lift_pk (p : ℕ) (hp : Nat.Prime p) (k : ℕ) (hk : k ≠ 0)
    (g h : Polynomial ℤ)
    (hcop : IsCoprime (Polynomial.map (Int.castRingHom (ZMod p)) g)
                       (Polynomial.map (Int.castRingHom (ZMod p)) h))
    : IsCoprime (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) g)
                 (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) h) := by
  obtain ⟨s_bar, t_bar, hbez_p⟩ := hcop
  obtain ⟨s0, hs0⟩ := Polynomial.map_surjective (Int.castRingHom (ZMod p)) ZMod.intCast_surjective s_bar
  obtain ⟨t0, ht0⟩ := Polynomial.map_surjective (Int.castRingHom (ZMod p)) ZMod.intCast_surjective t_bar
  let φk := Int.castRingHom (ZMod (p ^ k))
  set bez_k := Polynomial.map φk s0 * Polynomial.map φk g +
               Polynomial.map φk t0 * Polynomial.map φk h
  have hdvd_pk_p : p ∣ p ^ k := dvd_pow_self p hk
  have hbez_k_mod : Polynomial.map (ZMod.castHom hdvd_pk_p (ZMod p)) bez_k = 1 := by
    simp only [bez_k, Polynomial.map_add, Polynomial.map_mul, Polynomial.map_map]
    have : (ZMod.castHom hdvd_pk_p (ZMod p)).comp φk = Int.castRingHom (ZMod p) := by
      ext x; simp [ZMod.castHom_apply, ZMod.cast_intCast hdvd_pk_p]
    rw [this, hs0, ht0]; exact hbez_p
  haveI : NeZero (p ^ k) := ⟨pow_ne_zero k hp.ne_zero⟩
  have hunit : IsUnit bez_k := by
    rw [show bez_k = 1 + (bez_k - 1) from by ring]
    apply IsNilpotent.isUnit_one_add
    rw [Polynomial.isNilpotent_iff]
    intro i
    have hbk0 : Polynomial.map (ZMod.castHom hdvd_pk_p (ZMod p)) (bez_k - 1) = 0 := by
      rw [Polynomial.map_sub, Polynomial.map_one, hbez_k_mod, sub_self]
    have : (ZMod.castHom hdvd_pk_p (ZMod p)) ((bez_k - 1).coeff i) = 0 := by
      have h := congr_arg (fun r => r.coeff i) hbk0
      simpa only [Polynomial.coeff_map, Polynomial.coeff_zero] using h
    rw [ZMod.castHom_apply, ZMod.cast_eq_val, ZMod.natCast_eq_zero_iff] at this
    obtain ⟨c, hc⟩ := this
    refine ⟨k, ?_⟩
    rw [show (bez_k - 1).coeff i = ((p * c : ℕ) : ZMod (p ^ k)) from by
      rw [← hc, ZMod.natCast_zmod_val]]
    push_cast; rw [mul_pow]
    have : (p : ZMod (p ^ k)) ^ k = 0 := by exact_mod_cast ZMod.natCast_self (p ^ k)
    rw [this, zero_mul]
  obtain ⟨u, hu⟩ := hunit
  exact ⟨↑u⁻¹ * Polynomial.map φk s0, ↑u⁻¹ * Polynomial.map φk t0, by
    rw [show ↑u⁻¹ * Polynomial.map φk s0 * Polynomial.map φk g +
            ↑u⁻¹ * Polynomial.map φk t0 * Polynomial.map φk h =
            ↑u⁻¹ * bez_k from by ring, ← hu]; simp [Units.inv_mul]⟩

-- ============================================================
-- 4c. hensel_two_factor_deg: 2-factor with degree preservation
-- ============================================================

/-- 2-factor Hensel lifting with degree preservation.
    Requires h monic in ℤ[x]. Outputs include Monic h' and natDegree preservation. -/
theorem hensel_two_factor_deg
    (p : ℕ) (hp : Nat.Prime p) (k : ℕ) (hk : 0 < k)
    (f g h : Polynomial ℤ)
    (hprod : Polynomial.map (Int.castRingHom (ZMod p)) f =
             Polynomial.map (Int.castRingHom (ZMod p)) g *
             Polynomial.map (Int.castRingHom (ZMod p)) h)
    (hcop : IsCoprime (Polynomial.map (Int.castRingHom (ZMod p)) g)
                       (Polynomial.map (Int.castRingHom (ZMod p)) h))
    (hh_monic : Monic h)
    : ∃ g' h' : Polynomial ℤ,
        Polynomial.map (Int.castRingHom (ZMod (p ^ k))) f =
        Polynomial.map (Int.castRingHom (ZMod (p ^ k))) g' *
        Polynomial.map (Int.castRingHom (ZMod (p ^ k))) h'
        ∧ Polynomial.map (Int.castRingHom (ZMod p)) g' =
          Polynomial.map (Int.castRingHom (ZMod p)) g
        ∧ Polynomial.map (Int.castRingHom (ZMod p)) h' =
          Polynomial.map (Int.castRingHom (ZMod p)) h
        ∧ Monic h'
        ∧ h'.natDegree = h.natDegree := by
  induction k with
  | zero => omega
  | succ k ih =>
    by_cases hk' : k = 0
    · subst hk'; exact ⟨g, h, by rwa [pow_one], rfl, rfl, hh_monic, rfl⟩
    · have hk_pos : 0 < k := Nat.pos_of_ne_zero hk'
      obtain ⟨g_k, h_k, hprod_k, hg_k, hh_k, hh_k_monic, hh_k_deg⟩ := ih hk_pos
      -- Coprimality in ZMod(p^k)
      have hcop_k : IsCoprime (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) g_k)
                               (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) h_k) :=
        isCoprime_lift_pk p hp k hk' g_k h_k (by rw [hg_k, hh_k]; exact hcop)
      -- Monic and degree hypotheses for hensel_step_with_degree
      have hpk_gt : 1 < p ^ k := Nat.one_lt_pow hk' hp.one_lt
      have hh_monic_mk : Monic (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) h_k) :=
        hh_k_monic.map _
      haveI : Fact (1 < p ^ k) := ⟨hpk_gt⟩
      haveI : Nontrivial (ZMod (p ^ k)) := ZMod.nontrivial (p ^ k)
      have hh_deg_mk : h_k.natDegree =
          (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) h_k).natDegree :=
        (Polynomial.natDegree_map_of_leadingCoeff_ne_zero _ (by
          rw [hh_k_monic.leadingCoeff, map_one]; exact one_ne_zero)).symm
      -- Apply hensel_step_with_degree
      obtain ⟨g', h', hprod', hg', hh', hh'_deg_m2, hh'_deg, hh'_lc⟩ :=
        hensel_step_with_degree (p ^ k) hpk_gt _ g_k h_k hprod_k hcop_k hh_monic_mk hh_deg_mk
      -- h' is monic (from H6 + h_k monic)
      have hh'_monic : Monic h' := by
        rwa [Polynomial.Monic, hh'_lc, ← Polynomial.Monic]
      -- Project from p^{2k} to p^{k+1}
      have h_le : k + 1 ≤ 2 * k := by omega
      have hdvd_pk : p ^ (k + 1) ∣ p ^ (2 * k) := Nat.pow_dvd_pow p h_le
      have hdvd_pk2 : p ^ (2 * k) = (p ^ k) ^ 2 := by ring
      refine ⟨g', h', ?_, ?_, ?_, hh'_monic, by omega⟩
      · -- Product mod p^{k+1}: project from mod p^{2k}
        have hcast : ∀ q : Polynomial ℤ,
            Polynomial.map (Int.castRingHom (ZMod (p ^ (k + 1)))) q =
            Polynomial.map (ZMod.castHom hdvd_pk (ZMod (p ^ (k + 1))))
              (Polynomial.map (Int.castRingHom (ZMod (p ^ (2 * k)))) q) := by
          intro q; rw [Polynomial.map_map]; congr 1; ext x
          simp [ZMod.castHom_apply, ZMod.cast_intCast hdvd_pk]
        rw [hcast f, hcast g', hcast h', ← Polynomial.map_mul]; congr 1; rw [hdvd_pk2]; exact hprod'
      · -- mod p: g' ≡ g (via g' ≡ g_k mod p^k → g' ≡ g_k mod p → g' ≡ g mod p)
        have hcast : ∀ q : Polynomial ℤ,
            Polynomial.map (Int.castRingHom (ZMod p)) q =
            Polynomial.map (ZMod.castHom (dvd_pow_self p hk') (ZMod p))
              (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) q) := by
          intro q; rw [Polynomial.map_map]; congr 1; ext x
          simp [ZMod.castHom_apply, ZMod.cast_intCast (dvd_pow_self p hk')]
        have : Polynomial.map (Int.castRingHom (ZMod p)) g' =
               Polynomial.map (Int.castRingHom (ZMod p)) g_k := by
          rw [hcast g', hcast g_k, hg']
        rw [this, hg_k]
      · -- mod p: h' ≡ h
        have hcast : ∀ q : Polynomial ℤ,
            Polynomial.map (Int.castRingHom (ZMod p)) q =
            Polynomial.map (ZMod.castHom (dvd_pow_self p hk') (ZMod p))
              (Polynomial.map (Int.castRingHom (ZMod (p ^ k))) q) := by
          intro q; rw [Polynomial.map_map]; congr 1; ext x
          simp [ZMod.castHom_apply, ZMod.cast_intCast (dvd_pow_self p hk')]
        have : Polynomial.map (Int.castRingHom (ZMod p)) h' =
               Polynomial.map (Int.castRingHom (ZMod p)) h_k := by
          rw [hcast h', hcast h_k, hh']
        rw [this, hh_k]

-- ============================================================
-- 5. hensel_multifactor: multi-factor Hensel lifting
-- ============================================================

/-- Helper: IsCoprime with a list product (right side). -/
private lemma isCoprime_list_prod_right
    {R : Type*} [CommSemiring R]
    (a : R) (l : List R) (h : ∀ b ∈ l, IsCoprime a b) : IsCoprime a l.prod := by
  induction l with
  | nil => exact isCoprime_one_right
  | cons b rest ih =>
    rw [List.prod_cons]
    exact (h b (.head ..)).mul_right (ih (fun c hc => h c (.tail _ hc)))

/-- Multi-factor Hensel lifting: given pairwise coprime factors of f mod p,
    there exist ℤ[x] lifts whose images mod p^k multiply to map_{p^k}(f)
    and agree with the originals mod p. -/
theorem hensel_multifactor
    (p : ℕ) (hp : Nat.Prime p) (k : ℕ) (hk : 0 < k)
    (f : Polynomial ℤ) (factors : List (Polynomial ℤ))
    (hne : factors ≠ [])
    (hprod : Polynomial.map (Int.castRingHom (ZMod p)) f =
             (factors.map (Polynomial.map (Int.castRingHom (ZMod p)))).prod)
    (hcop : factors.Pairwise (fun a b =>
        IsCoprime (Polynomial.map (Int.castRingHom (ZMod p)) a)
                  (Polynomial.map (Int.castRingHom (ZMod p)) b)))
    : ∃ lifted : List (Polynomial ℤ),
        lifted.length = factors.length
        ∧ Polynomial.map (Int.castRingHom (ZMod (p ^ k))) f =
          (lifted.map (Polynomial.map (Int.castRingHom (ZMod (p ^ k))))).prod
        ∧ List.Forall₂ (fun g h => Polynomial.map (Int.castRingHom (ZMod p)) g =
                                     Polynomial.map (Int.castRingHom (ZMod p)) h)
            factors lifted := by
  induction factors generalizing f with
  | nil => exact absurd rfl hne
  | cons g rest ih =>
    by_cases hrest : rest = []
    · -- Singleton: factors = [g]. Take lifted = [f].
      subst hrest; simp only [List.map_cons, List.map_nil, List.prod_cons, List.prod_nil,
        mul_one] at hprod
      exact ⟨[f], rfl, by simp, .cons hprod.symm .nil⟩
    · -- Split: g :: rest with rest nonempty
      simp only [List.map_cons, List.prod_cons] at hprod
      have hcop1 := (List.pairwise_cons.mp hcop).1
      have hcop2 := (List.pairwise_cons.mp hcop).2
      -- Coprimality with product of rest
      have hcop_gr : IsCoprime
          (Polynomial.map (Int.castRingHom (ZMod p)) g)
          (Polynomial.map (Int.castRingHom (ZMod p)) rest.prod) := by
        rw [Polynomial.map_list_prod]
        exact isCoprime_list_prod_right _ _ (fun b hb => by
          obtain ⟨x, hx_mem, hx_eq⟩ := List.mem_map.mp hb
          rw [← hx_eq]; exact hcop1 x hx_mem)
      -- Product condition for hensel_two_factor
      have hprod_2 : Polynomial.map (Int.castRingHom (ZMod p)) f =
          Polynomial.map (Int.castRingHom (ZMod p)) g *
          Polynomial.map (Int.castRingHom (ZMod p)) rest.prod := by
        rw [Polynomial.map_list_prod]; exact hprod
      -- Apply hensel_two_factor
      obtain ⟨g', h', hprod_k, hg'_p, hh'_p⟩ :=
        hensel_two_factor p hp k hk f g rest.prod hprod_2 hcop_gr
      -- Recurse on rest with h'
      have hprod_rest : Polynomial.map (Int.castRingHom (ZMod p)) h' =
          (rest.map (Polynomial.map (Int.castRingHom (ZMod p)))).prod := by
        rw [hh'_p, Polynomial.map_list_prod]
      obtain ⟨lifted_rest, hlen, hprod_rest_k, hforall⟩ :=
        ih h' hrest hprod_rest hcop2
      -- Combine: g' :: lifted_rest
      refine ⟨g' :: lifted_rest, by simp [hlen], ?_, .cons hg'_p.symm hforall⟩
      simp only [List.map_cons, List.prod_cons]
      rw [← hprod_rest_k]; exact hprod_k
