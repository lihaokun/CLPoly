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
-- 5. hensel_multifactor + hensel_correct
-- ============================================================

-- TODO: hensel_multifactor (induction on list length) + HenselCorrect
