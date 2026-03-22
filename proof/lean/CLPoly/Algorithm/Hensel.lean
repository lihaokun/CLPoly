/-
  CLPoly/Algorithm/Hensel.lean — L2 Hensel 提升正确性证明

  Phase 3 T3.4: Hensel Lifting (ℤ[x] approach)
  对应 C++: polynomial_factorize_univar.hh:307-597

  证明策略：在 ℤ[x] 中构造 g' = g + C(m)*τ, h' = h + C(m)*σ，
  验证 g'*h' - f = C(m²)*(something)，投射到 ZMod(m²) 得 0。
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

/-- 若 map_m p = 0，则 p 的每个系数被 m 整除 -/
private lemma coeff_dvd_of_map_eq_zero (m : ℕ) (p : Polynomial ℤ)
    (hp : Polynomial.map (Int.castRingHom (ZMod m)) p = 0) :
    ∀ i, (m : ℤ) ∣ p.coeff i := by
  intro i
  have h1 : (p.coeff i : ZMod m) = 0 := by
    have h := congr_arg (fun q => q.coeff i) hp
    simpa [Polynomial.coeff_map] using h
  exact (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mp h1

/-- 若 map_m p = 0，则 ∃ q, p = C(↑m) * q（ℤ[x] 中整除提取） -/
private lemma exists_C_mul_of_map_eq_zero (m : ℕ) (hm : 0 < m) (p : Polynomial ℤ)
    (hp : Polynomial.map (Int.castRingHom (ZMod m)) p = 0) :
    ∃ q : Polynomial ℤ, p = Polynomial.C ((m : ℤ)) * q := by
  have hdvd := coeff_dvd_of_map_eq_zero m p hp
  refine ⟨Polynomial.ofFinsupp (p.toFinsupp.mapRange (· / (↑m : ℤ)) (by simp)), ?_⟩
  ext n
  simp only [Polynomial.coeff_ofFinsupp, Finsupp.mapRange_apply, coeff_C_mul]
  exact (Int.ediv_mul_cancel (hdvd n)).symm

/-- C(m) 在 ZMod m 中映射为 0 -/
private lemma map_C_m_eq_zero (m : ℕ) :
    Polynomial.map (Int.castRingHom (ZMod m)) (Polynomial.C (↑m : ℤ)) = 0 := by
  simp [Polynomial.map_C, ZMod.intCast_zmod_eq_zero_iff_dvd]

/-- C(m²) 在 ZMod(m²) 中映射为 0 -/
private lemma map_C_m_sq_eq_zero (m : ℕ) :
    Polynomial.map (Int.castRingHom (ZMod (m ^ 2)))
      (Polynomial.C (↑m * ↑m : ℤ)) = 0 := by
  simp [Polynomial.map_C, ZMod.intCast_zmod_eq_zero_iff_dvd]
  exact ⟨1, by ring⟩

-- ============================================================
-- 1. hensel_step: 2-factor 单步 (ℤ[x] 方法)
-- ============================================================

/-- Hensel 单步：mod m → mod m²。
    g' = g + C(m)*τ, h' = h + C(m)*σ。
    乘积验证：g'*h' - f = C(m²)*(something) → 0 mod m²。 -/
theorem hensel_step
    (m : ℕ) (hm : 1 < m)
    (f g h : Polynomial ℤ)
    (hprod : Polynomial.map (Int.castRingHom (ZMod m)) f =
             Polynomial.map (Int.castRingHom (ZMod m)) g *
             Polynomial.map (Int.castRingHom (ZMod m)) h)
    (hcop : IsCoprime (Polynomial.map (Int.castRingHom (ZMod m)) g)
                       (Polynomial.map (Int.castRingHom (ZMod m)) h))
    : ∃ g' h' : Polynomial ℤ,
        -- (H1) f ≡ g'·h' (mod m²)
        Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) f =
        Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) g' *
        Polynomial.map (Int.castRingHom (ZMod (m ^ 2))) h'
        -- (H2) g' ≡ g (mod m)
        ∧ Polynomial.map (Int.castRingHom (ZMod m)) g' =
          Polynomial.map (Int.castRingHom (ZMod m)) g
        -- (H3) h' ≡ h (mod m)
        ∧ Polynomial.map (Int.castRingHom (ZMod m)) h' =
          Polynomial.map (Int.castRingHom (ZMod m)) h := by
  -- Abbreviate
  let φ := Int.castRingHom (ZMod m)
  let φ2 := Int.castRingHom (ZMod (m ^ 2))
  -- Step A: f - g*h = C(m) * e_int
  have hfgh : Polynomial.map φ (f - g * h) = 0 := by
    simp only [Polynomial.map_sub, Polynomial.map_mul, hprod, sub_self]
  obtain ⟨e_int, he_int⟩ := exists_C_mul_of_map_eq_zero m (by omega) (f - g * h) hfgh
  -- So f = g*h + C(m) * e_int
  have hf_eq : f = g * h + C (↑m : ℤ) * e_int :=
    eq_add_of_sub_eq he_int
  -- Step B: Bézout lift to ℤ[x]
  obtain ⟨s_bar, t_bar, hbez⟩ := hcop
  obtain ⟨s, hs⟩ := Polynomial.map_surjective φ (ZMod.intCast_surjective) s_bar
  obtain ⟨t, ht⟩ := Polynomial.map_surjective φ (ZMod.intCast_surjective) t_bar
  -- s*g + t*h ≡ 1 (mod m)
  have hbez_m : Polynomial.map φ (s * g + t * h - 1) = 0 := by
    simp only [Polynomial.map_sub, Polynomial.map_add, Polynomial.map_mul,
               Polynomial.map_one, hs, ht, hbez, sub_self]
  obtain ⟨w, hw⟩ := exists_C_mul_of_map_eq_zero m (by omega) _ hbez_m
  -- s*g + t*h = 1 + C(m)*w
  have hbez_eq : s * g + t * h = 1 + C (↑m : ℤ) * w :=
    eq_add_of_sub_eq hw
  -- Step C: Bézout decomposition of e_int (mod m)
  -- s*e_int*g + t*e_int*h = e_int * (s*g + t*h) = e_int * (1 + C(m)*w) = e_int + C(m)*w*e_int
  -- So s*e_int*g + t*e_int*h - e_int = C(m)*w*e_int
  -- Define σ = s*e_int (mod m lift), τ = t*e_int (mod m lift)
  -- But for degree control, we need divmod in ZMod m[x]
  -- For now, use σ = s*e_int, τ = t*e_int (no degree control in this version)
  -- Step D: Define g', h'
  let τ_int := t * e_int
  let σ_int := s * e_int
  let g' := g + C (↑m : ℤ) * τ_int
  let h' := h + C (↑m : ℤ) * σ_int
  refine ⟨g', h', ?_, ?_, ?_⟩
  · -- (H1) map_{m²}(g'*h') = map_{m²}(f)
    -- Key: g'*h' - f = C(m²) * (something)
    -- g'*h' = (g + C(m)*τ)(h + C(m)*σ)
    --       = g*h + C(m)*(τ*h + g*σ) + C(m)*C(m)*τ*σ
    -- f = g*h + C(m)*e_int
    -- g'*h' - f = C(m)*(τ*h + g*σ - e_int) + C(m)²*τ*σ
    -- τ*h + g*σ - e_int = (t*e_int)*h + g*(s*e_int) - e_int
    --                    = e_int*(t*h + s*g) - e_int
    --                    = e_int*(1 + C(m)*w) - e_int  [by hbez_eq]
    --                    = C(m)*w*e_int
    -- So g'*h' - f = C(m)*C(m)*w*e_int + C(m)*C(m)*τ*σ = C(m²)*(w*e_int + τ*σ)
    have key : g' * h' - f = C (↑m : ℤ) * C (↑m : ℤ) * (w * e_int + τ_int * σ_int) := by
      simp only [g', h', τ_int, σ_int, hf_eq, hbez_eq]; ring
    -- map_{m²}(C(m)*C(m)*anything) = 0 since C(m)*C(m) = C(m²) and m² ≡ 0
    have : Polynomial.map φ2 (g' * h' - f) = 0 := by
      rw [key, Polynomial.map_mul, Polynomial.map_mul, Polynomial.map_C, Polynomial.map_C]
      simp only [Int.cast_natCast]
      have : (↑m : ZMod (m ^ 2)) * (↑m : ZMod (m ^ 2)) = 0 := by
        rw [← Nat.cast_mul, show m * m = m ^ 2 from by ring, ZMod.natCast_self]
      simp [this]
    -- map_{m²}(g'*h') - map_{m²}(f) = 0 → map_{m²}(g'*h') = map_{m²}(f)
    rwa [Polynomial.map_sub, sub_eq_zero, Polynomial.map_mul] at this
  · -- (H2) map_m(g') = map_m(g)
    show Polynomial.map φ g' = Polynomial.map φ g
    simp only [g', Polynomial.map_add, Polynomial.map_mul, Polynomial.map_C]
    simp [Int.cast_natCast, ZMod.natCast_self]
  · -- (H3) map_m(h') = map_m(h)
    show Polynomial.map φ h' = Polynomial.map φ h
    simp only [h', Polynomial.map_add, Polynomial.map_mul, Polynomial.map_C]
    simp [Int.cast_natCast, ZMod.natCast_self]

-- TODO: hensel_two_factor, hensel_multifactor, hensel_correct
-- These build on hensel_step via induction (structurally simpler)
