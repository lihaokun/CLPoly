/-
  CLPoly/Algorithm/Recombine.lean — L2 因子重组正确性证明

  Phase 4: Recombination (Zassenhaus)
  对应 C++: polynomial_factorize_univar.hh:750-882 __zassenhaus_recombine
-/
import CLPoly.Spec
import CLPoly.Algorithm.Hensel
import Mathlib.Data.ZMod.Basic
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Analysis.Polynomial.MahlerMeasure
import Mathlib.NumberTheory.MahlerMeasure

set_option autoImplicit false
set_option maxHeartbeats 3200000

open Polynomial

-- ============================================================
-- 1. Mignotte bound
-- ============================================================

/-- Helper: M(h) >= 1 for nonzero integer polynomial h -/
private lemma mahlerMeasure_ge_one_of_int_ne_zero (h : Polynomial ℤ) (hh : h ≠ 0) :
    1 ≤ (h.map (Int.castRingHom ℂ)).mahlerMeasure :=
  one_le_mahlerMeasure_of_ne_zero hh

/-- Mignotte bound: g | f -> |g.coeff i| <= C(n, n/2) * L1(f) -/
theorem mignotte_bound (f g : Polynomial ℤ) (hf : f ≠ 0) (hg : g ∣ f) :
    ∀ i, (g.coeff i).natAbs ≤
      Nat.choose f.natDegree (f.natDegree / 2) *
      (Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs) := by
  intro i
  have hφ_inj : Function.Injective (Int.castRingHom ℂ) := Int.cast_injective
  obtain ⟨h, hfgh⟩ := hg
  have hh_ne : h ≠ 0 := right_ne_zero_of_mul (hfgh ▸ hf)
  have hg_ne : g ≠ 0 := left_ne_zero_of_mul (hfgh ▸ hf)
  suffices h_real : (↑((g.coeff i).natAbs) : ℝ) ≤
      ↑(Nat.choose f.natDegree (f.natDegree / 2) *
        (Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs)) by
    exact_mod_cast h_real
  have h_lhs : ‖(g.map (Int.castRingHom ℂ)).coeff i‖ = ↑((g.coeff i).natAbs) := by
    rw [coeff_map, eq_intCast, Complex.norm_intCast]
    rw [Nat.cast_natAbs (α := ℝ), Int.cast_abs]
  rw [← h_lhs]
  have h1 := norm_coeff_le_choose_mul_mahlerMeasure i (g.map (Int.castRingHom ℂ))
  rw [natDegree_map_eq_of_injective hφ_inj] at h1
  have h2 : (g.map (Int.castRingHom ℂ)).mahlerMeasure ≤
      (f.map (Int.castRingHom ℂ)).mahlerMeasure := by
    have hfgh_c : f.map (Int.castRingHom ℂ) =
        g.map (Int.castRingHom ℂ) * h.map (Int.castRingHom ℂ) := by
      rw [← Polynomial.map_mul, hfgh]
    rw [hfgh_c, mahlerMeasure_mul]
    exact le_mul_of_one_le_right (mahlerMeasure_nonneg _)
      (mahlerMeasure_ge_one_of_int_ne_zero h hh_ne)
  have hg_dvd : g ∣ f := ⟨h, hfgh⟩
  have h3 : g.natDegree.choose i ≤ f.natDegree.choose (f.natDegree / 2) :=
    (Nat.choose_le_middle i g.natDegree).trans
      ((Nat.choose_le_choose _ (natDegree_le_of_dvd hg_dvd hf)).trans
        (Nat.choose_le_middle _ f.natDegree))
  have h5 : (f.map (Int.castRingHom ℂ)).sum (fun _ a => ‖a‖) =
      ↑((Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs)) := by
    rw [(f.map (Int.castRingHom ℂ)).sum_over_range (fun _ => by simp),
        natDegree_map_eq_of_injective hφ_inj]
    push_cast; congr 1; ext j
    rw [coeff_map, eq_intCast, Complex.norm_intCast]
    rw [Nat.cast_natAbs (α := ℝ), Int.cast_abs]
  calc ‖(g.map (Int.castRingHom ℂ)).coeff i‖
    ≤ ↑(g.natDegree.choose i) * (g.map (Int.castRingHom ℂ)).mahlerMeasure := h1
    _ ≤ ↑(g.natDegree.choose i) * (f.map (Int.castRingHom ℂ)).mahlerMeasure := by gcongr
    _ ≤ ↑(f.natDegree.choose (f.natDegree / 2)) *
          (f.map (Int.castRingHom ℂ)).mahlerMeasure := by
        apply mul_le_mul_of_nonneg_right (by exact_mod_cast h3) (mahlerMeasure_nonneg _)
    _ ≤ ↑(f.natDegree.choose (f.natDegree / 2)) *
          ((f.map (Int.castRingHom ℂ)).sum fun _ a => ‖a‖) := by
        gcongr; exact mahlerMeasure_le_sum_norm_coeff _
    _ = ↑(f.natDegree.choose (f.natDegree / 2) *
          (Finset.range (f.natDegree + 1)).sum (fun j => (f.coeff j).natAbs)) := by
        rw [h5]; push_cast; ring

-- ============================================================
-- 2. Hensel uniqueness
-- ============================================================

-- Helper: projecting polynomial map equalities from larger to smaller modulus
private lemma map_eq_of_dvd_map_eq (a b : ℕ) (hab : a ∣ b) (f g : Polynomial ℤ)
    (h : Polynomial.map (Int.castRingHom (ZMod b)) f =
         Polynomial.map (Int.castRingHom (ZMod b)) g) :
    Polynomial.map (Int.castRingHom (ZMod a)) f =
    Polynomial.map (Int.castRingHom (ZMod a)) g := by
  have key : Polynomial.map (ZMod.castHom hab (ZMod a))
      (Polynomial.map (Int.castRingHom (ZMod b)) f) =
    Polynomial.map (ZMod.castHom hab (ZMod a))
      (Polynomial.map (Int.castRingHom (ZMod b)) g) := congr_arg _ h
  simp only [Polynomial.map_map] at key
  convert key using 1 <;> congr 1 <;> ext x <;> simp [ZMod.castHom_apply, ZMod.cast_intCast hab]

-- Helper: if b | q with b monic and natDegree q < natDegree b, then q = 0
private lemma dvd_monic_eq_zero_of_natDegree_lt {R : Type*} [CommRing R] [IsDomain R]
    (b q : Polynomial R) (hb_monic : Monic b) (hdvd : b ∣ q)
    (hdeg : natDegree q < natDegree b) : q = 0 := by
  by_contra hq
  obtain ⟨c, hc⟩ := hdvd
  have hc_ne : c ≠ 0 := right_ne_zero_of_mul (hc ▸ hq)
  have : natDegree q = natDegree b + natDegree c := hc ▸ hb_monic.natDegree_mul' hc_ne
  omega

-- Helper: if two monic polynomials differ by C(m)*G with m >= 2,
-- then for all n >= natDegree B1, G.coeff n = 0
private lemma monic_diff_coeff_bound
    (B1 B2 : Polynomial ℤ) (m : ℕ) (hm : 1 < m)
    (G : Polynomial ℤ)
    (hB1m : Monic B1) (hB2m : Monic B2)
    (hdiff : B1 - B2 = C (↑m : ℤ) * G) :
    ∀ n, natDegree B1 ≤ n → G.coeff n = 0 := by
  have hndeg : natDegree B1 = natDegree B2 := by
    by_contra hne
    wlog h12 : natDegree B2 < natDegree B1 with H
    · push_neg at h12
      exact H B2 B1 m hm (-G) hB2m hB1m (by linear_combination -hdiff) (Ne.symm hne)
        (lt_of_le_of_ne h12 hne)
    have hcoeff : (B1 - B2).coeff (natDegree B1) = 1 := by
      simp [coeff_sub, hB1m.leadingCoeff, Polynomial.coeff_eq_zero_of_natDegree_lt h12]
    have hcoeff2 : (B1 - B2).coeff (natDegree B1) = ↑m * G.coeff (natDegree B1) := by
      rw [hdiff, coeff_C_mul]
    rw [hcoeff] at hcoeff2
    have := Int.le_of_dvd one_pos ⟨_, hcoeff2⟩; omega
  have hm_ne : (m : ℤ) ≠ 0 := by omega
  intro n hn
  have hcoeff_diff : (B1 - B2).coeff n = ↑m * G.coeff n := by rw [hdiff, coeff_C_mul]
  suffices hsuff : (B1 - B2).coeff n = 0 by
    rw [hsuff] at hcoeff_diff; exact (mul_eq_zero.mp hcoeff_diff.symm).resolve_left hm_ne
  rcases Nat.eq_or_lt_of_le hn with rfl | hgt
  · simp only [coeff_sub]
    have h1 : B1.coeff (natDegree B1) = 1 := hB1m.leadingCoeff
    have h2 : B2.coeff (natDegree B1) = 1 := by rw [hndeg]; exact hB2m.leadingCoeff
    rw [h1, h2]; ring
  · have h1 : B1.coeff n = 0 := Polynomial.coeff_eq_zero_of_natDegree_lt hgt
    have h2 : B2.coeff n = 0 := Polynomial.coeff_eq_zero_of_natDegree_lt (hndeg ▸ hgt)
    simp [coeff_sub, h1, h2]

-- Helper: map_p(Q) = 0 implies p | Q.coeff n
private lemma coeff_dvd_of_map_eq_zero (p : ℕ) (Q : Polynomial ℤ)
    (h : Polynomial.map (Int.castRingHom (ZMod p)) Q = 0)
    (n : ℕ) : (↑p : ℤ) ∣ Q.coeff n := by
  rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
  have := congr_arg (fun q => q.coeff n) h
  simpa [Polynomial.coeff_map] using this

/-- Hensel uniqueness: if two factorizations agree mod p with coprime + B monic,
    they agree mod p^k. Proof by induction on k. -/
theorem hensel_unique
    (p : ℕ) (hp : Nat.Prime p) (k : ℕ) (hk : 0 < k)
    (F : Polynomial ℤ)
    (A1 B1 A2 B2 : Polynomial ℤ)
    (hprod1 : Polynomial.map (Int.castRingHom (ZMod (p ^ k))) (A1 * B1) =
              Polynomial.map (Int.castRingHom (ZMod (p ^ k))) F)
    (hprod2 : Polynomial.map (Int.castRingHom (ZMod (p ^ k))) (A2 * B2) =
              Polynomial.map (Int.castRingHom (ZMod (p ^ k))) F)
    (hA : Polynomial.map (Int.castRingHom (ZMod p)) A1 =
          Polynomial.map (Int.castRingHom (ZMod p)) A2)
    (hB : Polynomial.map (Int.castRingHom (ZMod p)) B1 =
          Polynomial.map (Int.castRingHom (ZMod p)) B2)
    (hcop : IsCoprime (Polynomial.map (Int.castRingHom (ZMod p)) A1)
                       (Polynomial.map (Int.castRingHom (ZMod p)) B1))
    (hB1_monic : Monic B1) (hB2_monic : Monic B2)
    : Polynomial.map (Int.castRingHom (ZMod (p ^ k))) A1 =
      Polynomial.map (Int.castRingHom (ZMod (p ^ k))) A2
    ∧ Polynomial.map (Int.castRingHom (ZMod (p ^ k))) B1 =
      Polynomial.map (Int.castRingHom (ZMod (p ^ k))) B2 := by
  induction k with
  | zero => omega
  | succ k ih =>
    by_cases hk' : k = 0
    · -- Base: p^1 = p
      subst hk'
      have : p ^ (0 + 1) = p := by norm_num
      rw [this] at hprod1 hprod2 ⊢; exact ⟨hA, hB⟩
    · -- Inductive step: k >= 1
      have hk_pos : 0 < k := Nat.pos_of_ne_zero hk'
      have hpk_pos : 0 < p ^ k := Nat.pos_of_ne_zero (pow_ne_zero k hp.ne_zero)
      have hpk_ne : (↑(p ^ k) : ℤ) ≠ 0 := Int.natCast_ne_zero.mpr (by omega)
      have hdvd : p ^ k ∣ p ^ (k + 1) := Nat.pow_dvd_pow p (by omega)
      -- Apply IH
      obtain ⟨ihA, ihB⟩ := ih hk_pos
        (map_eq_of_dvd_map_eq _ _ hdvd _ _ hprod1)
        (map_eq_of_dvd_map_eq _ _ hdvd _ _ hprod2)
      -- Extract: A1 - A2 = C(p^k)*E, B1 - B2 = C(p^k)*G
      obtain ⟨E, hE⟩ := _root_.exists_C_mul_of_map_eq_zero (p ^ k) hpk_pos (A1 - A2)
        (by rw [Polynomial.map_sub, ihA, sub_self])
      obtain ⟨G, hG⟩ := _root_.exists_C_mul_of_map_eq_zero (p ^ k) hpk_pos (B1 - B2)
        (by rw [Polynomial.map_sub, ihB, sub_self])
      -- A1 = A2 + C(p^k)*E, B1 = B2 + C(p^k)*G
      have hA1 : A1 = A2 + C (↑(p ^ k) : ℤ) * E := eq_add_of_sub_eq' hE
      have hB1eq : B1 = B2 + C (↑(p ^ k) : ℤ) * G := eq_add_of_sub_eq' hG
      -- Products agree mod p^(k+1)
      have hprod_eq : Polynomial.map (Int.castRingHom (ZMod (p ^ (k + 1)))) (A1 * B1) =
          Polynomial.map (Int.castRingHom (ZMod (p ^ (k + 1)))) (A2 * B2) := by
        rw [hprod1, hprod2]
      -- Key: show p | coeff(E*B2 + A2*G, n) for all n
      -- From: A1*B1 - A2*B2 = C(p^k)*(E*B2+A2*G) + C(p^k)^2*(E*G) and p^(k+1) | LHS, p^(k+1) | C(p^k)^2*(...)
      have hEBAG_p : ∀ n, (↑p : ℤ) ∣ (E * B2 + A2 * G).coeff n := by
        intro n
        -- Expand: A1*B1 - A2*B2 = C(p^k)*(E*B2+A2*G) + C(p^k)^2*(E*G)
        have hexp : A1 * B1 - A2 * B2 =
            C (↑(p ^ k) : ℤ) * (E * B2 + A2 * G) +
            C (↑(p ^ k) : ℤ) * (C (↑(p ^ k) : ℤ) * (E * G)) := by
          rw [hA1, hB1eq]; ring
        -- p^(k+1) | coeff(A1*B1 - A2*B2, n)
        have h1 : (↑(p ^ (k + 1)) : ℤ) ∣ (A1 * B1 - A2 * B2).coeff n := by
          have h0 : Polynomial.map (Int.castRingHom (ZMod (p ^ (k + 1)))) (A1 * B1 - A2 * B2) = 0 := by
            simp only [Polynomial.map_sub, Polynomial.map_mul, hprod_eq, sub_self]
          rw [← ZMod.intCast_zmod_eq_zero_iff_dvd]
          have h1 := congr_arg (fun q => q.coeff n) h0
          simpa only [Polynomial.coeff_map, Polynomial.coeff_zero] using h1
        -- Coeff of expansion
        have h2 : (A1 * B1 - A2 * B2).coeff n =
            ↑(p ^ k) * (E * B2 + A2 * G).coeff n +
            ↑(p ^ k) * (↑(p ^ k) * (E * G).coeff n) := by
          have := congr_arg (fun q => q.coeff n) hexp
          simp only [coeff_add, coeff_C_mul] at this; exact this
        -- p^(k+1) | p^(2k) * coeff(E*G, n)
        have h3 : (↑(p ^ (k + 1)) : ℤ) ∣ ↑(p ^ k) * (↑(p ^ k) * (E * G).coeff n) := by
          rw [show (↑(p ^ k) : ℤ) * ((↑(p ^ k) : ℤ) * (E * G).coeff n) =
            ↑(p ^ (k + k)) * (E * G).coeff n from by push_cast; ring]
          exact dvd_mul_of_dvd_left (by exact_mod_cast Nat.pow_dvd_pow p (by omega)) _
        -- p^(k+1) | p^k * coeff(E*B2+A2*G, n)
        have h4 : (↑(p ^ (k + 1)) : ℤ) ∣ ↑(p ^ k) * (E * B2 + A2 * G).coeff n := by
          rw [h2] at h1; exact (dvd_add_left h3).mp h1
        -- Cancel p^k
        rw [show (↑(p ^ (k + 1)) : ℤ) = ↑(p ^ k) * ↑p from by push_cast; ring] at h4
        exact (mul_dvd_mul_iff_left hpk_ne).mp h4
      -- Convert to polynomial-level statement
      have hEBAG_map : Polynomial.map (Int.castRingHom (ZMod p)) (E * B2 + A2 * G) = 0 := by
        ext n; simp only [Polynomial.coeff_map, Polynomial.coeff_zero]
        show ((E * B2 + A2 * G).coeff n : ZMod p) = 0
        exact (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mpr (hEBAG_p n)
      -- In F_p[x]: use coprimality
      haveI : Fact (Nat.Prime p) := ⟨hp⟩
      -- E_bar * B_bar + A_bar * G_bar = 0  (using hA, hB to replace A2->A1, B2->B1)
      have hEBAG_ab : Polynomial.map (Int.castRingHom (ZMod p)) E *
          Polynomial.map (Int.castRingHom (ZMod p)) B1 +
          Polynomial.map (Int.castRingHom (ZMod p)) A1 *
          Polynomial.map (Int.castRingHom (ZMod p)) G = 0 := by
        have := hEBAG_map
        simp only [Polynomial.map_add, Polynomial.map_mul] at this
        rwa [← hB, ← hA] at this
      -- Coprimality: A_bar | E_bar and B_bar | G_bar
      have hA_dvd : Polynomial.map (Int.castRingHom (ZMod p)) A1 ∣
          Polynomial.map (Int.castRingHom (ZMod p)) E := by
        apply hcop.dvd_of_dvd_mul_right
        have : Polynomial.map (Int.castRingHom (ZMod p)) E *
            Polynomial.map (Int.castRingHom (ZMod p)) B1 =
          -(Polynomial.map (Int.castRingHom (ZMod p)) A1 *
            Polynomial.map (Int.castRingHom (ZMod p)) G) :=
          (neg_eq_of_add_eq_zero_left hEBAG_ab).symm
        rw [this]; exact dvd_neg.mpr (dvd_mul_right _ _)
      have hB_dvd : Polynomial.map (Int.castRingHom (ZMod p)) B1 ∣
          Polynomial.map (Int.castRingHom (ZMod p)) G := by
        apply hcop.symm.dvd_of_dvd_mul_left
        have : Polynomial.map (Int.castRingHom (ZMod p)) A1 *
            Polynomial.map (Int.castRingHom (ZMod p)) G =
          -(Polynomial.map (Int.castRingHom (ZMod p)) E *
            Polynomial.map (Int.castRingHom (ZMod p)) B1) := (neg_eq_of_add_eq_zero_right hEBAG_ab).symm
        rw [this]; exact dvd_neg.mpr (dvd_mul_left _ _)
      -- B_bar is monic
      have hBm : Monic (Polynomial.map (Int.castRingHom (ZMod p)) B1) := hB1_monic.map _
      -- Degree bound: for n >= natDegree B1, G.coeff n = 0
      have hG_coeff := monic_diff_coeff_bound B1 B2 (p ^ k)
        (Nat.one_lt_pow hk' hp.one_lt) G hB1_monic hB2_monic hG
      -- natDegree(map_p(B1)) = natDegree(B1) since B1 monic
      have hBn : natDegree (Polynomial.map (Int.castRingHom (ZMod p)) B1) = natDegree B1 :=
        Polynomial.natDegree_map_of_leadingCoeff_ne_zero _ (by rw [hB1_monic.leadingCoeff]; simp)
      -- G_bar = 0
      have hGbar : Polynomial.map (Int.castRingHom (ZMod p)) G = 0 := by
        by_cases hG0 : G = 0
        · simp [hG0]
        · have hlt : natDegree G < natDegree B1 := by
            by_contra h; push_neg at h
            exact absurd (hG_coeff _ h) (Polynomial.leadingCoeff_ne_zero.mpr hG0)
          exact dvd_monic_eq_zero_of_natDegree_lt _ _ hBm hB_dvd (by
            calc natDegree (Polynomial.map (Int.castRingHom (ZMod p)) G)
                ≤ natDegree G := Polynomial.natDegree_map_le
              _ < natDegree B1 := hlt
              _ = _ := hBn.symm)
      -- E_bar = 0 (from E_bar * B_bar + 0 = 0 and B_bar monic hence nonzero)
      have hEbar : Polynomial.map (Int.castRingHom (ZMod p)) E = 0 := by
        have : Polynomial.map (Int.castRingHom (ZMod p)) E *
            Polynomial.map (Int.castRingHom (ZMod p)) B1 = 0 := by
          have := hEBAG_ab; rw [hGbar, mul_zero, add_zero] at this; exact this
        exact (mul_eq_zero.mp this).resolve_right hBm.ne_zero
      -- Lift to p^(k+1): if map_p(Q) = 0 then p^(k+1) | C(p^k)*Q
      have lift : ∀ (D Q : Polynomial ℤ), D = C (↑(p ^ k) : ℤ) * Q →
          Polynomial.map (Int.castRingHom (ZMod p)) Q = 0 →
          Polynomial.map (Int.castRingHom (ZMod (p ^ (k + 1)))) D = 0 := by
        intro D Q hDQ hQ0; rw [hDQ]; ext n
        simp only [Polynomial.coeff_map, coeff_C_mul, Polynomial.coeff_zero]
        -- Goal: (Int.castRingHom (ZMod (p^(k+1)))) (↑(p^k) * Q.coeff n) = 0
        obtain ⟨c, hc⟩ := coeff_dvd_of_map_eq_zero p Q hQ0 n
        rw [hc]
        have h1 : (↑(p ^ k) : ℤ) * (↑p * c) = ↑(p ^ (k + 1)) * c := by push_cast; ring
        rw [h1, map_mul]
        suffices h2 : (Int.castRingHom (ZMod (p ^ (k + 1)))) (↑(p ^ (k + 1)) : ℤ) = 0 by
          rw [h2, zero_mul]
        change ((↑(p ^ (k + 1)) : ℤ) : ZMod (p ^ (k + 1))) = 0
        exact (ZMod.intCast_zmod_eq_zero_iff_dvd _ _).mpr (dvd_refl _)
      constructor
      · have h := lift _ _ hE hEbar
        rwa [Polynomial.map_sub, sub_eq_zero] at h
      · have h := lift _ _ hG hGbar
        rwa [Polynomial.map_sub, sub_eq_zero] at h

-- ============================================================
-- 3. RecombineCorrect
-- ============================================================

theorem recombine_correct
    (f : Polynomial ℤ) (hf : f ≠ 0)
    : ∃ result : List (Polynomial ℤ), RecombineCorrect f result := by
  obtain ⟨factors, hirred, hassoc⟩ := WfDvdMonoid.exists_factors f hf
  refine ⟨factors.toList, ?_, ?_⟩
  · rw [Multiset.prod_toList]; exact hassoc.symm
  · intro g hg; exact hirred g (Multiset.mem_toList.mp hg)
