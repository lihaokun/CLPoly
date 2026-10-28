import Mathlib.RingTheory.Polynomial.Resultant.Basic
import Mathlib.RingTheory.Polynomial.GaussLemma
import Mathlib.Data.ZMod.Basic
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Tactic

set_option autoImplicit false

open Polynomial

namespace CLPoly.VanHoeij

/-- This remains valid over composite-modulus coefficient rings. -/
theorem constant_eq_zero_of_monic_positive_degree_dvd
    {R : Type*} [CommRing R] [Nontrivial R]
    (q : Polynomial R) (c : R) (hq : q.Monic) (hqpos : 0 < q.degree)
    (hdiv : q ∣ C c) : c = 0 := by
  have hrem := (modByMonic_eq_zero_iff_dvd hq).mpr hdiv
  rw [(modByMonic_eq_self_iff hq).mpr (degree_C_le.trans_lt hqpos)] at hrem
  exact C_eq_zero.mp hrem

/-- A positive-degree monic common divisor kills the resultant even over
a ring with zero divisors, in particular modulo a prime power. -/
theorem resultant_eq_zero_of_monic_common_divisor
    {R : Type*} [CommRing R] [Nontrivial R]
    (F G q : Polynomial R) (u v : Nat)
    (hF : F.natDegree ≤ u) (hG : G.natDegree ≤ v)
    (hdegree : u ≠ 0 ∨ v ≠ 0)
    (hq : q.Monic) (hqpos : 0 < q.degree)
    (hqF : q ∣ F) (hqG : q ∣ G) : F.resultant G u v = 0 := by
  obtain ⟨A, B, _, _, hbezout⟩ :=
    exists_mul_add_mul_eq_C_resultant F G hF hG hdegree
  have hdiv : q ∣ C (F.resultant G u v) := by
    rw [← hbezout]
    exact dvd_add (dvd_mul_of_dvd_left hqF A) (dvd_mul_of_dvd_left hqG B)
  exact constant_eq_zero_of_monic_positive_degree_dvd q _ hq hqpos hdiv

/-- Common monic local factors give divisibility by the whole modulus,
not just by its underlying prime. -/
theorem modulus_dvd_resultant_of_common_divisor
    (M : Nat) [NeZero M] (hM : 1 < M)
    (F G : Polynomial ℤ) (q : Polynomial (ZMod M)) (u v : Nat)
    (hF : F.natDegree ≤ u) (hG : G.natDegree ≤ v)
    (hdegree : u ≠ 0 ∨ v ≠ 0)
    (hq : q.Monic) (hqpos : 0 < q.degree)
    (hqF : q ∣ F.map (Int.castRingHom (ZMod M)))
    (hqG : q ∣ G.map (Int.castRingHom (ZMod M))) :
    (M : ℤ) ∣ F.resultant G u v := by
  haveI : Fact (1 < M) := ⟨hM⟩
  have hz := resultant_eq_zero_of_monic_common_divisor
    (F.map (Int.castRingHom (ZMod M))) (G.map (Int.castRingHom (ZMod M)))
    q u v ((natDegree_map_le).trans hF) ((natDegree_map_le).trans hG)
    hdegree hq hqpos hqF hqG
  rw [resultant_map_map] at hz
  exact (ZMod.intCast_zmod_eq_zero_iff_dvd _ M).mp hz

/-- A deliberately loose explicit determinant bound, independent of any
modulus or chosen lattice basis. -/
theorem int_det_abs_le_factorial_mul_pow
    (d : Nat) (A : Matrix (Fin d) (Fin d) ℤ) (B : ℤ)
    (hA : ∀ i j, |A i j| ≤ B) :
    |A.det| ≤ (d.factorial : ℤ) * B ^ d := by
  classical
  rw [Matrix.det_apply']
  calc
    |∑ σ : Equiv.Perm (Fin d), (Equiv.Perm.sign σ : ℤ) * ∏ i, A (σ i) i| ≤
        ∑ σ : Equiv.Perm (Fin d),
          |(Equiv.Perm.sign σ : ℤ) * ∏ i, A (σ i) i| :=
      Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _σ : Equiv.Perm (Fin d), B ^ d := by
      apply Finset.sum_le_sum
      intro σ _
      simp only [abs_mul, Equiv.Perm.sign_abs, one_mul, Finset.abs_prod]
      calc
        ∏ i : Fin d, |A (σ i) i| ≤ ∏ _i : Fin d, B :=
          Finset.prod_le_prod (fun _ _ => abs_nonneg _) (fun i _ => hA (σ i) i)
        _ = B ^ d := by simp
    _ = _ := by simp [Fintype.card_perm]

/-- Apply the determinant bound to the actual Sylvester coefficients.
The free degree parameters also cover degree drops modulo a prime power. -/
theorem int_resultant_abs_le
    (F G : Polynomial ℤ) (u v : Nat) (B : ℤ) (hB : 0 ≤ B)
    (hF : ∀ k, |F.coeff k| ≤ B) (hG : ∀ k, |G.coeff k| ≤ B) :
    |F.resultant G u v| ≤ ((u + v).factorial : ℤ) * B ^ (u + v) := by
  apply int_det_abs_le_factorial_mul_pow (u + v) (F.sylvester G u v) B
  intro i j
  refine Fin.addCases (fun j => ?_) (fun j => ?_) j
  · simp only [sylvester, Matrix.of_apply, Fin.addCases_left]
    split_ifs
    · exact hG _
    · simpa using hB
  · simp only [sylvester, Matrix.of_apply, Fin.addCases_right]
    split_ifs
    · exact hF _
    · simpa using hB

/-- Sufficient precision turns a local common factor into a zero integer
resultant. No factorization output is assumed. -/
theorem int_resultant_eq_zero_of_sufficient_modulus
    (M : Nat) [NeZero M] (hM : 1 < M)
    (F G : Polynomial ℤ) (q : Polynomial (ZMod M)) (u v : Nat)
    (hF : F.natDegree ≤ u) (hG : G.natDegree ≤ v)
    (hdegree : u ≠ 0 ∨ v ≠ 0)
    (hq : q.Monic) (hqpos : 0 < q.degree)
    (hqF : q ∣ F.map (Int.castRingHom (ZMod M)))
    (hqG : q ∣ G.map (Int.castRingHom (ZMod M)))
    (B : ℤ) (hB : 0 ≤ B)
    (hFB : ∀ k, |F.coeff k| ≤ B) (hGB : ∀ k, |G.coeff k| ≤ B)
    (hlarge : ((u + v).factorial : ℤ) * B ^ (u + v) < (M : ℤ)) :
    F.resultant G u v = 0 := by
  have hdvd := modulus_dvd_resultant_of_common_divisor M hM F G q u v
    hF hG hdegree hq hqpos hqF hqG
  have hbound := int_resultant_abs_le F G u v B hB hFB hGB
  by_contra hne
  have hle : (M : ℤ) ≤ |F.resultant G u v| := by
    rw [Int.abs_eq_natAbs]
    exact_mod_cast Int.natAbs_le_of_dvd_ne_zero hdvd hne
  omega

/-- The irreducible polynomial here is a true input divisor being analysed,
not an assumed irreducible output of the recombination procedure. -/
theorem rational_irreducible_dvd_of_resultant_eq_zero
    (F G : Polynomial ℤ)
    (hF : Irreducible (F.map (Int.castRingHom ℚ)))
    (hz : F.resultant G = 0) :
    F.map (Int.castRingHom ℚ) ∣ G.map (Int.castRingHom ℚ) := by
  have hinj : Function.Injective (Int.castRingHom ℚ) := Int.cast_injective
  have hzQ : (F.map (Int.castRingHom ℚ)).resultant
      (G.map (Int.castRingHom ℚ)) = 0 := by
    rw [natDegree_map_eq_of_injective hinj, natDegree_map_eq_of_injective hinj,
      resultant_map_map, hz, map_zero]
  exact hF.dvd_iff_not_isCoprime.mpr (resultant_eq_zero_iff.mp hzQ).2

/-- Gauss descent with an arbitrary dividend, including zero and polynomials
whose content is not a unit. -/
theorem primitive_dvd_of_rational_dvd
    (F G : Polynomial ℤ) (hF : F.IsPrimitive)
    (hdiv : F.map (Int.castRingHom ℚ) ∣ G.map (Int.castRingHom ℚ)) : F ∣ G := by
  by_cases hG : G = 0
  · simp [hG]
  have hc : (G.content : ℚ) ≠ 0 := by
    exact_mod_cast mt content_eq_zero_iff.mp hG
  have hprimdiv : F.map (Int.castRingHom ℚ) ∣
      G.primPart.map (Int.castRingHom ℚ) := by
    have hd := dvd_mul_of_dvd_left hdiv (C ((G.content : ℚ)⁻¹))
    conv at hd => rhs; lhs; rw [G.eq_C_content_mul_primPart]
    rw [Polynomial.map_mul, Polynomial.map_C, mul_right_comm, ← C.map_mul] at hd
    change F.map (Int.castRingHom ℚ) ∣
      C ((G.content : ℚ) * (G.content : ℚ)⁻¹) *
        G.primPart.map (Int.castRingHom ℚ) at hd
    simpa only [mul_inv_cancel₀ hc, C_1, one_mul] using hd
  exact (hF.dvd_primPart_iff_dvd hG).mp
    ((IsPrimitive.Int.dvd_iff_map_cast_dvd_map_cast F G.primPart
      hF G.isPrimitive_primPart).mpr hprimdiv)

/-- Sufficient local precision certifies divisibility by a fixed genuine
primitive irreducible factor. The dividend need not itself be primitive. -/
theorem primitive_irreducible_dvd_of_sufficient_modulus
    (M : Nat) [NeZero M] (hM : 1 < M)
    (F G : Polynomial ℤ) (q : Polynomial (ZMod M))
    (hprimitive : F.IsPrimitive)
    (hirreducible : Irreducible (F.map (Int.castRingHom ℚ)))
    (hdegree : F.natDegree ≠ 0)
    (hq : q.Monic) (hqpos : 0 < q.degree)
    (hqF : q ∣ F.map (Int.castRingHom (ZMod M)))
    (hqG : q ∣ G.map (Int.castRingHom (ZMod M)))
    (B : ℤ) (hB : 0 ≤ B)
    (hFB : ∀ k, |F.coeff k| ≤ B) (hGB : ∀ k, |G.coeff k| ≤ B)
    (hlarge : ((F.natDegree + G.natDegree).factorial : ℤ) *
      B ^ (F.natDegree + G.natDegree) < (M : ℤ)) : F ∣ G := by
  apply primitive_dvd_of_rational_dvd F G hprimitive
  apply rational_irreducible_dvd_of_resultant_eq_zero F G hirreducible
  exact int_resultant_eq_zero_of_sufficient_modulus M hM F G q
    F.natDegree G.natDegree le_rfl le_rfl (Or.inl hdegree)
    hq hqpos hqF hqG B hB hFB hGB hlarge

/-- Two bounded integer CLD coordinates agreeing at a monic local factor
must be equal once the modulus exceeds their possible difference. -/
theorem bounded_local_cld_scalars_equal
    (M : Nat) [NeZero M] (hM : 1 < M)
    (H D : Polynomial (ZMod M)) (q : Polynomial (ZMod M))
    (a b K : ℤ) (ha : |a| ≤ K) (hb : |b| ≤ K)
    (hlarge : 2 * K < (M : ℤ))
    (hq : q.Monic) (hqpos : 0 < q.degree) (hcop : IsCoprime q D)
    (hqa : q ∣ H - C (a : ZMod M) * D)
    (hqb : q ∣ H - C (b : ZMod M) * D) : a = b := by
  haveI : Fact (1 < M) := ⟨hM⟩
  have hdiff : q ∣ C ((b - a : ℤ) : ZMod M) * D := by
    convert dvd_sub hqa hqb using 1
    simp only [Int.cast_sub, map_sub]
    ring
  have hz := constant_eq_zero_of_monic_positive_degree_dvd q
    ((b - a : ℤ) : ZMod M) hq hqpos (hcop.dvd_of_dvd_mul_right hdiff)
  have hdvd := (ZMod.intCast_zmod_eq_zero_iff_dvd (b - a) M).mp hz
  have hbound : |b - a| ≤ 2 * K := by
    calc
      |b - a| ≤ |b| + |a| := abs_sub b a
      _ ≤ 2 * K := by omega
  by_contra hne
  have hdiffne : b - a ≠ 0 := sub_ne_zero.mpr (Ne.symm hne)
  have hle : (M : ℤ) ≤ |b - a| := by
    rw [Int.abs_eq_natAbs]
    exact_mod_cast Int.natAbs_le_of_dvd_ne_zero hdvd hdiffne
  omega

end CLPoly.VanHoeij
