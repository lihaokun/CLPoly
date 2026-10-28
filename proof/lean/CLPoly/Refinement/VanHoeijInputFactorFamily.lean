import CLPoly.Refinement.VanHoeijTargetsFromGroups
import Mathlib.RingTheory.Polynomial.GaussLemma
import Mathlib.Algebra.Squarefree.Basic

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem primitive_squarefree_of_good_prime
    (source : Polynomial Int) (prime : Nat) [Fact (Nat.Prime prime)]
    (hprimitive : source.IsPrimitive)
    (hleading : (source.leadingCoeff : ZMod prime) ≠ 0)
    (hsquarefree : Squarefree (source.map (Int.castRingHom (ZMod prime)))) :
    Squarefree source := by
  intro factor hsquare
  have hd : factor ∣ source := (dvd_mul_right factor factor).trans hsquare
  have hp := Polynomial.isPrimitive_of_dvd hprimitive hd
  have hl := divisor_leading_mod_nonzero prime source factor hleading hd
  have hm : factor.map (Int.castRingHom (ZMod prime)) *
      factor.map (Int.castRingHom (ZMod prime)) ∣ source.map (Int.castRingHom (ZMod prime)) := by
    obtain ⟨quotient, hq⟩ := hsquare
    exact ⟨quotient.map (Int.castRingHom (ZMod prime)), by rw [hq, Polynomial.map_mul, Polynomial.map_mul]⟩
  have hu := hsquarefree _ hm
  have hz : factor.natDegree = 0 := by
    rw [← Polynomial.natDegree_map_of_leadingCoeff_ne_zero (Int.castRingHom (ZMod prime)) hl]
    exact Polynomial.natDegree_eq_zero_of_isUnit hu
  have he := Polynomial.eq_C_of_natDegree_eq_zero hz
  have hc : IsUnit (factor.coeff 0) :=
    (Polynomial.isPrimitive_iff_isUnit_of_C_dvd.mp hp) _ (by rw [he]; simp)
  rw [he]
  exact Polynomial.isUnit_C.mpr hc

theorem normalize_eq_self_of_positive_leading
    (source : Polynomial Int) (hpositive : (0 : Int) < source.leadingCoeff) :
    normalize source = source := by
  rw [normalize_apply, Polynomial.coe_normUnit, Int.normUnit_eq, if_pos hpositive.le]
  simp

/-- A source-derived factor family used to analyse the actual lattice and
validation runs. Its existence alone makes no claim about algorithm output. -/
theorem exists_positive_primitive_factor_family
    (source : Polynomial Int) (hprimitive : source.IsPrimitive)
    (hpositive : (0 : Int) < source.leadingCoeff) (hsquarefree : Squarefree source) :
    ∃ factors : Array (Polynomial Int),
      factors.toList.prod = source ∧ factors.toList.Nodup ∧
      (∀ i : Fin factors.size, factors[i.val].IsPrimitive) ∧
      (∀ i : Fin factors.size, (0 : Int) < factors[i.val].leadingCoeff) ∧
      (∀ i : Fin factors.size, 0 < factors[i.val].natDegree) ∧
      (∀ i : Fin factors.size, Irreducible (factors[i.val].map (Int.castRingHom QQ))) := by
  classical
  let factors := (UniqueFactorizationMonoid.normalizedFactors source).toList.toArray
  have hf : source ≠ 0 := hprimitive.ne_zero
  have hmem : ∀ i : Fin factors.size,
      factors[i.val] ∈ UniqueFactorizationMonoid.normalizedFactors source := by
    intro i
    have hi := Array.getElem_mem_toList (xs := factors) i.isLt
    simpa only [factors, List.toList_toArray, Multiset.mem_toList] using hi
  have hp : ∀ i : Fin factors.size, factors[i.val].IsPrimitive := by
    intro i
    exact Polynomial.isPrimitive_of_dvd hprimitive
      (UniqueFactorizationMonoid.dvd_of_mem_normalizedFactors (hmem i))
  have hirr : ∀ i : Fin factors.size,
      Irreducible (factors[i.val].map (Int.castRingHom QQ)) := by
    intro i
    exact (Polynomial.IsPrimitive.Int.irreducible_iff_irreducible_map_cast (hp i)).mp
      (UniqueFactorizationMonoid.irreducible_of_normalized_factor _ (hmem i))
  refine ⟨factors, ?_, ?_, hp, ?_, ?_, hirr⟩
  · simp only [factors, List.toList_toArray, Multiset.prod_toList]
    exact (UniqueFactorizationMonoid.prod_normalizedFactors_eq hf).trans
      (normalize_eq_self_of_positive_leading source hpositive)
  · apply Multiset.coe_nodup.mp
    simpa only [factors, List.toList_toArray, Multiset.coe_toList] using
      (UniqueFactorizationMonoid.squarefree_iff_nodup_normalizedFactors hf).mp hsquarefree
  · intro i
    have hn := UniqueFactorizationMonoid.normalize_normalized_factor _ (hmem i)
    have hl : normalize factors[i.val].leadingCoeff = factors[i.val].leadingCoeff := by
      rw [← Polynomial.leadingCoeff_normalize, hn]
    exact lt_of_le_of_ne (Int.nonneg_of_normalize_eq_self hl)
      (Ne.symm (Polynomial.leadingCoeff_ne_zero.mpr (hp i).ne_zero))
  · intro i
    have h := Irreducible.natDegree_pos (hirr i)
    rwa [Polynomial.natDegree_map_eq_of_injective (Int.castRingHom QQ).injective_int] at h

end Refinement.StrictRecombine
