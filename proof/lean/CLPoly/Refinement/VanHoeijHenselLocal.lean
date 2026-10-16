import CLPoly.Refinement.VanHoeijHenselDegree
import CLPoly.Refinement.VanHoeijLocalResidual
import Mathlib.FieldTheory.Perfect

set_option autoImplicit false

namespace Refinement.StrictRecombine

open CLPoly.Math

theorem hensel_output_modulus_factor_dvd
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime aTarget output)
    (modulus : Nat) (hmodulus : output.2 = (modulus : Int))
    (i : Nat) (hi : i < output.1.size) :
    StrictHensel.toPolyMod modulus output.1[i]! ∣ StrictHensel.toPolyMod modulus f := by
  obtain ⟨m, scale, hm, hu, hp⟩ :=
    StrictFactorZZ.selectionHenselFactors_normalized_product_eq_unit_mul_source
      hcount hfactors hleading hselection hentry
  have heq : m = modulus := by exact_mod_cast hm.symm.trans hmodulus
  subst m
  have hmem : StrictHensel.toPolyMod modulus output.1[i]! ∈
      output.1.toList.map (StrictHensel.toPolyMod modulus) := by
    apply List.mem_map.mpr
    exact ⟨output.1[i], Array.getElem_mem_toList hi, by rw [getElem!_pos output.1 i hi]⟩
  have hd := List.dvd_prod hmem
  rw [hp] at hd
  exact (Polynomial.isUnit_C.mpr hu).dvd_mul_left.mp hd

theorem hensel_output_modulus_pairwise_coprime
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime aTarget output)
    (modulus : Nat) (hmodulus : output.2 = (modulus : Int))
    (i j : Nat) (hi : i < output.1.size) (hj : j < output.1.size) (hne : i ≠ j) :
    IsCoprime (StrictHensel.toPolyMod modulus output.1[i]!)
      (StrictHensel.toPolyMod modulus output.1[j]!) := by
  have hc := StrictFactorZZ.henselFactors_mod_pairwise_coprime hcount
    (fun adjusted ha => StrictFactorZZ.selectionAdjustedFactors_pairwise_isCoprime
      hfactors hleading hselection ha) hentry
  have hp : IsCoprime (StrictHensel.toPolyMod selection.prime.toNat output.1[i]!)
      (StrictHensel.toPolyMod selection.prime.toNat output.1[j]!) := by
    rw [getElem!_pos output.1 i hi, getElem!_pos output.1 j hj]
    rcases lt_or_gt_of_ne hne with hij | hji
    · exact hc i j hi hj hij
    · exact (hc j i hj hi hji).symm
  obtain ⟨exponent, he, hm⟩ := hentry.outputModulus_eq_prime_pow
  have hpower : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmodulus.symm.trans hm
  subst modulus
  exact isCoprime_lift_pk selection.prime.toNat hselection.goodPrime.prime exponent
    (Nat.ne_of_gt he) (SparsePolyZZ.toPoly output.1[i]!) (SparsePolyZZ.toPoly output.1[j]!) hp

theorem hensel_output_modulus_derivative_coprime
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime aTarget output)
    (modulus : Nat) (hmodulus : output.2 = (modulus : Int))
    (i : Nat) (hi : i < output.1.size) :
    IsCoprime (StrictHensel.toPolyMod modulus output.1[i]!)
      (Polynomial.derivative (StrictHensel.toPolyMod modulus f)) := by
  have hproduct := hselection.productAssociated.trans
    (StrictFactorZZ.selectionHenselFactors_product_associated hcount hfactors hleading hselection hentry)
  have hmem : StrictHensel.toPolyMod selection.prime.toNat output.1[i]! ∈
      output.1.toList.map (StrictHensel.toPolyMod selection.prime.toNat) := by
    apply List.mem_map.mpr
    exact ⟨output.1[i], Array.getElem_mem_toList hi, by rw [getElem!_pos output.1 i hi]⟩
  have hd := (List.dvd_prod hmem).trans hproduct.symm.dvd
  have hsep := PerfectField.separable_iff_squarefree.mpr hselection.goodPrime.sqfree
  have hcop := hsep.of_isCoprime_of_dvd_left hd
  obtain ⟨exponent, he, hm⟩ := hentry.outputModulus_eq_prime_pow
  have hpower : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmodulus.symm.trans hm
  subst modulus
  have hl := isCoprime_lift_pk selection.prime.toNat hselection.goodPrime.prime exponent
    (Nat.ne_of_gt he) (SparsePolyZZ.toPoly output.1[i]!)
    (Polynomial.derivative (SparsePolyZZ.toPoly f)) (by
      simpa only [Polynomial.derivative_map] using hcop)
  simpa only [StrictHensel.toPolyMod, Polynomial.derivative_map] using hl

theorem hensel_output_modulus_degree_pos
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime aTarget output)
    (modulus : Nat) (hm : 1 < modulus) (i : Nat) (hi : i < output.1.size) :
    0 < (StrictHensel.toPolyMod modulus output.1[i]!).degree := by
  haveI : Fact (1 < modulus) := ⟨hm⟩
  have hirr := StrictFactorZZ.selectionHenselFactors_mod_irreducible
    hcount hfactors hleading hselection hentry i hi
  have hp := Polynomial.natDegree_pos_iff_degree_pos.mpr (Polynomial.degree_pos_of_irreducible hirr)
  have hmonic := hentry.outputToPolyMonic i hi
  rw [getElem!_pos output.1 i hi, ← Polynomial.natDegree_pos_iff_degree_pos]
  change 0 < ((SparsePolyZZ.toPoly output.1[i]).map (Int.castRingHom (ZMod modulus))).natDegree
  rw [hmonic.natDegree_map]
  change 0 < ((SparsePolyZZ.toPoly output.1[i]).map
    (Int.castRingHom (ZMod selection.prime.toNat))).natDegree at hp
  rwa [hmonic.natDegree_map] at hp

end Refinement.StrictRecombine
