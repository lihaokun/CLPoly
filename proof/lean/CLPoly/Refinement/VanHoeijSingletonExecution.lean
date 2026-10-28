import CLPoly.Refinement.VanHoeijSingleton
import CLPoly.Refinement.VanHoeijSingletonConditions

set_option autoImplicit false

namespace Refinement.StrictRecombine
open CLPoly.Math

theorem hensel_constructed_singletons_validation_irreducible
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : ∀ factor ∈ selection.factors.toList,
      SparsePolyZp.Canonical selection.prime.toNat factor)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) =
        (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f
      selection.factors selection.prime aTarget output)
    (modulus : Nat) (hmod : output.2 = (modulus : Int))
    (count : Int32) (hcountSize : count.toInt.toNat = output.1.size)
    (groups : Array (Array Int32))
    (hconstruct : Generated.StrictVanHoeij.singletonGroups count = .ok groups)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hdegree : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hvalidate : Generated.StrictVanHoeij.validatePartition f output.1 modulus groups before =
      .ok (true, after)) :
    factorArrayProduct after.factors = SparsePolyZZ.toPoly f ∧
      ∀ factor ∈ after.factors.toList, Irreducible (SparsePolyZZ.toPoly factor) := by
  refine ⟨validatePartition_product f output.1 modulus groups before after
    hcanonical hdegree hvalidate, ?_⟩
  obtain ⟨exponent, he, hpower⟩ := hentry.outputModulus_eq_prime_pow
  have hmodulus : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmod.symm.trans hpower
  have hm : 0 < modulus := by
    rw [hmodulus]
    exact pow_pos (Fact.out : Nat.Prime selection.prime.toNat).pos _
  have hd : selection.prime.toNat ∣ modulus := by
    rw [hmodulus]
    exact dvd_pow_self _ (Nat.ne_of_gt he)
  have hv : ∀ i (hi : i < groups.size),
      CandidateIndicesValid groups[i] (Array.replicate output.1.size false) := by
    simpa only [hcountSize] using singletonGroups_indices_valid count groups hconstruct
  exact validatePartition_singletons_irreducible f output.1 modulus selection.prime.toNat
    groups before after hm hd hv (singletonGroups_singleton count groups hconstruct)
    (StrictFactorZZ.selectionHenselFactors_mod_irreducible hcount hfactors hleadingSemantic
      hselection hentry) hcanonical hdegree hprimitive hselection.goodPrime.lc_nonzero hvalidate

end Refinement.StrictRecombine
