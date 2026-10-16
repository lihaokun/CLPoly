import CLPoly.Refinement.VanHoeijZeroWidthSupport
import CLPoly.Refinement.VanHoeijSingletonTargets
import CLPoly.Refinement.VanHoeijInputFactorFamily
import CLPoly.Refinement.VanHoeijHenselAtoms

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- At sufficient precision, zero truncated width forces the source entry's
actual singleton partition to validate successfully. No lattice execution or
assumed output irreducibility is used for this branch. -/
theorem hensel_zero_width_singleton_validation_complete
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {aTarget : Int32} {lifted : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime aTarget lifted)
    (modulus : Nat) (hmodulus : lifted.2 = (modulus : Int))
    (original truncated : Array SparsePolyZZ)
    (data : CldSeparationData f lifted.1 original modulus)
    (divisor : ZZ)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hwidth : Generated.StrictRecombine.cldSpiralWidth truncated = 0)
    (count : Int32) (groups : Array (Array Int32))
    (hcountValue : count.toInt.toNat = lifted.1.size)
    (hgroups : Generated.StrictVanHoeij.singletonGroups count = .ok groups)
    (hsourceCanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hsourcePrimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hsourcePositive : (0 : Int) < (SparsePolyZZ.toPoly f).leadingCoeff)
    (hsourceFits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63)
    (hresultant : separationResultantThreshold (SparsePolyZZ.toPoly f) lifted.1.size divisor < (modulus : Int))
    (hcoordinates : 2 * separationCoordinateBound lifted.1.size (SparsePolyZZ.toPoly f).natDegree < (modulus : Int))
    (before : Generated.StrictVanHoeij.VanHoeijResult) :
    ∃ after, Generated.StrictVanHoeij.validatePartition f lifted.1 modulus groups before = .ok (true, after) := by
  have hsq := primitive_squarefree_of_good_prime (SparsePolyZZ.toPoly f) selection.prime.toNat
    hsourcePrimitive hselection.goodPrime.lc_nonzero hselection.goodPrime.sqfree
  obtain ⟨factors, hproduct, hnodup, hprimitive, hpositive, hdegree, hirreducible⟩ :=
    exists_positive_primitive_factor_family (SparsePolyZZ.toPoly f) hsourcePrimitive hsourcePositive hsq
  have hexact : ∀ i : Fin factors.size, ∃ index : Fin lifted.1.size,
      factors[i.val].map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factors[i.val].leadingCoeff : ZMod modulus) *
          StrictHensel.toPolyMod modulus lifted.1[index.val]! := by
    intro i
    have hd : factors[i.val] ∣ SparsePolyZZ.toPoly f := by
      rw [← hproduct]; exact List.dvd_prod (Array.getElem_mem_toList i.isLt)
    obtain ⟨complement, hc⟩ := hd
    exact hensel_zero_width_divisor_singleton_support hcount hfactors hleading hselection hentry
      modulus hmodulus original truncated data divisor htruncate hwidth factors[i.val] complement hc
      (hprimitive i) (hirreducible i) (hdegree i) hresultant hcoordinates
  obtain ⟨hp, hseparate⟩ := hensel_physical_atoms_prime_separated hcount hfactors
    hleading hselection hentry
  have hatom : ∀ i : Fin lifted.1.size, StrictHensel.toPolyMod selection.prime.toNat lifted.1[i.val] ∣
      (SparsePolyZZ.toPoly f).map (Int.castRingHom (ZMod selection.prime.toNat)) := by
    intro i
    have ha := hselection.productAssociated.trans
      (StrictFactorZZ.selectionHenselFactors_product_associated hcount hfactors hleading hselection hentry)
    have hmem : StrictHensel.toPolyMod selection.prime.toNat lifted.1[i.val] ∈
        lifted.1.toList.map (StrictHensel.toPolyMod selection.prime.toNat) :=
      List.mem_map.mpr ⟨lifted.1[i.val], Array.getElem_mem_toList i.isLt, rfl⟩
    exact (List.dvd_prod hmem).trans ha.symm.dvd
  obtain ⟨exponent, he, hpower⟩ := hentry.outputModulus_eq_prime_pow
  have hpower' : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmodulus.symm.trans hpower
  have hd : selection.prime.toNat ∣ modulus := by
    rw [hpower']; exact dvd_pow_self _ (Nat.ne_of_gt he)
  have hmpos : 0 < modulus := lt_trans Nat.zero_lt_one data.modulus_gt_one
  have hdivisor : 0 ≤ divisor := (truncateCld_spec original truncated divisor htruncate).1.le
  obtain ⟨targets, ht⟩ := exactPartitionTargets_of_singleton_supports (SparsePolyZZ.toPoly f)
    factors hproduct hnodup hsourcePrimitive.ne_zero data.source_degree_pos hprimitive hpositive hdegree
    selection.prime.toNat modulus hmpos hd lifted.1 count groups hcountValue hgroups hp
    hseparate hatom hselection.goodPrime.lc_nonzero divisor hdivisor hresultant hexact
  have hnonempty : ∀ i, i < lifted.1.size → 0 < lifted.1[i]!.size := by
    intro i hi
    obtain ⟨head, tail, heq⟩ := hentry.outputOneHead i hi
    rw [getElem!_pos lifted.1 i hi]
    simpa only [Array.length_toList, List.length_cons] using
      (show 0 < lifted.1[i].toList.length by rw [heq]; simp)
  have hvalid := singletonGroups_indices_valid count groups hgroups
  rw [hcountValue] at hvalid
  exact validatePartition_complete_of_targets f selection.prime.toNat modulus
    hmpos hd lifted.1 groups targets ht hsourceCanonical data.source_degree_pos hsourceFits
    hselection.goodPrime.lc_nonzero hvalid
    (fun i hi => by simpa only [getElem!_pos lifted.1 i hi] using hentry.outputCanonical i hi)
    hnonempty
    (fun i hi => by simpa only [getElem!_pos lifted.1 i hi] using hentry.outputToPolyMonic i hi)
    divisor hdivisor hresultant before

end Refinement.StrictRecombine
