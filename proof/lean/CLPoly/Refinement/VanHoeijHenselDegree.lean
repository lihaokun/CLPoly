import CLPoly.Refinement.FactorZZ

set_option autoImplicit false

namespace Refinement.StrictRecombine

theorem hensel_output_degree_le_source
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {target : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f selection.factors
      selection.prime target output)
    (index : Nat) (hi : index < output.1.size) :
    (SparsePolyZZ.toPoly output.1[index]).natDegree ≤ (SparsePolyZZ.toPoly f).natDegree := by
  have hproduct := hselection.productAssociated.trans
    (StrictFactorZZ.selectionHenselFactors_product_associated hcount hfactors
      hleadingSemantic hselection hentry)
  have hmem : StrictHensel.toPolyMod selection.prime.toNat output.1[index] ∈
      output.1.toList.map (StrictHensel.toPolyMod selection.prime.toNat) := by
    exact List.mem_map.mpr ⟨output.1[index], Array.getElem_mem_toList hi, rfl⟩
  have hdvd := (List.dvd_prod hmem).trans hproduct.symm.dvd
  have hsourceNe : Polynomial.map (Int.castRingHom (ZMod selection.prime.toNat))
      (SparsePolyZZ.toPoly f) ≠ 0 := by
    intro hzero
    have hlc := Polynomial.leadingCoeff_map_of_leadingCoeff_ne_zero
      (Int.castRingHom (ZMod selection.prime.toNat)) hselection.goodPrime.lc_nonzero
    rw [hzero] at hlc
    exact hselection.goodPrime.lc_nonzero hlc.symm
  have hdegree := Polynomial.natDegree_le_of_dvd hdvd hsourceNe
  have hmonic := hentry.outputToPolyMonic index hi
  have hmap := hmonic.natDegree_map (Int.castRingHom (ZMod selection.prime.toNat))
  change (StrictHensel.toPolyMod selection.prime.toNat output.1[index]).natDegree =
    (SparsePolyZZ.toPoly output.1[index]).natDegree at hmap
  rw [hmap] at hdegree
  exact hdegree.trans Polynomial.natDegree_map_le

theorem hensel_output_head_degree_fits
    {termination : Generated.StrictHensel.DivmodTermination}
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    {target : Int32} {output : Array SparsePolyZZ × ZZ}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hcount : 2 ≤ selection.factors.size)
    (hfactors : StrictSelectPrime.SelectionPhysical selection)
    (hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hentry : StrictHensel.HenselLiftEntryCorrect termination f selection.factors
      selection.prime target output)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (index : Nat) (hi : index < output.1.size) :
    output.1[index][0]!.1.deg < 2 ^ 63 := by
  have hdegree := hensel_output_degree_le_source hcount hfactors hleadingSemantic
    hselection hentry index hi
  have hnonempty : 0 < output.1[index].size := by
    obtain ⟨head, tail, hlist⟩ := hentry.outputOneHead index hi
    have hlength := congrArg List.length hlist
    simp only [Array.length_toList, List.length_cons] at hlength
    omega
  rw [sparsePolyZZ_natDegree_eq_head _ (hentry.outputCanonical index hi) hnonempty] at hdegree
  rw [getElem!_pos _ 0 hnonempty]
  omega

end Refinement.StrictRecombine
