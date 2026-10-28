import CLPoly.Refinement.VanHoeijInputFactorFamily
import CLPoly.Refinement.VanHoeijActualExactGrouping
import CLPoly.Refinement.VanHoeijHenselAtoms

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- At sufficient precision, the same actual full-width lattice execution
produces groups whose actual partition validation succeeds. The mathematical
factor family and ordered targets are derived, not assumed. Reachability of
this precision by the production controller is a separate obligation. -/
theorem hensel_lll_validation_complete
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
    (hindexFits : lifted.1.size ≤ 2 ^ 31)
    (modulus : Nat) (hmodulus : lifted.2 = (modulus : Int))
    (original truncated : Array SparsePolyZZ)
    (data : CldSeparationData f lifted.1 original modulus)
    (scale divisor height bound : ZZ) (matrix : Generated.StrictRecombine.LLLMatrix)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hcanonical : CanonicalCldBasis truncated (Generated.StrictRecombine.cldSpiralWidth truncated)
      scale ((modulus : Int) / divisor) (Generated.StrictRecombine.cldSpiralWidth truncated) matrix)
    (hscale : 1 ≤ scale) (hm : 1 ≤ (modulus : Int) / divisor)
    (hinput : ConcreteLLLInputValid matrix) (hB : 0 ≤ bound)
    (hsize : matrix.size ≤ lifted.1.size + (SparsePolyZZ.toPoly f).natDegree)
    (hbound : bound ≤ separationCutoff lifted.1.size (SparsePolyZZ.toPoly f).natDegree)
    (output : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hrun : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution matrix hinput bound = .ok output)
    (hheight : height = ((SparsePolyZZ.toPoly f).natDegree *
      2 ^ ((SparsePolyZZ.toPoly f).natDegree - 1) *
      (Finset.range ((SparsePolyZZ.toPoly f).natDegree + 1)).sum
        (fun k => ((SparsePolyZZ.toPoly f).coeff k).natAbs) : Nat))
    (hh : height < (modulus : Int))
    (hcutoff : (lifted.1.size : ZZ) * scale * scale +
      (Generated.StrictRecombine.cldSpiralWidth truncated : ZZ) *
      ((height + divisor - 1) / divisor + 2 * (lifted.1.size : Int)) *
      ((height + divisor - 1) / divisor + 2 * (lifted.1.size : Int)) ≤ bound)
    (hsourceCanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hsourcePrimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hsourcePositive : (0 : Int) < (SparsePolyZZ.toPoly f).leadingCoeff)
    (hsourceFits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63)
    (hresultant : separationResultantThreshold (SparsePolyZZ.toPoly f) lifted.1.size divisor < (modulus : Int))
    (hcoordinates : 2 * separationCoordinateBound lifted.1.size (SparsePolyZZ.toPoly f).natDegree < (modulus : Int))
    (before : Generated.StrictVanHoeij.VanHoeijResult) :
    ∃ groups after,
      Generated.StrictRecombine.extractCandidates
        (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound)
        output.2.1 lifted.1.size = .ok groups ∧
      Generated.StrictVanHoeij.validatePartition f lifted.1 modulus groups before = .ok (true, after) := by
  have htruncSize := (truncateCld_spec original truncated divisor htruncate).2.1
  have hfactorCount : lifted.1.size ≤ matrix.size := by
    rw [hcanonical.size_eq, htruncSize, data.size_eq]
    omega
  let rows := Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound
  have hvalid := lllReduceWithNorms_candidate_columns_valid matrix hinput bound lifted.1.size
    hfactorCount output hrun
  obtain ⟨groups, hg⟩ := extractCandidates_succeeds output.2.1 rows lifted.1.size hvalid
  have hsq := primitive_squarefree_of_good_prime (SparsePolyZZ.toPoly f) selection.prime.toNat
    hsourcePrimitive hselection.goodPrime.lc_nonzero hselection.goodPrime.sqfree
  obtain ⟨factors, hproduct, hnodup, hprimitive, hpositive, hdegree, hirreducible⟩ :=
    exists_positive_primitive_factor_family (SparsePolyZZ.toPoly f) hsourcePrimitive hsourcePositive hsq
  have hexact : ∀ i : Fin factors.size, ∃ (indices : Array Nat) (classId : Nat),
      LegalCombination lifted.1.size indices.size indices ∧ classId < groups.size ∧
      factors[i.val].map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factors[i.val].leadingCoeff : ZMod modulus) *
          (indices.toList.map fun j => StrictHensel.toPolyMod modulus lifted.1[j]!).prod ∧
      factors[i.val].map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factors[i.val].leadingCoeff : ZMod modulus) *
          SelectedProductMod modulus groups[classId]! lifted.1 0 ∧
      ∀ index, index < lifted.1.size →
        (CandidateMember groups classId index.toUInt32.toInt32 ↔ index ∈ indices.toList) := by
    intro i
    have hd : factors[i.val] ∣ SparsePolyZZ.toPoly f := by
      rw [← hproduct]; exact List.dvd_prod (Array.getElem_mem_toList i.isLt)
    obtain ⟨complement, hc⟩ := hd
    obtain ⟨indices, returned, classId, hl, hs, he, hclass, hgroup, hmembers⟩ :=
      hensel_exact_support_group hcount hfactors hleading hselection hentry hindexFits modulus
        hmodulus original truncated data scale divisor height bound matrix htruncate hcanonical
        hscale hm hinput hB hsize hbound output hrun hheight hh hcutoff factors[i.val]
        complement hc (hprimitive i) (hirreducible i) (Nat.ne_of_gt (hdegree i)) hresultant hcoordinates
    have heq : returned = groups := Except.ok.inj (he.symm.trans hg)
    subst returned
    exact ⟨indices, classId, hl, hclass, hs, hgroup, hmembers⟩
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
  obtain ⟨targets, ht⟩ := exactPartitionTargets_of_individual_groups (SparsePolyZZ.toPoly f)
    factors hproduct hnodup hsourcePrimitive.ne_zero data.source_degree_pos hprimitive hpositive hdegree
    selection.prime.toNat modulus hmpos hd lifted.1 output.2.1 rows groups hindexFits hg hp
    hseparate hatom hselection.goodPrime.lc_nonzero divisor hdivisor hresultant hexact
  have hrows : rows.isEmpty = false := by
    cases hempty : rows.isEmpty with
    | false => rfl
    | true =>
      have hgroups : groups = #[] := by
        simpa [Generated.StrictRecombine.extractCandidates, hempty] using hg.symm
      have htargets : targets = #[] := Array.size_eq_zero_iff.mp (by simpa [hgroups] using ht.size_eq)
      have hsource := ht.product_eq
      simp only [htargets, Array.toList_empty, List.prod_nil] at hsource
      have hn := data.source_degree_pos
      rw [← hsource] at hn
      simp at hn
  have hnonempty : ∀ i, i < lifted.1.size → 0 < lifted.1[i]!.size := by
    intro i hi
    obtain ⟨head, tail, heq⟩ := hentry.outputOneHead i hi
    have hlen := congrArg List.length heq
    rw [getElem!_pos lifted.1 i hi]
    simpa only [Array.length_toList, List.length_cons] using
      (show 0 < lifted.1[i].toList.length by rw [heq]; simp)
  obtain ⟨after, ha⟩ := validatePartition_complete_of_targets f selection.prime.toNat modulus
    hmpos hd lifted.1 groups targets ht hsourceCanonical data.source_degree_pos hsourceFits
    hselection.goodPrime.lc_nonzero
    (extractCandidates_indices_valid rows output.2.1 lifted.1.size groups hindexFits hrows hg)
    (fun i hi => by simpa only [getElem!_pos lifted.1 i hi] using hentry.outputCanonical i hi)
    hnonempty
    (fun i hi => by simpa only [getElem!_pos lifted.1 i hi] using hentry.outputToPolyMonic i hi)
    divisor hdivisor hresultant before
  exact ⟨groups, after, hg, ha⟩

end Refinement.StrictRecombine
