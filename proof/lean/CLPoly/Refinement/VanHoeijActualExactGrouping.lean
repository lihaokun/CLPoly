import CLPoly.Refinement.VanHoeijGroupSupportProduct

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- The same actual Hensel support supplies both the preserved indicator and
the sufficient-precision separation argument. Their combination identifies
an exact returned group, without assuming grouping success or nonempty rows. -/
theorem hensel_exact_support_group
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
    (factor complement : Polynomial Int) (hfactor : SparsePolyZZ.toPoly f = factor * complement)
    (hprimitive : factor.IsPrimitive) (hirreducible : Irreducible (factor.map (Int.castRingHom QQ)))
    (hfactorDegree : factor.natDegree ≠ 0)
    (hresultant : separationResultantThreshold (SparsePolyZZ.toPoly f) lifted.1.size divisor < (modulus : Int))
    (hcoordinates : 2 * separationCoordinateBound lifted.1.size (SparsePolyZZ.toPoly f).natDegree < (modulus : Int)) :
    ∃ (indices : Array Nat) (groups : Array (Array Int32)) (classId : Nat),
      LegalCombination lifted.1.size indices.size indices ∧
      factor.map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          (indices.toList.map fun i => StrictHensel.toPolyMod modulus lifted.1[i]!).prod ∧
      Generated.StrictRecombine.extractCandidates
        (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound)
        output.2.1 lifted.1.size = .ok groups ∧
      classId < groups.size ∧
      factor.map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          SelectedProductMod modulus groups[classId]! lifted.1 0 ∧
      ∀ index, index < lifted.1.size →
        (CandidateMember groups classId index.toUInt32.toInt32 ↔ index ∈ indices.toList) := by
  obtain ⟨indices, coefficients, hlegal, hsupport, hindicator, hshort⟩ :=
    hensel_factor_truncated_short hcount hfactors hleading hselection hentry original truncated
      modulus hmodulus divisor height scale _ _ matrix data.size_eq
      (fun i hi => data.values i (by rw [data.size_eq]; exact hi)) htruncate hcanonical
      factor complement hfactor (by have h := data.modulus_gt_one; omega)
      (truncateCld_spec original truncated divisor htruncate).1 hh
      (Nat.pos_of_ne_zero hfactorDegree) hheight
  obtain ⟨exponent, he, hp⟩ := hentry.outputModulus_eq_prime_pow
  have hp' : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmodulus.symm.trans hp
  have hd : selection.prime.toNat ∣ modulus := by
    rw [hp']; exact dvd_pow_self _ (Nat.ne_of_gt he)
  have hnonempty := divisor_modular_support_nonempty selection.prime.toNat modulus hd
    lifted.1 indices hlegal factor (Nat.pos_of_ne_zero hfactorDegree)
    (divisor_leading_mod_nonzero selection.prime.toNat (SparsePolyZZ.toPoly f) factor
      hselection.goodPrime.lc_nonzero ⟨complement, hfactor⟩) hsupport
  have htruncSize := (truncateCld_spec original truncated divisor htruncate).2.1
  have hfactorCount : lifted.1.size ≤ matrix.size := by
    rw [hcanonical.size_eq, htruncSize, data.size_eq]
    omega
  have hshort' : dotProduct
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ))
      (fun column => ((Matrix.vecMul coefficients
        (basisPrefixMatrix matrix matrix.size matrix.size) column : ZZ) : QQ)) ≤ (bound : QQ) :=
    hshort.trans (by exact_mod_cast hcutoff)
  obtain ⟨groups, classId, hg, hc, hmembers⟩ :=
    lllReduceWithNorms_exact_support_group matrix hinput bound output hrun lifted.1.size
      hfactorCount hindexFits indices hlegal hnonempty coefficients hindicator hshort' (by
        intro row hkeep left right hl hr
        exact data.retained_row_equal_on_support truncated scale divisor bound matrix htruncate
          hcanonical hscale hm hinput hB hsize hbound output hrun row hkeep factor
          ⟨complement, hfactor⟩ hprimitive hirreducible hfactorDegree indices.toList _ hsupport
          ⟨left.val, lt_of_lt_of_le left.isLt hfactorCount⟩
          ⟨right.val, lt_of_lt_of_le right.isLt hfactorCount⟩
          (by rw [data.size_eq]; exact left.isLt)
          (by rw [data.size_eq]; exact right.isLt) hl hr hresultant hcoordinates)
  have hproduct := extractCandidates_exact_support_product output.2.1
    (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound) lifted.1 groups
    hindexFits hg classId hc indices hlegal hmembers modulus
  refine ⟨indices, groups, classId, hlegal, hsupport, hg, hc, ?_, hmembers⟩
  rw [hproduct]
  exact hsupport

end Refinement.StrictRecombine
