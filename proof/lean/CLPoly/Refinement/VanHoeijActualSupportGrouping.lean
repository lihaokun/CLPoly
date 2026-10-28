import CLPoly.Refinement.VanHoeijGroupingTotal
import CLPoly.Refinement.VanHoeijRetainedNonempty

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Sufficient precision puts members of each genuine divisor's actual
Hensel support in a common actual output group. The grouping run is produced
by the theorem, not assumed. -/
theorem hensel_grouping_complete_on_divisor
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
    (scale divisor bound : ZZ) (matrix : Generated.StrictRecombine.LLLMatrix)
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
    (height : ZZ)
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
    ∃ (indices : Array Nat) (groups : Array (Array Int32)),
      LegalCombination lifted.1.size indices.size indices ∧
      factor.map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          (indices.toList.map fun i => StrictHensel.toPolyMod modulus lifted.1[i]!).prod ∧
      Generated.StrictRecombine.extractCandidates
        (Generated.StrictVanHoeij.shortSpanFromNorms output.2.2.2 bound)
        output.2.1 lifted.1.size = .ok groups ∧
      ∀ left right : Fin lifted.1.size, left.val ∈ indices.toList → right.val ∈ indices.toList →
        ∃ classId, CandidateMember groups classId left.val.toUInt32.toInt32 ∧
          CandidateMember groups classId right.val.toUInt32.toInt32 := by
  obtain ⟨indices, hlegal, hsupport, hconstant⟩ :=
    hensel_retained_rows_constant_on_divisor hcount hfactors hleading hselection hentry
      modulus hmodulus original truncated data scale divisor bound matrix htruncate hcanonical
      hscale hm hinput hB hsize hbound output hrun factor complement hfactor
      hprimitive hirreducible hfactorDegree hresultant hcoordinates
  have htruncSize : truncated.size = original.size :=
    (truncateCld_spec original truncated divisor htruncate).2.1
  have hfactorCount : lifted.1.size ≤ matrix.size := by
    have hs := hcanonical.size_eq
    rw [htruncSize, data.size_eq] at hs
    omega
  have hnonempty := hensel_lll_retained_nonempty hcount hfactors hleading hselection hentry
    modulus hmodulus original truncated data scale divisor height bound _ _ matrix
    htruncate hcanonical hheight hh hcutoff hinput output hrun
  obtain ⟨groups, hg, hgroups⟩ := lllReduceWithNorms_grouping_complete matrix hinput bound
    output hrun lifted.1.size hfactorCount hnonempty
  refine ⟨indices, groups, hlegal, hsupport, hg, ?_⟩
  intro left right hl hr
  apply hgroups left right
  intro row hkeep
  exact hconstant row hkeep
    ⟨left.val, lt_of_lt_of_le left.isLt hfactorCount⟩
    ⟨right.val, lt_of_lt_of_le right.isLt hfactorCount⟩
    (by rw [data.size_eq]; exact left.isLt)
    (by rw [data.size_eq]; exact right.isLt) hl hr

end Refinement.StrictRecombine
