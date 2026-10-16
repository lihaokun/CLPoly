import CLPoly.Refinement.VanHoeijRowSeparation
import CLPoly.Refinement.VanHoeijHensel

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- At the explicit precision threshold, the support returned by the
actual Hensel support theorem is constant on every actual retained LLL row.
The support equation is derived here, not supplied as an oracle. -/
theorem hensel_retained_rows_constant_on_divisor
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
    (factor complement : Polynomial Int) (hfactor : SparsePolyZZ.toPoly f = factor * complement)
    (hprimitive : factor.IsPrimitive) (hirreducible : Irreducible (factor.map (Int.castRingHom QQ)))
    (hfactorDegree : factor.natDegree ≠ 0)
    (hresultant : separationResultantThreshold (SparsePolyZZ.toPoly f) lifted.1.size divisor < (modulus : Int))
    (hcoordinates : 2 * separationCoordinateBound lifted.1.size (SparsePolyZZ.toPoly f).natDegree < (modulus : Int)) :
    ∃ indices : Array Nat,
      LegalCombination lifted.1.size indices.size indices ∧
      factor.map (Int.castRingHom (ZMod modulus)) =
        Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          (indices.toList.map fun i => StrictHensel.toPolyMod modulus lifted.1[i]!).prod ∧
      ∀ row : Fin matrix.size,
        row.val < Generated.StrictVanHoeij.retainedPrefixCount output.2.2.2 (bound : QQ) output.2.2.2.size →
        ∀ left right : Fin matrix.size, left.val < original.size → right.val < original.size →
          left.val ∈ indices.toList → right.val ∈ indices.toList →
          (output.2.1[row.val]!)[left.val]! = (output.2.1[row.val]!)[right.val]! := by
  obtain ⟨indices, exponent, hlegal, he, hp, hs⟩ :=
    StrictFactorZZ.hensel_divisor_physical_support hcount hfactors hleading hselection hentry
      factor complement hfactor
  have hpower : modulus = selection.prime.toNat ^ exponent := by
    exact_mod_cast hmodulus.symm.trans hp
  have hsupport : factor.map (Int.castRingHom (ZMod modulus)) =
      Polynomial.C (factor.leadingCoeff : ZMod modulus) *
        (indices.toList.map fun i => StrictHensel.toPolyMod modulus lifted.1[i]!).prod := by
    rw [hpower]
    exact hs
  refine ⟨indices, hlegal, hsupport, ?_⟩
  intro row hkeep left right hl hr hls hrs
  exact data.retained_row_equal_on_support truncated scale divisor bound matrix htruncate hcanonical
    hscale hm hinput hB hsize hbound output hrun row hkeep factor ⟨complement, hfactor⟩
    hprimitive hirreducible hfactorDegree indices.toList _ hsupport left right hl hr hls hrs
    hresultant hcoordinates

end Refinement.StrictRecombine
