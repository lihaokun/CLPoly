import CLPoly.Refinement.VanHoeijGrouping

set_option autoImplicit false

namespace Refinement.StrictRecombine
open CLPoly.Math

/-- Actual same-class members cannot disagree on the support of any genuine
positive-degree integer divisor. One support witness works for every class. -/
theorem hensel_execution_group_support
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
    (original truncated : Array SparsePolyZZ)
    (modulus : Nat) (hmod : output.2 = (modulus : Int))
    (divisor height scale : ZZ) (target columns : Nat)
    (initial matrix : Generated.StrictRecombine.LLLMatrix)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hliftedDegree : ∀ index (hi : index < output.1.size),
      output.1[index][0]!.1.deg < 2 ^ 63)
    (hcldRun : Generated.StrictRecombine.cldPolys ⟨termination⟩ f output.1 modulus = .ok original)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hscale : scale ≠ 0)
    (hinit : Generated.StrictRecombine.makeInitialMatrix truncated.size scale = .ok initial)
    (hbuild : Generated.StrictVanHoeij.buildCldMatrix initial truncated 0 target
      ((modulus : Int) / divisor) = .ok (matrix, columns))
    (factor complement : Polynomial Int)
    (hfactor : SparsePolyZZ.toPoly f = factor * complement)
    (hm : 1 < modulus) (hd : 0 < divisor) (hh : height < (modulus : Int))
    (hdegree : 0 < factor.natDegree)
    (hheightRun : Generated.StrictVanHoeij.cldBound f = .ok height)
    (bound : ZZ)
    (hbound : bound = (output.1.size : ZZ) * scale * scale + (columns : ZZ) *
      ((height + divisor - 1) / divisor + 2 * (output.1.size : Int)) *
      ((height + divisor - 1) / divisor + 2 * (output.1.size : Int)))
    (hinput : ConcreteLLLInputValid matrix)
    (reduced : Generated.StrictRecombine.LLLMatrix ×
      Generated.StrictRecombine.LLLMatrix × Array Nat × Array QQ)
    (hlll : Generated.StrictVanHoeij.lllReduceWithNorms concreteLLLExecution
      matrix hinput bound = .ok reduced)
    (candidates : Array (Array Int32))
    (hextract : Generated.StrictRecombine.extractCandidates
      (Generated.StrictVanHoeij.shortSpanFromNorms reduced.2.2.2 bound)
      reduced.2.1 output.1.size = .ok candidates) :
    ∃ indices : Array Nat,
      LegalCombination output.1.size indices.size indices ∧
      Polynomial.map (Int.castRingHom (ZMod modulus)) factor =
        Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          (indices.toList.map fun index => StrictHensel.toPolyMod modulus output.1[index]!).prod ∧
      ∀ classId leftValue rightValue,
        CandidateMember candidates classId leftValue →
        CandidateMember candidates classId rightValue →
        ∃ left right : Nat, left < output.1.size ∧ right < output.1.size ∧
          leftValue = left.toUInt32.toInt32 ∧ rightValue = right.toUInt32.toInt32 ∧
          (left ∈ indices.toList ↔ right ∈ indices.toList) := by
  obtain ⟨indices, coefficients, hlegal, hsupport, hcoordinates, hshort⟩ :=
    hensel_execution_truncated_short hcount hfactors hleadingSemantic hselection hentry
      original truncated modulus hmod divisor height scale target columns initial matrix
      hcanonical hnonempty hfits hliftedDegree hcldRun htruncate hscale hinit hbuild
      factor complement hfactor hm hd hh hdegree hheightRun
  have hsize := cldPolys_size_of_success ⟨termination⟩ f output.1 original modulus hcldRun
  have htruncSize := (truncateCld_spec original truncated divisor htruncate).2.1
  have hbasis := buildCldMatrix_canonical truncated 0 target scale
    ((modulus : Int) / divisor) initial (matrix, columns)
    (makeInitialMatrix_canonical truncated _ scale ((modulus : Int) / divisor)
      initial hscale hinit) hbuild
  have hfactorCount : output.1.size ≤ matrix.size := by
    have hs := hbasis.size_eq
    dsimp only at hs
    omega
  refine ⟨indices, hlegal, hsupport, ?_⟩
  intro classId leftValue rightValue hleft hright
  apply lll_class_support_membership matrix hinput bound output.1.size hfactorCount
    reduced hlll indices coefficients hcoordinates
    (by simpa only [hbound] using hshort) candidates hextract
    classId leftValue rightValue hleft hright

end Refinement.StrictRecombine
