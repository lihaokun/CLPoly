import CLPoly.Refinement.VanHoeijExecution

set_option autoImplicit false

namespace Refinement.StrictRecombine
open CLPoly.Math

/-- The true-factor short-vector result for actual height, CLD, truncation and
matrix-construction runs. Physical machine-degree bounds remain explicit. -/
theorem hensel_execution_truncated_short
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
    (hheightRun : Generated.StrictVanHoeij.cldBound f = .ok height) :
    ∃ (indices : Array Nat) (coefficients : Fin matrix.size → ZZ),
      LegalCombination output.1.size indices.size indices ∧
      Polynomial.map (Int.castRingHom (ZMod modulus)) factor =
        Polynomial.C (factor.leadingCoeff : ZMod modulus) *
          (indices.toList.map fun index => StrictHensel.toPolyMod modulus output.1[index]!).prod ∧
      (∀ index, index.val < output.1.size →
        coefficients index = if index.val ∈ indices.toList then 1 else 0) ∧
      dotProduct
        (fun index => ((Matrix.vecMul coefficients
          (basisPrefixMatrix matrix matrix.size matrix.size) index : ZZ) : QQ))
        (fun index => ((Matrix.vecMul coefficients
          (basisPrefixMatrix matrix matrix.size matrix.size) index : ZZ) : QQ)) ≤
        ((output.1.size : ZZ) * scale * scale + (columns : ZZ) *
          ((height + divisor - 1) / divisor + 2 * (output.1.size : Int)) *
          ((height + divisor - 1) / divisor + 2 * (output.1.size : Int)) : ZZ) := by
  have hfBound : StrictHensel.DegreesBound f := by
    intro term ht
    have hhead := StrictHensel.canonical_headDominates f hnonempty hcanonical term ht
    rw [getElem!_pos f 0 hnonempty,
      ← sparsePolyZZ_natDegree_eq_head f hcanonical hnonempty] at hhead
    omega
  obtain ⟨hsize, hcld⟩ := hensel_cld_run_valid hentry modulus hm hfBound
    hliftedDegree original hcldRun
  have hbasis := buildCldMatrix_canonical truncated 0 target scale
    ((modulus : Int) / divisor) initial (matrix, columns)
    (makeInitialMatrix_canonical truncated _ scale ((modulus : Int) / divisor)
      initial hscale hinit) hbuild
  have hheight := cldBound_polynomial f height hcanonical hnonempty hfits hheightRun
  exact hensel_factor_truncated_short hcount hfactors hleadingSemantic hselection hentry
    original truncated modulus hmod divisor height scale
    (Generated.StrictRecombine.cldSpiralWidth truncated) columns matrix hsize hcld
    htruncate (by simpa only [Nat.zero_add] using hbasis) factor complement hfactor
    (by omega) hd hh hdegree hheight

end Refinement.StrictRecombine
