import CLPoly.Refinement.VanHoeijFullWidthRound
import CLPoly.Refinement.VanHoeijColumnLoopComplete

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- The actual positive-target column loop cannot normally return incomplete
above the input-fixed precision threshold. The full-width acceptance condition
is discharged from the same Hensel/CLD execution and source properties. -/
theorem columnLoop_above_input_threshold_complete
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
    (height : Int) (count width : Int32)
    (parameters : Generated.StrictVanHoeijParameters.Parameters)
    (hcountValue : count.toInt = (lifted.1.size : Int))
    (hselected : Generated.StrictVanHoeijParameters.selectParameters count height (modulus : Int) = some parameters)
    (htruncate : Generated.StrictVanHoeij.truncateCld original parameters.divisor = .ok truncated)
    (hwidth : width.toInt.toNat = Generated.StrictRecombine.cldSpiralWidth truncated)
    (countSize : count.toInt.toNat = truncated.size)
    (modulusNonzero : parameters.modulus ≠ 0)
    (hheight : height = ((SparsePolyZZ.toPoly f).natDegree *
      2 ^ ((SparsePolyZZ.toPoly f).natDegree - 1) *
      (Finset.range ((SparsePolyZZ.toPoly f).natDegree + 1)).sum
        (fun k => ((SparsePolyZZ.toPoly f).coeff k).natAbs) : Nat))
    (hsourceCanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hsourcePrimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hsourcePositive : (0 : Int) < (SparsePolyZZ.toPoly f).leadingCoeff)
    (hsourceFits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 63)
    (hprecision : vanHoeijInputThreshold (SparsePolyZZ.toPoly f) count height < (modulus : Int))
    (target : Int32) (targetBound : target.toInt ≤ width.toInt) (hpositive : 0 < target)
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hrun : Generated.StrictVanHoeijController.columnLoop concreteVanHoeijControllerExecution
      f lifted.1 truncated (modulus : Int) parameters.modulus 32 parameters.residual
      count width target countSize (by decide) modulusNonzero targetBound before = .ok after) :
    after.complete = true := by
  apply columnLoop_complete_of_full_width concreteVanHoeijControllerExecution f lifted.1
    truncated (modulus : Int) parameters.modulus 32 parameters.residual count width target
    countSize (by decide) modulusNonzero targetBound hpositive ?_ before after hrun
  intro first last hfail
  have h := columnRound_full_width_accepts hcount hfactors hleading hselection hentry
    hindexFits modulus hmodulus original truncated data height count width parameters
    hcountValue hselected htruncate hwidth countSize modulusNonzero hheight hsourceCanonical
    hsourcePrimitive hsourcePositive hsourceFits hprecision first last false hfail
  cases h

end Refinement.StrictRecombine
