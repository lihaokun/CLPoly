import CLPoly.Refinement.VanHoeijColumnLoop
import CLPoly.Refinement.VanHoeijColumnBounds
import CLPoly.Refinement.VanHoeijIrreducibility

set_option autoImplicit false

namespace Refinement.StrictRecombine
open CLPoly.Math

theorem hensel_columnLoop_validation_irreducible
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
    (divisor height scale : ZZ) (count width target : Int32)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hliftedDegree : ∀ index (hi : index < output.1.size),
      output.1[index][0]!.1.deg < 2 ^ 63)
    (hcldRun : Generated.StrictRecombine.cldPolys ⟨termination⟩ f output.1 modulus = .ok original)
    (htruncate : Generated.StrictVanHoeij.truncateCld original divisor = .ok truncated)
    (hscale : scale ≠ 0)
    (hmodulusNonzero : (modulus : Int) / divisor ≠ 0)
    (hm : 1 < modulus) (hd : 0 < divisor) (hh : height < (modulus : Int))
    (hheightRun : Generated.StrictVanHoeij.cldBound f = .ok height)
    (hcountNonnegative : 0 ≤ count)
    (hcountSize : count.toInt.toNat = truncated.size)
    (htargetBound : target.toInt ≤ width.toInt)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (before after : Generated.StrictVanHoeij.VanHoeijResult)
    (hbefore : before.complete = false) (hafter : after.complete = true)
    (hrun : Generated.StrictVanHoeijController.columnLoop concreteVanHoeijControllerExecution
      f output.1 truncated modulus ((modulus : Int) / divisor) scale
      ((height + divisor - 1) / divisor + 2 * count.toInt) count width target
      hcountSize hscale hmodulusNonzero htargetBound before = .ok after) :
    factorArrayProduct after.factors = SparsePolyZZ.toPoly f ∧
      ∀ factor ∈ after.factors.toList, Irreducible (SparsePolyZZ.toPoly factor) := by
  obtain ⟨acceptedTarget, acceptedBefore, ht, htb, hr⟩ :=
    columnLoop_complete_origin concreteVanHoeijControllerExecution f output.1 truncated modulus
      ((modulus : Int) / divisor) scale ((height + divisor - 1) / divisor + 2 * count.toInt)
      count width target hcountSize hscale hmodulusNonzero htargetBound before after
      hbefore hafter hrun
  obtain ⟨initial, built, hi, hb, ready, reduced, candidates, hl, hc, hv⟩ :=
    columnRound_trace concreteVanHoeijControllerExecution f output.1 truncated modulus
      ((modulus : Int) / divisor) scale ((height + divisor - 1) / divisor + 2 * count.toInt)
      count acceptedTarget hcountSize hscale hmodulusNonzero acceptedBefore after true hr
  have hsize := cldPolys_size_of_success ⟨termination⟩ f output.1 original modulus hcldRun
  have htruncSize : truncated.size = output.1.size :=
    (truncateCld_spec original truncated divisor htruncate).2.1.trans hsize
  have hcountLifted : count.toInt.toNat = output.1.size := hcountSize.trans htruncSize
  have hn := Int32.le_iff_toInt_le.mp hcountNonnegative
  change (0 : Int) ≤ count.toInt at hn
  have hcountInt : count.toInt = (output.1.size : Int) := by omega
  have hindexFits : output.1.size ≤ 2 ^ 31 := by
    have hlimit := Int32.toInt_lt count
    omega
  have hadded := buildCldMatrix_added_int32 truncated 0 acceptedTarget
    ((modulus : Int) / divisor) initial built hb
  rw [hcountSize] at hi
  rw [hcountLifted] at hc
  exact hensel_lattice_validation_irreducible hcount hfactors hleadingSemantic hselection hentry
    hindexFits original truncated modulus hmod divisor height scale acceptedTarget.toInt.toNat
    built.2 initial built.1 hcanonical hnonempty hfits hliftedDegree hcldRun htruncate
    hscale hi hb hm hd hh hheightRun
    (count.toInt * scale * scale + (Int32.ofNat built.2).toInt *
      ((height + divisor - 1) / divisor + 2 * count.toInt) *
      ((height + divisor - 1) / divisor + 2 * count.toInt))
    (by rw [hcountInt, hadded])
    ready reduced hl candidates hc hprimitive _ after hv

end Refinement.StrictRecombine
