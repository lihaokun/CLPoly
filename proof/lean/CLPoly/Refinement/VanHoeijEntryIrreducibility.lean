import CLPoly.Refinement.VanHoeijEntryTrace
import CLPoly.Refinement.VanHoeijSingletonExecution
import CLPoly.Refinement.VanHoeijColumnIrreducibility
import CLPoly.Refinement.VanHoeijHenselDegree

set_option autoImplicit false

namespace Refinement.StrictRecombine
open CLPoly.Math

theorem hensel_tryRecombine_complete_irreducible
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
    (modulus : Nat) (hmod : output.2 = (modulus : Int))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hcountFits : output.1.size < 2 ^ 31)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (after : Generated.StrictVanHoeij.VanHoeijResult)
    (hcomplete : after.complete = true)
    (hrun : Generated.StrictVanHoeijEntry.tryRecombine (concreteVanHoeijEntryExecution termination)
      f output.1 modulus hcountFits = .ok after) :
    factorArrayProduct after.factors = SparsePolyZZ.toPoly f ∧
      ∀ factor ∈ after.factors.toList, Irreducible (SparsePolyZZ.toPoly factor) := by
  obtain ⟨hguard, height, groups, hheight, hgroups, hbranch⟩ :=
    tryRecombine_complete_trace (concreteVanHoeijEntryExecution termination)
      f output.1 modulus hcountFits after hcomplete hrun
  have hcountSize : (Int32.ofNat output.1.size).toInt.toNat = output.1.size := by
    rw [Int32.toInt_ofNat_of_lt hcountFits]
    exact Int.toNat_natCast _
  rcases hbranch with hsingleton | hcolumns
  · exact hensel_constructed_singletons_validation_irreducible hcount hfactors hleadingSemantic
      hselection hentry modulus hmod (Int32.ofNat output.1.size) hcountSize groups hgroups
      hcanonical (by omega) hprimitive {} after hsingleton
  · obtain ⟨parameters, original, truncated, before, hp, hc, ht, hbefore,
      countSize, modulusNonzero, hloop⟩ := hcolumns
    obtain ⟨hh, hd, hdivisor, hm, hr, hlarge, hgap⟩ :=
      selectParameters_spec (Int32.ofNat output.1.size) height modulus parameters hp
    have hmodulusNonzero : (modulus : Int) / parameters.divisor ≠ 0 := by
      simpa only [hm] using modulusNonzero
    have hmodulus : 1 < modulus := by exact_mod_cast hguard.2
    have hnonnegative : 0 ≤ Int32.ofNat output.1.size := by
      apply Int32.le_iff_toInt_le.mpr
      rw [Int32.toInt_zero, Int32.toInt_ofNat_of_lt hcountFits]
      omega
    apply hensel_columnLoop_validation_irreducible hcount hfactors hleadingSemantic
      hselection hentry original truncated modulus hmod parameters.divisor height 32
      (Int32.ofNat output.1.size) (Generated.StrictVanHoeijParameters.cldWidth truncated)
      (min (2 : Int32) (Generated.StrictVanHoeijParameters.cldWidth truncated))
      hcanonical hnonempty hfits
      (hensel_output_head_degree_fits hcount hfactors hleadingSemantic hselection hentry hfits)
      hc ht (by decide) hmodulusNonzero
      hmodulus hd hh hheight hnonnegative countSize
      (initialCldTarget_bound _) hprimitive before after hbefore hcomplete
    simpa only [hm, hr, concreteVanHoeijEntryExecution] using hloop

end Refinement.StrictRecombine
