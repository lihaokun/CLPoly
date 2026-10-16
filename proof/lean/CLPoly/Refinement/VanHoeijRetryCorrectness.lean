import CLPoly.Refinement.VanHoeijRetryHensel
import CLPoly.Refinement.VanHoeijRetryTrace
import CLPoly.Refinement.VanHoeijEntryIrreducibility

set_option autoImplicit false

namespace Refinement.StrictRecombine

/-- Correctness of the new concrete retry entry on an actual successful run.
Source selection and recursive Hensel execution readiness remain explicit.
No returned factor's irreducibility, degree bound, count bound, or recovery
precision is assumed. This does not assert eventual pure van Hoeij success. -/
theorem concreteVanHoeijRetry_refines_FactorZZCorrect
    (this : DenseUPolyZp) (termination : Generated.StrictHensel.DivmodTermination)
    (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hprime : this._p = selection.prime)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hphysical : StrictSelectPrime.SelectionPhysical selection)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hdegree : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hruntime : ∀ target current,
      Generated.StrictHensel.__hensel_lift_upoly_raw_ir
        (StrictHensel.strictHenselRawOps termination)
        (StrictHensel.strictHenselTreeBuildRawOps this mulProvider)
        f selection.factors selection.prime target
        (Fact.out : Nat.Prime selection.prime.toNat).two_le = .ok current →
      StrictHensel.HenselLiftRuntimeReadiness this termination mulProvider f selection.factors target)
    (output : Array SparsePolyZZ)
    (hrun : Generated.StrictVanHoeijRetry.lllFactorize
      (concreteVanHoeijRetryExecution this termination mulProvider)
      f selection.factors selection.prime (Fact.out : Nat.Prime selection.prime.toNat).two_le = .ok output) :
    FactorZZCorrect (SparsePolyZZ.toPoly f) (output.toList.map SparsePolyZZ.toPoly) := by
  let execution := concreteVanHoeijRetryExecution this termination mulProvider
  let hp := (Fact.out : Nat.Prime selection.prime.toNat).two_le
  have hleadingSemantic : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff := by
    intro leading hleading
    have hfront : f[0] = leading := by
      rw [Array.getElem?_eq_getElem hnonempty] at hleading
      exact Option.some.inj hleading
    rw [sparsePolyZZ_leadingCoeff_eq_head f hcanonical hnonempty, hfront]
  have hheadDegree : f[0].1.deg < 2 ^ 63 := by
    rw [sparsePolyZZ_natDegree_eq_head f hcanonical hnonempty] at hdegree
    omega
  have hcert : ∀ target current,
      Generated.StrictHensel.__hensel_lift_upoly_raw_ir execution.step execution.tree
        f selection.factors selection.prime target hp = .ok current →
      StrictHensel.HenselLiftEntryCorrect termination f selection.factors selection.prime target current := by
    intro target current hl
    exact concreteRetry_hensel_success this termination mulProvider hprime hcfg
      hphysical hselection hcanonical hnonempty hdegree target (hruntime target current hl) current hl
  have horigin := lllFactorize_success_origin execution f selection.factors selection.prime hp output hrun
  obtain ⟨target, current, hl, after, hr, hc, hout⟩ := horigin
  have hentry := hcert target current hl
  have hcount := StrictHensel.__hensel_lift_upoly_raw_ir_factorCount_of_success
    execution.step execution.tree f selection.factors selection.prime target hp current hl
  have hcountFits := (StrictHensel.hensel_entry_success_count
    execution.step execution.tree f selection.factors selection.prime target hp current hl).2
  obtain ⟨exponent, hexponent, hmod⟩ := hentry.outputModulus_eq_prime_pow
  have hr' : Generated.StrictVanHoeijEntry.tryRecombine (concreteVanHoeijEntryExecution termination)
      f current.1 (selection.prime.toNat ^ exponent) hcountFits = .ok after := by
    simpa only [execution, concreteVanHoeijRetryExecution, hmod] using hr
  obtain ⟨hproduct, hirreducible⟩ := hensel_tryRecombine_complete_irreducible
    hcount hphysical hleadingSemantic hselection hentry (selection.prime.toNat ^ exponent)
    hmod hcanonical hnonempty hdegree hcountFits hprimitive after hc hr'
  rw [← hout]
  constructor
  · simpa only [factorArrayProduct] using Associated.of_eq hproduct.symm
  · intro factor hfactor
    obtain ⟨physical, hphysical, rfl⟩ := List.mem_map.mp hfactor
    exact hirreducible physical hphysical

end Refinement.StrictRecombine
