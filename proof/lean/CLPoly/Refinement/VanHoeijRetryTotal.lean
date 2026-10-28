import CLPoly.Refinement.VanHoeijRetryCorrectness
import CLPoly.Refinement.VanHoeijRetryArithmetic
import CLPoly.Refinement.VanHoeijRetryPrecision
import CLPoly.Refinement.VanHoeijEntryTotal
import CLPoly.Refinement.VanHoeijHeuristicComplete

set_option autoImplicit false

namespace Refinement.StrictRecombine
open _root_.Generated.StrictVanHoeijRetry

/-- Actual certificate-only retry succeeds when a sufficient exponent fits
the machine interface. Hensel readiness is the existing callee contract;
no recombination run or output irreducibility is assumed. -/
theorem concrete_retryLoop_succeeds
    (this : DenseUPolyZp) (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hprime : this._p = selection.prime)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hphysical : StrictSelectPrime.SelectionPhysical selection)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hcount : 2 ≤ selection.factors.size) (hcountFits : selection.factors.size < 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : 0 < (SparsePolyZZ.toPoly f).natDegree)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hpositive : (0 : Int) < (SparsePolyZZ.toPoly f).leadingCoeff)
    (hruntime : ∀ target, 0 < target → StrictHensel.HenselLiftRuntimeReadiness this
      StrictHensel.concreteDivmodTermination mulProvider f selection.factors target)
    (height : Int)
    (hheight : height = ((SparsePolyZZ.toPoly f).natDegree *
      2 ^ ((SparsePolyZZ.toPoly f).natDegree - 1) *
      (Finset.range ((SparsePolyZZ.toPoly f).natDegree + 1)).sum
        (fun k => ((SparsePolyZZ.toPoly f).coeff k).natAbs) : Nat))
    (hfeasible : vanHoeijInputThreshold (SparsePolyZZ.toPoly f)
      (Int32.ofNat selection.factors.size) height <
        ((selection.prime.toNat ^ 2147483647 : Nat) : Int))
    (full precision : Int32) (hprecision : 0 < precision) :
    ∃ output, retryLoop (concreteVanHoeijRetryExecution this
      StrictHensel.concreteDivmodTermination mulProvider) f selection.factors selection.prime
      (Fact.out : Nat.Prime selection.prime.toNat).two_le full precision = .ok output := by
  let execution := concreteVanHoeijRetryExecution this StrictHensel.concreteDivmodTermination mulProvider
  let hp := (Fact.out : Nat.Prime selection.prime.toNat).two_le
  have hleading : ∀ leading, f[0]? = some leading →
      (leading.2 : ZMod selection.prime.toNat) = (SparsePolyZZ.toPoly f).leadingCoeff := by
    intro leading hl
    have he : f[0] = leading := by
      rw [Array.getElem?_eq_getElem hnonempty] at hl
      exact Option.some.inj hl
    rw [sparsePolyZZ_leadingCoeff_eq_head f hcanonical hnonempty, he]
  change ∃ output, retryLoop execution f selection.factors selection.prime hp full precision = .ok output
  induction hmeasure : (2147483647 - precision.toInt).toNat using Nat.strong_induction_on
      generalizing precision with
  | h measure ih =>
    obtain ⟨current, hl, hentry⟩ := concreteRetry_hensel_complete this
      StrictHensel.concreteDivmodTermination mulProvider hprime hcfg hphysical hselection
      hcanonical hnonempty hfits precision (hruntime precision hprecision) hcount hcountFits hprecision
    have hs := (StrictFactorZZ.selectionHenselFactors_pointwise_associated hcount hphysical
      hleading hselection hentry).1
    have hcFits : current.1.size < 2 ^ 31 := by omega
    obtain ⟨exponent, hexponent, hmod⟩ := hentry.outputModulus_eq_prime_pow
    obtain ⟨after, ha⟩ := hensel_tryRecombine_succeeds hcount hphysical hleading hselection hentry
      hcFits (selection.prime.toNat ^ exponent) hmod hcanonical hnonempty hdegree hfits
    have ha' : Generated.StrictVanHoeijEntry.tryRecombine execution.entry f current.1 current.2
        (execution.liftedCountFits f selection.factors selection.prime precision hp current hl) = .ok after := by
      simpa only [execution, concreteVanHoeijRetryExecution, hmod] using ha
    have hmaxComplete (he : precision = 2147483647) : after.complete = true := by
      have hpower := hensel_output_modulus_ge_requested hprecision hentry
      rw [he] at hpower
      change ((selection.prime.toNat ^ 2147483647 : Nat) : Int) ≤ current.2 at hpower
      have hthreshold : vanHoeijInputThreshold (SparsePolyZZ.toPoly f)
          (Int32.ofNat current.1.size) height < ((selection.prime.toNat ^ exponent : Nat) : Int) := by
        rw [← hs, ← hmod]
        exact hfeasible.trans_le hpower
      exact tryRecombine_above_input_threshold_complete hcount hphysical hleading hselection hentry
        hcFits (selection.prime.toNat ^ exponent) hmod hcanonical hnonempty hdegree hfits
        hprimitive hpositive height hheight hthreshold after ha
    rw [retryLoop]
    split
    next fault hbad => have he := hl.symm.trans hbad; contradiction
    next current' hl' =>
      have he := Except.ok.inj (hl'.symm.trans hl)
      subst current'
      split
      next fault hbad => have he := ha'.symm.trans hbad; contradiction
      next after' ha'' =>
        have he := Except.ok.inj (ha''.symm.trans ha')
        subst after'
        split
        · exact ⟨after.factors, rfl⟩
        next hnotComplete =>
          have hnotMax : precision ≠ 2147483647 := by
            intro he
            exact hnotComplete (hmaxComplete he)
          obtain ⟨next, hn⟩ := nextPrecision_succeeds precision full hnotMax
          have hg := nextPrecision_increases precision full next
            (Int32.le_iff_toInt_le.mpr (le_of_lt (Int32.lt_iff_toInt_lt.mp hprecision))) hn
          have hnextPositive : 0 < next := by
            apply Int32.lt_iff_toInt_lt.mpr
            have hp' := Int32.lt_iff_toInt_lt.mp hprecision
            change (0 : Int) < precision.toInt at hp'
            change (0 : Int) < next.toInt
            omega
          have hu := Int32.toInt_lt next
          obtain ⟨output, ho⟩ := ih _ (by omega) next hnextPositive rfl
          refine ⟨output, ?_⟩
          split
          next fault hbad => rw [hn] at hbad; contradiction
          next next' hn' =>
            have he := Except.ok.inj (hn'.symm.trans hn)
            subst next'
            exact ho

/-- End-to-end success and irreducibility after the actual heuristic prefix
has returned a positive machine exponent. The prefix's own execution and
signed-range obligations remain explicit here. -/
theorem concrete_lllFactorize_succeeds_of_heuristic
    (this : DenseUPolyZp) (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hprime : this._p = selection.prime)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hphysical : StrictSelectPrime.SelectionPhysical selection)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hcount : 2 ≤ selection.factors.size) (hcountFits : selection.factors.size < 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : 0 < (SparsePolyZZ.toPoly f).natDegree)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hpositive : (0 : Int) < (SparsePolyZZ.toPoly f).leadingCoeff)
    (hruntime : ∀ target, StrictHensel.HenselLiftRuntimeReadiness this
      StrictHensel.concreteDivmodTermination mulProvider f selection.factors target)
    (height : Int)
    (hheight : height = ((SparsePolyZZ.toPoly f).natDegree *
      2 ^ ((SparsePolyZZ.toPoly f).natDegree - 1) *
      (Finset.range ((SparsePolyZZ.toPoly f).natDegree + 1)).sum
        (fun k => ((SparsePolyZZ.toPoly f).coeff k).natAbs) : Nat))
    (hfeasible : vanHoeijInputThreshold (SparsePolyZZ.toPoly f)
      (Int32.ofNat selection.factors.size) height <
        ((selection.prime.toNat ^ 2147483647 : Nat) : Int))
    (full precision : Int32) (hprecision : 0 < precision)
    (hheuristic : Generated.StrictVanHoeijHeuristic.startingPrecision
      f (Int32.ofNat selection.factors.size) selection.prime = .ok (precision, full)) :
    ∃ output, lllFactorize (concreteVanHoeijRetryExecution this
      StrictHensel.concreteDivmodTermination mulProvider) f selection.factors selection.prime
      (Fact.out : Nat.Prime selection.prime.toNat).two_le = .ok output ∧
      FactorZZCorrect (SparsePolyZZ.toPoly f) (output.toList.map SparsePolyZZ.toPoly) := by
  obtain ⟨output, ho⟩ := concrete_retryLoop_succeeds this mulProvider hprime hcfg hphysical hselection
    hcount hcountFits hcanonical hnonempty hdegree hfits hprimitive hpositive
    (fun target _ => hruntime target) height hheight hfeasible full precision hprecision
  have hrun : lllFactorize (concreteVanHoeijRetryExecution this
      StrictHensel.concreteDivmodTermination mulProvider) f selection.factors selection.prime
      (Fact.out : Nat.Prime selection.prime.toNat).two_le = .ok output := by
    unfold lllFactorize
    simp only [hheuristic, bind, Except.bind]
    exact ho
  exact ⟨output, hrun, concreteVanHoeijRetry_refines_FactorZZCorrect this
    StrictHensel.concreteDivmodTermination mulProvider hprime hcfg hphysical hselection
    hcanonical hnonempty hfits hprimitive (fun target _ _ => hruntime target) output hrun⟩

/-- Actual checked heuristic prefix followed by the certificate-only retry.
Both representability premises depend only on the input. -/
theorem concrete_lllFactorize_succeeds
    (this : DenseUPolyZp) (mulProvider : StrictDDF.RawMulWorkspaceProvider this)
    {f : SparsePolyZZ} {selection : PrimeSelectionResult}
    [Fact (Nat.Prime selection.prime.toNat)]
    (hprime : this._p = selection.prime)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hphysical : StrictSelectPrime.SelectionPhysical selection)
    (hselection : StrictSelectPrime.SelectionCorrect (SparsePolyZZ.toPoly f) selection)
    (hcount : 2 ≤ selection.factors.size) (hcountFits : selection.factors.size < 2 ^ 31)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (hdegree : 0 < (SparsePolyZZ.toPoly f).natDegree)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hpositive : (0 : Int) < (SparsePolyZZ.toPoly f).leadingCoeff)
    (hruntime : ∀ target, StrictHensel.HenselLiftRuntimeReadiness this
      StrictHensel.concreteDivmodTermination mulProvider f selection.factors target)
    (height : Int)
    (hheight : height = ((SparsePolyZZ.toPoly f).natDegree *
      2 ^ ((SparsePolyZZ.toPoly f).natDegree - 1) *
      (Finset.range ((SparsePolyZZ.toPoly f).natDegree + 1)).sum
        (fun k => ((SparsePolyZZ.toPoly f).coeff k).natAbs) : Nat))
    (hfeasible : vanHoeijInputThreshold (SparsePolyZZ.toPoly f)
      (Int32.ofNat selection.factors.size) height <
        ((selection.prime.toNat ^ 2147483647 : Nat) : Int))
    (hrecoveryFeasible : heuristicRecoveryTarget f <
      (selection.prime.toNat : Int) ^ 2147483647) :
    ∃ output, lllFactorize (concreteVanHoeijRetryExecution this
      StrictHensel.concreteDivmodTermination mulProvider) f selection.factors selection.prime
      (Fact.out : Nat.Prime selection.prime.toNat).two_le = .ok output ∧
      FactorZZCorrect (SparsePolyZZ.toPoly f) (output.toList.map SparsePolyZZ.toPoly) := by
  obtain ⟨precision, full, hh, hp⟩ := startingPrecision_succeeds_positive f
    (Int32.ofNat selection.factors.size) selection.prime
    (Fact.out : Nat.Prime selection.prime.toNat).two_le hcanonical hnonempty hfits
    hpositive hrecoveryFeasible
  exact concrete_lllFactorize_succeeds_of_heuristic this mulProvider hprime hcfg
    hphysical hselection hcount hcountFits hcanonical hnonempty hdegree hfits
    hprimitive hpositive hruntime height hheight hfeasible full precision hp hh

end Refinement.StrictRecombine
