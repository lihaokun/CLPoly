import CLPoly.Generated.StrictVanHoeijFactorZZ
import CLPoly.Refinement.VanHoeijRetryCorrectness

set_option autoImplicit false

namespace Refinement.StrictRecombine

private theorem source_front_semantics (f : SparsePolyZZ)
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size) (p : UInt64) (_hp : Nat.Prime p.toNat) :
    ((SparsePolyZZ.front! f).2 : ZMod p.toNat) =
      ((SparsePolyZZ.toPoly f).leadingCoeff : ZMod p.toNat) := by
  rw [SparsePolyZZ.front!, getElem!_pos f 0 hnonempty,
    sparsePolyZZ_leadingCoeff_eq_head f hcanonical hnonempty]

/-- The retry operations are constructed at exactly the selected prime,
from the same physical candidate provider used by strict modular factoring. -/
noncomputable def concreteRetryAtPrime {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (p : UInt64) (hp : Nat.Prime p.toNat) : Generated.StrictVanHoeijRetry.RetryExecution :=
  let candidate := provider.physical p hp
  concreteVanHoeijRetryExecution candidate.dense StrictHensel.concreteDivmodTermination
    candidate.providers.mul

theorem concreteSelectedPrime {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (initialRng : State) (useLargePrime : Bool) (f : SparsePolyZZ)
    (hinitial : Nat.Prime (if useLargePrime then
      ((18446744073709551615 : UInt64) - 58).toNat else (2 : UInt64).toNat))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (selection : PrimeSelectionResult)
    (hrun : StrictFactorZZ.concreteSelectPrime engine provider initialRng f useLargePrime = .ok selection) :
    Nat.Prime selection.prime.toNat := by
  exact (StrictFactorZZ.concreteSelectPrime_success engine provider initialRng useLargePrime f
    hinitial hcanonical hnonempty hdegree (by omega)
    (source_front_semantics f hcanonical hnonempty) selection hrun).1.goodPrime.prime

/-- Executable control flow is the generated outer caller, with actual strict
selection and retry operations. Input proofs are erased before execution. -/
noncomputable def concreteFactorSquarefree {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (initialRng : State) (useLargePrime : Bool) (f : SparsePolyZZ)
    (hinitial : Nat.Prime (if useLargePrime then
      ((18446744073709551615 : UInt64) - 58).toNat else (2 : UInt64).toNat))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31) : RawExec (Array SparsePolyZZ) :=
  Generated.StrictVanHoeijFactorZZ.factorSquarefree
    (StrictFactorZZ.concreteSelectPrime engine provider initialRng)
    (concreteRetryAtPrime engine provider) useLargePrime f
    (concreteSelectedPrime engine provider initialRng useLargePrime f hinitial
      hcanonical hnonempty hdegree hfits)

/-- New outer caller correctness. Both the single-factor and retry branches
refer to the actual selected candidate and the same returned physical array. -/
theorem concreteFactorSquarefree_refines_FactorZZCorrect {State : Type}
    (engine : Generated.StrictEDF.RandomEngine State)
    (provider : StrictSelectPrime.CandidateRuntimeProvider engine)
    (initialRng : State) (useLargePrime : Bool) (f : SparsePolyZZ)
    (hinitial : Nat.Prime (if useLargePrime then
      ((18446744073709551615 : UInt64) - 58).toNat else (2 : UInt64).toNat))
    (hcanonical : StrictPolynomialMod.SparsePolyZZCanonical f)
    (hnonempty : 0 < f.size)
    (hdegree : 2 ≤ (SparsePolyZZ.toPoly f).natDegree)
    (hfits : (SparsePolyZZ.toPoly f).natDegree < 2 ^ 31)
    (hprimitive : (SparsePolyZZ.toPoly f).IsPrimitive)
    (hruntime : ∀ (selection : PrimeSelectionResult)
      (hp : Nat.Prime selection.prime.toNat) (target : Int32),
      let candidate := provider.physical selection.prime hp
      StrictHensel.HenselLiftRuntimeReadiness candidate.dense
        StrictHensel.concreteDivmodTermination candidate.providers.mul f selection.factors target)
    (output : Array SparsePolyZZ)
    (hrun : concreteFactorSquarefree engine provider initialRng useLargePrime f
      hinitial hcanonical hnonempty hdegree hfits = .ok output) :
    FactorZZCorrect (SparsePolyZZ.toPoly f) (output.toList.map SparsePolyZZ.toPoly) := by
  unfold concreteFactorSquarefree Generated.StrictVanHoeijFactorZZ.factorSquarefree at hrun
  split at hrun
  next hguard => contradiction
  next hguard =>
    split at hrun
    next fault hs => contradiction
    next selection hs =>
      have hselected := StrictFactorZZ.concreteSelectPrime_success
        engine provider initialRng useLargePrime f hinitial hcanonical hnonempty
        hdegree (by omega) (source_front_semantics f hcanonical hnonempty) selection hs
      let hp := hselected.1.goodPrime.prime
      letI selectionPrime : Fact (Nat.Prime selection.prime.toNat) := ⟨hp⟩
      split at hrun
      next hsingle =>
        have hout := Except.ok.inj hrun
        rw [← hout]
        have hcount : selection.factors.size ≤ 1 := by
          have hbranches : selection.irreducible = true ∨ selection.factors.size ≤ 1 := by
            simpa only [Bool.or_eq_true, decide_eq_true_eq] using hsingle
          rcases hbranches with hirred | hcount
          · exact StrictFactorZZ.concreteSelectPrime_irreducible_size engine provider
              initialRng useLargePrime f selection hs hirred
          · exact hcount
        simpa using StrictFactorZZ.selection_atMostOne_refines_singleton_FactorZZCorrect
          hselected.1 hprimitive hdegree hcount
      next hmultiple =>
        let candidate := provider.physical selection.prime hp
        apply concreteVanHoeijRetry_refines_FactorZZCorrect candidate.dense
          StrictHensel.concreteDivmodTermination candidate.providers.mul rfl
          candidate.configured hselected.2 hselected.1 hcanonical hnonempty hfits hprimitive
          (fun target current _ => hruntime selection hp target) output
        simpa only [concreteRetryAtPrime, candidate] using hrun

end Refinement.StrictRecombine
