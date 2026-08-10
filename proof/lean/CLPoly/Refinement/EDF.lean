/-
  Strict EDF refinement is intentionally not exported yet.

  The former contents of this module selected an L2 existence witness when a
  bounded generated run failed or when its result did not satisfy
  `EDFCorrect`.  That construction was a certified specification fallback, not
  a refinement of the C++ `__edf_Zp` execution, and has therefore been removed.

  A future theorem in this module must start from a cpp2lean-generated EDF trace
  (including RNG state transitions), prove every successful split and recursive
  output invariant, and expose failure/nontermination explicitly.  It must not
  manufacture an L2 result on an unproved execution branch.
-/
import CLPoly.Algorithm.EDF
import CLPoly.Generated.StrictEDF
import CLPoly.Refinement.Basic
import CLPoly.Refinement.DDF
import CLPoly.Refinement.EDFRandom
import CLPoly.Refinement.EDFSubtractOne

set_option autoImplicit false

open Polynomial
open CLPoly.Math

namespace Refinement

namespace StrictEDF

/-- Concrete invariant carried by every recursive C++ EDF call.  It combines
the sparse representation obligations needed by raw operations with the L2
equal-degree hypotheses needed to prove the returned factor list correct. -/
structure EDFEntryInvariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (f : SparsePolyZp) (d : UInt64) : Prop where
  canonical : SparsePolyZp.Canonical this._p.toNat f
  primeMatches : 0 < f.size → f[0]!.2.prime = this._p
  degreeBound : ∀ term ∈ f.toList, term.1.deg < 2 ^ 62
  monic : (SparsePolyZp.toPoly this._p.toNat f).Monic
  degreePositive : 0 < (SparsePolyZp.toPoly this._p.toNat f).natDegree
  dPositive : 0 < d.toNat
  squarefree : Squarefree (SparsePolyZp.toPoly this._p.toNat f)
  equalDegree : ∀ q : Polynomial (ZMod this._p.toNat),
    Irreducible q → q ∣ SparsePolyZp.toPoly this._p.toNat f →
      q.natDegree = d.toNat

theorem EDFEntryInvariant.nonempty
    {this : DenseUPolyZp} [Fact (Nat.Prime this._p.toNat)]
    {f : SparsePolyZp} {d : UInt64} (h : EDFEntryInvariant this f d) :
    0 < f.size := by
  exact Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
    this._p.toNat f h.monic.ne_zero

theorem EDFEntryInvariant.natDegree_lt
    {this : DenseUPolyZp} [Fact (Nat.Prime this._p.toNat)]
    {f : SparsePolyZp} {d : UInt64} (h : EDFEntryInvariant this f d) :
    (SparsePolyZp.toPoly this._p.toNat f).natDegree < 2 ^ 63 := by
  exact lt_trans
    (StrictDDF.canonical_natDegree_lt_of_terms_lt this._p.toNat f
      h.canonical h.monic.ne_zero (2 ^ 62) h.degreeBound)
    (by norm_num)

/-- Decode the concrete C++ EDF accumulator without manufacturing or replacing
any factor. -/
noncomputable def edfResultToL2 (p : Nat) (result : Array SparsePolyZp) :
    List (Polynomial (ZMod p)) :=
  result.toList.map (SparsePolyZp.toPoly p)

@[simp] theorem edfResultToL2_empty (p : Nat) :
    edfResultToL2 p (#[] : Array SparsePolyZp) = [] := by
  simp [edfResultToL2]

@[simp] theorem edfResultToL2_push (p : Nat) (result : Array SparsePolyZp)
    (factor : SparsePolyZp) :
    edfResultToL2 p (result.push factor) =
      edfResultToL2 p result ++ [SparsePolyZp.toPoly p factor] := by
  simp [edfResultToL2]

private theorem certifyRawExec_ok {α : Type} (run : RawExec α) (output : α)
    (hrun : run = .ok output) :
    ∃ certified, Generated.StrictEDF.certifyRawExec run = .ok certified ∧
      certified.val = output := by
  cases run with
  | error fault => simp at hrun
  | ok value =>
      cases hrun
      exact ⟨⟨output, rfl⟩, rfl, rfl⟩

private theorem certifyRawExec_ok_eq {α : Type} (run : RawExec α)
    (output : α) (hrun : run = .ok output) :
    Generated.StrictEDF.certifyRawExec run = .ok ⟨output, hrun⟩ := by
  cases run with
  | error fault => simp at hrun
  | ok value =>
      cases hrun
      rfl

/-- Exact execution of the generated C++ base branch.  The theorem refers to
the generated state machine itself and to the concrete `makeMonic` run; it is
not an L2 execution fallback. -/
theorem rawState_base_run
    {State : Type} (ops : Generated.StrictEDF.EDFRawOps State)
    (termination : Generated.StrictEDF.EDFTermination ops)
    (state : Generated.StrictEDF.EDFState ops)
    (hdegree : ((get_deg state.f).toUInt64 == state.d) = true)
    (hmonic : ops.makeMonic state.f = .ok state.f) :
    Generated.StrictEDF.__edf_Zp_raw_ir_state ops termination state =
      .ok (state.result.push state.f, state.rng) := by
  rw [Generated.StrictEDF.__edf_Zp_raw_ir_state.eq_1]
  simp only [hdegree, ↓reduceIte]
  rcases certifyRawExec_ok _ _ hmonic with ⟨certified, hcertified, hvalue⟩
  rw [hcertified]
  simp [hvalue]

/-- L2 semantic fact used for the same C++ base branch: under the real EDF
preconditions, the factor appended by that branch is irreducible of degree
`d`. -/
theorem base_factor_correct
    {p : Nat} [Fact (Nat.Prime p)]
    (f : Polynomial (ZMod p)) (d : Nat)
    (hmonic : f.Monic) (hd : 0 < d) (hdegree : f.natDegree = d)
    (hfactors : ∀ q, Irreducible q → q ∣ f → q.natDegree = d) :
    EDFCorrect f d [f] := by
  have hirreducible := edf_base_irred f d hmonic hd hdegree hfactors
  exact ⟨by simp, by simp [hirreducible, hmonic, hdegree]⟩

/-- Exact raw-to-safe bridge for the odd-characteristic candidate branch.
The theorem requires, in source order, the actual powmod execution, the actual
generated subtract-one execution, and the actual GCD execution.  No
intermediate result can be replaced by an L2 witness. -/
theorem candidateRun_odd_run
    {State : Type} (ops : Generated.StrictEDF.EDFRawOps State)
    (f : SparsePolyZp) (d : UInt64) (r hpow hminus factor : SparsePolyZp)
    (hbudget : 0 < d.toNat ∧ d.toNat < UInt64.size)
    (hodd : (f[0]!.2.prime == 2) = false)
    (hpowRun : ops.powmod r
      (((f[0]!.2.prime.toNat : Int) ^ d.toNat - 1) / 2) f = .ok hpow)
    (hminusRun : Generated.StrictEDF.__upoly_subtract_one_raw_ir hpow
      f[0]!.2.prime = .ok hminus)
    (hgcdRun : ops.gcd hminus f = .ok factor) :
    Generated.StrictEDF.candidateRun ops f d r hbudget = .ok factor := by
  unfold Generated.StrictEDF.candidateRun
  simp only [hodd, Bool.false_eq_true, ↓reduceIte]
  rw [certifyRawExec_ok_eq _ _ hpowRun]
  simp only
  rw [certifyRawExec_ok_eq _ _ hminusRun]
  simp [hgcdRun]

/-- A finite exact retry trace is the raw-to-safe termination bridge for the
C++ `while (true)` loop.  The returned split retains the successful random
and candidate executions; the theorem does not synthesize a factor. -/
theorem retryLoop_terminates
    {State : Type} (ops : Generated.StrictEDF.EDFRawOps State)
    (f : SparsePolyZp) (d : UInt64) (rng : State)
    (trace : Generated.StrictEDF.RetryTrace ops f d rng) :
    ∃ split,
      Generated.StrictEDF.retryLoop ops f d rng trace = .ok split := by
  induction trace with
  | empty rng r rngNext randomRun isEmpty next ih =>
      simpa [Generated.StrictEDF.retryLoop] using ih
  | failed rng r rngNext candidate hbudget randomRun randomNonempty
      candidateExec notProper next ih =>
      simpa [Generated.StrictEDF.retryLoop] using ih
  | success rng r rngNext candidate hbudget randomRun randomNonempty
      candidateExec candidateNonempty proper =>
      exact ⟨⟨candidate, rng, rngNext, r, randomRun,
        ⟨hbudget, candidateExec⟩, candidateNonempty, proper⟩, rfl⟩

/-- Concrete odd-characteristic candidate pipeline assembled exclusively from
strict C++ raw boundaries already used by DDF, plus the generated EDF
subtract-one entry. -/
noncomputable def strictOddCandidateIR
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (f : SparsePolyZp) (d : UInt64) (r : SparsePolyZp) :
    RawExec SparsePolyZp := do
  let exponent := (f[0]!.2.prime.toNat ^ d.toNat - 1) / 2
  let hpow ← StrictDDF.strictPowmodIR this r exponent f providers.mul
    (providers.mod f)
  let hminus ← Generated.StrictEDF.__upoly_subtract_one_raw_ir hpow
    f[0]!.2.prime
  StrictDDF.strictDDFGCDIR this hminus f (providers.gcd hminus f)

/-- End-to-end semantics of the concrete odd-characteristic candidate
pipeline.  The returned candidate is the normalized GCD of the actual raw
powmod result minus one with the input polynomial. -/
theorem strictOddCandidateIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (f : SparsePolyZp) (d : UInt64) (r : SparsePolyZp)
    (hfCanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hrCanonical : SparsePolyZp.Canonical this._p.toNat r)
    (hfNonempty : 0 < f.size)
    (hfMonic : (SparsePolyZp.toPoly this._p.toNat f).Monic)
    (hfDegree : 0 < (SparsePolyZp.toPoly this._p.toNat f).natDegree) :
    ∃ factor,
      strictOddCandidateIR this providers f d r = .ok factor ∧
      SparsePolyZp.Canonical this._p.toNat factor ∧
      SparsePolyZp.toPoly this._p.toNat factor = normalize
        (EuclideanDomain.gcd
          ((SparsePolyZp.toPoly this._p.toNat r ^
              ((this._p.toNat ^ d.toNat - 1) / 2) %ₘ
                SparsePolyZp.toPoly this._p.toNat f) - 1)
          (SparsePolyZp.toPoly this._p.toNat f)) := by
  have hfHeadMem : f[0] ∈ f.toList := by
    simpa using Array.getElem_mem f 0 hfNonempty
  have hhead : f[0]! = f[0] := by
    rw [getElem!_def, getElem?_def]
    simp [hfNonempty]
  have hprimeNat : f[0]!.2.prime.toNat = this._p.toNat := by
    rw [hhead]
    exact (hfCanonical.1 f[0] hfHeadMem).1
  have hprime : f[0]!.2.prime = this._p := UInt64.toNat_inj.mp hprimeNat
  let exponent := (f[0]!.2.prime.toNat ^ d.toNat - 1) / 2
  rcases StrictDDF.strictPowmodIR_refines this providers.hcfg r exponent f
      providers.mul (providers.mod f) hrCanonical hfCanonical hfNonempty
      hfMonic hfDegree with
    ⟨hpow, hpowRun, hpowCanonical, hpowSemantic⟩
  rcases __upoly_subtract_one_raw_ir_certified hpow this._p
      providers.h2p hpowCanonical with
    ⟨hminus, hminusRun, hminusSemantic, hminusCanonical⟩
  have hminusRunHead :
      Generated.StrictEDF.__upoly_subtract_one_raw_ir hpow
        f[0]!.2.prime = .ok hminus := by
    simpa [hprime] using hminusRun
  rcases StrictDDF.strictDDFGCDIR_refines this hminus f
      (providers.gcd hminus f) hminusCanonical
      hfCanonical hfNonempty hfMonic with
    ⟨factor, hgcdRun, hfactorCanonical, hfactorSemantic⟩
  refine ⟨factor, ?_, hfactorCanonical, ?_⟩
  · unfold strictOddCandidateIR
    dsimp only
    rw [hpowRun]
    change (do
      let hminus ← Generated.StrictEDF.__upoly_subtract_one_raw_ir hpow
        f[0]!.2.prime
      StrictDDF.strictDDFGCDIR this hminus f
        (providers.gcd hminus f)) = .ok factor
    rw [hminusRunHead]
    simpa using hgcdRun
  · rw [hfactorSemantic, hminusSemantic, hpowSemantic]
    simp only [exponent, hprimeNat]

end StrictEDF

-- The public L1→L2 EDF theorem remains deliberately absent until the retry
-- and recursive split simulations below this base layer are closed.

end Refinement
