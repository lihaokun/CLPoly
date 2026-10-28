/-
  Strict EDF refinement.

  This module proves the cpp2lean-generated recursive execution directly,
  including RNG transitions.  The former bounded safe wrapper and its L2
  existence-witness fallback remain removed; the public checked contract is
  generated in `CLPoly.Refinement.Generated`.
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

/-- Pure L2 facts needed after a successful concrete EDF split. -/
structure EDFPolynomialSplit {p : Nat} [Fact (Nat.Prime p)]
    (f g h : Polynomial (ZMod p)) (d : Nat) : Prop where
  product : g * h = f
  gMonic : g.Monic
  hMonic : h.Monic
  gSquarefree : Squarefree g
  hSquarefree : Squarefree h
  gDegreePositive : 0 < g.natDegree
  hDegreePositive : 0 < h.natDegree
  gDegreeLt : g.natDegree < f.natDegree
  hDegreeLt : h.natDegree < f.natDegree
  gEqualDegree : ∀ q, Irreducible q → q ∣ g → q.natDegree = d
  hEqualDegree : ∀ q, Irreducible q → q ∣ h → q.natDegree = d

/-- A proper monic divisor and the normalized exact quotient inherit every
mathematical invariant required by the two recursive EDF calls. -/
theorem edfPolynomialSplit_of_properDivisor
    {p : Nat} [Fact (Nat.Prime p)]
    (f g : Polynomial (ZMod p)) (d : Nat)
    (hfMonic : f.Monic) (hfSquarefree : Squarefree f)
    (hfEqualDegree : ∀ q, Irreducible q → q ∣ f → q.natDegree = d)
    (hgMonic : g.Monic) (hgDivides : g ∣ f)
    (hgPositive : 0 < g.natDegree)
    (hgProper : g.natDegree < f.natDegree) :
    EDFPolynomialSplit f g (normalize (f /ₘ g)) d := by
  let quotient := f /ₘ g
  have hmod : f %ₘ g = 0 := (modByMonic_eq_zero_iff_dvd hgMonic).mpr hgDivides
  have hproductRaw : g * quotient = f := by
    have hdivision := Polynomial.modByMonic_add_div f g
    rw [hmod, zero_add] at hdivision
    exact hdivision
  have hquotientDegree : quotient.natDegree = f.natDegree - g.natDegree :=
    Polynomial.natDegree_divByMonic f hgMonic
  have hquotientPositive : 0 < quotient.natDegree := by omega
  have hquotientNonzero : quotient ≠ 0 := by
    intro hzero
    rw [hzero] at hquotientPositive
    simp at hquotientPositive
  have hnormalizeMonic : (normalize quotient).Monic :=
    Polynomial.monic_normalize hquotientNonzero
  have hquotientMonic : quotient.Monic := by
    rw [Polynomial.Monic]
    have hleading := congrArg Polynomial.leadingCoeff hproductRaw
    rw [Polynomial.leadingCoeff_mul, hgMonic.leadingCoeff,
      hfMonic.leadingCoeff, one_mul] at hleading
    exact hleading
  have hnormalizeEq : normalize quotient = quotient :=
    hquotientMonic.normalize_eq_self
  have hnormalizeDegree : (normalize quotient).natDegree = quotient.natDegree := by
    rw [hnormalizeEq]
  have hnormalizeDivides : normalize quotient ∣ f := by
    rw [normalize_dvd_iff]
    exact ⟨g, by simpa [mul_comm] using hproductRaw.symm⟩
  have hnormalizeProduct : g * normalize quotient = f := by
    rw [hnormalizeEq]
    exact hproductRaw
  refine ⟨hnormalizeProduct, hgMonic, hnormalizeMonic,
    Squarefree.squarefree_of_dvd hgDivides hfSquarefree,
    Squarefree.squarefree_of_dvd hnormalizeDivides hfSquarefree,
    hgPositive, ?_, hgProper, ?_, ?_, ?_⟩
  · rw [hnormalizeDegree]
    exact hquotientPositive
  · rw [hnormalizeDegree, hquotientDegree]
    omega
  · intro q hq hqg
    exact hfEqualDegree q hq (dvd_trans hqg hgDivides)
  · intro q hq hqh
    exact hfEqualDegree q hq (dvd_trans hqh hnormalizeDivides)

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
not an L2 execution substitution. -/
theorem rawState_base_run
    {State : Type} (ops : Generated.StrictEDF.EDFRawOps State)
    (splitLaw : Generated.StrictEDF.EDFSplitLaw ops)
    (termination : Generated.StrictEDF.EDFTermination ops)
    (state : Generated.StrictEDF.EDFState ops)
    (hdegree : ((get_deg state.f).toUInt64 == state.d) = true)
    (hmonic : ops.makeMonic state.f = .ok state.f) :
    Generated.StrictEDF.__edf_Zp_raw_ir_state ops splitLaw termination state =
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

/-- Recursive-split consequences of the exact odd-characteristic candidate
pipeline.  The divisor is the actual raw GCD output returned above. -/
theorem strictOddCandidateIR_factor
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
      (SparsePolyZp.toPoly this._p.toNat factor).Monic ∧
      SparsePolyZp.toPoly this._p.toNat factor ∣
        SparsePolyZp.toPoly this._p.toNat f := by
  rcases strictOddCandidateIR_refines this providers f d r hfCanonical
      hrCanonical hfNonempty hfMonic hfDegree with
    ⟨factor, hrun, hcanonical, hsemantic⟩
  let candidateBase :=
    (SparsePolyZp.toPoly this._p.toNat r ^
        ((this._p.toNat ^ d.toNat - 1) / 2) %ₘ
          SparsePolyZp.toPoly this._p.toNat f) - 1
  let gcdResult := EuclideanDomain.gcd candidateBase
    (SparsePolyZp.toPoly this._p.toNat f)
  have hgcdDivides : gcdResult ∣ SparsePolyZp.toPoly this._p.toNat f :=
    EuclideanDomain.gcd_dvd_right candidateBase
      (SparsePolyZp.toPoly this._p.toNat f)
  have hgcdNonzero : gcdResult ≠ 0 := by
    intro hzero
    have := hgcdDivides
    rw [hzero, zero_dvd_iff] at this
    exact hfMonic.ne_zero this
  have hfactorMonic :
      (SparsePolyZp.toPoly this._p.toNat factor).Monic := by
    rw [hsemantic]
    exact Polynomial.monic_normalize hgcdNonzero
  have hfactorDivides : SparsePolyZp.toPoly this._p.toNat factor ∣
      SparsePolyZp.toPoly this._p.toNat f := by
    rw [hsemantic, normalize_dvd_iff]
    exact hgcdDivides
  exact ⟨factor, hrun, hcanonical, hfactorMonic, hfactorDivides⟩

/-- Exact characteristic-two trace-map iteration: execute the verified raw
multiplication, the sparse C++-level addition, then the verified raw modular
reduction, in source order. -/
def strictSquareAddModIR
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (g r f : SparsePolyZp) : RawExec SparsePolyZp := do
  let square ← StrictDDF.strictMulIR this g g providers.mul
  let sum := square + r
  StrictDDF.strictModIR this sum f ((providers.mod f).workspace sum)

/-- Semantic refinement of one concrete characteristic-two trace-map
iteration. -/
theorem strictSquareAddModIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (g r f : SparsePolyZp)
    (hgCanonical : SparsePolyZp.Canonical this._p.toNat g)
    (hrCanonical : SparsePolyZp.Canonical this._p.toNat r)
    (hfCanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hfNonempty : 0 < f.size)
    (hfMonic : (SparsePolyZp.toPoly this._p.toNat f).Monic) :
    ∃ next,
      strictSquareAddModIR this providers g r f = .ok next ∧
      SparsePolyZp.Canonical this._p.toNat next ∧
      SparsePolyZp.toPoly this._p.toNat next =
        (SparsePolyZp.toPoly this._p.toNat g ^ 2 +
          SparsePolyZp.toPoly this._p.toNat r) %ₘ
            SparsePolyZp.toPoly this._p.toNat f := by
  rcases StrictDDF.strictMulIR_refines_mul this providers.hcfg g g
      providers.mul hgCanonical hgCanonical with
    ⟨square, hsquareRun, hsquareCanonical, hsquareSemantic⟩
  let sum := square + r
  have hsumCanonical : SparsePolyZp.Canonical this._p.toNat sum :=
    SparsePolyZp.Canonical.add this._p.toNat square r
      hsquareCanonical hrCanonical
  have hsumSemantic : SparsePolyZp.toPoly this._p.toNat sum =
      SparsePolyZp.toPoly this._p.toNat square +
        SparsePolyZp.toPoly this._p.toNat r := by
    exact SparsePolyZp.toPoly_add this._p.toNat
      (Fact.out : Nat.Prime this._p.toNat).pos square r
      hsquareCanonical.1 hrCanonical.1
  rcases StrictDDF.strictModIR_refines_modByMonic this providers.hcfg sum f
      ((providers.mod f).workspace sum) hsumCanonical hfCanonical hfNonempty
      hfMonic with
    ⟨next, hmodRun, hnextCanonical, hnextSemantic⟩
  refine ⟨next, ?_, hnextCanonical, ?_⟩
  · unfold strictSquareAddModIR
    rw [hsquareRun]
    exact hmodRun
  · rw [hnextSemantic, hsumSemantic, hsquareSemantic]
    ring_nf

/-- The generated characteristic-two trace loop terminates with a canonical
array whenever its concrete square/add/mod operation does.  The induction is
on the source loop distance `d-i`. -/
theorem traceLoop_terminates_canonical
    {State : Type} (ops : Generated.StrictEDF.EDFRawOps State)
    (d : UInt64) (f r : SparsePolyZp) (i : UInt64) (g : SparsePolyZp)
    (hbudget : i.toNat ≤ d.toNat ∧ d.toNat < UInt64.size)
    (hgCanonical : SparsePolyZp.Canonical f[0]!.2.prime.toNat g)
    (hsquare : ∀ current,
      SparsePolyZp.Canonical f[0]!.2.prime.toNat current →
      ∃ next,
        ops.squareAddMod current r f = .ok next ∧
        SparsePolyZp.Canonical f[0]!.2.prime.toNat next) :
    ∃ trace,
      Generated.StrictEDF.traceLoop ops d f r i g hbudget = .ok trace ∧
      SparsePolyZp.Canonical f[0]!.2.prime.toNat trace := by
  generalize hdistance : d.toNat - i.toNat = distance
  induction distance using Nat.strong_induction_on generalizing i g with
  | h distance ih =>
      rw [Generated.StrictEDF.traceLoop]
      split
      next hmore =>
        rcases hsquare g hgCanonical with ⟨next, hnextRun, hnextCanonical⟩
        rw [certifyRawExec_ok_eq _ _ hnextRun]
        apply ih (d.toNat - (i + 1).toNat)
        · have hi : i.toNat < d.toNat := by simpa using hmore
          have hiSize : i.toNat + 1 < UInt64.size :=
            Nat.lt_of_le_of_lt (Nat.succ_le_of_lt hi) hbudget.2
          simp [UInt64.toNat_add, Nat.mod_eq_of_lt hiSize]
          rw [← hdistance]
          omega
        · exact hnextCanonical
        · rfl
      next hdone =>
        exact ⟨g, rfl, hgCanonical⟩

/-- The unique concrete operation record used by the strict generated EDF
shell.  Every computational field is an executable L1 boundary; proof
obligations live separately in `EDFSplitLaw`. -/
noncomputable def strictEDFRawOps
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this) :
    Generated.StrictEDF.EDFRawOps State where
  random := Generated.StrictEDF.__upoly_random_raw_ir engine
  modPoly := fun dividend modulus =>
    StrictDDF.strictModIR this dividend modulus
      ((providers.mod modulus).workspace dividend)
  squareAddMod := strictSquareAddModIR this providers
  powmod := fun base exponent modulus =>
    StrictDDF.strictPowmodIR this base exponent.toNat modulus providers.mul
      (providers.mod modulus)
  gcd := fun left right =>
    StrictDDF.strictDDFGCDIR this left right (providers.gcd left right)
  exactDiv := StrictDDF.strictExactDivIR this
  makeMonic := StrictDDF.strictMakeMonicIR this
  EntryInvariant := EDFEntryInvariant this

/-- The generated characteristic-two candidate branch terminates through the
actual mod/trace/GCD calls and returns a canonical monic divisor of the input.
The factor is the concrete raw GCD output; no L2 witness is substituted. -/
theorem strictCandidateRun_charTwo_factor
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (f : SparsePolyZp) (d : UInt64) (r : SparsePolyZp)
    (hinvariant : EDFEntryInvariant this f d)
    (hrCanonical : SparsePolyZp.Canonical this._p.toNat r)
    (hbudget : 0 < d.toNat ∧ d.toNat < UInt64.size)
    (hcharTwo : (f[0]!.2.prime == 2) = true) :
    ∃ factor,
      Generated.StrictEDF.candidateRun
          (strictEDFRawOps engine this providers) f d r hbudget = .ok factor ∧
      SparsePolyZp.Canonical this._p.toNat factor ∧
      (SparsePolyZp.toPoly this._p.toNat factor).Monic ∧
      SparsePolyZp.toPoly this._p.toNat factor ∣
        SparsePolyZp.toPoly this._p.toNat f := by
  have hfNonempty := hinvariant.nonempty
  have hprime : f[0]!.2.prime = this._p :=
    hinvariant.primeMatches hfNonempty
  rcases StrictDDF.strictModIR_refines_modByMonic this providers.hcfg r f
      ((providers.mod f).workspace r) hrCanonical hinvariant.canonical
      hfNonempty hinvariant.monic with
    ⟨g0, hg0Run, hg0Canonical, _hg0Semantic⟩
  have hsquare : ∀ current,
      SparsePolyZp.Canonical f[0]!.2.prime.toNat current →
      ∃ next,
        (strictEDFRawOps engine this providers).squareAddMod current r f =
            .ok next ∧
        SparsePolyZp.Canonical f[0]!.2.prime.toNat next := by
    intro current hcurrent
    rcases strictSquareAddModIR_refines this providers current r f
        (by simpa [hprime] using hcurrent) hrCanonical hinvariant.canonical
        hfNonempty hinvariant.monic with
      ⟨next, hnextRun, hnextCanonical, _hnextSemantic⟩
    exact ⟨next, hnextRun, by simpa [hprime] using hnextCanonical⟩
  rcases traceLoop_terminates_canonical
      (strictEDFRawOps engine this providers) d f r 1 g0
      ⟨by simpa using hbudget.1, hbudget.2⟩
      (by simpa [hprime] using hg0Canonical) hsquare with
    ⟨trace, htraceRun, htraceCanonical⟩
  rcases StrictDDF.strictDDFGCDIR_refines this trace f
      (providers.gcd trace f) (by simpa [hprime] using htraceCanonical)
      hinvariant.canonical hfNonempty hinvariant.monic with
    ⟨factor, hgcdRun, hfactorCanonical, hfactorSemantic⟩
  have hgcdDivides :
      EuclideanDomain.gcd
          (SparsePolyZp.toPoly this._p.toNat trace)
          (SparsePolyZp.toPoly this._p.toNat f) ∣
        SparsePolyZp.toPoly this._p.toNat f :=
    EuclideanDomain.gcd_dvd_right _ _
  have hgcdNonzero :
      EuclideanDomain.gcd
          (SparsePolyZp.toPoly this._p.toNat trace)
          (SparsePolyZp.toPoly this._p.toNat f) ≠ 0 := by
    intro hzero
    rw [hzero, zero_dvd_iff] at hgcdDivides
    exact hinvariant.monic.ne_zero hgcdDivides
  refine ⟨factor, ?_, hfactorCanonical, ?_, ?_⟩
  · unfold Generated.StrictEDF.candidateRun
    simp only [hcharTwo, ↓reduceIte]
    have hg0Ops :
        (strictEDFRawOps engine this providers).modPoly r f = .ok g0 := by
      simpa [strictEDFRawOps] using hg0Run
    rw [certifyRawExec_ok_eq _ _ hg0Ops]
    simp only
    rw [htraceRun]
    simpa [strictEDFRawOps] using hgcdRun
  · rw [hfactorSemantic]
    exact Polynomial.monic_normalize hgcdNonzero
  · rw [hfactorSemantic, normalize_dvd_iff]
    exact hgcdDivides

/-- The generated odd-characteristic branch is exactly the previously
certified strict powmod/subtract-one/GCD pipeline.  The only representation
bridge is the equality between the source `Int` exponent converted by the raw
operation record and the corresponding nonnegative `Nat` exponent. -/
theorem strictCandidateRun_odd_factor
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (f : SparsePolyZp) (d : UInt64) (r : SparsePolyZp)
    (hinvariant : EDFEntryInvariant this f d)
    (hrCanonical : SparsePolyZp.Canonical this._p.toNat r)
    (hbudget : 0 < d.toNat ∧ d.toNat < UInt64.size)
    (hodd : (f[0]!.2.prime == 2) = false) :
    ∃ factor,
      Generated.StrictEDF.candidateRun
          (strictEDFRawOps engine this providers) f d r hbudget = .ok factor ∧
      SparsePolyZp.Canonical this._p.toNat factor ∧
      (SparsePolyZp.toPoly this._p.toNat factor).Monic ∧
      SparsePolyZp.toPoly this._p.toNat factor ∣
        SparsePolyZp.toPoly this._p.toNat f := by
  have hfNonempty := hinvariant.nonempty
  have hprime : f[0]!.2.prime = this._p :=
    hinvariant.primeMatches hfNonempty
  have honeLePow : 1 ≤ this._p.toNat ^ d.toNat :=
    Nat.one_le_pow d.toNat this._p.toNat
      ((Fact.out : Nat.Prime this._p.toNat).pos)
  have hexponentInt :
      (((this._p.toNat : Int) ^ d.toNat - 1) / 2) =
        (((this._p.toNat ^ d.toNat - 1) / 2 : Nat) : Int) := by
    rw [← Int.natCast_pow]
    calc
      ((this._p.toNat ^ d.toNat : Int) - 1) / 2 =
          ((this._p.toNat ^ d.toNat - 1 : Nat) : Int) / 2 := by
            rw [Int.ofNat_sub honeLePow]
            norm_num
      _ = (((this._p.toNat ^ d.toNat - 1) / 2 : Nat) : Int) :=
        (Int.natCast_ediv _ _).symm
  have hexponentToNat :
      ((((f[0]!.2.prime.toNat : Int) ^ d.toNat - 1) / 2).toNat) =
        (f[0]!.2.prime.toNat ^ d.toNat - 1) / 2 := by
    rw [hprime]
    rw [hexponentInt]
    exact Int.toNat_natCast _
  rcases strictOddCandidateIR_factor this providers f d r
      hinvariant.canonical hrCanonical hfNonempty hinvariant.monic
      hinvariant.degreePositive with
    ⟨factor, hrun, hfactorCanonical, hfactorMonic, hfactorDivides⟩
  unfold strictOddCandidateIR at hrun
  dsimp only at hrun
  generalize hpowRun :
      StrictDDF.strictPowmodIR this r
          ((f[0]!.2.prime.toNat ^ d.toNat - 1) / 2) f providers.mul
          (providers.mod f) = powRun at hrun
  cases powRun with
  | error fault =>
      simp only [Bind.bind, Except.bind] at hrun
      contradiction
  | ok hpow =>
      simp only [Bind.bind, Except.bind] at hrun
      generalize hminusRun :
          Generated.StrictEDF.__upoly_subtract_one_raw_ir hpow
            f[0]!.2.prime = minusRun at hrun
      cases minusRun with
      | error fault =>
          contradiction
      | ok hminus =>
          have hpowOps :
              (strictEDFRawOps engine this providers).powmod r
                  (((f[0]!.2.prime.toNat : Int) ^ d.toNat - 1) / 2) f =
                .ok hpow := by
            simpa [strictEDFRawOps, hexponentToNat] using hpowRun
          have hgcdOps :
              (strictEDFRawOps engine this providers).gcd hminus f =
                .ok factor := by
            simpa [strictEDFRawOps] using hrun
          refine ⟨factor, ?_, hfactorCanonical, hfactorMonic,
            hfactorDivides⟩
          exact candidateRun_odd_run
            (strictEDFRawOps engine this providers) f d r hpow hminus factor
            hbudget hodd hpowOps hminusRun hgcdOps

/-- Uniform correctness of the exact generated candidate execution.  The
proof follows the source characteristic test and delegates to the two strict
raw branches above. -/
theorem strictCandidateRun_factor
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (f : SparsePolyZp) (d : UInt64) (r : SparsePolyZp)
    (hinvariant : EDFEntryInvariant this f d)
    (hrCanonical : SparsePolyZp.Canonical this._p.toNat r)
    (hbudget : 0 < d.toNat ∧ d.toNat < UInt64.size) :
    ∃ factor,
      Generated.StrictEDF.candidateRun
          (strictEDFRawOps engine this providers) f d r hbudget = .ok factor ∧
      SparsePolyZp.Canonical this._p.toNat factor ∧
      (SparsePolyZp.toPoly this._p.toNat factor).Monic ∧
      SparsePolyZp.toPoly this._p.toNat factor ∣
        SparsePolyZp.toPoly this._p.toNat f := by
  cases hcharacteristic : (f[0]!.2.prime == 2) with
  | false =>
      exact strictCandidateRun_odd_factor engine this providers f d r
        hinvariant hrCanonical hbudget hcharacteristic
  | true =>
      exact strictCandidateRun_charTwo_factor engine this providers f d r
        hinvariant hrCanonical hbudget hcharacteristic

/-- Invert the successful random and candidate executions retained by a retry
trace.  This identifies the specific returned candidate with the strict branch
certificate, rather than merely proving that some candidate exists. -/
theorem strictSuccessfulCandidate_factor
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (f : SparsePolyZp) (d : UInt64)
    (rngBefore rngAfter : State) (r g : SparsePolyZp)
    (hinvariant : EDFEntryInvariant this f d)
    (hrandom : (strictEDFRawOps engine this providers).random
      (get_deg f) f[0]!.2.prime rngBefore = .ok (r, rngAfter))
    (hbudget : 0 < d.toNat ∧ d.toNat < UInt64.size)
    (hcandidate : Generated.StrictEDF.candidateRun
      (strictEDFRawOps engine this providers) f d r hbudget = .ok g) :
    SparsePolyZp.Canonical this._p.toNat g ∧
      (SparsePolyZp.toPoly this._p.toNat g).Monic ∧
      SparsePolyZp.toPoly this._p.toNat g ∣
        SparsePolyZp.toPoly this._p.toNat f := by
  have hfNonempty := hinvariant.nonempty
  have hprime : f[0]!.2.prime = this._p :=
    hinvariant.primeMatches hfNonempty
  have hpWord : 0 < f[0]!.2.prime := by
    simpa [hprime, UInt64.lt_iff_toNat_lt] using
      (Fact.out : Nat.Prime this._p.toNat).pos
  have hrandomStrict :
      Generated.StrictEDF.__upoly_random_raw_ir engine (get_deg f)
        f[0]!.2.prime rngBefore = .ok (r, rngAfter) := by
    simpa [strictEDFRawOps] using hrandom
  rcases __upoly_random_raw_ir_canonical engine (get_deg f)
      f[0]!.2.prime rngBefore hpWord with
    ⟨randomOutput, randomNext, hrandomCanonicalRun,
      hrandomCanonical, _hrandomDegree⟩
  rw [hrandomStrict] at hrandomCanonicalRun
  have hrandomPair := Except.ok.inj hrandomCanonicalRun
  have hrandomOutputEq : r = randomOutput := congrArg Prod.fst hrandomPair
  subst randomOutput
  have hrCanonical : SparsePolyZp.Canonical this._p.toNat r := by
    simpa [hprime] using hrandomCanonical
  rcases strictCandidateRun_factor engine this providers f d r hinvariant
      hrCanonical hbudget with
    ⟨candidate, hcandidateStrict, hcanonical, hmonic, hdivides⟩
  rw [hcandidate] at hcandidateStrict
  injection hcandidateStrict with hcandidateEq
  subst candidate
  exact ⟨hcanonical, hmonic, hdivides⟩

/-- The exact division and two make-monic calls used by a successful split
preserve the concrete factor and quotient, whose L2 product is the input. -/
theorem strictSplitCalls_product
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (f : SparsePolyZp) (d : UInt64)
    (rngBefore rngAfter : State) (r g hRaw gMonic hMonic : SparsePolyZp)
    (hinvariant : EDFEntryInvariant this f d)
    (hrandom : (strictEDFRawOps engine this providers).random
      (get_deg f) f[0]!.2.prime rngBefore = .ok (r, rngAfter))
    (hbudget : 0 < d.toNat ∧ d.toNat < UInt64.size)
    (hcandidate : Generated.StrictEDF.candidateRun
      (strictEDFRawOps engine this providers) f d r hbudget = .ok g)
    (hdivRun : (strictEDFRawOps engine this providers).exactDiv f g =
      .ok hRaw)
    (hgMonicRun : (strictEDFRawOps engine this providers).makeMonic g =
      .ok gMonic)
    (hhMonicRun : (strictEDFRawOps engine this providers).makeMonic
      hRaw.normalization = .ok hMonic) :
    SparsePolyZp.toPoly this._p.toNat gMonic *
        SparsePolyZp.toPoly this._p.toNat hMonic =
      SparsePolyZp.toPoly this._p.toNat f := by
  rcases strictSuccessfulCandidate_factor engine this providers f d rngBefore
      rngAfter r g hinvariant hrandom hbudget hcandidate with
    ⟨hgCanonical, hgMonic, hgDivides⟩
  have hgNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat g hgMonic.ne_zero
  rcases StrictDDF.strictExactDivIR_refines this f g providers.hcfg
      hinvariant.canonical hgCanonical hgNonempty hgMonic hgDivides with
    ⟨quotient, hquotientRun, hquotientCanonical, hquotientSemantic⟩
  have hdivRunStrict : StrictDDF.strictExactDivIR this f g = .ok hRaw := by
    simpa [strictEDFRawOps] using hdivRun
  rw [hdivRunStrict] at hquotientRun
  injection hquotientRun with hquotientEq
  subst quotient
  have hquotientMonic :
      (SparsePolyZp.toPoly this._p.toNat hRaw).Monic := by
    rw [hquotientSemantic]
    exact StrictSquarefreeZp.divByMonic_monic_of_monic_of_dvd
      (SparsePolyZp.toPoly this._p.toNat f)
      (SparsePolyZp.toPoly this._p.toNat g) hinvariant.monic hgMonic
      hgDivides
  have hhNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat hRaw hquotientMonic.ne_zero
  have hgMonicExact := StrictDDF.strictMakeMonicIR_eq_of_monic this g
    hgCanonical hgNonempty hgMonic
  have hgMonicRunStrict : StrictDDF.strictMakeMonicIR this g =
      .ok gMonic := by
    simpa [strictEDFRawOps] using hgMonicRun
  rw [hgMonicExact] at hgMonicRunStrict
  injection hgMonicRunStrict with hgMonicEq
  subst gMonic
  have hnormalizeRaw : hRaw.normalization = hRaw :=
    StrictSquarefreeZp.sparsePolyZp_normalization_eq_of_canonical
      this._p.toNat hRaw hquotientCanonical
  have hhMonicExact := StrictDDF.strictMakeMonicIR_eq_of_monic this hRaw
    hquotientCanonical hhNonempty hquotientMonic
  have hhMonicRunStrict : StrictDDF.strictMakeMonicIR this
      hRaw.normalization = .ok hMonic := by
    simpa [strictEDFRawOps] using hhMonicRun
  rw [hnormalizeRaw, hhMonicExact] at hhMonicRunStrict
  injection hhMonicRunStrict with hhMonicEq
  subst hMonic
  rw [hquotientSemantic]
  have hmod : SparsePolyZp.toPoly this._p.toNat f %ₘ
      SparsePolyZp.toPoly this._p.toNat g = 0 :=
    (Polynomial.modByMonic_eq_zero_iff_dvd hgMonic).mpr hgDivides
  have hdivision := Polynomial.modByMonic_add_div
    (SparsePolyZp.toPoly this._p.toNat f)
    (SparsePolyZp.toPoly this._p.toNat g)
  rw [hmod, zero_add] at hdivision
  exact hdivision

/-- All raw calls following a successful retry terminate on the concrete
strict implementations.  The witnesses are precisely the exact quotient and
the already-monic factor/quotient consumed by the generated state machine. -/
theorem strictSplitCalls_terminate
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (f : SparsePolyZp) (d : UInt64)
    (split : Generated.StrictEDF.SplitState
      (strictEDFRawOps engine this providers) f d)
    (hinvariant : EDFEntryInvariant this f d) :
    ∃ hRaw gMonic hMonic,
      (strictEDFRawOps engine this providers).exactDiv f split.factor =
          .ok hRaw ∧
      (strictEDFRawOps engine this providers).makeMonic split.factor =
          .ok gMonic ∧
      (strictEDFRawOps engine this providers).makeMonic
          hRaw.normalization = .ok hMonic := by
  rcases split.candidateRun with ⟨hbudget, hcandidate⟩
  rcases strictSuccessfulCandidate_factor engine this providers f d
      split.rngBefore split.rng split.randomPoly split.factor hinvariant
      split.randomRun hbudget hcandidate with
    ⟨hgCanonical, hgMonic, hgDivides⟩
  have hgNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat split.factor hgMonic.ne_zero
  rcases StrictDDF.strictExactDivIR_refines this f split.factor providers.hcfg
      hinvariant.canonical hgCanonical hgNonempty hgMonic hgDivides with
    ⟨hRaw, hdivRun, hhCanonical, hhSemantic⟩
  have hhMonic : (SparsePolyZp.toPoly this._p.toNat hRaw).Monic := by
    rw [hhSemantic]
    exact StrictSquarefreeZp.divByMonic_monic_of_monic_of_dvd
      (SparsePolyZp.toPoly this._p.toNat f)
      (SparsePolyZp.toPoly this._p.toNat split.factor)
      hinvariant.monic hgMonic hgDivides
  have hhNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat hRaw hhMonic.ne_zero
  have hgMonicRun := StrictDDF.strictMakeMonicIR_eq_of_monic this
    split.factor hgCanonical hgNonempty hgMonic
  have hnormalizeRaw : hRaw.normalization = hRaw :=
    StrictSquarefreeZp.sparsePolyZp_normalization_eq_of_canonical
      this._p.toNat hRaw hhCanonical
  have hhMonicRun := StrictDDF.strictMakeMonicIR_eq_of_monic this hRaw
    hhCanonical hhNonempty hhMonic
  refine ⟨hRaw, split.factor, hRaw, ?_, ?_, ?_⟩
  · simpa [strictEDFRawOps] using hdivRun
  · simpa [strictEDFRawOps] using hgMonicRun
  · simpa [strictEDFRawOps, hnormalizeRaw] using hhMonicRun

/-- On canonical nonempty inputs, the generated recursive measure is exactly
one more than the L2 natural degree. -/
theorem edfMeasure_eq_natDegree_succ
    (p : Nat) [Fact (Nat.Prime p)] (f : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical p f) (hnonempty : 0 < f.size) :
    Generated.StrictEDF.edfMeasure f =
      (SparsePolyZp.toPoly p f).natDegree + 1 := by
  have hdegree :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
      p f hcanonical hnonempty
  have hnatDegree : (SparsePolyZp.toPoly p f).natDegree = f[0].1.deg :=
    Polynomial.natDegree_eq_of_degree_eq_some hdegree
  have hnotEmpty : f.isEmpty = false := by
    simp [Array.isEmpty, Nat.ne_of_gt hnonempty]
  simp [Generated.StrictEDF.edfMeasure, hnotEmpty,
    getElem!_pos f 0 hnonempty, hnatDegree]

/-- A positive signed C++ degree comparison is an exact L2 natural-degree
comparison when both canonical degrees fit the signed word range. -/
theorem strict_get_deg_lt_implies_natDegree_lt
    (p : Nat) [Fact (Nat.Prime p)] (g f : SparsePolyZp)
    (hgCanonical : SparsePolyZp.Canonical p g)
    (hfCanonical : SparsePolyZp.Canonical p f)
    (hgNonempty : 0 < g.size) (hfNonempty : 0 < f.size)
    (hgBound : (SparsePolyZp.toPoly p g).natDegree < 2 ^ 63)
    (hfBound : (SparsePolyZp.toPoly p f).natDegree < 2 ^ 63)
    (hgPositive : get_deg g > 0) (hlt : get_deg g < get_deg f) :
    (SparsePolyZp.toPoly p g).natDegree <
      (SparsePolyZp.toPoly p f).natDegree := by
  have hltInt : (get_deg g).toInt < (get_deg f).toInt :=
    Int64.lt_iff_toInt_lt.mp hlt
  have hgPosInt : 0 < (get_deg g).toInt := by
    simpa [Int64.lt_iff_toInt_lt] using hgPositive
  have hfPosInt : 0 < (get_deg f).toInt := lt_trans hgPosInt hltInt
  have hclampLt : (get_deg g).toNatClampNeg <
      (get_deg f).toNatClampNeg := by
    change (get_deg g).toInt.toNat < (get_deg f).toInt.toNat
    have hgNat := Int.toNat_of_nonneg (le_of_lt hgPosInt)
    have hfNat := Int.toNat_of_nonneg (le_of_lt hfPosInt)
    omega
  rw [StrictDDF.strict_get_deg_toNatClampNeg p g hgCanonical hgNonempty
      hgBound,
    StrictDDF.strict_get_deg_toNatClampNeg p f hfCanonical hfNonempty
      hfBound] at hclampLt
  exact hclampLt

/-- Concrete split law for the generated strict EDF shell.  It derives the
random polynomial from the actual engine transition, the divisor from the
actual candidate execution, and both recursive children from the supplied
exact-div and make-monic executions. -/
noncomputable def strictEDFSplitLaw
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this) :
    Generated.StrictEDF.EDFSplitLaw
      (strictEDFRawOps engine this providers) where
  splitStep := by
    intro f d rngBefore rngAfter r g hRaw gMonic hMonic
      hinvariant hrandom hcandidate hproper hdivRun hgMonicRun hhMonicRun
    change EDFEntryInvariant this f d at hinvariant
    have hfNonempty := hinvariant.nonempty
    have hprime : f[0]!.2.prime = this._p :=
      hinvariant.primeMatches hfNonempty
    have hpWord : 0 < f[0]!.2.prime := by
      simpa [hprime, UInt64.lt_iff_toNat_lt] using
        (Fact.out : Nat.Prime this._p.toNat).pos
    have hrandomStrict :
        Generated.StrictEDF.__upoly_random_raw_ir engine (get_deg f)
          f[0]!.2.prime rngBefore = .ok (r, rngAfter) := by
      simpa [strictEDFRawOps] using hrandom
    rcases __upoly_random_raw_ir_canonical engine (get_deg f)
        f[0]!.2.prime rngBefore hpWord with
      ⟨randomOutput, randomNext, hrandomCanonicalRun,
        hrandomCanonical, _hrandomDegree⟩
    rw [hrandomStrict] at hrandomCanonicalRun
    have hrandomPair := Except.ok.inj hrandomCanonicalRun
    have hrandomOutputEq : r = randomOutput :=
      congrArg Prod.fst hrandomPair
    have hrandomNextEq : rngAfter = randomNext :=
      congrArg Prod.snd hrandomPair
    subst randomOutput
    subst randomNext
    have hrCanonical : SparsePolyZp.Canonical this._p.toNat r := by
      simpa [hprime] using hrandomCanonical
    rcases hcandidate with ⟨hbudget, hcandidateRun⟩
    rcases strictCandidateRun_factor engine this providers f d r hinvariant
        hrCanonical hbudget with
      ⟨candidate, hcandidateStrict, hgCanonical, hgMonic, hgDivides⟩
    rw [hcandidateRun] at hcandidateStrict
    injection hcandidateStrict with hcandidateEq
    subst candidate
    have hgNonempty :=
      Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
        this._p.toNat g hgMonic.ne_zero
    have hfDegreeBound :
        (SparsePolyZp.toPoly this._p.toNat f).natDegree < 2 ^ 62 :=
      StrictDDF.canonical_natDegree_lt_of_terms_lt this._p.toNat f
        hinvariant.canonical hinvariant.monic.ne_zero (2 ^ 62)
        hinvariant.degreeBound
    have hgDegreeLe : (SparsePolyZp.toPoly this._p.toNat g).natDegree ≤
        (SparsePolyZp.toPoly this._p.toNat f).natDegree :=
      Polynomial.natDegree_le_of_dvd hgDivides hinvariant.monic.ne_zero
    have hgDegreeBound :
        (SparsePolyZp.toPoly this._p.toNat g).natDegree < 2 ^ 63 := by
      omega
    have hfDegreeBound63 :
        (SparsePolyZp.toPoly this._p.toNat f).natDegree < 2 ^ 63 := by
      omega
    have hgDegreePositive :
        0 < (SparsePolyZp.toPoly this._p.toNat g).natDegree :=
      (StrictDDF.strict_get_deg_pos_iff_natDegree_pos this._p.toNat g
        hgCanonical hgDegreeBound).mp hproper.1
    have hgDegreeLt : (SparsePolyZp.toPoly this._p.toNat g).natDegree <
        (SparsePolyZp.toPoly this._p.toNat f).natDegree :=
      strict_get_deg_lt_implies_natDegree_lt this._p.toNat g f hgCanonical
        hinvariant.canonical hgNonempty hfNonempty hgDegreeBound
        hfDegreeBound63 hproper.1 hproper.2
    rcases StrictDDF.strictExactDivIR_refines this f g providers.hcfg
        hinvariant.canonical hgCanonical hgNonempty hgMonic hgDivides with
      ⟨quotient, hquotientRun, hquotientCanonical, hquotientSemantic⟩
    have hdivRunStrict : StrictDDF.strictExactDivIR this f g = .ok hRaw := by
      simpa [strictEDFRawOps] using hdivRun
    rw [hdivRunStrict] at hquotientRun
    injection hquotientRun with hquotientEq
    subst quotient
    have hdivByMonicMonic :
        (SparsePolyZp.toPoly this._p.toNat f /ₘ
          SparsePolyZp.toPoly this._p.toNat g).Monic :=
      StrictSquarefreeZp.divByMonic_monic_of_monic_of_dvd
        (SparsePolyZp.toPoly this._p.toNat f)
        (SparsePolyZp.toPoly this._p.toNat g) hinvariant.monic hgMonic
        hgDivides
    have hquotientMonic :
        (SparsePolyZp.toPoly this._p.toNat hRaw).Monic := by
      rw [hquotientSemantic]
      exact hdivByMonicMonic
    have hhNonempty :=
      Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
        this._p.toNat hRaw hquotientMonic.ne_zero
    have hgMonicExact :=
      StrictDDF.strictMakeMonicIR_eq_of_monic this g hgCanonical hgNonempty
        hgMonic
    have hgMonicRunStrict : StrictDDF.strictMakeMonicIR this g =
        .ok gMonic := by
      simpa [strictEDFRawOps] using hgMonicRun
    rw [hgMonicExact] at hgMonicRunStrict
    injection hgMonicRunStrict with hgMonicEq
    subst gMonic
    have hnormalizeRaw : hRaw.normalization = hRaw :=
      StrictSquarefreeZp.sparsePolyZp_normalization_eq_of_canonical
        this._p.toNat hRaw hquotientCanonical
    have hhMonicExact :=
      StrictDDF.strictMakeMonicIR_eq_of_monic this hRaw hquotientCanonical
        hhNonempty hquotientMonic
    have hhMonicRunStrict :
        StrictDDF.strictMakeMonicIR this hRaw.normalization = .ok hMonic := by
      simpa [strictEDFRawOps] using hhMonicRun
    rw [hnormalizeRaw, hhMonicExact] at hhMonicRunStrict
    injection hhMonicRunStrict with hhMonicEq
    subst hMonic
    have hpolySplit := edfPolynomialSplit_of_properDivisor
      (SparsePolyZp.toPoly this._p.toNat f)
      (SparsePolyZp.toPoly this._p.toNat g) d.toNat hinvariant.monic
      hinvariant.squarefree hinvariant.equalDegree hgMonic hgDivides
      hgDegreePositive hgDegreeLt
    have hnormalizeQuotient : normalize
        (SparsePolyZp.toPoly this._p.toNat f /ₘ
          SparsePolyZp.toPoly this._p.toNat g) =
        SparsePolyZp.toPoly this._p.toNat hRaw := by
      rw [hquotientSemantic]
      exact hdivByMonicMonic.normalize_eq_self
    rw [hnormalizeQuotient] at hpolySplit
    have hgDegreeTerms : ∀ term ∈ g.toList, term.1.deg < 2 ^ 62 := by
      intro term hterm
      have := StrictDDF.canonical_term_degree_le_natDegree
        this._p.toNat g hgCanonical term hterm
      omega
    have hhDegreeTerms : ∀ term ∈ hRaw.toList,
        term.1.deg < 2 ^ 62 := by
      intro term hterm
      have htermDegree := StrictDDF.canonical_term_degree_le_natDegree
        this._p.toNat hRaw hquotientCanonical term hterm
      exact lt_of_le_of_lt htermDegree
        (lt_trans hpolySplit.hDegreeLt hfDegreeBound)
    have hgPrimeMatches : 0 < g.size → g[0]!.2.prime = this._p := by
      intro hnonempty
      have hheadMem : g[0] ∈ g.toList := by
        simpa using Array.getElem_mem g 0 hnonempty
      have hprimeNat := (hgCanonical.1 g[0] hheadMem).1
      rw [getElem!_pos g 0 hnonempty]
      exact UInt64.toNat_inj.mp hprimeNat
    have hhPrimeMatches : 0 < hRaw.size →
        hRaw[0]!.2.prime = this._p := by
      intro hnonempty
      have hheadMem : hRaw[0] ∈ hRaw.toList := by
        simpa using Array.getElem_mem hRaw 0 hnonempty
      have hprimeNat := (hquotientCanonical.1 hRaw[0] hheadMem).1
      rw [getElem!_pos hRaw 0 hnonempty]
      exact UInt64.toNat_inj.mp hprimeNat
    let hgInvariant : EDFEntryInvariant this g d :=
      ⟨hgCanonical, hgPrimeMatches, hgDegreeTerms, hpolySplit.gMonic,
        hpolySplit.gDegreePositive, hinvariant.dPositive,
        hpolySplit.gSquarefree, hpolySplit.gEqualDegree⟩
    let hhInvariant : EDFEntryInvariant this hRaw d :=
      ⟨hquotientCanonical, hhPrimeMatches, hhDegreeTerms,
        hpolySplit.hMonic, hpolySplit.hDegreePositive,
        hinvariant.dPositive, hpolySplit.hSquarefree,
        hpolySplit.hEqualDegree⟩
    refine ⟨hgInvariant, hhInvariant, ?_, ?_⟩
    · rw [edfMeasure_eq_natDegree_succ this._p.toNat g hgCanonical
          hgNonempty,
        edfMeasure_eq_natDegree_succ this._p.toNat f hinvariant.canonical
          hfNonempty]
      exact Nat.add_lt_add_right hpolySplit.gDegreeLt 1
    · rw [edfMeasure_eq_natDegree_succ this._p.toNat hRaw
          hquotientCanonical hhNonempty,
        edfMeasure_eq_natDegree_succ this._p.toNat f hinvariant.canonical
          hfNonempty]
      exact Nat.add_lt_add_right hpolySplit.hDegreeLt 1

/-- The well-founded generated EDF state machine terminates along its exact
retry traces and appends a genuine `EDFCorrect` factorization of the current
input to the existing raw accumulator. -/
theorem strictEDFState_refines
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (strictEDFRawOps engine this providers))
    (state : Generated.StrictEDF.EDFState
      (strictEDFRawOps engine this providers))
    (hresultCanonical : ∀ factor ∈ state.result.toList,
      SparsePolyZp.Canonical this._p.toNat factor) :
    ∃ output rng factors,
      Generated.StrictEDF.__edf_Zp_raw_ir_state
          (strictEDFRawOps engine this providers)
          (strictEDFSplitLaw engine this providers) termination state =
        .ok (output, rng) ∧
      edfResultToL2 this._p.toNat output =
        edfResultToL2 this._p.toNat state.result ++ factors ∧
      EDFCorrect (SparsePolyZp.toPoly this._p.toNat state.f)
        state.d.toNat factors ∧
      ∀ factor ∈ output.toList,
        SparsePolyZp.Canonical this._p.toNat factor := by
  generalize hmeasure : Generated.StrictEDF.edfMeasure state.f = measure
  induction measure using Nat.strong_induction_on generalizing state with
  | h measure ih =>
      have hvalid : EDFEntryInvariant this state.f state.d := state.valid
      have hfNonempty := hvalid.nonempty
      have hfDegreeBound := hvalid.natDegree_lt
      rw [Generated.StrictEDF.__edf_Zp_raw_ir_state.eq_1]
      split
      next hbase =>
        have hdegreeWord := StrictDDF.strict_get_deg_toUInt64_toNat
          this._p.toNat state.f hvalid.canonical hfNonempty
          hfDegreeBound
        have hwordEq : (get_deg state.f).toUInt64 = state.d := by
          simpa using hbase
        have hdegreeEq :
            (SparsePolyZp.toPoly this._p.toNat state.f).natDegree =
              state.d.toNat := by
          have := congrArg UInt64.toNat hwordEq
          omega
        have hmake := StrictDDF.strictMakeMonicIR_eq_of_monic this state.f
          hvalid.canonical hfNonempty hvalid.monic
        have hmakeOps :
            (strictEDFRawOps engine this providers).makeMonic state.f =
              .ok state.f := by
          simpa [strictEDFRawOps] using hmake
        rw [certifyRawExec_ok_eq _ _ hmakeOps]
        refine ⟨state.result.push state.f, state.rng,
          [SparsePolyZp.toPoly this._p.toNat state.f], rfl, ?_, ?_, ?_⟩
        · simp
        · exact base_factor_correct
            (SparsePolyZp.toPoly this._p.toNat state.f) state.d.toNat
            hvalid.monic hvalid.dPositive hdegreeEq hvalid.equalDegree
        · intro factor hfactor
          rw [Array.toList_push, List.mem_append] at hfactor
          rcases hfactor with hfactor | hfactor
          · exact hresultCanonical factor hfactor
          · have heq : factor = state.f := by simpa using hfactor
            subst factor
            exact hvalid.canonical
      next hnotBase =>
        split
        next hnonpositive =>
          have hpositive : get_deg state.f > 0 :=
            (StrictDDF.strict_get_deg_pos_iff_natDegree_pos
              this._p.toNat state.f hvalid.canonical
              hfDegreeBound).mpr hvalid.degreePositive
          have hpositiveInt : (0 : Int64).toInt < (get_deg state.f).toInt :=
            Int64.lt_iff_toInt_lt.mp hpositive
          have hnonpositiveInt : (get_deg state.f).toInt ≤ (0 : Int64).toInt :=
            Int64.le_iff_toInt_le.mp hnonpositive
          omega
        next hpositiveBranch =>
          rcases retryLoop_terminates
              (strictEDFRawOps engine this providers) state.f state.d
              state.rng
              (termination.retryTrace state.f state.d state.rng hvalid) with
            ⟨splitState, hretryRun⟩
          rw [hretryRun]
          rcases strictSplitCalls_terminate engine this providers state.f
              state.d splitState hvalid with
            ⟨hRaw, gMonic, hMonic, hdivRun, hgMonicRun, hhMonicRun⟩
          simp only
          rw [certifyRawExec_ok_eq _ _ hdivRun]
          simp only
          rw [certifyRawExec_ok_eq _ _ hgMonicRun]
          simp only
          rw [certifyRawExec_ok_eq _ _ hhMonicRun]
          simp only
          have hstep := (strictEDFSplitLaw engine this providers).splitStep
            state.f state.d splitState.rngBefore splitState.rng
            splitState.randomPoly splitState.factor hRaw gMonic hMonic
            hvalid splitState.randomRun splitState.candidateRun
            splitState.proper hdivRun hgMonicRun hhMonicRun
          let leftState : Generated.StrictEDF.EDFState
              (strictEDFRawOps engine this providers) :=
            ⟨state.result, gMonic, state.d, splitState.rng, hstep.1⟩
          have hleftMeasure :
              Generated.StrictEDF.edfMeasure leftState.f < measure := by
            dsimp [leftState]
            rw [← hmeasure]
            exact hstep.2.2.1
          rcases ih _ hleftMeasure leftState (by simpa [leftState] using hresultCanonical) rfl with
            ⟨leftOutput, leftRng, leftFactors, hleftRun,
              hleftDecode, hleftCorrect, hleftCanonical⟩
          rw [hleftRun]
          let rightState : Generated.StrictEDF.EDFState
              (strictEDFRawOps engine this providers) :=
            ⟨leftOutput, hMonic, state.d, leftRng, hstep.2.1⟩
          have hrightMeasure :
              Generated.StrictEDF.edfMeasure rightState.f < measure := by
            dsimp [rightState]
            rw [← hmeasure]
            exact hstep.2.2.2
          rcases ih _ hrightMeasure rightState (by
              simpa [rightState] using hleftCanonical) rfl with
            ⟨output, rng, rightFactors, hrightRun,
              hrightDecode, hrightCorrect, hrightCanonical⟩
          have hproduct := strictSplitCalls_product engine this providers
            state.f state.d splitState.rngBefore splitState.rng
            splitState.randomPoly splitState.factor hRaw gMonic hMonic
            hvalid splitState.randomRun splitState.candidateRun.choose
            splitState.candidateRun.choose_spec hdivRun hgMonicRun hhMonicRun
          have hassociated : Associated
              (SparsePolyZp.toPoly this._p.toNat state.f)
              (SparsePolyZp.toPoly this._p.toNat gMonic *
                SparsePolyZp.toPoly this._p.toNat hMonic) := by
            rw [hproduct]
          have hcorrect := edf_combine
            (SparsePolyZp.toPoly this._p.toNat state.f)
            (SparsePolyZp.toPoly this._p.toNat gMonic)
            (SparsePolyZp.toPoly this._p.toNat hMonic) state.d.toNat
            hassociated leftFactors rightFactors hleftCorrect hrightCorrect
          refine ⟨output, rng, leftFactors ++ rightFactors, hrightRun, ?_,
            hcorrect, hrightCanonical⟩
          rw [hrightDecode, hleftDecode]
          simp [rightState, leftState, List.append_assoc]

/-- Exact entry-point refinement for the cpp2lean-generated C++ `__edf_Zp`
semantics.  The concrete accumulator and RNG transition are preserved, while
the newly appended concrete factors decode to an `EDFCorrect` L2 result. -/
theorem strictEDFEntryIR_refines_edf
    {State : Type} (engine : Generated.StrictEDF.RandomEngine State)
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : StrictDDF.DDFRawProviders this)
    (termination : Generated.StrictEDF.EDFTermination
      (strictEDFRawOps engine this providers))
    (result : Array SparsePolyZp) (f : SparsePolyZp) (d : UInt64)
    (rng : State) (hinvariant : EDFEntryInvariant this f d)
    (hresultCanonical : ∀ factor ∈ result.toList,
      SparsePolyZp.Canonical this._p.toNat factor) :
    ∃ output rng' factors,
      Generated.StrictEDF.__edf_Zp_raw_ir
          (strictEDFRawOps engine this providers)
          (strictEDFSplitLaw engine this providers) termination
          result f d rng hinvariant = .ok (output, rng') ∧
      edfResultToL2 this._p.toNat output =
        edfResultToL2 this._p.toNat result ++ factors ∧
      EDFCorrect (SparsePolyZp.toPoly this._p.toNat f) d.toNat factors ∧
      ∀ factor ∈ output.toList,
        SparsePolyZp.Canonical this._p.toNat factor := by
  simpa [Generated.StrictEDF.__edf_Zp_raw_ir] using
    strictEDFState_refines engine this providers termination
      (⟨result, f, d, rng, hinvariant⟩ : Generated.StrictEDF.EDFState
        (strictEDFRawOps engine this providers)) hresultCanonical

end StrictEDF

end Refinement
