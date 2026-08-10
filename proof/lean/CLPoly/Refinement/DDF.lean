/-
  Strict DDF refinement boundary.

  The previous contents proved modular exponentiation by importing the legacy
  `Generated.__upoly_mod_ir`, whose body dispatches directly to the hand-written
  `SparsePolyZp.divmod`.  That chain did not establish execution of the C++
  polynomial-division loops, so it is no longer exported as L1→L2 refinement.

  DDF will be reopened only after the RawHeap `_poly_divrem` theorem proves the
  quotient/remainder identity and remainder-degree bound from the four
  cpp2lean-generated loops.  Powmod and the DDF loop must consume that theorem
  directly.
-/
import CLPoly.Algorithm.DDF
import CLPoly.Generated.StrictDDF
import CLPoly.Impl.StrictDivremRefinement
import CLPoly.Impl.StrictPolynomialGCDRefinement
import CLPoly.Impl.StrictMulDispatchRefinement
import CLPoly.Refinement.DDFSubtractX
import CLPoly.Refinement.SquarefreeZp

set_option autoImplicit false

namespace Refinement

namespace StrictDDF

set_option maxHeartbeats 0

open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictPolynomialGCDRefinement
open CLPoly.Impl.StrictEuclidRefinement
open Generated.StrictDivrem
open CLPoly.Math

/-- Physical buffers for one generated `__upoly_mod` call.  This package
contains allocation/layout facts only; in particular, it cannot choose a
remainder or any L2 polynomial result. -/
structure RawModWorkspace (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) where
  heap : RawHeap
  dividendPtr : RawPtr UInt64
  divisorPtr : RawPtr UInt64
  quotientPtr : RawPtr UInt64
  remainderPtr : RawPtr UInt64
  scratchPtr : RawPtr Word3
  dividendValid : heap.ValidU64Slice dividendPtr
    (sparseDenseLength dividend)
  divisorValid : heap.ValidU64Slice divisorPtr
    (sparseDenseLength divisor)
  quotientValid : heap.ValidU64Slice quotientPtr
    (sparseDenseLength dividend - (sparseDenseLength divisor - 1))
  remainderValid : heap.ValidU64Slice remainderPtr
    (Nat.min (sparseDenseLength dividend) (sparseDenseLength divisor - 1))
  scratchValid : heap.ValidWord3Slice scratchPtr
    (sparseDenseLength dividend)
  inputsDisjoint : CLPoly.Impl.StrictMulRefinement.U64SlicesDisjoint
    divisorPtr (sparseDenseLength divisor) dividendPtr
      (sparseDenseLength dividend)
  remainderDividendDisjoint : remainderPtr.region ≠ dividendPtr.region
  scratchDividendDisjoint : scratchPtr.region ≠ dividendPtr.region
  scratchDivisorDisjoint : scratchPtr.region ≠ divisorPtr.region
  quotientDivisorDisjoint : quotientPtr.region ≠ divisorPtr.region
  quotientScratchDisjoint : quotientPtr.region ≠ scratchPtr.region
  remainderScratchDisjoint : remainderPtr.region ≠ scratchPtr.region
  remainderQuotientDisjoint : remainderPtr.region ≠ quotientPtr.region
  remainderDivisorDisjoint : remainderPtr.region ≠ divisorPtr.region
  quotientCapacity : sparseDenseLength dividend -
      (sparseDenseLength divisor - 1) <
      CLPoly.Impl.StrictWordArithmetic.limbBase

/-- Concrete raw storage and readiness obligations for one nonempty public
C++ `polynomial_GCD` call.  This package contains no polynomial result and no
semantic oracle: every output is obtained by executing the generated dense
GCD/HGCD/Euclid path. -/
structure RawGCDWorkspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (left right : SparsePolyZp) where
  M : HgcdMat
  hM : M.Valid
  resultPtr : RawPtr UInt64
  leftPtr : RawPtr UInt64
  rightPtr : RawPtr UInt64
  aBuf : RawPtr UInt64
  bBuf : RawPtr UInt64
  J : RawPtr UInt64
  Q : RawPtr UInt64
  R : RawPtr UInt64
  W3 : RawPtr Word3
  W : RawPtr UInt64
  scratch : RawPtr UInt64
  euclidQ : RawPtr UInt64
  euclidR : RawPtr UInt64
  euclidW3 : RawPtr Word3
  heap : RawHeap
  hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this
  hp : 1 < this._p.toNat
  loopDecrease : Generated.StrictGCDHGCD.HgcdGcdLoopLengthDecreases
    this M hM W scratch
  leftValid : heap.ValidU64Slice leftPtr (sparseDenseLength left)
  rightValid : heap.ValidU64Slice rightPtr (sparseDenseLength right)
  leftBound : sparseDenseLength left ≤ 2 ^ 63
  rightBound : sparseDenseLength right ≤ 2 ^ 63
  inputsDisjoint : CLPoly.Impl.StrictMulRefinement.U64SlicesDisjoint rightPtr
    (sparseDenseLength right) leftPtr (sparseDenseLength left)
  readyAB : ∀ finalHeap,
    RawDensePolyRep this finalHeap leftPtr (sparseDenseLength left)
        (SparsePolyZp.toPoly this._p.toNat left) →
    RawDensePolyRep this finalHeap rightPtr (sparseDenseLength right)
        (SparsePolyZp.toPoly this._p.toNat right) →
    DenseGcdOrderedReady this hcfg hp M hM resultPtr leftPtr rightPtr aBuf
      bBuf J Q R W3 W scratch euclidQ euclidR euclidW3
      (sparseDenseLength left) (sparseDenseLength right) finalHeap
      (SparsePolyZp.toPoly this._p.toNat left)
      (SparsePolyZp.toPoly this._p.toNat right)
  readyBA : ∀ finalHeap,
    RawDensePolyRep this finalHeap leftPtr (sparseDenseLength left)
        (SparsePolyZp.toPoly this._p.toNat left) →
    RawDensePolyRep this finalHeap rightPtr (sparseDenseLength right)
        (SparsePolyZp.toPoly this._p.toNat right) →
    DenseGcdOrderedReady this hcfg hp M hM resultPtr rightPtr leftPtr aBuf
      bBuf J Q R W3 W scratch euclidQ euclidR euclidW3
      (sparseDenseLength right) (sparseDenseLength left) finalHeap
      (SparsePolyZp.toPoly this._p.toNat right)
      (SparsePolyZp.toPoly this._p.toNat left)

/-- Actual raw execution of one DDF modular reduction: construct both dense
inputs, execute generated `_poly_divrem`, then run the concrete reverse
dense-to-sparse scan over the returned remainder slice. -/
def strictModIR (this : DenseUPolyZp) (dividend divisor : SparsePolyZp)
    (workspace : RawModWorkspace this dividend divisor) :
    RawExec SparsePolyZp :=
  match sparse_upoly_zp_dense_constructor_raw_ir workspace.dividendPtr
      dividend workspace.heap with
  | .error fault => .error fault
  | .ok dividendHeap =>
    match sparse_upoly_zp_dense_constructor_raw_ir workspace.divisorPtr
        divisor dividendHeap with
    | .error fault => .error fault
    | .ok inputHeap =>
      match dense_upoly_zp__poly_divrem_ir this workspace.quotientPtr
          workspace.remainderPtr workspace.dividendPtr
          (sparseDenseLength dividend) workspace.divisorPtr
          (sparseDenseLength divisor) workspace.scratchPtr inputHeap with
      | .error fault => .error fault
      | .ok (outputHeap, _, remainderLength) =>
        dense_upoly_zp_to_upoly_raw_ir this outputHeap workspace.remainderPtr
          remainderLength

/-- Total semantic refinement of the raw modular-reduction execution. -/
theorem strictModIR_refines_modByMonic
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (dividend divisor : SparsePolyZp)
    (workspace : RawModWorkspace this dividend divisor)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hdivisorCanonical : SparsePolyZp.Canonical this._p.toNat divisor)
    (hdivisorNonempty : 0 < divisor.size)
    (hdivisorMonic : (SparsePolyZp.toPoly this._p.toNat divisor).Monic) :
    ∃ remainder,
      strictModIR this dividend divisor workspace = .ok remainder ∧
      SparsePolyZp.Canonical this._p.toNat remainder ∧
      SparsePolyZp.toPoly this._p.toNat remainder =
        SparsePolyZp.toPoly this._p.toNat dividend %ₘ
          SparsePolyZp.toPoly this._p.toNat divisor := by
  have hpNonzero : this._p ≠ 0 := by
    intro hzero
    apply (Fact.out : Nat.Prime this._p.toNat).ne_zero
    simpa using congrArg UInt64.toNat hzero
  rcases sparse_upoly_zp_dense_constructor_raw_ir_result this
      workspace.dividendPtr dividend workspace.heap hpNonzero
      workspace.dividendValid hdividendCanonical with
    ⟨dividendHeap, hdividendRun, hdividendLayout, hdividendRep⟩
  have hdivisorValid' : dividendHeap.ValidU64Slice workspace.divisorPtr
      (sparseDenseLength divisor) :=
    (hdividendLayout workspace.divisorPtr _).mp workspace.divisorValid
  rcases sparse_upoly_zp_dense_constructor_raw_ir_result this
      workspace.divisorPtr divisor dividendHeap hpNonzero hdivisorValid'
      hdivisorCanonical with
    ⟨inputHeap, hdivisorRun, hdivisorLayout, hdivisorRep⟩
  have hdividendDense :=
    sparse_upoly_zp_dense_constructor_raw_ir_preserves_rawDense this
      workspace.divisorPtr workspace.dividendPtr divisor dividendHeap inputHeap
      (sparseDenseLength dividend)
      (SparsePolyZp.toPoly this._p.toNat dividend) hpNonzero hdivisorValid'
      hdivisorCanonical workspace.inputsDisjoint hdividendRep.dense
      hdivisorRun
  have hlayout : RawHeap.SameLayout workspace.heap inputHeap :=
    fun ptr length => (hdividendLayout ptr length).trans
      (hdivisorLayout ptr length)
  have hquotientValid :=
    (hlayout workspace.quotientPtr _).mp workspace.quotientValid
  have hremainderValid :=
    (hlayout workspace.remainderPtr _).mp workspace.remainderValid
  have hscratchValid :=
    (hlayout (RawPtr.reinterpret workspace.scratchPtr) _).mp
      workspace.scratchValid
  have hdivisorLength : 0 < sparseDenseLength divisor := by
    simp [sparseDenseLength, hdivisorNonempty]
  rcases polyDivrem_refines this workspace.quotientPtr workspace.remainderPtr
      workspace.dividendPtr workspace.divisorPtr
      (sparseDenseLength dividend) (sparseDenseLength divisor)
      workspace.scratchPtr inputHeap
      (SparsePolyZp.toPoly this._p.toNat dividend)
      (SparsePolyZp.toPoly this._p.toNat divisor) hdivisorLength
      hdividendDense.2.1 hdivisorRep.dense.1 hquotientValid hremainderValid
      hscratchValid hdividendDense.2.2.1 hdivisorRep.dense.2.1
      hdividendDense.2.2.2.1 hdivisorRep.dense.2.2.1
      hdividendDense.2.2.2.2 hdivisorRep.dense.2.2.2
      workspace.quotientCapacity workspace.remainderDividendDisjoint
      workspace.scratchDividendDisjoint workspace.scratchDivisorDisjoint
      workspace.quotientDivisorDisjoint workspace.quotientScratchDisjoint
      workspace.remainderScratchDisjoint workspace.remainderQuotientDisjoint
      workspace.remainderDivisorDisjoint hcfg
      (Fact.out : Nat.Prime this._p.toNat) with
    ⟨outputHeap, quotientLength, remainderLength, quotientPoly, remainderPoly,
      hdivRun, _, _, _, hremainderRep, hremainderCanonicalRaw,
      hremainderNorm, houtputLayout, _, halgebra, hremainderDegree, _,
      hremainderLengthLe,
      hremainderLengthLt⟩
  have hremainderValidInput : inputHeap.ValidU64Slice workspace.remainderPtr
      remainderLength := inputHeap.validU64Slice_mono workspace.remainderPtr _ _
    hremainderValid hremainderLengthLe
  have hremainderValidOut : outputHeap.ValidU64Slice workspace.remainderPtr
      remainderLength :=
    (houtputLayout workspace.remainderPtr remainderLength).mp
      hremainderValidInput
  have hremainderRaw : RawDensePolyRep this outputHeap
      workspace.remainderPtr remainderLength remainderPoly :=
    ⟨hremainderValidOut, hremainderCanonicalRaw, hremainderRep,
      hremainderNorm⟩
  rcases dense_upoly_zp_to_upoly_raw_ir_refines this outputHeap
      workspace.remainderPtr remainderLength remainderPoly hremainderRaw with
    ⟨remainder, hremainderRun, hremainderSemantic⟩
  rcases dense_upoly_zp_to_upoly_raw_ir_canonical this outputHeap
      workspace.remainderPtr remainderLength hremainderValidOut
      hremainderCanonicalRaw with
    ⟨canonicalRemainder, hcanonicalRun, hcanonical⟩
  have hsame : canonicalRemainder = remainder :=
    Except.ok.inj (hcanonicalRun.symm.trans hremainderRun)
  subst canonicalRemainder
  have hdegree : Polynomial.degree remainderPoly <
      Polynomial.degree (SparsePolyZp.toPoly this._p.toNat divisor) := by
    by_cases hzero : remainderPoly = 0
    · rw [hzero, Polynomial.degree_zero, bot_lt_iff_ne_bot,
        Polynomial.degree_ne_bot]
      exact hdivisorMonic.ne_zero
    · have hlt := hremainderDegree.resolve_left hzero
      rw [Polynomial.degree_eq_natDegree hzero,
        Polynomial.degree_eq_natDegree hdivisorMonic.ne_zero]
      exact_mod_cast hlt
  have hmod : SparsePolyZp.toPoly this._p.toNat dividend %ₘ
      SparsePolyZp.toPoly this._p.toNat divisor = remainderPoly := by
    exact (Polynomial.div_modByMonic_unique quotientPoly remainderPoly
      hdivisorMonic ⟨by rw [halgebra]; ring, hdegree⟩).2
  refine ⟨remainder, ?_, hcanonical, ?_⟩
  · simp [strictModIR, hdividendRun, hdivisorRun, hdivRun, hremainderRun]
  · rw [hremainderSemantic, ← hmod]

/-- Execute the complete nonempty public C++ GCD boundary and project its
concrete sparse result. -/
def strictGCDIR (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (left right : SparsePolyZp)
    (workspace : RawGCDWorkspace this left right) : RawExec SparsePolyZp :=
  match polynomialGCDPublicNonemptyRawIR this workspace.M workspace.hM
      workspace.resultPtr workspace.leftPtr workspace.rightPtr
      workspace.aBuf workspace.bBuf workspace.J workspace.Q workspace.R
      workspace.W3 workspace.W workspace.scratch workspace.euclidQ
      workspace.euclidR workspace.euclidW3 left right workspace.loopDecrease
      workspace.heap with
  | .error fault => .error fault
  | .ok final =>
    match final.output with
    | none => .error .assertionFailure
    | some result => .ok result

/-- End-to-end raw execution and semantics of the concrete public C++ GCD.
No `SparsePolyZp.gcd` implementation occurs in this execution path. -/
theorem strictGCDIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (left right : SparsePolyZp)
    (workspace : RawGCDWorkspace this left right)
    (hleftCanonical : SparsePolyZp.Canonical this._p.toNat left)
    (hrightCanonical : SparsePolyZp.Canonical this._p.toNat right)
    (hleftNonzero : SparsePolyZp.toPoly this._p.toNat left ≠ 0)
    (hrightNonzero : SparsePolyZp.toPoly this._p.toNat right ≠ 0) :
    ∃ result,
      strictGCDIR this left right workspace = .ok result ∧
      SparsePolyZp.Canonical this._p.toNat result ∧
      SparsePolyZp.toPoly this._p.toNat result = normalize
        (EuclideanDomain.gcd (SparsePolyZp.toPoly this._p.toNat left)
          (SparsePolyZp.toPoly this._p.toNat right)) := by
  rcases polynomialGCDPublicNonemptyRawIR_refines this workspace.M
      workspace.hM workspace.resultPtr workspace.leftPtr workspace.rightPtr
      workspace.aBuf workspace.bBuf workspace.J workspace.Q workspace.R
      workspace.W3 workspace.W workspace.scratch workspace.euclidQ
      workspace.euclidR workspace.euclidW3 left right workspace.loopDecrease
      workspace.heap workspace.hcfg workspace.hp workspace.leftValid
      workspace.rightValid hleftCanonical hrightCanonical hleftNonzero
      hrightNonzero workspace.leftBound workspace.rightBound
      workspace.inputsDisjoint workspace.readyAB workspace.readyBA with
    ⟨final, result, hrun, houtput, hcanonical, hsemantic⟩
  refine ⟨result, ?_, hcanonical, hsemantic⟩
  simp [strictGCDIR, hrun, houtput]

/-- Exact source `pair_vec_div` call used by the DDF factor-removal branch.
The implementation is the already verified generated VHC/single-term C++
control flow, not polynomial division from the L2 model. -/
def strictExactDivIR (this : DenseUPolyZp)
    (dividend divisor : SparsePolyZp) : RawExec SparsePolyZp :=
  Refinement.StrictSquarefreeZp.pairVecDivIR this dividend divisor

/-- Project the mutated polynomial from the checked C++ make-monic entry. -/
def strictMakeMonicIR (this : DenseUPolyZp) (f : SparsePolyZp) :
    RawExec SparsePolyZp :=
  match Refinement.StrictSquarefreeZp.upolyMakeMonicIR this f with
  | .error fault => .error fault
  | .ok (_, monic) => .ok monic

/-- A canonical monic input takes the actual source early-return branch of
`__upoly_make_monic` and is returned unchanged. -/
theorem strictMakeMonicIR_eq_of_monic
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (f : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hnonempty : 0 < f.size)
    (hmonic : (SparsePolyZp.toPoly this._p.toNat f).Monic) :
    strictMakeMonicIR this f = .ok f := by
  have hrun := Refinement.StrictSquarefreeZp.upolyMakeMonicIR_eq_of_monic
    this f hcanonical hnonempty hmonic
  simp [strictMakeMonicIR, hrun]

/-- The exact two public C++ `polynomial_GCD` branches needed by DDF.  A zero
left operand takes the source wrapper's copy-and-monic path; otherwise the
verified dense Euclidean implementation is executed. -/
def strictDDFGCDIR (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (left right : SparsePolyZp)
    (workspace : RawGCDWorkspace this left right) : RawExec SparsePolyZp :=
  if left.isEmpty then strictMakeMonicIR this right
  else strictGCDIR this left right workspace

/-- Semantic refinement of the public GCD call in `__ddf_Zp`, including the
source-level zero-left branch that the nonempty dense GCD entry cannot take. -/
theorem strictDDFGCDIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (left right : SparsePolyZp)
    (workspace : RawGCDWorkspace this left right)
    (hleftCanonical : SparsePolyZp.Canonical this._p.toNat left)
    (hrightCanonical : SparsePolyZp.Canonical this._p.toNat right)
    (hrightNonempty : 0 < right.size)
    (hrightMonic : (SparsePolyZp.toPoly this._p.toNat right).Monic) :
    ∃ result,
      strictDDFGCDIR this left right workspace = .ok result ∧
      SparsePolyZp.Canonical this._p.toNat result ∧
      SparsePolyZp.toPoly this._p.toNat result = normalize
        (EuclideanDomain.gcd (SparsePolyZp.toPoly this._p.toNat left)
          (SparsePolyZp.toPoly this._p.toNat right)) := by
  by_cases hleftEmpty : left.isEmpty
  · have hleftSize : left.size = 0 := by
      simpa [Array.isEmpty] using hleftEmpty
    have hleftEq : left = #[] := Array.size_eq_zero_iff.mp hleftSize
    have hrun := Refinement.StrictSquarefreeZp.upolyMakeMonicIR_eq_of_monic
      this right hrightCanonical hrightNonempty hrightMonic
    refine ⟨right, ?_, hrightCanonical, ?_⟩
    · simp [strictDDFGCDIR, hleftEmpty, strictMakeMonicIR, hrun]
    · subst left
      simpa using hrightMonic.normalize_eq_self.symm
  · have hleftNonempty : 0 < left.size := by
      apply Nat.pos_of_ne_zero
      simpa [Array.isEmpty] using hleftEmpty
    have hleftNonzero : SparsePolyZp.toPoly this._p.toNat left ≠ 0 := by
      intro hzero
      have hdegree :=
        Refinement.StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
          this._p.toNat left hleftCanonical hleftNonempty
      rw [hzero] at hdegree
      simp at hdegree
    rcases strictGCDIR_refines this left right workspace hleftCanonical
        hrightCanonical hleftNonzero hrightMonic.ne_zero with
      ⟨result, hrun, hcanonical, hsemantic⟩
    exact ⟨result, by simp [strictDDFGCDIR, hleftEmpty, hrun],
      hcanonical, hsemantic⟩

theorem strictExactDivIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (dividend divisor : SparsePolyZp)
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (hdividendCanonical : SparsePolyZp.Canonical this._p.toNat dividend)
    (hdivisorCanonical : SparsePolyZp.Canonical this._p.toNat divisor)
    (hdivisorNonempty : 0 < divisor.size)
    (hdivisorMonic :
      (SparsePolyZp.toPoly this._p.toNat divisor).Monic)
    (hdivides : SparsePolyZp.toPoly this._p.toNat divisor ∣
      SparsePolyZp.toPoly this._p.toNat dividend) :
    ∃ quotient,
      strictExactDivIR this dividend divisor = .ok quotient ∧
      SparsePolyZp.Canonical this._p.toNat quotient ∧
      SparsePolyZp.toPoly this._p.toNat quotient =
        SparsePolyZp.toPoly this._p.toNat dividend /ₘ
          SparsePolyZp.toPoly this._p.toNat divisor := by
  exact Refinement.StrictSquarefreeZp.pairVecDivIR_refines_divByMonic
    this dividend divisor hcfg hdividendCanonical hdivisorCanonical
    hdivisorNonempty hdivisorMonic hdivides

/-- Physical buffers for an ordered sparse multiplication.  As with raw
modular reduction, no semantic product is stored in the workspace. -/
structure RawMulWorkspace (this : DenseUPolyZp)
    (left right : SparsePolyZp) where
  heap : RawHeap
  leftPtr : RawPtr UInt64
  rightPtr : RawPtr UInt64
  outputPtr : RawPtr UInt64
  scratchPtr : RawPtr UInt64
  leftValid : heap.ValidU64Slice leftPtr (sparseDenseLength left)
  rightValid : heap.ValidU64Slice rightPtr (sparseDenseLength right)
  outputValid : heap.ValidU64Slice outputPtr
    (2 * sparseDenseLength left - 1)
  scratchValid : heap.ValidU64Slice scratchPtr
    (8 * sparseDenseLength left)
  inputsDisjoint : CLPoly.Impl.StrictMulRefinement.U64SlicesDisjoint
    rightPtr (sparseDenseLength right) leftPtr (sparseDenseLength left)
  outputLeftDisjoint : CLPoly.Impl.StrictMulRefinement.U64SlicesDisjoint
    outputPtr (2 * sparseDenseLength left - 1) leftPtr
      (sparseDenseLength left)
  outputRightDisjoint : CLPoly.Impl.StrictMulRefinement.U64SlicesDisjoint
    outputPtr (2 * sparseDenseLength left - 1) rightPtr
      (sparseDenseLength right)
  outputScratchDisjoint : CLPoly.Impl.StrictMulRefinement.U64SlicesDisjoint
    outputPtr (2 * sparseDenseLength left - 1) scratchPtr
      (8 * sparseDenseLength left)
  scratchLeftDisjoint : CLPoly.Impl.StrictMulRefinement.U64SlicesDisjoint
    scratchPtr (8 * sparseDenseLength left) leftPtr
      (sparseDenseLength left)
  scratchRightDisjoint : CLPoly.Impl.StrictMulRefinement.U64SlicesDisjoint
    scratchPtr (8 * sparseDenseLength left) rightPtr
      (sparseDenseLength right)
  leftLengthWord : sparseDenseLength left <
    CLPoly.Impl.StrictWordArithmetic.limbBase

/-- Actual sparse→dense, generated `_mul`, dense→sparse execution.  Callers
order operands by dense length before choosing this entry. -/
def strictMulOrderedIR (this : DenseUPolyZp)
    (left right : SparsePolyZp) (workspace : RawMulWorkspace this left right) :
    RawExec SparsePolyZp :=
  match sparse_upoly_zp_dense_constructor_raw_ir workspace.leftPtr left
      workspace.heap with
  | .error fault => .error fault
  | .ok leftHeap =>
    match sparse_upoly_zp_dense_constructor_raw_ir workspace.rightPtr right
        leftHeap with
    | .error fault => .error fault
    | .ok inputHeap =>
      match Generated.StrictMul.dense_upoly_zp__mul_ir this
          workspace.outputPtr workspace.leftPtr (sparseDenseLength left)
          workspace.rightPtr (sparseDenseLength right) workspace.scratchPtr
          inputHeap with
      | .error fault => .error fault
      | .ok outputHeap =>
        dense_upoly_zp_to_upoly_raw_ir this outputHeap workspace.outputPtr
          (sparseDenseLength left + sparseDenseLength right - 1)

/-- Total semantic refinement of the ordered raw multiplication execution. -/
theorem strictMulOrderedIR_refines_mul
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (left right : SparsePolyZp)
    (workspace : RawMulWorkspace this left right)
    (hleftCanonical : SparsePolyZp.Canonical this._p.toNat left)
    (hrightCanonical : SparsePolyZp.Canonical this._p.toNat right)
    (hleftNonempty : 0 < left.size)
    (hrightNonempty : 0 < right.size)
    (horder : sparseDenseLength right ≤ sparseDenseLength left) :
    ∃ product,
      strictMulOrderedIR this left right workspace = .ok product ∧
      SparsePolyZp.Canonical this._p.toNat product ∧
      SparsePolyZp.toPoly this._p.toNat product =
        SparsePolyZp.toPoly this._p.toNat left *
          SparsePolyZp.toPoly this._p.toNat right := by
  have hpNonzero : this._p ≠ 0 := by
    intro hzero
    apply (Fact.out : Nat.Prime this._p.toNat).ne_zero
    simpa using congrArg UInt64.toNat hzero
  rcases sparse_upoly_zp_dense_constructor_raw_ir_result this
      workspace.leftPtr left workspace.heap hpNonzero workspace.leftValid
      hleftCanonical with ⟨leftHeap, hleftRun, hleftLayout, hleftRep⟩
  have hrightValid' : leftHeap.ValidU64Slice workspace.rightPtr
      (sparseDenseLength right) :=
    (hleftLayout workspace.rightPtr _).mp workspace.rightValid
  rcases sparse_upoly_zp_dense_constructor_raw_ir_result this
      workspace.rightPtr right leftHeap hpNonzero hrightValid'
      hrightCanonical with
    ⟨inputHeap, hrightRun, hrightLayout, hrightRep⟩
  have hleftDense :=
    sparse_upoly_zp_dense_constructor_raw_ir_preserves_rawDense this
      workspace.rightPtr workspace.leftPtr right leftHeap inputHeap
      (sparseDenseLength left) (SparsePolyZp.toPoly this._p.toNat left)
      hpNonzero hrightValid' hrightCanonical workspace.inputsDisjoint
      hleftRep.dense hrightRun
  have hlayout : RawHeap.SameLayout workspace.heap inputHeap :=
    fun ptr length => (hleftLayout ptr length).trans
      (hrightLayout ptr length)
  have houtputValid :=
    (hlayout workspace.outputPtr _).mp workspace.outputValid
  have hscratchValid :=
    (hlayout workspace.scratchPtr _).mp workspace.scratchValid
  have hleftLengthPositive : 0 < sparseDenseLength left := by
    simp [sparseDenseLength, hleftNonempty]
  have hrightLengthPositive : 0 < sparseDenseLength right := by
    simp [sparseDenseLength, hrightNonempty]
  rcases CLPoly.Impl.StrictMulRefinement.mul_refines_rawDense this
      workspace.outputPtr workspace.leftPtr (sparseDenseLength left)
      workspace.rightPtr (sparseDenseLength right) workspace.scratchPtr
      inputHeap (SparsePolyZp.toPoly this._p.toNat left)
      (SparsePolyZp.toPoly this._p.toNat right) hcfg
      (Fact.out : Nat.Prime this._p.toNat).one_lt hleftLengthPositive
      hrightLengthPositive horder workspace.leftLengthWord houtputValid
      hscratchValid workspace.outputLeftDisjoint
      workspace.outputRightDisjoint workspace.outputScratchDisjoint
      workspace.scratchLeftDisjoint workspace.scratchRightDisjoint
      hleftDense.2 hrightRep.dense with
    ⟨outputHeap, hmulRun, hmulLayout, hproductRep⟩
  rcases dense_upoly_zp_to_upoly_raw_ir_refines this outputHeap
      workspace.outputPtr
      (sparseDenseLength left + sparseDenseLength right - 1)
      (SparsePolyZp.toPoly this._p.toNat left *
        SparsePolyZp.toPoly this._p.toNat right) hproductRep with
    ⟨product, hproductRun, hproductSemantic⟩
  rcases dense_upoly_zp_to_upoly_raw_ir_canonical this outputHeap
      workspace.outputPtr
      (sparseDenseLength left + sparseDenseLength right - 1)
      hproductRep.1 hproductRep.2.1 with
    ⟨canonicalProduct, hcanonicalRun, hcanonical⟩
  have hsame : canonicalProduct = product :=
    Except.ok.inj (hcanonicalRun.symm.trans hproductRun)
  subst canonicalProduct
  refine ⟨product, ?_, hcanonical, hproductSemantic⟩
  simp [strictMulOrderedIR, hleftRun, hrightRun, hmulRun, hproductRun]

/-- A physical allocator for multiplication calls whose operands are known
nonempty.  It cannot inspect or prescribe the semantic product. -/
structure RawMulWorkspaceProvider (this : DenseUPolyZp) where
  workspace : ∀ (left right : SparsePolyZp), 0 < left.size →
    0 < right.size → RawMulWorkspace this left right

/-- Total sparse multiplication boundary used by strict powmod.  The empty
branches are the representation-level short circuit of sparse multiplication;
every nonempty branch executes the generated dense `_mul`. -/
def strictMulIR (this : DenseUPolyZp) (left right : SparsePolyZp)
    (provider : RawMulWorkspaceProvider this) : RawExec SparsePolyZp :=
  if hleft : 0 < left.size then
    if hright : 0 < right.size then
      if sparseDenseLength right ≤ sparseDenseLength left then
        strictMulOrderedIR this left right
          (provider.workspace left right hleft hright)
      else
        strictMulOrderedIR this right left
          (provider.workspace right left hright hleft)
    else
      .ok #[]
  else
    .ok #[]

/-- Total semantic refinement of the representation-level multiplication
dispatcher. -/
theorem strictMulIR_refines_mul
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (left right : SparsePolyZp) (provider : RawMulWorkspaceProvider this)
    (hleftCanonical : SparsePolyZp.Canonical this._p.toNat left)
    (hrightCanonical : SparsePolyZp.Canonical this._p.toNat right) :
    ∃ product,
      strictMulIR this left right provider = .ok product ∧
      SparsePolyZp.Canonical this._p.toNat product ∧
      SparsePolyZp.toPoly this._p.toNat product =
        SparsePolyZp.toPoly this._p.toNat left *
          SparsePolyZp.toPoly this._p.toNat right := by
  by_cases hleft : 0 < left.size
  · by_cases hright : 0 < right.size
    · by_cases horder : sparseDenseLength right ≤ sparseDenseLength left
      · rcases strictMulOrderedIR_refines_mul this hcfg left right
          (provider.workspace left right hleft hright) hleftCanonical
          hrightCanonical hleft hright horder with
          ⟨product, hrun, hcanonical, hsemantic⟩
        exact ⟨product, by simp [strictMulIR, hleft, hright, horder, hrun],
          hcanonical, hsemantic⟩
      · have hreverse : sparseDenseLength left ≤ sparseDenseLength right :=
          Nat.le_of_not_ge horder
        rcases strictMulOrderedIR_refines_mul this hcfg right left
            (provider.workspace right left hright hleft) hrightCanonical
            hleftCanonical hright hleft hreverse with
          ⟨product, hrun, hcanonical, hsemantic⟩
        refine ⟨product, by simp [strictMulIR, hleft, hright, horder, hrun],
          hcanonical, ?_⟩
        simpa [mul_comm] using hsemantic
    · have hrightSize : right.size = 0 := Nat.eq_zero_of_not_pos hright
      have hrightEmpty : right = #[] := Array.size_eq_zero_iff.mp hrightSize
      refine ⟨#[], by simp [strictMulIR, hleft, hright], ?_, ?_⟩
      · refine ⟨?_, ?_, ?_⟩
        · intro x hx; simp at hx
        · simp
        · intro x hx; simp at hx
      · subst right
        simp [SparsePolyZp.toPoly]
  · have hleftSize : left.size = 0 := Nat.eq_zero_of_not_pos hleft
    have hleftEmpty : left = #[] := Array.size_eq_zero_iff.mp hleftSize
    refine ⟨#[], by simp [strictMulIR, hleft], ?_, ?_⟩
    · refine ⟨?_, ?_, ?_⟩
      · intro x hx; simp at hx
      · simp
      · intro x hx; simp at hx
    · subst left
      simp [SparsePolyZp.toPoly]

/-- Physical allocation for repeated reductions by one fixed modulus. -/
structure RawModWorkspaceProvider (this : DenseUPolyZp)
    (modulus : SparsePolyZp) where
  workspace : ∀ dividend : SparsePolyZp,
    RawModWorkspace this dividend modulus

/-- One actual multiply-then-reduce step. -/
def strictMulModIR (this : DenseUPolyZp) (left right modulus : SparsePolyZp)
    (mulProvider : RawMulWorkspaceProvider this)
    (modProvider : RawModWorkspaceProvider this modulus) :
    RawExec SparsePolyZp := do
  let product ← strictMulIR this left right mulProvider
  strictModIR this product modulus (modProvider.workspace product)

/-- Semantic composition of one actual raw multiply and one actual raw
division/remainder call. -/
theorem strictMulModIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (left right modulus : SparsePolyZp)
    (mulProvider : RawMulWorkspaceProvider this)
    (modProvider : RawModWorkspaceProvider this modulus)
    (hleftCanonical : SparsePolyZp.Canonical this._p.toNat left)
    (hrightCanonical : SparsePolyZp.Canonical this._p.toNat right)
    (hmodulusCanonical : SparsePolyZp.Canonical this._p.toNat modulus)
    (hmodulusNonempty : 0 < modulus.size)
    (hmodulusMonic : (SparsePolyZp.toPoly this._p.toNat modulus).Monic) :
    ∃ reduced,
      strictMulModIR this left right modulus mulProvider modProvider =
        .ok reduced ∧
      SparsePolyZp.Canonical this._p.toNat reduced ∧
      SparsePolyZp.toPoly this._p.toNat reduced =
        (SparsePolyZp.toPoly this._p.toNat left *
          SparsePolyZp.toPoly this._p.toNat right) %ₘ
            SparsePolyZp.toPoly this._p.toNat modulus := by
  rcases strictMulIR_refines_mul this hcfg left right mulProvider
      hleftCanonical hrightCanonical with
    ⟨product, hmul, hproductCanonical, hproductSemantic⟩
  rcases strictModIR_refines_modByMonic this hcfg product modulus
      (modProvider.workspace product) hproductCanonical hmodulusCanonical
      hmodulusNonempty hmodulusMonic with
    ⟨reduced, hmod, hreducedCanonical, hreducedSemantic⟩
  refine ⟨reduced, ?_, hreducedCanonical, ?_⟩
  · rw [strictMulModIR, hmul]
    exact hmod
  · rw [hreducedSemantic, hproductSemantic]

/-- The source binary-powmod loop, expressed as well-founded recursion on the
natural exponent.  Its only arithmetic calls are the strict raw multiplication
and raw division/remainder boundaries above. -/
def strictPowmodLoopIR (this : DenseUPolyZp) (modulus : SparsePolyZp)
    (mulProvider : RawMulWorkspaceProvider this)
    (modProvider : RawModWorkspaceProvider this modulus)
    (e : Nat) (base result : SparsePolyZp) : RawExec SparsePolyZp :=
  if hzero : e = 0 then
    .ok result
  else do
    let result' ← if e % 2 ≠ 0 then
        strictMulModIR this result base modulus mulProvider modProvider
      else
        .ok result
    let e' := e / 2
    let base' ← if 0 < e' then
        strictMulModIR this base base modulus mulProvider modProvider
      else
        .ok base
    strictPowmodLoopIR this modulus mulProvider modProvider e' base' result'
termination_by e
decreasing_by
  exact Nat.div_lt_self (Nat.pos_of_ne_zero hzero) (by omega)

/-- Actual `__upoly_powmod` entry specialized to a natural exponent, as used
by DDF.  The singleton one uses the generated `__make_zp` conversion and the
prime stored in the nonempty source modulus. -/
def strictPowmodIR (this : DenseUPolyZp) (base : SparsePolyZp) (e : Nat)
    (modulus : SparsePolyZp) (mulProvider : RawMulWorkspaceProvider this)
    (modProvider : RawModWorkspaceProvider this modulus) :
    RawExec SparsePolyZp :=
  if hmodulus : 0 < modulus.size then
    let prime := modulus[0].2.prime
    let one : SparsePolyZp := #[(UMonomial.mk 0,
      Generated.StrictDDF.__make_zp_ir 1 prime)]
    match strictModIR this base modulus (modProvider.workspace base) with
    | .error fault => .error fault
    | .ok reducedBase => strictPowmodLoopIR this modulus mulProvider
        modProvider e reducedBase one
  else
    .error .assertionFailure

/-- Physical providers and termination invariants used to instantiate the
cpp2lean-generated effectful DDF shell.  The fields select raw workspaces;
they do not contain factorization results. -/
structure DDFRawProviders (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] where
  hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this
  h2p : 2 * this._p.toNat ≤ UInt64.size
  mul : RawMulWorkspaceProvider this
  mod : ∀ modulus, RawModWorkspaceProvider this modulus
  gcd : ∀ left right, RawGCDWorkspace this left right

/-- Decode the concrete C++ accumulator into the L2 list consumed by
`ddfLoop_correct`. -/
noncomputable def ddfResultToL2 (p : Nat)
    (result : Array (SparsePolyZp × UInt64)) :
    List (Polynomial (ZMod p) × Nat) :=
  result.toList.map fun item =>
    (SparsePolyZp.toPoly p item.1, item.2.toNat)

/-- Every stored term of a canonical sparse polynomial occurs as a nonzero
coefficient of its L2 denotation, hence lies below its natural degree. -/
theorem canonical_term_degree_le_natDegree (p : Nat) [Fact (Nat.Prime p)]
    (poly : SparsePolyZp) (hcanonical : SparsePolyZp.Canonical p poly) :
    ∀ term ∈ poly.toList,
      term.1.deg ≤ (SparsePolyZp.toPoly p poly).natDegree := by
  intro term hterm
  have hcoefficient :
      (SparsePolyZp.toPoly p poly).coeff term.1.deg ≠ 0 := by
    unfold SparsePolyZp.toPoly
    rw [Refinement.StrictSquarefreeZp.listSum_coeff_of_mem_chain
      p poly.toList term
      hcanonical.2.1 hterm]
    exact Zp.toZMod_ne_zero_of_val_ne_zero p term.2
      (hcanonical.1 term hterm) (hcanonical.2.2 term hterm)
  exact Polynomial.le_natDegree_of_ne_zero hcoefficient

/-- A nonzero canonical sparse polynomial inherits a uniform stored-degree
bound as an L2 natural-degree bound. -/
theorem canonical_natDegree_lt_of_terms_lt (p : Nat) [Fact (Nat.Prime p)]
    (poly : SparsePolyZp) (hcanonical : SparsePolyZp.Canonical p poly)
    (hnonzero : SparsePolyZp.toPoly p poly ≠ 0) (bound : Nat)
    (hterms : ∀ term ∈ poly.toList, term.1.deg < bound) :
    (SparsePolyZp.toPoly p poly).natDegree < bound := by
  have hnonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      p poly hnonzero
  have hheadMem : poly[0] ∈ poly.toList := by
    simpa using Array.getElem_mem poly 0 hnonempty
  have hdegree :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
      p poly hcanonical hnonempty
  have hnatDegree : (SparsePolyZp.toPoly p poly).natDegree =
      poly[0].1.deg := Polynomial.natDegree_eq_of_degree_eq_some hdegree
  rw [hnatDegree]
  exact hterms poly[0] hheadMem

/-- On a canonical sparse polynomial whose degree fits signed 64 bits, the
C++ `get_deg` positivity test is exactly L2 positive natural degree. -/
theorem strict_get_deg_pos_iff_natDegree_pos (p : Nat) [Fact (Nat.Prime p)]
    (poly : SparsePolyZp) (hcanonical : SparsePolyZp.Canonical p poly)
    (hdegree : (SparsePolyZp.toPoly p poly).natDegree < 2 ^ 63) :
    get_deg poly > 0 ↔ 0 < (SparsePolyZp.toPoly p poly).natDegree := by
  by_cases hempty : poly.isEmpty
  · have hsize : poly.size = 0 := by simpa [Array.isEmpty] using hempty
    have hpoly : poly = #[] := Array.size_eq_zero_iff.mp hsize
    subst poly
    simp [get_deg]
  · have hnonempty : 0 < poly.size := by
      apply Nat.pos_of_ne_zero
      simpa [Array.isEmpty] using hempty
    have hdegreeHead :=
      Refinement.StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
        p poly hcanonical hnonempty
    have hpolyNonzero : SparsePolyZp.toPoly p poly ≠ 0 := by
      intro hzero
      rw [hzero] at hdegreeHead
      simp at hdegreeHead
    have hnatDegree : (SparsePolyZp.toPoly p poly).natDegree =
        poly[0].1.deg := by
      rw [Polynomial.degree_eq_natDegree hpolyNonzero] at hdegreeHead
      exact WithBot.coe_eq_coe.mp hdegreeHead
    rw [hnatDegree] at hdegree
    exact Refinement.StrictSquarefreeZp.generated_getDeg_pos_iff_natDegree_pos
      p poly hcanonical hnonempty hdegree

/-- Exact natural value of the signed C++ degree word on a nonempty canonical
polynomial. -/
theorem strict_get_deg_toNatClampNeg (p : Nat) [Fact (Nat.Prime p)]
    (poly : SparsePolyZp) (hcanonical : SparsePolyZp.Canonical p poly)
    (hnonempty : 0 < poly.size)
    (hdegree : (SparsePolyZp.toPoly p poly).natDegree < 2 ^ 63) :
    (get_deg poly).toNatClampNeg =
      (SparsePolyZp.toPoly p poly).natDegree := by
  have hne : poly ≠ #[] := by
    intro hempty
    subst poly
    simp at hnonempty
  have hdegreeEq :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_toPoly_degree_eq_head
      p poly hcanonical hnonempty
  have hnatDegree : (SparsePolyZp.toPoly p poly).natDegree =
      poly[0].1.deg := Polynomial.natDegree_eq_of_degree_eq_some hdegreeEq
  have hwordNat : poly[0].1.deg.toUInt64.toNat = poly[0].1.deg := by
    change (OfNat.ofNat poly[0].1.deg : UInt64).toNat = poly[0].1.deg
    rw [UInt64.toNat_ofNat, Nat.mod_eq_of_lt]
    omega
  have hheadBound : poly[0].1.deg < 2 ^ 63 := by
    rw [← hnatDegree]
    exact hdegree
  have hsignedNat : poly[0].1.deg.toUInt64.toInt64.toNatClampNeg =
      poly[0].1.deg := by
    rw [UInt64_toInt64_toNatClampNeg_eq_toNat_of_lt (by
      simpa [hwordNat] using hheadBound), hwordNat]
  have hget : get_deg poly = poly[0].1.deg.toUInt64.toInt64 := by
    simp [get_deg, Array.isEmpty_iff, hne, getElem!_pos poly 0 hnonempty]
  rw [hget, hsignedNat, hnatDegree]

/-- The source expression `(2 * d).toInt64` denotes mathematical `2*d`
under the DDF loop's explicit no-overflow budget. -/
theorem strict_twice_degree_toNatClampNeg (d : UInt64)
    (hbound : d.toNat < 2 ^ 62) :
    ((2 * d).toInt64).toNatClampNeg = 2 * d.toNat := by
  have hproduct64 : 2 * d.toNat < 2 ^ 64 := by omega
  have hproduct63 : 2 * d.toNat < 2 ^ 63 := by omega
  have htoNat : (2 * d).toNat = 2 * d.toNat := by
    rw [UInt64.toNat_mul]
    norm_num
    exact Nat.mod_eq_of_lt hproduct64
  rw [UInt64_toInt64_toNatClampNeg_eq_toNat_of_lt (by
    rw [htoNat]
    exact hproduct63), htoNat]

/-- Taking the generated recursive branch implies the corresponding L2
degree guard is false. -/
theorem strict_ddf_active_degree_ge (p : Nat) [Fact (Nat.Prime p)]
    (d : UInt64) (fStar : SparsePolyZp)
    (hd : d.toNat < 2 ^ 62)
    (hdPos : 0 < d.toNat)
    (hfCanonical : SparsePolyZp.Canonical p fStar)
    (hfNonempty : 0 < fStar.size)
    (hfDegree : (SparsePolyZp.toPoly p fStar).natDegree < 2 ^ 63)
    (hactive : ¬get_deg fStar < (2 * d).toInt64) :
    2 * d.toNat ≤ (SparsePolyZp.toPoly p fStar).natDegree := by
  have hleft := strict_get_deg_toNatClampNeg p fStar hfCanonical
    hfNonempty hfDegree
  have hright := strict_twice_degree_toNatClampNeg d hd
  have hactiveInt : ¬(get_deg fStar).toInt < (2 * d).toInt64.toInt := by
    simpa only [Int64.lt_iff_toInt_lt] using hactive
  have hleInt : (2 * d).toInt64.toInt ≤ (get_deg fStar).toInt :=
    Int.not_lt.mp hactiveInt
  change (get_deg fStar).toInt.toNat =
    (SparsePolyZp.toPoly p fStar).natDegree at hleft
  change (2 * d).toInt64.toInt.toNat = 2 * d.toNat at hright
  have hrightIntNonneg : 0 ≤ (2 * d).toInt64.toInt := by
    have hrightNatPos : 0 < (2 * d).toInt64.toInt.toNat := by
      rw [hright]
      omega
    omega
  have hleftIntNonneg : 0 ≤ (get_deg fStar).toInt :=
    le_trans hrightIntNonneg hleInt
  omega

/-- Taking the generated termination branch implies exactly the corresponding
L2 degree bound.  This is the converse machine-boundary bridge needed to align
the two well-founded loop equations; it does not appeal to either algorithm's
result. -/
theorem strict_ddf_term_degree_lt (p : Nat) [Fact (Nat.Prime p)]
    (d : UInt64) (fStar : SparsePolyZp)
    (hd : d.toNat < 2 ^ 62)
    (hdPos : 0 < d.toNat)
    (hfCanonical : SparsePolyZp.Canonical p fStar)
    (hfNonempty : 0 < fStar.size)
    (hfDegree : (SparsePolyZp.toPoly p fStar).natDegree < 2 ^ 63)
    (hterm : get_deg fStar < (2 * d).toInt64) :
    (SparsePolyZp.toPoly p fStar).natDegree < 2 * d.toNat := by
  have hleft := strict_get_deg_toNatClampNeg p fStar hfCanonical
    hfNonempty hfDegree
  have hright := strict_twice_degree_toNatClampNeg d hd
  have htermInt : (get_deg fStar).toInt < (2 * d).toInt64.toInt := by
    simpa only [Int64.lt_iff_toInt_lt] using hterm
  change (get_deg fStar).toInt.toNat =
    (SparsePolyZp.toPoly p fStar).natDegree at hleft
  change (2 * d).toInt64.toInt.toNat = 2 * d.toNat at hright
  have hrightIntPos : 0 < (2 * d).toInt64.toInt := by
    have hrightNatPos : 0 < (2 * d).toInt64.toInt.toNat := by
      rw [hright]
      omega
    omega
  omega

/-- Full representation-and-mathematical invariant for the generated DDF
loop.  The final seven fields are exactly P0–P6 of `ddfLoop_correct`; the
preceding fields justify every machine-word and sparse-array operation. -/
structure DDFLoopInvariant (this : DenseUPolyZp)
    (original : Polynomial (ZMod this._p.toNat))
    (d : UInt64) (fStar h : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64)) (prime : UInt64) : Prop where
  prime_eq : prime = this._p
  d_bound : d.toNat < 2 ^ 62
  fStar_canonical : SparsePolyZp.Canonical this._p.toNat fStar
  h_canonical : SparsePolyZp.Canonical this._p.toNat h
  fStar_degree_bound : ∀ term ∈ fStar.toList, term.1.deg < 2 ^ 62
  h_degree_bound : ∀ term ∈ h.toList, term.1.deg < 2 ^ 62
  result_canonical : ∀ item ∈ result.toList,
    SparsePolyZp.Canonical this._p.toNat item.1
  p0 : 1 ≤ d.toNat
  p1 : original = SparsePolyZp.toPoly this._p.toNat fStar *
    ((ddfResultToL2 this._p.toNat result).map Prod.fst).prod
  p2 : SparsePolyZp.toPoly this._p.toNat fStar ∣
    (SparsePolyZp.toPoly this._p.toNat h -
      Polynomial.X ^ (this._p.toNat ^ (d.toNat - 1)))
  p3 : (SparsePolyZp.toPoly this._p.toNat fStar).Monic
  p4 : Squarefree (SparsePolyZp.toPoly this._p.toNat fStar)
  p5 : ∀ q : Polynomial (ZMod this._p.toNat),
    Irreducible q → q ∣ SparsePolyZp.toPoly this._p.toNat fStar →
      d.toNat ≤ q.natDegree
  p6 : ∀ item ∈ ddfResultToL2 this._p.toNat result,
    item.1 ∣ original ∧ item.1.Monic ∧
      (∀ q : Polynomial (ZMod this._p.toNat),
        Irreducible q → q ∣ item.1 → q.natDegree = item.2)

/-- The exact C++ initial state `d = 1`, `fStar = f`, `h = X`, empty
accumulator satisfies every representation condition and P0–P6. -/
theorem DDFLoopInvariant.initial
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (f : SparsePolyZp) (q : UInt64) (hq : q = this._p)
    (hfCanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hfDegree : ∀ term ∈ f.toList, term.1.deg < 2 ^ 62)
    (hfMonic : (SparsePolyZp.toPoly this._p.toNat f).Monic)
    (hfSquarefree : Squarefree (SparsePolyZp.toPoly this._p.toNat f)) :
    DDFLoopInvariant this (SparsePolyZp.toPoly this._p.toNat f) 1 f
      #[(UMonomial.mk (1 : Int64),
        Generated.StrictDDF.__make_zp_ir (1 : Int64) q)] #[] q := by
  have hqNat : q.toNat = this._p.toNat := by rw [hq]
  have hx := Refinement.strict_singleton_x_data
    (p := this._p.toNat) q hqNat
  have hxCanonical : SparsePolyZp.Canonical this._p.toNat
      (#[(UMonomial.mk (1 : Int64),
        Generated.StrictDDF.__make_zp_ir (1 : Int64) q)] : SparsePolyZp) := by
    apply (canonicalRep_iff_canonical _).mp
    simpa [Generated.StrictDDF.__make_zp_ir] using hx.1
  have hxSemantic : SparsePolyZp.toPoly this._p.toNat
      (#[(UMonomial.mk (1 : Int64),
        Generated.StrictDDF.__make_zp_ir (1 : Int64) q)] : SparsePolyZp) =
      Polynomial.X := by
    simpa [Generated.StrictDDF.__make_zp_ir] using hx.2
  refine ⟨hq, ?_, hfCanonical, hxCanonical, hfDegree, ?_, ?_,
    ?_, ?_, ?_, hfMonic, hfSquarefree, ?_, ?_⟩
  · change (1 : Nat) < 2 ^ 62
    norm_num
  · intro term hterm
    have htermEq : term =
        (UMonomial.mk (1 : Int64),
          Generated.StrictDDF.__make_zp_ir (1 : Int64) q) := by
      simpa using hterm
    subst term
    change (1 : Nat) < 2 ^ 62
    norm_num
  · intro item hitem
    simp at hitem
  · change (1 : Nat) ≤ 1
    exact le_rfl
  · simp [ddfResultToL2]
  · rw [hxSemantic]
    simp
  · intro irreducible hirreducible _
    exact Irreducible.natDegree_pos hirreducible
  · intro item hitem
    simp [ddfResultToL2] at hitem

/-- The concrete no-split recursive call advances the unsigned counter
without wraparound and strictly decreases the generated raw measure. -/
theorem strict_ddf_noSplit_machine_step
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (original : Polynomial (ZMod this._p.toNat))
    (d : UInt64) (fStar h : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64)) (prime : UInt64)
    (hinvariant : DDFLoopInvariant this original d fStar h result prime)
    (hactive : ¬get_deg fStar < (2 * d).toInt64) :
    (d + 1).toNat = d.toNat + 1 ∧
      (d + 1).toNat < 2 ^ 62 ∧
      Generated.StrictDDF.ddfRawMeasure fStar (d + 1) <
        Generated.StrictDDF.ddfRawMeasure fStar d := by
  have hfNonzero : SparsePolyZp.toPoly this._p.toNat fStar ≠ 0 :=
    hinvariant.p3.ne_zero
  have hfNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat fStar hfNonzero
  have hfDegree : (SparsePolyZp.toPoly this._p.toNat fStar).natDegree <
      2 ^ 62 := canonical_natDegree_lt_of_terms_lt this._p.toNat fStar
        hinvariant.fStar_canonical hfNonzero (2 ^ 62)
        hinvariant.fStar_degree_bound
  have hactiveNat := strict_ddf_active_degree_ge this._p.toNat d fStar
    hinvariant.d_bound (lt_of_lt_of_le hinvariant.p0 (Nat.le_refl _))
    hinvariant.fStar_canonical hfNonempty (lt_trans hfDegree (by norm_num))
    hactive
  have hdIncrement :=
    Refinement.StrictSquarefreeZp.UInt64_toNat_add_one_of_lt d (by
      have : d.toNat + 1 < UInt64.size := by
        simpa [UInt64.size] using
          (show d.toNat + 1 < 2 ^ 64 by omega)
      exact this)
  have hdNextBound : (d + 1).toNat < 2 ^ 62 := by
    rw [hdIncrement]
    omega
  have hfNotEmpty : fStar.isEmpty = false := by
    simp [Array.isEmpty, Nat.ne_of_gt hfNonempty]
  have hgetDegree := strict_get_deg_toNatClampNeg this._p.toNat fStar
    hinvariant.fStar_canonical hfNonempty (lt_trans hfDegree (by norm_num))
  refine ⟨hdIncrement, hdNextBound, ?_⟩
  simp only [Generated.StrictDDF.ddfRawMeasure, hfNotEmpty,
    Bool.false_eq_true, ↓reduceIte, Int64.toNat]
  rw [hgetDegree, hdIncrement]
  omega

/-- Full no-split transition of the generated DDF loop.  Its hypotheses are
exactly the successful raw-call observations supplied by `DDFRawOps`; the
conclusion re-establishes the concrete representation conditions and L2
P0–P6, together with the generated well-founded decrease. -/
theorem DDFLoopInvariant.noSplit
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (original : Polynomial (ZMod this._p.toNat))
    (d : UInt64) (fStar h : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64)) (prime : UInt64)
    (hPow gdRaw : SparsePolyZp)
    (hinvariant : DDFLoopInvariant this original d fStar h result prime)
    (hactive : ¬get_deg fStar < (2 * d).toInt64)
    (hPowCanonical : SparsePolyZp.Canonical this._p.toNat hPow)
    (hPowSemantic : SparsePolyZp.toPoly this._p.toNat hPow =
      SparsePolyZp.toPoly this._p.toNat h ^ this._p.toNat %ₘ
        SparsePolyZp.toPoly this._p.toNat fStar)
    (hgdCanonical : SparsePolyZp.Canonical this._p.toNat gdRaw)
    (hgdSemantic : SparsePolyZp.toPoly this._p.toNat gdRaw = normalize
      (EuclideanDomain.gcd
        (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
        (SparsePolyZp.toPoly this._p.toNat fStar)))
    (hnoSplit : (!gdRaw.isEmpty && get_deg gdRaw > 0) = false) :
    DDFLoopInvariant this original (d + 1) fStar hPow result prime ∧
      Generated.StrictDDF.ddfRawMeasure fStar (d + 1) <
        Generated.StrictDDF.ddfRawMeasure fStar d := by
  rcases strict_ddf_noSplit_machine_step this original d fStar h result
      prime hinvariant hactive with ⟨hdIncrement, hdNextBound, hmeasure⟩
  have hfNonzero : SparsePolyZp.toPoly this._p.toNat fStar ≠ 0 :=
    hinvariant.p3.ne_zero
  have hfDegree : (SparsePolyZp.toPoly this._p.toNat fStar).natDegree <
      2 ^ 62 := canonical_natDegree_lt_of_terms_lt this._p.toNat fStar
        hinvariant.fStar_canonical hfNonzero (2 ^ 62)
        hinvariant.fStar_degree_bound
  have hPowDegreeLe :
      (SparsePolyZp.toPoly this._p.toNat hPow).natDegree ≤
        (SparsePolyZp.toPoly this._p.toNat fStar).natDegree := by
    rw [hPowSemantic]
    exact Polynomial.natDegree_modByMonic_le _ hinvariant.p3
  have hPowDegree :
      (SparsePolyZp.toPoly this._p.toNat hPow).natDegree < 2 ^ 62 :=
    lt_of_le_of_lt hPowDegreeLe hfDegree
  have hPowTerms : ∀ term ∈ hPow.toList, term.1.deg < 2 ^ 62 := by
    intro term hterm
    exact lt_of_le_of_lt
      (canonical_term_degree_le_natDegree this._p.toNat hPow
        hPowCanonical term hterm) hPowDegree
  have hcongNext : SparsePolyZp.toPoly this._p.toNat fStar ∣
      (SparsePolyZp.toPoly this._p.toNat hPow -
        Polynomial.X ^ (this._p.toNat ^ ((d + 1).toNat - 1))) := by
    rw [hdIncrement, Nat.add_sub_cancel]
    rw [hPowSemantic]
    exact ddf_h_cong_step
      (SparsePolyZp.toPoly this._p.toNat h)
      (SparsePolyZp.toPoly this._p.toNat fStar) d.toNat
      hinvariant.p0 hinvariant.p3 hinvariant.p2
  have hgcdRawNonzero : EuclideanDomain.gcd
      (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
      (SparsePolyZp.toPoly this._p.toNat fStar) ≠ 0 := by
    intro hzero
    have hdiv := EuclideanDomain.gcd_dvd_right
      (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
      (SparsePolyZp.toPoly this._p.toNat fStar)
    rw [hzero] at hdiv
    exact hfNonzero (zero_dvd_iff.mp hdiv)
  have hgdNonzero : SparsePolyZp.toPoly this._p.toNat gdRaw ≠ 0 := by
    rw [hgdSemantic]
    exact mt normalize_eq_zero.mp hgcdRawNonzero
  have hgdDvd : SparsePolyZp.toPoly this._p.toNat gdRaw ∣
      SparsePolyZp.toPoly this._p.toNat fStar := by
    rw [hgdSemantic]
    exact normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right _ _)
  have hgdDegreeLe :
      (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree ≤
        (SparsePolyZp.toPoly this._p.toNat fStar).natDegree :=
    Polynomial.natDegree_le_of_dvd hgdDvd hfNonzero
  have hgdDegree :
      (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree < 2 ^ 63 :=
    lt_of_le_of_lt hgdDegreeLe (lt_trans hfDegree (by norm_num))
  have hgdDegreeZero :
      (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree = 0 := by
    apply Nat.eq_zero_of_not_pos
    intro hpositive
    have hgetPositive :=
      (strict_get_deg_pos_iff_natDegree_pos this._p.toNat gdRaw
        hgdCanonical hgdDegree).mpr hpositive
    have hgdNonempty :=
      Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
        this._p.toNat gdRaw hgdNonzero
    have hnotEmpty : gdRaw.isEmpty = false := by
      simp [Array.isEmpty, Nat.ne_of_gt hgdNonempty]
    simp [hnotEmpty, hgetPositive] at hnoSplit
  have hp5Next : ∀ q : Polynomial (ZMod this._p.toNat),
      Irreducible q → q ∣ SparsePolyZp.toPoly this._p.toNat fStar →
        (d + 1).toNat ≤ q.natDegree := by
    intro q hq hqfs
    rw [hdIncrement]
    by_contra hnot
    have hdegreeEq : q.natDegree = d.toNat := by
      have := hinvariant.p5 q hq hqfs
      omega
    have hcongD : SparsePolyZp.toPoly this._p.toNat fStar ∣
        (SparsePolyZp.toPoly this._p.toNat hPow -
          Polynomial.X ^ (this._p.toNat ^ d.toNat)) := by
      simpa [hdIncrement] using hcongNext
    have hqgd := (ddf_gd_irred_characterization
      (SparsePolyZp.toPoly this._p.toNat fStar)
      (SparsePolyZp.toPoly this._p.toNat hPow) d.toNat hinvariant.p0
      hinvariant.p3 hinvariant.p4 hcongD hinvariant.p5 q hq).mpr
        ⟨hqfs, hdegreeEq⟩
    have hqgdRaw : q ∣ SparsePolyZp.toPoly this._p.toNat gdRaw := by
      rw [hgdSemantic]
      exact hqgd
    have hle := Polynomial.natDegree_le_of_dvd hqgdRaw hgdNonzero
    rw [hgdDegreeZero] at hle
    have hqPositive := Irreducible.natDegree_pos hq
    omega
  refine ⟨⟨hinvariant.prime_eq, hdNextBound,
    hinvariant.fStar_canonical, hPowCanonical,
    hinvariant.fStar_degree_bound, hPowTerms, hinvariant.result_canonical,
    ?_, hinvariant.p1, hcongNext, hinvariant.p3, hinvariant.p4,
    hp5Next, hinvariant.p6⟩, hmeasure⟩
  rw [hdIncrement]
  omega

/-- Full split transition of the generated DDF loop after the concrete GCD,
make-monic, exact-division, normalization, and modular-reduction calls have
returned.  No L2 operation is used as an implementation; the semantic
equalities below are precisely the refinement results of those raw calls. -/
theorem DDFLoopInvariant.split
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (original : Polynomial (ZMod this._p.toNat))
    (d : UInt64) (fStar h : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64)) (prime : UInt64)
    (hPow gd fNext hNext : SparsePolyZp)
    (hinvariant : DDFLoopInvariant this original d fStar h result prime)
    (hactive : ¬get_deg fStar < (2 * d).toInt64)
    (hPowCanonical : SparsePolyZp.Canonical this._p.toNat hPow)
    (hPowSemantic : SparsePolyZp.toPoly this._p.toNat hPow =
      SparsePolyZp.toPoly this._p.toNat h ^ this._p.toNat %ₘ
        SparsePolyZp.toPoly this._p.toNat fStar)
    (hgdCanonical : SparsePolyZp.Canonical this._p.toNat gd)
    (hgdSemantic : SparsePolyZp.toPoly this._p.toNat gd = normalize
      (EuclideanDomain.gcd
        (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
        (SparsePolyZp.toPoly this._p.toNat fStar)))
    (hgdPositive : 0 < (SparsePolyZp.toPoly this._p.toNat gd).natDegree)
    (hfNextCanonical : SparsePolyZp.Canonical this._p.toNat fNext)
    (hfNextSemantic : SparsePolyZp.toPoly this._p.toNat fNext =
      SparsePolyZp.toPoly this._p.toNat fStar /ₘ
        SparsePolyZp.toPoly this._p.toNat gd)
    (hNextCanonical : SparsePolyZp.Canonical this._p.toNat hNext)
    (hNextSemantic : SparsePolyZp.toPoly this._p.toNat hNext =
      SparsePolyZp.toPoly this._p.toNat hPow %ₘ
        SparsePolyZp.toPoly this._p.toNat fNext) :
    DDFLoopInvariant this original (d + 1) fNext hNext
        (result.push (gd, d)) prime ∧
      Generated.StrictDDF.ddfRawMeasure fNext (d + 1) <
        Generated.StrictDDF.ddfRawMeasure fStar d := by
  have hfNonzero : SparsePolyZp.toPoly this._p.toNat fStar ≠ 0 :=
    hinvariant.p3.ne_zero
  have hfNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat fStar hfNonzero
  have hfDegree : (SparsePolyZp.toPoly this._p.toNat fStar).natDegree <
      2 ^ 62 := canonical_natDegree_lt_of_terms_lt this._p.toNat fStar
        hinvariant.fStar_canonical hfNonzero (2 ^ 62)
        hinvariant.fStar_degree_bound
  have hactiveNat := strict_ddf_active_degree_ge this._p.toNat d fStar
    hinvariant.d_bound (lt_of_lt_of_le hinvariant.p0 (Nat.le_refl _))
    hinvariant.fStar_canonical hfNonempty (lt_trans hfDegree (by norm_num))
    hactive
  have hdIncrement :=
    Refinement.StrictSquarefreeZp.UInt64_toNat_add_one_of_lt d (by
      simpa [UInt64.size] using
        (show d.toNat + 1 < 2 ^ 64 by omega))
  have hdNextBound : (d + 1).toNat < 2 ^ 62 := by
    rw [hdIncrement]
    omega
  have hcongD : SparsePolyZp.toPoly this._p.toNat fStar ∣
      (SparsePolyZp.toPoly this._p.toNat hPow -
        Polynomial.X ^ (this._p.toNat ^ d.toNat)) := by
    rw [hPowSemantic]
    exact ddf_h_cong_step
      (SparsePolyZp.toPoly this._p.toNat h)
      (SparsePolyZp.toPoly this._p.toNat fStar) d.toNat
      hinvariant.p0 hinvariant.p3 hinvariant.p2
  have hgcdRawNonzero : EuclideanDomain.gcd
      (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
      (SparsePolyZp.toPoly this._p.toNat fStar) ≠ 0 := by
    intro hzero
    have hdiv := EuclideanDomain.gcd_dvd_right
      (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
      (SparsePolyZp.toPoly this._p.toNat fStar)
    rw [hzero] at hdiv
    exact hfNonzero (zero_dvd_iff.mp hdiv)
  have hgdMonic : (SparsePolyZp.toPoly this._p.toNat gd).Monic := by
    rw [hgdSemantic]
    exact Polynomial.monic_normalize hgcdRawNonzero
  have hgdDvd : SparsePolyZp.toPoly this._p.toNat gd ∣
      SparsePolyZp.toPoly this._p.toNat fStar := by
    rw [hgdSemantic]
    exact normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right _ _)
  have hF3 : SparsePolyZp.toPoly this._p.toNat fStar =
      SparsePolyZp.toPoly this._p.toNat gd *
        SparsePolyZp.toPoly this._p.toNat fNext := by
    rw [hfNextSemantic]
    exact (ddf_monic_divByMonic_mul_eq
      (SparsePolyZp.toPoly this._p.toNat fStar)
      (SparsePolyZp.toPoly this._p.toNat gd) hgdMonic hgdDvd).symm
  have hfNextMonic :
      (SparsePolyZp.toPoly this._p.toNat fNext).Monic := by
    rw [hfNextSemantic]
    exact Refinement.StrictSquarefreeZp.divByMonic_monic_of_monic_of_dvd
      (SparsePolyZp.toPoly this._p.toNat fStar)
      (SparsePolyZp.toPoly this._p.toNat gd) hinvariant.p3 hgdMonic hgdDvd
  have hfNextNonzero : SparsePolyZp.toPoly this._p.toNat fNext ≠ 0 :=
    hfNextMonic.ne_zero
  have hfNextDegreeLt :
      (SparsePolyZp.toPoly this._p.toNat fNext).natDegree <
        (SparsePolyZp.toPoly this._p.toNat fStar).natDegree := by
    rw [hfNextSemantic]
    exact ddf_natDegree_divByMonic_lt
      (SparsePolyZp.toPoly this._p.toNat fStar)
      (SparsePolyZp.toPoly this._p.toNat gd) hgdMonic hgdDvd
      hgdPositive hfNonzero
  have hfNextDegree :
      (SparsePolyZp.toPoly this._p.toNat fNext).natDegree < 2 ^ 62 :=
    lt_trans hfNextDegreeLt hfDegree
  have hfNextTerms : ∀ term ∈ fNext.toList, term.1.deg < 2 ^ 62 := by
    intro term hterm
    exact lt_of_le_of_lt
      (canonical_term_degree_le_natDegree this._p.toNat fNext
        hfNextCanonical term hterm) hfNextDegree
  have hNextDegreeLe :
      (SparsePolyZp.toPoly this._p.toNat hNext).natDegree ≤
        (SparsePolyZp.toPoly this._p.toNat fNext).natDegree := by
    rw [hNextSemantic]
    exact Polynomial.natDegree_modByMonic_le _ hfNextMonic
  have hNextDegree :
      (SparsePolyZp.toPoly this._p.toNat hNext).natDegree < 2 ^ 62 :=
    lt_of_le_of_lt hNextDegreeLe hfNextDegree
  have hNextTerms : ∀ term ∈ hNext.toList, term.1.deg < 2 ^ 62 := by
    intro term hterm
    exact lt_of_le_of_lt
      (canonical_term_degree_le_natDegree this._p.toNat hNext
        hNextCanonical term hterm) hNextDegree
  have hprodNext : original = SparsePolyZp.toPoly this._p.toNat fNext *
      ((ddfResultToL2 this._p.toNat (result.push (gd, d))).map
        Prod.fst).prod := by
    have hdecodePush :
        ddfResultToL2 this._p.toNat (result.push (gd, d)) =
          ddfResultToL2 this._p.toNat result ++
            [(SparsePolyZp.toPoly this._p.toNat gd, d.toNat)] := by
      simp [ddfResultToL2, Array.toList_push]
    rw [hdecodePush]
    simp only [List.map_append, List.map_singleton, List.prod_append,
      List.prod_singleton]
    rw [hinvariant.p1, hF3]
    ring
  have hfNextDvdF : SparsePolyZp.toPoly this._p.toNat fNext ∣
      SparsePolyZp.toPoly this._p.toNat fStar := by
    exact ⟨SparsePolyZp.toPoly this._p.toNat gd, hF3.trans (mul_comm _ _)⟩
  have hcongNext : SparsePolyZp.toPoly this._p.toNat fNext ∣
      (SparsePolyZp.toPoly this._p.toNat hNext -
        Polynomial.X ^ (this._p.toNat ^ ((d + 1).toNat - 1))) := by
    rw [hdIncrement, Nat.add_sub_cancel, hNextSemantic]
    have hdivCong := dvd_trans hfNextDvdF hcongD
    have hmodDvd : SparsePolyZp.toPoly this._p.toNat fNext ∣
        (SparsePolyZp.toPoly this._p.toNat hPow -
          SparsePolyZp.toPoly this._p.toNat hPow %ₘ
            SparsePolyZp.toPoly this._p.toNat fNext) := by
      refine ⟨SparsePolyZp.toPoly this._p.toNat hPow /ₘ
        SparsePolyZp.toPoly this._p.toNat fNext, ?_⟩
      have := (Polynomial.modByMonic_add_div
        (SparsePolyZp.toPoly this._p.toNat hPow)
        (SparsePolyZp.toPoly this._p.toNat fNext)).symm
      rw [add_comm] at this
      exact sub_eq_of_eq_add this
    have hdecomp :
        (SparsePolyZp.toPoly this._p.toNat hPow %ₘ
            SparsePolyZp.toPoly this._p.toNat fNext -
          Polynomial.X ^ (this._p.toNat ^ d.toNat)) =
        -(SparsePolyZp.toPoly this._p.toNat hPow -
            SparsePolyZp.toPoly this._p.toNat hPow %ₘ
              SparsePolyZp.toPoly this._p.toNat fNext) +
          (SparsePolyZp.toPoly this._p.toNat hPow -
            Polynomial.X ^ (this._p.toNat ^ d.toNat)) := by ring
    rw [hdecomp]
    exact dvd_add (dvd_neg.mpr hmodDvd) hdivCong
  have hfNextSquarefree :
      Squarefree (SparsePolyZp.toPoly this._p.toNat fNext) :=
    Squarefree.squarefree_of_dvd hfNextDvdF hinvariant.p4
  have hp5Next : ∀ q : Polynomial (ZMod this._p.toNat),
      Irreducible q → q ∣ SparsePolyZp.toPoly this._p.toNat fNext →
        (d + 1).toNat ≤ q.natDegree := by
    intro q hq hqfn
    rw [hdIncrement]
    have hqfs := dvd_trans hqfn hfNextDvdF
    by_contra hnot
    have hdegreeEq : q.natDegree = d.toNat := by
      have := hinvariant.p5 q hq hqfs
      omega
    have hqgd := (ddf_gd_irred_characterization
      (SparsePolyZp.toPoly this._p.toNat fStar)
      (SparsePolyZp.toPoly this._p.toNat hPow) d.toNat hinvariant.p0
      hinvariant.p3 hinvariant.p4 hcongD hinvariant.p5 q hq).mpr
        ⟨hqfs, hdegreeEq⟩
    have hqgdRaw : q ∣ SparsePolyZp.toPoly this._p.toNat gd := by
      rw [hgdSemantic]
      exact hqgd
    have hqSquare : q * q ∣ SparsePolyZp.toPoly this._p.toNat fStar := by
      rw [hF3]
      exact mul_dvd_mul hqgdRaw hqfn
    exact absurd (hinvariant.p4 q hqSquare) hq.1
  have hp6Next : ∀ item ∈
      ddfResultToL2 this._p.toNat (result.push (gd, d)),
      item.1 ∣ original ∧ item.1.Monic ∧
        (∀ q : Polynomial (ZMod this._p.toNat),
          Irreducible q → q ∣ item.1 → q.natDegree = item.2) := by
    intro item hitem
    simp only [ddfResultToL2, Array.toList_push, List.map_append,
      List.map_singleton, List.mem_append, List.mem_singleton] at hitem
    rcases hitem with hOld | rfl
    · exact hinvariant.p6 item hOld
    · refine ⟨dvd_trans hgdDvd ⟨_, hinvariant.p1⟩, hgdMonic, ?_⟩
      intro q hq hqgd
      exact ((ddf_gd_irred_characterization
        (SparsePolyZp.toPoly this._p.toNat fStar)
        (SparsePolyZp.toPoly this._p.toNat hPow) d.toNat hinvariant.p0
        hinvariant.p3 hinvariant.p4 hcongD hinvariant.p5 q hq).mp (by
          rw [← hgdSemantic]
          exact hqgd)).2
  have hfNextNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat fNext hfNextNonzero
  have hfNotEmpty : fStar.isEmpty = false := by
    simp [Array.isEmpty, Nat.ne_of_gt hfNonempty]
  have hfNextNotEmpty : fNext.isEmpty = false := by
    simp [Array.isEmpty, Nat.ne_of_gt hfNextNonempty]
  have hgetF := strict_get_deg_toNatClampNeg this._p.toNat fStar
    hinvariant.fStar_canonical hfNonempty (lt_trans hfDegree (by norm_num))
  have hgetNext := strict_get_deg_toNatClampNeg this._p.toNat fNext
    hfNextCanonical hfNextNonempty (lt_trans hfNextDegree (by norm_num))
  have hmeasure : Generated.StrictDDF.ddfRawMeasure fNext (d + 1) <
      Generated.StrictDDF.ddfRawMeasure fStar d := by
    simp only [Generated.StrictDDF.ddfRawMeasure, hfNotEmpty,
      hfNextNotEmpty, Bool.false_eq_true, ↓reduceIte, Int64.toNat]
    rw [hgetF, hgetNext, hdIncrement]
    omega
  refine ⟨⟨hinvariant.prime_eq, hdNextBound, hfNextCanonical,
    hNextCanonical, hfNextTerms, hNextTerms, ?_, ?_, hprodNext,
    hcongNext, hfNextMonic, hfNextSquarefree, hp5Next, hp6Next⟩,
    hmeasure⟩
  · intro item hitem
    simp only [Array.toList_push, List.mem_append, List.mem_singleton] at hitem
    rcases hitem with hOld | rfl
    · exact hinvariant.result_canonical item hOld
    · exact hgdCanonical
  · rw [hdIncrement]
    omega

private lemma modByMonic_idem {p : Nat} [Fact (Nat.Prime p)]
    (a m : Polynomial (ZMod p))
    (hm : m.Monic) : (a %ₘ m) %ₘ m = a %ₘ m :=
  (Polynomial.modByMonic_eq_self_iff hm).mpr
    (Polynomial.degree_modByMonic_lt a hm)

private lemma mul_modByMonic_congr {p : Nat} [Fact (Nat.Prime p)]
    {a a' b b' m : Polynomial (ZMod p)}
    (ha : a %ₘ m = a' %ₘ m) (hb : b %ₘ m = b' %ₘ m) :
    (a * b) %ₘ m = (a' * b') %ₘ m := by
  calc
    (a * b) %ₘ m = (a %ₘ m) * (b %ₘ m) %ₘ m :=
      Polynomial.mul_modByMonic _ _ _
    _ = (a' %ₘ m) * (b' %ₘ m) %ₘ m := by rw [ha, hb]
    _ = (a' * b') %ₘ m := (Polynomial.mul_modByMonic _ _ _).symm

private lemma pow_modByMonic_congr {p : Nat} [Fact (Nat.Prime p)]
    {a b m : Polynomial (ZMod p)} (h : a %ₘ m = b %ₘ m) :
    ∀ n : Nat, (a ^ n) %ₘ m = (b ^ n) %ₘ m := by
  intro n
  induction n with
  | zero => simp
  | succ n ih =>
      rw [pow_succ, pow_succ]
      exact mul_modByMonic_congr ih h

/-- Strong-induction invariant for the actual well-founded binary-powmod
execution. -/
theorem strictPowmodLoopIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (modulus : SparsePolyZp) (mulProvider : RawMulWorkspaceProvider this)
    (modProvider : RawModWorkspaceProvider this modulus)
    (hmodulusCanonical : SparsePolyZp.Canonical this._p.toNat modulus)
    (hmodulusNonempty : 0 < modulus.size)
    (hmodulusMonic : (SparsePolyZp.toPoly this._p.toNat modulus).Monic) :
    ∀ (e : Nat) (base result : SparsePolyZp),
      SparsePolyZp.Canonical this._p.toNat base →
      SparsePolyZp.Canonical this._p.toNat result →
      SparsePolyZp.toPoly this._p.toNat base %ₘ
          SparsePolyZp.toPoly this._p.toNat modulus =
        SparsePolyZp.toPoly this._p.toNat base →
      SparsePolyZp.toPoly this._p.toNat result %ₘ
          SparsePolyZp.toPoly this._p.toNat modulus =
        SparsePolyZp.toPoly this._p.toNat result →
      ∃ output,
        strictPowmodLoopIR this modulus mulProvider modProvider e base result =
          .ok output ∧
        SparsePolyZp.Canonical this._p.toNat output ∧
        SparsePolyZp.toPoly this._p.toNat output =
          (SparsePolyZp.toPoly this._p.toNat result *
            SparsePolyZp.toPoly this._p.toNat base ^ e) %ₘ
              SparsePolyZp.toPoly this._p.toNat modulus := by
  intro e
  induction e using Nat.strongRecOn with
  | ind e ih =>
      intro base result hbaseCanonical hresultCanonical hbaseReduced
        hresultReduced
      by_cases hezero : e = 0
      · subst e
        refine ⟨result, by simp [strictPowmodLoopIR], hresultCanonical, ?_⟩
        simpa using hresultReduced.symm
      · have hepos : 0 < e := Nat.pos_of_ne_zero hezero
        let k := e / 2
        have hklt : k < e := Nat.div_lt_self hepos (by omega)
        by_cases hodd : e % 2 ≠ 0
        · have hoddOne : e % 2 = 1 := by omega
          rcases strictMulModIR_refines this hcfg result base modulus
              mulProvider modProvider hresultCanonical hbaseCanonical
              hmodulusCanonical hmodulusNonempty hmodulusMonic with
            ⟨result', hresultRun, hresult'Canonical, hresult'Semantic⟩
          have hresult'Reduced :
              SparsePolyZp.toPoly this._p.toNat result' %ₘ
                  SparsePolyZp.toPoly this._p.toNat modulus =
                SparsePolyZp.toPoly this._p.toNat result' := by
            rw [hresult'Semantic]
            exact modByMonic_idem _ _ hmodulusMonic
          by_cases hkpos : 0 < k
          · rcases strictMulModIR_refines this hcfg base base modulus
                mulProvider modProvider hbaseCanonical hbaseCanonical
                hmodulusCanonical hmodulusNonempty hmodulusMonic with
              ⟨base', hbaseRun, hbase'Canonical, hbase'Semantic⟩
            have hbase'Reduced :
                SparsePolyZp.toPoly this._p.toNat base' %ₘ
                    SparsePolyZp.toPoly this._p.toNat modulus =
                  SparsePolyZp.toPoly this._p.toNat base' := by
              rw [hbase'Semantic]
              exact modByMonic_idem _ _ hmodulusMonic
            rcases ih k hklt base' result' hbase'Canonical hresult'Canonical
                hbase'Reduced hresult'Reduced with
              ⟨output, hrecursiveRun, houtputCanonical, houtputSemantic⟩
            refine ⟨output, ?_, houtputCanonical, houtputSemantic.trans ?_⟩
            · rw [strictPowmodLoopIR]
              simp [hezero, hodd, hresultRun, k, hkpos, hbaseRun,
                ]
              simpa [k] using hrecursiveRun
            · have heform : e = 2 * k + 1 := by omega
              calc
                (SparsePolyZp.toPoly this._p.toNat result' *
                    SparsePolyZp.toPoly this._p.toNat base' ^ k) %ₘ
                      SparsePolyZp.toPoly this._p.toNat modulus =
                    ((SparsePolyZp.toPoly this._p.toNat result *
                        SparsePolyZp.toPoly this._p.toNat base) *
                      (SparsePolyZp.toPoly this._p.toNat base *
                        SparsePolyZp.toPoly this._p.toNat base) ^ k) %ₘ
                      SparsePolyZp.toPoly this._p.toNat modulus :=
                    mul_modByMonic_congr
                      (by rw [hresult'Semantic,
                        modByMonic_idem _ _ hmodulusMonic])
                      (pow_modByMonic_congr
                        (by rw [hbase'Semantic,
                          modByMonic_idem _ _ hmodulusMonic]) k)
                _ = (SparsePolyZp.toPoly this._p.toNat result *
                      SparsePolyZp.toPoly this._p.toNat base ^ e) %ₘ
                        SparsePolyZp.toPoly this._p.toNat modulus := by
                    rw [heform]
                    congr 1 <;> ring
          · have hkzero : k = 0 := Nat.eq_zero_of_not_pos hkpos
            have heone : e = 1 := by omega
            subst e
            refine ⟨result', ?_, hresult'Canonical, ?_⟩
            · rw [strictPowmodLoopIR]
              simp only [↓reduceIte, Nat.one_mod, one_ne_zero, hresultRun,
                Nat.reduceDiv, lt_self_iff_false]
              change strictPowmodLoopIR this modulus mulProvider modProvider
                0 base result' = .ok result'
              rw [strictPowmodLoopIR]
              rfl
            · simpa using hresult'Semantic
        · have heven : e % 2 = 0 := by omega
          have hkpos : 0 < k := by dsimp [k]; omega
          rcases strictMulModIR_refines this hcfg base base modulus
              mulProvider modProvider hbaseCanonical hbaseCanonical
              hmodulusCanonical hmodulusNonempty hmodulusMonic with
            ⟨base', hbaseRun, hbase'Canonical, hbase'Semantic⟩
          have hbase'Reduced :
              SparsePolyZp.toPoly this._p.toNat base' %ₘ
                  SparsePolyZp.toPoly this._p.toNat modulus =
                SparsePolyZp.toPoly this._p.toNat base' := by
            rw [hbase'Semantic]
            exact modByMonic_idem _ _ hmodulusMonic
          rcases ih k hklt base' result hbase'Canonical hresultCanonical
              hbase'Reduced hresultReduced with
            ⟨output, hrecursiveRun, houtputCanonical, houtputSemantic⟩
          refine ⟨output, ?_, houtputCanonical, houtputSemantic.trans ?_⟩
          · rw [strictPowmodLoopIR]
            simp [hezero, hodd, k, hkpos, hbaseRun]
            simpa [k] using hrecursiveRun
          · have heform : e = 2 * k := by omega
            calc
              (SparsePolyZp.toPoly this._p.toNat result *
                  SparsePolyZp.toPoly this._p.toNat base' ^ k) %ₘ
                    SparsePolyZp.toPoly this._p.toNat modulus =
                  (SparsePolyZp.toPoly this._p.toNat result *
                    (SparsePolyZp.toPoly this._p.toNat base *
                      SparsePolyZp.toPoly this._p.toNat base) ^ k) %ₘ
                    SparsePolyZp.toPoly this._p.toNat modulus :=
                  mul_modByMonic_congr rfl
                    (pow_modByMonic_congr
                      (by rw [hbase'Semantic,
                        modByMonic_idem _ _ hmodulusMonic]) k)
              _ = (SparsePolyZp.toPoly this._p.toNat result *
                    SparsePolyZp.toPoly this._p.toNat base ^ e) %ₘ
                      SparsePolyZp.toPoly this._p.toNat modulus := by
                  rw [heform]
                  congr 1 <;> ring

private lemma singletonOne_refines (q : UInt64)
    [Fact (Nat.Prime q.toNat)] :
    SparsePolyZp.Canonical q.toNat
        (#[(UMonomial.mk 0,
          Generated.StrictDDF.__make_zp_ir 1 q)] : SparsePolyZp) ∧
      SparsePolyZp.toPoly q.toNat
        (#[(UMonomial.mk 0,
          Generated.StrictDDF.__make_zp_ir 1 q)] : SparsePolyZp) = 1 := by
  have hqgt : 1 < q.toNat := (Fact.out : Nat.Prime q.toNat).one_lt
  have hmod : (1 : Int) % (q.toNat : Int) = 1 := by
    apply Int.emod_eq_of_lt
    · norm_num
    · exact_mod_cast hqgt
  simp [SparsePolyZp.Canonical, SparsePolyZp.WellFormed_arr,
    SparsePolyZp.AllReduced, SparsePolyZp.toPoly, listSum,
    Generated.StrictDDF.__make_zp_ir, Zp.Reduced, Zp.ofInt, Zp.toZMod,
    hmod, hqgt]

/-- End-to-end semantic refinement of the actual strict powmod entry. -/
theorem strictPowmodIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : CLPoly.Impl.StrictWordArithmetic.DensePreinvConfigured this)
    (base : SparsePolyZp) (e : Nat) (modulus : SparsePolyZp)
    (mulProvider : RawMulWorkspaceProvider this)
    (modProvider : RawModWorkspaceProvider this modulus)
    (hbaseCanonical : SparsePolyZp.Canonical this._p.toNat base)
    (hmodulusCanonical : SparsePolyZp.Canonical this._p.toNat modulus)
    (hmodulusNonempty : 0 < modulus.size)
    (hmodulusMonic : (SparsePolyZp.toPoly this._p.toNat modulus).Monic)
    (hmodulusDegree : 0 <
      (SparsePolyZp.toPoly this._p.toNat modulus).natDegree) :
    ∃ output,
      strictPowmodIR this base e modulus mulProvider modProvider = .ok output ∧
      SparsePolyZp.Canonical this._p.toNat output ∧
      SparsePolyZp.toPoly this._p.toNat output =
        SparsePolyZp.toPoly this._p.toNat base ^ e %ₘ
          SparsePolyZp.toPoly this._p.toNat modulus := by
  have hindex : 0 < modulus.size := hmodulusNonempty
  have hmember : modulus[0] ∈ modulus.toList := by
    simpa using Array.getElem_mem modulus 0 hindex
  have hprime : modulus[0].2.prime.toNat = this._p.toNat :=
    (hmodulusCanonical.1 modulus[0] hmember).1
  have hprimeWord : modulus[0].2.prime = this._p := by
    exact UInt64.toNat_inj.mp hprime
  rcases strictModIR_refines_modByMonic this hcfg base modulus
      (modProvider.workspace base) hbaseCanonical hmodulusCanonical
      hmodulusNonempty hmodulusMonic with
    ⟨reducedBase, hbaseRun, hreducedBaseCanonical,
      hreducedBaseSemantic⟩
  have hreducedBaseReduced :
      SparsePolyZp.toPoly this._p.toNat reducedBase %ₘ
          SparsePolyZp.toPoly this._p.toNat modulus =
        SparsePolyZp.toPoly this._p.toNat reducedBase := by
    rw [hreducedBaseSemantic]
    exact modByMonic_idem _ _ hmodulusMonic
  let one : SparsePolyZp := #[(UMonomial.mk 0,
    Generated.StrictDDF.__make_zp_ir 1 modulus[0].2.prime)]
  have honeData : SparsePolyZp.Canonical this._p.toNat one ∧
      SparsePolyZp.toPoly this._p.toNat one = 1 := by
    dsimp [one]
    rw [hprimeWord]
    exact singletonOne_refines this._p
  have honeReduced : SparsePolyZp.toPoly this._p.toNat one %ₘ
      SparsePolyZp.toPoly this._p.toNat modulus =
        SparsePolyZp.toPoly this._p.toNat one := by
    rw [honeData.2]
    apply (Polynomial.modByMonic_eq_self_iff hmodulusMonic).mpr
    have hmodulusNonzero :
        SparsePolyZp.toPoly this._p.toNat modulus ≠ 0 := by
      intro hzero
      simp [hzero] at hmodulusDegree
    simp [Polynomial.degree_eq_natDegree hmodulusNonzero, hmodulusDegree]
  rcases strictPowmodLoopIR_refines this hcfg modulus mulProvider modProvider
      hmodulusCanonical hmodulusNonempty hmodulusMonic e reducedBase one
      hreducedBaseCanonical honeData.1 hreducedBaseReduced honeReduced with
    ⟨output, hloopRun, houtputCanonical, houtputSemantic⟩
  refine ⟨output, ?_, houtputCanonical, houtputSemantic.trans ?_⟩
  · simp [strictPowmodIR, hmodulusNonempty, hbaseRun, one, hloopRun]
  · rw [honeData.2, one_mul]
    exact pow_modByMonic_congr
      (by rw [hreducedBaseSemantic,
        modByMonic_idem _ _ hmodulusMonic]) e

/-- Execute and refine the common powmod/subtract-X/public-GCD prefix of one
active generated DDF iteration.  Every returned witness is identified by an
actual `.ok` execution equation. -/
theorem strictDDFRawPrefix
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : DDFRawProviders this)
    (original : Polynomial (ZMod this._p.toNat))
    (d : UInt64) (fStar h : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64)) (prime : UInt64)
    (hinvariant : DDFLoopInvariant this original d fStar h result prime)
    (hactive : ¬get_deg fStar < (2 * d).toInt64) :
    ∃ hPow gdRaw,
      strictPowmodIR this h prime.toNat fStar providers.mul
          (providers.mod fStar) = .ok hPow ∧
      strictDDFGCDIR this
          (Generated.StrictDDF.__upoly_subtract_x_ir hPow prime) fStar
          (providers.gcd
            (Generated.StrictDDF.__upoly_subtract_x_ir hPow prime) fStar) =
        .ok gdRaw ∧
      SparsePolyZp.Canonical this._p.toNat hPow ∧
      SparsePolyZp.toPoly this._p.toNat hPow =
        SparsePolyZp.toPoly this._p.toNat h ^ this._p.toNat %ₘ
          SparsePolyZp.toPoly this._p.toNat fStar ∧
      SparsePolyZp.Canonical this._p.toNat gdRaw ∧
      SparsePolyZp.toPoly this._p.toNat gdRaw = normalize
        (EuclideanDomain.gcd
          (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
          (SparsePolyZp.toPoly this._p.toNat fStar)) := by
  have hprime : prime = this._p := hinvariant.prime_eq
  subst prime
  have hfNonzero : SparsePolyZp.toPoly this._p.toNat fStar ≠ 0 :=
    hinvariant.p3.ne_zero
  have hfNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat fStar hfNonzero
  have hfDegree : (SparsePolyZp.toPoly this._p.toNat fStar).natDegree <
      2 ^ 62 := canonical_natDegree_lt_of_terms_lt this._p.toNat fStar
        hinvariant.fStar_canonical hfNonzero (2 ^ 62)
        hinvariant.fStar_degree_bound
  have hactiveNat := strict_ddf_active_degree_ge this._p.toNat d fStar
    hinvariant.d_bound (lt_of_lt_of_le hinvariant.p0 (Nat.le_refl _))
    hinvariant.fStar_canonical hfNonempty (lt_trans hfDegree (by norm_num))
    hactive
  have hfPositive : 0 <
      (SparsePolyZp.toPoly this._p.toNat fStar).natDegree := by
    have := hinvariant.p0
    omega
  rcases strictPowmodIR_refines this providers.hcfg h this._p.toNat fStar
      providers.mul (providers.mod fStar) hinvariant.h_canonical
      hinvariant.fStar_canonical hfNonempty hinvariant.p3 hfPositive with
    ⟨hPow, hpowRun, hPowCanonical, hPowSemanticRaw⟩
  have hPowSemantic : SparsePolyZp.toPoly this._p.toNat hPow =
      SparsePolyZp.toPoly this._p.toNat h ^ this._p.toNat %ₘ
        SparsePolyZp.toPoly this._p.toNat fStar := by
    exact hPowSemanticRaw
  have hPowDegreeLe :
      (SparsePolyZp.toPoly this._p.toNat hPow).natDegree ≤
        (SparsePolyZp.toPoly this._p.toNat fStar).natDegree := by
    rw [hPowSemantic]
    exact Polynomial.natDegree_modByMonic_le _ hinvariant.p3
  have hPowDegree :
      (SparsePolyZp.toPoly this._p.toNat hPow).natDegree < 2 ^ 62 :=
    lt_of_le_of_lt hPowDegreeLe hfDegree
  have hPowTerms : ∀ term ∈ hPow.toList, term.1.deg < 2 ^ 63 := by
    intro term hterm
    exact lt_trans (lt_of_le_of_lt
      (canonical_term_degree_le_natDegree this._p.toNat hPow
        hPowCanonical term hterm) hPowDegree) (by norm_num)
  have hsubtractCanonical := __upoly_subtract_x_ir_canonical hPow this._p
    rfl hPowCanonical hPowTerms
  have hsubtractSemantic := strict_upoly_subtract_x_refines providers.h2p
    hPow this._p rfl ((canonicalRep_iff_canonical hPow).mpr
      hPowCanonical) hPowTerms
  rcases strictDDFGCDIR_refines this
      (Generated.StrictDDF.__upoly_subtract_x_ir hPow this._p) fStar
      (providers.gcd
        (Generated.StrictDDF.__upoly_subtract_x_ir hPow this._p) fStar)
      hsubtractCanonical
      hinvariant.fStar_canonical hfNonempty hinvariant.p3 with
    ⟨gdRaw, hgcdRun, hgdCanonical, hgdSemanticRaw⟩
  refine ⟨hPow, gdRaw, hpowRun, hgcdRun, hPowCanonical,
    hPowSemantic, hgdCanonical, ?_⟩
  rw [hgdSemanticRaw]
  simpa using congrArg (fun poly => normalize
    (EuclideanDomain.gcd poly
      (SparsePolyZp.toPoly this._p.toNat fStar))) hsubtractSemantic

/-- Execute the concrete suffix of a split DDF iteration.  In particular the
factor is not supplied by a specification: `makeMonic` takes its proved
early-return branch, and quotient/remainder witnesses come from the raw
division loops. -/
theorem strictDDFRawSplitSuffix
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : DDFRawProviders this)
    (original : Polynomial (ZMod this._p.toNat))
    (d : UInt64) (fStar h : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64)) (prime : UInt64)
    (hPow gdRaw : SparsePolyZp)
    (hinvariant : DDFLoopInvariant this original d fStar h result prime)
    (hactive : ¬get_deg fStar < (2 * d).toInt64)
    (hPowCanonical : SparsePolyZp.Canonical this._p.toNat hPow)
    (hPowSemantic : SparsePolyZp.toPoly this._p.toNat hPow =
      SparsePolyZp.toPoly this._p.toNat h ^ this._p.toNat %ₘ
        SparsePolyZp.toPoly this._p.toNat fStar)
    (hgdCanonical : SparsePolyZp.Canonical this._p.toNat gdRaw)
    (hgdSemantic : SparsePolyZp.toPoly this._p.toNat gdRaw = normalize
      (EuclideanDomain.gcd
        (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
        (SparsePolyZp.toPoly this._p.toNat fStar)))
    (hsplit : (!gdRaw.isEmpty && get_deg gdRaw > 0) = true) :
    ∃ quotient hNext,
      strictMakeMonicIR this gdRaw = .ok gdRaw ∧
      strictExactDivIR this fStar gdRaw = .ok quotient ∧
      SparsePolyZp.normalization quotient = quotient ∧
      strictModIR this hPow (SparsePolyZp.normalization quotient)
          ((providers.mod (SparsePolyZp.normalization quotient)).workspace
            hPow) = .ok hNext ∧
      SparsePolyZp.Canonical this._p.toNat quotient ∧
      SparsePolyZp.toPoly this._p.toNat quotient =
        SparsePolyZp.toPoly this._p.toNat fStar /ₘ
          SparsePolyZp.toPoly this._p.toNat gdRaw ∧
      SparsePolyZp.Canonical this._p.toNat hNext ∧
      SparsePolyZp.toPoly this._p.toNat hNext =
        SparsePolyZp.toPoly this._p.toNat hPow %ₘ
          SparsePolyZp.toPoly this._p.toNat quotient ∧
      DDFLoopInvariant this original (d + 1) quotient hNext
        (result.push (gdRaw, d)) prime ∧
      Generated.StrictDDF.ddfRawMeasure quotient (d + 1) <
        Generated.StrictDDF.ddfRawMeasure fStar d := by
  have hfNonzero : SparsePolyZp.toPoly this._p.toNat fStar ≠ 0 :=
    hinvariant.p3.ne_zero
  have hfDegree : (SparsePolyZp.toPoly this._p.toNat fStar).natDegree <
      2 ^ 62 := canonical_natDegree_lt_of_terms_lt this._p.toNat fStar
        hinvariant.fStar_canonical hfNonzero (2 ^ 62)
        hinvariant.fStar_degree_bound
  have hgcdRawNonzero : EuclideanDomain.gcd
      (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
      (SparsePolyZp.toPoly this._p.toNat fStar) ≠ 0 := by
    intro hzero
    have hright := EuclideanDomain.gcd_dvd_right
      (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
      (SparsePolyZp.toPoly this._p.toNat fStar)
    rw [hzero] at hright
    exact hfNonzero (zero_dvd_iff.mp hright)
  have hgdNonzero : SparsePolyZp.toPoly this._p.toNat gdRaw ≠ 0 := by
    rw [hgdSemantic]
    exact mt normalize_eq_zero.mp hgcdRawNonzero
  have hgdMonic : (SparsePolyZp.toPoly this._p.toNat gdRaw).Monic := by
    rw [hgdSemantic]
    exact Polynomial.monic_normalize hgcdRawNonzero
  have hgdDvd : SparsePolyZp.toPoly this._p.toNat gdRaw ∣
      SparsePolyZp.toPoly this._p.toNat fStar := by
    rw [hgdSemantic]
    exact normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right _ _)
  have hgdDegreeLe :
      (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree ≤
        (SparsePolyZp.toPoly this._p.toNat fStar).natDegree :=
    Polynomial.natDegree_le_of_dvd hgdDvd hfNonzero
  have hgdDegree :
      (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree < 2 ^ 63 :=
    lt_of_le_of_lt hgdDegreeLe (lt_trans hfDegree (by norm_num))
  have hsplitParts : gdRaw.isEmpty = false ∧ get_deg gdRaw > 0 := by
    simpa using hsplit
  have hgdPositive :
      0 < (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree :=
    (strict_get_deg_pos_iff_natDegree_pos this._p.toNat gdRaw
      hgdCanonical hgdDegree).mp hsplitParts.2
  have hgdNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat gdRaw hgdNonzero
  have hmonicRun := strictMakeMonicIR_eq_of_monic this gdRaw
    hgdCanonical hgdNonempty hgdMonic
  rcases strictExactDivIR_refines this fStar gdRaw providers.hcfg
      hinvariant.fStar_canonical hgdCanonical hgdNonempty hgdMonic hgdDvd with
    ⟨quotient, hdivRun, hquotientCanonical, hquotientSemantic⟩
  have hnormalization : SparsePolyZp.normalization quotient = quotient :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_normalization_eq_of_canonical
      this._p.toNat quotient hquotientCanonical
  have hfNextMonic :
      (SparsePolyZp.toPoly this._p.toNat quotient).Monic := by
    rw [hquotientSemantic]
    exact Refinement.StrictSquarefreeZp.divByMonic_monic_of_monic_of_dvd
      (SparsePolyZp.toPoly this._p.toNat fStar)
      (SparsePolyZp.toPoly this._p.toNat gdRaw) hinvariant.p3
      hgdMonic hgdDvd
  have hfNextNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat quotient hfNextMonic.ne_zero
  rcases strictModIR_refines_modByMonic this providers.hcfg hPow quotient
      ((providers.mod quotient).workspace hPow) hPowCanonical
      hquotientCanonical hfNextNonempty hfNextMonic with
    ⟨hNext, hmodRun, hNextCanonical, hNextSemantic⟩
  have hmodGenerated : strictModIR this hPow
      (SparsePolyZp.normalization quotient)
      ((providers.mod (SparsePolyZp.normalization quotient)).workspace hPow) =
        .ok hNext := by
    rw [hnormalization]
    exact hmodRun
  have hstep := DDFLoopInvariant.split this original d fStar h result prime
    hPow gdRaw quotient hNext hinvariant hactive hPowCanonical hPowSemantic
    hgdCanonical hgdSemantic hgdPositive hquotientCanonical
    hquotientSemantic hNextCanonical hNextSemantic
  exact ⟨quotient, hNext, hmonicRun, hdivRun, hnormalization,
    hmodGenerated, hquotientCanonical, hquotientSemantic, hNextCanonical,
    hNextSemantic, hstep.1, hstep.2⟩

/-- The generated boolean no-split test is the L2 zero-degree branch for the
concrete normalized GCD result. -/
theorem strictDDFRawNoSplit_degree_zero
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (original : Polynomial (ZMod this._p.toNat))
    (d : UInt64) (fStar h : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64)) (prime : UInt64)
    (hPow gdRaw : SparsePolyZp)
    (hinvariant : DDFLoopInvariant this original d fStar h result prime)
    (hgdCanonical : SparsePolyZp.Canonical this._p.toNat gdRaw)
    (hgdSemantic : SparsePolyZp.toPoly this._p.toNat gdRaw = normalize
      (EuclideanDomain.gcd
        (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
        (SparsePolyZp.toPoly this._p.toNat fStar)))
    (hnoSplit : (!gdRaw.isEmpty && get_deg gdRaw > 0) = false) :
    (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree = 0 := by
  have hfNonzero : SparsePolyZp.toPoly this._p.toNat fStar ≠ 0 :=
    hinvariant.p3.ne_zero
  have hfDegree : (SparsePolyZp.toPoly this._p.toNat fStar).natDegree <
      2 ^ 62 := canonical_natDegree_lt_of_terms_lt this._p.toNat fStar
        hinvariant.fStar_canonical hfNonzero (2 ^ 62)
        hinvariant.fStar_degree_bound
  have hgcdRawNonzero : EuclideanDomain.gcd
      (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
      (SparsePolyZp.toPoly this._p.toNat fStar) ≠ 0 := by
    intro hzero
    have hright := EuclideanDomain.gcd_dvd_right
      (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
      (SparsePolyZp.toPoly this._p.toNat fStar)
    rw [hzero] at hright
    exact hfNonzero (zero_dvd_iff.mp hright)
  have hgdNonzero : SparsePolyZp.toPoly this._p.toNat gdRaw ≠ 0 := by
    rw [hgdSemantic]
    exact mt normalize_eq_zero.mp hgcdRawNonzero
  have hgdDvd : SparsePolyZp.toPoly this._p.toNat gdRaw ∣
      SparsePolyZp.toPoly this._p.toNat fStar := by
    rw [hgdSemantic]
    exact normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right _ _)
  have hgdDegreeLe :
      (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree ≤
        (SparsePolyZp.toPoly this._p.toNat fStar).natDegree :=
    Polynomial.natDegree_le_of_dvd hgdDvd hfNonzero
  have hgdDegree :
      (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree < 2 ^ 63 :=
    lt_of_le_of_lt hgdDegreeLe (lt_trans hfDegree (by norm_num))
  apply Nat.eq_zero_of_not_pos
  intro hpositive
  have hgetPositive :=
    (strict_get_deg_pos_iff_natDegree_pos this._p.toNat gdRaw
      hgdCanonical hgdDegree).mpr hpositive
  have hgdNonempty :=
    Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
      this._p.toNat gdRaw hgdNonzero
  have hnotEmpty : gdRaw.isEmpty = false := by
    simp [Array.isEmpty, Nat.ne_of_gt hgdNonempty]
  simp [hnotEmpty, hgetPositive] at hnoSplit

/-- The unique concrete operation record supplied to the generated DDF shell.
All computational fields execute raw C++ refinements; both proof fields are
discharged from those executions and `DDFLoopInvariant`, not supplied by a
caller. -/
noncomputable def strictDDFRawOps
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : DDFRawProviders this)
    (original : Polynomial (ZMod this._p.toNat)) :
    Generated.StrictDDF.DDFRawOps where
  powmod := fun base e modulus =>
    strictPowmodIR this base e modulus providers.mul (providers.mod modulus)
  gcd := fun left right =>
    strictDDFGCDIR this left right (providers.gcd left right)
  makeMonic := strictMakeMonicIR this
  exactDiv := strictExactDivIR this
  mod := fun dividend divisor =>
    strictModIR this dividend divisor
      ((providers.mod divisor).workspace dividend)
  Invariant := DDFLoopInvariant this original
  splitStep := by
    intro d fStar h result prime hPow gdRaw gd quotient hNext
    intro hinvariant hactive hpowRun hgcdRun hsplit hmonicRun hdivRun hmodRun
    rcases strictDDFRawPrefix this providers original d fStar h result prime
        hinvariant hactive with
      ⟨hPow', gdRaw', hpowRun', hgcdRun', hPowCanonical,
        hPowSemantic, hgdCanonical, hgdSemantic⟩
    have hPowEq : hPow' = hPow :=
      Except.ok.inj (hpowRun'.symm.trans hpowRun)
    subst hPow'
    have hgdEq : gdRaw' = gdRaw :=
      Except.ok.inj (hgcdRun'.symm.trans hgcdRun)
    subst gdRaw'
    have hfNonzero : SparsePolyZp.toPoly this._p.toNat fStar ≠ 0 :=
      hinvariant.p3.ne_zero
    have hfDegree : (SparsePolyZp.toPoly this._p.toNat fStar).natDegree <
        2 ^ 62 := canonical_natDegree_lt_of_terms_lt this._p.toNat fStar
          hinvariant.fStar_canonical hfNonzero (2 ^ 62)
          hinvariant.fStar_degree_bound
    have hgcdRawNonzero : EuclideanDomain.gcd
        (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
        (SparsePolyZp.toPoly this._p.toNat fStar) ≠ 0 := by
      intro hzero
      have hright := EuclideanDomain.gcd_dvd_right
        (SparsePolyZp.toPoly this._p.toNat hPow - Polynomial.X)
        (SparsePolyZp.toPoly this._p.toNat fStar)
      rw [hzero] at hright
      exact hfNonzero (zero_dvd_iff.mp hright)
    have hgdNonzero : SparsePolyZp.toPoly this._p.toNat gdRaw ≠ 0 := by
      rw [hgdSemantic]
      exact mt normalize_eq_zero.mp hgcdRawNonzero
    have hgdMonic : (SparsePolyZp.toPoly this._p.toNat gdRaw).Monic := by
      rw [hgdSemantic]
      exact Polynomial.monic_normalize hgcdRawNonzero
    have hgdDvd : SparsePolyZp.toPoly this._p.toNat gdRaw ∣
        SparsePolyZp.toPoly this._p.toNat fStar := by
      rw [hgdSemantic]
      exact normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right _ _)
    have hgdDegreeLe :
        (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree ≤
          (SparsePolyZp.toPoly this._p.toNat fStar).natDegree :=
      Polynomial.natDegree_le_of_dvd hgdDvd hfNonzero
    have hgdDegree :
        (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree < 2 ^ 63 :=
      lt_of_le_of_lt hgdDegreeLe (lt_trans hfDegree (by norm_num))
    have hsplitParts : gdRaw.isEmpty = false ∧ get_deg gdRaw > 0 := by
      simpa using hsplit
    have hgdPositive :
        0 < (SparsePolyZp.toPoly this._p.toNat gdRaw).natDegree :=
      (strict_get_deg_pos_iff_natDegree_pos this._p.toNat gdRaw
        hgdCanonical hgdDegree).mp hsplitParts.2
    have hgdNonempty :=
      Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
        this._p.toNat gdRaw hgdNonzero
    have hmonicActual := strictMakeMonicIR_eq_of_monic this gdRaw
      hgdCanonical hgdNonempty hgdMonic
    have hgdOutEq : gd = gdRaw :=
      Except.ok.inj (hmonicRun.symm.trans hmonicActual)
    subst gd
    rcases strictExactDivIR_refines this fStar gdRaw providers.hcfg
        hinvariant.fStar_canonical hgdCanonical hgdNonempty hgdMonic hgdDvd with
      ⟨quotient', hdivRun', hquotientCanonical, hquotientSemantic⟩
    have hquotientEq : quotient' = quotient :=
      Except.ok.inj (hdivRun'.symm.trans hdivRun)
    subst quotient'
    have hnormalization : SparsePolyZp.normalization quotient = quotient :=
      Refinement.StrictSquarefreeZp.sparsePolyZp_normalization_eq_of_canonical
        this._p.toNat quotient hquotientCanonical
    rw [hnormalization] at hmodRun
    have hfNextMonic :
        (SparsePolyZp.toPoly this._p.toNat quotient).Monic := by
      rw [hquotientSemantic]
      exact Refinement.StrictSquarefreeZp.divByMonic_monic_of_monic_of_dvd
        (SparsePolyZp.toPoly this._p.toNat fStar)
        (SparsePolyZp.toPoly this._p.toNat gdRaw) hinvariant.p3
        hgdMonic hgdDvd
    have hfNextNonempty :=
      Refinement.StrictSquarefreeZp.sparsePolyZp_size_pos_of_toPoly_ne_zero
        this._p.toNat quotient hfNextMonic.ne_zero
    rcases strictModIR_refines_modByMonic this providers.hcfg hPow quotient
        ((providers.mod quotient).workspace hPow) hPowCanonical
        hquotientCanonical hfNextNonempty hfNextMonic with
      ⟨hNext', hmodRun', hNextCanonical, hNextSemantic⟩
    have hNextEq : hNext' = hNext :=
      Except.ok.inj (hmodRun'.symm.trans hmodRun)
    subst hNext'
    simpa [hnormalization] using
      DDFLoopInvariant.split this original d fStar h result prime hPow
        gdRaw quotient hNext hinvariant hactive hPowCanonical hPowSemantic
        hgdCanonical hgdSemantic hgdPositive hquotientCanonical
        hquotientSemantic hNextCanonical hNextSemantic
  noSplitStep := by
    intro d fStar h result prime hPow gdRaw
    intro hinvariant hactive hpowRun hgcdRun hnoSplit
    rcases strictDDFRawPrefix this providers original d fStar h result prime
        hinvariant hactive with
      ⟨hPow', gdRaw', hpowRun', hgcdRun', hPowCanonical,
        hPowSemantic, hgdCanonical, hgdSemantic⟩
    have hPowEq : hPow' = hPow :=
      Except.ok.inj (hpowRun'.symm.trans hpowRun)
    subst hPow'
    have hgdEq : gdRaw' = gdRaw :=
      Except.ok.inj (hgcdRun'.symm.trans hgcdRun)
    subst gdRaw'
    exact DDFLoopInvariant.noSplit this original d fStar h result prime
      hPow gdRaw hinvariant hactive hPowCanonical hPowSemantic
      hgdCanonical hgdSemantic hnoSplit

/-- The concrete generated DDF entry with every allocation-owning call bound
to its raw implementation and every recursive obligation discharged by the
proved invariant. -/
noncomputable def strictDDFIR
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : DDFRawProviders this)
    (original : Polynomial (ZMod this._p.toNat)) (f : SparsePolyZp)
    (hinitial : ¬f.isEmpty → DDFLoopInvariant this original 1 f
      #[(UMonomial.mk (1 : Int64), Generated.StrictDDF.__make_zp_ir
        (1 : Int64) f[0]!.2.prime)] #[] f[0]!.2.prime) :
    RawExec (Array (SparsePolyZp × UInt64)) :=
  Generated.StrictDDF.__ddf_Zp_raw_ir
    (strictDDFRawOps this providers original) f hinitial

/-- L2 decoding of the generated DDF loop's epilogue.  The raw loop returns
the residual polynomial separately; the public C++ entry appends it exactly
when it has positive degree. -/
noncomputable def ddfRawFinishToL2 (p : Nat) (fStar : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64)) :
    List (Polynomial (ZMod p) × Nat) :=
  if 0 < (SparsePolyZp.toPoly p fStar).natDegree then
    ddfResultToL2 p result ++
      [(SparsePolyZp.toPoly p fStar,
        (SparsePolyZp.toPoly p fStar).natDegree)]
  else
    ddfResultToL2 p result

@[simp] theorem ddfResultToL2_push (p : Nat)
    (result : Array (SparsePolyZp × UInt64))
    (factor : SparsePolyZp) (d : UInt64) :
    ddfResultToL2 p (result.push (factor, d)) =
      ddfResultToL2 p result ++
        [(SparsePolyZp.toPoly p factor, d.toNat)] := by
  simp [ddfResultToL2, Array.toList_push]

/-- Public refinement contract named after the exact cpp2lean-generated C++
entry point.  It proves both representation safety and L2 polynomial
semantics for the generated implementation. -/
theorem __upoly_subtract_x_ir_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (h2p : 2 * this._p.toNat ≤ UInt64.size) (h : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat h)
    (hdegree : ∀ x ∈ h.toList, x.1.deg < 2 ^ 63) :
    SparsePolyZp.Canonical this._p.toNat
        (Generated.StrictDDF.__upoly_subtract_x_ir h this._p) ∧
      SparsePolyZp.toPoly this._p.toNat
          (Generated.StrictDDF.__upoly_subtract_x_ir h this._p) =
        SparsePolyZp.toPoly this._p.toNat h - Polynomial.X := by
  refine ⟨__upoly_subtract_x_ir_canonical h this._p rfl hcanonical hdegree,
    ?_⟩
  exact strict_upoly_subtract_x_refines h2p h this._p rfl
    ((canonicalRep_iff_canonical h).mpr hcanonical) hdegree

/-- Compatibility projection of `__upoly_subtract_x_ir_refines`. -/
theorem strictSubtractXIR_refines
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (h2p : 2 * this._p.toNat ≤ UInt64.size) (h : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat h)
    (hdegree : ∀ x ∈ h.toList, x.1.deg < 2 ^ 63) :
    SparsePolyZp.toPoly this._p.toNat
        (Generated.StrictDDF.__upoly_subtract_x_ir h this._p) =
      SparsePolyZp.toPoly this._p.toNat h - Polynomial.X := by
  exact (__upoly_subtract_x_ir_refines this h2p h hcanonical hdegree).2

/-- Representation half of the generated subtract-X refinement.  This is the
coefficient invariant consumed by the raw sparse-to-dense GCD constructor. -/
theorem strictSubtractXIR_allReduced
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (h : SparsePolyZp)
    (hcanonical : SparsePolyZp.Canonical this._p.toNat h)
    (hdegree : ∀ x ∈ h.toList, x.1.deg < 2 ^ 31) :
    SparsePolyZp.AllReduced this._p.toNat
      (Generated.StrictDDF.__upoly_subtract_x_ir h this._p).toList := by
  have hdegree64 : ∀ x ∈ h.toList, x.1.deg < 2 ^ 63 := by
    intro x hx
    exact lt_trans (hdegree x hx) (by norm_num)
  exact (__upoly_subtract_x_ir_canonical h this._p rfl hcanonical hdegree64).1

end StrictDDF

end Refinement
