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

set_option autoImplicit false

namespace Refinement

namespace StrictDDF

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

end StrictDDF

end Refinement
