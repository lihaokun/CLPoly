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

set_option autoImplicit false

namespace Refinement

namespace StrictDDF

open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictPolynomialGCDRefinement
open Generated.StrictDivrem

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

end StrictDDF

end Refinement
