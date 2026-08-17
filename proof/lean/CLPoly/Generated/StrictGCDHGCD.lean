import CLPoly.Generated.StrictHGCDChecked
import CLPoly.Generated.StrictEuclidGCD

set_option autoImplicit false

namespace Generated.StrictGCDHGCD

open Generated.StrictDivrem
open Generated.StrictEuclidGCD
open Generated.StrictHGCD

/-- Observable return of `_gcd_hgcd`; `lenG` models the source reference
output and the polynomial itself remains in the caller-provided `G` region. -/
structure GcdHgcdRawResult where
  heap : RawHeap
  lenG : Nat

/-- Proof-only termination contract for the HGCD main loop.  It refers to the
actual successful checked recursive execution and cannot select or replace
its result. -/
def HgcdGcdLoopLengthDecreases (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid) (W scratch : RawPtr UInt64) : Prop :=
  ∀ (heap : RawHeap) (G J R : RawPtr UInt64)
    (lenG lenJ lenR : Nat) (result : HgcdRecursiveResult),
    0 < lenR → lenR < lenJ →
    hgcdRecursiveCallChecked this M hM false G J J R lenJ lenR W scratch
        heap = .ok result →
    result.lenB < lenJ

/-- Exact raw lowering of the `while (len_j != 0)` block in `_gcd_hgcd`.
The Euclid helper's five local vectors are exposed as fixed raw allocations.
The recursion measure is the source variable `len_j`. -/
def gcdHgcdLoop (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid)
    (G J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (W scratch : RawPtr UInt64)
    (euclidA euclidB euclidQ euclidR : RawPtr UInt64)
    (euclidW3 : RawPtr Word3)
    (euclidDecrease : DivremLengthDecreases this euclidQ euclidW3)
    (hdecrease : HgcdGcdLoopLengthDecreases this M hM W scratch) :
    (heap : RawHeap) → (lenG lenJ : Nat) → RawExec GcdHgcdRawResult
  | heap, lenG, 0 => .ok ⟨heap, lenG⟩
  | heap, lenG, lenJ + 1 =>
      match hdiv : dense_upoly_zp__poly_divrem_ir this Q R G lenG J
          (lenJ + 1) W3 heap with
      | .error fault => .error fault
      | .ok (heap1, lenQ, lenR) =>
        if lenR = 0 then
          match heap1.copyU64 G J (lenJ + 1) with
          | .error fault => .error fault
          | .ok heap2 => .ok ⟨heap2, lenJ + 1⟩
        else if lenJ + 1 < hgcdRecursiveCutoff then
          match dense_upoly_zp__gcd_euclid_raw_ir this G J R euclidA euclidB
              euclidQ euclidR euclidW3 (lenJ + 1) lenR euclidDecrease heap1 with
          | .error fault => .error fault
          | .ok euclid => .ok ⟨euclid.heap, euclid.lenG⟩
        else
          match hcall : hgcdRecursiveCallChecked this M hM false G J J R
              (lenJ + 1) lenR W scratch heap1 with
          | .error fault => .error fault
          | .ok result =>
            gcdHgcdLoop this M hM G J Q R W3 W scratch euclidA euclidB
              euclidQ euclidR euclidW3 euclidDecrease hdecrease result.heap
              result.lenA result.lenB
termination_by heap lenG lenJ => lenJ
decreasing_by
  apply hdecrease heap1 G J R result.lenA (lenJ + 1) lenR result
  · omega
  · exact polyDivrem_remainder_lt this Q R G lenG J (lenJ + 1) W3 heap
      heap1 lenQ lenR hdiv
  · exact hcall

/-- Exact raw execution of C++ `_gcd_hgcd` after exposing its workspace and
the nested Euclid helper's local vector allocations. -/
def dense_upoly_zp__gcd_hgcd_raw_ir (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid)
    (G A B J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (W scratch : RawPtr UInt64)
    (euclidA euclidB euclidQ euclidR : RawPtr UInt64)
    (euclidW3 : RawPtr Word3)
    (lenA lenB : Nat)
    (euclidDecrease : DivremLengthDecreases this euclidQ euclidW3)
    (loopDecrease : HgcdGcdLoopLengthDecreases this M hM W scratch)
    (heap : RawHeap) : RawExec GcdHgcdRawResult :=
  match dense_upoly_zp__poly_divrem_ir this Q R A lenA B lenB W3 heap with
  | .error fault => .error fault
  | .ok (heap1, _, lenR) =>
    if lenR = 0 then
      match heap1.copyU64 G B lenB with
      | .error fault => .error fault
      | .ok heap2 => .ok ⟨heap2, lenB⟩
    else
      match hgcdRecursiveCallChecked this M hM false G J B R lenB lenR W
          scratch heap1 with
      | .error fault => .error fault
      | .ok first =>
        gcdHgcdLoop this M hM G J Q R W3 W scratch euclidA euclidB euclidQ
          euclidR euclidW3 euclidDecrease loopDecrease first.heap first.lenA
          first.lenB

end Generated.StrictGCDHGCD
