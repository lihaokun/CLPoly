import CLPoly.Generated.StrictDivrem

namespace Generated.StrictEuclidGCD

open Generated.StrictDivrem

/-- Proof-only termination contract for the generated division call used by
the source Euclid loop. -/
def DivremLengthDecreases (this : DenseUPolyZp) (Q : RawPtr UInt64)
    (W3 : RawPtr Word3) : Prop :=
  ∀ (heap : RawHeap) (A B R : RawPtr UInt64) (lenA lenB : Nat)
    (heap' : RawHeap) (lenQ lenR : Nat),
    0 < lenB →
    dense_upoly_zp__poly_divrem_ir this Q R A lenA B lenB W3 heap =
      .ok (heap', lenQ, lenR) →
    lenR < lenB

/-- Raw specialization of the C++ Euclid loop.  `Q` is the persistent local
quotient object.  The source moves `b` into `a` and `r` into `b`; therefore
the old `a` allocation becomes the next remainder scratch buffer. -/
def euclidLoop (this : DenseUPolyZp) (Q : RawPtr UInt64)
    (W3 : RawPtr Word3) (hdec : DivremLengthDecreases this Q W3) :
    (heap : RawHeap) → (A : RawPtr UInt64) → (lenA : Nat) →
      (B : RawPtr UInt64) → (lenB : Nat) → (R : RawPtr UInt64) →
      RawExec (RawHeap × RawPtr UInt64 × Nat)
  | heap, A, lenA, _, 0, _ => .ok (heap, A, lenA)
  | heap, A, lenA, B, lenB + 1, R =>
      match hrun : dense_upoly_zp__poly_divrem_ir this Q R A lenA B
          (lenB + 1) W3 heap with
      | .error fault => .error fault
      | .ok (heap', _, lenR) =>
          euclidLoop this Q W3 hdec heap' B (lenB + 1) R lenR A
termination_by heap A lenA B lenB R => lenB
decreasing_by
  exact hdec heap A B R lenA (lenB + 1) _ _ _ (by omega) hrun

/-- Observable result of the C++ `_gcd_euclid` raw helper.  The output pointer
is fixed by the call; only the final heap and logical output length are
returned in place of the source `size_t& len_G`. -/
structure GcdEuclidRawResult where
  heap : RawHeap
  lenG : Nat

/-- Exact raw execution order of `dense_upoly_zp::_gcd_euclid` after exposing
the five vector allocations as caller-provided regions: copy `A`, copy `B`,
run the generated Euclid rotation, then copy the last live dividend to `G`.
All allocation validity and separation remain semantic preconditions rather
than default-valued array operations. -/
def dense_upoly_zp__gcd_euclid_raw_ir (this : DenseUPolyZp)
    (G A B aBuf bBuf Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenA lenB : Nat) (hdec : DivremLengthDecreases this Q W3)
    (heap : RawHeap) : RawExec GcdEuclidRawResult :=
  match heap.copyU64 aBuf A lenA with
  | .error fault => .error fault
  | .ok heap1 =>
    match heap1.copyU64 bBuf B lenB with
    | .error fault => .error fault
    | .ok heap2 =>
      match euclidLoop this Q W3 hdec heap2 aBuf lenA bBuf lenB R with
      | .error fault => .error fault
      | .ok (heap3, out, outLen) =>
        match heap3.copyU64 G out outLen with
        | .error fault => .error fault
        | .ok heap4 => .ok ⟨heap4, outLen⟩

end Generated.StrictEuclidGCD
