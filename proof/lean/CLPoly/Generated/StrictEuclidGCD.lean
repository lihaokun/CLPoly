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

end Generated.StrictEuclidGCD
