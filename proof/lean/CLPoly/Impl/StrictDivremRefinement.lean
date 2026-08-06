import CLPoly.Generated.StrictDivrem

namespace CLPoly.Impl.StrictDivremRefinement

open Generated.StrictDivrem

/-- Natural-language proof outline:

At an iteration with `i < lenA`, validity of `A[0..lenA)` gives a successful
read of `A[i]`; validity of the `lenA`-element W3 slice gives a successful
three-limb write at `W3[i]`.  `writeWord3_preserves_valid` shows that this
write changes no allocation sizes, so both slice invariants hold for the
recursive heap.  The recursive measure `lenA - (i+1)` is smaller than
`lenA-i`.  When `i ≥ lenA`, the source loop exits without a heap access. -/
theorem initW3Loop_ok (heap : RawHeap) (A : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenA i : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hW3 : heap.ValidWord3Slice W3 lenA) (hi : i ≤ lenA) :
    ∃ heap', initW3Loop heap A W3 lenA i = .ok heap' ∧
      heap'.ValidU64Slice A lenA ∧ heap'.ValidWord3Slice W3 lenA := by
  rw [initW3Loop]
  split
  next hlt =>
    rcases heap.readU64_of_valid A lenA i hA hlt with ⟨value, hread⟩
    simp only [hread]
    rcases heap.writeWord3_of_valid W3 lenA i
      { lo := value, mid := 0, hi := 0 } hW3 hlt with ⟨heap1, hwrite⟩
    simp only [hwrite]
    have hA1 : heap1.ValidU64Slice A lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i
        { lo := value, mid := 0, hi := 0 } hwrite A lenA).mp hA
    have hW31 : heap1.ValidWord3Slice W3 lenA :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 i
        { lo := value, mid := 0, hi := 0 } hwrite
        (RawPtr.reinterpret W3) (3 * lenA)).mp hW3
    exact initW3Loop_ok heap1 A W3 lenA (i + 1) hA1 hW31 (by omega)
  next hnot =>
    exact ⟨heap, rfl, hA, hW3⟩
termination_by lenA - i
decreasing_by omega

end CLPoly.Impl.StrictDivremRefinement
