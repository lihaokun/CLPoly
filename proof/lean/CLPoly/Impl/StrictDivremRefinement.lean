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

/-- Natural-language proof outline:

For `j ≤ d`, the divisor invariant makes `B[j]` readable.  The bound
`i+d < lenW3` implies `i+j < lenW3`, so the accumulator is readable and
writable.  The multiplication and carry routines are total word operations.
The W3 write preserves both B and W3 allocation invariants; recurse at `j+1`
with the strictly smaller measure `d+1-j`.  At `j>d`, the source loop exits. -/
theorem addMulLoop_ok (heap : RawHeap) (B : RawPtr UInt64)
    (W3 : RawPtr Word3) (lenW3 i d j : Nat) (c : UInt64)
    (other : RawPtr UInt64) (otherLen : Nat)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hOther : heap.ValidU64Slice other otherLen)
    (htop : i + d < lenW3) (hj : j ≤ d + 1) :
    ∃ heap', addMulLoop heap B W3 i d j c = .ok heap' ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 ∧
      heap'.ValidU64Slice other otherLen := by
  rw [addMulLoop]
  split
  next hle =>
    have hjB : j < d + 1 := by omega
    have hijW : i + j < lenW3 := by omega
    rcases heap.readU64_of_valid B (d + 1) j hB hjB with ⟨bj, hreadB⟩
    simp only [hreadB]
    rcases heap.readWord3_of_valid W3 lenW3 (i + j) hW3 hijW with
      ⟨accum, hreadW⟩
    simp only [hreadW]
    let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 c bj
    let accum' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
      accum product.fst product.snd
    rcases heap.writeWord3_of_valid W3 lenW3 (i + j) accum' hW3 hijW with
      ⟨heap1, hwrite⟩
    dsimp [product, accum'] at hwrite ⊢
    simp only [hwrite]
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    have hOther1 : heap1.ValidU64Slice other otherLen :=
      (RawHeap.writeWord3_preserves_valid heap heap1 W3 (i + j) accum'
        hwrite other otherLen).mp hOther
    exact addMulLoop_ok heap1 B W3 lenW3 i d (j + 1) c other otherLen
      hB1 hW31 hOther1 htop (by omega)
  next hnot =>
    exact ⟨heap, rfl, hB, hW3, hOther⟩
termination_by d + 1 - j
decreasing_by omega

/-- The descending quotient loop cannot fault when Q, B and W3 have the
capacities documented by the C++ raw API.  Its induction variable is exactly
the source `ii`: the successor case processes coefficient `i = ii-1`, then
recurses on the predecessor. -/
theorem quotientLoop_ok (this : DenseUPolyZp) (Q B : RawPtr UInt64)
    (W3 : RawPtr Word3) (qLen d lenW3 : Nat) (invLc : UInt64)
    (heap : RawHeap) (ii : Nat)
    (hQ : heap.ValidU64Slice Q qLen)
    (hB : heap.ValidU64Slice B (d + 1))
    (hW3 : heap.ValidWord3Slice W3 lenW3)
    (hii : ii ≤ qLen) (hspan : qLen + d ≤ lenW3) :
    ∃ heap', quotientLoop this Q B W3 d invLc heap ii = .ok heap' ∧
      heap'.ValidU64Slice Q qLen ∧
      heap'.ValidU64Slice B (d + 1) ∧
      heap'.ValidWord3Slice W3 lenW3 := by
  cases ii with
  | zero => exact ⟨heap, rfl, hQ, hB, hW3⟩
  | succ i =>
    have hiQ : i < qLen := by omega
    have hiW : i + d < lenW3 := by omega
    simp only [quotientLoop]
    rcases heap.readWord3_of_valid W3 lenW3 (i + d) hW3 hiW with
      ⟨accum, hread⟩
    simp only [hread]
    let r := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      accum.hi accum.mid accum.lo this._p this._ninv this._norm
    let qi := Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this r invLc
    rcases heap.writeU64_of_valid Q qLen i qi hQ hiQ with ⟨heap1, hwrite⟩
    dsimp [r, qi] at hwrite ⊢
    simp only [hwrite]
    have hQ1 : heap1.ValidU64Slice Q qLen :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite Q qLen).mp hQ
    have hB1 : heap1.ValidU64Slice B (d + 1) :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite B (d + 1)).mp hB
    have hW31 : heap1.ValidWord3Slice W3 lenW3 :=
      (RawHeap.writeU64_preserves_valid heap heap1 Q i _ hwrite
        (RawPtr.reinterpret W3) (3 * lenW3)).mp hW3
    split
    next hnonzero =>
      rcases addMulLoop_ok heap1 B W3 lenW3 i d 0 (this._p -
          Generated.StrictGCD.dense_upoly_zp_nmod_mul_ir this
            (Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
              accum.hi accum.mid accum.lo this._p this._ninv this._norm) invLc)
          Q qLen hB1 hW31 hQ1 hiW (by omega) with
        ⟨heap2, hadd, hB2, hW32, hQ2⟩
      simp only [hadd]
      exact quotientLoop_ok this Q B W3 qLen d lenW3 invLc heap2 i
        hQ2 hB2 hW32 (by omega) hspan
    next hzero =>
      exact quotientLoop_ok this Q B W3 qLen d lenW3 invLc heap1 i
        hQ1 hB1 hW31 (by omega) hspan

end CLPoly.Impl.StrictDivremRefinement
