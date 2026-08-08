import CLPoly.Generated.StrictPolyAddSub
import CLPoly.Impl.StrictDivremRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictPolyAddSubRefinement

open Generated.StrictPolyAddSub
open CLPoly.Impl.StrictDivremRefinement

theorem addCommonLoop_ok (this : DenseUPolyZp) (C A B : RawPtr UInt64)
    (limit i : Nat) (heap : RawHeap)
    (hC : heap.ValidU64Slice C limit)
    (hA : heap.ValidU64Slice A limit)
    (hB : heap.ValidU64Slice B limit) (hi : i ≤ limit) :
    ∃ heap', addCommonLoop this C A B limit i heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  unfold addCommonLoop
  split
  next hlt =>
    rcases heap.readU64_of_valid A limit i hA hlt with ⟨a, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid B limit i hB hlt with ⟨b, hb⟩
    simp only [hb]
    rcases heap.writeU64_of_valid C limit i
      (dense_upoly_zp_nmod_add_ir this a b) hC hlt with ⟨heap1, hw⟩
    simp only [hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C i _ hw
    have hC1 := (hlayout1 C limit).mp hC
    have hA1 := (hlayout1 A limit).mp hA
    have hB1 := (hlayout1 B limit).mp hB
    rcases addCommonLoop_ok this C A B limit (i + 1) heap1
      hC1 hA1 hB1 (by omega) with ⟨heap2, hrun, hlayout2⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next hnot =>
    exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by limit - i
decreasing_by omega

theorem subCommonLoop_ok (this : DenseUPolyZp) (C A B : RawPtr UInt64)
    (limit i : Nat) (heap : RawHeap)
    (hC : heap.ValidU64Slice C limit)
    (hA : heap.ValidU64Slice A limit)
    (hB : heap.ValidU64Slice B limit) (hi : i ≤ limit) :
    ∃ heap', subCommonLoop this C A B limit i heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  unfold subCommonLoop
  split
  next hlt =>
    rcases heap.readU64_of_valid A limit i hA hlt with ⟨a, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid B limit i hB hlt with ⟨b, hb⟩
    simp only [hb]
    rcases heap.writeU64_of_valid C limit i
      (dense_upoly_zp_nmod_sub_ir this a b) hC hlt with ⟨heap1, hw⟩
    simp only [hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C i _ hw
    rcases subCommonLoop_ok this C A B limit (i + 1) heap1
      ((hlayout1 C limit).mp hC) ((hlayout1 A limit).mp hA)
      ((hlayout1 B limit).mp hB) (by omega) with
      ⟨heap2, hrun, hlayout2⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next hnot =>
    exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by limit - i
decreasing_by omega

theorem subNegTailLoop_ok (this : DenseUPolyZp) (C B : RawPtr UInt64)
    (limit i : Nat) (heap : RawHeap)
    (hC : heap.ValidU64Slice C limit)
    (hB : heap.ValidU64Slice B limit) (hi : i ≤ limit) :
    ∃ heap', subNegTailLoop this C B limit i heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  unfold subNegTailLoop
  split
  next hlt =>
    rcases heap.readU64_of_valid B limit i hB hlt with ⟨b, hb⟩
    simp only [hb]
    let value := if b = 0 then 0 else this._p - b
    rcases heap.writeU64_of_valid C limit i value hC hlt with ⟨heap1, hw⟩
    rw [show (if b = 0 then 0 else this._p - b) = value by rfl, hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C i value hw
    rcases subNegTailLoop_ok this C B limit (i + 1) heap1
      ((hlayout1 C limit).mp hC) ((hlayout1 B limit).mp hB) (by omega) with
      ⟨heap2, hrun, hlayout2⟩
    simp only
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next hnot =>
    exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by limit - i
decreasing_by omega

theorem polyAdd_ok (this : DenseUPolyZp) (C A : RawPtr UInt64)
    (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat) (heap : RawHeap)
    (hC : heap.ValidU64Slice C (max lenA lenB))
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB) :
    ∃ heap' length,
      dense_upoly_zp__poly_add_ir this C A lenA B lenB heap =
        .ok (heap', length) ∧
      RawHeap.SameLayout heap heap' ∧ length ≤ max lenA lenB := by
  let minLen := min lenA lenB
  let maxLen := max lenA lenB
  have hMinA : minLen ≤ lenA := by simp [minLen]
  have hMinB : minLen ≤ lenB := by simp [minLen]
  have hMinMax : minLen ≤ maxLen := by simp [minLen, maxLen]
  have hCMin := heap.validU64Slice_mono C maxLen minLen hC hMinMax
  have hAMin := heap.validU64Slice_mono A lenA minLen hA hMinA
  have hBMin := heap.validU64Slice_mono B lenB minLen hB hMinB
  rcases addCommonLoop_ok this C A B minLen 0 heap
    hCMin hAMin hBMin (by omega) with ⟨heap1, hloop, hlayout1⟩
  have hC1 : heap1.ValidU64Slice C maxLen := (hlayout1 C maxLen).mp hC
  have hA1 : heap1.ValidU64Slice A lenA := (hlayout1 A lenA).mp hA
  have hB1 : heap1.ValidU64Slice B lenB := (hlayout1 B lenB).mp hB
  have htail : ∃ heap2,
      (if lenA > lenB then
          if C.sameAddress A then .ok heap1
          else heap1.copyU64 (C.add minLen) (A.add minLen) (lenA - minLen)
        else if lenB > lenA then
          if C.sameAddress B then .ok heap1
          else heap1.copyU64 (C.add minLen) (B.add minLen) (lenB - minLen)
        else .ok heap1) = .ok heap2 ∧ RawHeap.SameLayout heap1 heap2 := by
    split
    next hlongA =>
      split
      next => exact ⟨heap1, rfl, fun _ _ => Iff.rfl⟩
      next =>
        have hDst := heap1.validU64Slice_add C maxLen minLen
          (lenA - minLen) hC1 (by simp [maxLen, minLen])
        have hSrc := heap1.validU64Slice_add A lenA minLen
          (lenA - minLen) hA1 (by omega)
        exact copyU64_ok heap1 (C.add minLen) (A.add minLen)
          (lenA - minLen) hDst hSrc
    next hnotA =>
      split
      next hlongB =>
        split
        next => exact ⟨heap1, rfl, fun _ _ => Iff.rfl⟩
        next =>
          have hDst := heap1.validU64Slice_add C maxLen minLen
            (lenB - minLen) hC1 (by simp [maxLen, minLen])
          have hSrc := heap1.validU64Slice_add B lenB minLen
            (lenB - minLen) hB1 (by omega)
          exact copyU64_ok heap1 (C.add minLen) (B.add minLen)
            (lenB - minLen) hDst hSrc
      next => exact ⟨heap1, rfl, fun _ _ => Iff.rfl⟩
  rcases htail with ⟨heap2, htailRun, hlayout2⟩
  have hC2 : heap2.ValidU64Slice C maxLen := (hlayout2 C maxLen).mp hC1
  rcases normaliseU64_ok heap2 C maxLen hC2 with
    ⟨length, hnorm, hlength⟩
  refine ⟨heap2, length, ?_, ?_, ?_⟩
  · simp only [dense_upoly_zp__poly_add_ir, minLen, maxLen, hloop,
      htailRun, hnorm]
  · exact fun ptr count =>
      (hlayout1 ptr count).trans (hlayout2 ptr count)
  · simpa [maxLen] using hlength

theorem polySub_ok (this : DenseUPolyZp) (C A : RawPtr UInt64)
    (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat) (heap : RawHeap)
    (hC : heap.ValidU64Slice C (max lenA lenB))
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB) :
    ∃ heap' length,
      dense_upoly_zp__poly_sub_ir this C A lenA B lenB heap =
        .ok (heap', length) ∧
      RawHeap.SameLayout heap heap' ∧ length ≤ max lenA lenB := by
  let minLen := min lenA lenB
  let maxLen := max lenA lenB
  have hMinA : minLen ≤ lenA := by simp [minLen]
  have hMinB : minLen ≤ lenB := by simp [minLen]
  have hMinMax : minLen ≤ maxLen := by simp [minLen, maxLen]
  have hCMin := heap.validU64Slice_mono C maxLen minLen hC hMinMax
  have hAMin := heap.validU64Slice_mono A lenA minLen hA hMinA
  have hBMin := heap.validU64Slice_mono B lenB minLen hB hMinB
  rcases subCommonLoop_ok this C A B minLen 0 heap
    hCMin hAMin hBMin (by omega) with ⟨heap1, hloop, hlayout1⟩
  have hC1 : heap1.ValidU64Slice C maxLen := (hlayout1 C maxLen).mp hC
  have hA1 : heap1.ValidU64Slice A lenA := (hlayout1 A lenA).mp hA
  have hB1 : heap1.ValidU64Slice B lenB := (hlayout1 B lenB).mp hB
  have htail : ∃ heap2,
      (if lenA > lenB then
          if C.sameAddress A then .ok heap1
          else heap1.copyU64 (C.add minLen) (A.add minLen) (lenA - minLen)
        else if lenB > lenA then
          subNegTailLoop this C B lenB minLen heap1
        else .ok heap1) = .ok heap2 ∧ RawHeap.SameLayout heap1 heap2 := by
    split
    next hlongA =>
      split
      next => exact ⟨heap1, rfl, fun _ _ => Iff.rfl⟩
      next =>
        have hDst := heap1.validU64Slice_add C maxLen minLen
          (lenA - minLen) hC1 (by simp [maxLen, minLen])
        have hSrc := heap1.validU64Slice_add A lenA minLen
          (lenA - minLen) hA1 (by omega)
        exact copyU64_ok heap1 (C.add minLen) (A.add minLen)
          (lenA - minLen) hDst hSrc
    next hnotA =>
      split
      next hlongB =>
        exact subNegTailLoop_ok this C B lenB minLen heap1
          (heap1.validU64Slice_mono C maxLen lenB hC1 (by
            simp [maxLen])) hB1 (by simp [minLen])
      next => exact ⟨heap1, rfl, fun _ _ => Iff.rfl⟩
  rcases htail with ⟨heap2, htailRun, hlayout2⟩
  have hC2 : heap2.ValidU64Slice C maxLen := (hlayout2 C maxLen).mp hC1
  rcases normaliseU64_ok heap2 C maxLen hC2 with
    ⟨length, hnorm, hlength⟩
  refine ⟨heap2, length, ?_, ?_, ?_⟩
  · simp only [dense_upoly_zp__poly_sub_ir, minLen, maxLen, hloop,
      htailRun, hnorm]
  · exact fun ptr count =>
      (hlayout1 ptr count).trans (hlayout2 ptr count)
  · simpa [maxLen] using hlength

end CLPoly.Impl.StrictPolyAddSubRefinement
