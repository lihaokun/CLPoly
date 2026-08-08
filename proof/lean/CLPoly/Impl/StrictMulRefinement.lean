import CLPoly.Generated.StrictMul
import CLPoly.Impl.StrictDivremRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictMulRefinement

open Generated.StrictMul

theorem classicalDotLoop_ok (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (acc : Word3)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB) :
    ∃ result, classicalDotLoop heap A B k stop j acc = .ok result := by
  unfold classicalDotLoop
  split
  next hle =>
    have hjA := hAIndex j (Nat.le_refl _) hle
    have hjB := hBIndex j (Nat.le_refl _) hle
    rcases heap.readU64_of_valid A lenA j hA hjA with ⟨a, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid B lenB (k - j) hB hjB with ⟨b, hb⟩
    simp only [hb]
    let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 a b
    let acc' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
      acc product.fst product.snd
    apply classicalDotLoop_ok heap A B lenA lenB k stop (j + 1) acc'
      hA hB
    · intro t hjt hts
      exact hAIndex t (by omega) hts
    · intro t hjt hts
      exact hBIndex t (by omega) hts
  next => exact ⟨acc, rfl⟩
termination_by stop + 1 - j
decreasing_by omega

theorem classical_index_bounds (lenA lenB k j : Nat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hk : k < lenA + lenB - 1)
    (hjMin : (if k ≥ lenB then k - lenB + 1 else 0) ≤ j)
    (hjMax : j ≤ (if k < lenA then k else lenA - 1)) :
    j < lenA ∧ k - j < lenB := by
  split at hjMax
  next hkA =>
    constructor
    · omega
    · split at hjMin <;> omega
  next hkNotA =>
    constructor
    · omega
    · split at hjMin <;> omega

theorem classicalOuterLoop_ok (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (lenA lenB lenC k : Nat) (heap : RawHeap)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hlenC : lenC = lenA + lenB - 1)
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB) :
    ∃ heap', classicalOuterLoop this C A B lenA lenB lenC k heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' := by
  unfold classicalOuterLoop
  split
  next hk =>
    let jMin := if k ≥ lenB then k - lenB + 1 else 0
    let jMax := if k < lenA then k else lenA - 1
    have hrange : jMin ≤ jMax := by
      dsimp [jMin, jMax]
      split <;> split <;> omega
    rcases classicalDotLoop_ok heap A B lenA lenB k jMax jMin
      { lo := 0, mid := 0, hi := 0 } hA hB
      (by
        intro t hjt hts
        exact (classical_index_bounds lenA lenB k t hApos hBpos
          (by omega) hjt hts).1)
      (by
        intro t hjt hts
        exact (classical_index_bounds lenA lenB k t hApos hBpos
          (by omega) hjt hts).2) with ⟨acc, hdot⟩
    simp only [jMin, jMax, hdot]
    let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      acc.hi acc.mid acc.lo this._p this._ninv this._norm
    rcases heap.writeU64_of_valid C lenC k value hC hk with ⟨heap1, hw⟩
    simp only [value, hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C k value hw
    rcases classicalOuterLoop_ok this C A B lenA lenB lenC (k + 1) heap1
      hApos hBpos hlenC ((hlayout1 C lenC).mp hC)
      ((hlayout1 A lenA).mp hA) ((hlayout1 B lenB).mp hB) with
      ⟨heap2, hrun, hlayout2⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next => exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by lenC - k
decreasing_by omega

theorem classicalMul_ok (this : DenseUPolyZp)
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hC : heap.ValidU64Slice C (lenA + lenB - 1))
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB) :
    ∃ heap', dense_upoly_zp__classical_mul_ir this C A lenA B lenB heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' := by
  rcases classicalOuterLoop_ok this C A B lenA lenB
    (lenA + lenB - 1) 0 heap hApos hBpos rfl hC hA hB with
    ⟨heap', hrun, hlayout⟩
  refine ⟨heap', ?_, hlayout⟩
  simp [dense_upoly_zp__classical_mul_ir, Nat.ne_of_gt hApos,
    Nat.ne_of_gt hBpos, hrun]

end CLPoly.Impl.StrictMulRefinement
