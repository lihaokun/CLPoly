import CLPoly.Generated.StrictMul
import CLPoly.Impl.StrictDivremRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictMulRefinement

open Generated.StrictMul
open CLPoly.Impl.StrictWordArithmetic

/-- Mathematical value of the exact raw cells visited by the C++ dot loop.
It shares the loop's reads and failure behavior, but performs unbounded natural
addition so that machine accumulation can be related to it explicitly. -/
def classicalDotNat (heap : RawHeap) (A B : RawPtr UInt64)
    (k stop j : Nat) : RawExec Nat :=
  if h : j ≤ stop then
    match heap.readU64 A j with
    | .error fault => .error fault
    | .ok a =>
      match heap.readU64 B (k - j) with
      | .error fault => .error fault
      | .ok b =>
        match classicalDotNat heap A B k stop (j + 1) with
        | .error fault => .error fault
        | .ok tail => .ok (a.toNat * b.toNat + tail)
  else
    .ok 0
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotLoop_modEq (heap : RawHeap) (A B : RawPtr UInt64)
    (k stop j : Nat) (acc result : Word3) (sum : Nat)
    (hrun : classicalDotLoop heap A B k stop j acc = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    Nat.ModEq (limbBase ^ 3) (word3Value result)
      (word3Value acc + sum) := by
  unfold classicalDotLoop at hrun
  unfold classicalDotNat at hsum
  split at hrun
  next hle =>
    simp only [hle, ↓reduceDIte] at hsum
    cases ha : heap.readU64 A j with
    | error fault => simp [ha] at hrun
    | ok a =>
      simp only [ha] at hrun hsum
      cases hb : heap.readU64 B (k - j) with
      | error fault => simp [hb] at hrun
      | ok b =>
        simp only [hb] at hrun hsum
        cases ht : classicalDotNat heap A B k stop (j + 1) with
        | error fault => simp [ht] at hsum
        | ok tail =>
          simp only [ht] at hsum
          have hsumEq : sum = a.toNat * b.toNat + tail :=
            Except.ok.inj hsum.symm
          let product := Generated.StrictGCD.dense_upoly_zp__umul128_ir 0 0 a b
          let acc' := Generated.StrictGCD.dense_upoly_zp__add_carry3_ir
            acc product.fst product.snd
          have hrec := classicalDotLoop_modEq heap A B k stop (j + 1)
            acc' result tail hrun ht
          have hstep : Nat.ModEq (limbBase ^ 3) (word3Value acc')
              (word3Value acc + a.toNat * b.toNat) := by
            simpa [product, acc'] using addMulWord3_modEq acc a b
          have htotal := hrec.trans (hstep.add_right tail)
          simpa [hsumEq, Nat.add_assoc] using htotal
  next hnot =>
    simp only [hnot, ↓reduceDIte] at hsum
    have hresult : result = acc := Except.ok.inj hrun.symm
    have hzero : sum = 0 := Except.ok.inj hsum.symm
    subst result
    subst sum
    simpa using (Nat.ModEq.refl (word3Value acc) :
      Nat.ModEq (limbBase ^ 3) (word3Value acc) (word3Value acc))
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotNat_ok (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB) :
    ∃ sum, classicalDotNat heap A B k stop j = .ok sum := by
  unfold classicalDotNat
  split
  next hle =>
    rcases heap.readU64_of_valid A lenA j hA
      (hAIndex j (Nat.le_refl _) hle) with ⟨a, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid B lenB (k - j) hB
      (hBIndex j (Nat.le_refl _) hle) with ⟨b, hb⟩
    simp only [hb]
    rcases classicalDotNat_ok heap A B lenA lenB k stop (j + 1)
      hA hB
      (by intro t hjt hts; exact hAIndex t (by omega) hts)
      (by intro t hjt hts; exact hBIndex t (by omega) hts) with
      ⟨tail, htail⟩
    rw [htail]
    exact ⟨a.toNat * b.toNat + tail, rfl⟩
  next => exact ⟨0, rfl⟩
termination_by stop + 1 - j
decreasing_by omega

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

theorem classicalDotLoop_raw_sum (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (acc : Word3)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB) :
    ∃ result sum,
      classicalDotLoop heap A B k stop j acc = .ok result ∧
      classicalDotNat heap A B k stop j = .ok sum ∧
      Nat.ModEq (limbBase ^ 3) (word3Value result)
        (word3Value acc + sum) := by
  rcases classicalDotLoop_ok heap A B lenA lenB k stop j acc hA hB
    hAIndex hBIndex with ⟨result, hrun⟩
  rcases classicalDotNat_ok heap A B lenA lenB k stop j hA hB
    hAIndex hBIndex with ⟨sum, hsum⟩
  exact ⟨result, sum, hrun, hsum,
    classicalDotLoop_modEq heap A B k stop j acc result sum hrun hsum⟩

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
