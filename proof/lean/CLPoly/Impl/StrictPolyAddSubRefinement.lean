import CLPoly.Generated.StrictPolyAddSub
import CLPoly.Impl.StrictDivremRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictPolyAddSubRefinement

open Generated.StrictPolyAddSub
open CLPoly.Impl.StrictDivremRefinement

theorem nmodAdd_toNat (this : DenseUPolyZp) (a b : UInt64)
    (hp : this._p ≠ 0) (ha : a < this._p) (hb : b < this._p) :
    (dense_upoly_zp_nmod_add_ir this a b).toNat =
      (a.toNat + b.toNat) % this._p.toNat := by
  have hpNat : 0 < this._p.toNat := by
    exact Nat.pos_of_ne_zero (fun h => hp (UInt64.toNat_inj.mp (by simpa using h)))
  have haNat : a.toNat < this._p.toNat := by
    simpa [UInt64.lt_iff_toNat_lt] using ha
  have hbNat : b.toNat < this._p.toNat := by
    simpa [UInt64.lt_iff_toNat_lt] using hb
  have haLe : a ≤ this._p := by
    simpa [UInt64.le_iff_toNat_le] using Nat.le_of_lt haNat
  have hsub : (this._p - a).toNat = this._p.toNat - a.toNat :=
    UInt64.toNat_sub_of_le _ _ haLe
  simp only [dense_upoly_zp_nmod_add_ir]
  split
  next hlt =>
    have hltNat : b.toNat < this._p.toNat - a.toNat := by
      simpa [UInt64.lt_iff_toNat_lt, hsub] using hlt
    have hsumP : a.toNat + b.toNat < this._p.toNat := by omega
    have hsum64 : a.toNat + b.toNat < UInt64.size :=
      lt_trans hsumP (UInt64.toNat_lt this._p)
    rw [UInt64.toNat_add, Nat.mod_eq_of_lt hsum64,
      Nat.mod_eq_of_lt hsumP]
  next hnot =>
    have hleNat : this._p.toNat - a.toNat ≤ b.toNat := by
      have : ¬b.toNat < (this._p - a).toNat := by
        simpa [UInt64.lt_iff_toNat_lt] using hnot
      omega
    have hleWord : this._p - a ≤ b := by
      simpa [UInt64.le_iff_toNat_le, hsub] using hleNat
    rw [UInt64.toNat_sub_of_le _ _ hleWord, hsub]
    rw [Nat.mod_eq_sub_mod (by omega), Nat.mod_eq_of_lt (by omega)]
    omega

theorem nmodSub_toNat (this : DenseUPolyZp) (a b : UInt64)
    (hp : this._p ≠ 0) (ha : a < this._p) (hb : b < this._p) :
    (dense_upoly_zp_nmod_sub_ir this a b).toNat =
      (a.toNat + this._p.toNat - b.toNat) % this._p.toNat := by
  have hpNat : 0 < this._p.toNat := by
    exact Nat.pos_of_ne_zero (fun h => hp (UInt64.toNat_inj.mp (by simpa using h)))
  have haNat : a.toNat < this._p.toNat := by
    simpa [UInt64.lt_iff_toNat_lt] using ha
  have hbNat : b.toNat < this._p.toNat := by
    simpa [UInt64.lt_iff_toNat_lt] using hb
  simp only [dense_upoly_zp_nmod_sub_ir]
  split
  next hle =>
    have hleNat : b.toNat ≤ a.toNat := by
      simpa [UInt64.le_iff_toNat_le] using hle
    rw [UInt64.toNat_sub_of_le _ _ hle]
    have heq : a.toNat + this._p.toNat - b.toNat =
        this._p.toNat + (a.toNat - b.toNat) := by omega
    have hdiffLt : a.toNat - b.toNat < this._p.toNat :=
      lt_of_le_of_lt (Nat.sub_le _ _) haNat
    rw [heq, Nat.add_mod, Nat.mod_self, zero_add]
    simp [Nat.mod_eq_of_lt hdiffLt]
  next hnot =>
    have hltNat : a.toNat < b.toNat := by
      have : ¬b.toNat ≤ a.toNat := by
        simpa [UInt64.le_iff_toNat_le] using hnot
      omega
    have hbLe : b ≤ this._p := by
      simpa [UInt64.le_iff_toNat_le] using Nat.le_of_lt hbNat
    have hsub : (this._p - b).toNat = this._p.toNat - b.toNat :=
      UInt64.toNat_sub_of_le _ _ hbLe
    have hsumP : (this._p.toNat - b.toNat) + a.toNat < this._p.toNat := by
      omega
    have hsum64 : (this._p.toNat - b.toNat) + a.toNat < UInt64.size :=
      lt_trans hsumP (UInt64.toNat_lt this._p)
    rw [UInt64.toNat_add, hsub, Nat.mod_eq_of_lt hsum64]
    have heq : a.toNat + this._p.toNat - b.toNat =
        (this._p.toNat - b.toNat) + a.toNat := by omega
    rw [heq, Nat.mod_eq_of_lt hsumP]

theorem nmodAdd_lt (this : DenseUPolyZp) (a b : UInt64)
    (hp : this._p ≠ 0) (ha : a < this._p) (hb : b < this._p) :
    dense_upoly_zp_nmod_add_ir this a b < this._p := by
  rw [UInt64.lt_iff_toNat_lt, nmodAdd_toNat this a b hp ha hb]
  exact Nat.mod_lt _ (Nat.pos_of_ne_zero
    (fun h => hp (UInt64.toNat_inj.mp (by simpa using h))))

theorem nmodSub_lt (this : DenseUPolyZp) (a b : UInt64)
    (hp : this._p ≠ 0) (ha : a < this._p) (hb : b < this._p) :
    dense_upoly_zp_nmod_sub_ir this a b < this._p := by
  rw [UInt64.lt_iff_toNat_lt, nmodSub_toNat this a b hp ha hb]
  exact Nat.mod_lt _ (Nat.pos_of_ne_zero
    (fun h => hp (UInt64.toNat_inj.mp (by simpa using h))))

theorem nmodAdd_cast (this : DenseUPolyZp) (a b : UInt64)
    (hp : this._p ≠ 0) (ha : a < this._p) (hb : b < this._p) :
    ((dense_upoly_zp_nmod_add_ir this a b).toNat : ZMod this._p.toNat) =
      (a.toNat : ZMod this._p.toNat) + (b.toNat : ZMod this._p.toNat) := by
  rw [nmodAdd_toNat this a b hp ha hb]
  rw [ZMod.natCast_mod]
  push_cast
  rfl

theorem nmodSub_cast (this : DenseUPolyZp) (a b : UInt64)
    (hp : this._p ≠ 0) (ha : a < this._p) (hb : b < this._p) :
    ((dense_upoly_zp_nmod_sub_ir this a b).toNat : ZMod this._p.toNat) =
      (a.toNat : ZMod this._p.toNat) - (b.toNat : ZMod this._p.toNat) := by
  rw [nmodSub_toNat this a b hp ha hb]
  rw [ZMod.natCast_mod, Nat.cast_sub (by
    have hbNat : b.toNat < this._p.toNat := by
      simpa [UInt64.lt_iff_toNat_lt] using hb
    omega), Nat.cast_add, ZMod.natCast_self, add_zero]

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
