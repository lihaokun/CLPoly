import CLPoly.Generated.StrictPolyAddSub
import CLPoly.Impl.StrictDivremRefinement
import CLPoly.Impl.StrictEuclidRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictPolyAddSubRefinement

open Generated.StrictPolyAddSub
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictEuclidRefinement

/-- The only pointer relationships supported by the source API: the output
starts at exactly the input address, or belongs to a different allocation. -/
def ExactOrDisjoint (out input : RawPtr UInt64) : Prop :=
  out = input ∨ out.region ≠ input.region

theorem address_ne_of_exactOrDisjoint {out input : RawPtr UInt64}
    {writeIndex readIndex : Nat} (h : ExactOrDisjoint out input)
    (hne : writeIndex ≠ readIndex) :
    out.region ≠ input.region ∨
      out.limbOffset + writeIndex ≠ input.limbOffset + readIndex := by
  rcases h with rfl | hregions
  · right
    omega
  · exact Or.inl hregions

theorem sameAddress_eq_true_iff (left right : RawPtr UInt64) :
    left.sameAddress right = true ↔ left = right := by
  cases left
  cases right
  simp [RawPtr.sameAddress]

theorem region_ne_of_exactOrDisjoint_not_sameAddress
    {out input : RawPtr UInt64} (h : ExactOrDisjoint out input)
    (hne : out.sameAddress input = false) : out.region ≠ input.region := by
  rcases h with heq | hregions
  · subst input
    simp [RawPtr.sameAddress] at hne
  · exact hregions

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

def nmodNeg (this : DenseUPolyZp) (b : UInt64) : UInt64 :=
  if b = 0 then 0 else this._p - b

theorem nmodNeg_lt (this : DenseUPolyZp) (b : UInt64)
    (hp : this._p ≠ 0) (hb : b < this._p) : nmodNeg this b < this._p := by
  unfold nmodNeg
  split
  next heq => simpa [heq] using (show (0 : UInt64) < this._p from UInt64.pos_iff_ne_zero.mpr hp)
  next hne =>
    rw [UInt64.lt_iff_toNat_lt]
    have hbNat : b.toNat < this._p.toNat := by
      simpa [UInt64.lt_iff_toNat_lt] using hb
    have hbLe : b ≤ this._p := by
      simpa [UInt64.le_iff_toNat_le] using Nat.le_of_lt hbNat
    rw [UInt64.toNat_sub_of_le _ _ hbLe]
    have hbPos : 0 < b.toNat := Nat.pos_of_ne_zero
      (fun hz => hne (UInt64.toNat_inj.mp (by simpa using hz)))
    omega

theorem nmodNeg_cast (this : DenseUPolyZp) (b : UInt64)
    (hp : this._p ≠ 0) (hb : b < this._p) :
    ((nmodNeg this b).toNat : ZMod this._p.toNat) =
      -(b.toNat : ZMod this._p.toNat) := by
  unfold nmodNeg
  split
  next heq => simp [heq]
  next hne =>
    have hbNat : b.toNat < this._p.toNat := by
      simpa [UInt64.lt_iff_toNat_lt] using hb
    have hbLe : b ≤ this._p := by
      simpa [UInt64.le_iff_toNat_le] using Nat.le_of_lt hbNat
    rw [UInt64.toNat_sub_of_le _ _ hbLe, Nat.cast_sub (Nat.le_of_lt hbNat),
      ZMod.natCast_self, zero_sub]

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

/-- The generated add loop writes only `C[i..limit)`.  This is the memory
fact needed to justify both exact in-place aliases admitted by the C++ API. -/
theorem addCommonLoop_preserves_outside (this : DenseUPolyZp)
    (C A B guard : RawPtr UInt64) (limit i readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hC : heap.ValidU64Slice C limit)
    (hA : heap.ValidU64Slice A limit)
    (hB : heap.ValidU64Slice B limit)
    (hread : heap.readU64 guard readIndex = .ok old)
    (houtside : ∀ k, i ≤ k → k < limit →
      C.region ≠ guard.region ∨
        C.limbOffset + k ≠ guard.limbOffset + readIndex)
    (hrun : addCommonLoop this C A B limit i heap = .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold addCommonLoop at hrun
  split at hrun
  next hlt =>
    rcases heap.readU64_of_valid A limit i hA hlt with ⟨a, ha⟩
    simp only [ha] at hrun
    rcases heap.readU64_of_valid B limit i hB hlt with ⟨b, hb⟩
    simp only [hb] at hrun
    rcases heap.writeU64_of_valid C limit i
      (dense_upoly_zp_nmod_add_ir this a b) hC hlt with ⟨heap1, hw⟩
    simp only [hw] at hrun
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 C guard
      i readIndex _ old hw hread (houtside i (Nat.le_refl _) hlt)
    have hlayout := RawHeap.writeU64_sameLayout heap heap1 C i _ hw
    apply addCommonLoop_preserves_outside this C A B guard limit (i + 1)
      readIndex heap1 heap' old
      ((hlayout C limit).mp hC) ((hlayout A limit).mp hA)
      ((hlayout B limit).mp hB) hread1
    · intro k hik hkl
      exact houtside k (by omega) hkl
    · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by limit - i
decreasing_by omega

theorem addCommonLoop_value (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (limit i k : Nat) (heap heap' : RawHeap)
    (a b : UInt64)
    (hC : heap.ValidU64Slice C limit)
    (hA : heap.ValidU64Slice A limit)
    (hB : heap.ValidU64Slice B limit)
    (hAliasA : ExactOrDisjoint C A) (hAliasB : ExactOrDisjoint C B)
    (hik : i ≤ k) (hkl : k < limit)
    (ha : heap.readU64 A k = .ok a)
    (hb : heap.readU64 B k = .ok b)
    (hrun : addCommonLoop this C A B limit i heap = .ok heap') :
    heap'.readU64 C k =
      .ok (dense_upoly_zp_nmod_add_ir this a b) := by
  unfold addCommonLoop at hrun
  split at hrun
  next hlt =>
    rcases heap.readU64_of_valid A limit i hA hlt with ⟨ai, hai⟩
    simp only [hai] at hrun
    rcases heap.readU64_of_valid B limit i hB hlt with ⟨bi, hbi⟩
    simp only [hbi] at hrun
    rcases heap.writeU64_of_valid C limit i
      (dense_upoly_zp_nmod_add_ir this ai bi) hC hlt with ⟨heap1, hw⟩
    simp only [hw] at hrun
    have hlayout := RawHeap.writeU64_sameLayout heap heap1 C i _ hw
    by_cases heq : k = i
    · subst k
      have haiEq : ai = a := Except.ok.inj (hai.symm.trans ha)
      have hbiEq : bi = b := Except.ok.inj (hbi.symm.trans hb)
      subst ai
      subst bi
      have hnow := RawHeap.readU64_writeU64_same heap heap1 C i
        (dense_upoly_zp_nmod_add_ir this a b) hw
      apply addCommonLoop_preserves_outside this C A B C limit (i + 1)
        i heap1 heap' (dense_upoly_zp_nmod_add_ir this a b)
        ((hlayout C limit).mp hC) ((hlayout A limit).mp hA)
        ((hlayout B limit).mp hB) hnow
      · intro j hij hjl
        right
        omega
      · exact hrun
    · have hik' : i + 1 ≤ k := by omega
      have ha1 := RawHeap.readU64_writeU64_ne heap heap1 C A i k _ a
        hw ha (address_ne_of_exactOrDisjoint hAliasA (Ne.symm heq))
      have hb1 := RawHeap.readU64_writeU64_ne heap heap1 C B i k _ b
        hw hb (address_ne_of_exactOrDisjoint hAliasB (Ne.symm heq))
      exact addCommonLoop_value this C A B limit (i + 1) k heap1 heap'
        a b ((hlayout C limit).mp hC) ((hlayout A limit).mp hA)
        ((hlayout B limit).mp hB) hAliasA hAliasB hik' hkl ha1 hb1 hrun
  next hnot => omega
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

theorem subCommonLoop_preserves_outside (this : DenseUPolyZp)
    (C A B guard : RawPtr UInt64) (limit i readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hC : heap.ValidU64Slice C limit)
    (hA : heap.ValidU64Slice A limit)
    (hB : heap.ValidU64Slice B limit)
    (hread : heap.readU64 guard readIndex = .ok old)
    (houtside : ∀ k, i ≤ k → k < limit →
      C.region ≠ guard.region ∨
        C.limbOffset + k ≠ guard.limbOffset + readIndex)
    (hrun : subCommonLoop this C A B limit i heap = .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold subCommonLoop at hrun
  split at hrun
  next hlt =>
    rcases heap.readU64_of_valid A limit i hA hlt with ⟨a, ha⟩
    simp only [ha] at hrun
    rcases heap.readU64_of_valid B limit i hB hlt with ⟨b, hb⟩
    simp only [hb] at hrun
    rcases heap.writeU64_of_valid C limit i
      (dense_upoly_zp_nmod_sub_ir this a b) hC hlt with ⟨heap1, hw⟩
    simp only [hw] at hrun
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 C guard
      i readIndex _ old hw hread (houtside i (Nat.le_refl _) hlt)
    have hlayout := RawHeap.writeU64_sameLayout heap heap1 C i _ hw
    apply subCommonLoop_preserves_outside this C A B guard limit (i + 1)
      readIndex heap1 heap' old
      ((hlayout C limit).mp hC) ((hlayout A limit).mp hA)
      ((hlayout B limit).mp hB) hread1
    · intro k hik hkl
      exact houtside k (by omega) hkl
    · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by limit - i
decreasing_by omega

theorem subCommonLoop_value (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (limit i k : Nat) (heap heap' : RawHeap)
    (a b : UInt64)
    (hC : heap.ValidU64Slice C limit)
    (hA : heap.ValidU64Slice A limit)
    (hB : heap.ValidU64Slice B limit)
    (hAliasA : ExactOrDisjoint C A) (hAliasB : ExactOrDisjoint C B)
    (hik : i ≤ k) (hkl : k < limit)
    (ha : heap.readU64 A k = .ok a)
    (hb : heap.readU64 B k = .ok b)
    (hrun : subCommonLoop this C A B limit i heap = .ok heap') :
    heap'.readU64 C k =
      .ok (dense_upoly_zp_nmod_sub_ir this a b) := by
  unfold subCommonLoop at hrun
  split at hrun
  next hlt =>
    rcases heap.readU64_of_valid A limit i hA hlt with ⟨ai, hai⟩
    simp only [hai] at hrun
    rcases heap.readU64_of_valid B limit i hB hlt with ⟨bi, hbi⟩
    simp only [hbi] at hrun
    rcases heap.writeU64_of_valid C limit i
      (dense_upoly_zp_nmod_sub_ir this ai bi) hC hlt with ⟨heap1, hw⟩
    simp only [hw] at hrun
    have hlayout := RawHeap.writeU64_sameLayout heap heap1 C i _ hw
    by_cases heq : k = i
    · subst k
      have haiEq : ai = a := Except.ok.inj (hai.symm.trans ha)
      have hbiEq : bi = b := Except.ok.inj (hbi.symm.trans hb)
      subst ai
      subst bi
      have hnow := RawHeap.readU64_writeU64_same heap heap1 C i
        (dense_upoly_zp_nmod_sub_ir this a b) hw
      apply subCommonLoop_preserves_outside this C A B C limit (i + 1)
        i heap1 heap' (dense_upoly_zp_nmod_sub_ir this a b)
        ((hlayout C limit).mp hC) ((hlayout A limit).mp hA)
        ((hlayout B limit).mp hB) hnow
      · intro j hij hjl
        right
        omega
      · exact hrun
    · have hik' : i + 1 ≤ k := by omega
      have ha1 := RawHeap.readU64_writeU64_ne heap heap1 C A i k _ a
        hw ha (address_ne_of_exactOrDisjoint hAliasA (Ne.symm heq))
      have hb1 := RawHeap.readU64_writeU64_ne heap heap1 C B i k _ b
        hw hb (address_ne_of_exactOrDisjoint hAliasB (Ne.symm heq))
      exact subCommonLoop_value this C A B limit (i + 1) k heap1 heap'
        a b ((hlayout C limit).mp hC) ((hlayout A limit).mp hA)
        ((hlayout B limit).mp hB) hAliasA hAliasB hik' hkl ha1 hb1 hrun
  next hnot => omega
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

theorem subNegTailLoop_preserves_outside (this : DenseUPolyZp)
    (C B guard : RawPtr UInt64) (limit i readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hC : heap.ValidU64Slice C limit)
    (hB : heap.ValidU64Slice B limit)
    (hread : heap.readU64 guard readIndex = .ok old)
    (houtside : ∀ k, i ≤ k → k < limit →
      C.region ≠ guard.region ∨
        C.limbOffset + k ≠ guard.limbOffset + readIndex)
    (hrun : subNegTailLoop this C B limit i heap = .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold subNegTailLoop at hrun
  split at hrun
  next hlt =>
    rcases heap.readU64_of_valid B limit i hB hlt with ⟨b, hb⟩
    simp only [hb] at hrun
    let value := if b = 0 then 0 else this._p - b
    rcases heap.writeU64_of_valid C limit i value hC hlt with ⟨heap1, hw⟩
    rw [show (if b = 0 then 0 else this._p - b) = value by rfl, hw] at hrun
    simp only at hrun
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 C guard
      i readIndex value old hw hread (houtside i (Nat.le_refl _) hlt)
    have hlayout := RawHeap.writeU64_sameLayout heap heap1 C i value hw
    apply subNegTailLoop_preserves_outside this C B guard limit (i + 1)
      readIndex heap1 heap' old
      ((hlayout C limit).mp hC) ((hlayout B limit).mp hB) hread1
    · intro k hik hkl
      exact houtside k (by omega) hkl
    · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by limit - i
decreasing_by omega

theorem subNegTailLoop_value (this : DenseUPolyZp)
    (C B : RawPtr UInt64) (limit i k : Nat) (heap heap' : RawHeap)
    (b : UInt64)
    (hC : heap.ValidU64Slice C limit)
    (hB : heap.ValidU64Slice B limit)
    (hAliasB : ExactOrDisjoint C B)
    (hik : i ≤ k) (hkl : k < limit)
    (hb : heap.readU64 B k = .ok b)
    (hrun : subNegTailLoop this C B limit i heap = .ok heap') :
    heap'.readU64 C k = .ok (nmodNeg this b) := by
  unfold subNegTailLoop at hrun
  split at hrun
  next hlt =>
    rcases heap.readU64_of_valid B limit i hB hlt with ⟨bi, hbi⟩
    simp only [hbi] at hrun
    let value := if bi = 0 then 0 else this._p - bi
    rcases heap.writeU64_of_valid C limit i value hC hlt with ⟨heap1, hw⟩
    rw [show (if bi = 0 then 0 else this._p - bi) = value by rfl, hw] at hrun
    simp only at hrun
    have hlayout := RawHeap.writeU64_sameLayout heap heap1 C i value hw
    by_cases heq : k = i
    · subst k
      have hbiEq : bi = b := Except.ok.inj (hbi.symm.trans hb)
      subst bi
      have hnow := RawHeap.readU64_writeU64_same heap heap1 C i value hw
      apply subNegTailLoop_preserves_outside this C B C limit (i + 1)
        i heap1 heap' (nmodNeg this b)
        ((hlayout C limit).mp hC) ((hlayout B limit).mp hB)
      · simpa [nmodNeg, value] using hnow
      · intro j hij hjl
        right
        omega
      · exact hrun
    · have hik' : i + 1 ≤ k := by omega
      have hb1 := RawHeap.readU64_writeU64_ne heap heap1 C B i k value b
        hw hb (address_ne_of_exactOrDisjoint hAliasB (Ne.symm heq))
      exact subNegTailLoop_value this C B limit (i + 1) k heap1 heap'
        b ((hlayout C limit).mp hC) ((hlayout B limit).mp hB)
        hAliasB hik' hkl hb1 hrun
  next hnot => omega
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

/-- Complete equal-length branch of the actual C++ `_poly_add`, including
in-place execution and the source normalization call. -/
theorem polyAdd_equalLength_refines (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (length : Nat) (heap heap' : RawHeap)
    (outLen : Nat) (left right : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hC : heap.ValidU64Slice C length)
    (hLeft : RawDensePolyRep this heap A length left)
    (hRight : RawDensePolyRep this heap B length right)
    (hAliasA : ExactOrDisjoint C A) (hAliasB : ExactOrDisjoint C B)
    (hrun : dense_upoly_zp__poly_add_ir this C A length B length heap =
      .ok (heap', outLen)) :
    RawDensePolyRep this heap' C outLen (left + right) := by
  rcases hLeft with ⟨hA, hcanonA, hrepA, hnormA⟩
  rcases hRight with ⟨hB, hcanonB, hrepB, hnormB⟩
  rcases addCommonLoop_ok this C A B length 0 heap hC hA hB (by omega) with
    ⟨heap1, hloop, hlayout⟩
  have hC1 : heap1.ValidU64Slice C length := (hlayout C length).mp hC
  rcases normaliseU64_ok heap1 C length hC1 with
    ⟨result, hnorm, hresultLe⟩
  have hrun' := hrun
  simp only [dense_upoly_zp__poly_add_ir, min_self, max_self, hloop,
    gt_iff_lt, lt_self_iff_false, ↓reduceIte, hnorm] at hrun'
  have hheap : heap' = heap1 := by
    exact congrArg Prod.fst (Except.ok.inj hrun').symm
  have hlen : outLen = result := by
    exact congrArg Prod.snd (Except.ok.inj hrun').symm
  subst heap'
  subst outLen
  rcases slicePolyRep_exists_unique heap1 C length this._p.toNat hC1 with
    ⟨output, hrepOutput, _⟩
  have houtput : output = left + right := by
    ext degree
    by_cases hd : degree < length
    · rcases slicePolyRep_coeff heap A length this._p.toNat left hrepA
        degree hd with ⟨a, ha, hcoeffA⟩
      rcases slicePolyRep_coeff heap B length this._p.toNat right hrepB
        degree hd with ⟨b, hb, hcoeffB⟩
      rcases slicePolyRep_coeff heap1 C length this._p.toNat output
        hrepOutput degree hd with ⟨c, hc, hcoeffC⟩
      have hvalue := addCommonLoop_value this C A B length 0 degree
        heap heap1 a b hC hA hB hAliasA hAliasB (by omega) hd ha hb hloop
      have hcEq : c = dense_upoly_zp_nmod_add_ir this a b :=
        Except.ok.inj (hc.symm.trans hvalue)
      have haLt : a < this._p := by
        simpa [UInt64.lt_iff_toNat_lt] using hcanonA degree a hd ha
      have hbLt : b < this._p := by
        simpa [UInt64.lt_iff_toNat_lt] using hcanonB degree b hd hb
      rw [Polynomial.coeff_add, hcoeffC, hcEq,
        nmodAdd_cast this a b hp haLt hbLt, hcoeffA, hcoeffB]
    · rw [slicePolyRep_coeff_zero_of_length_le heap1 C length
          this._p.toNat output hrepOutput degree (by omega),
        Polynomial.coeff_add,
        slicePolyRep_coeff_zero_of_length_le heap A length
          this._p.toNat left hrepA degree (by omega),
        slicePolyRep_coeff_zero_of_length_le heap B length
          this._p.toNat right hrepB degree (by omega), add_zero]
  subst output
  have hcanonicalFull : CanonicalU64Prefix heap1 C length this._p := by
    intro k value hk hread
    rcases slicePolyRep_coeff heap A length this._p.toNat left hrepA k hk with
      ⟨a, ha, _⟩
    rcases slicePolyRep_coeff heap B length this._p.toNat right hrepB k hk with
      ⟨b, hb, _⟩
    have hvalue := addCommonLoop_value this C A B length 0 k heap heap1
      a b hC hA hB hAliasA hAliasB (by omega) hk ha hb hloop
    have heq : value = dense_upoly_zp_nmod_add_ir this a b :=
      Except.ok.inj (hread.symm.trans hvalue)
    subst value
    apply (show dense_upoly_zp_nmod_add_ir this a b < this._p from
      nmodAdd_lt this a b hp
        (by simpa [UInt64.lt_iff_toNat_lt] using hcanonA k a hk ha)
        (by simpa [UInt64.lt_iff_toNat_lt] using hcanonB k b hk hb))
  have hvalidResult := heap1.validU64Slice_mono C length result hC1 hresultLe
  refine ⟨hvalidResult, ?_, ?_, ?_⟩
  · intro k value hk hread
    exact hcanonicalFull k value (lt_of_lt_of_le hk hresultLe) hread
  · exact slicePolyRep_of_normaliseU64 heap1 C length this._p.toNat result
      (left + right) hC1 hrepOutput hnorm
  · exact normaliseU64_result_fixed heap1 C length result hC1 hnorm

theorem polySub_equalLength_refines (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (length : Nat) (heap heap' : RawHeap)
    (outLen : Nat) (left right : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hC : heap.ValidU64Slice C length)
    (hLeft : RawDensePolyRep this heap A length left)
    (hRight : RawDensePolyRep this heap B length right)
    (hAliasA : ExactOrDisjoint C A) (hAliasB : ExactOrDisjoint C B)
    (hrun : dense_upoly_zp__poly_sub_ir this C A length B length heap =
      .ok (heap', outLen)) :
    RawDensePolyRep this heap' C outLen (left - right) := by
  rcases hLeft with ⟨hA, hcanonA, hrepA, hnormA⟩
  rcases hRight with ⟨hB, hcanonB, hrepB, hnormB⟩
  rcases subCommonLoop_ok this C A B length 0 heap hC hA hB (by omega) with
    ⟨heap1, hloop, hlayout⟩
  have hC1 : heap1.ValidU64Slice C length := (hlayout C length).mp hC
  rcases normaliseU64_ok heap1 C length hC1 with
    ⟨result, hnorm, hresultLe⟩
  have hrun' := hrun
  simp only [dense_upoly_zp__poly_sub_ir, min_self, max_self, hloop,
    gt_iff_lt, lt_self_iff_false, ↓reduceIte, hnorm] at hrun'
  have hheap : heap' = heap1 := by
    exact congrArg Prod.fst (Except.ok.inj hrun').symm
  have hlen : outLen = result := by
    exact congrArg Prod.snd (Except.ok.inj hrun').symm
  subst heap'
  subst outLen
  rcases slicePolyRep_exists_unique heap1 C length this._p.toNat hC1 with
    ⟨output, hrepOutput, _⟩
  have houtput : output = left - right := by
    ext degree
    by_cases hd : degree < length
    · rcases slicePolyRep_coeff heap A length this._p.toNat left hrepA
        degree hd with ⟨a, ha, hcoeffA⟩
      rcases slicePolyRep_coeff heap B length this._p.toNat right hrepB
        degree hd with ⟨b, hb, hcoeffB⟩
      rcases slicePolyRep_coeff heap1 C length this._p.toNat output
        hrepOutput degree hd with ⟨c, hc, hcoeffC⟩
      have hvalue := subCommonLoop_value this C A B length 0 degree
        heap heap1 a b hC hA hB hAliasA hAliasB (by omega) hd ha hb hloop
      have hcEq : c = dense_upoly_zp_nmod_sub_ir this a b :=
        Except.ok.inj (hc.symm.trans hvalue)
      have haLt : a < this._p := by
        simpa [UInt64.lt_iff_toNat_lt] using hcanonA degree a hd ha
      have hbLt : b < this._p := by
        simpa [UInt64.lt_iff_toNat_lt] using hcanonB degree b hd hb
      rw [Polynomial.coeff_sub, hcoeffC, hcEq,
        nmodSub_cast this a b hp haLt hbLt, hcoeffA, hcoeffB]
    · rw [slicePolyRep_coeff_zero_of_length_le heap1 C length
          this._p.toNat output hrepOutput degree (by omega),
        Polynomial.coeff_sub,
        slicePolyRep_coeff_zero_of_length_le heap A length
          this._p.toNat left hrepA degree (by omega),
        slicePolyRep_coeff_zero_of_length_le heap B length
          this._p.toNat right hrepB degree (by omega), sub_zero]
  subst output
  have hcanonicalFull : CanonicalU64Prefix heap1 C length this._p := by
    intro k value hk hread
    rcases slicePolyRep_coeff heap A length this._p.toNat left hrepA k hk with
      ⟨a, ha, _⟩
    rcases slicePolyRep_coeff heap B length this._p.toNat right hrepB k hk with
      ⟨b, hb, _⟩
    have hvalue := subCommonLoop_value this C A B length 0 k heap heap1
      a b hC hA hB hAliasA hAliasB (by omega) hk ha hb hloop
    have heq : value = dense_upoly_zp_nmod_sub_ir this a b :=
      Except.ok.inj (hread.symm.trans hvalue)
    subst value
    apply (show dense_upoly_zp_nmod_sub_ir this a b < this._p from
      nmodSub_lt this a b hp
        (by simpa [UInt64.lt_iff_toNat_lt] using hcanonA k a hk ha)
        (by simpa [UInt64.lt_iff_toNat_lt] using hcanonB k b hk hb))
  have hvalidResult := heap1.validU64Slice_mono C length result hC1 hresultLe
  refine ⟨hvalidResult, ?_, ?_, ?_⟩
  · intro k value hk hread
    exact hcanonicalFull k value (lt_of_lt_of_le hk hresultLe) hread
  · exact slicePolyRep_of_normaliseU64 heap1 C length this._p.toNat result
      (left - right) hC1 hrepOutput hnorm
  · exact normaliseU64_result_fixed heap1 C length result hC1 hnorm

end CLPoly.Impl.StrictPolyAddSubRefinement
