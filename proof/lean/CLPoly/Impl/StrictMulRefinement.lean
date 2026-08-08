import CLPoly.Generated.StrictMul
import CLPoly.Impl.StrictDivremRefinement
import CLPoly.Impl.StrictEuclidRefinement
import CLPoly.Impl.StrictPolyAddSubRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictMulRefinement

open Generated.StrictMul
open Generated.StrictPolyAddSub
open CLPoly.Impl.StrictWordArithmetic
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictPolyAddSubRefinement

theorem slicePolyRep_prefix_exists (heap : RawHeap) (ptr : RawPtr UInt64)
    (length prefixLength p : Nat) (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice ptr length)
    (hle : prefixLength ≤ length)
    (hrep : SlicePolyRep heap ptr length p poly) :
    ∃ prefixPoly : Polynomial (ZMod p),
      SlicePolyRep heap ptr prefixLength p prefixPoly ∧
      ∀ i, i < prefixLength → prefixPoly.coeff i = poly.coeff i := by
  have hvalidPrefix := heap.validU64Slice_mono ptr length prefixLength hvalid hle
  rcases slicePolyRep_exists_unique heap ptr prefixLength p hvalidPrefix with
    ⟨prefixPoly, hprefix, _⟩
  refine ⟨prefixPoly, hprefix, ?_⟩
  intro i hi
  rcases slicePolyRep_coeff heap ptr prefixLength p prefixPoly hprefix i hi with
    ⟨prefixValue, hreadPrefix, hcoeffPrefix⟩
  rcases slicePolyRep_coeff heap ptr length p poly hrep i (by omega) with
    ⟨fullValue, hreadFull, hcoeffFull⟩
  have heq : prefixValue = fullValue :=
    Except.ok.inj (hreadPrefix.symm.trans hreadFull)
  rw [hcoeffPrefix, hcoeffFull, heq]

theorem karAddHalvesLoop_ok (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m i : Nat) (heap : RawHeap)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m) :
    ∃ heap', karAddHalvesLoop this A B t1 t2 m i heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  unfold karAddHalvesLoop
  split
  next hi =>
    rcases heap.readU64_of_valid A (2 * m) i hA (by omega) with ⟨alo, halo⟩
    simp only [halo]
    rcases heap.readU64_of_valid A (2 * m) (m + i) hA (by omega) with
      ⟨ahi, hahi⟩
    simp only [hahi]
    rcases heap.readU64_of_valid B (2 * m) i hB (by omega) with ⟨blo, hblo⟩
    simp only [hblo]
    rcases heap.readU64_of_valid B (2 * m) (m + i) hB (by omega) with
      ⟨bhi, hbhi⟩
    simp only [hbhi]
    let av := dense_upoly_zp_nmod_add_ir this alo ahi
    rcases heap.writeU64_of_valid t1 m i av hT1 hi with ⟨heap1, hw1⟩
    simp only [av, hw1]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 i av hw1
    let bv := dense_upoly_zp_nmod_add_ir this blo bhi
    rcases heap1.writeU64_of_valid t2 m i bv ((hlayout1 t2 m).mp hT2) hi with
      ⟨heap2, hw2⟩
    simp only [bv, hw2]
    have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 t2 i bv hw2
    rcases karAddHalvesLoop_ok this A B t1 t2 m (i + 1) heap2
      ((hlayout2 A (2 * m)).mp ((hlayout1 A (2 * m)).mp hA))
      ((hlayout2 B (2 * m)).mp ((hlayout1 B (2 * m)).mp hB))
      ((hlayout2 t1 m).mp ((hlayout1 t1 m).mp hT1))
      ((hlayout2 t2 m).mp ((hlayout1 t2 m).mp hT2)) with
      ⟨heap3, hrun, hlayout3⟩
    rw [hrun]
    exact ⟨heap3, rfl, fun ptr length =>
      (hlayout1 ptr length).trans
        ((hlayout2 ptr length).trans (hlayout3 ptr length))⟩
  next => exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by m - i
decreasing_by omega

theorem karSubLoop_ok (this : DenseUPolyZp)
    (dst sub : RawPtr UInt64) (count i : Nat) (heap : RawHeap)
    (hDst : heap.ValidU64Slice dst count)
    (hSub : heap.ValidU64Slice sub count) :
    ∃ heap', karSubLoop this dst sub count i heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  unfold karSubLoop
  split
  next hi =>
    rcases heap.readU64_of_valid dst count i hDst hi with ⟨a, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid sub count i hSub hi with ⟨b, hb⟩
    simp only [hb]
    let value := dense_upoly_zp_nmod_sub_ir this a b
    rcases heap.writeU64_of_valid dst count i value hDst hi with ⟨heap1, hw⟩
    simp only [value, hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 dst i value hw
    rcases karSubLoop_ok this dst sub count (i + 1) heap1
      ((hlayout1 dst count).mp hDst) ((hlayout1 sub count).mp hSub) with
      ⟨heap2, hrun, hlayout2⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next => exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by count - i
decreasing_by omega

theorem karAssembleLoop_ok (this : DenseUPolyZp)
    (C sP1 : RawPtr UInt64) (m count i : Nat) (heap : RawHeap)
    (hC : heap.ValidU64Slice C (m + count))
    (hP1 : heap.ValidU64Slice sP1 count) :
    ∃ heap', karAssembleLoop this C sP1 m count i heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  unfold karAssembleLoop
  split
  next hi =>
    rcases heap.readU64_of_valid C (m + count) (m + i) hC (by omega) with
      ⟨base, hbase⟩
    simp only [hbase]
    rcases heap.readU64_of_valid sP1 count i hP1 hi with ⟨cross, hcross⟩
    simp only [hcross]
    let value := dense_upoly_zp_nmod_add_ir this base cross
    rcases heap.writeU64_of_valid C (m + count) (m + i) value hC (by omega) with
      ⟨heap1, hw⟩
    simp only [value, hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C (m + i) value hw
    rcases karAssembleLoop_ok this C sP1 m count (i + 1) heap1
      ((hlayout1 C (m + count)).mp hC) ((hlayout1 sP1 count).mp hP1) with
      ⟨heap2, hrun, hlayout2⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next => exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩
termination_by count - i
decreasing_by omega

theorem karOddTail_ok (A B t1 t2 : RawPtr UInt64) (m h : Nat)
    (heap : RawHeap)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h) :
    ∃ heap', karOddTail A B t1 t2 m h heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  unfold karOddTail
  split
  next hodd =>
    rcases heap.readU64_of_valid A (m + h) (m + m) hA (by omega) with
      ⟨aTail, ha⟩
    simp only [ha]
    rcases heap.readU64_of_valid B (m + h) (m + m) hB (by omega) with
      ⟨bTail, hb⟩
    simp only [hb]
    rcases heap.writeU64_of_valid t1 h m aTail hT1 hodd with ⟨heap1, hw1⟩
    simp only [hw1]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 m aTail hw1
    rcases heap1.writeU64_of_valid t2 h m bTail ((hlayout1 t2 h).mp hT2)
      hodd with ⟨heap2, hw2⟩
    simp only [hw2]
    have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 t2 m bTail hw2
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  next => exact ⟨heap, rfl, fun _ _ => Iff.rfl⟩

theorem karOddTail_preserves_region_ne (A B t1 t2 guard : RawPtr UInt64)
    (m h guardLen : Nat) (heap heap' : RawHeap)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1Guard : t1.region ≠ guard.region)
    (hT2Guard : t2.region ≠ guard.region)
    (hrun : karOddTail A B t1 t2 m h heap = .ok heap') :
    SameU64Prefix heap heap' guard guardLen := by
  intro i old hi hread
  unfold karOddTail at hrun
  split at hrun
  next hodd =>
    rcases heap.readU64_of_valid A (m + h) (m + m) hA (by omega) with
      ⟨aTail, ha⟩
    simp only [ha] at hrun
    rcases heap.readU64_of_valid B (m + h) (m + m) hB (by omega) with
      ⟨bTail, hb⟩
    simp only [hb] at hrun
    rcases heap.writeU64_of_valid t1 h m aTail hT1 hodd with ⟨heap1, hw1⟩
    simp only [hw1] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 m aTail hw1
    rcases heap1.writeU64_of_valid t2 h m bTail ((hlayout1 t2 h).mp hT2)
      hodd with ⟨heap2, hw2⟩
    simp only [hw2] at hrun
    have heq : heap' = heap2 := Except.ok.inj hrun.symm
    subst heap'
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 guard m i
      aTail old hw1 hread (Or.inl hT1Guard)
    exact RawHeap.readU64_writeU64_ne heap1 heap2 t2 guard m i bTail old
      hw2 hread1 (Or.inl hT2Guard)
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread

theorem karPrepareHalves_ok (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m h : Nat) (heap : RawHeap)
    (hmh : m ≤ h)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h) :
    ∃ heap', karPrepareHalves this A B t1 t2 m h heap = .ok heap' ∧
      RawHeap.SameLayout heap heap' := by
  have hA2m := heap.validU64Slice_mono A (m + h) (2 * m) hA (by omega)
  have hB2m := heap.validU64Slice_mono B (m + h) (2 * m) hB (by omega)
  have hT1m := heap.validU64Slice_mono t1 h m hT1 hmh
  have hT2m := heap.validU64Slice_mono t2 h m hT2 hmh
  rcases karAddHalvesLoop_ok this A B t1 t2 m 0 heap hA2m hB2m hT1m
      hT2m with ⟨heap1, hadd, hlayout1⟩
  have hA1 := (hlayout1 A (m + h)).mp hA
  have hB1 := (hlayout1 B (m + h)).mp hB
  have hT11 := (hlayout1 t1 h).mp hT1
  have hT21 := (hlayout1 t2 h).mp hT2
  rcases karOddTail_ok A B t1 t2 m h heap1 hA1 hB1 hT11 hT21 with
    ⟨heap2, htail, hlayout2⟩
  refine ⟨heap2, ?_, fun ptr length =>
    (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
  simp [karPrepareHalves, hadd, htail]

theorem karOddTail_values (A B t1 t2 : RawPtr UInt64) (m h : Nat)
    (heap heap' : RawHeap)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1T2 : t1.region ≠ t2.region)
    (hodd : h > m)
    (hrun : karOddTail A B t1 t2 m h heap = .ok heap') :
    ∃ aTail bTail,
      heap.readU64 A (m + m) = .ok aTail ∧
      heap.readU64 B (m + m) = .ok bTail ∧
      heap'.readU64 t1 m = .ok aTail ∧
      heap'.readU64 t2 m = .ok bTail := by
  unfold karOddTail at hrun
  simp only [hodd, ↓reduceIte] at hrun
  rcases heap.readU64_of_valid A (m + h) (m + m) hA (by omega) with
    ⟨aTail, ha⟩
  simp only [ha] at hrun
  rcases heap.readU64_of_valid B (m + h) (m + m) hB (by omega) with
    ⟨bTail, hb⟩
  simp only [hb] at hrun
  rcases heap.writeU64_of_valid t1 h m aTail hT1 hodd with ⟨heap1, hw1⟩
  simp only [hw1] at hrun
  have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 m aTail hw1
  rcases heap1.writeU64_of_valid t2 h m bTail ((hlayout1 t2 h).mp hT2)
    hodd with ⟨heap2, hw2⟩
  simp only [hw2] at hrun
  have heq : heap' = heap2 := Except.ok.inj hrun.symm
  subst heap'
  have hread1 := RawHeap.readU64_writeU64_same heap heap1 t1 m aTail hw1
  have hread1' := RawHeap.readU64_writeU64_ne heap1 heap2 t2 t1 m m
    bTail aTail hw2 hread1 (Or.inl (Ne.symm hT1T2))
  have hread2 := RawHeap.readU64_writeU64_same heap1 heap2 t2 m bTail hw2
  exact ⟨aTail, bTail, ha, hb, hread1', hread2⟩

theorem karOddTail_preserves_own_prefixes
    (A B t1 t2 : RawPtr UInt64) (m h : Nat) (heap heap' : RawHeap)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1T2 : t1.region ≠ t2.region)
    (hrun : karOddTail A B t1 t2 m h heap = .ok heap') :
    SameU64Prefix heap heap' t1 m ∧ SameU64Prefix heap heap' t2 m := by
  constructor <;> intro i old hi hread
  · unfold karOddTail at hrun
    split at hrun
    next hodd =>
      rcases heap.readU64_of_valid A (m + h) (m + m) hA (by omega) with
        ⟨aTail, ha⟩
      simp only [ha] at hrun
      rcases heap.readU64_of_valid B (m + h) (m + m) hB (by omega) with
        ⟨bTail, hb⟩
      simp only [hb] at hrun
      rcases heap.writeU64_of_valid t1 h m aTail hT1 hodd with ⟨heap1, hw1⟩
      simp only [hw1] at hrun
      have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 m aTail hw1
      rcases heap1.writeU64_of_valid t2 h m bTail ((hlayout1 t2 h).mp hT2)
        hodd with ⟨heap2, hw2⟩
      simp only [hw2] at hrun
      have heq : heap' = heap2 := Except.ok.inj hrun.symm
      subst heap'
      have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 t1 m i
        aTail old hw1 hread (Or.inr (by omega))
      exact RawHeap.readU64_writeU64_ne heap1 heap2 t2 t1 m i bTail old
        hw2 hread1 (Or.inl (Ne.symm hT1T2))
    next hnot =>
      have heq : heap' = heap := Except.ok.inj hrun.symm
      simpa [heq] using hread
  · unfold karOddTail at hrun
    split at hrun
    next hodd =>
      rcases heap.readU64_of_valid A (m + h) (m + m) hA (by omega) with
        ⟨aTail, ha⟩
      simp only [ha] at hrun
      rcases heap.readU64_of_valid B (m + h) (m + m) hB (by omega) with
        ⟨bTail, hb⟩
      simp only [hb] at hrun
      rcases heap.writeU64_of_valid t1 h m aTail hT1 hodd with ⟨heap1, hw1⟩
      simp only [hw1] at hrun
      have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 m aTail hw1
      rcases heap1.writeU64_of_valid t2 h m bTail ((hlayout1 t2 h).mp hT2)
        hodd with ⟨heap2, hw2⟩
      simp only [hw2] at hrun
      have heq : heap' = heap2 := Except.ok.inj hrun.symm
      subst heap'
      have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 t2 m i
        aTail old hw1 hread (Or.inl hT1T2)
      exact RawHeap.readU64_writeU64_ne heap1 heap2 t2 t2 m i bTail old
        hw2 hread1 (Or.inr (by omega))
    next hnot =>
      have heq : heap' = heap := Except.ok.inj hrun.symm
      simpa [heq] using hread

theorem karOddTail_coeffs (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m h : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1T2 : t1.region ≠ t2.region)
    (hCanonicalA : CanonicalU64Prefix heap A (m + h) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (m + h) this._p)
    (hRepA : SlicePolyRep heap A (m + h) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (m + h) this._p.toNat right)
    (hodd : h > m)
    (hrun : karOddTail A B t1 t2 m h heap = .ok heap') :
    ∃ aTail bTail,
      heap'.readU64 t1 m = .ok aTail ∧
      heap'.readU64 t2 m = .ok bTail ∧
      (aTail.toNat : ZMod this._p.toNat) = left.coeff (m + m) ∧
      (bTail.toNat : ZMod this._p.toNat) = right.coeff (m + m) ∧
      aTail.toNat < this._p.toNat ∧ bTail.toNat < this._p.toNat := by
  rcases karOddTail_values A B t1 t2 m h heap heap' hA hB hT1 hT2
      hT1T2 hodd hrun with ⟨aTail, bTail, ha, hb, ht1, ht2⟩
  rcases slicePolyRep_coeff heap A (m + h) this._p.toNat left hRepA
      (m + m) (by omega) with ⟨a', ha', hcoeffA⟩
  have haEq : a' = aTail := Except.ok.inj (ha'.symm.trans ha)
  subst a'
  rcases slicePolyRep_coeff heap B (m + h) this._p.toNat right hRepB
      (m + m) (by omega) with ⟨b', hb', hcoeffB⟩
  have hbEq : b' = bTail := Except.ok.inj (hb'.symm.trans hb)
  subst b'
  exact ⟨aTail, bTail, ht1, ht2, hcoeffA.symm, hcoeffB.symm,
    hCanonicalA (m + m) aTail (by omega) ha,
    hCanonicalB (m + m) bTail (by omega) hb⟩

theorem karAddHalvesLoop_preserves_outside (this : DenseUPolyZp)
    (A B t1 t2 guard : RawPtr UInt64) (m i readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m)
    (hread : heap.readU64 guard readIndex = .ok old)
    (hout1 : ∀ j, i ≤ j → j < m → t1.region ≠ guard.region ∨
      t1.limbOffset + j ≠ guard.limbOffset + readIndex)
    (hout2 : ∀ j, i ≤ j → j < m → t2.region ≠ guard.region ∨
      t2.limbOffset + j ≠ guard.limbOffset + readIndex)
    (hrun : karAddHalvesLoop this A B t1 t2 m i heap = .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold karAddHalvesLoop at hrun
  split at hrun
  next hi =>
    rcases heap.readU64_of_valid A (2 * m) i hA (by omega) with ⟨alo, halo⟩
    simp only [halo] at hrun
    rcases heap.readU64_of_valid A (2 * m) (m + i) hA (by omega) with
      ⟨ahi, hahi⟩
    simp only [hahi] at hrun
    rcases heap.readU64_of_valid B (2 * m) i hB (by omega) with ⟨blo, hblo⟩
    simp only [hblo] at hrun
    rcases heap.readU64_of_valid B (2 * m) (m + i) hB (by omega) with
      ⟨bhi, hbhi⟩
    simp only [hbhi] at hrun
    let av := dense_upoly_zp_nmod_add_ir this alo ahi
    rcases heap.writeU64_of_valid t1 m i av hT1 hi with ⟨heap1, hw1⟩
    simp only [av, hw1] at hrun
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 guard i
      readIndex av old hw1 hread (hout1 i (Nat.le_refl _) hi)
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 i av hw1
    let bv := dense_upoly_zp_nmod_add_ir this blo bhi
    rcases heap1.writeU64_of_valid t2 m i bv ((hlayout1 t2 m).mp hT2) hi with
      ⟨heap2, hw2⟩
    simp only [bv, hw2] at hrun
    have hread2 := RawHeap.readU64_writeU64_ne heap1 heap2 t2 guard i
      readIndex bv old hw2 hread1 (hout2 i (Nat.le_refl _) hi)
    have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 t2 i bv hw2
    apply karAddHalvesLoop_preserves_outside this A B t1 t2 guard m (i + 1)
      readIndex heap2 heap' old
      ((hlayout2 A (2 * m)).mp ((hlayout1 A (2 * m)).mp hA))
      ((hlayout2 B (2 * m)).mp ((hlayout1 B (2 * m)).mp hB))
      ((hlayout2 t1 m).mp ((hlayout1 t1 m).mp hT1))
      ((hlayout2 t2 m).mp ((hlayout1 t2 m).mp hT2)) hread2
    · intro j hij hjm
      exact hout1 j (by omega) hjm
    · intro j hij hjm
      exact hout2 j (by omega) hjm
    · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by m - i
decreasing_by omega

theorem karAddHalvesLoop_current_values (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m i : Nat) (heap heap' : RawHeap)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m)
    (hT1T2 : t1.region ≠ t2.region)
    (hi : i < m)
    (hrun : karAddHalvesLoop this A B t1 t2 m i heap = .ok heap') :
    ∃ alo ahi blo bhi,
      heap.readU64 A i = .ok alo ∧
      heap.readU64 A (m + i) = .ok ahi ∧
      heap.readU64 B i = .ok blo ∧
      heap.readU64 B (m + i) = .ok bhi ∧
      heap'.readU64 t1 i =
        .ok (dense_upoly_zp_nmod_add_ir this alo ahi) ∧
      heap'.readU64 t2 i =
        .ok (dense_upoly_zp_nmod_add_ir this blo bhi) := by
  unfold karAddHalvesLoop at hrun
  simp only [hi, ↓reduceDIte] at hrun
  rcases heap.readU64_of_valid A (2 * m) i hA (by omega) with ⟨alo, halo⟩
  simp only [halo] at hrun
  rcases heap.readU64_of_valid A (2 * m) (m + i) hA (by omega) with
    ⟨ahi, hahi⟩
  simp only [hahi] at hrun
  rcases heap.readU64_of_valid B (2 * m) i hB (by omega) with ⟨blo, hblo⟩
  simp only [hblo] at hrun
  rcases heap.readU64_of_valid B (2 * m) (m + i) hB (by omega) with
    ⟨bhi, hbhi⟩
  simp only [hbhi] at hrun
  let av := dense_upoly_zp_nmod_add_ir this alo ahi
  rcases heap.writeU64_of_valid t1 m i av hT1 hi with ⟨heap1, hw1⟩
  simp only [av, hw1] at hrun
  have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 i av hw1
  let bv := dense_upoly_zp_nmod_add_ir this blo bhi
  rcases heap1.writeU64_of_valid t2 m i bv ((hlayout1 t2 m).mp hT2) hi with
    ⟨heap2, hw2⟩
  simp only [bv, hw2] at hrun
  have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 t2 i bv hw2
  have hav1 := RawHeap.readU64_writeU64_same heap heap1 t1 i av hw1
  have hav2 := RawHeap.readU64_writeU64_ne heap1 heap2 t2 t1 i i bv av
    hw2 hav1 (Or.inl (Ne.symm hT1T2))
  have hbv2 := RawHeap.readU64_writeU64_same heap1 heap2 t2 i bv hw2
  have havFinal := karAddHalvesLoop_preserves_outside this A B t1 t2 t1 m
    (i + 1) i heap2 heap' av
    ((hlayout2 A (2 * m)).mp ((hlayout1 A (2 * m)).mp hA))
    ((hlayout2 B (2 * m)).mp ((hlayout1 B (2 * m)).mp hB))
    ((hlayout2 t1 m).mp ((hlayout1 t1 m).mp hT1))
    ((hlayout2 t2 m).mp ((hlayout1 t2 m).mp hT2)) hav2
    (by intro j _ _; exact Or.inr (by omega))
    (by intro _ _ _; exact Or.inl (Ne.symm hT1T2)) hrun
  have hbvFinal := karAddHalvesLoop_preserves_outside this A B t1 t2 t2 m
    (i + 1) i heap2 heap' bv
    ((hlayout2 A (2 * m)).mp ((hlayout1 A (2 * m)).mp hA))
    ((hlayout2 B (2 * m)).mp ((hlayout1 B (2 * m)).mp hB))
    ((hlayout2 t1 m).mp ((hlayout1 t1 m).mp hT1))
    ((hlayout2 t2 m).mp ((hlayout1 t2 m).mp hT2)) hbv2
    (by intro _ _ _; exact Or.inl hT1T2)
    (by intro j _ _; exact Or.inr (by omega)) hrun
  exact ⟨alo, ahi, blo, bhi, halo, hahi, hblo, hbhi,
    by simpa [av] using havFinal, by simpa [bv] using hbvFinal⟩

theorem karAddHalvesLoop_current_coeffs (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m i : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m)
    (hT1T2 : t1.region ≠ t2.region)
    (hCanonicalA : CanonicalU64Prefix heap A (2 * m) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (2 * m) this._p)
    (hRepA : SlicePolyRep heap A (2 * m) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (2 * m) this._p.toNat right)
    (hi : i < m)
    (hrun : karAddHalvesLoop this A B t1 t2 m i heap = .ok heap') :
    ∃ value1 value2,
      heap'.readU64 t1 i = .ok value1 ∧
      heap'.readU64 t2 i = .ok value2 ∧
      (value1.toNat : ZMod this._p.toNat) =
        left.coeff i + left.coeff (m + i) ∧
      (value2.toNat : ZMod this._p.toNat) =
        right.coeff i + right.coeff (m + i) ∧
      value1.toNat < this._p.toNat ∧ value2.toNat < this._p.toNat := by
  rcases karAddHalvesLoop_current_values this A B t1 t2 m i heap heap'
      hA hB hT1 hT2 hT1T2 hi hrun with
    ⟨alo, ahi, blo, bhi, halo, hahi, hblo, hbhi, ht1, ht2⟩
  have haloLt : alo < this._p := by
    simpa [UInt64.lt_iff_toNat_lt] using hCanonicalA i alo (by omega) halo
  have hahiLt : ahi < this._p := by
    simpa [UInt64.lt_iff_toNat_lt] using hCanonicalA (m + i) ahi (by omega) hahi
  have hbloLt : blo < this._p := by
    simpa [UInt64.lt_iff_toNat_lt] using hCanonicalB i blo (by omega) hblo
  have hbhiLt : bhi < this._p := by
    simpa [UInt64.lt_iff_toNat_lt] using hCanonicalB (m + i) bhi (by omega) hbhi
  have hpWord : this._p ≠ 0 := by
    intro hzero
    have hzeroNat := congrArg UInt64.toNat hzero
    simp at hzeroNat
    omega
  rcases slicePolyRep_coeff heap A (2 * m) this._p.toNat left hRepA i
      (by omega) with ⟨alo', halo', hcoeffAlo⟩
  have haloEq : alo' = alo := Except.ok.inj (halo'.symm.trans halo)
  subst alo'
  rcases slicePolyRep_coeff heap A (2 * m) this._p.toNat left hRepA (m + i)
      (by omega) with ⟨ahi', hahi', hcoeffAhi⟩
  have hahiEq : ahi' = ahi := Except.ok.inj (hahi'.symm.trans hahi)
  subst ahi'
  rcases slicePolyRep_coeff heap B (2 * m) this._p.toNat right hRepB i
      (by omega) with ⟨blo', hblo', hcoeffBlo⟩
  have hbloEq : blo' = blo := Except.ok.inj (hblo'.symm.trans hblo)
  subst blo'
  rcases slicePolyRep_coeff heap B (2 * m) this._p.toNat right hRepB (m + i)
      (by omega) with ⟨bhi', hbhi', hcoeffBhi⟩
  have hbhiEq : bhi' = bhi := Except.ok.inj (hbhi'.symm.trans hbhi)
  subst bhi'
  let value1 := dense_upoly_zp_nmod_add_ir this alo ahi
  let value2 := dense_upoly_zp_nmod_add_ir this blo bhi
  refine ⟨value1, value2, by simpa [value1] using ht1,
    by simpa [value2] using ht2, ?_, ?_, ?_, ?_⟩
  · simpa [value1, hcoeffAlo, hcoeffAhi] using
      nmodAdd_cast this alo ahi hpWord haloLt hahiLt
  · simpa [value2, hcoeffBlo, hcoeffBhi] using
      nmodAdd_cast this blo bhi hpWord hbloLt hbhiLt
  · simpa [value1, UInt64.lt_iff_toNat_lt] using
      nmodAdd_lt this alo ahi hpWord haloLt hahiLt
  · simpa [value2, UInt64.lt_iff_toNat_lt] using
      nmodAdd_lt this blo bhi hpWord hbloLt hbhiLt

theorem karAddHalvesLoop_coeffs (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m i : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m)
    (hT1T2 : t1.region ≠ t2.region)
    (hT1A : t1.region ≠ A.region) (hT1B : t1.region ≠ B.region)
    (hT2A : t2.region ≠ A.region) (hT2B : t2.region ≠ B.region)
    (hCanonicalA : CanonicalU64Prefix heap A (2 * m) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (2 * m) this._p)
    (hRepA : SlicePolyRep heap A (2 * m) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (2 * m) this._p.toNat right)
    (hrun : karAddHalvesLoop this A B t1 t2 m i heap = .ok heap') :
    ∀ k, i ≤ k → k < m → ∃ value1 value2,
      heap'.readU64 t1 k = .ok value1 ∧
      heap'.readU64 t2 k = .ok value2 ∧
      (value1.toNat : ZMod this._p.toNat) =
        left.coeff k + left.coeff (m + k) ∧
      (value2.toNat : ZMod this._p.toNat) =
        right.coeff k + right.coeff (m + k) ∧
      value1.toNat < this._p.toNat ∧ value2.toNat < this._p.toNat := by
  intro k hik hkm
  by_cases heq : k = i
  · subst k
    exact karAddHalvesLoop_current_coeffs this A B t1 t2 m i heap heap'
      left right hp hA hB hT1 hT2 hT1T2 hCanonicalA hCanonicalB hRepA
      hRepB hkm hrun
  · have hik' : i + 1 ≤ k := by omega
    unfold karAddHalvesLoop at hrun
    have hi : i < m := by omega
    simp only [hi, ↓reduceDIte] at hrun
    rcases heap.readU64_of_valid A (2 * m) i hA (by omega) with ⟨alo, halo⟩
    simp only [halo] at hrun
    rcases heap.readU64_of_valid A (2 * m) (m + i) hA (by omega) with
      ⟨ahi, hahi⟩
    simp only [hahi] at hrun
    rcases heap.readU64_of_valid B (2 * m) i hB (by omega) with ⟨blo, hblo⟩
    simp only [hblo] at hrun
    rcases heap.readU64_of_valid B (2 * m) (m + i) hB (by omega) with
      ⟨bhi, hbhi⟩
    simp only [hbhi] at hrun
    let av := dense_upoly_zp_nmod_add_ir this alo ahi
    rcases heap.writeU64_of_valid t1 m i av hT1 hi with ⟨heap1, hw1⟩
    simp only [av, hw1] at hrun
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 t1 i av hw1
    let bv := dense_upoly_zp_nmod_add_ir this blo bhi
    rcases heap1.writeU64_of_valid t2 m i bv ((hlayout1 t2 m).mp hT2) hi with
      ⟨heap2, hw2⟩
    simp only [bv, hw2] at hrun
    have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 t2 i bv hw2
    have hA1 := (hlayout1 A (2 * m)).mp hA
    have hB1 := (hlayout1 B (2 * m)).mp hB
    have hA2 := (hlayout2 A (2 * m)).mp hA1
    have hB2 := (hlayout2 B (2 * m)).mp hB1
    have hT12 := (hlayout2 t1 m).mp ((hlayout1 t1 m).mp hT1)
    have hT22 := (hlayout2 t2 m).mp ((hlayout1 t2 m).mp hT2)
    have hsameA : SameU64Prefix heap heap2 A (2 * m) := by
      intro j old hj hread
      have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 A i j av old
        hw1 hread (Or.inl hT1A)
      exact RawHeap.readU64_writeU64_ne heap1 heap2 t2 A i j bv old
        hw2 hread1 (Or.inl hT2A)
    have hsameB : SameU64Prefix heap heap2 B (2 * m) := by
      intro j old hj hread
      have hread1 := RawHeap.readU64_writeU64_ne heap heap1 t1 B i j av old
        hw1 hread (Or.inl hT1B)
      exact RawHeap.readU64_writeU64_ne heap1 heap2 t2 B i j bv old
        hw2 hread1 (Or.inl hT2B)
    have hCanonicalA2 : CanonicalU64Prefix heap2 A (2 * m) this._p := by
      intro j value hj hread2
      rcases heap.readU64_of_valid A (2 * m) j hA hj with ⟨old, hread⟩
      have hpreserved := hsameA j old hj hread
      have heq : value = old := Except.ok.inj (hread2.symm.trans hpreserved)
      subst value
      exact hCanonicalA j old hj hread
    have hCanonicalB2 : CanonicalU64Prefix heap2 B (2 * m) this._p := by
      intro j value hj hread2
      rcases heap.readU64_of_valid B (2 * m) j hB hj with ⟨old, hread⟩
      have hpreserved := hsameB j old hj hread
      have heq : value = old := Except.ok.inj (hread2.symm.trans hpreserved)
      subst value
      exact hCanonicalB j old hj hread
    have hRepA2 := slicePolyRep_of_same_prefix heap heap2 A (2 * m)
      this._p.toNat left hA hA2 hsameA hRepA
    have hRepB2 := slicePolyRep_of_same_prefix heap heap2 B (2 * m)
      this._p.toNat right hB hB2 hsameB hRepB
    exact karAddHalvesLoop_coeffs this A B t1 t2 m (i + 1) heap2 heap'
      left right hp hA2 hB2 hT12 hT22 hT1T2 hT1A hT1B hT2A hT2B
      hCanonicalA2 hCanonicalB2 hRepA2 hRepB2 hrun k hik' hkm
termination_by m - i
decreasing_by omega

theorem karSubLoop_preserves_outside (this : DenseUPolyZp)
    (dst sub guard : RawPtr UInt64) (count i readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hDst : heap.ValidU64Slice dst count)
    (hSub : heap.ValidU64Slice sub count)
    (hread : heap.readU64 guard readIndex = .ok old)
    (hout : ∀ j, i ≤ j → j < count → dst.region ≠ guard.region ∨
      dst.limbOffset + j ≠ guard.limbOffset + readIndex)
    (hrun : karSubLoop this dst sub count i heap = .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold karSubLoop at hrun
  split at hrun
  next hi =>
    rcases heap.readU64_of_valid dst count i hDst hi with ⟨a, ha⟩
    simp only [ha] at hrun
    rcases heap.readU64_of_valid sub count i hSub hi with ⟨b, hb⟩
    simp only [hb] at hrun
    let value := dense_upoly_zp_nmod_sub_ir this a b
    rcases heap.writeU64_of_valid dst count i value hDst hi with ⟨heap1, hw⟩
    simp only [value, hw] at hrun
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 dst guard i
      readIndex value old hw hread (hout i (Nat.le_refl _) hi)
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 dst i value hw
    apply karSubLoop_preserves_outside this dst sub guard count (i + 1)
      readIndex heap1 heap' old ((hlayout1 dst count).mp hDst)
      ((hlayout1 sub count).mp hSub) hread1
    · intro j hij hjc
      exact hout j (by omega) hjc
    · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by count - i
decreasing_by omega

theorem karAssembleLoop_preserves_outside (this : DenseUPolyZp)
    (C sP1 guard : RawPtr UInt64) (m count i readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hC : heap.ValidU64Slice C (m + count))
    (hP1 : heap.ValidU64Slice sP1 count)
    (hread : heap.readU64 guard readIndex = .ok old)
    (hout : ∀ j, i ≤ j → j < count → C.region ≠ guard.region ∨
      C.limbOffset + (m + j) ≠ guard.limbOffset + readIndex)
    (hrun : karAssembleLoop this C sP1 m count i heap = .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold karAssembleLoop at hrun
  split at hrun
  next hi =>
    rcases heap.readU64_of_valid C (m + count) (m + i) hC (by omega) with
      ⟨base, hbase⟩
    simp only [hbase] at hrun
    rcases heap.readU64_of_valid sP1 count i hP1 hi with ⟨cross, hcross⟩
    simp only [hcross] at hrun
    let value := dense_upoly_zp_nmod_add_ir this base cross
    rcases heap.writeU64_of_valid C (m + count) (m + i) value hC (by omega) with
      ⟨heap1, hw⟩
    simp only [value, hw] at hrun
    have hread1 := RawHeap.readU64_writeU64_ne heap heap1 C guard (m + i)
      readIndex value old hw hread (hout i (Nat.le_refl _) hi)
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C (m + i) value hw
    apply karAssembleLoop_preserves_outside this C sP1 guard m count (i + 1)
      readIndex heap1 heap' old ((hlayout1 C (m + count)).mp hC)
      ((hlayout1 sP1 count).mp hP1) hread1
    · intro j hij hjc
      exact hout j (by omega) hjc
    · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by count - i
decreasing_by omega

theorem karPrepareHalves_preserves_region_ne (this : DenseUPolyZp)
    (A B t1 t2 guard : RawPtr UInt64) (m h guardLen : Nat)
    (heap heap' : RawHeap)
    (hmh : m ≤ h)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1Guard : t1.region ≠ guard.region)
    (hT2Guard : t2.region ≠ guard.region)
    (hrun : karPrepareHalves this A B t1 t2 m h heap = .ok heap') :
    SameU64Prefix heap heap' guard guardLen := by
  have hA2m := heap.validU64Slice_mono A (m + h) (2 * m) hA (by omega)
  have hB2m := heap.validU64Slice_mono B (m + h) (2 * m) hB (by omega)
  have hT1m := heap.validU64Slice_mono t1 h m hT1 hmh
  have hT2m := heap.validU64Slice_mono t2 h m hT2 hmh
  rcases karAddHalvesLoop_ok this A B t1 t2 m 0 heap hA2m hB2m hT1m
      hT2m with ⟨heap1, hadd, hlayout1⟩
  have htail : karOddTail A B t1 t2 m h heap1 = .ok heap' := by
    simpa [karPrepareHalves, hadd] using hrun
  have hsameAdd : SameU64Prefix heap heap1 guard guardLen := by
    intro i old hi hread
    apply karAddHalvesLoop_preserves_outside this A B t1 t2 guard m 0 i
      heap heap1 old hA2m hB2m hT1m hT2m hread
    · intro _ _ _
      exact Or.inl hT1Guard
    · intro _ _ _
      exact Or.inl hT2Guard
    · exact hadd
  have hsameTail := karOddTail_preserves_region_ne A B t1 t2 guard m h
    guardLen heap1 heap'
    ((hlayout1 A (m + h)).mp hA) ((hlayout1 B (m + h)).mp hB)
    ((hlayout1 t1 h).mp hT1) ((hlayout1 t2 h).mp hT2)
    hT1Guard hT2Guard htail
  intro i old hi hread
  exact hsameTail i old hi (hsameAdd i old hi hread)

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

def classicalDotPoly {p : Nat} (left right : Polynomial (ZMod p))
    (k stop j : Nat) : ZMod p :=
  if h : j ≤ stop then
    left.coeff j * right.coeff (k - j) +
      classicalDotPoly left right k stop (j + 1)
  else
    0
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotPoly_eq_sum_Icc {p : Nat}
    (left right : Polynomial (ZMod p)) (k stop j : Nat) :
    classicalDotPoly left right k stop j =
      ∑ t ∈ Finset.Icc j stop, left.coeff t * right.coeff (k - t) := by
  rw [classicalDotPoly]
  split
  next hle =>
    rw [← Finset.insert_Icc_succ_left_eq_Icc hle]
    simp
    exact classicalDotPoly_eq_sum_Icc left right k stop (j + 1)
  next hnot =>
    symm
    apply Finset.sum_eq_zero
    intro t ht
    simp only [Finset.mem_Icc] at ht
    omega
termination_by stop + 1 - j
decreasing_by omega

theorem classicalDotPoly_source_eq_coeff {p : Nat}
    (heap : RawHeap) (A B : RawPtr UInt64) (lenA lenB k : Nat)
    (left right : Polynomial (ZMod p))
    (hLenA : 0 < lenA)
    (hRepA : SlicePolyRep heap A lenA p left)
    (hRepB : SlicePolyRep heap B lenB p right) :
    classicalDotPoly left right k
      (if k < lenA then k else lenA - 1)
      (if k ≥ lenB then k - lenB + 1 else 0) =
        (left * right).coeff k := by
  rw [classicalDotPoly_eq_sum_Icc]
  rw [Polynomial.coeff_mul,
    Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  by_cases hkA : k < lenA <;> by_cases hkB : k ≥ lenB
  all_goals simp only [hkA, hkB, if_true, if_false]
  all_goals
  apply Finset.sum_subset
  · intro t ht
    simp only [Finset.mem_Icc] at ht
    simp only [Finset.mem_range]
    omega
  · intro t htRange htIcc
    simp only [Finset.mem_range] at htRange
    simp only [Finset.mem_Icc, not_and_or, not_le] at htIcc
    rcases htIcc with hBelow | hAbove
    · by_cases hkt : lenB ≤ k
      · have hzero := slicePolyRep_coeff_zero_of_length_le heap B lenB p
          right hRepB (k - t) (by omega)
        rw [hzero, mul_zero]
      · omega
    · by_cases hkt : k < lenA
      · omega
      · have hzero := slicePolyRep_coeff_zero_of_length_le heap A lenA p
          left hRepA t (by omega)
        rw [hzero, zero_mul]

theorem classicalDotNat_cast_eq_poly (heap : RawHeap)
    (A B : RawPtr UInt64) (lenA lenB p k stop j sum : Nat)
    (left right : Polynomial (ZMod p))
    (hRepA : SlicePolyRep heap A lenA p left)
    (hRepB : SlicePolyRep heap B lenB p right)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotNat heap A B k stop j = .ok sum) :
    (sum : ZMod p) = classicalDotPoly left right k stop j := by
  unfold classicalDotNat at hrun
  split at hrun
  next hle =>
    rw [classicalDotPoly, dif_pos hle]
    have hjA := hAIndex j (Nat.le_refl _) hle
    have hjB := hBIndex j (Nat.le_refl _) hle
    rcases slicePolyRep_coeff heap A lenA p left hRepA j hjA with
      ⟨a, ha, hcoeffA⟩
    simp only [ha] at hrun
    rcases slicePolyRep_coeff heap B lenB p right hRepB (k - j) hjB with
      ⟨b, hb, hcoeffB⟩
    simp only [hb] at hrun
    cases ht : classicalDotNat heap A B k stop (j + 1) with
    | error fault => simp [ht] at hrun
    | ok tail =>
      simp only [ht] at hrun
      have hsum : sum = a.toNat * b.toNat + tail :=
        Except.ok.inj hrun.symm
      have htail := classicalDotNat_cast_eq_poly heap A B lenA lenB p k
        stop (j + 1) tail left right hRepA hRepB
        (by intro t hjt hts; exact hAIndex t (by omega) hts)
        (by intro t hjt hts; exact hBIndex t (by omega) hts) ht
      rw [hsum, Nat.cast_add, Nat.cast_mul, hcoeffA, hcoeffB, htail]
  next hnot =>
    rw [classicalDotPoly, dif_neg hnot]
    have hzero : sum = 0 := Except.ok.inj hrun.symm
    subst sum
    simp
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

theorem classicalDotNat_bound (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (p : UInt64) (sum : Nat)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotNat heap A B k stop j = .ok sum) :
    sum ≤ (stop + 1 - j) * (p.toNat - 1) ^ 2 := by
  unfold classicalDotNat at hrun
  split at hrun
  next hle =>
    have hjA := hAIndex j (Nat.le_refl _) hle
    have hjB := hBIndex j (Nat.le_refl _) hle
    rcases heap.readU64_of_valid A lenA j hA hjA with ⟨a, ha⟩
    simp only [ha] at hrun
    rcases heap.readU64_of_valid B lenB (k - j) hB hjB with ⟨b, hb⟩
    simp only [hb] at hrun
    cases ht : classicalDotNat heap A B k stop (j + 1) with
    | error fault => simp [ht] at hrun
    | ok tail =>
      simp only [ht] at hrun
      have hsum : sum = a.toNat * b.toNat + tail :=
        Except.ok.inj hrun.symm
      have haLe : a.toNat ≤ p.toNat - 1 := by
        have := hCanonicalA j a hjA ha
        omega
      have hbLe : b.toNat ≤ p.toNat - 1 := by
        have := hCanonicalB (k - j) b hjB hb
        omega
      have hprod : a.toNat * b.toNat ≤ (p.toNat - 1) ^ 2 := by
        rw [pow_two]
        exact Nat.mul_le_mul haLe hbLe
      have htail := classicalDotNat_bound heap A B lenA lenB k stop (j + 1)
        p tail hA hB hCanonicalA hCanonicalB
        (by intro t hjt hts; exact hAIndex t (by omega) hts)
        (by intro t hjt hts; exact hBIndex t (by omega) hts) ht
      rw [hsum]
      calc
        a.toNat * b.toNat + tail ≤
            (p.toNat - 1) ^ 2 +
              (stop + 1 - (j + 1)) * (p.toNat - 1) ^ 2 :=
          Nat.add_le_add hprod htail
        _ = (stop + 1 - j) * (p.toNat - 1) ^ 2 := by
          have : stop + 1 - j = (stop + 1 - (j + 1)) + 1 := by omega
          rw [this]
          ring
  next hnot =>
    have hzero : sum = 0 := Except.ok.inj hrun.symm
    subst sum
    exact Nat.zero_le _
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

theorem classicalDotLoop_exact_zero (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (p : UInt64) (result : Word3) (sum : Nat)
    (hp : 1 < p.toNat)
    (hcount : stop + 1 - j < limbBase)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotLoop heap A B k stop j
      { lo := 0, mid := 0, hi := 0 } = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    word3Value result = sum := by
  have hbound := classicalDotNat_bound heap A B lenA lenB k stop j p sum
    hA hB hCanonicalA hCanonicalB hAIndex hBIndex hsum
  have hlazy := lazyAccumulation_word3_budget p (stop + 1 - j) 0
    hp hcount (by omega)
  have hpB : p.toNat < limbBase := by
    simpa [limbBase] using UInt64.toNat_lt p
  have hsumLt : sum < limbBase ^ 3 := by
    calc
      sum ≤ (stop + 1 - j) * (p.toNat - 1) ^ 2 := hbound
      _ < p.toNat * limbBase ^ 2 := by simpa using hlazy
      _ < limbBase ^ 3 := by
        have hpowPos : 0 < limbBase ^ 2 :=
          pow_pos (by norm_num [limbBase]) 2
        have hmul : p.toNat * limbBase ^ 2 < limbBase * limbBase ^ 2 :=
          Nat.mul_lt_mul_of_pos_right hpB hpowPos
        simpa [pow_succ, Nat.mul_comm, Nat.mul_left_comm,
          Nat.mul_assoc] using hmul
  have hmod := classicalDotLoop_modEq heap A B k stop j
    { lo := 0, mid := 0, hi := 0 } result sum hrun hsum
  have hmod' : Nat.ModEq (limbBase ^ 3) (word3Value result) sum := by
    simpa [word3Value] using hmod
  exact hmod'.eq_of_lt_of_lt (word3Value_lt result) hsumLt

theorem classicalDotReduced_toNat (this : DenseUPolyZp)
    (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (result : Word3) (sum : Nat)
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hcount : stop + 1 - j < limbBase)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotLoop heap A B k stop j
      { lo := 0, mid := 0, hi := 0 } = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    (Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      result.hi result.mid result.lo this._p this._ninv this._norm).toNat =
      sum % this._p.toNat := by
  have hexact := classicalDotLoop_exact_zero heap A B lenA lenB k stop j
    this._p result sum hp hcount hA hB hCanonicalA hCanonicalB
    hAIndex hBIndex hrun hsum
  have hsumBound := classicalDotNat_bound heap A B lenA lenB k stop j
    this._p sum hA hB hCanonicalA hCanonicalB hAIndex hBIndex hsum
  have hlazy := lazyAccumulation_word3_budget this._p (stop + 1 - j) 0
    hp hcount (by omega)
  have hvalue : word3Value result < this._p.toNat * limbBase ^ 2 := by
    rw [hexact]
    exact lt_of_le_of_lt hsumBound (by simpa using hlazy)
  have hhi : result.hi.toNat < this._p.toNat :=
    word3_hi_lt_of_value_lt result this._p hvalue
  rw [lll_mod_preinv_ir_correct_of_configured this result.hi result.mid
    result.lo hcfg hhi, hexact]

theorem classicalDotReduced_cast (this : DenseUPolyZp)
    (heap : RawHeap) (A B : RawPtr UInt64)
    (lenA lenB k stop j : Nat) (result : Word3) (sum : Nat)
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hcount : stop + 1 - j < limbBase)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hAIndex : ∀ t, j ≤ t → t ≤ stop → t < lenA)
    (hBIndex : ∀ t, j ≤ t → t ≤ stop → k - t < lenB)
    (hrun : classicalDotLoop heap A B k stop j
      { lo := 0, mid := 0, hi := 0 } = .ok result)
    (hsum : classicalDotNat heap A B k stop j = .ok sum) :
    ((Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      result.hi result.mid result.lo this._p this._ninv this._norm).toNat :
        ZMod this._p.toNat) = (sum : ZMod this._p.toNat) := by
  rw [classicalDotReduced_toNat this heap A B lenA lenB k stop j result sum
    hcfg hp hcount hA hB hCanonicalA hCanonicalB hAIndex hBIndex hrun hsum,
    ZMod.natCast_mod]

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

theorem classicalReduced_source_eq_coeff (this : DenseUPolyZp)
    (heap : RawHeap) (A B : RawPtr UInt64) (lenA lenB lenC k : Nat)
    (left right : Polynomial (ZMod this._p.toNat)) (acc : Word3)
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hLenAWord : lenA < limbBase)
    (hlenC : lenC = lenA + lenB - 1) (hk : k < lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hRepA : SlicePolyRep heap A lenA this._p.toNat left)
    (hRepB : SlicePolyRep heap B lenB this._p.toNat right)
    (hdot : classicalDotLoop heap A B k
      (if k < lenA then k else lenA - 1)
      (if k ≥ lenB then k - lenB + 1 else 0)
      { lo := 0, mid := 0, hi := 0 } = .ok acc) :
    ((Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
      acc.hi acc.mid acc.lo this._p this._ninv this._norm).toNat :
        ZMod this._p.toNat) = (left * right).coeff k ∧
      (Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        acc.hi acc.mid acc.lo this._p this._ninv this._norm).toNat <
          this._p.toNat := by
  let jMin := if k ≥ lenB then k - lenB + 1 else 0
  let jMax := if k < lenA then k else lenA - 1
  have hkSource : k < lenA + lenB - 1 := by omega
  have hrange : jMin ≤ jMax := by
    dsimp [jMin, jMax]
    split <;> split <;> omega
  have hAIndex : ∀ t, jMin ≤ t → t ≤ jMax → t < lenA := by
    intro t hjt hts
    exact (classical_index_bounds lenA lenB k t hApos hBpos hkSource
      hjt hts).1
  have hBIndex : ∀ t, jMin ≤ t → t ≤ jMax → k - t < lenB := by
    intro t hjt hts
    exact (classical_index_bounds lenA lenB k t hApos hBpos hkSource
      hjt hts).2
  rcases classicalDotNat_ok heap A B lenA lenB k jMax jMin hA hB
    hAIndex hBIndex with ⟨sum, hsum⟩
  have hcount : jMax + 1 - jMin < limbBase := by
    have hjMaxA := hAIndex jMax hrange (Nat.le_refl _)
    omega
  constructor
  · calc
      ((Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        acc.hi acc.mid acc.lo this._p this._ninv this._norm).toNat :
          ZMod this._p.toNat) = (sum : ZMod this._p.toNat) :=
        classicalDotReduced_cast this heap A B lenA lenB k jMax jMin acc sum
          hcfg hp hcount hA hB hCanonicalA hCanonicalB hAIndex hBIndex
          (by simpa [jMin, jMax] using hdot) hsum
      _ = classicalDotPoly left right k jMax jMin :=
        classicalDotNat_cast_eq_poly heap A B lenA lenB this._p.toNat k
          jMax jMin sum left right hRepA hRepB hAIndex hBIndex hsum
      _ = (left * right).coeff k := by
        simpa [jMin, jMax] using classicalDotPoly_source_eq_coeff heap A B
          lenA lenB k left right hApos hRepA hRepB
  · rw [classicalDotReduced_toNat this heap A B lenA lenB k jMax jMin acc sum
      hcfg hp hcount hA hB hCanonicalA hCanonicalB hAIndex hBIndex
      (by simpa [jMin, jMax] using hdot) hsum]
    exact Nat.mod_lt sum (by omega)

/-- The generated schoolbook outer loop writes only `C[k..lenC)`.  This
memory fact keeps both source buffers and already-produced coefficients tied
to their original raw cells throughout the actual C++ loop. -/
theorem classicalOuterLoop_preserves_outside (this : DenseUPolyZp)
    (C A B guard : RawPtr UInt64) (lenA lenB lenC k readIndex : Nat)
    (heap heap' : RawHeap) (old : UInt64)
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hread : heap.readU64 guard readIndex = .ok old)
    (houtside : ∀ i, k ≤ i → i < lenC →
      C.region ≠ guard.region ∨
        C.limbOffset + i ≠ guard.limbOffset + readIndex)
    (hrun : classicalOuterLoop this C A B lenA lenB lenC k heap =
      .ok heap') :
    heap'.readU64 guard readIndex = .ok old := by
  unfold classicalOuterLoop at hrun
  split at hrun
  next hk =>
    let jMin := if k ≥ lenB then k - lenB + 1 else 0
    let jMax := if k < lenA then k else lenA - 1
    cases hdot : classicalDotLoop heap A B k jMax jMin
        { lo := 0, mid := 0, hi := 0 } with
    | error fault => simp [jMin, jMax, hdot] at hrun
    | ok acc =>
      simp only [jMin, jMax, hdot] at hrun
      let value := Generated.StrictGCD.dense_upoly_zp__lll_mod_preinv_ir
        acc.hi acc.mid acc.lo this._p this._ninv this._norm
      rcases heap.writeU64_of_valid C lenC k value hC hk with ⟨heap1, hw⟩
      simp only [value, hw] at hrun
      have hread1 := RawHeap.readU64_writeU64_ne heap heap1 C guard
        k readIndex value old hw hread (houtside k (Nat.le_refl _) hk)
      have hlayout := RawHeap.writeU64_sameLayout heap heap1 C k value hw
      apply classicalOuterLoop_preserves_outside this C A B guard lenA lenB
        lenC (k + 1) readIndex heap1 heap' old
        ((hlayout C lenC).mp hC) ((hlayout A lenA).mp hA)
        ((hlayout B lenB).mp hB) hread1
      · intro i hik hil
        exact houtside i (by omega) hil
      · exact hrun
  next hnot =>
    have heq : heap' = heap := Except.ok.inj hrun.symm
    simpa [heq] using hread
termination_by lenC - k
decreasing_by omega

theorem classicalOuterLoop_same_prefix_region_ne (this : DenseUPolyZp)
    (C A B guard : RawPtr UInt64) (lenA lenB lenC k guardLen : Nat)
    (heap heap' : RawHeap)
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hregions : C.region ≠ guard.region)
    (hrun : classicalOuterLoop this C A B lenA lenB lenC k heap =
      .ok heap') :
    SameU64Prefix heap heap' guard guardLen := by
  intro i old hi hread
  apply classicalOuterLoop_preserves_outside this C A B guard lenA lenB
    lenC k i heap heap' old hC hA hB hread
  · intro _ _ _
    exact Or.inl hregions
  · exact hrun

theorem canonicalU64Prefix_of_same_prefix (before after : RawHeap)
    (ptr : RawPtr UInt64) (length : Nat) (p : UInt64)
    (hvalid : before.ValidU64Slice ptr length)
    (hsame : SameU64Prefix before after ptr length)
    (hcanonical : CanonicalU64Prefix before ptr length p) :
    CanonicalU64Prefix after ptr length p := by
  intro k value hk hreadAfter
  rcases before.readU64_of_valid ptr length k
      hvalid hk with ⟨old, hreadBefore⟩
  have hpreserved := hsame k old hk hreadBefore
  have hvalue : value = old := Except.ok.inj (hreadAfter.symm.trans hpreserved)
  subst value
  exact hcanonical k old hk hreadBefore

/-- Coefficient-level output invariant for the prefix already produced by
the C++ schoolbook outer loop. -/
def ClassicalCoeffPrefix {p : Nat} (heap : RawHeap)
    (C : RawPtr UInt64) (upto : Nat) (poly : Polynomial (ZMod p)) : Prop :=
  ∀ i, i < upto → ∃ value : UInt64,
    heap.readU64 C i = .ok value ∧
      (value.toNat : ZMod p) = poly.coeff i ∧ value.toNat < p

theorem slicePolyRep_of_classicalCoeffPrefix {p : Nat}
    (heap : RawHeap) (C : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod p))
    (hvalid : heap.ValidU64Slice C length)
    (hprefix : ClassicalCoeffPrefix heap C length poly)
    (hzero : ∀ i, length ≤ i → poly.coeff i = 0) :
    SlicePolyRep heap C length p poly := by
  rcases slicePolyRep_exists_unique heap C length p hvalid with
    ⟨observed, hObserved, _⟩
  have heq : observed = poly := by
    ext i
    by_cases hi : i < length
    · rcases hprefix i hi with ⟨value, hread, hcoeff, _⟩
      rcases slicePolyRep_coeff heap C length p observed hObserved i hi with
        ⟨observedValue, hObservedRead, hObservedCoeff⟩
      have hvalue : observedValue = value :=
        Except.ok.inj (hObservedRead.symm.trans hread)
      rw [hObservedCoeff, hvalue, hcoeff]
    · rw [slicePolyRep_coeff_zero_of_length_le heap C length p observed
          hObserved i (by omega),
        hzero i (by omega)]
  rw [heq] at hObserved
  exact hObserved

theorem canonicalU64Prefix_of_classicalCoeffPrefix
    (heap : RawHeap) (C : RawPtr UInt64) (length : Nat)
    (modulus : UInt64) (poly : Polynomial (ZMod modulus.toNat))
    (hprefix : ClassicalCoeffPrefix heap C length poly) :
    CanonicalU64Prefix heap C length modulus := by
  intro i value hi hread
  rcases hprefix i hi with ⟨stored, hstored, _, hstoredLt⟩
  have heq : value = stored := Except.ok.inj (hread.symm.trans hstored)
  subst value
  exact hstoredLt

noncomputable def karHalfSumPoly {p : Nat}
    (poly : Polynomial (ZMod p)) (m : Nat) : Polynomial (ZMod p) :=
  ∑ i : Fin m,
    Polynomial.monomial i.val (poly.coeff i.val + poly.coeff (m + i.val))

theorem coeff_karHalfSumPoly {p : Nat}
    (poly : Polynomial (ZMod p)) (m degree : Nat) :
    (karHalfSumPoly poly m).coeff degree =
      if degree < m then poly.coeff degree + poly.coeff (m + degree) else 0 := by
  classical
  unfold karHalfSumPoly
  rw [Polynomial.finset_sum_coeff]
  by_cases hdegree : degree < m
  · rw [if_pos hdegree, Finset.sum_eq_single ⟨degree, hdegree⟩]
    · simp
    · intro index _ hne
      have hval : index.val ≠ degree := by
        intro heq
        apply hne
        exact Fin.ext heq
      simp [Polynomial.coeff_monomial, hval]
    · simp
  · rw [if_neg hdegree]
    apply Finset.sum_eq_zero
    intro index _
    have hval : index.val ≠ degree := by
      intro heq
      apply hdegree
      simpa [heq] using index.isLt
    simp [Polynomial.coeff_monomial, hval]

theorem karHalfSumPoly_congr_of_prefix {p : Nat}
    (prefixPoly fullPoly : Polynomial (ZMod p)) (m : Nat)
    (hcoeff : ∀ i, i < 2 * m → prefixPoly.coeff i = fullPoly.coeff i) :
    karHalfSumPoly prefixPoly m = karHalfSumPoly fullPoly m := by
  ext degree
  rw [coeff_karHalfSumPoly, coeff_karHalfSumPoly]
  by_cases hdegree : degree < m
  · rw [if_pos hdegree, if_pos hdegree,
      hcoeff degree (by omega), hcoeff (m + degree) (by omega)]
  · rw [if_neg hdegree, if_neg hdegree]

theorem karAddHalvesLoop_refines_slices (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat)
    (hA : heap.ValidU64Slice A (2 * m))
    (hB : heap.ValidU64Slice B (2 * m))
    (hT1 : heap.ValidU64Slice t1 m)
    (hT2 : heap.ValidU64Slice t2 m)
    (hT1T2 : t1.region ≠ t2.region)
    (hT1A : t1.region ≠ A.region) (hT1B : t1.region ≠ B.region)
    (hT2A : t2.region ≠ A.region) (hT2B : t2.region ≠ B.region)
    (hCanonicalA : CanonicalU64Prefix heap A (2 * m) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (2 * m) this._p)
    (hRepA : SlicePolyRep heap A (2 * m) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (2 * m) this._p.toNat right)
    (hrun : karAddHalvesLoop this A B t1 t2 m 0 heap = .ok heap') :
    SlicePolyRep heap' t1 m this._p.toNat (karHalfSumPoly left m) ∧
      SlicePolyRep heap' t2 m this._p.toNat (karHalfSumPoly right m) ∧
      CanonicalU64Prefix heap' t1 m this._p ∧
      CanonicalU64Prefix heap' t2 m this._p := by
  have hvalues := karAddHalvesLoop_coeffs this A B t1 t2 m 0 heap heap'
    left right hp hA hB hT1 hT2 hT1T2 hT1A hT1B hT2A hT2B
    hCanonicalA hCanonicalB hRepA hRepB hrun
  rcases karAddHalvesLoop_ok this A B t1 t2 m 0 heap hA hB hT1 hT2 with
    ⟨okHeap, hok, hlayout⟩
  have hokHeap : okHeap = heap' := Except.ok.inj (hok.symm.trans hrun)
  subst okHeap
  have hT1' := (hlayout t1 m).mp hT1
  have hT2' := (hlayout t2 m).mp hT2
  have hprefix1 : ClassicalCoeffPrefix heap' t1 m (karHalfSumPoly left m) := by
    intro k hk
    rcases hvalues k (Nat.zero_le _) hk with
      ⟨value1, value2, hread1, _, hcoeff1, _, hlt1, _⟩
    refine ⟨value1, hread1, ?_, hlt1⟩
    simpa [coeff_karHalfSumPoly, hk] using hcoeff1
  have hprefix2 : ClassicalCoeffPrefix heap' t2 m (karHalfSumPoly right m) := by
    intro k hk
    rcases hvalues k (Nat.zero_le _) hk with
      ⟨value1, value2, _, hread2, _, hcoeff2, _, hlt2⟩
    refine ⟨value2, hread2, ?_, hlt2⟩
    simpa [coeff_karHalfSumPoly, hk] using hcoeff2
  refine ⟨slicePolyRep_of_classicalCoeffPrefix heap' t1 m
      (karHalfSumPoly left m) hT1' hprefix1 ?_,
    slicePolyRep_of_classicalCoeffPrefix heap' t2 m
      (karHalfSumPoly right m) hT2' hprefix2 ?_,
    canonicalU64Prefix_of_classicalCoeffPrefix heap' t1 m this._p
      (karHalfSumPoly left m) hprefix1,
    canonicalU64Prefix_of_classicalCoeffPrefix heap' t2 m this._p
      (karHalfSumPoly right m) hprefix2⟩
  · intro degree hdegree
    rw [coeff_karHalfSumPoly, if_neg (by omega)]
  · intro degree hdegree
    rw [coeff_karHalfSumPoly, if_neg (by omega)]

theorem karOddTail_refines_slices_odd (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m h : Nat) (heap heap' : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hsucc : h = m + 1)
    (hA : heap.ValidU64Slice A (m + h))
    (hB : heap.ValidU64Slice B (m + h))
    (hT1 : heap.ValidU64Slice t1 h)
    (hT2 : heap.ValidU64Slice t2 h)
    (hT1T2 : t1.region ≠ t2.region)
    (hCanonicalA : CanonicalU64Prefix heap A (m + h) this._p)
    (hCanonicalB : CanonicalU64Prefix heap B (m + h) this._p)
    (hCanonicalT1 : CanonicalU64Prefix heap t1 m this._p)
    (hCanonicalT2 : CanonicalU64Prefix heap t2 m this._p)
    (hRepA : SlicePolyRep heap A (m + h) this._p.toNat left)
    (hRepB : SlicePolyRep heap B (m + h) this._p.toNat right)
    (hRepT1 : SlicePolyRep heap t1 m this._p.toNat (karHalfSumPoly left m))
    (hRepT2 : SlicePolyRep heap t2 m this._p.toNat (karHalfSumPoly right m))
    (hrun : karOddTail A B t1 t2 m h heap = .ok heap') :
    SlicePolyRep heap' t1 h this._p.toNat
        (karHalfSumPoly left m +
          Polynomial.monomial m (left.coeff (m + m))) ∧
      SlicePolyRep heap' t2 h this._p.toNat
        (karHalfSumPoly right m +
          Polynomial.monomial m (right.coeff (m + m))) ∧
      CanonicalU64Prefix heap' t1 h this._p ∧
      CanonicalU64Prefix heap' t2 h this._p := by
  have hodd : h > m := by omega
  rcases karOddTail_coeffs this A B t1 t2 m h heap heap' left right
      hA hB hT1 hT2 hT1T2 hCanonicalA hCanonicalB hRepA hRepB hodd hrun with
    ⟨aTail, bTail, ht1, ht2, hcoeffA, hcoeffB, haLt, hbLt⟩
  rcases karOddTail_preserves_own_prefixes A B t1 t2 m h heap heap'
      hA hB hT1 hT2 hT1T2 hrun with ⟨hsameT1, hsameT2⟩
  rcases karOddTail_ok A B t1 t2 m h heap hA hB hT1 hT2 with
    ⟨okHeap, hok, hlayout⟩
  have hokHeap : okHeap = heap' := Except.ok.inj (hok.symm.trans hrun)
  subst okHeap
  have hT1m := heap.validU64Slice_mono t1 h m hT1 (by omega)
  have hT2m := heap.validU64Slice_mono t2 h m hT2 (by omega)
  have hT1' := (hlayout t1 h).mp hT1
  have hT2' := (hlayout t2 h).mp hT2
  have hT1m' := (hlayout t1 m).mp hT1m
  have hT2m' := (hlayout t2 m).mp hT2m
  have hRepT1' := slicePolyRep_of_same_prefix heap heap' t1 m
    this._p.toNat (karHalfSumPoly left m) hT1m hT1m' hsameT1 hRepT1
  have hRepT2' := slicePolyRep_of_same_prefix heap heap' t2 m
    this._p.toNat (karHalfSumPoly right m) hT2m hT2m' hsameT2 hRepT2
  rcases slicePolyRep_extend_exists heap' t1 m this._p.toNat aTail
      (karHalfSumPoly left m) (by simpa [hsucc] using hT1') hRepT1' ht1 with
    ⟨full1, hfull1, heq1⟩
  rcases slicePolyRep_extend_exists heap' t2 m this._p.toNat bTail
      (karHalfSumPoly right m) (by simpa [hsucc] using hT2') hRepT2' ht2 with
    ⟨full2, hfull2, heq2⟩
  have htarget1 : full1 = karHalfSumPoly left m +
      Polynomial.monomial m (left.coeff (m + m)) := by
    simpa [hcoeffA] using heq1
  have htarget2 : full2 = karHalfSumPoly right m +
      Polynomial.monomial m (right.coeff (m + m)) := by
    simpa [hcoeffB] using heq2
  rw [htarget1] at hfull1
  rw [htarget2] at hfull2
  refine ⟨by simpa [hsucc] using hfull1, by simpa [hsucc] using hfull2, ?_, ?_⟩
  · intro k value hk hread
    by_cases hkm : k < m
    · rcases heap.readU64_of_valid t1 m k hT1m hkm with ⟨old, hold⟩
      have hpreserved := hsameT1 k old hkm hold
      have heq : value = old := Except.ok.inj (hread.symm.trans hpreserved)
      subst value
      exact hCanonicalT1 k old hkm hold
    · have hkEq : k = m := by omega
      subst k
      have heq : value = aTail := Except.ok.inj (hread.symm.trans ht1)
      subst value
      exact haLt
  · intro k value hk hread
    by_cases hkm : k < m
    · rcases heap.readU64_of_valid t2 m k hT2m hkm with ⟨old, hold⟩
      have hpreserved := hsameT2 k old hkm hold
      have heq : value = old := Except.ok.inj (hread.symm.trans hpreserved)
      subst value
      exact hCanonicalT2 k old hkm hold
    · have hkEq : k = m := by omega
      subst k
      have heq : value = bTail := Except.ok.inj (hread.symm.trans ht2)
      subst value
      exact hbLt

theorem normaliseU64_eq_length_of_classicalCoeffPrefix {p : Nat}
    (heap : RawHeap) (C : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod p))
    (hpos : 0 < length)
    (hprefix : ClassicalCoeffPrefix heap C length poly)
    (hlast : poly.coeff (length - 1) ≠ 0) :
    heap.normaliseU64 C length = .ok length := by
  cases length with
  | zero => omega
  | succ n =>
    rcases hprefix n (by omega) with ⟨value, hread, hcoeff, _⟩
    have hvalue : value ≠ 0 := by
      intro hzero
      subst value
      apply hlast
      simpa using hcoeff.symm
    simp [RawHeap.normaliseU64, hread, hvalue]

theorem mul_coeff_zero_of_slice_lengths {p : Nat}
    (heap : RawHeap) (A B : RawPtr UInt64) (lenA lenB degree : Nat)
    (left right : Polynomial (ZMod p))
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hRepA : SlicePolyRep heap A lenA p left)
    (hRepB : SlicePolyRep heap B lenB p right)
    (hdegree : lenA + lenB - 1 ≤ degree) :
    (left * right).coeff degree = 0 := by
  rw [Polynomial.coeff_mul]
  apply Finset.sum_eq_zero
  intro pair hpair
  simp only [Finset.mem_antidiagonal] at hpair
  by_cases hleft : lenA ≤ pair.1
  · rw [slicePolyRep_coeff_zero_of_length_le heap A lenA p left hRepA
      pair.1 hleft, zero_mul]
  · have hright : lenB ≤ pair.2 := by omega
    rw [slicePolyRep_coeff_zero_of_length_le heap B lenB p right hRepB
      pair.2 hright, mul_zero]

theorem mul_last_coeff_ne_zero_of_rawDense (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (heap : RawHeap) (A B : RawPtr UInt64) (lenA lenB : Nat)
    (left right : Polynomial (ZMod this._p.toNat))
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hLeft : RawDensePolyRep this heap A lenA left)
    (hRight : RawDensePolyRep this heap B lenB right) :
    (left * right).coeff (lenA + lenB - 1 - 1) ≠ 0 := by
  rcases hLeft with ⟨hA, hCanonicalA, hRepA, hNormA⟩
  rcases hRight with ⟨hB, hCanonicalB, hRepB, hNormB⟩
  have hleftLast := normaliseU64_poly_last_coeff_ne_zero heap A lenA
    this._p.toNat lenA left hA hRepA hCanonicalA hNormA
    (Nat.ne_of_gt hApos)
  have hrightLast := normaliseU64_poly_last_coeff_ne_zero heap B lenB
    this._p.toNat lenB right hB hRepB hCanonicalB hNormB
    (Nat.ne_of_gt hBpos)
  have hleftDegree := normaliseU64_poly_natDegree_eq heap A lenA
    this._p.toNat lenA left hA hRepA hCanonicalA hNormA
    (Nat.ne_of_gt hApos)
  have hrightDegree := normaliseU64_poly_natDegree_eq heap B lenB
    this._p.toNat lenB right hB hRepB hCanonicalB hNormB
    (Nat.ne_of_gt hBpos)
  have hcoeff := Polynomial.coeff_mul_add_eq_of_natDegree_le
    (show left.natDegree ≤ lenA - 1 by omega)
    (show right.natDegree ≤ lenB - 1 by omega)
  rw [show lenA + lenB - 1 - 1 = (lenA - 1) + (lenB - 1) by omega,
    hcoeff]
  exact mul_ne_zero hleftLast hrightLast

theorem classicalCoeffPrefix_succ_of_write {p : Nat}
    (before after : RawHeap) (C : RawPtr UInt64) (upto : Nat)
    (poly : Polynomial (ZMod p)) (value : UInt64)
    (hprefix : ClassicalCoeffPrefix before C upto poly)
    (hwrite : before.writeU64 C upto value = .ok after)
    (hvalue : (value.toNat : ZMod p) = poly.coeff upto)
    (hvalueLt : value.toNat < p) :
    ClassicalCoeffPrefix after C (upto + 1) poly := by
  intro i hi
  by_cases heq : i = upto
  · subst i
    exact ⟨value, RawHeap.readU64_writeU64_same before after C upto value
      hwrite, hvalue, hvalueLt⟩
  · have hiOld : i < upto := by omega
    rcases hprefix i hiOld with ⟨old, hread, hold, holdLt⟩
    refine ⟨old, ?_, hold, holdLt⟩
    exact RawHeap.readU64_writeU64_ne before after C C upto i value old
      hwrite hread (Or.inr (by omega))

theorem classicalOuterLoop_preserves_coeff_prefix {p : Nat}
    (this : DenseUPolyZp) (C A B : RawPtr UInt64)
    (lenA lenB lenC k : Nat) (heap heap' : RawHeap)
    (poly : Polynomial (ZMod p))
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hprefix : ClassicalCoeffPrefix heap C k poly)
    (hrun : classicalOuterLoop this C A B lenA lenB lenC k heap =
      .ok heap') :
    ClassicalCoeffPrefix heap' C k poly := by
  intro i hi
  rcases hprefix i hi with ⟨old, hread, hold, holdLt⟩
  refine ⟨old, ?_, hold, holdLt⟩
  apply classicalOuterLoop_preserves_outside this C A B C lenA lenB lenC
    k i heap heap' old hC hA hB hread
  · intro target hktarget _
    exact Or.inr (by omega)
  · exact hrun

theorem classicalOuterLoop_refines_coeff_prefix (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (lenA lenB lenC k : Nat) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hLenAWord : lenA < limbBase)
    (hlenC : lenC = lenA + lenB - 1)
    (hC : heap.ValidU64Slice C lenC)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCA : C.region ≠ A.region) (hCB : C.region ≠ B.region)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hRepA : SlicePolyRep heap A lenA this._p.toNat left)
    (hRepB : SlicePolyRep heap B lenB this._p.toNat right)
    (hprefix : ClassicalCoeffPrefix heap C k (left * right)) :
    ∃ heap', classicalOuterLoop this C A B lenA lenB lenC k heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' ∧
      ClassicalCoeffPrefix heap' C lenC (left * right) := by
  unfold classicalOuterLoop
  split
  next hk =>
    let jMin := if k ≥ lenB then k - lenB + 1 else 0
    let jMax := if k < lenA then k else lenA - 1
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
    have hvalueResult : (value.toNat : ZMod this._p.toNat) =
          (left * right).coeff k ∧ value.toNat < this._p.toNat := by
      dsimp [value]
      exact classicalReduced_source_eq_coeff this heap A B lenA lenB lenC k
        left right acc hcfg hp hApos hBpos hLenAWord hlenC hk hA hB
        hCanonicalA hCanonicalB hRepA hRepB (by simpa [jMin, jMax] using hdot)
    have hvalue := hvalueResult.1
    have hvalueLt := hvalueResult.2
    rcases heap.writeU64_of_valid C lenC k value hC hk with ⟨heap1, hw⟩
    simp only [value, hw]
    have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 C k value hw
    have hC1 := (hlayout1 C lenC).mp hC
    have hA1 := (hlayout1 A lenA).mp hA
    have hB1 := (hlayout1 B lenB).mp hB
    have hsameA : SameU64Prefix heap heap1 A lenA := by
      intro i old hi hread
      exact RawHeap.readU64_writeU64_ne heap heap1 C A k i value old hw
        hread (Or.inl hCA)
    have hsameB : SameU64Prefix heap heap1 B lenB := by
      intro i old hi hread
      exact RawHeap.readU64_writeU64_ne heap heap1 C B k i value old hw
        hread (Or.inl hCB)
    have hCanonicalA1 := canonicalU64Prefix_of_same_prefix heap heap1 A lenA
      this._p hA hsameA hCanonicalA
    have hCanonicalB1 := canonicalU64Prefix_of_same_prefix heap heap1 B lenB
      this._p hB hsameB hCanonicalB
    have hRepA1 := slicePolyRep_of_same_prefix heap heap1 A lenA
      this._p.toNat left hA hA1 hsameA hRepA
    have hRepB1 := slicePolyRep_of_same_prefix heap heap1 B lenB
      this._p.toNat right hB hB1 hsameB hRepB
    have hprefix1 := classicalCoeffPrefix_succ_of_write heap heap1 C k
      (left * right) value hprefix hw hvalue hvalueLt
    rcases classicalOuterLoop_refines_coeff_prefix this C A B lenA lenB lenC
      (k + 1) heap1 left right hcfg hp hApos hBpos hLenAWord hlenC hC1
      hA1 hB1 hCA hCB hCanonicalA1 hCanonicalB1 hRepA1 hRepB1 hprefix1 with
      ⟨heap2, hrun, hlayout2, hfull⟩
    rw [hrun]
    exact ⟨heap2, rfl, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length), hfull⟩
  next hnot =>
    refine ⟨heap, rfl, fun _ _ => Iff.rfl, ?_⟩
    intro i hi
    exact hprefix i (by omega)
termination_by lenC - k
decreasing_by omega

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

theorem classicalMul_refines_slice (this : DenseUPolyZp)
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hLenAWord : lenA < limbBase)
    (hC : heap.ValidU64Slice C (lenA + lenB - 1))
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hCA : C.region ≠ A.region) (hCB : C.region ≠ B.region)
    (hCanonicalA : CanonicalU64Prefix heap A lenA this._p)
    (hCanonicalB : CanonicalU64Prefix heap B lenB this._p)
    (hRepA : SlicePolyRep heap A lenA this._p.toNat left)
    (hRepB : SlicePolyRep heap B lenB this._p.toNat right) :
    ∃ heap', dense_upoly_zp__classical_mul_ir this C A lenA B lenB heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' ∧
      SlicePolyRep heap' C (lenA + lenB - 1) this._p.toNat
        (left * right) := by
  have hempty : ClassicalCoeffPrefix heap C 0 (left * right) := by
    intro _ hi
    omega
  rcases classicalOuterLoop_refines_coeff_prefix this C A B lenA lenB
      (lenA + lenB - 1) 0 heap left right hcfg hp hApos hBpos hLenAWord
      rfl hC hA hB hCA hCB hCanonicalA hCanonicalB hRepA hRepB hempty with
    ⟨heap', hrun, hlayout, hprefix⟩
  have hvalid' := (hlayout C (lenA + lenB - 1)).mp hC
  have hslice := slicePolyRep_of_classicalCoeffPrefix heap' C
    (lenA + lenB - 1) (left * right) hvalid' hprefix
    (by
      intro degree hdegree
      exact mul_coeff_zero_of_slice_lengths heap A B lenA lenB degree
        left right hApos hBpos hRepA hRepB hdegree)
  refine ⟨heap', ?_, hlayout, hslice⟩
  simp [dense_upoly_zp__classical_mul_ir, Nat.ne_of_gt hApos,
    Nat.ne_of_gt hBpos, hrun]

theorem classicalMul_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (hApos : 0 < lenA) (hBpos : 0 < lenB)
    (hLenAWord : lenA < limbBase)
    (hC : heap.ValidU64Slice C (lenA + lenB - 1))
    (hCA : C.region ≠ A.region) (hCB : C.region ≠ B.region)
    (hLeft : RawDensePolyRep this heap A lenA left)
    (hRight : RawDensePolyRep this heap B lenB right) :
    ∃ heap', dense_upoly_zp__classical_mul_ir this C A lenA B lenB heap =
        .ok heap' ∧ RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' C (lenA + lenB - 1) (left * right) := by
  have hempty : ClassicalCoeffPrefix heap C 0 (left * right) := by
    intro _ hi
    omega
  rcases classicalOuterLoop_refines_coeff_prefix this C A B lenA lenB
      (lenA + lenB - 1) 0 heap left right hcfg hp hApos hBpos hLenAWord
      rfl hC hLeft.1 hRight.1 hCA hCB hLeft.2.1 hRight.2.1
      hLeft.2.2.1 hRight.2.2.1 hempty with
    ⟨heap', hrun, hlayout, hprefix⟩
  have hvalid' := (hlayout C (lenA + lenB - 1)).mp hC
  have hcanonical' := canonicalU64Prefix_of_classicalCoeffPrefix heap' C
    (lenA + lenB - 1) this._p (left * right) hprefix
  have hslice' := slicePolyRep_of_classicalCoeffPrefix heap' C
    (lenA + lenB - 1) (left * right) hvalid' hprefix
    (by
      intro degree hdegree
      exact mul_coeff_zero_of_slice_lengths heap A B lenA lenB degree
        left right hApos hBpos hLeft.2.2.1 hRight.2.2.1 hdegree)
  have hlast := mul_last_coeff_ne_zero_of_rawDense this heap A B lenA lenB
    left right hApos hBpos hLeft hRight
  have hnorm' := normaliseU64_eq_length_of_classicalCoeffPrefix heap' C
    (lenA + lenB - 1) (left * right) (by omega) hprefix hlast
  refine ⟨heap', ?_, hlayout, hvalid', hcanonical', hslice', hnorm'⟩
  simp [dense_upoly_zp__classical_mul_ir, Nat.ne_of_gt hApos,
    Nat.ne_of_gt hBpos, hrun]

end CLPoly.Impl.StrictMulRefinement
