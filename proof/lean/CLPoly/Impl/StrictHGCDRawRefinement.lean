import CLPoly.Generated.StrictHGCD
import CLPoly.Impl.StrictMulDispatchRefinement
import CLPoly.Impl.StrictPolyAddSubRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictHGCDRawRefinement

open Generated.StrictHGCD
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictPolyAddSubRefinement
open CLPoly.Impl.StrictMulRefinement
open CLPoly.Impl.StrictWordArithmetic

/-- Bounds-safe accessors justified by the C array layout invariant. -/
def hgcdMatPtr (M : HgcdMat) (hM : M.Valid) (i : Fin 4) : RawPtr UInt64 :=
  M.poly[i.val]'(by rw [hM.1]; exact i.isLt)

def hgcdMatLen (M : HgcdMat) (hM : M.Valid) (i : Fin 4) : Nat :=
  M.len[i.val]'(by rw [hM.2]; exact i.isLt)

/-- Raw representation of the four polynomial entries of an HGCD matrix. -/
def HgcdMatPolyRep (heap : RawHeap) (M : HgcdMat) (p : Nat)
    (entries : Fin 4 → Polynomial (ZMod p)) (hM : M.Valid) : Prop :=
  ∀ i : Fin 4, SlicePolyRep heap (hgcdMatPtr M hM i)
    (hgcdMatLen M hM i) p (entries i)

/-- Pointwise heap-frame transport for all four matrix entries. -/
theorem hgcdMatPolyRep_of_same_prefixes (before after : RawHeap)
    (M : HgcdMat) (p : Nat) (entries : Fin 4 → Polynomial (ZMod p))
    (hM : M.Valid) (hlayout : RawHeap.SameLayout before after)
    (hvalid : ∀ i : Fin 4, before.ValidU64Slice
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hsame : ∀ i : Fin 4, SameU64Prefix before after
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hrep : HgcdMatPolyRep before M p entries hM) :
    HgcdMatPolyRep after M p entries hM := by
  intro i
  exact slicePolyRep_of_same_prefix before after (hgcdMatPtr M hM i)
    (hgcdMatLen M hM i) p (entries i) (hvalid i)
    ((hlayout (hgcdMatPtr M hM i) (hgcdMatLen M hM i)).mp (hvalid i))
    (hsame i) (hrep i)

/-- A successful generated recursive memcpy preserves a represented HGCD
matrix whenever its destination is disjoint from every live matrix slice. -/
theorem copyU64_preserves_hgcdMatPolyRep (heap heap' : RawHeap)
    (dst src : RawPtr UInt64) (count : Nat) (M : HgcdMat) (p : Nat)
    (entries : Fin 4 → Polynomial (ZMod p)) (hM : M.Valid)
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hValidMatrix : ∀ i : Fin 4, heap.ValidU64Slice
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hMatrix : ∀ i : Fin 4, U64SlicesDisjoint dst count
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hcopy : heap.copyU64 dst src count = .ok heap')
    (hrep : HgcdMatPolyRep heap M p entries hM) :
    RawHeap.SameLayout heap heap' ∧
      HgcdMatPolyRep heap' M p entries hM := by
  rcases copyU64_ok heap dst src count hDst hSrc with
    ⟨copyHeap, hcopy', hlayout⟩
  have heq : copyHeap = heap' := Except.ok.inj (hcopy'.symm.trans hcopy)
  subst copyHeap
  refine ⟨hlayout, hgcdMatPolyRep_of_same_prefixes heap heap' M p entries
    hM hlayout hValidMatrix ?_ hrep⟩
  intro i
  exact copyU64_preserves_prefix heap heap' dst src (hgcdMatPtr M hM i)
    count (hgcdMatLen M hM i) hDst hSrc (hMatrix i) hcopy

noncomputable def identityEntries (p : Nat) : Fin 4 → Polynomial (ZMod p)
  | ⟨0, _⟩ => 1
  | ⟨1, _⟩ => 0
  | ⟨2, _⟩ => 0
  | ⟨3, _⟩ => 1

def identityEntryLen : Fin 4 → Nat
  | ⟨0, _⟩ => 1
  | ⟨1, _⟩ => 0
  | ⟨2, _⟩ => 0
  | ⟨3, _⟩ => 1

/-- A normalized raw polynomial survives any execution that preserves both
its allocation layout and every cell in its declared prefix. -/
theorem rawDensePolyRep_of_same_prefix (this : DenseUPolyZp)
    (before after : RawHeap) (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hlayout : RawHeap.SameLayout before after)
    (hsame : SameU64Prefix before after ptr length)
    (hrep : RawDensePolyRep this before ptr length poly) :
    RawDensePolyRep this after ptr length poly := by
  have hvalidAfter := (hlayout ptr length).mp hrep.1
  refine ⟨hvalidAfter,
    canonicalU64Prefix_of_same_prefix before after ptr length this._p
      hrep.1 hsame hrep.2.1,
    slicePolyRep_of_same_prefix before after ptr length this._p.toNat poly
      hrep.1 hvalidAfter hsame hrep.2.2.1, ?_⟩
  have hnorm := normaliseU64_eq_of_prefix_map before after ptr ptr length
    hrep.1 hsame
  exact hnorm.symm.trans hrep.2.2.2

/-- The actual recursive RawHeap memcpy transports the complete normalized
polynomial representation to a disjoint destination. -/
theorem copyU64_refines_rawDense (this : DenseUPolyZp)
    (heap : RawHeap) (dst src : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hDst : heap.ValidU64Slice dst length)
    (hdisjoint : U64SlicesDisjoint dst length src length)
    (hrep : RawDensePolyRep this heap src length poly) :
    ∃ heap', heap.copyU64 dst src length = .ok heap' ∧
      RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' dst length poly := by
  rcases copyU64_refines_disjoint heap dst src length hDst hrep.1 hdisjoint
      with ⟨heap', hcopy, hlayout, hcontents⟩
  rcases copyU64_refines_slice_canonical heap dst src length
      this._p.toNat poly this._p hDst hrep.1 hdisjoint hrep.2.2.1
      hrep.2.1 with ⟨repHeap, hcopyRep, _, hslice, hcanonical⟩
  have heq : repHeap = heap' := Except.ok.inj (hcopyRep.symm.trans hcopy)
  subst repHeap
  have hnorm := normaliseU64_eq_of_prefix_map heap heap' src dst length
    hrep.1 hcontents
  exact ⟨heap', hcopy, hlayout, (hlayout dst length).mp hDst, hcanonical,
    hslice, hnorm.symm.trans hrep.2.2.2⟩

theorem slicePolyRep_zero_length_any (heap : RawHeap) (ptr : RawPtr UInt64)
    (p : Nat) : SlicePolyRep heap ptr 0 p 0 := by
  refine ⟨#[], rfl, rfl, ?_⟩
  ext degree
  rw [coeff_coeffArrayPoly]
  simp

theorem slicePolyRep_one_of_read_one (heap : RawHeap) (ptr : RawPtr UInt64)
    (p : Nat) (hvalid : heap.ValidU64Slice ptr 1)
    (hread : heap.readU64 ptr 0 = .ok 1) :
    SlicePolyRep heap ptr 1 p 1 := by
  rcases slicePolyRep_exists_unique heap ptr 1 p hvalid with
    ⟨poly, hrep, _⟩
  have heq : poly = 1 := by
    ext degree
    by_cases hd : degree = 0
    · subst degree
      rcases slicePolyRep_coeff heap ptr 1 p poly hrep 0 (by omega) with
        ⟨value, hvalue, hcoeff⟩
      have hv : value = 1 := Except.ok.inj (hvalue.symm.trans hread)
      subst value
      simpa using hcoeff
    · rw [slicePolyRep_coeff_zero_of_length_le heap ptr 1 p poly hrep
          degree (by omega)]
      simp [Polynomial.coeff_one, hd]
  simpa [heq] using hrep

theorem matOne_refines (M : HgcdMat) (heap : RawHeap) (p : Nat)
    (hM : M.Valid)
    (h0 : heap.ValidU64Slice (hgcdMatPtr M hM ⟨0, by omega⟩) 1)
    (h3 : heap.ValidU64Slice (hgcdMatPtr M hM ⟨3, by omega⟩) 1)
    (h03 : CLPoly.Impl.StrictMulRefinement.U64SlicesDisjoint
      (hgcdMatPtr M hM ⟨0, by omega⟩) 1
      (hgcdMatPtr M hM ⟨3, by omega⟩) 1) :
    ∃ heap' M', dense_upoly_zp__mat_one_ir M heap = .ok (heap', M') ∧
      RawHeap.SameLayout heap heap' ∧
      M'.poly = M.poly ∧ M'.len = #[1, 0, 0, 1] ∧
      ∃ hM' : M'.Valid,
        HgcdMatPolyRep heap' M' p (identityEntries p) hM' := by
  have hvalid : M.poly.size = 4 ∧ M.len.size = 4 := by
    simpa [HgcdMat.Valid] using hM
  let p0 := hgcdMatPtr M hM ⟨0, by omega⟩
  let p3 := hgcdMatPtr M hM ⟨3, by omega⟩
  rcases heap.writeU64_of_valid p0 1 0 1 (by simpa [p0] using h0)
      (by omega) with ⟨heap1, hwrite0⟩
  have hlayout1 := RawHeap.writeU64_sameLayout heap heap1 p0 0 1 hwrite0
  have h3one := (hlayout1 p3 1).mp (by simpa [p3] using h3)
  rcases heap1.writeU64_of_valid p3 1 0 1 h3one (by omega) with
    ⟨heap2, hwrite3⟩
  have hlayout2 := RawHeap.writeU64_sameLayout heap1 heap2 p3 0 1 hwrite3
  let M' : HgcdMat := { M with len := #[1, 0, 0, 1] }
  have hwrite0' : heap.writeU64 (M.poly[0]'(by omega)) 0 1 = .ok heap1 := by
    simpa [p0, hgcdMatPtr] using hwrite0
  have hwrite3' : heap1.writeU64 (M.poly[3]'(by omega)) 0 1 = .ok heap2 := by
    simpa [p3, hgcdMatPtr] using hwrite3
  have hrun : dense_upoly_zp__mat_one_ir M heap = .ok (heap2, M') := by
    simp [dense_upoly_zp__mat_one_ir, hvalid, hwrite0', hwrite3', M']
  refine ⟨heap2, M', hrun, fun ptr length =>
    (hlayout1 ptr length).trans (hlayout2 ptr length), rfl, rfl, ?_⟩
  have hsame0 := CLPoly.Impl.StrictMulRefinement.writeU64_preserves_prefix
    heap1 heap2 p3 p0 1 1 0 1
    (by simpa [p0, p3] using
      CLPoly.Impl.StrictMulRefinement.u64SlicesDisjoint_symm h03)
    (by omega) hwrite3
  have hread0Heap1 := RawHeap.readU64_writeU64_same heap heap1 p0 0 1 hwrite0
  have hread0 : heap2.readU64 p0 0 = .ok 1 :=
    hsame0 0 1 (by omega) hread0Heap1
  have hread3 : heap2.readU64 p3 0 = .ok 1 :=
    RawHeap.readU64_writeU64_same heap1 heap2 p3 0 1 hwrite3
  have hM' : M'.Valid := by simp [M', HgcdMat.Valid, hvalid.1]
  refine ⟨hM', ?_⟩
  intro i
  fin_cases i
  · simpa [HgcdMatPolyRep, hgcdMatPtr, hgcdMatLen, M', identityEntries,
      p0] using
      slicePolyRep_one_of_read_one heap2 p0 p
        ((hlayout2 p0 1).mp ((hlayout1 p0 1).mp (by simpa [p0] using h0)))
        hread0
  · simpa [HgcdMatPolyRep, hgcdMatPtr, hgcdMatLen, M', identityEntries] using
      slicePolyRep_zero_length_any heap2
        (hgcdMatPtr M hM ⟨1, by omega⟩) p
  · simpa [HgcdMatPolyRep, hgcdMatPtr, hgcdMatLen, M', identityEntries] using
      slicePolyRep_zero_length_any heap2
        (hgcdMatPtr M hM ⟨2, by omega⟩) p
  · simpa [HgcdMatPolyRep, hgcdMatPtr, hgcdMatLen, M', identityEntries,
      p3] using
      slicePolyRep_one_of_read_one heap2 p3 p ((hlayout2 p3 1).mp h3one)
        hread3

/-- The two real writes in `_mat_one` preserve every declared prefix outside
the two matrix-entry cells. -/
theorem matOne_preserves_prefix (M : HgcdMat) (heap heap' : RawHeap)
    (M' : HgcdMat) (guard : RawPtr UInt64) (guardLen : Nat)
    (hM : M.Valid)
    (h0Guard : U64SlicesDisjoint
      (hgcdMatPtr M hM ⟨0, by omega⟩) 1 guard guardLen)
    (h3Guard : U64SlicesDisjoint
      (hgcdMatPtr M hM ⟨3, by omega⟩) 1 guard guardLen)
    (hrun : dense_upoly_zp__mat_one_ir M heap = .ok (heap', M')) :
    SameU64Prefix heap heap' guard guardLen := by
  have hvalid : M.poly.size = 4 ∧ M.len.size = 4 := by
    simpa [HgcdMat.Valid] using hM
  let p0 := hgcdMatPtr M hM ⟨0, by omega⟩
  let p3 := hgcdMatPtr M hM ⟨3, by omega⟩
  cases hwrite0 : heap.writeU64 p0 0 1 with
  | error fault =>
      have hwrite0' : heap.writeU64 (M.poly[0]'(by omega)) 0 1 =
          .error fault := by simpa [p0, hgcdMatPtr] using hwrite0
      simp [dense_upoly_zp__mat_one_ir, hvalid, hwrite0'] at hrun
  | ok heap1 =>
      have hwrite0' : heap.writeU64 (M.poly[0]'(by omega)) 0 1 =
          .ok heap1 := by simpa [p0, hgcdMatPtr] using hwrite0
      cases hwrite3 : heap1.writeU64 p3 0 1 with
      | error fault =>
          have hwrite3' : heap1.writeU64 (M.poly[3]'(by omega)) 0 1 =
              .error fault := by simpa [p3, hgcdMatPtr] using hwrite3
          simp [dense_upoly_zp__mat_one_ir, hvalid, hwrite0', hwrite3'] at hrun
      | ok heap2 =>
          have hwrite3' : heap1.writeU64 (M.poly[3]'(by omega)) 0 1 =
              .ok heap2 := by simpa [p3, hgcdMatPtr] using hwrite3
          have heq : heap' = heap2 := by
            have hrun' : (.ok (heap2, { M with len := #[1, 0, 0, 1] }) :
                RawExec (RawHeap × HgcdMat)) = .ok (heap', M') := by
              simpa [dense_upoly_zp__mat_one_ir, hvalid, hwrite0', hwrite3']
                using hrun
            exact (congrArg Prod.fst (Except.ok.inj hrun')).symm
          subst heap'
          have hsame0 := writeU64_preserves_prefix heap heap1 p0 guard 1
            guardLen 0 1 h0Guard (by omega) hwrite0
          have hsame3 := writeU64_preserves_prefix heap1 heap2 p3 guard 1
            guardLen 0 1 h3Guard (by omega) hwrite3
          exact sameU64Prefix_trans hsame0 hsame3

/-- End-to-end refinement of the exact initialization prefix of C++
`_hgcd_iter`: identity matrix construction followed by the ordered A and B
copies.  Every alias restriction below is a physical L1 memory condition. -/
theorem hgcdIterInit_refines (this : DenseUPolyZp)
    (M : HgcdMat) (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (left right : Polynomial (ZMod this._p.toNat))
    (hM : M.Valid)
    (h0 : heap.ValidU64Slice (hgcdMatPtr M hM ⟨0, by omega⟩) 1)
    (h3 : heap.ValidU64Slice (hgcdMatPtr M hM ⟨3, by omega⟩) 1)
    (h03 : U64SlicesDisjoint (hgcdMatPtr M hM ⟨0, by omega⟩) 1
      (hgcdMatPtr M hM ⟨3, by omega⟩) 1)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAa : U64SlicesDisjoint A lenA a lenA)
    (hBb : U64SlicesDisjoint B lenB b lenB)
    (hAb : U64SlicesDisjoint A lenA b lenB)
    (hBA : U64SlicesDisjoint B lenB A lenA)
    (h0a : U64SlicesDisjoint (hgcdMatPtr M hM ⟨0, by omega⟩) 1 a lenA)
    (h3a : U64SlicesDisjoint (hgcdMatPtr M hM ⟨3, by omega⟩) 1 a lenA)
    (h0b : U64SlicesDisjoint (hgcdMatPtr M hM ⟨0, by omega⟩) 1 b lenB)
    (h3b : U64SlicesDisjoint (hgcdMatPtr M hM ⟨3, by omega⟩) 1 b lenB)
    (hAMatrix : ∀ i : Fin 4, U64SlicesDisjoint A lenA
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hBMatrix : ∀ i : Fin 4, U64SlicesDisjoint B lenB
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hMatrixValid : ∀ i : Fin 4, heap.ValidU64Slice
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hLeft : RawDensePolyRep this heap a lenA left)
    (hRight : RawDensePolyRep this heap b lenB right) :
    ∃ initial, hgcdIterInit M A B T t lenT a lenA b lenB heap =
        .ok initial ∧
      initial.A = A ∧ initial.lenA = lenA ∧
      initial.B = B ∧ initial.lenB = lenB ∧
      initial.T = T ∧ initial.lenT = lenT ∧ initial.t = t ∧
      initial.sgn = 1 ∧
      ∃ hInitialM : initial.matrix.Valid,
        HgcdMatPolyRep initial.heap initial.matrix this._p.toNat
          (identityEntries this._p.toNat) hInitialM ∧
        RawDensePolyRep this initial.heap initial.A initial.lenA left ∧
        RawDensePolyRep this initial.heap initial.B initial.lenB right := by
  rcases matOne_refines M heap this._p.toNat hM h0 h3 h03 with
    ⟨heap1, M1, hone, hlayout1, hpoly1, hlen1, hM1, hMatrix1⟩
  have hsameA := matOne_preserves_prefix M heap heap1 M1 a lenA hM
    h0a h3a hone
  have hsameB := matOne_preserves_prefix M heap heap1 M1 b lenB hM
    h0b h3b hone
  have hLeft1 := rawDensePolyRep_of_same_prefix this heap heap1 a lenA left
    hlayout1 hsameA hLeft
  have hRight1 := rawDensePolyRep_of_same_prefix this heap heap1 b lenB right
    hlayout1 hsameB hRight
  have hPtr1 : ∀ i : Fin 4,
      hgcdMatPtr M1 hM1 i = hgcdMatPtr M hM i := by
    intro i
    fin_cases i <;> simp [hgcdMatPtr, hpoly1]
  have hLen1 : ∀ i : Fin 4, hgcdMatLen M1 hM1 i = identityEntryLen i := by
    intro i
    fin_cases i <;> simp [hgcdMatLen, hlen1, identityEntryLen]
  have hValidMatrix1 : ∀ i : Fin 4, heap1.ValidU64Slice
      (hgcdMatPtr M1 hM1 i) (hgcdMatLen M1 hM1 i) := by
    intro i
    rw [hPtr1 i, hLen1 i]
    exact (hlayout1 _ _).mp (hMatrixValid i)
  have hA1 := (hlayout1 A lenA).mp hA
  rcases copyU64_refines_rawDense this heap1 A a lenA left hA1 hAa
      hLeft1 with ⟨heap2, hcopyA, hlayout2, hA2⟩
  have hsameB2 := copyU64_preserves_prefix heap1 heap2 A a b lenA lenB
    hA1 hLeft1.1 hAb hcopyA
  have hRight2 := rawDensePolyRep_of_same_prefix this heap1 heap2 b lenB
    right hlayout2 hsameB2 hRight1
  have hAMatrix1 : ∀ i : Fin 4, U64SlicesDisjoint A lenA
      (hgcdMatPtr M1 hM1 i) (hgcdMatLen M1 hM1 i) := by
    intro i
    simpa [hPtr1 i, hLen1 i] using hAMatrix i
  rcases copyU64_preserves_hgcdMatPolyRep heap1 heap2 A a lenA M1
      this._p.toNat (identityEntries this._p.toNat) hM1 hA1 hLeft1.1
      hValidMatrix1 hAMatrix1 hcopyA hMatrix1 with
    ⟨_, hMatrix2⟩
  have hB2 := (hlayout2 B lenB).mp ((hlayout1 B lenB).mp hB)
  rcases copyU64_refines_rawDense this heap2 B b lenB right hB2 hBb
      hRight2 with ⟨heap3, hcopyB, hlayout3, hB3⟩
  have hsameA3 := copyU64_preserves_prefix heap2 heap3 B b A lenB lenA
    hB2 hRight2.1 hBA hcopyB
  have hA3 := rawDensePolyRep_of_same_prefix this heap2 heap3 A lenA left
    hlayout3 hsameA3 hA2
  have hValidMatrix2 : ∀ i : Fin 4, heap2.ValidU64Slice
      (hgcdMatPtr M1 hM1 i) (hgcdMatLen M1 hM1 i) := by
    intro i
    exact (hlayout2 _ _).mp (hValidMatrix1 i)
  have hBMatrix1 : ∀ i : Fin 4, U64SlicesDisjoint B lenB
      (hgcdMatPtr M1 hM1 i) (hgcdMatLen M1 hM1 i) := by
    intro i
    simpa [hPtr1 i, hLen1 i] using hBMatrix i
  rcases copyU64_preserves_hgcdMatPolyRep heap2 heap3 B b lenB M1
      this._p.toNat (identityEntries this._p.toNat) hM1 hB2 hRight2.1
      hValidMatrix2 hBMatrix1 hcopyB hMatrix2 with
    ⟨_, hMatrix3⟩
  let initial : HgcdIterState := {
    heap := heap3, matrix := M1, A := A, lenA := lenA, B := B,
    lenB := lenB, T := T, lenT := lenT, t := t, sgn := 1 }
  refine ⟨initial, ?_, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl,
    hM1, ?_, ?_, ?_⟩
  · simp [hgcdIterInit, hone, hcopyA, hcopyB, initial]
  · simpa [initial] using hMatrix3
  · simpa [initial] using hA3
  · simpa [initial] using hB3

/-- The source's zero-quotient/zero-entry branch performs exactly the two
matrix-entry swaps and no heap access.  This exposes the real descriptor
state consumed by the next HGCD iteration. -/
theorem matRowUpdate_zero_exec (this : DenseUPolyZp) (M : HgcdMat)
    (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (hM : M.Valid)
    (hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0) :
    ∃ result,
      dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T lenT t
          scratch heap = .ok result ∧
      result.heap = heap ∧ result.T = T ∧ result.lenT = lenT ∧
      result.t = t ∧ result.matrix.Valid := by
  have hvalid : M.poly.size = 4 ∧ M.len.size = 4 := by
    simpa [HgcdMat.Valid] using hM
  let p0 := hgcdMatPtr M hM i0
  let p1 := hgcdMatPtr M hM i1
  let l0 := hgcdMatLen M hM i0
  let l1 := hgcdMatLen M hM i1
  let poly' := (M.poly.set i1.val p0 (by omega)).set i0.val p1 (by simp; omega)
  let len' := (M.len.set i1.val l0 (by omega)).set i0.val l1 (by simp; omega)
  let result := MatRowUpdateResult.mk heap
    ({ poly := poly', len := len' } : HgcdMat) T lenT t
  have hbranch : ¬(lenQ ≠ 0 ∧ l0 ≠ 0) := by
    intro h
    rcases hzero with hq | hentry
    · exact h.1 hq
    · exact h.2 (by simpa [l0, hgcdMatLen] using hentry)
  have hbranch' : ¬(lenQ ≠ 0 ∧ M.len[i0.val]'(by omega) ≠ 0) := by
    simpa [l0, hgcdMatLen] using hbranch
  have hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result := by
    simp [dense_upoly_zp__mat_row_update_ir, hvalid, hbranch', result,
      poly', len', p0, p1, l0, l1, hgcdMatPtr, hgcdMatLen]
  refine ⟨result, hrun, rfl, rfl, rfl, rfl, ?_⟩
  simp [result, HgcdMat.Valid, poly', len', hvalid]

/-- Successful execution of the nonzero source branch can only arise from
the actual `_mul` followed by the actual `_poly_add`; the returned pointer
state is exactly the three source swaps. -/
theorem matRowUpdate_nonzero_success_shape (this : DenseUPolyZp)
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (hne : i0 ≠ i1) (hQ : lenQ ≠ 0)
    (hEntry : hgcdMatLen M hM i0 ≠ 0)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    ∃ heap1 heap2 sumLen,
      Generated.StrictMul.dense_upoly_zp__mul_ir this T
        (if lenQ ≥ hgcdMatLen M hM i0 then Q else hgcdMatPtr M hM i0)
        (if lenQ ≥ hgcdMatLen M hM i0 then lenQ else hgcdMatLen M hM i0)
        (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatPtr M hM i0 else Q)
        (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatLen M hM i0 else lenQ)
        scratch heap = .ok heap1 ∧
      Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir this t
        (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1) T
        (lenQ + hgcdMatLen M hM i0 - 1) heap1 = .ok (heap2, sumLen) ∧
      result.heap = heap2 ∧ result.T = T ∧
      result.lenT = lenQ + hgcdMatLen M hM i0 - 1 ∧
      result.t = hgcdMatPtr M hM i1 ∧
      ∃ hResult : result.matrix.Valid,
        hgcdMatPtr result.matrix hResult i0 = t ∧
        hgcdMatLen result.matrix hResult i0 = sumLen ∧
        hgcdMatPtr result.matrix hResult i1 = hgcdMatPtr M hM i0 ∧
        hgcdMatLen result.matrix hResult i1 = hgcdMatLen M hM i0 := by
  have hvalid : M.poly.size = 4 ∧ M.len.size = 4 := by
    simpa [HgcdMat.Valid] using hM
  have hneVal : i0.val ≠ i1.val := by
    intro heq
    apply hne
    exact Fin.ext heq
  let p0 := hgcdMatPtr M hM i0
  let p1 := hgcdMatPtr M hM i1
  let l0 := hgcdMatLen M hM i0
  let l1 := hgcdMatLen M hM i1
  let left := if lenQ ≥ l0 then Q else p0
  let leftLen := if lenQ ≥ l0 then lenQ else l0
  let right := if lenQ ≥ l0 then p0 else Q
  let rightLen := if lenQ ≥ l0 then l0 else lenQ
  have hEntry' : M.len[i0.val]'(by omega) ≠ 0 := by
    simpa [hgcdMatLen] using hEntry
  cases hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this T left
      leftLen right rightLen scratch heap with
  | error fault =>
      have hmul' : Generated.StrictMul.dense_upoly_zp__mul_ir this T
          (if lenQ ≥ M.len[i0.val]'(by omega) then Q else M.poly[i0.val]'(by omega))
          (if lenQ ≥ M.len[i0.val]'(by omega) then lenQ else M.len[i0.val]'(by omega))
          (if lenQ ≥ M.len[i0.val]'(by omega) then M.poly[i0.val]'(by omega) else Q)
          (if lenQ ≥ M.len[i0.val]'(by omega) then M.len[i0.val]'(by omega) else lenQ)
          scratch heap = .error fault := by
        simpa [left, leftLen, right, rightLen, p0, l0, hgcdMatPtr,
          hgcdMatLen] using hmul
      simp [dense_upoly_zp__mat_row_update_ir, hvalid, hQ, hEntry', hmul'] at hrun
  | ok heap1 =>
      have hmul' : Generated.StrictMul.dense_upoly_zp__mul_ir this T
          (if lenQ ≥ M.len[i0.val]'(by omega) then Q else M.poly[i0.val]'(by omega))
          (if lenQ ≥ M.len[i0.val]'(by omega) then lenQ else M.len[i0.val]'(by omega))
          (if lenQ ≥ M.len[i0.val]'(by omega) then M.poly[i0.val]'(by omega) else Q)
          (if lenQ ≥ M.len[i0.val]'(by omega) then M.len[i0.val]'(by omega) else lenQ)
          scratch heap = .ok heap1 := by
        simpa [left, leftLen, right, rightLen, p0, l0, hgcdMatPtr,
          hgcdMatLen] using hmul
      cases hadd : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir
          this t p1 l1 T (lenQ + l0 - 1) heap1 with
      | error fault =>
          have hadd' : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir
              this t (M.poly[i1.val]'(by omega)) (M.len[i1.val]'(by omega)) T
              (lenQ + M.len[i0.val]'(by omega) - 1) heap1 = .error fault := by
            simpa [p1, l0, l1, hgcdMatPtr, hgcdMatLen] using hadd
          simp [dense_upoly_zp__mat_row_update_ir, hvalid, hQ, hEntry',
            hmul', hadd'] at hrun
      | ok pair =>
          rcases pair with ⟨heap2, sumLen⟩
          have hadd' : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir
              this t (M.poly[i1.val]'(by omega)) (M.len[i1.val]'(by omega)) T
              (lenQ + M.len[i0.val]'(by omega) - 1) heap1 =
                .ok (heap2, sumLen) := by
            simpa [p1, l0, l1, hgcdMatPtr, hgcdMatLen] using hadd
          have heq : result = MatRowUpdateResult.mk heap2
              ({
                poly := (M.poly.set i1.val p0 (by omega)).set i0.val t
                  (by simp; omega)
                len := (M.len.set i1.val l0 (by omega)).set i0.val sumLen
                  (by simp; omega)
              } : HgcdMat) T (lenQ + l0 - 1) p1 := by
            have hrun' := hrun
            simp [dense_upoly_zp__mat_row_update_ir, hvalid, hQ, hEntry',
              hmul', hadd'] at hrun'
            simpa [p0, p1, l0, hgcdMatPtr, hgcdMatLen] using hrun'.symm
          subst result
          refine ⟨heap1, heap2, sumLen, ?_, ?_, rfl, rfl, ?_, ?_, ?_⟩
          · simpa [left, leftLen, right, rightLen, p0, l0,
              hgcdMatPtr, hgcdMatLen] using hmul
          · simpa [p1, l0, l1, hgcdMatPtr, hgcdMatLen] using hadd
          · simp [l0, hgcdMatLen]
          · simp [p1, hgcdMatPtr]
          · have hResult :
                ({
                  poly := (M.poly.set i1.val p0 (by omega)).set i0.val t
                    (by simp; omega)
                  len := (M.len.set i1.val l0 (by omega)).set i0.val sumLen
                    (by simp; omega)
                } : HgcdMat).Valid := by
              simp [HgcdMat.Valid, hvalid]
            refine ⟨hResult, ?_, ?_, ?_, ?_⟩
            · simp [hgcdMatPtr]
            · simp [hgcdMatLen]
            · simp only [hgcdMatPtr]
              rw [Array.getElem_set_ne
                (by simpa using (show i0.val < M.poly.size by
                  rw [hvalid.1]; exact i0.isLt))
                (by simpa using (show i1.val < M.poly.size by
                  rw [hvalid.1]; exact i1.isLt)) hneVal,
                Array.getElem_set_self]
              rfl
            · simp only [hgcdMatLen]
              rw [Array.getElem_set_ne
                (by simpa using (show i0.val < M.len.size by
                  rw [hvalid.2]; exact i0.isLt))
                (by simpa using (show i1.val < M.len.size by
                  rw [hvalid.2]; exact i1.isLt)) hneVal,
                Array.getElem_set_self]
              rfl

/-- The multiplication call exposed by the real nonzero row-update branch
computes exactly the quotient times the old `i0` entry.  The conditional
operand order is the one used by C++ to satisfy `_mul`'s length contract;
commutativity is used only after refining that actual dispatcher call. -/
theorem matRowUpdate_mul_result (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (hM : M.Valid) (i0 : Fin 4)
    (Q T scratch : RawPtr UInt64) (lenQ : Nat) (heap heap1 : RawHeap)
    (quotient entry0 : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hQpos : 0 < lenQ) (hEntryPos : 0 < hgcdMatLen M hM i0)
    (hLenWord : max lenQ (hgcdMatLen M hM i0) < limbBase)
    (hT : heap.ValidU64Slice T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1))
    (hScratch : heap.ValidU64Slice scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hTQ : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) Q lenQ)
    (hTEntry : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hTScratch : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hScratchQ : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) Q lenQ)
    (hScratchEntry : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntryRep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this T
      (if lenQ ≥ hgcdMatLen M hM i0 then Q else hgcdMatPtr M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then lenQ else hgcdMatLen M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatPtr M hM i0 else Q)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatLen M hM i0 else lenQ)
      scratch heap = .ok heap1) :
    RawHeap.SameLayout heap heap1 ∧
      RawDensePolyRep this heap1 T
        (lenQ + hgcdMatLen M hM i0 - 1) (quotient * entry0) := by
  by_cases horder : lenQ ≥ hgcdMatLen M hM i0
  · rcases mul_refines_rawDense this T Q lenQ (hgcdMatPtr M hM i0)
        (hgcdMatLen M hM i0) scratch heap quotient entry0 hcfg hp hQpos
        hEntryPos horder
        (by simpa [max_eq_left horder] using hLenWord)
        (by simpa [max_eq_left horder] using hT)
        (by simpa [max_eq_left horder] using hScratch)
        (by simpa [max_eq_left horder] using hTQ)
        (by simpa [max_eq_left horder] using hTEntry)
        (by simpa [max_eq_left horder] using hTScratch)
        (by simpa [max_eq_left horder] using hScratchQ)
        (by simpa [max_eq_left horder] using hScratchEntry)
        hQRep hEntryRep with ⟨heap', hrun, hlayout, hrep⟩
    have heq : heap' = heap1 := Except.ok.inj (hrun.symm.trans (by
      simpa [horder] using hmul))
    subst heap'
    exact ⟨hlayout, hrep⟩
  · have hle : lenQ ≤ hgcdMatLen M hM i0 := by omega
    rcases mul_refines_rawDense this T (hgcdMatPtr M hM i0)
        (hgcdMatLen M hM i0) Q lenQ scratch heap entry0 quotient hcfg hp
        hEntryPos hQpos hle
        (by simpa [max_eq_right hle] using hLenWord)
        (by simpa [max_eq_right hle] using hT)
        (by simpa [max_eq_right hle] using hScratch)
        (by simpa [max_eq_right hle] using hTEntry)
        (by simpa [max_eq_right hle] using hTQ)
        (by simpa [max_eq_right hle] using hTScratch)
        (by simpa [max_eq_right hle] using hScratchEntry)
        (by simpa [max_eq_right hle] using hScratchQ)
        hEntryRep hQRep with ⟨heap', hrun, hlayout, hrep⟩
    have heq : heap' = heap1 := Except.ok.inj (hrun.symm.trans (by
      simpa [horder] using hmul))
    subst heap'
    exact ⟨hlayout, by simpa [mul_comm, Nat.add_comm] using hrep⟩

/-- The same generated multiplication call leaves the old `i1` entry
byte-for-byte unchanged when that entry is outside both writable areas. -/
theorem matRowUpdate_mul_preserves_entry1 (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid) (i0 i1 : Fin 4)
    (Q T scratch : RawPtr UInt64) (lenQ : Nat) (heap heap1 : RawHeap)
    (entry1 : Polynomial (ZMod this._p.toNat))
    (hQpos : 0 < lenQ) (hEntryPos : 0 < hgcdMatLen M hM i0)
    (hT : heap.ValidU64Slice T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1))
    (hScratch : heap.ValidU64Slice scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hQValid : heap.ValidU64Slice Q lenQ)
    (hEntryValid : heap.ValidU64Slice (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0))
    (hScratchQ : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) Q lenQ)
    (hScratchEntry : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hTEntry1 : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hScratchEntry1 : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hlayout : RawHeap.SameLayout heap heap1)
    (hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this T
      (if lenQ ≥ hgcdMatLen M hM i0 then Q else hgcdMatPtr M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then lenQ else hgcdMatLen M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatPtr M hM i0 else Q)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatLen M hM i0 else lenQ)
      scratch heap = .ok heap1) :
    RawDensePolyRep this heap1 (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1 := by
  by_cases horder : lenQ ≥ hgcdMatLen M hM i0
  · have hsame := mul_preserves_prefix this T Q lenQ
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0) scratch
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1) heap heap1 hQpos
      hEntryPos horder
      (by simpa [max_eq_left horder] using hT) hQValid hEntryValid
      (by simpa [max_eq_left horder] using hScratch)
      (by simpa [max_eq_left horder] using hScratchEntry)
      (by simpa [max_eq_left horder] using hTEntry1)
      (by simpa [max_eq_left horder] using hScratchEntry1)
      (by simpa [horder] using hmul)
    exact rawDensePolyRep_of_same_prefix this heap heap1
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1) entry1 hlayout hsame
      hEntry1Rep
  · have hle : lenQ ≤ hgcdMatLen M hM i0 := by omega
    have hsame := mul_preserves_prefix this T (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) Q lenQ scratch (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) heap heap1 hEntryPos hQpos hle
      (by simpa [max_eq_right hle] using hT) hEntryValid hQValid
      (by simpa [max_eq_right hle] using hScratch)
      (by simpa [max_eq_right hle] using hScratchQ)
      (by simpa [max_eq_right hle] using hTEntry1)
      (by simpa [max_eq_right hle] using hScratchEntry1)
      (by simpa [horder] using hmul)
    exact rawDensePolyRep_of_same_prefix this heap heap1
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1) entry1 hlayout hsame
      hEntry1Rep

/-- The generated row addition writes only `t`; therefore the old `i0`
buffer installed as the new `i1` entry retains its normalized polynomial. -/
theorem matRowUpdate_add_preserves_entry0 (this : DenseUPolyZp)
    (t p1 T p0 : RawPtr UInt64) (l1 productLen l0 sumLen : Nat)
    (heap1 heap2 : RawHeap)
    (entry0 entry1 product : Polynomial (ZMod this._p.toNat))
    (hOutput : heap1.ValidU64Slice t (max l1 productLen))
    (hEntry1 : RawDensePolyRep this heap1 p1 l1 entry1)
    (hProduct : RawDensePolyRep this heap1 T productLen product)
    (hEntry0 : RawDensePolyRep this heap1 p0 l0 entry0)
    (htp0 : t.region ≠ p0.region)
    (hadd : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir this t
      p1 l1 T productLen heap1 = .ok (heap2, sumLen)) :
    RawDensePolyRep this heap2 p0 l0 entry0 := by
  rcases polyAdd_ok this t p1 l1 T productLen heap1 hOutput hEntry1.1
      hProduct.1 with ⟨heap', length, hrun, hlayout, _⟩
  have heq : (heap', length) = (heap2, sumLen) :=
    Except.ok.inj (hrun.symm.trans hadd)
  have hheap : heap' = heap2 := congrArg Prod.fst heq
  subst heap'
  have hsame := polyAdd_preserves_prefix_region_ne this t p1 l1 T p0
    productLen l0 heap1 heap2 sumLen hOutput hEntry1.1 hProduct.1 htp0
    hadd
  exact rawDensePolyRep_of_same_prefix this heap1 heap2 p0 l0 entry0
    hlayout hsame hEntry0

/-- Algebraic result of the actual add call exposed by the nonzero row-update
branch.  The product premise is supplied by strict `_mul`; this theorem then
binds the generated `_poly_add` result to the descriptor installed in M[i0]. -/
theorem matRowUpdate_nonzero_sum_rep (this : DenseUPolyZp)
    (M' : HgcdMat) (hM' : M'.Valid) (i0 : Fin 4)
    (t p1 T : RawPtr UInt64) (l1 productLen sumLen : Nat)
    (heap1 heap2 : RawHeap)
    (entry1 product : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hOutput : heap1.ValidU64Slice t (max l1 productLen))
    (hEntry1 : RawDensePolyRep this heap1 p1 l1 entry1)
    (hProduct : RawDensePolyRep this heap1 T productLen product)
    (hAliasEntry : ExactOrDisjoint t p1)
    (hAliasProduct : ExactOrDisjoint t T)
    (hadd : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir this t
      p1 l1 T productLen heap1 = .ok (heap2, sumLen))
    (hptr : hgcdMatPtr M' hM' i0 = t)
    (hlen : hgcdMatLen M' hM' i0 = sumLen) :
    RawDensePolyRep this heap2 (hgcdMatPtr M' hM' i0)
      (hgcdMatLen M' hM' i0) (entry1 + product) := by
  have hsum := polyAdd_refines this t p1 l1 T productLen heap1 heap2
    sumLen entry1 product hp hOutput hEntry1 hProduct hAliasEntry
    hAliasProduct hadd
  simpa [hptr, hlen] using hsum

/-- Complete semantic refinement of the successful nonzero C++ row update.
Both installed matrix entries are tied to the same generated `_mul` and
`_poly_add` executions exposed by `hrun`; no specification execution is
substituted for either call. -/
theorem matRowUpdate_nonzero_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1) (hQpos : 0 < lenQ)
    (hEntryPos : 0 < hgcdMatLen M hM i0)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hLenWord : max lenQ (hgcdMatLen M hM i0) < limbBase)
    (hT : heap.ValidU64Slice T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1))
    (hScratch : heap.ValidU64Slice scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hAddOutput : heap.ValidU64Slice t
      (max (hgcdMatLen M hM i1)
        (lenQ + hgcdMatLen M hM i0 - 1)))
    (hTQ : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) Q lenQ)
    (hTEntry0 : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hTEntry1 : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1)
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hTScratch : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) scratch
      (8 * max lenQ (hgcdMatLen M hM i0)))
    (hScratchQ : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) Q lenQ)
    (hScratchEntry0 : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0))
    (hScratchEntry1 : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0))
      (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1))
    (hAliasEntry1 : ExactOrDisjoint t (hgcdMatPtr M hM i1))
    (hAliasProduct : ExactOrDisjoint t T)
    (htEntry0 : t.region ≠ (hgcdMatPtr M hM i0).region)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      RawDensePolyRep this result.heap
        (hgcdMatPtr result.matrix hResult i0)
        (hgcdMatLen result.matrix hResult i0)
        (entry1 + quotient * entry0) ∧
      RawDensePolyRep this result.heap
        (hgcdMatPtr result.matrix hResult i1)
        (hgcdMatLen result.matrix hResult i1) entry0 := by
  rcases matRowUpdate_nonzero_success_shape this M i0 i1 Q lenQ T lenT t
      scratch heap result hM hne (by omega) (by omega) hrun with
    ⟨heap1, heap2, sumLen, hmul, hadd, hResultHeap, _, _, _,
      hResult, hptr0, hlen0, hptr1, hlen1⟩
  have hmulSemantic := matRowUpdate_mul_result this M hM i0 Q T scratch
    lenQ heap heap1 quotient entry0 hcfg hp hQpos hEntryPos hLenWord hT
    hScratch hTQ hTEntry0 hTScratch hScratchQ hScratchEntry0 hQRep
    hEntry0Rep hmul
  rcases hmulSemantic with ⟨hlayoutMul, hProduct1⟩
  have hEntry1_1 := matRowUpdate_mul_preserves_entry1 this M hM i0 i1 Q T
    scratch lenQ heap heap1 entry1 hQpos hEntryPos hT hScratch hQRep.1
    hEntry0Rep.1 hScratchQ hScratchEntry0 hTEntry1 hScratchEntry1
    hEntry1Rep hlayoutMul hmul
  have hEntry0_1 := matRowUpdate_mul_preserves_entry1 this M hM i0 i0 Q T
    scratch lenQ heap heap1 entry0 hQpos hEntryPos hT hScratch hQRep.1
    hEntry0Rep.1 hScratchQ hScratchEntry0 hTEntry0 hScratchEntry0
    hEntry0Rep hlayoutMul hmul
  have hAddOutput1 := (hlayoutMul t
    (max (hgcdMatLen M hM i1)
      (lenQ + hgcdMatLen M hM i0 - 1))).mp hAddOutput
  have hpWord : this._p ≠ 0 := by
    intro hzero
    have := congrArg UInt64.toNat hzero
    simp at this
    omega
  have hNew0 := matRowUpdate_nonzero_sum_rep this result.matrix hResult i0
    t (hgcdMatPtr M hM i1) T (hgcdMatLen M hM i1)
    (lenQ + hgcdMatLen M hM i0 - 1) sumLen heap1 heap2 entry1
    (quotient * entry0) hpWord hAddOutput1 hEntry1_1 hProduct1
    hAliasEntry1 hAliasProduct hadd hptr0 hlen0
  have hOld0_2 := matRowUpdate_add_preserves_entry0 this t
    (hgcdMatPtr M hM i1) T (hgcdMatPtr M hM i0)
    (hgcdMatLen M hM i1) (lenQ + hgcdMatLen M hM i0 - 1)
    (hgcdMatLen M hM i0) sumLen heap1 heap2 entry0 entry1
    (quotient * entry0) hAddOutput1 hEntry1_1 hProduct1 hEntry0_1
    htEntry0 hadd
  refine ⟨hResult, ?_, ?_⟩
  · simpa [hResultHeap] using hNew0
  · simpa [hResultHeap, hptr1, hlen1] using hOld0_2

end CLPoly.Impl.StrictHGCDRawRefinement
