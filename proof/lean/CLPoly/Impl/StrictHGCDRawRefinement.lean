import CLPoly.Generated.StrictHGCD
import CLPoly.Impl.StrictMulDispatchRefinement
import CLPoly.Impl.StrictPolyAddSubRefinement
import CLPoly.Impl.StrictHGCDRefinement

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

/-- Normalized raw representation of all four live HGCD matrix entries. -/
def HgcdMatRawDenseRep (this : DenseUPolyZp) (heap : RawHeap)
    (M : HgcdMat) (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid) : Prop :=
  ∀ i : Fin 4, RawDensePolyRep this heap (hgcdMatPtr M hM i)
    (hgcdMatLen M hM i) (entries i)

/-- Pointwise frame transport for the normalized raw matrix invariant. -/
theorem hgcdMatRawDenseRep_of_same_prefixes (this : DenseUPolyZp)
    (before after : RawHeap) (M : HgcdMat)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid) (hlayout : RawHeap.SameLayout before after)
    (hsame : ∀ i : Fin 4, SameU64Prefix before after
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hrep : HgcdMatRawDenseRep this before M entries hM) :
    HgcdMatRawDenseRep this after M entries hM := by
  intro i
  exact rawDensePolyRep_of_same_prefix this before after
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i) (entries i) hlayout
    (hsame i) (hrep i)

/-- A real recursive memcpy preserves the normalized raw matrix invariant
when its destination is disjoint from every live entry. -/
theorem copyU64_preserves_hgcdMatRawDenseRep (this : DenseUPolyZp)
    (heap heap' : RawHeap) (dst src : RawPtr UInt64) (count : Nat)
    (M : HgcdMat) (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid) (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hMatrix : ∀ i : Fin 4, U64SlicesDisjoint dst count
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i))
    (hcopy : heap.copyU64 dst src count = .ok heap')
    (hrep : HgcdMatRawDenseRep this heap M entries hM) :
    RawHeap.SameLayout heap heap' ∧
      HgcdMatRawDenseRep this heap' M entries hM := by
  rcases copyU64_ok heap dst src count hDst hSrc with
    ⟨copyHeap, hcopy', hlayout⟩
  have heq : copyHeap = heap' := Except.ok.inj (hcopy'.symm.trans hcopy)
  subst copyHeap
  refine ⟨hlayout, hgcdMatRawDenseRep_of_same_prefixes this heap heap' M
    entries hM hlayout ?_ hrep⟩
  intro i
  exact copyU64_preserves_prefix heap heap' dst src (hgcdMatPtr M hM i)
    count (hgcdMatLen M hM i) hDst hSrc (hMatrix i) hcopy

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

/-- A zero-length C++ coefficient slice is the normalized raw
representation of the zero polynomial. -/
theorem rawDensePolyRep_zero_length (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64)
    (hvalid : heap.ValidU64Slice ptr 0) :
    RawDensePolyRep this heap ptr 0 0 := by
  refine ⟨hvalid, ?_, slicePolyRep_zero_length_any heap ptr this._p.toNat, ?_⟩
  · intro k value hk _
    omega
  · simp [RawHeap.normaliseU64]

/-- A readable limb `1` is the normalized raw representation of the constant
one whenever the C++ modulus has at least two residues. -/
theorem rawDensePolyRep_one_of_read_one (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64)
    (hp : 1 < this._p.toNat)
    (hvalid : heap.ValidU64Slice ptr 1)
    (hread : heap.readU64 ptr 0 = .ok 1) :
    RawDensePolyRep this heap ptr 1 1 := by
  refine ⟨hvalid, ?_, slicePolyRep_one_of_read_one heap ptr
    this._p.toNat hvalid hread, ?_⟩
  · intro k value hk hvalue
    have hk0 : k = 0 := by omega
    subst k
    have hv : value = 1 := Except.ok.inj (hvalue.symm.trans hread)
    subst value
    simpa using hp
  · simp [RawHeap.normaliseU64, hread]

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
      heap'.readU64 (hgcdMatPtr M hM ⟨0, by omega⟩) 0 = .ok 1 ∧
      heap'.readU64 (hgcdMatPtr M hM ⟨3, by omega⟩) 0 = .ok 1 ∧
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
  refine ⟨heap2, M', hrun, fun ptr length =>
    (hlayout1 ptr length).trans (hlayout2 ptr length), rfl, rfl,
    (by simpa [p0]), (by simpa [p3]), hM', ?_⟩
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
    (hp : 1 < this._p.toNat)
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
        HgcdMatRawDenseRep this initial.heap initial.matrix
          (identityEntries this._p.toNat) hInitialM ∧
        RawDensePolyRep this initial.heap initial.A initial.lenA left ∧
        RawDensePolyRep this initial.heap initial.B initial.lenB right ∧
        CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right left right
          (identityEntries this._p.toNat ⟨0, by omega⟩)
          (identityEntries this._p.toNat ⟨1, by omega⟩)
          (identityEntries this._p.toNat ⟨2, by omega⟩)
          (identityEntries this._p.toNat ⟨3, by omega⟩) ∧
        CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet initial.sgn
          (identityEntries this._p.toNat ⟨0, by omega⟩)
          (identityEntries this._p.toNat ⟨1, by omega⟩)
          (identityEntries this._p.toNat ⟨2, by omega⟩)
          (identityEntries this._p.toNat ⟨3, by omega⟩) := by
  rcases matOne_refines M heap this._p.toNat hM h0 h3 h03 with
    ⟨heap1, M1, hone, hlayout1, hpoly1, hlen1, hread0, hread3,
      hM1, hMatrix1⟩
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
  have hRawMatrix1 : HgcdMatRawDenseRep this heap1 M1
      (identityEntries this._p.toNat) hM1 := by
    intro i
    fin_cases i
    · simp only [identityEntries]
      rw [hLen1 ⟨0, by omega⟩]
      apply rawDensePolyRep_one_of_read_one this heap1 _ hp
      · convert hValidMatrix1 ⟨0, by omega⟩ using 1
        simpa [identityEntryLen] using (hLen1 ⟨0, by omega⟩).symm
      · rw [hPtr1 ⟨0, by omega⟩]
        exact hread0
    · simp only [identityEntries]
      rw [hLen1 ⟨1, by omega⟩]
      apply rawDensePolyRep_zero_length this heap1
      convert hValidMatrix1 ⟨1, by omega⟩ using 1
      simpa [identityEntryLen] using (hLen1 ⟨1, by omega⟩).symm
    · simp only [identityEntries]
      rw [hLen1 ⟨2, by omega⟩]
      apply rawDensePolyRep_zero_length this heap1
      convert hValidMatrix1 ⟨2, by omega⟩ using 1
      simpa [identityEntryLen] using (hLen1 ⟨2, by omega⟩).symm
    · simp only [identityEntries]
      rw [hLen1 ⟨3, by omega⟩]
      apply rawDensePolyRep_one_of_read_one this heap1 _ hp
      · convert hValidMatrix1 ⟨3, by omega⟩ using 1
        simpa [identityEntryLen] using (hLen1 ⟨3, by omega⟩).symm
      · rw [hPtr1 ⟨3, by omega⟩]
        exact hread3
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
  rcases copyU64_preserves_hgcdMatRawDenseRep this heap1 heap2 A a lenA
      M1 (identityEntries this._p.toNat) hM1 hA1 hLeft1.1 hAMatrix1
      hcopyA hRawMatrix1 with ⟨_, hMatrix2⟩
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
  rcases copyU64_preserves_hgcdMatRawDenseRep this heap2 heap3 B b lenB
      M1 (identityEntries this._p.toNat) hM1 hB2 hRight2.1 hBMatrix1
      hcopyB hMatrix2 with
    ⟨_, hMatrix3⟩
  let initial : HgcdIterState := {
    heap := heap3, matrix := M1, A := A, lenA := lenA, B := B,
    lenB := lenB, T := T, lenT := lenT, t := t, sgn := 1 }
  refine ⟨initial, ?_, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl,
    hM1, ?_, ?_, ?_, ?_, ?_⟩
  · simp [hgcdIterInit, hone, hcopyA, hcopyB, initial]
  · simpa [initial] using hMatrix3
  · simpa [initial] using hA3
  · simpa [initial] using hB3
  · simp [CLPoly.Impl.StrictHGCDRefinement.HgcdTransform, identityEntries]
  · simp [initial, CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet,
      identityEntries]

/-- A stopped generated HGCD loop returns its state unchanged. -/
theorem hgcdIterLoop_stop (this : DenseUPolyZp) (m : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (state : HgcdIterState) (hstop : state.lenB < m + 1) :
    hgcdIterLoop this m Q W3 scratch state = .ok state := by
  rw [hgcdIterLoop]
  simp [show ¬state.lenB ≥ m + 1 by omega]

/-- Every successful nonterminal iteration is exactly one generated divrem,
the source pointer rotation, two generated row updates, and the recursive
tail call on the shorter remainder. -/
theorem hgcdIterLoop_step_shape (this : DenseUPolyZp) (m : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (state final : HgcdIterState) (hguard : state.lenB ≥ m + 1)
    (hrun : hgcdIterLoop this m Q W3 scratch state = .ok final) :
    ∃ heap1 lenQ lenR row23 row01,
      Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
        state.A state.lenA state.B state.lenB W3 state.heap =
          .ok (heap1, lenQ, lenR) ∧
      dense_upoly_zp__mat_row_update_ir this state.matrix
        ⟨2, by omega⟩ ⟨3, by omega⟩ Q lenQ state.A state.lenT state.t
        scratch heap1 = .ok row23 ∧
      dense_upoly_zp__mat_row_update_ir this row23.matrix
        ⟨0, by omega⟩ ⟨1, by omega⟩ Q lenQ row23.T row23.lenT row23.t
        scratch row23.heap = .ok row01 ∧
      hgcdIterLoop this m Q W3 scratch {
        heap := row01.heap
        matrix := row01.matrix
        A := state.B
        lenA := state.lenB
        B := state.T
        lenB := lenR
        T := row01.T
        lenT := row01.lenT
        t := row01.t
        sgn := -state.sgn
      } = .ok final ∧ lenR < state.lenB := by
  rw [hgcdIterLoop] at hrun
  simp only [if_pos hguard] at hrun
  split at hrun
  next fault hdiv => simp at hrun
  next =>
    rename_i heap1 lenQ lenR hdiv
    split at hrun
    next fault hrow23 => simp at hrun
    next row23 hrow23 =>
      split at hrun
      next fault hrow01 => simp at hrun
      next row01 hrow01 =>
        have hlt := Generated.StrictDivrem.polyDivrem_remainder_lt this Q
          state.T state.A
          state.lenA state.B state.lenB W3 state.heap heap1 lenQ lenR hdiv
        exact ⟨heap1, lenQ, lenR, row23, row01, hdiv, hrow23,
          hrow01, hrun, hlt⟩

/-- Semantic divrem closure for one successful nonterminal HGCD iteration.
All polynomial facts below concern the exact divrem result exposed by
`hgcdIterLoop_step_shape`; no second or specification-level execution is
substituted. -/
theorem hgcdIterLoop_step_divrem_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (m : Nat) (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch : RawPtr UInt64) (state final : HgcdIterState)
    (dividend divisor : Polynomial (ZMod this._p.toNat))
    (hguard : state.lenB ≥ m + 1)
    (hARep : RawDensePolyRep this state.heap state.A state.lenA dividend)
    (hBRep : RawDensePolyRep this state.heap state.B state.lenB divisor)
    (hQ : state.heap.ValidU64Slice Q
      (state.lenA - (state.lenB - 1)))
    (hR : state.heap.ValidU64Slice state.T
      (Nat.min state.lenA (state.lenB - 1)))
    (hW3 : state.heap.ValidWord3Slice W3 state.lenA)
    (hqCapacity : state.lenA - (state.lenB - 1) < limbBase)
    (hRA : state.T.region ≠ state.A.region)
    (hWA : W3.region ≠ state.A.region)
    (hWB : W3.region ≠ state.B.region)
    (hQB : Q.region ≠ state.B.region)
    (hQW : Q.region ≠ W3.region)
    (hRW : state.T.region ≠ W3.region)
    (hRQ : state.T.region ≠ Q.region)
    (hRB : state.T.region ≠ state.B.region)
    (hcfg : DensePreinvConfigured this)
    (hrun : hgcdIterLoop this m Q W3 scratch state = .ok final) :
    ∃ heap1 lenQ lenR quotient remainder row23 row01,
      Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
        state.A state.lenA state.B state.lenB W3 state.heap =
          .ok (heap1, lenQ, lenR) ∧
      RawDensePolyRep this heap1 Q lenQ quotient ∧
      RawDensePolyRep this heap1 state.B state.lenB divisor ∧
      RawDensePolyRep this heap1 state.T lenR remainder ∧
      normalize (EuclideanDomain.gcd dividend divisor) =
        normalize (EuclideanDomain.gcd divisor remainder) ∧
      dense_upoly_zp__mat_row_update_ir this state.matrix
        ⟨2, by omega⟩ ⟨3, by omega⟩ Q lenQ state.A state.lenT state.t
        scratch heap1 = .ok row23 ∧
      dense_upoly_zp__mat_row_update_ir this row23.matrix
        ⟨0, by omega⟩ ⟨1, by omega⟩ Q lenQ row23.T row23.lenT row23.t
        scratch row23.heap = .ok row01 ∧
      hgcdIterLoop this m Q W3 scratch {
        heap := row01.heap
        matrix := row01.matrix
        A := state.B
        lenA := state.lenB
        B := state.T
        lenB := lenR
        T := row01.T
        lenT := row01.lenT
        t := row01.t
        sgn := -state.sgn
      } = .ok final ∧
      lenR < state.lenB := by
  have hlenB : 0 < state.lenB := by omega
  rcases hgcdIterLoop_step_shape this m Q W3 scratch state final hguard hrun
      with ⟨heap1, lenQ, lenR, row23, row01, hdiv, hrow23, hrow01,
        htail, hlt⟩
  rcases polyDivrem_next_state this Q state.T state.A state.B state.lenA
      state.lenB W3 state.heap dividend divisor hlenB hARep hBRep hQ hR
      hW3 hqCapacity hRA hWA hWB hQB hQW hRW hRQ hRB hcfg with
    ⟨semanticHeap, semanticLenQ, semanticLenR, quotient, remainder,
      hsemantic, hQRep, hBRep1, hRRep, hgcd, _, _, _, _⟩
  have heq : (semanticHeap, semanticLenQ, semanticLenR) =
      (heap1, lenQ, lenR) := Except.ok.inj (hsemantic.symm.trans hdiv)
  cases heq
  exact ⟨heap1, lenQ, lenR, quotient, remainder, row23, row01, hdiv,
    hQRep, hBRep1, hRRep, hgcd, hrow23, hrow01, htail, hlt⟩

/-- The source's zero-quotient/zero-entry branch performs exactly the two
matrix-entry swaps and no heap access.  This exposes the real descriptor
state consumed by the next HGCD iteration. -/
theorem matRowUpdate_zero_exec (this : DenseUPolyZp) (M : HgcdMat)
    (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (hM : M.Valid)
    (hne : i0 ≠ i1)
    (hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0) :
    ∃ result,
      dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T lenT t
          scratch heap = .ok result ∧
      result.heap = heap ∧ result.T = T ∧ result.lenT = lenT ∧
      result.t = t ∧
      ∃ hResult : result.matrix.Valid,
        hgcdMatPtr result.matrix hResult i0 = hgcdMatPtr M hM i1 ∧
        hgcdMatLen result.matrix hResult i0 = hgcdMatLen M hM i1 ∧
        hgcdMatPtr result.matrix hResult i1 = hgcdMatPtr M hM i0 ∧
        hgcdMatLen result.matrix hResult i1 = hgcdMatLen M hM i0 ∧
        ∀ j : Fin 4, j ≠ i0 → j ≠ i1 →
          hgcdMatPtr result.matrix hResult j = hgcdMatPtr M hM j ∧
          hgcdMatLen result.matrix hResult j = hgcdMatLen M hM j := by
  have hvalid : M.poly.size = 4 ∧ M.len.size = 4 := by
    simpa [HgcdMat.Valid] using hM
  have hneVal : i0.val ≠ i1.val := by
    intro heq
    exact hne (Fin.ext heq)
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
  have hResult : result.matrix.Valid := by
    simp [result, HgcdMat.Valid, poly', len', hvalid]
  refine ⟨result, hrun, rfl, rfl, rfl, rfl, hResult, ?_, ?_, ?_, ?_, ?_⟩
  · simp [result, poly', p1, hgcdMatPtr]
  · simp [result, len', l1, hgcdMatLen]
  · simp only [result, poly', hgcdMatPtr]
    rw [Array.getElem_set_ne
      (by simpa using (show i0.val < M.poly.size by
        rw [hvalid.1]; exact i0.isLt))
      (by simpa using (show i1.val < M.poly.size by
        rw [hvalid.1]; exact i1.isLt)) hneVal,
      Array.getElem_set_self]
    simp [p0, hgcdMatPtr]
  · simp only [result, len', hgcdMatLen]
    rw [Array.getElem_set_ne
      (by simpa using (show i0.val < M.len.size by
        rw [hvalid.2]; exact i0.isLt))
      (by simpa using (show i1.val < M.len.size by
        rw [hvalid.2]; exact i1.isLt)) hneVal,
      Array.getElem_set_self]
    simp [l0, hgcdMatLen]
  · intro j hji0 hji1
    have hi0j : i0.val ≠ j.val := by
      intro heq
      exact hji0 (Fin.ext heq.symm)
    have hi1j : i1.val ≠ j.val := by
      intro heq
      exact hji1 (Fin.ext heq.symm)
    have hi0PolySet : i0.val < (M.poly.set i1.val p0 (by omega)).size := by
      simpa [hvalid.1] using i0.isLt
    have hjPolySet : j.val < (M.poly.set i1.val p0 (by omega)).size := by
      simpa [hvalid.1] using j.isLt
    have hi1Poly : i1.val < M.poly.size := by simpa [hvalid.1] using i1.isLt
    have hjPoly : j.val < M.poly.size := by simpa [hvalid.1] using j.isLt
    have hi0LenSet : i0.val < (M.len.set i1.val l0 (by omega)).size := by
      simpa [hvalid.2] using i0.isLt
    have hjLenSet : j.val < (M.len.set i1.val l0 (by omega)).size := by
      simpa [hvalid.2] using j.isLt
    have hi1Len : i1.val < M.len.size := by simpa [hvalid.2] using i1.isLt
    have hjLen : j.val < M.len.size := by simpa [hvalid.2] using j.isLt
    constructor
    · simp only [result, poly', hgcdMatPtr]
      rw [Array.getElem_set_ne hi0PolySet hjPolySet hi0j,
        Array.getElem_set_ne hi1Poly hjPoly hi1j]
    · simp only [result, len', hgcdMatLen]
      rw [Array.getElem_set_ne hi0LenSet hjLenSet hi0j,
        Array.getElem_set_ne hi1Len hjLen hi1j]

/-- Complete semantic refinement of the source's zero row-update branch.
The physical descriptor swap represents `[entry1 + quotient*entry0,
entry0]`; the vanishing factor is derived from the actual zero-length raw
representation, not assumed at L2. -/
theorem matRowUpdate_zero_refines (this : DenseUPolyZp)
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1)
    (hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0)
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
  rcases matRowUpdate_zero_exec this M i0 i1 Q lenQ T lenT t scratch
      heap hM hne hzero with
    ⟨actual, hactual, hheap, _, _, _, hActualM, hptr0, hlen0,
      hptr1, hlen1, _⟩
  have heq : actual = result := Except.ok.inj (hactual.symm.trans hrun)
  subst actual
  have hvanish : quotient = 0 ∨ entry0 = 0 := by
    rcases hzero with hq | he
    · left
      subst lenQ
      exact slicePolyRep_zero_length heap Q this._p.toNat quotient
        hQRep.2.2.1
    · right
      have hrep0 : SlicePolyRep heap (hgcdMatPtr M hM i0) 0
          this._p.toNat entry0 := by
        simpa [he] using hEntry0Rep.2.2.1
      exact slicePolyRep_zero_length heap (hgcdMatPtr M hM i0)
        this._p.toNat entry0 hrep0
  refine ⟨hActualM, ?_, ?_⟩
  · rcases hvanish with hq | he
    · simpa [hheap, hptr0, hlen0, hq] using hEntry1Rep
    · simpa [hheap, hptr0, hlen0, he] using hEntry1Rep
  · simpa [hheap, hptr1, hlen1] using hEntry0Rep

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
        hgcdMatLen result.matrix hResult i1 = hgcdMatLen M hM i0 ∧
        ∀ j : Fin 4, j ≠ i0 → j ≠ i1 →
          hgcdMatPtr result.matrix hResult j = hgcdMatPtr M hM j ∧
          hgcdMatLen result.matrix hResult j = hgcdMatLen M hM j := by
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
            refine ⟨hResult, ?_, ?_, ?_, ?_, ?_⟩
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
            · intro j hji0 hji1
              have hi0j : i0.val ≠ j.val := by
                intro heq
                exact hji0 (Fin.ext heq.symm)
              have hi1j : i1.val ≠ j.val := by
                intro heq
                exact hji1 (Fin.ext heq.symm)
              have hi0PolySet : i0.val <
                  (M.poly.set i1.val p0 (by omega)).size := by
                simpa [hvalid.1] using i0.isLt
              have hjPolySet : j.val <
                  (M.poly.set i1.val p0 (by omega)).size := by
                simpa [hvalid.1] using j.isLt
              have hi1Poly : i1.val < M.poly.size := by
                simpa [hvalid.1] using i1.isLt
              have hjPoly : j.val < M.poly.size := by
                simpa [hvalid.1] using j.isLt
              have hi0LenSet : i0.val <
                  (M.len.set i1.val l0 (by omega)).size := by
                simpa [hvalid.2] using i0.isLt
              have hjLenSet : j.val <
                  (M.len.set i1.val l0 (by omega)).size := by
                simpa [hvalid.2] using j.isLt
              have hi1Len : i1.val < M.len.size := by
                simpa [hvalid.2] using i1.isLt
              have hjLen : j.val < M.len.size := by
                simpa [hvalid.2] using j.isLt
              constructor
              · simp only [hgcdMatPtr]
                rw [Array.getElem_set_ne hi0PolySet hjPolySet hi0j,
                  Array.getElem_set_ne hi1Poly hjPoly hi1j]
              · simp only [hgcdMatLen]
                rw [Array.getElem_set_ne hi0LenSet hjLenSet hi0j,
                  Array.getElem_set_ne hi1Len hjLen hi1j]

/-- Branch-complete descriptor frame for `_mat_row_update`: indices other
than the selected pair retain exactly their source pointer and length. -/
theorem matRowUpdate_preserves_other_descriptor (this : DenseUPolyZp)
    (M : HgcdMat) (i0 i1 j : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (hne : i0 ≠ i1) (hji0 : j ≠ i0) (hji1 : j ≠ i1)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      hgcdMatPtr result.matrix hResult j = hgcdMatPtr M hM j ∧
      hgcdMatLen result.matrix hResult j = hgcdMatLen M hM j := by
  by_cases hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0
  · rcases matRowUpdate_zero_exec this M i0 i1 Q lenQ T lenT t scratch
      heap hM hne hzero with
    ⟨actual, hactual, _, _, _, _, hActualM, _, _, _, _, hother⟩
    have heq : actual = result := Except.ok.inj (hactual.symm.trans hrun)
    subst actual
    exact ⟨hActualM, hother j hji0 hji1⟩
  · rcases matRowUpdate_nonzero_success_shape this M i0 i1 Q lenQ T lenT
      t scratch heap result hM hne (by omega) (by omega) hrun with
    ⟨_, _, _, _, _, _, _, _, _, hResult, _, _, _, _, hother⟩
    exact ⟨hResult, hother j hji0 hji1⟩

/-- Descriptor frame specialized to the two source row updates in one HGCD
iteration: `(2,3)` leaves `(0,1)` unchanged, then `(0,1)` leaves `(2,3)`
unchanged. -/
theorem hgcdTwoRowUpdates_descriptor_frame (this : DenseUPolyZp)
    (M : HgcdMat) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (row23 row01 : MatRowUpdateResult) (hM : M.Valid)
    (hrow23 : dense_upoly_zp__mat_row_update_ir this M
      ⟨2, by omega⟩ ⟨3, by omega⟩ Q lenQ T lenT t scratch heap = .ok row23)
    (hrow01 : dense_upoly_zp__mat_row_update_ir this row23.matrix
      ⟨0, by omega⟩ ⟨1, by omega⟩ Q lenQ row23.T row23.lenT row23.t
      scratch row23.heap = .ok row01) :
    ∃ h23 : row23.matrix.Valid, ∃ h01 : row01.matrix.Valid,
      (∀ j : Fin 2,
        hgcdMatPtr row23.matrix h23 ⟨j.val, by omega⟩ =
          hgcdMatPtr M hM ⟨j.val, by omega⟩ ∧
        hgcdMatLen row23.matrix h23 ⟨j.val, by omega⟩ =
          hgcdMatLen M hM ⟨j.val, by omega⟩) ∧
      (∀ j : Fin 2,
        hgcdMatPtr row01.matrix h01 ⟨j.val + 2, by omega⟩ =
          hgcdMatPtr row23.matrix h23 ⟨j.val + 2, by omega⟩ ∧
        hgcdMatLen row01.matrix h01 ⟨j.val + 2, by omega⟩ =
          hgcdMatLen row23.matrix h23 ⟨j.val + 2, by omega⟩) := by
  rcases matRowUpdate_preserves_other_descriptor this M
      (2 : Fin 4) (3 : Fin 4) (0 : Fin 4) Q lenQ T lenT t
      scratch heap row23 hM (by decide) (by decide) (by decide)
      hrow23 with
    ⟨h23, hzero⟩
  rcases matRowUpdate_preserves_other_descriptor this row23.matrix
      (0 : Fin 4) (1 : Fin 4) (2 : Fin 4) Q lenQ row23.T
      row23.lenT row23.t scratch row23.heap row01 h23 (by decide)
      (by decide) (by decide) hrow01 with ⟨h01, htwo⟩
  refine ⟨h23, h01, ?_, ?_⟩
  · intro j
    fin_cases j
    · exact hzero
    · rcases matRowUpdate_preserves_other_descriptor this M
        (2 : Fin 4) (3 : Fin 4) (1 : Fin 4) Q lenQ T lenT t
        scratch heap row23 hM (by decide) (by decide) (by decide)
        hrow23 with
        ⟨h23', hone⟩
      simpa only using hone
  · intro j
    fin_cases j
    · exact htwo
    · rcases matRowUpdate_preserves_other_descriptor this row23.matrix
        (0 : Fin 4) (1 : Fin 4) (3 : Fin 4) Q lenQ row23.T
        row23.lenT row23.t scratch row23.heap row01 h23 (by decide)
        (by decide) (by decide) hrow01 with ⟨h01', hthree⟩
      simpa only using hthree

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

/-- Generic frame form of the generated multiplication inside a row update.
Any raw polynomial slice disjoint from both physical write areas (`T` and
`scratch`) is preserved. -/
theorem matRowUpdate_mul_preserves_guard (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid) (i0 : Fin 4)
    (Q T scratch guard : RawPtr UInt64) (lenQ guardLen : Nat)
    (heap heap1 : RawHeap) (guardPoly : Polynomial (ZMod this._p.toNat))
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
    (hTGuard : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) guard guardLen)
    (hScratchGuard : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) guard guardLen)
    (hGuardRep : RawDensePolyRep this heap guard guardLen guardPoly)
    (hlayout : RawHeap.SameLayout heap heap1)
    (hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this T
      (if lenQ ≥ hgcdMatLen M hM i0 then Q else hgcdMatPtr M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then lenQ else hgcdMatLen M hM i0)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatPtr M hM i0 else Q)
      (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatLen M hM i0 else lenQ)
      scratch heap = .ok heap1) :
    RawDensePolyRep this heap1 guard guardLen guardPoly := by
  by_cases horder : lenQ ≥ hgcdMatLen M hM i0
  · have hsame := mul_preserves_prefix this T Q lenQ
      (hgcdMatPtr M hM i0) (hgcdMatLen M hM i0) scratch guard guardLen
      heap heap1 hQpos hEntryPos horder
      (by simpa [max_eq_left horder] using hT) hQValid hEntryValid
      (by simpa [max_eq_left horder] using hScratch)
      (by simpa [max_eq_left horder] using hScratchEntry)
      (by simpa [max_eq_left horder] using hTGuard)
      (by simpa [max_eq_left horder] using hScratchGuard)
      (by simpa [horder] using hmul)
    exact rawDensePolyRep_of_same_prefix this heap heap1 guard guardLen
      guardPoly hlayout hsame hGuardRep
  · have hle : lenQ ≤ hgcdMatLen M hM i0 := by omega
    have hsame := mul_preserves_prefix this T (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) Q lenQ scratch guard guardLen heap heap1
      hEntryPos hQpos hle
      (by simpa [max_eq_right hle] using hT) hEntryValid hQValid
      (by simpa [max_eq_right hle] using hScratch)
      (by simpa [max_eq_right hle] using hScratchQ)
      (by simpa [max_eq_right hle] using hTGuard)
      (by simpa [max_eq_right hle] using hScratchGuard)
      (by simpa [horder] using hmul)
    exact rawDensePolyRep_of_same_prefix this heap heap1 guard guardLen
      guardPoly hlayout hsame hGuardRep

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
      hResult, hptr0, hlen0, hptr1, hlen1, _⟩
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

/-- Branch-complete semantic refinement of the generated C++ row update.
The same conclusion is obtained from either the physical swap-only branch or
the real `_mul`/`_poly_add` branch. -/
theorem matRowUpdate_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1)
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
  by_cases hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0
  · exact matRowUpdate_zero_refines this M i0 i1 Q lenQ T lenT t scratch
      heap result hM quotient entry0 entry1 hne hzero hQRep hEntry0Rep
      hEntry1Rep hrun
  · have hQpos : 0 < lenQ := Nat.pos_of_ne_zero (by
      intro heq
      exact hzero (Or.inl heq))
    have hEntryPos : 0 < hgcdMatLen M hM i0 := Nat.pos_of_ne_zero (by
      intro heq
      exact hzero (Or.inr heq))
    exact matRowUpdate_nonzero_refines this M i0 i1 Q lenQ T lenT t
      scratch heap result hM quotient entry0 entry1 hne hQpos hEntryPos
      hcfg hp hLenWord hT hScratch hAddOutput hTQ hTEntry0 hTEntry1
      hTScratch hScratchQ hScratchEntry0 hScratchEntry1 hAliasEntry1
      hAliasProduct htEntry0 hQRep hEntry0Rep hEntry1Rep hrun

/-- Branch-complete frame theorem for a generated row update.  A normalized
raw polynomial disjoint from `T`, `scratch`, and the add destination `t`
survives the exact `_mul`/`_poly_add` execution (or the swap-only branch). -/
theorem matRowUpdate_preserves_guard (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch guard : RawPtr UInt64)
    (guardLen : Nat) (heap : RawHeap) (result : MatRowUpdateResult)
    (hM : M.Valid)
    (quotient entry0 entry1 guardPoly : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1)
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
    (hTGuard : U64SlicesDisjoint T
      (2 * max lenQ (hgcdMatLen M hM i0) - 1) guard guardLen)
    (hScratchGuard : U64SlicesDisjoint scratch
      (8 * max lenQ (hgcdMatLen M hM i0)) guard guardLen)
    (htGuard : t.region ≠ guard.region)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hGuardRep : RawDensePolyRep this heap guard guardLen guardPoly)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    RawDensePolyRep this result.heap guard guardLen guardPoly := by
  by_cases hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0
  · rcases matRowUpdate_zero_exec this M i0 i1 Q lenQ T lenT t scratch
      heap hM hne hzero with ⟨actual, hactual, hheap, _⟩
    have heq : actual = result := Except.ok.inj (hactual.symm.trans hrun)
    subst actual
    simpa [hheap] using hGuardRep
  · have hQpos : 0 < lenQ := Nat.pos_of_ne_zero (by
      intro heq
      exact hzero (Or.inl heq))
    have hEntryPos : 0 < hgcdMatLen M hM i0 := Nat.pos_of_ne_zero (by
      intro heq
      exact hzero (Or.inr heq))
    rcases matRowUpdate_nonzero_success_shape this M i0 i1 Q lenQ T lenT
        t scratch heap result hM hne (by omega) (by omega) hrun with
      ⟨heap1, heap2, sumLen, hmul, hadd, hResultHeap, _, _, _, _, _, _, _,
        _, _⟩
    rcases matRowUpdate_mul_result this M hM i0 Q T scratch lenQ heap
        heap1 quotient entry0 hcfg hp hQpos hEntryPos hLenWord hT hScratch
        hTQ hTEntry0 hTScratch hScratchQ hScratchEntry0 hQRep hEntry0Rep
        hmul with ⟨hlayoutMul, hProduct1⟩
    have hGuard1 := matRowUpdate_mul_preserves_guard this M hM i0 Q T
      scratch guard lenQ guardLen heap heap1 guardPoly hQpos hEntryPos hT
      hScratch hQRep.1 hEntry0Rep.1 hScratchQ hScratchEntry0 hTGuard
      hScratchGuard hGuardRep hlayoutMul hmul
    have hEntry1_1 := matRowUpdate_mul_preserves_guard this M hM i0 Q T
      scratch (hgcdMatPtr M hM i1) lenQ (hgcdMatLen M hM i1) heap heap1
      entry1 hQpos hEntryPos hT hScratch hQRep.1 hEntry0Rep.1 hScratchQ
      hScratchEntry0 hTEntry1 hScratchEntry1
      hEntry1Rep hlayoutMul hmul
    have hAddOutput1 := (hlayoutMul t
      (max (hgcdMatLen M hM i1)
        (lenQ + hgcdMatLen M hM i0 - 1))).mp hAddOutput
    have hGuard2 := matRowUpdate_add_preserves_entry0 this t
      (hgcdMatPtr M hM i1) T guard (hgcdMatLen M hM i1)
      (lenQ + hgcdMatLen M hM i0 - 1) guardLen sumLen heap1 heap2
      guardPoly entry1 (quotient * entry0) hAddOutput1 hEntry1_1
      hProduct1 hGuard1 htGuard hadd
    simpa [hResultHeap] using hGuard2

/-- Physical L1 workspace invariant required by one generated
`_mat_row_update` call.  It records only allocation/aliasing facts about the
actual quotient, matrix entries, multiplication buffer, scratch buffer, and
addition destination; no L2 polynomial result is assumed. -/
structure MatRowUpdateWorkspace (M : HgcdMat) (i0 i1 : Fin 4)
    (Q : RawPtr UInt64) (lenQ : Nat) (T : RawPtr UInt64)
    (t scratch : RawPtr UInt64) (heap : RawHeap) (hM : M.Valid) : Prop where
  lenWord : max lenQ (hgcdMatLen M hM i0) < limbBase
  validT : heap.ValidU64Slice T
    (2 * max lenQ (hgcdMatLen M hM i0) - 1)
  validScratch : heap.ValidU64Slice scratch
    (8 * max lenQ (hgcdMatLen M hM i0))
  validAddOutput : heap.ValidU64Slice t
    (max (hgcdMatLen M hM i1)
      (lenQ + hgcdMatLen M hM i0 - 1))
  disjointTQ : U64SlicesDisjoint T
    (2 * max lenQ (hgcdMatLen M hM i0) - 1) Q lenQ
  disjointTMatrix : ∀ i : Fin 4, U64SlicesDisjoint T
    (2 * max lenQ (hgcdMatLen M hM i0) - 1)
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i)
  disjointTScratch : U64SlicesDisjoint T
    (2 * max lenQ (hgcdMatLen M hM i0) - 1) scratch
    (8 * max lenQ (hgcdMatLen M hM i0))
  disjointScratchQ : U64SlicesDisjoint scratch
    (8 * max lenQ (hgcdMatLen M hM i0)) Q lenQ
  disjointScratchMatrix : ∀ i : Fin 4, U64SlicesDisjoint scratch
    (8 * max lenQ (hgcdMatLen M hM i0))
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i)
  aliasEntry1 : ExactOrDisjoint t (hgcdMatPtr M hM i1)
  aliasProduct : ExactOrDisjoint t T
  tDisjointQ : t.region ≠ Q.region
  tDisjointOther : ∀ i : Fin 4, i ≠ i1 →
    t.region ≠ (hgcdMatPtr M hM i).region

/-- Apply the branch-complete semantic row-update theorem from the packaged
L1 workspace invariant. -/
theorem matRowUpdate_refines_of_workspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (workspace : MatRowUpdateWorkspace M i0 i1 Q lenQ T t scratch heap hM)
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
  exact matRowUpdate_refines this M i0 i1 Q lenQ T lenT t scratch heap
    result hM quotient entry0 entry1 hne hcfg hp workspace.lenWord
    workspace.validT workspace.validScratch workspace.validAddOutput
    workspace.disjointTQ (workspace.disjointTMatrix i0)
    (workspace.disjointTMatrix i1) workspace.disjointTScratch
    workspace.disjointScratchQ (workspace.disjointScratchMatrix i0)
    (workspace.disjointScratchMatrix i1) workspace.aliasEntry1
    workspace.aliasProduct (workspace.tDisjointOther i0 hne)
    hQRep hEntry0Rep hEntry1Rep hrun

/-- Any non-target raw polynomial covered by the same physical workspace
survives the exact generated row-update execution. -/
theorem matRowUpdate_preserves_matrix_entry_of_workspace
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 guardIndex : Fin 4)
    (Q : RawPtr UInt64) (lenQ : Nat) (T : RawPtr UInt64) (lenT : Nat)
    (t scratch : RawPtr UInt64) (heap : RawHeap)
    (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 guardPoly : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1) (hguard : guardIndex ≠ i1)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (workspace : MatRowUpdateWorkspace M i0 i1 Q lenQ T t scratch heap hM)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hGuardRep : RawDensePolyRep this heap (hgcdMatPtr M hM guardIndex)
      (hgcdMatLen M hM guardIndex) guardPoly)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    RawDensePolyRep this result.heap (hgcdMatPtr M hM guardIndex)
      (hgcdMatLen M hM guardIndex) guardPoly := by
  exact matRowUpdate_preserves_guard this M i0 i1 Q lenQ T lenT t scratch
    (hgcdMatPtr M hM guardIndex) (hgcdMatLen M hM guardIndex) heap result
    hM quotient entry0 entry1 guardPoly hne hcfg hp workspace.lenWord
    workspace.validT workspace.validScratch workspace.validAddOutput
    workspace.disjointTQ (workspace.disjointTMatrix i0)
    (workspace.disjointTMatrix i1) workspace.disjointTScratch
    workspace.disjointScratchQ (workspace.disjointScratchMatrix i0)
    (workspace.disjointScratchMatrix i1)
    (workspace.disjointTMatrix guardIndex)
    (workspace.disjointScratchMatrix guardIndex)
    (workspace.tDisjointOther guardIndex hguard) hQRep hEntry0Rep
    hEntry1Rep hGuardRep hrun

/-- The quotient itself survives one generated row update, which is needed
because the source reuses the same quotient for the second matrix row. -/
theorem matRowUpdate_preserves_quotient_of_workspace
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (workspace : MatRowUpdateWorkspace M i0 i1 Q lenQ T t scratch heap hM)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    RawDensePolyRep this result.heap Q lenQ quotient := by
  exact matRowUpdate_preserves_guard this M i0 i1 Q lenQ T lenT t scratch
    Q lenQ heap result hM quotient entry0 entry1 quotient hne hcfg hp
    workspace.lenWord workspace.validT workspace.validScratch
    workspace.validAddOutput workspace.disjointTQ
    (workspace.disjointTMatrix i0) (workspace.disjointTMatrix i1)
    workspace.disjointTScratch workspace.disjointScratchQ
    (workspace.disjointScratchMatrix i0)
    (workspace.disjointScratchMatrix i1) workspace.disjointTQ
    workspace.disjointScratchQ workspace.tDisjointQ hQRep hEntry0Rep
    hEntry1Rep hQRep hrun

/-- The two actual generated `_mat_row_update` calls refine the complete
four-entry mathematical HGCD matrix step.  In particular, this theorem
executes and frames both physical calls; it does not replace them with the
L2 matrix formula. -/
theorem hgcdTwoRowUpdates_refine_matrix (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (row23 row01 : MatRowUpdateResult) (hM : M.Valid)
    (h23 : row23.matrix.Valid)
    (quotient : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (workspace23 : MatRowUpdateWorkspace M (2 : Fin 4) (3 : Fin 4)
      Q lenQ T t scratch heap hM)
    (workspace01 : MatRowUpdateWorkspace row23.matrix (0 : Fin 4)
      (1 : Fin 4) Q lenQ row23.T row23.t scratch row23.heap h23)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hrow23 : dense_upoly_zp__mat_row_update_ir this M
      (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch heap = .ok row23)
    (hrow01 : dense_upoly_zp__mat_row_update_ir this row23.matrix
      (0 : Fin 4) (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t
      scratch row23.heap = .ok row01) :
    ∃ h01 : row01.matrix.Valid,
      HgcdMatRawDenseRep this row01.heap row01.matrix
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
        h01 := by
  have hE0 := hMatrix (0 : Fin 4)
  have hE1 := hMatrix (1 : Fin 4)
  have hE2 := hMatrix (2 : Fin 4)
  have hE3 := hMatrix (3 : Fin 4)
  rcases matRowUpdate_refines_of_workspace this M (2 : Fin 4)
      (3 : Fin 4) Q lenQ T lenT t scratch heap row23 hM quotient
      (entries 2) (entries 3) (by decide) hcfg hp workspace23 hQRep hE2
      hE3 hrow23 with ⟨h23', hNew2', hNew3'⟩
  have hh23 : h23' = h23 := Subsingleton.elim _ _
  subst h23'
  have hQ23 := matRowUpdate_preserves_quotient_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch heap row23 hM
    quotient (entries 2) (entries 3) (by decide) hcfg hp workspace23
    hQRep hE2 hE3 hrow23
  have hE0Old23 := matRowUpdate_preserves_matrix_entry_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) (0 : Fin 4) Q lenQ T lenT t scratch heap
    row23 hM quotient (entries 2) (entries 3) (entries 0) (by decide)
    (by decide) hcfg hp workspace23 hQRep hE2 hE3 hE0 hrow23
  have hE1Old23 := matRowUpdate_preserves_matrix_entry_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) (1 : Fin 4) Q lenQ T lenT t scratch heap
    row23 hM quotient (entries 2) (entries 3) (entries 1) (by decide)
    (by decide) hcfg hp workspace23 hQRep hE2 hE3 hE1 hrow23
  rcases hgcdTwoRowUpdates_descriptor_frame this M Q lenQ T lenT t scratch
      heap row23 row01 hM hrow23 hrow01 with
    ⟨h23Frame, h01Frame, hframe23, hframe01⟩
  have hh23Frame : h23Frame = h23 := Subsingleton.elim _ _
  subst h23Frame
  have hE0_23 : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (0 : Fin 4))
      (hgcdMatLen row23.matrix h23 (0 : Fin 4)) (entries 0) := by
    have hptr : hgcdMatPtr row23.matrix h23 (0 : Fin 4) =
        hgcdMatPtr M hM (0 : Fin 4) := by
      simpa only using (hframe23 (0 : Fin 2)).1
    have hlen : hgcdMatLen row23.matrix h23 (0 : Fin 4) =
        hgcdMatLen M hM (0 : Fin 4) := by
      simpa only using (hframe23 (0 : Fin 2)).2
    rw [hptr, hlen]
    exact hE0Old23
  have hE1_23 : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (1 : Fin 4))
      (hgcdMatLen row23.matrix h23 (1 : Fin 4)) (entries 1) := by
    have hptr : hgcdMatPtr row23.matrix h23 (1 : Fin 4) =
        hgcdMatPtr M hM (1 : Fin 4) := by
      simpa only using (hframe23 (1 : Fin 2)).1
    have hlen : hgcdMatLen row23.matrix h23 (1 : Fin 4) =
        hgcdMatLen M hM (1 : Fin 4) := by
      simpa only using (hframe23 (1 : Fin 2)).2
    rw [hptr, hlen]
    exact hE1Old23
  rcases matRowUpdate_refines_of_workspace this row23.matrix (0 : Fin 4)
      (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t scratch row23.heap
      row01 h23 quotient (entries 0) (entries 1) (by decide) hcfg hp
      workspace01 hQ23 hE0_23 hE1_23 hrow01 with
    ⟨h01, hNew0, hNew1⟩
  have hNew2 := matRowUpdate_preserves_matrix_entry_of_workspace this
    row23.matrix (0 : Fin 4) (1 : Fin 4) (2 : Fin 4) Q lenQ row23.T
    row23.lenT row23.t scratch row23.heap row01 h23 quotient
    (entries 0) (entries 1) (entries 3 + quotient * entries 2)
    (by decide) (by decide) hcfg hp workspace01 hQ23 hE0_23 hE1_23
    hNew2' hrow01
  have hNew3 := matRowUpdate_preserves_matrix_entry_of_workspace this
    row23.matrix (0 : Fin 4) (1 : Fin 4) (3 : Fin 4) Q lenQ row23.T
    row23.lenT row23.t scratch row23.heap row01 h23 quotient
    (entries 0) (entries 1) (entries 2) (by decide) (by decide) hcfg hp
    workspace01 hQ23 hE0_23 hE1_23 hNew3' hrow01
  have hh01Frame : h01Frame = h01 := Subsingleton.elim _ _
  subst h01Frame
  refine ⟨h01, ?_⟩
  intro i
  fin_cases i
  · simpa [CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries] using hNew0
  · simpa [CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries] using hNew1
  · have hptr : hgcdMatPtr row01.matrix h01 (2 : Fin 4) =
        hgcdMatPtr row23.matrix h23 (2 : Fin 4) := by
      simpa only using (hframe01 (0 : Fin 2)).1
    have hlen : hgcdMatLen row01.matrix h01 (2 : Fin 4) =
        hgcdMatLen row23.matrix h23 (2 : Fin 4) := by
      simpa only using (hframe01 (0 : Fin 2)).2
    change RawDensePolyRep this row01.heap
      (hgcdMatPtr row01.matrix h01 (2 : Fin 4))
      (hgcdMatLen row01.matrix h01 (2 : Fin 4))
      (entries 3 + quotient * entries 2)
    rw [hptr, hlen]
    exact hNew2
  · have hptr : hgcdMatPtr row01.matrix h01 (3 : Fin 4) =
        hgcdMatPtr row23.matrix h23 (3 : Fin 4) := by
      simpa only using (hframe01 (1 : Fin 2)).1
    have hlen : hgcdMatLen row01.matrix h01 (3 : Fin 4) =
        hgcdMatLen row23.matrix h23 (3 : Fin 4) := by
      simpa only using (hframe01 (1 : Fin 2)).2
    change RawDensePolyRep this row01.heap
      (hgcdMatPtr row01.matrix h01 (3 : Fin 4))
      (hgcdMatLen row01.matrix h01 (3 : Fin 4)) (entries 2)
    rw [hptr, hlen]
    exact hNew3

end CLPoly.Impl.StrictHGCDRawRefinement
