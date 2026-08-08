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

theorem array_getElem_eq_of_eq {α : Type} (left right : Array α)
    (h : left = right) (i : Nat) (hLeft : i < left.size)
    (hRight : i < right.size) :
    left[i]'hLeft = right[i]'hRight := by
  subst right
  rfl

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

/-- Purely physical obligations of the two C++ matrix-stabilization loops.
No L2 polynomial value or expected algorithm result is stored here. -/
structure HgcdMatStabilizeWorkspace (heap : RawHeap)
    (original current : HgcdMat) (hOriginal : original.Valid)
    (hCurrent : current.Valid) (stage : RawPtr UInt64) : Prop where
  stageValid : heap.ValidU64Slice stage
    (hgcdMatStageSize current hCurrent)
  stageCurrentDisjoint : ∀ i : Fin 4,
    U64SlicesDisjoint stage (hgcdMatStageSize current hCurrent)
      (hgcdMatPtr current hCurrent i) (hgcdMatLen current hCurrent i)
  originalValid : ∀ i : Fin 4, heap.ValidU64Slice
    (hgcdMatPtr original hOriginal i) (hgcdMatLen current hCurrent i)
  originalStageDisjoint : ∀ i : Fin 4,
    U64SlicesDisjoint (hgcdMatPtr original hOriginal i)
      (hgcdMatLen current hCurrent i) stage
      (hgcdMatStageSize current hCurrent)
  originalPairwiseDisjoint : ∀ i j : Fin 4, i ≠ j →
    U64SlicesDisjoint (hgcdMatPtr original hOriginal i)
      (hgcdMatLen current hCurrent i) (hgcdMatPtr original hOriginal j)
      (hgcdMatLen current hCurrent j)

theorem validU64Slice_stage_entry (heap : RawHeap) (stage : RawPtr UInt64)
    (M : HgcdMat) (hM : M.Valid) (i : Fin 4)
    (hvalid : heap.ValidU64Slice stage (hgcdMatStageSize M hM)) :
    heap.ValidU64Slice (stage.add (hgcdMatStageOffset M hM i.val))
      (hgcdMatLen M hM i) := by
  exact heap.validU64Slice_add stage (hgcdMatStageSize M hM)
    (hgcdMatStageOffset M hM i.val) (hgcdMatLen M hM i) hvalid
    (hgcdMatStageOffset_entry_le_size M hM i.val i.isLt)

theorem u64SlicesDisjoint_stage_entry_left (stage other : RawPtr UInt64)
    (total offset length otherLength : Nat)
    (hbound : offset + length ≤ total)
    (hdisjoint : U64SlicesDisjoint stage total other otherLength) :
    U64SlicesDisjoint (stage.add offset) length other otherLength := by
  intro i hi j hj
  rcases hdisjoint (offset + i) (by omega) j hj with hregion | hoff
  · exact Or.inl hregion
  · right
    have hwidth : RawLimbWidth.width UInt64 = 1 := rfl
    simpa [RawPtr.add, hwidth, Nat.add_assoc, Nat.add_left_comm, Nat.add_comm]
      using hoff

theorem u64SlicesDisjoint_stage_entry_right (other stage : RawPtr UInt64)
    (otherLength total offset length : Nat)
    (hbound : offset + length ≤ total)
    (hdisjoint : U64SlicesDisjoint other otherLength stage total) :
    U64SlicesDisjoint other otherLength (stage.add offset) length := by
  exact u64SlicesDisjoint_symm
    (u64SlicesDisjoint_stage_entry_left stage other total offset length
      otherLength hbound (u64SlicesDisjoint_symm hdisjoint))

theorem hgcdMatStabilizeWorkspace_of_sameLayout (before after : RawHeap)
    (original current : HgcdMat) (hOriginal : original.Valid)
    (hCurrent : current.Valid) (stage : RawPtr UInt64)
    (hlayout : RawHeap.SameLayout before after)
    (hws : HgcdMatStabilizeWorkspace before original current hOriginal
      hCurrent stage) :
    HgcdMatStabilizeWorkspace after original current hOriginal hCurrent stage := by
  exact {
    stageValid := (hlayout stage (hgcdMatStageSize current hCurrent)).mp
      hws.stageValid
    stageCurrentDisjoint := hws.stageCurrentDisjoint
    originalValid := fun i =>
      (hlayout (hgcdMatPtr original hOriginal i)
        (hgcdMatLen current hCurrent i)).mp (hws.originalValid i)
    originalStageDisjoint := hws.originalStageDisjoint
    originalPairwiseDisjoint := hws.originalPairwiseDisjoint }

/-- Raw representation of the four source-order entries after the first
stabilization loop has copied them into its contiguous stage buffer. -/
def HgcdMatStageRawDenseRep (this : DenseUPolyZp) (heap : RawHeap)
    (stage : RawPtr UInt64) (M : HgcdMat)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid) : Prop :=
  ∀ i : Fin 4, RawDensePolyRep this heap
    (stage.add (hgcdMatStageOffset M hM i.val))
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

/-- A represented raw polynomial disjoint from the memcpy destination is a
frame of that actual generated copy. -/
theorem copyU64_preserves_rawDenseRep (this : DenseUPolyZp)
    (heap heap' : RawHeap) (dst src : RawPtr UInt64) (count : Nat)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hdisjoint : U64SlicesDisjoint dst count ptr length)
    (hcopy : heap.copyU64 dst src count = .ok heap')
    (hrep : RawDensePolyRep this heap ptr length poly) :
    RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' ptr length poly := by
  rcases copyU64_ok heap dst src count hDst hSrc with
    ⟨copyHeap, hcopy', hlayout⟩
  have heq : copyHeap = heap' := Except.ok.inj (hcopy'.symm.trans hcopy)
  subst copyHeap
  exact ⟨hlayout, rawDensePolyRep_of_same_prefix this heap heap' ptr length
    poly hlayout (copyU64_preserves_prefix heap heap' dst src ptr count length
      hDst hSrc hdisjoint hcopy) hrep⟩

/-- Semantic refinement of the first actual stabilization loop.  The induction
follows its generated `i` recursion and accumulates the already copied stage
entries; the workspace contributes only validity and non-aliasing. -/
theorem hgcdMatStageLoop_refines (this : DenseUPolyZp)
    (M original : HgcdMat) (hM : M.Valid) (hOriginal : original.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (hi : i ≤ 4) (heap : RawHeap)
    (hws : HgcdMatStabilizeWorkspace heap original M hOriginal hM stage)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hPrior : ∀ j : Fin 4, j.val < i → RawDensePolyRep this heap
      (stage.add (hgcdMatStageOffset M hM j.val))
      (hgcdMatLen M hM j) (entries j)) :
    ∃ result, hgcdMatStageLoop M hM stage i
        (hgcdMatStageOffset M hM i) heap = .ok result ∧
      result.off = hgcdMatStageSize M hM ∧
      RawHeap.SameLayout heap result.heap ∧
      HgcdMatRawDenseRep this result.heap M entries hM ∧
      HgcdMatStageRawDenseRep this result.heap stage M entries hM := by
  rw [hgcdMatStageLoop]
  split
  next hlt =>
    dsimp only
    let index : Fin 4 := ⟨i, hlt⟩
    let dst := stage.add (hgcdMatStageOffset M hM i)
    let length := hgcdMatLen M hM index
    have hDst : heap.ValidU64Slice dst length :=
      validU64Slice_stage_entry heap stage M hM index hws.stageValid
    have hSrc := hMatrix index
    have hDstSrc : U64SlicesDisjoint dst length
        (hgcdMatPtr M hM index) length := by
      exact u64SlicesDisjoint_stage_entry_left stage (hgcdMatPtr M hM index)
        (hgcdMatStageSize M hM) (hgcdMatStageOffset M hM i) length length
        (hgcdMatStageOffset_entry_le_size M hM i hlt)
        (hws.stageCurrentDisjoint index)
    rcases copyU64_refines_rawDense this heap dst (hgcdMatPtr M hM index)
        length (entries index) hDst hDstSrc hSrc with
      ⟨heap1, hcopy, hlayout, hStageI⟩
    have hMatrix1 : HgcdMatRawDenseRep this heap1 M entries hM := by
      apply (copyU64_preserves_hgcdMatRawDenseRep this heap heap1 dst
        (hgcdMatPtr M hM index) length M entries hM hDst hSrc.1 ?_
        hcopy hMatrix).2
      intro k
      exact u64SlicesDisjoint_stage_entry_left stage (hgcdMatPtr M hM k)
        (hgcdMatStageSize M hM) (hgcdMatStageOffset M hM i) length
        (hgcdMatLen M hM k)
        (hgcdMatStageOffset_entry_le_size M hM i hlt)
        (hws.stageCurrentDisjoint k)
    have hPrior1 : ∀ j : Fin 4, j.val < i + 1 →
        RawDensePolyRep this heap1
          (stage.add (hgcdMatStageOffset M hM j.val))
          (hgcdMatLen M hM j) (entries j) := by
      intro j hj
      by_cases hji : j.val = i
      · have hjeq : j = index := Fin.ext hji
        subst j
        simpa [dst, length, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
          hgcdMatLenRaw] using hStageI
      · have hjlt : j.val < i := by omega
        have hOld := hPrior j hjlt
        have hBefore :
            hgcdMatStageOffset M hM j.val + hgcdMatLen M hM j ≤
              hgcdMatStageOffset M hM i := by
          exact hgcdMatStageOffset_entry_le_later M hM j.val i hjlt hlt j.isLt
        have hDisjoint : U64SlicesDisjoint dst length
            (stage.add (hgcdMatStageOffset M hM j.val))
            (hgcdMatLen M hM j) := by
          exact u64SlicesDisjoint_symm
            (u64SlicesDisjoint_add_of_le stage
              (hgcdMatStageOffset M hM j.val) (hgcdMatLen M hM j)
              (hgcdMatStageOffset M hM i) length hBefore)
        exact (copyU64_preserves_rawDenseRep this heap heap1 dst
          (hgcdMatPtr M hM index) length
          (stage.add (hgcdMatStageOffset M hM j.val))
          (hgcdMatLen M hM j) (entries j) hDst hSrc.1 hDisjoint hcopy hOld).2
    have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap heap1 original M
      hOriginal hM stage hlayout hws
    have hcopyLoop : heap.copyU64
        (stage.add (hgcdMatStageOffset M hM i))
        (hgcdMatPtrRaw M hM index) (hgcdMatLenRaw M hM index) = .ok heap1 := by
      simpa [dst, length, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
        hgcdMatLenRaw] using hcopy
    have hrec := hgcdMatStageLoop_refines this M original hM hOriginal stage
      entries (i + 1) (by omega) heap1 hws1 hMatrix1 hPrior1
    rw [hcopyLoop]
    rcases hrec with ⟨final, hrunFinal, hoffFinal, hlayoutFinal,
      hMatrixFinal, hStageFinal⟩
    refine ⟨final, ?_, hoffFinal, ?_, hMatrixFinal, hStageFinal⟩
    · simpa [index, dst, length, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
        hgcdMatLenRaw, hgcdMatStageOffset_step M hM i hlt] using hrunFinal
    · intro ptr count
      exact (hlayout ptr count).trans (hlayoutFinal ptr count)
  next hstop =>
    have hi4 : i = 4 := by omega
    subst i
    refine ⟨{ heap := heap, off := hgcdMatStageOffset M hM 4 }, rfl,
      rfl, (fun _ _ => Iff.rfl), hMatrix, ?_⟩
    intro j
    exact hPrior j (by omega)
termination_by 4 - i
decreasing_by omega

/-- Every concrete copy in the staging loop is a frame for a represented
live polynomial allocated outside the staging region.  This follows the
generated loop itself; no semantic result is supplied by the caller. -/
theorem hgcdMatStageLoop_preserves_rawDenseRep (this : DenseUPolyZp)
    (M original : HgcdMat) (hM : M.Valid) (hOriginal : original.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (hi : i ≤ 4) (heap : RawHeap)
    (hws : HgcdMatStabilizeWorkspace heap original M hOriginal hM stage)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hStagePtr : stage.region ≠ ptr.region)
    (hRep : RawDensePolyRep this heap ptr length poly)
    (result : HgcdMatStageResult)
    (hrun : hgcdMatStageLoop M hM stage i
      (hgcdMatStageOffset M hM i) heap = .ok result) :
    RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap ptr length poly := by
  rw [hgcdMatStageLoop] at hrun
  split at hrun
  next hlt =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hlt⟩
      let dst := stage.add (hgcdMatStageOffset M hM i)
      let count := hgcdMatLen M hM index
      have hDst : heap.ValidU64Slice dst count :=
        validU64Slice_stage_entry heap stage M hM index hws.stageValid
      have hSrc := hMatrix index
      have hLiveDisjoint : U64SlicesDisjoint dst count ptr length := by
        apply u64SlicesDisjoint_of_region_ne
        simpa [dst, RawPtr.add] using hStagePtr
      have hFrame := copyU64_preserves_rawDenseRep this heap heap1 dst
        (hgcdMatPtr M hM index) count ptr length poly hDst hSrc.1
        hLiveDisjoint (by
          simpa [dst, count, index, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
            hgcdMatLenRaw] using hcopy) hRep
      have hMatrix1 : HgcdMatRawDenseRep this heap1 M entries hM := by
        apply (copyU64_preserves_hgcdMatRawDenseRep this heap heap1 dst
          (hgcdMatPtr M hM index) count M entries hM hDst hSrc.1 ?_ (by
            simpa [dst, count, index, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
              hgcdMatLenRaw] using hcopy) hMatrix).2
        intro k
        exact u64SlicesDisjoint_stage_entry_left stage (hgcdMatPtr M hM k)
          (hgcdMatStageSize M hM) (hgcdMatStageOffset M hM i) count
          (hgcdMatLen M hM k)
          (hgcdMatStageOffset_entry_le_size M hM i hlt)
          (hws.stageCurrentDisjoint k)
      have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap heap1
        original M hOriginal hM stage hFrame.1 hws
      have hrec := hgcdMatStageLoop_preserves_rawDenseRep this M original hM
        hOriginal stage entries (i + 1) (by omega) heap1 hws1 hMatrix1 ptr
        length poly hStagePtr hFrame.2 result
        (by simpa [index, dst, count, hgcdMatPtr, hgcdMatPtrRaw, hgcdMatLen,
          hgcdMatLenRaw, hgcdMatStageOffset_step M hM i hlt] using hrun)
      exact ⟨fun p n => (hFrame.1 p n).trans (hrec.1 p n), hrec.2⟩
  next hstop =>
    have heq := Except.ok.inj hrun
    subst result
    exact ⟨fun _ _ => Iff.rfl, hRep⟩
termination_by 4 - i
decreasing_by omega

/-- Entry-point refinement of the first generated stabilization loop. -/
theorem hgcdMatStageLoop_zero_refines (this : DenseUPolyZp)
    (M original : HgcdMat) (hM : M.Valid) (hOriginal : original.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap)
    (hws : HgcdMatStabilizeWorkspace heap original M hOriginal hM stage)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM) :
    ∃ result, hgcdMatStageLoop M hM stage 0 0 heap = .ok result ∧
      result.off = hgcdMatStageSize M hM ∧
      RawHeap.SameLayout heap result.heap ∧
      HgcdMatRawDenseRep this result.heap M entries hM ∧
      HgcdMatStageRawDenseRep this result.heap stage M entries hM := by
  simpa [hgcdMatStageOffset] using
    (hgcdMatStageLoop_refines this M original hM hOriginal stage entries 0
      (by omega) heap hws hMatrix (by intro j hj; omega))

/-- Semantic refinement of every suffix of the second generated stabilization
loop.  `base` fixes the iterator-produced lengths while `current` follows the
real descriptor pointer updates performed by the loop. -/
theorem hgcdMatRestoreLoop_refines (this : DenseUPolyZp)
    (original base current : HgcdMat)
    (hOriginal : original.Valid) (hBase : base.Valid)
    (hCurrent : current.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (hi : i ≤ 4) (heap : RawHeap)
    (result : HgcdMatRestoreResult)
    (hLen : current.len = base.len)
    (hws : HgcdMatStabilizeWorkspace heap original base hOriginal hBase stage)
    (hStage : HgcdMatStageRawDenseRep this heap stage base entries hBase)
    (hPrior : ∀ j : Fin 4, j.val < i → RawDensePolyRep this heap
      (hgcdMatPtr original hOriginal j) (hgcdMatLen base hBase j)
      (entries j))
    (hrun : hgcdMatRestoreLoop original current hOriginal hCurrent stage i
      (hgcdMatStageOffset base hBase i) heap = .ok result) :
    HgcdMatStageRawDenseRep this result.heap stage base entries hBase ∧
      ∀ j : Fin 4, RawDensePolyRep this result.heap
        (hgcdMatPtr original hOriginal j) (hgcdMatLen base hBase j)
        (entries j) := by
  rw [hgcdMatRestoreLoop] at hrun
  split at hrun
  next hlt =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hlt⟩
      let dst := hgcdMatPtr original hOriginal index
      let src := stage.add (hgcdMatStageOffset base hBase i)
      let length := hgcdMatLen base hBase index
      have hLength : hgcdMatLen current hCurrent index = length := by
        exact array_getElem_eq_of_eq current.len base.len hLen index.val
          (by rw [hCurrent.2]; exact index.isLt)
          (by rw [hBase.2]; exact index.isLt)
      have hLengthRaw : hgcdMatLenRaw current hCurrent index = length := by
        simpa [hgcdMatLenRaw, hgcdMatLen] using hLength
      have hDst : heap.ValidU64Slice dst length := hws.originalValid index
      have hSrc := hStage index
      have hDstSrc : U64SlicesDisjoint dst length src length := by
        exact u64SlicesDisjoint_stage_entry_right dst stage length
          (hgcdMatStageSize base hBase) (hgcdMatStageOffset base hBase i)
          length (hgcdMatStageOffset_entry_le_size base hBase i hlt)
          (hws.originalStageDisjoint index)
      rcases copyU64_refines_rawDense this heap dst src length (entries index)
          hDst hDstSrc hSrc with
        ⟨semanticHeap, hsemantic, hlayout, hRestoredI⟩
      have hcopyNorm : heap.copyU64 dst src length = .ok heap1 := by
        simpa [dst, src, length, index, hgcdMatPtr, hgcdMatPtrRaw,
          hLengthRaw] using hcopy
      have heq : semanticHeap = heap1 :=
        Except.ok.inj (hsemantic.symm.trans hcopyNorm)
      subst semanticHeap
      have hStage1 : HgcdMatStageRawDenseRep this heap1 stage base entries
          hBase := by
        intro k
        have hDisjoint : U64SlicesDisjoint dst length
            (stage.add (hgcdMatStageOffset base hBase k.val))
            (hgcdMatLen base hBase k) := by
          exact u64SlicesDisjoint_stage_entry_right dst stage length
            (hgcdMatStageSize base hBase)
            (hgcdMatStageOffset base hBase k.val) (hgcdMatLen base hBase k)
            (hgcdMatStageOffset_entry_le_size base hBase k.val k.isLt)
            (hws.originalStageDisjoint index)
        exact (copyU64_preserves_rawDenseRep this heap heap1 dst src length
          (stage.add (hgcdMatStageOffset base hBase k.val))
          (hgcdMatLen base hBase k) (entries k) hDst hSrc.1 hDisjoint
          hcopyNorm (hStage k)).2
      have hPrior1 : ∀ j : Fin 4, j.val < i + 1 →
          RawDensePolyRep this heap1 (hgcdMatPtr original hOriginal j)
            (hgcdMatLen base hBase j) (entries j) := by
        intro j hj
        by_cases hji : j.val = i
        · have hjeq : j = index := Fin.ext hji
          subst j
          exact hRestoredI
        · have hjlt : j.val < i := by omega
          have hDisjoint : U64SlicesDisjoint dst length
              (hgcdMatPtr original hOriginal j)
              (hgcdMatLen base hBase j) := by
            exact hws.originalPairwiseDisjoint index j (by
              intro heqIndex
              exact hji (congrArg Fin.val heqIndex).symm)
          exact (copyU64_preserves_rawDenseRep this heap heap1 dst src length
            (hgcdMatPtr original hOriginal j) (hgcdMatLen base hBase j)
            (entries j) hDst hSrc.1 hDisjoint hcopyNorm
            (hPrior j hjlt)).2
      let poly' := current.poly.set i (hgcdMatPtrRaw original hOriginal index)
        (by rw [hCurrent.1]; omega)
      let next : HgcdMat := { current with poly := poly' }
      have hNext : next.Valid := by
        exact ⟨by simp [next, poly', hCurrent.1], hCurrent.2⟩
      have hLenNext : next.len = base.len := by
        exact hLen
      have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap heap1 original
        base hOriginal hBase stage hlayout hws
      have htail : hgcdMatRestoreLoop original next hOriginal hNext stage
          (i + 1) (hgcdMatStageOffset base hBase (i + 1)) heap1 =
            .ok result := by
        rw [hgcdMatStageOffset_step base hBase i hlt]
        simpa [next, poly', index, hLengthRaw] using hrun
      exact hgcdMatRestoreLoop_refines this original base next hOriginal hBase
        hNext stage entries (i + 1) (by omega) heap1 result hLenNext hws1
        hStage1 hPrior1 htail
  next hstop =>
    have hi4 : i = 4 := by omega
    have heq := Except.ok.inj hrun
    subst result
    constructor
    · exact hStage
    · intro j
      exact hPrior j (by omega)
termination_by 4 - i
decreasing_by omega

/-- Every concrete copy in the restore loop is a frame for a represented
live polynomial outside all four saved matrix-entry regions. -/
theorem hgcdMatRestoreLoop_preserves_rawDenseRep (this : DenseUPolyZp)
    (original base current : HgcdMat)
    (hOriginal : original.Valid) (hBase : base.Valid)
    (hCurrent : current.Valid) (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (hi : i ≤ 4) (heap : RawHeap)
    (hLen : current.len = base.len)
    (hws : HgcdMatStabilizeWorkspace heap original base hOriginal hBase stage)
    (hStage : HgcdMatStageRawDenseRep this heap stage base entries hBase)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hOriginalPtr : ∀ j : Fin 4,
      (hgcdMatPtr original hOriginal j).region ≠ ptr.region)
    (hRep : RawDensePolyRep this heap ptr length poly)
    (result : HgcdMatRestoreResult)
    (hrun : hgcdMatRestoreLoop original current hOriginal hCurrent stage i
      (hgcdMatStageOffset base hBase i) heap = .ok result) :
    RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap ptr length poly := by
  rw [hgcdMatRestoreLoop] at hrun
  split at hrun
  next hlt =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hlt⟩
      let dst := hgcdMatPtr original hOriginal index
      let src := stage.add (hgcdMatStageOffset base hBase i)
      let count := hgcdMatLen base hBase index
      have hLength : hgcdMatLen current hCurrent index = count := by
        exact array_getElem_eq_of_eq current.len base.len hLen index.val
          (by rw [hCurrent.2]; exact index.isLt)
          (by rw [hBase.2]; exact index.isLt)
      have hLengthRaw : hgcdMatLenRaw current hCurrent index = count := by
        simpa [hgcdMatLenRaw, hgcdMatLen] using hLength
      have hDst : heap.ValidU64Slice dst count := hws.originalValid index
      have hSrc := hStage index
      have hLiveDisjoint : U64SlicesDisjoint dst count ptr length := by
        exact u64SlicesDisjoint_of_region_ne (hOriginalPtr index)
      have hcopyNorm : heap.copyU64 dst src count = .ok heap1 := by
        simpa [dst, src, count, index, hgcdMatPtr, hgcdMatPtrRaw,
          hLengthRaw] using hcopy
      have hFrame := copyU64_preserves_rawDenseRep this heap heap1 dst src
        count ptr length poly hDst hSrc.1 hLiveDisjoint hcopyNorm hRep
      have hStage1 : HgcdMatStageRawDenseRep this heap1 stage base entries
          hBase := by
        intro k
        have hDisjoint : U64SlicesDisjoint dst count
            (stage.add (hgcdMatStageOffset base hBase k.val))
            (hgcdMatLen base hBase k) := by
          exact u64SlicesDisjoint_stage_entry_right dst stage count
            (hgcdMatStageSize base hBase)
            (hgcdMatStageOffset base hBase k.val) (hgcdMatLen base hBase k)
            (hgcdMatStageOffset_entry_le_size base hBase k.val k.isLt)
            (hws.originalStageDisjoint index)
        exact (copyU64_preserves_rawDenseRep this heap heap1 dst src count
          (stage.add (hgcdMatStageOffset base hBase k.val))
          (hgcdMatLen base hBase k) (entries k) hDst hSrc.1 hDisjoint
          hcopyNorm (hStage k)).2
      let poly' := current.poly.set i
        (hgcdMatPtrRaw original hOriginal index)
        (by rw [hCurrent.1]; omega)
      let next : HgcdMat := { current with poly := poly' }
      have hNext : next.Valid := by
        exact ⟨by simp [next, poly', hCurrent.1], hCurrent.2⟩
      have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap heap1
        original base hOriginal hBase stage hFrame.1 hws
      have hrec := hgcdMatRestoreLoop_preserves_rawDenseRep this original base next
        hOriginal hBase hNext stage entries (i + 1) (by omega) heap1
        (by exact hLen) hws1 hStage1 ptr length poly hOriginalPtr hFrame.2
        result (by
          rw [hgcdMatStageOffset_step base hBase i hlt]
          simpa [next, poly', index, hLengthRaw] using hrun)
      exact ⟨fun p n => (hFrame.1 p n).trans (hrec.1 p n), hrec.2⟩
  next hstop =>
    have heq := Except.ok.inj hrun
    subst result
    exact ⟨fun _ _ => Iff.rfl, hRep⟩
termination_by 4 - i
decreasing_by omega

/-- Entry-point refinement of the second generated loop: staged polynomial
contents are copied back to all four saved pointers, and the returned matrix
descriptor represents those same entries at the iterator-produced lengths. -/
theorem hgcdMatRestoreLoop_zero_refines (this : DenseUPolyZp)
    (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap) (result : HgcdMatRestoreResult)
    (hws : HgcdMatStabilizeWorkspace heap original current hOriginal
      hCurrent stage)
    (hStage : HgcdMatStageRawDenseRep this heap stage current entries hCurrent)
    (hrun : hgcdMatRestoreLoop original current hOriginal hCurrent stage
      0 0 heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      HgcdMatStageRawDenseRep this result.heap stage current entries hCurrent ∧
      HgcdMatRawDenseRep this result.heap result.matrix entries hResult := by
  have hrefine := hgcdMatRestoreLoop_refines this original current current
    hOriginal hCurrent hCurrent stage entries 0 (by omega) heap result rfl hws
    hStage (by intro j hj; omega) (by
      simpa [hgcdMatStageOffset] using hrun)
  have hdesc := hgcdMatRestoreLoop_zero_descriptors original current hOriginal
    hCurrent stage heap result hrun
  refine ⟨hdesc.1, hrefine.1, ?_⟩
  intro j
  have hPtr : hgcdMatPtr result.matrix hdesc.1 j =
      hgcdMatPtr original hOriginal j := by
    exact array_getElem_eq_of_eq result.matrix.poly original.poly hdesc.2.1
      j.val (by rw [hdesc.1.1]; exact j.isLt)
      (by rw [hOriginal.1]; exact j.isLt)
  have hLength : hgcdMatLen result.matrix hdesc.1 j =
      hgcdMatLen current hCurrent j := by
    exact array_getElem_eq_of_eq result.matrix.len current.len hdesc.2.2
      j.val (by rw [hdesc.1.2]; exact j.isLt)
      (by rw [hCurrent.2]; exact j.isLt)
  rw [hPtr, hLength]
  exact hrefine.2 j

/-- The physical stabilization workspace is sufficient to make every suffix
of the generated restore loop terminate successfully. -/
theorem hgcdMatRestoreLoop_terminates (this : DenseUPolyZp)
    (original base current : HgcdMat)
    (hOriginal : original.Valid) (hBase : base.Valid)
    (hCurrent : current.Valid) (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (hi : i ≤ 4) (heap : RawHeap)
    (hLen : current.len = base.len)
    (hws : HgcdMatStabilizeWorkspace heap original base hOriginal hBase stage)
    (hStage : HgcdMatStageRawDenseRep this heap stage base entries hBase) :
    ∃ result, hgcdMatRestoreLoop original current hOriginal hCurrent stage i
      (hgcdMatStageOffset base hBase i) heap = .ok result := by
  rw [hgcdMatRestoreLoop]
  split
  next hlt =>
    let index : Fin 4 := ⟨i, hlt⟩
    let dst := hgcdMatPtr original hOriginal index
    let src := stage.add (hgcdMatStageOffset base hBase i)
    let length := hgcdMatLen base hBase index
    have hLength : hgcdMatLen current hCurrent index = length := by
      exact array_getElem_eq_of_eq current.len base.len hLen index.val
        (by rw [hCurrent.2]; exact index.isLt)
        (by rw [hBase.2]; exact index.isLt)
    have hLengthRaw : hgcdMatLenRaw current hCurrent index = length := by
      simpa [hgcdMatLenRaw, hgcdMatLen] using hLength
    have hBaseLengthRaw : hgcdMatLenRaw base hBase index = length := by
      rfl
    have hDst : heap.ValidU64Slice dst length := hws.originalValid index
    have hSrc := hStage index
    have hDstSrc : U64SlicesDisjoint dst length src length := by
      exact u64SlicesDisjoint_stage_entry_right dst stage length
        (hgcdMatStageSize base hBase) (hgcdMatStageOffset base hBase i) length
        (hgcdMatStageOffset_entry_le_size base hBase i hlt)
        (hws.originalStageDisjoint index)
    rcases copyU64_refines_rawDense this heap dst src length (entries index)
        hDst hDstSrc hSrc with ⟨heap1, hcopy, hlayout, _⟩
    have hStage1 : HgcdMatStageRawDenseRep this heap1 stage base entries
        hBase := by
      intro k
      have hDisjoint : U64SlicesDisjoint dst length
          (stage.add (hgcdMatStageOffset base hBase k.val))
          (hgcdMatLen base hBase k) := by
        exact u64SlicesDisjoint_stage_entry_right dst stage length
          (hgcdMatStageSize base hBase)
          (hgcdMatStageOffset base hBase k.val) (hgcdMatLen base hBase k)
          (hgcdMatStageOffset_entry_le_size base hBase k.val k.isLt)
          (hws.originalStageDisjoint index)
      exact (copyU64_preserves_rawDenseRep this heap heap1 dst src length
        (stage.add (hgcdMatStageOffset base hBase k.val))
        (hgcdMatLen base hBase k) (entries k) hDst hSrc.1 hDisjoint hcopy
        (hStage k)).2
    let poly' := current.poly.set i (hgcdMatPtrRaw original hOriginal index)
      (by rw [hCurrent.1]; omega)
    let next : HgcdMat := { current with poly := poly' }
    have hNext : next.Valid := by
      exact ⟨by simp [next, poly', hCurrent.1], hCurrent.2⟩
    have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap heap1 original
      base hOriginal hBase stage hlayout hws
    rcases hgcdMatRestoreLoop_terminates this original base next hOriginal hBase
        hNext stage entries (i + 1) (by omega) heap1 hLen hws1 hStage1 with
      ⟨result, htail⟩
    refine ⟨result, ?_⟩
    dsimp only
    rw [show heap.copyU64 (hgcdMatPtrRaw original hOriginal index) src
        (hgcdMatLenRaw current hCurrent index) = .ok heap1 by
      simpa [dst, length, hLengthRaw] using hcopy]
    simp only
    rw [hLengthRaw, ← hBaseLengthRaw,
      ← hgcdMatStageOffset_step base hBase i hlt]
    simpa [next, poly', index, hLengthRaw] using htail
  next hstop =>
    exact ⟨HgcdMatRestoreResult.mk heap current
      (hgcdMatStageOffset base hBase i), rfl⟩
termination_by 4 - i
decreasing_by omega

/-- Complete raw semantic refinement of the two-loop source stabilization
block, including successful termination under its physical workspace. -/
theorem hgcdMatStabilize_refines (this : DenseUPolyZp)
    (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap)
    (hws : HgcdMatStabilizeWorkspace heap original current hOriginal
      hCurrent stage)
    (hMatrix : HgcdMatRawDenseRep this heap current entries hCurrent) :
    ∃ result, hgcdMatStabilize original current hOriginal hCurrent stage heap =
        .ok result ∧
      ∃ hResult : result.matrix.Valid,
        HgcdMatRawDenseRep this result.heap result.matrix entries hResult := by
  rcases hgcdMatStageLoop_zero_refines this current original hCurrent hOriginal
      stage entries heap hws hMatrix with
    ⟨staged, hstage, hoff, hlayout, _, hStage⟩
  have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap staged.heap original
    current hOriginal hCurrent stage hlayout hws
  rcases hgcdMatRestoreLoop_terminates this original current current hOriginal
      hCurrent hCurrent stage entries 0 (by omega) staged.heap rfl hws1 hStage with
    ⟨result, hrestore⟩
  have hrestore0 : hgcdMatRestoreLoop original current hOriginal hCurrent stage
      0 0 staged.heap = .ok result := by
    simpa [hgcdMatStageOffset] using hrestore
  rcases hgcdMatRestoreLoop_zero_refines this original current hOriginal hCurrent
      stage entries staged.heap result hws1 hStage hrestore0 with
    ⟨hResult, _, hResultRep⟩
  refine ⟨result, ?_, hResult, hResultRep⟩
  simp [hgcdMatStabilize, hstage, hrestore0]

/-- The complete generated matrix stabilization is a frame for a live raw
polynomial outside both the staging allocation and the four saved matrix
entry allocations.  The result also has the same allocation layout. -/
theorem hgcdMatStabilize_preserves_rawDenseRep (this : DenseUPolyZp)
    (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap)
    (hws : HgcdMatStabilizeWorkspace heap original current hOriginal
      hCurrent stage)
    (hMatrix : HgcdMatRawDenseRep this heap current entries hCurrent)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hStagePtr : stage.region ≠ ptr.region)
    (hOriginalPtr : ∀ j : Fin 4,
      (hgcdMatPtr original hOriginal j).region ≠ ptr.region)
    (hRep : RawDensePolyRep this heap ptr length poly)
    (result : HgcdMatRestoreResult)
    (hrun : hgcdMatStabilize original current hOriginal hCurrent stage heap =
      .ok result) :
    RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap ptr length poly := by
  simp only [hgcdMatStabilize] at hrun
  generalize hstage : hgcdMatStageLoop current hCurrent stage 0 0 heap = staged
    at hrun
  cases staged with
  | error fault => simp at hrun
  | ok staged =>
    rcases hgcdMatStageLoop_zero_refines this current original hCurrent
        hOriginal stage entries heap hws hMatrix with
      ⟨semantic, hsemantic, _, hlayoutStage, _, hStage⟩
    have hStagedEq : semantic = staged :=
      Except.ok.inj (hsemantic.symm.trans hstage)
    subst semantic
    have hLiveStage := hgcdMatStageLoop_preserves_rawDenseRep this current
      original hCurrent hOriginal stage entries 0 (by omega) heap hws
      hMatrix ptr length poly hStagePtr hRep staged (by
        simpa [hgcdMatStageOffset] using hstage)
    have hws1 := hgcdMatStabilizeWorkspace_of_sameLayout heap staged.heap
      original current hOriginal hCurrent stage hlayoutStage hws
    have hrestore : hgcdMatRestoreLoop original current hOriginal hCurrent
        stage 0 (hgcdMatStageOffset current hCurrent 0) staged.heap =
          .ok result := by
      simpa [hgcdMatStageOffset] using hrun
    have hLiveRestore := hgcdMatRestoreLoop_preserves_rawDenseRep this
      original current current hOriginal hCurrent hCurrent stage entries 0
      (by omega) staged.heap rfl hws1 hStage ptr length poly hOriginalPtr
      hLiveStage.2 result hrestore
    exact ⟨fun p n => (hLiveStage.1 p n).trans (hLiveRestore.1 p n),
      hLiveRestore.2⟩

/-- The source's alias-protection branch really preserves both iterator
outputs: `pB` is copied away first, the `pA` representation is framed across
that write, and only then is `a3` overwritten from `pA`. -/
theorem hgcdRecursiveStoreIterOutputs_cross_refines (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = true)
    (hA3 : heap.ValidU64Slice a3 lenA3)
    (hB3 : heap.ValidU64Slice b3 lenB3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB)
    (hB3PB : U64SlicesDisjoint b3 lenB3 pB lenB3)
    (hB3PA : U64SlicesDisjoint b3 lenB3 pA lenA3)
    (hA3PA : U64SlicesDisjoint a3 lenA3 pA lenA3)
    (hA3B3 : U64SlicesDisjoint a3 lenA3 b3 lenB3) :
    ∃ heap', hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
        .ok heap' ∧
      RawDensePolyRep this heap' a3 lenA3 polyA ∧
      RawDensePolyRep this heap' b3 lenB3 polyB := by
  rcases copyU64_refines_rawDense this heap b3 pB lenB3 polyB hB3 hB3PB
      hBRep with ⟨heap1, hcopyB, hlayout1, hB3Rep⟩
  have hARep1 := (copyU64_preserves_rawDenseRep this heap heap1 b3 pB lenB3
    pA lenA3 polyA hB3 hBRep.1 hB3PA hcopyB hARep).2
  have hA31 : heap1.ValidU64Slice a3 lenA3 :=
    (hlayout1 a3 lenA3).mp hA3
  rcases copyU64_refines_rawDense this heap1 a3 pA lenA3 polyA hA31 hA3PA
      hARep1 with ⟨heap2, hcopyA, _, hA3Rep⟩
  have hB3Rep2 := (copyU64_preserves_rawDenseRep this heap1 heap2 a3 pA lenA3
    b3 lenB3 polyB hA31 hARep1.1 hA3B3 hcopyA hB3Rep).2
  exact ⟨heap2,
    hgcdRecursiveStoreIterOutputs_cross_exec a3 b3 pA pB lenA3 lenB3
      heap heap1 heap2 hcross hcopyB hcopyA,
    hA3Rep, hB3Rep2⟩

/-- In the regular branch where both pointers need normalization, the two
source-order copies also preserve both represented iterator outputs. -/
theorem hgcdRecursiveStoreIterOutputs_regular_both_refines
    (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = false)
    (hPACopy : (pA == a3) = false) (hPBCopy : (pB == b3) = false)
    (hA3 : heap.ValidU64Slice a3 lenA3)
    (hB3 : heap.ValidU64Slice b3 lenB3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB)
    (hA3PA : U64SlicesDisjoint a3 lenA3 pA lenA3)
    (hA3PB : U64SlicesDisjoint a3 lenA3 pB lenB3)
    (hB3PB : U64SlicesDisjoint b3 lenB3 pB lenB3)
    (hB3A3 : U64SlicesDisjoint b3 lenB3 a3 lenA3) :
    ∃ heap', hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
        .ok heap' ∧
      RawDensePolyRep this heap' a3 lenA3 polyA ∧
      RawDensePolyRep this heap' b3 lenB3 polyB := by
  rcases copyU64_refines_rawDense this heap a3 pA lenA3 polyA hA3 hA3PA
      hARep with ⟨heap1, hcopyA, hlayout1, hA3Rep⟩
  have hBRep1 := (copyU64_preserves_rawDenseRep this heap heap1 a3 pA lenA3
    pB lenB3 polyB hA3 hARep.1 hA3PB hcopyA hBRep).2
  have hB31 : heap1.ValidU64Slice b3 lenB3 :=
    (hlayout1 b3 lenB3).mp hB3
  rcases copyU64_refines_rawDense this heap1 b3 pB lenB3 polyB hB31 hB3PB
      hBRep1 with ⟨heap2, hcopyB, _, hB3Rep⟩
  have hA3Rep2 := (copyU64_preserves_rawDenseRep this heap1 heap2 b3 pB lenB3
    a3 lenA3 polyA hB31 hBRep1.1 hB3A3 hcopyB hA3Rep).2
  have hfirst : (if (pA == a3) = false then heap.copyU64 a3 pA lenA3
      else .ok heap) = .ok heap1 := by
    simp [hPACopy, hcopyA]
  have hsecond : (if (pB == b3) = false then heap1.copyU64 b3 pB lenB3
      else .ok heap1) = .ok heap2 := by
    simp [hPBCopy, hcopyB]
  exact ⟨heap2,
    hgcdRecursiveStoreIterOutputs_regular_exec a3 b3 pA pB lenA3 lenB3
      heap heap1 heap2 hcross hfirst hsecond,
    hA3Rep2, hB3Rep⟩

theorem hgcdRecursiveStoreIterOutputs_regular_skip_a_refines
    (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = false)
    (hPASkip : (pA == a3) = true) (hpA : pA = a3)
    (hPBCopy : (pB == b3) = false)
    (hB3 : heap.ValidU64Slice b3 lenB3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB)
    (hB3PB : U64SlicesDisjoint b3 lenB3 pB lenB3)
    (hB3A3 : U64SlicesDisjoint b3 lenB3 a3 lenA3) :
    ∃ heap', hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
        .ok heap' ∧
      RawDensePolyRep this heap' a3 lenA3 polyA ∧
      RawDensePolyRep this heap' b3 lenB3 polyB := by
  have hA3Rep : RawDensePolyRep this heap a3 lenA3 polyA := by
    simpa [hpA] using hARep
  rcases copyU64_refines_rawDense this heap b3 pB lenB3 polyB hB3 hB3PB
      hBRep with ⟨heap1, hcopyB, _, hB3Rep⟩
  have hA3Rep1 := (copyU64_preserves_rawDenseRep this heap heap1 b3 pB lenB3
    a3 lenA3 polyA hB3 hBRep.1 hB3A3 hcopyB hA3Rep).2
  have hfirst : (if (pA == a3) = false then heap.copyU64 a3 pA lenA3
      else .ok heap) = .ok heap := by simp [hPASkip]
  have hsecond : (if (pB == b3) = false then heap.copyU64 b3 pB lenB3
      else .ok heap) = .ok heap1 := by simp [hPBCopy, hcopyB]
  exact ⟨heap1,
    hgcdRecursiveStoreIterOutputs_regular_exec a3 b3 pA pB lenA3 lenB3
      heap heap heap1 hcross hfirst hsecond,
    hA3Rep1, hB3Rep⟩

theorem hgcdRecursiveStoreIterOutputs_regular_skip_b_refines
    (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = false)
    (hPACopy : (pA == a3) = false)
    (hPBSkip : (pB == b3) = true) (hpB : pB = b3)
    (hA3 : heap.ValidU64Slice a3 lenA3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB)
    (hA3PA : U64SlicesDisjoint a3 lenA3 pA lenA3)
    (hA3B3 : U64SlicesDisjoint a3 lenA3 b3 lenB3) :
    ∃ heap', hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
        .ok heap' ∧
      RawDensePolyRep this heap' a3 lenA3 polyA ∧
      RawDensePolyRep this heap' b3 lenB3 polyB := by
  have hB3Rep : RawDensePolyRep this heap b3 lenB3 polyB := by
    simpa [hpB] using hBRep
  rcases copyU64_refines_rawDense this heap a3 pA lenA3 polyA hA3 hA3PA
      hARep with ⟨heap1, hcopyA, _, hA3Rep⟩
  have hB3Rep1 := (copyU64_preserves_rawDenseRep this heap heap1 a3 pA lenA3
    b3 lenB3 polyB hA3 hARep.1 hA3B3 hcopyA hB3Rep).2
  have hfirst : (if (pA == a3) = false then heap.copyU64 a3 pA lenA3
      else .ok heap) = .ok heap1 := by simp [hPACopy, hcopyA]
  have hsecond : (if (pB == b3) = false then heap1.copyU64 b3 pB lenB3
      else .ok heap1) = .ok heap1 := by simp [hPBSkip]
  exact ⟨heap1,
    hgcdRecursiveStoreIterOutputs_regular_exec a3 b3 pA pB lenA3 lenB3
      heap heap1 heap1 hcross hfirst hsecond,
    hA3Rep, hB3Rep1⟩

theorem hgcdRecursiveStoreIterOutputs_regular_skip_both_refines
    (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = false)
    (hPASkip : (pA == a3) = true) (hpA : pA = a3)
    (hPBSkip : (pB == b3) = true) (hpB : pB = b3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB) :
    hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap = .ok heap ∧
      RawDensePolyRep this heap a3 lenA3 polyA ∧
      RawDensePolyRep this heap b3 lenB3 polyB := by
  have hfirst : (if (pA == a3) = false then heap.copyU64 a3 pA lenA3
      else .ok heap) = .ok heap := by simp [hPASkip]
  have hsecond : (if (pB == b3) = false then heap.copyU64 b3 pB lenB3
      else .ok heap) = .ok heap := by simp [hPBSkip]
  exact ⟨hgcdRecursiveStoreIterOutputs_regular_exec a3 b3 pA pB lenA3 lenB3
      heap heap heap hcross hfirst hsecond,
    by simpa [hpA] using hARep, by simpa [hpB] using hBRep⟩

/-- Complete raw refinement of all pointer-comparison branches in the source
output-normalization block.  Equality bridges are facts about the concrete C++
pointer comparison; the remaining hypotheses are only capacities and aliasing. -/
theorem hgcdRecursiveStoreIterOutputs_refines (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (polyA polyB : Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hPAEq : (pA == a3) = true → pA = a3)
    (hPBEq : (pB == b3) = true → pB = b3)
    (hA3 : heap.ValidU64Slice a3 lenA3)
    (hB3 : heap.ValidU64Slice b3 lenB3)
    (hARep : RawDensePolyRep this heap pA lenA3 polyA)
    (hBRep : RawDensePolyRep this heap pB lenB3 polyB)
    (hB3PB : U64SlicesDisjoint b3 lenB3 pB lenB3)
    (hB3PA : U64SlicesDisjoint b3 lenB3 pA lenA3)
    (hA3PA : U64SlicesDisjoint a3 lenA3 pA lenA3)
    (hA3PB : U64SlicesDisjoint a3 lenA3 pB lenB3)
    (hA3B3 : U64SlicesDisjoint a3 lenA3 b3 lenB3) :
    ∃ heap', hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
        .ok heap' ∧
      RawDensePolyRep this heap' a3 lenA3 polyA ∧
      RawDensePolyRep this heap' b3 lenB3 polyB := by
  by_cases hcross : (!(pA == a3) && (pB == a3)) = true
  · exact hgcdRecursiveStoreIterOutputs_cross_refines this a3 b3 pA pB
      lenA3 lenB3 polyA polyB heap hcross hA3 hB3 hARep hBRep hB3PB
      hB3PA hA3PA hA3B3
  · have hcrossFalse : (!(pA == a3) && (pB == a3)) = false := by
      cases hvalue : (!(pA == a3) && (pB == a3)) <;> simp_all
    cases hpa : pA == a3 with
    | false =>
      cases hpb : pB == b3 with
      | false =>
        exact hgcdRecursiveStoreIterOutputs_regular_both_refines this a3 b3
          pA pB lenA3 lenB3 polyA polyB heap hcrossFalse hpa hpb hA3 hB3
          hARep hBRep hA3PA hA3PB hB3PB (u64SlicesDisjoint_symm hA3B3)
      | true =>
        exact hgcdRecursiveStoreIterOutputs_regular_skip_b_refines this a3 b3
          pA pB lenA3 lenB3 polyA polyB heap hcrossFalse hpa hpb (hPBEq hpb)
          hA3 hARep hBRep hA3PA hA3B3
    | true =>
      cases hpb : pB == b3 with
      | false =>
        exact hgcdRecursiveStoreIterOutputs_regular_skip_a_refines this a3 b3
          pA pB lenA3 lenB3 polyA polyB heap hcrossFalse hpa (hPAEq hpa) hpb
          hB3 hARep hBRep hB3PB (u64SlicesDisjoint_symm hA3B3)
      | true =>
        have hboth := hgcdRecursiveStoreIterOutputs_regular_skip_both_refines
          this a3 b3 pA pB lenA3 lenB3 polyA polyB heap hcrossFalse hpa
          (hPAEq hpa) hpb (hPBEq hpb) hARep hBRep
        exact ⟨heap, hboth.1, hboth.2.1, hboth.2.2⟩

/-- Matrix frame for one source-conditional memcpy. -/
theorem optionalCopy_preserves_hgcdMatRawDenseRep (this : DenseUPolyZp)
    (condition : Bool) (heap heap' : RawHeap)
    (dst src : RawPtr UInt64) (count : Nat)
    (M : HgcdMat) (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid) (hDst : heap.ValidU64Slice dst count)
    (hSrc : heap.ValidU64Slice src count)
    (hframe : ∀ j : Fin 4, U64SlicesDisjoint dst count
      (hgcdMatPtr M hM j) (hgcdMatLen M hM j))
    (hrep : HgcdMatRawDenseRep this heap M entries hM)
    (hrun : (if condition then heap.copyU64 dst src count else .ok heap) =
      .ok heap') :
    RawHeap.SameLayout heap heap' ∧
      HgcdMatRawDenseRep this heap' M entries hM := by
  cases hcondition : condition with
  | false =>
    have heq : heap = heap' := by simpa [hcondition] using hrun
    subst heap'
    exact ⟨fun _ _ => Iff.rfl, hrep⟩
  | true =>
    have hcopy : heap.copyU64 dst src count = .ok heap' := by
      simpa [hcondition] using hrun
    exact copyU64_preserves_hgcdMatRawDenseRep this heap heap' dst src count
      M entries hM hDst hSrc hframe hcopy hrep

/-- The complete alias-sensitive output normalization is a frame for every
stable matrix entry when both possible destinations are disjoint from it. -/
theorem hgcdRecursiveStoreIterOutputs_preserves_matrix
    (this : DenseUPolyZp)
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (heap heap' : RawHeap)
    (M : HgcdMat) (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hM : M.Valid)
    (hA3 : heap.ValidU64Slice a3 lenA3)
    (hB3 : heap.ValidU64Slice b3 lenB3)
    (hPA : heap.ValidU64Slice pA lenA3)
    (hPB : heap.ValidU64Slice pB lenB3)
    (hAFrame : ∀ j : Fin 4, U64SlicesDisjoint a3 lenA3
      (hgcdMatPtr M hM j) (hgcdMatLen M hM j))
    (hBFrame : ∀ j : Fin 4, U64SlicesDisjoint b3 lenB3
      (hgcdMatPtr M hM j) (hgcdMatLen M hM j))
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hrun : hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
      .ok heap') :
    RawHeap.SameLayout heap heap' ∧
      HgcdMatRawDenseRep this heap' M entries hM := by
  simp only [hgcdRecursiveStoreIterOutputs] at hrun
  split at hrun
  next hcross =>
    split at hrun
    next fault hcopyB => simp at hrun
    next heap1 hcopyB =>
      have hframe1 := copyU64_preserves_hgcdMatRawDenseRep this heap heap1
        b3 pB lenB3 M entries hM hB3 hPB hBFrame hcopyB hMatrix
      have hA31 := (hframe1.1 a3 lenA3).mp hA3
      have hPA1 := (hframe1.1 pA lenA3).mp hPA
      have hframe2 := copyU64_preserves_hgcdMatRawDenseRep this heap1 heap'
        a3 pA lenA3 M entries hM hA31 hPA1 hAFrame hrun hframe1.2
      exact ⟨fun ptr count =>
        (hframe1.1 ptr count).trans (hframe2.1 ptr count), hframe2.2⟩
  next hregular =>
    generalize hfirst :
      (if !(pA == a3) then heap.copyU64 a3 pA lenA3 else .ok heap) = first
      at hrun
    cases first with
    | error fault => simp at hrun
    | ok heap1 =>
      have hframe1 := optionalCopy_preserves_hgcdMatRawDenseRep this
        (!(pA == a3)) heap heap1 a3 pA lenA3 M entries hM hA3 hPA
        hAFrame hMatrix hfirst
      have hB31 := (hframe1.1 b3 lenB3).mp hB3
      have hPB1 := (hframe1.1 pB lenB3).mp hPB
      have hsecond :
          (if !(pB == b3) then heap1.copyU64 b3 pB lenB3 else .ok heap1) =
            .ok heap' := by
        exact hrun
      have hframe2 := optionalCopy_preserves_hgcdMatRawDenseRep this
        (!(pB == b3)) heap1 heap' b3 pB lenB3 M entries hM hB31 hPB1
        hBFrame hframe1.2 hsecond
      exact ⟨fun ptr count =>
        (hframe1.1 ptr count).trans (hframe2.1 ptr count), hframe2.2⟩

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

/-- Physical buffer obligations for one guarded multiplication in the
low-half reconstruction of `_hgcd_recursive`.  The maximum operand length
matches the source's longer-first dispatch and does not encode an L2 value. -/
structure HgcdMulTermWorkspace (heap : RawHeap)
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat)
    (scratch : RawPtr UInt64) : Prop where
  lengthFits : max lenLeft lenRight < limbBase
  dstValid : heap.ValidU64Slice dst (2 * max lenLeft lenRight - 1)
  scratchValid : heap.ValidU64Slice scratch (8 * max lenLeft lenRight)
  dstLeft : U64SlicesDisjoint dst (2 * max lenLeft lenRight - 1)
    left lenLeft
  dstRight : U64SlicesDisjoint dst (2 * max lenLeft lenRight - 1)
    right lenRight
  dstScratch : U64SlicesDisjoint dst (2 * max lenLeft lenRight - 1)
    scratch (8 * max lenLeft lenRight)
  scratchLeft : U64SlicesDisjoint scratch (8 * max lenLeft lenRight)
    left lenLeft
  scratchRight : U64SlicesDisjoint scratch (8 * max lenLeft lenRight)
    right lenRight

/-- The exact guarded C++ multiplication block computes the polynomial
product.  In particular, a zero-length branch is justified directly from
the raw input representation. -/
theorem hgcdRecursiveMulTerm_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat)
    (scratch : RawPtr UInt64) (heap : RawHeap)
    (leftPoly rightPoly : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hwork : HgcdMulTermWorkspace heap dst left lenLeft right lenRight scratch)
    (hLeft : RawDensePolyRep this heap left lenLeft leftPoly)
    (hRight : RawDensePolyRep this heap right lenRight rightPoly) :
    ∃ result, hgcdRecursiveMulTerm this dst left lenLeft right lenRight
        scratch heap = .ok result ∧
      RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap dst result.length
        (leftPoly * rightPoly) := by
  by_cases hLeftPos : 0 < lenLeft
  · by_cases hRightPos : 0 < lenRight
    · by_cases horder : lenRight ≤ lenLeft
      · have hmax : max lenLeft lenRight = lenLeft := Nat.max_eq_left horder
        rcases mul_refines_rawDense this dst left lenLeft right lenRight
            scratch heap leftPoly rightPoly hcfg hp hLeftPos hRightPos horder
            (by simpa [hmax] using hwork.lengthFits)
            (by simpa [hmax] using hwork.dstValid)
            (by simpa [hmax] using hwork.scratchValid)
            (by simpa [hmax] using hwork.dstLeft)
            (by simpa [hmax] using hwork.dstRight)
            (by simpa [hmax] using hwork.dstScratch)
            (by simpa [hmax] using hwork.scratchLeft)
            (by simpa [hmax] using hwork.scratchRight)
            hLeft hRight with ⟨heap1, hrun, hlayout, hrep⟩
        refine ⟨{ heap := heap1, length := lenLeft + lenRight - 1 }, ?_,
          hlayout, hrep⟩
        simp [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hrun]
      · have hreverse : lenLeft ≤ lenRight := by omega
        have hmax : max lenLeft lenRight = lenRight := Nat.max_eq_right hreverse
        rcases mul_refines_rawDense this dst right lenRight left lenLeft
            scratch heap rightPoly leftPoly hcfg hp hRightPos hLeftPos hreverse
            (by simpa [hmax] using hwork.lengthFits)
            (by simpa [hmax] using hwork.dstValid)
            (by simpa [hmax] using hwork.scratchValid)
            (by simpa [hmax] using hwork.dstRight)
            (by simpa [hmax] using hwork.dstLeft)
            (by simpa [hmax] using hwork.dstScratch)
            (by simpa [hmax] using hwork.scratchRight)
            (by simpa [hmax] using hwork.scratchLeft)
            hRight hLeft with ⟨heap1, hrun, hlayout, hrep⟩
        refine ⟨{ heap := heap1, length := lenLeft + lenRight - 1 }, ?_,
          hlayout, ?_⟩
        · simp [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hrun,
            Nat.add_comm]
        · simpa [Nat.add_comm, mul_comm] using hrep
    · have hlenRight : lenRight = 0 := by omega
      subst lenRight
      have hzero : rightPoly = 0 :=
        slicePolyRep_zero_length heap right this._p.toNat rightPoly
          hRight.2.2.1
      subst rightPoly
      refine ⟨{ heap := heap, length := 0 }, ?_, (fun _ _ => Iff.rfl),
        ?_⟩
      · simp [hgcdRecursiveMulTerm]
      · simpa using rawDensePolyRep_zero_length this heap dst
          (heap.validU64Slice_mono dst (2 * max lenLeft 0 - 1) 0
            hwork.dstValid (by omega))
  · have hlenLeft : lenLeft = 0 := by omega
    subst lenLeft
    have hzero : leftPoly = 0 :=
      slicePolyRep_zero_length heap left this._p.toNat leftPoly
        hLeft.2.2.1
    subst leftPoly
    refine ⟨{ heap := heap, length := 0 }, ?_, (fun _ _ => Iff.rfl),
      ?_⟩
    · simp [hgcdRecursiveMulTerm]
    · simpa using rawDensePolyRep_zero_length this heap dst
        (heap.validU64Slice_mono dst (2 * max 0 lenRight - 1) 0
          hwork.dstValid (by omega))

/-- Frame rule for the exact guarded multiplication block.  It exposes that
the generated execution writes only its destination and scratch areas, which
is needed to sequence the two products in each reconstruction block. -/
theorem hgcdRecursiveMulTerm_preserves_guard (this : DenseUPolyZp)
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat)
    (scratch guard : RawPtr UInt64) (guardLen : Nat)
    (heap : RawHeap) (result : HgcdMulTermResult)
    (hwork : HgcdMulTermWorkspace heap dst left lenLeft right lenRight scratch)
    (hLeft : heap.ValidU64Slice left lenLeft)
    (hRight : heap.ValidU64Slice right lenRight)
    (hDstGuard : U64SlicesDisjoint dst (2 * max lenLeft lenRight - 1)
      guard guardLen)
    (hScratchGuard : U64SlicesDisjoint scratch
      (8 * max lenLeft lenRight) guard guardLen)
    (hrun : hgcdRecursiveMulTerm this dst left lenLeft right lenRight
      scratch heap = .ok result) :
    SameU64Prefix heap result.heap guard guardLen := by
  by_cases hLeftPos : 0 < lenLeft
  · by_cases hRightPos : 0 < lenRight
    · by_cases horder : lenRight ≤ lenLeft
      · have hmax : max lenLeft lenRight = lenLeft := Nat.max_eq_left horder
        cases hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this dst left lenLeft right
            lenRight scratch heap with
        | error fault =>
          simp [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hmul]
            at hrun
        | ok heap1 =>
          have hactual :
              (Except.ok (HgcdMulTermResult.mk heap1
                (lenLeft + lenRight - 1)) : RawExec HgcdMulTermResult) =
              Except.ok result := by
            simpa [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hmul]
              using hrun
          have heq : result =
              HgcdMulTermResult.mk heap1 (lenLeft + lenRight - 1) :=
            (Except.ok.inj hactual).symm
          subst result
          exact mul_preserves_prefix this dst left lenLeft right lenRight
            scratch guard guardLen heap heap1 hLeftPos hRightPos horder
            (by simpa [hmax] using hwork.dstValid) hLeft hRight
            (by simpa [hmax] using hwork.scratchValid)
            (by simpa [hmax] using hwork.scratchRight)
            (by simpa [hmax] using hDstGuard)
            (by simpa [hmax] using hScratchGuard) hmul
      · have hreverse : lenLeft ≤ lenRight := by omega
        have hmax : max lenLeft lenRight = lenRight := Nat.max_eq_right hreverse
        cases hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this dst right lenRight left
            lenLeft scratch heap with
        | error fault =>
          simp [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hmul]
            at hrun
        | ok heap1 =>
          have hactual :
              (Except.ok (HgcdMulTermResult.mk heap1
                (lenLeft + lenRight - 1)) : RawExec HgcdMulTermResult) =
              Except.ok result := by
            simpa [hgcdRecursiveMulTerm, hLeftPos, hRightPos, horder, hmul]
              using hrun
          have heq : result =
              HgcdMulTermResult.mk heap1 (lenLeft + lenRight - 1) :=
            (Except.ok.inj hactual).symm
          subst result
          exact mul_preserves_prefix this dst right lenRight left lenLeft
            scratch guard guardLen heap heap1 hRightPos hLeftPos hreverse
            (by simpa [hmax] using hwork.dstValid) hRight hLeft
            (by simpa [hmax] using hwork.scratchValid)
            (by simpa [hmax] using hwork.scratchLeft)
            (by simpa [hmax] using hDstGuard)
            (by simpa [hmax] using hScratchGuard) hmul
    · have hlenRight : lenRight = 0 := by omega
      subst lenRight
      have hactual :
          (Except.ok (HgcdMulTermResult.mk heap 0) :
            RawExec HgcdMulTermResult) = Except.ok result := by
        simpa [hgcdRecursiveMulTerm] using hrun
      have heq : result = HgcdMulTermResult.mk heap 0 :=
        (Except.ok.inj hactual).symm
      subst result
      exact fun _ _ _ hread => hread
  · have hlenLeft : lenLeft = 0 := by omega
    subst lenLeft
    have hactual :
        (Except.ok (HgcdMulTermResult.mk heap 0) :
          RawExec HgcdMulTermResult) = Except.ok result := by
      simpa [hgcdRecursiveMulTerm] using hrun
    have heq : result = HgcdMulTermResult.mk heap 0 :=
      (Except.ok.inj hactual).symm
    subst result
    exact fun _ _ _ hread => hread

def hgcdMulCapacity (leftLength rightLength : Nat) : Nat :=
  2 * max leftLength rightLength - 1

/-- Physical obligations for the two generated products and their in-place
subtraction.  It contains capacities and aliasing only, never an expected L2
polynomial or a supplied execution result. -/
structure HgcdReconstructWorkspace (heap : RawHeap)
    (out temp left1 : RawPtr UInt64) (lenLeft1 : Nat)
    (right1 : RawPtr UInt64) (lenRight1 : Nat)
    (left2 : RawPtr UInt64) (lenLeft2 : Nat)
    (right2 : RawPtr UInt64) (lenRight2 : Nat)
    (scratch : RawPtr UInt64) : Prop where
  first : HgcdMulTermWorkspace heap out left1 lenLeft1 right1 lenRight1 scratch
  second : HgcdMulTermWorkspace heap temp left2 lenLeft2 right2 lenRight2 scratch
  subDstValid : heap.ValidU64Slice out
    (max (hgcdMulCapacity lenLeft1 lenRight1)
      (hgcdMulCapacity lenLeft2 lenRight2))
  firstDstLeft2 : U64SlicesDisjoint out
    (hgcdMulCapacity lenLeft1 lenRight1) left2 lenLeft2
  firstScratchLeft2 : U64SlicesDisjoint scratch
    (8 * max lenLeft1 lenRight1) left2 lenLeft2
  firstDstRight2 : U64SlicesDisjoint out
    (hgcdMulCapacity lenLeft1 lenRight1) right2 lenRight2
  firstScratchRight2 : U64SlicesDisjoint scratch
    (8 * max lenLeft1 lenRight1) right2 lenRight2
  secondDstFirst : U64SlicesDisjoint temp
    (hgcdMulCapacity lenLeft2 lenRight2) out
    (hgcdMulCapacity lenLeft1 lenRight1)
  secondScratchFirst : U64SlicesDisjoint scratch
    (8 * max lenLeft2 lenRight2) out
    (hgcdMulCapacity lenLeft1 lenRight1)
  subAliasTemp : ExactOrDisjoint out temp

theorem hgcdMulTermWorkspace_of_sameLayout (heap heap' : RawHeap)
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat) (scratch : RawPtr UInt64)
    (hlayout : RawHeap.SameLayout heap heap')
    (hwork : HgcdMulTermWorkspace heap dst left lenLeft right lenRight scratch) :
    HgcdMulTermWorkspace heap' dst left lenLeft right lenRight scratch := by
  exact ⟨hwork.lengthFits,
    (hlayout dst (2 * max lenLeft lenRight - 1)).mp hwork.dstValid,
    (hlayout scratch (8 * max lenLeft lenRight)).mp hwork.scratchValid,
    hwork.dstLeft, hwork.dstRight, hwork.dstScratch,
    hwork.scratchLeft, hwork.scratchRight⟩

/-- Semantic refinement of the exact `b2` low reconstruction block.  Both
products and the final sign-selected subtraction are actual generated raw
executions. -/
theorem hgcdRecursiveReconstructB_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (b2 T0 r2 r0 aLo bLo scratch : RawPtr UInt64)
    (lenR2 lenR0 lenALo lenBLo : Nat) (sgn : Int) (heap : RawHeap)
    (polyR2 polyR0 polyALo polyBLo :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hwork : HgcdReconstructWorkspace heap b2 T0 r2 lenR2 aLo lenALo
      r0 lenR0 bLo lenBLo scratch)
    (hR2 : RawDensePolyRep this heap r2 lenR2 polyR2)
    (hR0 : RawDensePolyRep this heap r0 lenR0 polyR0)
    (hALo : RawDensePolyRep this heap aLo lenALo polyALo)
    (hBLo : RawDensePolyRep this heap bLo lenBLo polyBLo) :
    ∃ heap' length,
      hgcdRecursiveReconstructB this b2 T0 r2 r0 aLo bLo scratch
        lenR2 lenR0 lenALo lenBLo sgn heap = .ok (heap', length) ∧
      RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' b2 length
        (if sgn < 0 then polyR2 * polyALo - polyR0 * polyBLo
         else polyR0 * polyBLo - polyR2 * polyALo) ∧
      length ≤ max (lenR2 + lenALo) (lenR0 + lenBLo) := by
  rcases hgcdRecursiveMulTerm_refines this b2 r2 lenR2 aLo lenALo
      scratch heap polyR2 polyALo hcfg hp hwork.first hR2 hALo with
    ⟨term1, hrun1, hlayout1, hTerm1⟩
  have hR01 : RawDensePolyRep this term1.heap r0 lenR0 polyR0 := by
    apply rawDensePolyRep_of_same_prefix this heap term1.heap r0 lenR0
      polyR0 hlayout1
    · exact hgcdRecursiveMulTerm_preserves_guard this b2 r2 lenR2 aLo
        lenALo scratch r0 lenR0 heap term1 hwork.first hR2.1 hALo.1
        hwork.firstDstLeft2 hwork.firstScratchLeft2 hrun1
    · exact hR0
  have hBLo1 : RawDensePolyRep this term1.heap bLo lenBLo polyBLo := by
    apply rawDensePolyRep_of_same_prefix this heap term1.heap bLo lenBLo
      polyBLo hlayout1
    · exact hgcdRecursiveMulTerm_preserves_guard this b2 r2 lenR2 aLo
        lenALo scratch bLo lenBLo heap term1 hwork.first hR2.1 hALo.1
        hwork.firstDstRight2 hwork.firstScratchRight2 hrun1
    · exact hBLo
  have hSecond1 := hgcdMulTermWorkspace_of_sameLayout heap term1.heap T0 r0
    lenR0 bLo lenBLo scratch hlayout1 hwork.second
  rcases hgcdRecursiveMulTerm_refines this T0 r0 lenR0 bLo lenBLo
      scratch term1.heap polyR0 polyBLo hcfg hp hSecond1 hR01 hBLo1 with
    ⟨term2, hrun2, hlayout2, hTerm2⟩
  have hLen1 := hgcdRecursiveMulTerm_length_le this b2 r2 lenR2 aLo
    lenALo scratch heap term1 hrun1
  have hLen2 := hgcdRecursiveMulTerm_length_le this T0 r0 lenR0 bLo
    lenBLo scratch term1.heap term2 hrun2
  have hLen1Sum := hgcdRecursiveMulTerm_length_le_sum this b2 r2 lenR2 aLo
    lenALo scratch heap term1 hrun1
  have hLen2Sum := hgcdRecursiveMulTerm_length_le_sum this T0 r0 lenR0 bLo
    lenBLo scratch term1.heap term2 hrun2
  have hTerm1Final : RawDensePolyRep this term2.heap b2 term1.length
      (polyR2 * polyALo) := by
    apply rawDensePolyRep_of_same_prefix this term1.heap term2.heap b2
      term1.length (polyR2 * polyALo) hlayout2
    · apply hgcdRecursiveMulTerm_preserves_guard this T0 r0 lenR0 bLo
        lenBLo scratch b2 term1.length term1.heap term2 hSecond1 hR01.1
        hBLo1.1
      · exact u64SlicesDisjoint_mono hwork.secondDstFirst
          (Nat.le_refl _) hLen1
      · exact u64SlicesDisjoint_mono hwork.secondScratchFirst
          (Nat.le_refl _) hLen1
      · exact hrun2
    · exact hTerm1
  have hSubValid0 : heap.ValidU64Slice b2
      (max term1.length term2.length) :=
    heap.validU64Slice_mono b2
      (max (hgcdMulCapacity lenR2 lenALo)
        (hgcdMulCapacity lenR0 lenBLo))
      (max term1.length term2.length) hwork.subDstValid (by
        apply max_le
        · exact le_max_of_le_left hLen1
        · exact le_max_of_le_right hLen2)
  have hSubValid : term2.heap.ValidU64Slice b2
      (max term1.length term2.length) :=
    (hlayout2 b2 _).mp ((hlayout1 b2 _).mp hSubValid0)
  have hpWord : this._p ≠ 0 := by
    intro hzero
    have hzeroNat := congrArg UInt64.toNat hzero
    simp at hzeroNat
    omega
  by_cases hsgn : sgn < 0
  · rcases polySub_ok this b2 b2 term1.length T0 term2.length
        term2.heap hSubValid hTerm1Final.1 hTerm2.1 with
      ⟨heap3, length, hsub, hlayout3, hlength⟩
    have hrep := polySub_refines this b2 b2 term1.length T0 term2.length
      term2.heap heap3 length (polyR2 * polyALo) (polyR0 * polyBLo)
      hpWord hSubValid hTerm1Final hTerm2 (Or.inl rfl)
      hwork.subAliasTemp hsub
    refine ⟨heap3, length, ?_, ?_, ?_, ?_⟩
    · simp [hgcdRecursiveReconstructB, hrun1, hrun2, hsgn, hsub]
    · exact fun ptr count =>
        (hlayout1 ptr count).trans
          ((hlayout2 ptr count).trans (hlayout3 ptr count))
    · simpa [hsgn] using hrep
    · exact hlength.trans (max_le_max hLen1Sum hLen2Sum)

  · rcases polySub_ok this b2 T0 term2.length b2 term1.length
        term2.heap (by simpa [max_comm] using hSubValid) hTerm2.1
        hTerm1Final.1 with
      ⟨heap3, length, hsub, hlayout3, hlength⟩
    have hrep := polySub_refines this b2 T0 term2.length b2 term1.length
      term2.heap heap3 length (polyR0 * polyBLo) (polyR2 * polyALo)
      hpWord (by simpa [max_comm] using hSubValid) hTerm2 hTerm1Final
      hwork.subAliasTemp
      (Or.inl rfl) hsub
    refine ⟨heap3, length, ?_, ?_, ?_, ?_⟩
    · simp [hgcdRecursiveReconstructB, hrun1, hrun2, hsgn, hsub]
    · exact fun ptr count =>
        (hlayout1 ptr count).trans
          ((hlayout2 ptr count).trans (hlayout3 ptr count))
    · simpa [hsgn] using hrep
    · simpa [max_comm] using hlength.trans
        (max_le_max hLen2Sum hLen1Sum)

/-- Semantic refinement of the exact `a2` reconstruction block.  Its source
differs from `b2` only by reversing the sign-selected subtraction, so the
already verified physical execution is reused with an explicitly flipped
branch selector. -/
theorem hgcdRecursiveReconstructA_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (a2 T0 r3 r1 aLo bLo scratch : RawPtr UInt64)
    (lenR3 lenR1 lenALo lenBLo : Nat) (sgn : Int) (heap : RawHeap)
    (polyR3 polyR1 polyALo polyBLo :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hwork : HgcdReconstructWorkspace heap a2 T0 r3 lenR3 aLo lenALo
      r1 lenR1 bLo lenBLo scratch)
    (hR3 : RawDensePolyRep this heap r3 lenR3 polyR3)
    (hR1 : RawDensePolyRep this heap r1 lenR1 polyR1)
    (hALo : RawDensePolyRep this heap aLo lenALo polyALo)
    (hBLo : RawDensePolyRep this heap bLo lenBLo polyBLo) :
    ∃ heap' length,
      hgcdRecursiveReconstructA this a2 T0 r3 r1 aLo bLo scratch
        lenR3 lenR1 lenALo lenBLo sgn heap = .ok (heap', length) ∧
      RawHeap.SameLayout heap heap' ∧
      RawDensePolyRep this heap' a2 length
        (if sgn < 0 then polyR1 * polyBLo - polyR3 * polyALo
         else polyR3 * polyALo - polyR1 * polyBLo) ∧
      length ≤ max (lenR3 + lenALo) (lenR1 + lenBLo) := by
  let flippedSign : Int := if sgn < 0 then 0 else -1
  rcases hgcdRecursiveReconstructB_refines this a2 T0 r3 r1 aLo bLo
      scratch lenR3 lenR1 lenALo lenBLo flippedSign heap polyR3 polyR1
      polyALo polyBLo hcfg hp hwork hR3 hR1 hALo hBLo with
    ⟨heap', length, hrun, hlayout, hrep, hlength⟩
  refine ⟨heap', length, ?_, hlayout, ?_, hlength⟩
  · have hfunctions :
        hgcdRecursiveReconstructA this a2 T0 r3 r1 aLo bLo scratch
            lenR3 lenR1 lenALo lenBLo sgn heap =
          hgcdRecursiveReconstructB this a2 T0 r3 r1 aLo bLo scratch
            lenR3 lenR1 lenALo lenBLo flippedSign heap := by
      by_cases hsgn : sgn < 0
      · simp [hgcdRecursiveReconstructA, hgcdRecursiveReconstructB,
          flippedSign, hsgn]
      · simp [hgcdRecursiveReconstructA, hgcdRecursiveReconstructB,
          flippedSign, hsgn]
    exact hfunctions.trans hrun
  · by_cases hsgn : sgn < 0
    · simpa [flippedSign, hsgn] using hrep
    · simpa [flippedSign, hsgn] using hrep

/-- Purely physical buffers for the generated shifted-high reconstruction.
The output covers both the existing low representation and the inserted high
half; the alias condition is the one accepted by the real in-place add. -/
structure HgcdLiftHighWorkspace (heap : RawHeap)
    (out high : RawPtr UInt64) (lowLength m highLength : Nat) : Prop where
  outValid : heap.ValidU64Slice out (max (m + highLength) lowLength)
  highValid : heap.ValidU64Slice high highLength
  zeroHighDisjoint : U64SlicesDisjoint out
    (max (m + highLength) lowLength) high highLength
  addAliasHigh : ExactOrDisjoint (out.add m) high

/-- Splitting a normalized raw polynomial at an in-range offset leaves a
normalized suffix.  This is what permits the source's in-place `_poly_add`
on `out + m` without inventing a normalized padded whole buffer. -/
theorem rawDensePolyRep_split_suffix (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (length m : Nat)
    (poly : Polynomial (ZMod this._p.toNat)) (hm : m ≤ length)
    (hrep : RawDensePolyRep this heap ptr length poly) :
    ∃ low high : Polynomial (ZMod this._p.toNat),
      SlicePolyRep heap ptr m this._p.toNat low ∧
      RawDensePolyRep this heap (ptr.add m) (length - m) high ∧
      poly = low + Polynomial.X ^ m * high := by
  rcases slicePolyRep_split_exists heap ptr m (length - m) this._p.toNat
      poly (by simpa [Nat.add_sub_of_le hm] using hrep.1)
      (by simpa [Nat.add_sub_of_le hm] using hrep.2.2.1) with
    ⟨low, high, hLow, hHigh, hsplit⟩
  have hcanonical := canonicalU64Prefix_split heap ptr m (length - m)
    this._p (by simpa [Nat.add_sub_of_le hm] using hrep.2.1)
  have hvalidHigh := heap.validU64Slice_add ptr length m (length - m)
    hrep.1 (by omega)
  by_cases hstrict : m < length
  · have hnormHigh : heap.normaliseU64 (ptr.add m) (length - m) =
        .ok (length - m) := by
      rcases normaliseU64_spec heap ptr length hrep.1 with
        ⟨observed, hobserved, _, _, hlast⟩
      have heq : observed = length :=
        Except.ok.inj (hobserved.symm.trans hrep.2.2.2)
      subst observed
      rcases hlast with hzero | ⟨value, hread, hne⟩
      · omega
      · have hreadHigh : heap.readU64 (ptr.add m) (length - m - 1) =
            .ok value := by
          rw [RawHeap.readU64_add]
          have hindex : m + (length - m - 1) = length - 1 := by omega
          simpa [hindex] using hread
        rw [show length - m = (length - m - 1) + 1 by omega]
        simp [RawHeap.normaliseU64, hreadHigh, hne]
    exact ⟨low, high, hLow,
      ⟨hvalidHigh, hcanonical.2, hHigh, hnormHigh⟩, hsplit⟩
  · have heq : m = length := by omega
    subst m
    have hhighZero : high = 0 :=
      slicePolyRep_zero_length heap (ptr.add length) this._p.toNat high
        (by simpa using hHigh)
    subst high
    have hzeroRep : RawDensePolyRep this heap (ptr.add length)
        (length - length) 0 := by
      simpa using rawDensePolyRep_zero_length this heap (ptr.add length)
        (by simpa using hvalidHigh)
    exact ⟨low, 0, hLow, hzeroRep, by simpa using hsplit⟩

/-- Semantic result of the conditional source `memset`: the original
normalized prefix is framed, while the enlarged physical prefix represents
the same polynomial with canonical zero limbs appended. -/
theorem hgcdLiftHigh_zero_refines (this : DenseUPolyZp)
    (out : RawPtr UInt64) (lowLength m highLength : Nat)
    (heap : RawHeap) (lowPoly : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hOut : heap.ValidU64Slice out (max (m + highLength) lowLength))
    (hLow : RawDensePolyRep this heap out lowLength lowPoly) :
    ∃ heap1,
      (if lowLength < m + highLength then
        Generated.StrictMul.mulZeroPadLoop out lowLength
          (m + highLength - lowLength) 0 heap
       else (Except.ok heap : RawExec RawHeap)) = Except.ok heap1 ∧
      RawHeap.SameLayout heap heap1 ∧
      RawDensePolyRep this heap1 out lowLength lowPoly ∧
      SlicePolyRep heap1 out (max (m + highLength) lowLength)
        this._p.toNat lowPoly ∧
      CanonicalU64Prefix heap1 out (max (m + highLength) lowLength)
        this._p := by
  by_cases hpad : lowLength < m + highLength
  · have hsum : lowLength + (m + highLength - lowLength) =
        m + highLength := by omega
    have hmax : max (m + highLength) lowLength = m + highLength :=
      Nat.max_eq_left (Nat.le_of_lt hpad)
    have hvalid : heap.ValidU64Slice out
        (lowLength + (m + highLength - lowLength)) := by
      simpa [hsum, hmax] using hOut
    rcases mulZeroPadLoop_refines out lowLength
        (m + highLength - lowLength) 0 this._p.toNat heap lowPoly this._p
        (Nat.zero_le _) hp hvalid (by simpa using hLow.2.2.1)
        (by simpa using hLow.2.1) with
      ⟨heap1, hrun, hlayout, hslice, hcanonical⟩
    have hsame := mulZeroPadLoop_preserves_before_start out lowLength
      (m + highLength - lowLength) 0 lowLength heap heap1
      (Nat.le_refl _) hvalid hrun
    have hLow1 := rawDensePolyRep_of_same_prefix this heap heap1 out
      lowLength lowPoly hlayout hsame hLow
    exact ⟨heap1, by simpa [hpad] using hrun, hlayout, hLow1,
      by simpa [hsum, hmax] using hslice,
      by simpa [hsum, hmax] using hcanonical⟩
  · have hmax : max (m + highLength) lowLength = lowLength :=
      Nat.max_eq_right (by omega)
    exact ⟨heap, by simp [hpad], fun _ _ => Iff.rfl, hLow,
      by simpa [hmax] using hLow.2.2.1,
      by simpa [hmax] using hLow.2.1⟩

/-- Re-expand a normalized raw result to the full prefix scanned by the
actual normalization call.  The discarded suffix consists of physical zero
limbs, so both the polynomial and canonical-residue properties extend. -/
theorem rawDensePolyRep_extend_to_normalise_input (this : DenseUPolyZp)
    (heap : RawHeap) (ptr : RawPtr UInt64) (fullLength resultLength : Nat)
    (poly : Polynomial (ZMod this._p.toNat)) (hp : this._p ≠ 0)
    (hFull : heap.ValidU64Slice ptr fullLength)
    (hResult : RawDensePolyRep this heap ptr resultLength poly)
    (hnorm : heap.normaliseU64 ptr fullLength = .ok resultLength) :
    SlicePolyRep heap ptr fullLength this._p.toNat poly ∧
      CanonicalU64Prefix heap ptr fullLength this._p := by
  rcases normaliseU64_spec heap ptr fullLength hFull with
    ⟨observed, hobserved, hle, hzeros, _⟩
  have heq : observed = resultLength :=
    Except.ok.inj (hobserved.symm.trans hnorm)
  subst observed
  rcases slicePolyRep_exists_unique heap ptr fullLength this._p.toNat hFull with
    ⟨actual, hActual, _⟩
  have hpoly : actual = poly := by
    ext degree
    by_cases hshort : degree < resultLength
    · rcases slicePolyRep_coeff heap ptr fullLength this._p.toNat actual
          hActual degree (by omega) with ⟨v1, hr1, hc1⟩
      rcases slicePolyRep_coeff heap ptr resultLength this._p.toNat poly
          hResult.2.2.1 degree hshort with ⟨v2, hr2, hc2⟩
      have hv : v1 = v2 := Except.ok.inj (hr1.symm.trans hr2)
      simpa [hc1, hc2, hv]
    · have hpolyZero := slicePolyRep_coeff_zero_of_length_le heap ptr
          resultLength this._p.toNat poly hResult.2.2.1 degree (by omega)
      by_cases hfull : degree < fullLength
      · rcases slicePolyRep_coeff heap ptr fullLength this._p.toNat actual
            hActual degree hfull with ⟨value, hread, hcoeff⟩
        have hreadZero := hzeros degree (by omega) hfull
        have hv : value = 0 := Except.ok.inj (hread.symm.trans hreadZero)
        subst value
        simpa [hcoeff] using hpolyZero.symm
      · rw [slicePolyRep_coeff_zero_of_length_le heap ptr fullLength
            this._p.toNat actual hActual degree (by omega), hpolyZero]
  have hSlice : SlicePolyRep heap ptr fullLength this._p.toNat poly := by
    simpa [hpoly] using hActual
  refine ⟨hSlice, ?_⟩
  intro k value hk hread
  by_cases hshort : k < resultLength
  · exact hResult.2.1 k value hshort hread
  · have hreadZero := hzeros k (by omega) hk
    have hv : value = 0 := Except.ok.inj (hread.symm.trans hreadZero)
    subst value
    have hpPos : 0 < this._p.toNat := UInt64.pos_iff_ne_zero.mpr hp
    simpa using hpPos

/-- Total raw execution bridge for zero-fill, shifted in-place addition, and
whole-buffer normalization.  Every step is the corresponding generated L1
operation; no expected L2 output is an input to the execution. -/
theorem hgcdRecursiveLiftHigh_terminates (this : DenseUPolyZp)
    (out high : RawPtr UInt64) (lowLength m highLength : Nat)
    (heap : RawHeap) (lowPoly : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hwork : HgcdLiftHighWorkspace heap out high lowLength m highLength)
    (hLow : RawDensePolyRep this heap out lowLength lowPoly) :
    ∃ result,
      hgcdRecursiveLiftHigh this out high lowLength m highLength heap =
        .ok result ∧
      RawHeap.SameLayout heap result.heap ∧
      result.length ≤ max (m + highLength) lowLength := by
  let required := m + highLength
  let fullLength := max required lowLength
  have hzero : ∃ heap1,
      (if lowLength < required then
        Generated.StrictMul.mulZeroPadLoop out lowLength
          (required - lowLength) 0 heap
       else (Except.ok heap : RawExec RawHeap)) = Except.ok heap1 ∧
      RawHeap.SameLayout heap heap1 := by
    by_cases hpad : lowLength < required
    · have hsum : lowLength + (required - lowLength) = required := by omega
      have hfull : fullLength = required := by
        simp [fullLength, Nat.max_eq_left (Nat.le_of_lt hpad)]
      have hvalidFull : heap.ValidU64Slice out fullLength := by
        simpa [fullLength, required] using hwork.outValid
      have hvalid : heap.ValidU64Slice out
          (lowLength + (required - lowLength)) := by
        simpa [hsum, hfull] using hvalidFull
      rcases mulZeroPadLoop_refines out lowLength (required - lowLength) 0
          this._p.toNat heap lowPoly this._p (Nat.zero_le _) hp hvalid
          (by simpa using hLow.2.2.1) (by simpa using hLow.2.1) with
        ⟨heap1, hrun, hlayout, _, _⟩
      exact ⟨heap1, by simpa [hpad] using hrun, hlayout⟩
    · exact ⟨heap, by simp [hpad], fun _ _ => Iff.rfl⟩
  rcases hzero with ⟨heap1, hzeroRun, hlayout1⟩
  let oldHighLength := if m ≤ lowLength then lowLength - m else 0
  have hOut1 : heap1.ValidU64Slice out fullLength :=
    (hlayout1 out fullLength).mp (by simpa [fullLength, required] using
      hwork.outValid)
  have hOutHigh1 : heap1.ValidU64Slice (out.add m) oldHighLength := by
    apply heap1.validU64Slice_add out fullLength m oldHighLength hOut1
    dsimp [oldHighLength, fullLength, required]
    split <;> omega
  have hHigh1 : heap1.ValidU64Slice high highLength :=
    (hlayout1 high highLength).mp hwork.highValid
  have hAddOut1 : heap1.ValidU64Slice (out.add m)
      (max oldHighLength highLength) := by
    apply heap1.validU64Slice_add out fullLength m
      (max oldHighLength highLength) hOut1
    dsimp [oldHighLength, fullLength, required]
    split <;> omega
  rcases polyAdd_ok this (out.add m) (out.add m) oldHighLength high
      highLength heap1 hAddOut1 hOutHigh1 hHigh1 with
    ⟨heap2, ignoredLength, hadd, hlayout2, _⟩
  have hOut2 : heap2.ValidU64Slice out fullLength :=
    (hlayout2 out fullLength).mp hOut1
  rcases normaliseU64_ok heap2 out fullLength hOut2 with
    ⟨length, hnorm, hlength⟩
  refine ⟨HgcdLiftHighResult.mk heap2 length, ?_, ?_, ?_⟩
  · simp [hgcdRecursiveLiftHigh, required, fullLength, oldHighLength,
      hzeroRun, hadd, hnorm]
  · exact fun ptr count =>
      (hlayout1 ptr count).trans (hlayout2 ptr count)
  · simpa [fullLength, required] using hlength

/-- Full semantic refinement of the generated shifted-high block.  The
temporary padded whole is treated only as a slice representation; the result
becomes normalized solely through the source's final `_poly_normalise`. -/
theorem hgcdRecursiveLiftHigh_refines (this : DenseUPolyZp)
    (out high : RawPtr UInt64) (lowLength m highLength : Nat)
    (heap : RawHeap)
    (lowPoly highPoly : Polynomial (ZMod this._p.toNat))
    (hp : this._p ≠ 0)
    (hwork : HgcdLiftHighWorkspace heap out high lowLength m highLength)
    (hLow : RawDensePolyRep this heap out lowLength lowPoly)
    (hHigh : RawDensePolyRep this heap high highLength highPoly) :
    ∃ result,
      hgcdRecursiveLiftHigh this out high lowLength m highLength heap =
        .ok result ∧
      RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap out result.length
        (lowPoly + Polynomial.X ^ m * highPoly) := by
  let fullLength := max (m + highLength) lowLength
  let oldHighLength := if m ≤ lowLength then lowLength - m else 0
  rcases hgcdLiftHigh_zero_refines this out lowLength m highLength heap
      lowPoly hp hwork.outValid hLow with
    ⟨heap1, hzero, hlayout1, hLow1, hFullSlice1, hFullCanonical1⟩
  have hHigh1 : RawDensePolyRep this heap1 high highLength highPoly := by
    by_cases hpad : lowLength < m + highLength
    · have hrun : Generated.StrictMul.mulZeroPadLoop out lowLength
          (m + highLength - lowLength) 0 heap = .ok heap1 := by
        simpa [hpad] using hzero
      have hvalid : heap.ValidU64Slice out
          (lowLength + (m + highLength - lowLength)) := by
        have hmax : max (m + highLength) lowLength = m + highLength :=
          Nat.max_eq_left (Nat.le_of_lt hpad)
        simpa [hmax, Nat.add_sub_of_le (Nat.le_of_lt hpad)] using
          hwork.outValid
      have hsame := mulZeroPadLoop_preserves_prefix out high lowLength
        (m + highLength - lowLength) 0 highLength heap heap1 hvalid
        (by
          simpa [Nat.add_sub_of_le (Nat.le_of_lt hpad)] using
            u64SlicesDisjoint_mono hwork.zeroHighDisjoint
              (by omega) (Nat.le_refl _)) hrun
      exact rawDensePolyRep_of_same_prefix this heap heap1 high highLength
        highPoly hlayout1 hsame hHigh
    · have heq : heap1 = heap := by
        have : (Except.ok heap : RawExec RawHeap) = Except.ok heap1 := by
          simpa [hpad] using hzero
        exact (Except.ok.inj this).symm
      subst heap1
      exact hHigh
  have hparts : ∃ lowPart existingHigh : Polynomial (ZMod this._p.toNat),
      SlicePolyRep heap1 out m this._p.toNat lowPart ∧
      RawDensePolyRep this heap1 (out.add m) oldHighLength existingHigh ∧
      lowPoly = lowPart + Polynomial.X ^ m * existingHigh := by
    by_cases hm : m ≤ lowLength
    · simpa [oldHighLength, hm] using
        rawDensePolyRep_split_suffix this heap1 out lowLength m lowPoly hm hLow1
    · have hm' : lowLength < m := by omega
      have hmFull : m ≤ fullLength := by
        dsimp [fullLength]
        omega
      have hLowPrefix : SlicePolyRep heap1 out m this._p.toNat lowPoly :=
        slicePolyRep_prefix_of_coeff_zero heap1 out fullLength m
          this._p.toNat lowPoly
          (by simpa [fullLength] using
            (hlayout1 out fullLength).mp (by simpa [fullLength] using
              hwork.outValid)) hmFull
          (by simpa [fullLength] using hFullSlice1) (by
            intro degree hdegree
            exact slicePolyRep_coeff_zero_of_length_le heap1 out lowLength
              this._p.toNat lowPoly hLow1.2.2.1 degree (by omega))
      have hzeroValid : heap1.ValidU64Slice (out.add m) 0 := by
        exact heap1.validU64Slice_add out fullLength m 0
          ((hlayout1 out fullLength).mp (by simpa [fullLength] using
            hwork.outValid)) (by omega)
      refine ⟨lowPoly, 0, hLowPrefix, ?_, by simp⟩
      simpa [oldHighLength, hm] using
        rawDensePolyRep_zero_length this heap1 (out.add m) hzeroValid
  rcases hparts with ⟨lowPart, existingHigh, hLowPart, hExisting, hsplit⟩
  have hOut1 : heap1.ValidU64Slice out fullLength :=
    (hlayout1 out fullLength).mp (by simpa [fullLength] using hwork.outValid)
  have hAddOut1 : heap1.ValidU64Slice (out.add m)
      (max oldHighLength highLength) := by
    apply heap1.validU64Slice_add out fullLength m
      (max oldHighLength highLength) hOut1
    dsimp [fullLength, oldHighLength]
    split <;> omega
  rcases polyAdd_ok this (out.add m) (out.add m) oldHighLength high
      highLength heap1 hAddOut1 hExisting.1 hHigh1.1 with
    ⟨heap2, addLength, hadd, hlayout2, _⟩
  have hAdded := polyAdd_refines this (out.add m) (out.add m)
    oldHighLength high highLength heap1 heap2 addLength existingHigh highPoly
    hp hAddOut1 hExisting hHigh1 (Or.inl rfl) hwork.addAliasHigh hadd
  have hAddNorm := polyAdd_result_normalise this (out.add m) (out.add m)
    oldHighLength high highLength heap1 heap2 addLength hadd
  have hAddOut2 : heap2.ValidU64Slice (out.add m)
      (max oldHighLength highLength) :=
    (hlayout2 (out.add m) _).mp hAddOut1
  have hAddedFull := rawDensePolyRep_extend_to_normalise_input this heap2
    (out.add m) (max oldHighLength highLength) addLength
    (existingHigh + highPoly) hp hAddOut2 hAdded hAddNorm
  have hLowSame : SameU64Prefix heap1 heap2 out m :=
    polyAdd_preserves_prefix_disjoint this (out.add m) (out.add m)
      oldHighLength high out highLength m heap1 heap2 addLength hAddOut1
      hExisting.1 hHigh1.1 (by
        intro writeIndex _ readIndex _
        right
        have hwidth : RawLimbWidth.width UInt64 = 1 := rfl
        simp [RawPtr.add, hwidth]
        omega) hadd
  have hLowPart2 : SlicePolyRep heap2 out m this._p.toNat lowPart :=
    slicePolyRep_of_same_prefix heap1 heap2 out m this._p.toNat lowPart
      (heap1.validU64Slice_mono out fullLength m hOut1 (by omega))
      ((hlayout2 out m).mp
        (heap1.validU64Slice_mono out fullLength m hOut1 (by omega)))
      hLowSame hLowPart
  have hCanonicalLow1 := canonicalU64Prefix_mono heap1 out fullLength m
    this._p (by omega) (by simpa [fullLength] using hFullCanonical1)
  have hCanonicalLow2 := canonicalU64Prefix_of_same_prefix heap1 heap2 out
    m this._p (heap1.validU64Slice_mono out fullLength m hOut1 (by omega))
    hLowSame hCanonicalLow1
  have hlengthEq : m + max oldHighLength highLength = fullLength := by
    dsimp [oldHighLength, fullLength]
    split <;> omega
  have hOut2 : heap2.ValidU64Slice out fullLength :=
    (hlayout2 out fullLength).mp hOut1
  have hWholeSlice : SlicePolyRep heap2 out fullLength this._p.toNat
      (lowPoly + Polynomial.X ^ m * highPoly) := by
    have hjoin := slicePolyRep_join heap2 out m
      (max oldHighLength highLength) this._p.toNat lowPart
      (existingHigh + highPoly) (by simpa [hlengthEq] using hOut2)
      hLowPart2 hAddedFull.1
    have hpoly : lowPart + Polynomial.X ^ m *
          (existingHigh + highPoly) =
        lowPoly + Polynomial.X ^ m * highPoly := by
      rw [hsplit]
      ring
    simpa [hlengthEq, hpoly] using hjoin
  have hWholeCanonical : CanonicalU64Prefix heap2 out fullLength this._p := by
    have hjoin := canonicalU64Prefix_join heap2 out m
      (max oldHighLength highLength) this._p
      (by simpa [hlengthEq] using hOut2) hCanonicalLow2 hAddedFull.2
    simpa [hlengthEq] using hjoin
  rcases normaliseU64_ok heap2 out fullLength hOut2 with
    ⟨length, hnorm, hlength⟩
  let result := HgcdLiftHighResult.mk heap2 length
  refine ⟨result, ?_, ?_, ?_⟩
  · simp [result, hgcdRecursiveLiftHigh, fullLength, oldHighLength,
      hzero, hadd, hnorm]
  · exact fun ptr count =>
      (hlayout1 ptr count).trans (hlayout2 ptr count)
  · refine ⟨heap2.validU64Slice_mono out fullLength length hOut2 hlength,
      canonicalU64Prefix_mono heap2 out fullLength length this._p hlength
        hWholeCanonical,
      slicePolyRep_of_normaliseU64 heap2 out fullLength this._p.toNat length
        (lowPoly + Polynomial.X ^ m * highPoly) hOut2 hWholeSlice hnorm,
      normaliseU64_result_fixed heap2 out fullLength length hOut2 hnorm⟩

/-- Physical capacities needed by the optional four-entry early-return
matrix copy.  Source lengths come from `R`, exactly as in C++. -/
structure HgcdEarlyMatrixWorkspace (heap : RawHeap) (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid) : Prop where
  targetValid : ∀ i : Fin 4, heap.ValidU64Slice
    (hgcdMatPtrRaw M hM i) (hgcdMatLenRaw R hR i)
  sourceValid : ∀ i : Fin 4, heap.ValidU64Slice
    (hgcdMatPtrRaw R hR i) (hgcdMatLenRaw R hR i)

/-- Non-aliasing obligations that make the four sequential C++ matrix copies
semantic, rather than merely executable.  Later destinations cannot overwrite
an earlier result or any still-live source entry. -/
structure HgcdEarlyMatrixRefineWorkspace (heap : RawHeap) (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid) extends
    HgcdEarlyMatrixWorkspace heap M R hM hR where
  targetSourceDisjoint : ∀ i j : Fin 4, U64SlicesDisjoint
    (hgcdMatPtrRaw M hM i) (hgcdMatLenRaw R hR i)
    (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j)
  targetPairwiseDisjoint : ∀ i j : Fin 4, i ≠ j → U64SlicesDisjoint
    (hgcdMatPtrRaw M hM i) (hgcdMatLenRaw R hR i)
    (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j)

theorem hgcdEarlyMatrixRefineWorkspace_of_sameLayout
    (before after : RawHeap) (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid)
    (hlayout : RawHeap.SameLayout before after)
    (hwork : HgcdEarlyMatrixRefineWorkspace before M R hM hR) :
    HgcdEarlyMatrixRefineWorkspace after M R hM hR := by
  exact {
    targetValid := fun i => (hlayout _ _).mp (hwork.targetValid i)
    sourceValid := fun i => (hlayout _ _).mp (hwork.sourceValid i)
    targetSourceDisjoint := hwork.targetSourceDisjoint
    targetPairwiseDisjoint := hwork.targetPairwiseDisjoint }

theorem hgcdEarlyMatrixLoop_terminates (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid) (i : Nat) (heap : RawHeap)
    (hwork : HgcdEarlyMatrixWorkspace heap M R hM hR) :
    ∃ result, hgcdEarlyMatrixLoop M R hM hR i heap = .ok result ∧
      RawHeap.SameLayout heap result.heap ∧ result.matrix.Valid := by
  rw [hgcdEarlyMatrixLoop]
  split
  next hi =>
    let index : Fin 4 := ⟨i, hi⟩
    rcases copyU64_ok heap (hgcdMatPtrRaw M hM index)
        (hgcdMatPtrRaw R hR index) (hgcdMatLenRaw R hR index)
        (hwork.targetValid index) (hwork.sourceValid index) with
      ⟨heap1, hcopy, hlayout1⟩
    simp only [hcopy]
    let nextLen := M.len.set i (hgcdMatLenRaw R hR index)
      (by rw [hM.2]; exact hi)
    let next : HgcdMat := { M with len := nextLen }
    have hNext : next.Valid := by
      exact ⟨hM.1, by simp [next, nextLen, hM.2]⟩
    have hwork1 : HgcdEarlyMatrixWorkspace heap1 next R hNext hR := by
      constructor
      · intro j
        apply (hlayout1 (hgcdMatPtrRaw next hNext j)
          (hgcdMatLenRaw R hR j)).mp
        simpa [hgcdMatPtrRaw, next] using hwork.targetValid j
      · intro j
        exact (hlayout1 (hgcdMatPtrRaw R hR j)
          (hgcdMatLenRaw R hR j)).mp (hwork.sourceValid j)
    rcases hgcdEarlyMatrixLoop_terminates next R hNext hR (i + 1)
        heap1 hwork1 with ⟨result, hrun, hlayout2, hvalid⟩
    refine ⟨result, ?_, fun ptr count =>
      (hlayout1 ptr count).trans (hlayout2 ptr count), hvalid⟩
    simpa [index, next, nextLen, hcopy] using hrun
  next hi =>
    exact ⟨HgcdEarlyMatrixResult.mk heap M, rfl, fun _ _ => Iff.rfl, hM⟩
termination_by 4 - i
decreasing_by omega

/-- Content semantics of the actual four-copy early-return loop.  The
induction follows the generated loop and records exactly which destination
slots have already received the corresponding normalized source polynomial. -/
theorem hgcdEarlyMatrixLoop_copies (this : DenseUPolyZp)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (i : Nat) (heap : RawHeap) (result : HgcdEarlyMatrixResult)
    (hwork : HgcdEarlyMatrixRefineWorkspace heap M R hM hR)
    (hsource : ∀ j : Fin 4, RawDensePolyRep this heap
      (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j))
    (hdone : ∀ j : Fin 4, j.val < i → RawDensePolyRep this heap
      (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j) (entries j))
    (hrun : hgcdEarlyMatrixLoop M R hM hR i heap = .ok result) :
    (∀ j : Fin 4, RawDensePolyRep this result.heap
      (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j)) ∧
    ∀ j : Fin 4, RawDensePolyRep this result.heap
      (hgcdMatPtrRaw result.matrix
        (hgcdEarlyMatrixLoop_result_valid M R hM hR i heap result hrun) j)
      (hgcdMatLenRaw R hR j) (entries j) := by
  rw [hgcdEarlyMatrixLoop] at hrun
  split at hrun
  next hi =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hi⟩
      let dst := hgcdMatPtrRaw M hM index
      let src := hgcdMatPtrRaw R hR index
      let length := hgcdMatLenRaw R hR index
      rcases copyU64_refines_rawDense this heap dst src length (entries index)
          (hwork.targetValid index) (hwork.targetSourceDisjoint index index)
          (hsource index) with ⟨semanticHeap, hsemantic, hlayout, hcopied⟩
      have heq : semanticHeap = heap1 :=
        Except.ok.inj (hsemantic.symm.trans hcopy)
      subst semanticHeap
      have hsource1 : ∀ j : Fin 4, RawDensePolyRep this heap1
          (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j) := by
        intro j
        exact (copyU64_preserves_rawDenseRep this heap heap1 dst src length
          (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j)
          (hwork.targetValid index) (hwork.sourceValid index)
          (hwork.targetSourceDisjoint index j) hcopy (hsource j)).2
      let nextLen := M.len.set i (hgcdMatLenRaw R hR index)
        (by rw [hM.2]; exact hi)
      let next : HgcdMat := { M with len := nextLen }
      have hNext : next.Valid := by
        exact ⟨hM.1, by simp [next, nextLen, hM.2]⟩
      have hwork1 : HgcdEarlyMatrixRefineWorkspace heap1 next R hNext hR := by
        have transported := hgcdEarlyMatrixRefineWorkspace_of_sameLayout
          heap heap1 M R hM hR hlayout hwork
        exact {
          targetValid := fun j => by
            simpa [next, hgcdMatPtrRaw] using transported.targetValid j
          sourceValid := transported.sourceValid
          targetSourceDisjoint := fun j k => by
            simpa [next, hgcdMatPtrRaw] using
              transported.targetSourceDisjoint j k
          targetPairwiseDisjoint := fun j k hjk => by
            simpa [next, hgcdMatPtrRaw] using
              transported.targetPairwiseDisjoint j k hjk }
      have hdone1 : ∀ j : Fin 4, j.val < i + 1 →
          RawDensePolyRep this heap1 (hgcdMatPtrRaw next hNext j)
            (hgcdMatLenRaw R hR j) (entries j) := by
        intro j hj
        by_cases hji : j = index
        · subst j
          simpa [dst, next, hgcdMatPtrRaw] using hcopied
        · have hjlt : j.val < i := by
            have : j.val ≠ i := by
              intro hval
              exact hji (Fin.ext hval)
            omega
          have hpreserved := (copyU64_preserves_rawDenseRep this heap heap1
            dst src length (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j)
            (entries j) (hwork.targetValid index) (hwork.sourceValid index)
            (hwork.targetPairwiseDisjoint index j (Ne.symm hji)) hcopy
            (hdone j hjlt)).2
          simpa [next, hgcdMatPtrRaw] using hpreserved
      exact hgcdEarlyMatrixLoop_copies this next R hNext hR entries (i + 1)
        heap1 result hwork1 hsource1 hdone1 hrun
  next hi =>
    have heq : result = HgcdEarlyMatrixResult.mk heap M :=
      (Except.ok.inj hrun).symm
    subst result
    exact ⟨hsource, fun j => hdone j (by omega)⟩
termination_by 4 - i
decreasing_by omega

/-- Entry-point refinement of the optional matrix branch: the returned
descriptor has `R`'s four lengths and its original destination pointers now
represent exactly `R`'s four L2 polynomial entries. -/
theorem hgcdEarlyMatrixLoop_zero_refines (this : DenseUPolyZp)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap) (result : HgcdEarlyMatrixResult)
    (hwork : HgcdEarlyMatrixRefineWorkspace heap M R hM hR)
    (hsource : HgcdMatRawDenseRep this heap R entries hR)
    (hrun : hgcdEarlyMatrixLoop M R hM hR 0 heap = .ok result) :
    ∃ hResult : result.matrix.Valid,
      result.matrix.len = R.len ∧
      HgcdMatRawDenseRep this result.heap result.matrix entries hResult := by
  have hvalid := hgcdEarlyMatrixLoop_result_valid M R hM hR 0 heap result hrun
  have hlen := hgcdEarlyMatrixLoop_lengths M R hM hR 0 heap result
    (by intro j hj; omega) hrun
  have hcopies := hgcdEarlyMatrixLoop_copies this M R hM hR entries 0 heap
    result hwork (by
      intro j
      simpa [HgcdMatRawDenseRep, hgcdMatPtr, hgcdMatLen,
        hgcdMatPtrRaw, hgcdMatLenRaw] using hsource j)
    (by intro j hj; omega) hrun
  refine ⟨hvalid, hlen, ?_⟩
  intro j
  have hLength : hgcdMatLen result.matrix hvalid j =
      hgcdMatLenRaw R hR j := by
    exact array_getElem_eq_of_eq result.matrix.len R.len hlen j.val
      (by rw [hvalid.2]; exact j.isLt)
      (by rw [hR.2]; exact j.isLt)
  rw [hLength]
  simpa [hgcdMatPtr, hgcdMatPtrRaw] using hcopies.2 j

/-- A normalized polynomial outside all four matrix destinations is framed by
the complete generated early-copy loop. -/
theorem hgcdEarlyMatrixLoop_preserves_rawDenseRep (this : DenseUPolyZp)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (i : Nat) (heap : RawHeap) (result : HgcdEarlyMatrixResult)
    (ptr : RawPtr UInt64) (length : Nat)
    (poly : Polynomial (ZMod this._p.toNat))
    (hwork : HgcdEarlyMatrixWorkspace heap M R hM hR)
    (hframe : ∀ j : Fin 4, U64SlicesDisjoint
      (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j) ptr length)
    (hrep : RawDensePolyRep this heap ptr length poly)
    (hrun : hgcdEarlyMatrixLoop M R hM hR i heap = .ok result) :
    RawDensePolyRep this result.heap ptr length poly := by
  rw [hgcdEarlyMatrixLoop] at hrun
  split at hrun
  next hi =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hi⟩
      have hcopyFrame := copyU64_preserves_rawDenseRep this heap heap1
        (hgcdMatPtrRaw M hM index) (hgcdMatPtrRaw R hR index)
        (hgcdMatLenRaw R hR index) ptr length poly
        (hwork.targetValid index) (hwork.sourceValid index) (hframe index)
        hcopy hrep
      have hpreserved := hcopyFrame.2
      let nextLen := M.len.set i (hgcdMatLenRaw R hR index)
        (by rw [hM.2]; exact hi)
      let next : HgcdMat := { M with len := nextLen }
      have hNext : next.Valid := by
        exact ⟨hM.1, by simp [next, nextLen, hM.2]⟩
      have hwork1 : HgcdEarlyMatrixWorkspace heap1 next R hNext hR := by
        exact {
          targetValid := fun j => by
            apply (hcopyFrame.1 _ _).mp
            simpa [next, hgcdMatPtrRaw] using hwork.targetValid j
          sourceValid := fun j => by
            apply (hcopyFrame.1 _ _).mp
            exact hwork.sourceValid j }
      have hframe1 : ∀ j : Fin 4, U64SlicesDisjoint
          (hgcdMatPtrRaw next hNext j) (hgcdMatLenRaw R hR j) ptr length := by
        intro j
        simpa [next, hgcdMatPtrRaw] using hframe j
      exact hgcdEarlyMatrixLoop_preserves_rawDenseRep this next R hNext hR
        (i + 1) heap1 result ptr length poly hwork1 hframe1 hpreserved hrun
  next hi =>
    have heq : result = HgcdEarlyMatrixResult.mk heap M :=
      (Except.ok.inj hrun).symm
    subst result
    exact hrep
termination_by 4 - i
decreasing_by omega

structure HgcdEarlyReturnWorkspace (heap : RawHeap)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (A B a2 b2 : RawPtr UInt64) (lenA2 lenB2 : Nat) : Prop where
  AValid : heap.ValidU64Slice A lenA2
  BValid : heap.ValidU64Slice B lenB2
  a2Valid : heap.ValidU64Slice a2 lenA2
  b2Valid : heap.ValidU64Slice b2 lenB2
  matrix : HgcdEarlyMatrixWorkspace heap M R hM hR

/-- Full separation contract for the real `_hgcd_recursive` early return.
It mentions only physical slices; expected L2 output values remain theorem
inputs and are obtained from the actual memcpy executions. -/
structure HgcdEarlyReturnRefineWorkspace (heap : RawHeap)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (A B a2 b2 : RawPtr UInt64) (lenA2 lenB2 : Nat) : Prop where
  AValid : heap.ValidU64Slice A lenA2
  BValid : heap.ValidU64Slice B lenB2
  matrix : HgcdEarlyMatrixRefineWorkspace heap M R hM hR
  Aa2Disjoint : U64SlicesDisjoint A lenA2 a2 lenA2
  Ab2Disjoint : U64SlicesDisjoint A lenA2 b2 lenB2
  Bb2Disjoint : U64SlicesDisjoint B lenB2 b2 lenB2
  BADisjoint : U64SlicesDisjoint B lenB2 A lenA2
  ARDisjoint : ∀ j : Fin 4, U64SlicesDisjoint A lenA2
    (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j)
  BRDisjoint : ∀ j : Fin 4, U64SlicesDisjoint B lenB2
    (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j)
  matrixADisjoint : ∀ j : Fin 4, U64SlicesDisjoint
    (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j) A lenA2
  matrixBDisjoint : ∀ j : Fin 4, U64SlicesDisjoint
    (hgcdMatPtrRaw M hM j) (hgcdMatLenRaw R hR j) B lenB2

/-- Total execution bridge for the generated early-return branch, including
the optional matrix copy. -/
theorem hgcdRecursiveEarlyReturn_terminates (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid) (computeM : Bool)
    (A B a2 b2 : RawPtr UInt64) (lenA2 lenB2 : Nat) (sgn : Int)
    (heap : RawHeap)
    (hwork : HgcdEarlyReturnWorkspace heap M R hM hR A B a2 b2
      lenA2 lenB2) :
    ∃ result,
      hgcdRecursiveEarlyReturn M R hM hR computeM A B a2 b2
        lenA2 lenB2 sgn heap = .ok result ∧
      RawHeap.SameLayout heap result.heap ∧
      result.matrix.Valid ∧ result.lenA = lenA2 ∧
      result.lenB = lenB2 ∧ result.sgn = sgn := by
  rcases copyU64_ok heap A a2 lenA2 hwork.AValid hwork.a2Valid with
    ⟨heap1, hcopyA, hlayout1⟩
  have hB1 := (hlayout1 B lenB2).mp hwork.BValid
  have hb21 := (hlayout1 b2 lenB2).mp hwork.b2Valid
  rcases copyU64_ok heap1 B b2 lenB2 hB1 hb21 with
    ⟨heap2, hcopyB, hlayout2⟩
  by_cases hcompute : computeM = true
  · have hmatrix2 : HgcdEarlyMatrixWorkspace heap2 M R hM hR := by
      constructor <;> intro j
      · exact (hlayout2 _ _).mp ((hlayout1 _ _).mp
          (hwork.matrix.targetValid j))
      · exact (hlayout2 _ _).mp ((hlayout1 _ _).mp
          (hwork.matrix.sourceValid j))
    rcases hgcdEarlyMatrixLoop_terminates M R hM hR 0 heap2 hmatrix2 with
      ⟨matrixResult, hmatrix, hlayout3, hvalid⟩
    refine ⟨HgcdRecursiveEarlyResult.mk matrixResult.heap
        matrixResult.matrix lenA2 lenB2 sgn, ?_, ?_, hvalid,
      rfl, rfl, rfl⟩
    · simp [hgcdRecursiveEarlyReturn, hcopyA, hcopyB, hcompute, hmatrix]
    · exact fun ptr count =>
        (hlayout1 ptr count).trans
          ((hlayout2 ptr count).trans (hlayout3 ptr count))
  · have hfalse : computeM = false := by cases computeM <;> simp_all
    refine ⟨HgcdRecursiveEarlyResult.mk heap2 M lenA2 lenB2 sgn,
      ?_, ?_, hM, rfl, rfl, rfl⟩
    · simp [hgcdRecursiveEarlyReturn, hcopyA, hcopyB, hfalse]
    · exact fun ptr count =>
        (hlayout1 ptr count).trans (hlayout2 ptr count)

/-- End-to-end semantic refinement of the generated early-return block.  It
executes both output memcpys and, when requested, the four matrix memcpys;
the returned A/B polynomials and matrix are therefore consequences of those
same RawHeap executions. -/
theorem hgcdRecursiveEarlyReturn_refines (this : DenseUPolyZp)
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid) (computeM : Bool)
    (A B a2 b2 : RawPtr UInt64) (lenA2 lenB2 : Nat) (sgn : Int)
    (left right : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (heap : RawHeap)
    (hwork : HgcdEarlyReturnRefineWorkspace heap M R hM hR
      A B a2 b2 lenA2 lenB2)
    (hLeft : RawDensePolyRep this heap a2 lenA2 left)
    (hRight : RawDensePolyRep this heap b2 lenB2 right)
    (hMatrix : HgcdMatRawDenseRep this heap R entries hR) :
    ∃ result,
      hgcdRecursiveEarlyReturn M R hM hR computeM A B a2 b2
        lenA2 lenB2 sgn heap = .ok result ∧
      RawHeap.SameLayout heap result.heap ∧
      result.lenA = lenA2 ∧ result.lenB = lenB2 ∧ result.sgn = sgn ∧
      RawDensePolyRep this result.heap A lenA2 left ∧
      RawDensePolyRep this result.heap B lenB2 right ∧
      ∃ hResult : result.matrix.Valid,
        (computeM = true →
          result.matrix.len = R.len ∧
          HgcdMatRawDenseRep this result.heap result.matrix entries hResult) ∧
        (computeM = false → result.matrix = M) := by
  rcases copyU64_refines_rawDense this heap A a2 lenA2 left hwork.AValid
      hwork.Aa2Disjoint hLeft with ⟨heap1, hcopyA, hlayout1, hA1⟩
  have hRight1 := (copyU64_preserves_rawDenseRep this heap heap1 A a2 lenA2
    b2 lenB2 right hwork.AValid hLeft.1 hwork.Ab2Disjoint hcopyA hRight).2
  have hMatrix1 : HgcdMatRawDenseRep this heap1 R entries hR := by
    intro j
    exact (copyU64_preserves_rawDenseRep this heap heap1 A a2 lenA2
      (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j)
      hwork.AValid hLeft.1 (hwork.ARDisjoint j) hcopyA (hMatrix j)).2
  have hB1 := (hlayout1 B lenB2).mp hwork.BValid
  rcases copyU64_refines_rawDense this heap1 B b2 lenB2 right hB1
      hwork.Bb2Disjoint hRight1 with ⟨heap2, hcopyB, hlayout2, hB2⟩
  have hA2 := (copyU64_preserves_rawDenseRep this heap1 heap2 B b2 lenB2
    A lenA2 left hB1 hRight1.1 hwork.BADisjoint hcopyB hA1).2
  have hMatrix2 : HgcdMatRawDenseRep this heap2 R entries hR := by
    intro j
    exact (copyU64_preserves_rawDenseRep this heap1 heap2 B b2 lenB2
      (hgcdMatPtrRaw R hR j) (hgcdMatLenRaw R hR j) (entries j)
      hB1 hRight1.1 (hwork.BRDisjoint j) hcopyB (hMatrix1 j)).2
  have hMatrixWork1 := hgcdEarlyMatrixRefineWorkspace_of_sameLayout heap
    heap1 M R hM hR hlayout1 hwork.matrix
  have hMatrixWork2 := hgcdEarlyMatrixRefineWorkspace_of_sameLayout heap1
    heap2 M R hM hR hlayout2 hMatrixWork1
  by_cases hcompute : computeM = true
  · rcases hgcdEarlyMatrixLoop_terminates M R hM hR 0 heap2
        hMatrixWork2.toHgcdEarlyMatrixWorkspace with
      ⟨matrixResult, hmatrixRun, hlayout3, hmatrixValid⟩
    rcases hgcdEarlyMatrixLoop_zero_refines this M R hM hR entries heap2
        matrixResult hMatrixWork2 hMatrix2 hmatrixRun with
      ⟨hResult, hResultLen, hResultRep⟩
    have hAResult := hgcdEarlyMatrixLoop_preserves_rawDenseRep this M R hM hR
      0 heap2 matrixResult A lenA2 left
      hMatrixWork2.toHgcdEarlyMatrixWorkspace hwork.matrixADisjoint hA2
      hmatrixRun
    have hBResult := hgcdEarlyMatrixLoop_preserves_rawDenseRep this M R hM hR
      0 heap2 matrixResult B lenB2 right
      hMatrixWork2.toHgcdEarlyMatrixWorkspace hwork.matrixBDisjoint hB2
      hmatrixRun
    let result : HgcdRecursiveEarlyResult := {
      heap := matrixResult.heap
      matrix := matrixResult.matrix
      lenA := lenA2
      lenB := lenB2
      sgn := sgn }
    refine ⟨result, ?_, ?_, rfl, rfl, rfl, hAResult, hBResult,
      hResult, ?_, ?_⟩
    · simp [hgcdRecursiveEarlyReturn, hcopyA, hcopyB, hcompute, hmatrixRun,
        result]
    · exact fun ptr count => (hlayout1 ptr count).trans
        ((hlayout2 ptr count).trans (hlayout3 ptr count))
    · intro _
      exact ⟨hResultLen, hResultRep⟩
    · intro hfalse
      simp_all
  · have hfalse : computeM = false := by cases computeM <;> simp_all
    let result : HgcdRecursiveEarlyResult := {
      heap := heap2
      matrix := M
      lenA := lenA2
      lenB := lenB2
      sgn := sgn }
    refine ⟨result, ?_, ?_, rfl, rfl, rfl, hA2, hB2, hM, ?_, ?_⟩
    · simp [hgcdRecursiveEarlyReturn, hcopyA, hcopyB, hfalse, result]
    · exact fun ptr count =>
        (hlayout1 ptr count).trans (hlayout2 ptr count)
    · intro htrue
      simp_all
    · intro _
      rfl

/-- Strict semantic bridge for the middle source `_poly_divrem`.  Its result
record also exposes the exact `k`, `c0`, and `d0` passed to the second HGCD
call, while quotient and remainder come from the same raw execution. -/
theorem hgcdRecursiveMiddle_refines (this : DenseUPolyZp)
    (q d a2 b2 : RawPtr UInt64) (lenA2 lenB2 m lenA : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (dividend divisor : Polynomial (ZMod this._p.toNat))
    (hlenB : 0 < lenB2)
    (hA : RawDensePolyRep this heap a2 lenA2 dividend)
    (hB : RawDensePolyRep this heap b2 lenB2 divisor)
    (hQ : heap.ValidU64Slice q (lenA2 - (lenB2 - 1)))
    (hD : heap.ValidU64Slice d (Nat.min lenA2 (lenB2 - 1)))
    (hW3 : heap.ValidWord3Slice W3 lenA2)
    (hqCapacity : lenA2 - (lenB2 - 1) < limbBase)
    (hDA : d.region ≠ a2.region)
    (hWA : W3.region ≠ a2.region) (hWB : W3.region ≠ b2.region)
    (hQB : q.region ≠ b2.region) (hQW : q.region ≠ W3.region)
    (hDW : d.region ≠ W3.region) (hDQ : d.region ≠ q.region)
    (hDB : d.region ≠ b2.region)
    (hcfg : DensePreinvConfigured this)
    (hprime : Nat.Prime this._p.toNat)
    (hlenA : 0 < lenA) (hlenB2Bound : lenB2 ≤ lenA)
    (hm : 0 < m) (hnonearly : m + 1 ≤ lenB2) :
    ∃ result quotient remainder,
      hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3 heap =
        .ok result ∧
      SlicePolyRep result.heap q result.lenQ this._p.toNat quotient ∧
      CanonicalU64Prefix result.heap q result.lenQ this._p ∧
      result.heap.normaliseU64 q result.lenQ = .ok result.lenQ ∧
      SlicePolyRep result.heap d result.lenD this._p.toNat remainder ∧
      CanonicalU64Prefix result.heap d result.lenD this._p ∧
      result.heap.normaliseU64 d result.lenD = .ok result.lenD ∧
      RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap q result.lenQ quotient ∧
      RawDensePolyRep this result.heap d result.lenD remainder ∧
      RawDensePolyRep this result.heap b2 lenB2 divisor ∧
      dividend = quotient * divisor + remainder ∧
      (remainder = 0 ∨ remainder.natDegree < divisor.natDegree) ∧
      result.lenD < lenB2 ∧ result.lenC0 < lenA ∧
      0 < result.lenC0 ∧ result.lenD0 < result.lenC0 ∧
      result.k = 2 * m - lenB2 + 1 ∧
      result.c0 = b2.add result.k ∧
      result.lenC0 = (if lenB2 ≥ result.k then lenB2 - result.k else 0) ∧
      result.d0 = d.add result.k ∧
      result.lenD0 = (if result.lenD ≥ result.k then
        result.lenD - result.k else 0) := by
  rcases polyDivrem_refines this q d a2 b2 lenA2 lenB2 W3 heap
      dividend divisor hlenB hA.1 hB.1 hQ hD hW3 hA.2.1 hB.2.1
      hA.2.2.1 hB.2.2.1 hA.2.2.2 hB.2.2.2 hqCapacity hDA hWA hWB
      hQB hQW hDW hDQ hDB hcfg hprime with
    ⟨heap1, lenQ, lenD, quotient, remainder, hdiv, hQRep, hQCanonical,
      hQNorm, hDRep, hDCanonical, hDNorm, hlayout, hsameB, hidentity,
      hdegree, hlenQ, hlenD, hlt⟩
  let result : HgcdRecursiveMiddleResult := {
    heap := heap1
    lenQ := lenQ
    lenD := lenD
    k := 2 * m - lenB2 + 1
    c0 := b2.add (2 * m - lenB2 + 1)
    lenC0 := if lenB2 ≥ 2 * m - lenB2 + 1 then
      lenB2 - (2 * m - lenB2 + 1) else 0
    d0 := d.add (2 * m - lenB2 + 1)
    lenD0 := if lenD ≥ 2 * m - lenB2 + 1 then
      lenD - (2 * m - lenB2 + 1) else 0 }
  have hmiddle : hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3
      heap = .ok result := by
    simp [hgcdRecursiveMiddle, hdiv, result]
  have hQValid : heap1.ValidU64Slice q lenQ :=
    heap1.validU64Slice_mono q (lenA2 - (lenB2 - 1)) lenQ
      ((hlayout q _).mp hQ) hlenQ
  have hDValid : heap1.ValidU64Slice d lenD :=
    heap1.validU64Slice_mono d (Nat.min lenA2 (lenB2 - 1)) lenD
      ((hlayout d _).mp hD) hlenD
  have hQDense : RawDensePolyRep this heap1 q lenQ quotient :=
    ⟨hQValid, hQCanonical, hQRep, hQNorm⟩
  have hDDense : RawDensePolyRep this heap1 d lenD remainder :=
    ⟨hDValid, hDCanonical, hDRep, hDNorm⟩
  have hBDense := rawDensePolyRep_of_same_prefix this heap heap1 b2 lenB2
    divisor hlayout hsameB hB
  rcases hgcdRecursiveMiddle_second_call_bounds this q d a2 b2 lenA2 lenB2
      m lenA W3 heap result hlenA hm hlenB2Bound hnonearly hlt hmiddle with
    ⟨hcpos, hdorder, hdecrease⟩
  refine ⟨result, quotient, remainder, hmiddle, hQRep, hQCanonical, hQNorm,
    hDRep, hDCanonical, hDNorm, hlayout, hQDense, hDDense, hBDense,
    hidentity, hdegree, hlt,
    hdecrease, hcpos, hdorder, rfl, rfl, rfl, rfl, rfl⟩

/-- The nonempty second-call input is obtained by slicing the actual divisor
and remainder buffers at the generated `k`.  No heap write or L2 evaluation
is used to create `c0` or `d0`; a zero remainder suffix uses only its physical
zero-length validity. -/
theorem hgcdRecursiveMiddle_suffix_reps (this : DenseUPolyZp)
    (q d a2 b2 : RawPtr UInt64) (lenA2 lenB2 m : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (result : HgcdRecursiveMiddleResult)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (hDivisor : RawDensePolyRep this result.heap b2 lenB2 divisor)
    (hRemainder : RawDensePolyRep this result.heap d result.lenD remainder)
    (hD0Valid : result.heap.ValidU64Slice result.d0 result.lenD0)
    (hC0Pos : 0 < result.lenC0)
    (hrun : hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3 heap =
      .ok result) :
    ∃ cPoly dPoly : Polynomial (ZMod this._p.toNat),
      RawDensePolyRep this result.heap result.c0 result.lenC0 cPoly ∧
      RawDensePolyRep this result.heap result.d0 result.lenD0 dPoly := by
  have hlayout := hgcdRecursiveMiddle_layout this q d a2 b2 lenA2 lenB2 m
    W3 heap result hrun
  have hkB : result.k ≤ lenB2 := by
    rw [hlayout.2.2.2.1] at hC0Pos
    split at hC0Pos
    next hk => exact hk
    next hk => simp at hC0Pos
  rcases rawDensePolyRep_split_suffix this result.heap b2 lenB2 result.k
      divisor hkB hDivisor with ⟨cLow, cPoly, hCLow, hCRep, hCSplit⟩
  have hCResult : RawDensePolyRep this result.heap result.c0 result.lenC0
      cPoly := by
    rw [hlayout.2.2.1, hlayout.2.2.2.1]
    simp [hkB]
    exact hCRep
  by_cases hkD : result.k ≤ result.lenD
  · rcases rawDensePolyRep_split_suffix this result.heap d result.lenD
        result.k remainder hkD hRemainder with
      ⟨dLow, dPoly, hDLow, hDRep, hDSplit⟩
    refine ⟨cPoly, dPoly, hCResult, ?_⟩
    rw [hlayout.2.2.2.2.1, hlayout.2.2.2.2.2]
    simp [hkD]
    exact hDRep
  · have hLenD0 : result.lenD0 = 0 := by
      rw [hlayout.2.2.2.2.2]
      simp [hkD]
    have hZero := rawDensePolyRep_zero_length this result.heap result.d0
      (by simpa [hLenD0] using hD0Valid)
    exact ⟨cPoly, 0, hCResult, by simpa [hLenD0] using hZero⟩

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
      dividend = quotient * divisor + remainder ∧
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
      hsemantic, hQRep, hBRep1, hRRep, hdivision, hgcd, _, _, _, _⟩
  have heq : (semanticHeap, semanticLenQ, semanticLenR) =
      (heap1, lenQ, lenR) := Except.ok.inj (hsemantic.symm.trans hdiv)
  cases heq
  exact ⟨heap1, lenQ, lenR, quotient, remainder, row23, row01, hdiv,
    hQRep, hBRep1, hRRep, hdivision, hgcd, hrow23, hrow01, htail, hlt⟩

/-- The actual generated divrem call preserves all four live HGCD matrix
entries when its three write allocations are physically distinct from every
matrix allocation. -/
theorem polyDivrem_preserves_hgcdMatRawDenseRep (this : DenseUPolyZp)
    (M : HgcdMat) (Q R A B : RawPtr UInt64) (lenA lenB : Nat)
    (W3 : RawPtr Word3) (heap heap' : RawHeap) (lenQ lenR : Nat)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat)) (hM : M.Valid)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hQ : heap.ValidU64Slice Q (lenA - (lenB - 1)))
    (hR : heap.ValidU64Slice R (Nat.min lenA (lenB - 1)))
    (hW3 : heap.ValidWord3Slice W3 lenA)
    (hQMatrix : ∀ i : Fin 4,
      Q.region ≠ (hgcdMatPtr M hM i).region)
    (hRMatrix : ∀ i : Fin 4,
      R.region ≠ (hgcdMatPtr M hM i).region)
    (hW3Matrix : ∀ i : Fin 4,
      W3.region ≠ (hgcdMatPtr M hM i).region)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hrun : Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q R
      A lenA B lenB W3 heap = .ok (heap', lenQ, lenR)) :
    HgcdMatRawDenseRep this heap' M entries hM := by
  intro i
  rcases polyDivrem_preserves_u64_region_ne this Q R A B
      (hgcdMatPtr M hM i) lenA lenB (hgcdMatLen M hM i) W3 heap heap'
      lenQ lenR hA hB hQ hR hW3 (hMatrix i).1 (hQMatrix i)
      (hRMatrix i) (hW3Matrix i) hrun with ⟨hlayout, hsame⟩
  exact rawDensePolyRep_of_same_prefix this heap heap'
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i) (entries i) hlayout hsame
    (hMatrix i)

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

/-- Extra physical separation facts for one live polynomial outside the
matrix.  These facts let the exact multiplication/addition execution frame
that polynomial without assigning it any specification-level behavior. -/
structure MatRowUpdateGuardWorkspace (M : HgcdMat) (i0 : Fin 4)
    (Q : RawPtr UInt64) (lenQ : Nat) (T t scratch guard : RawPtr UInt64)
    (guardLen : Nat) (hM : M.Valid) : Prop where
  disjointT : U64SlicesDisjoint T
    (2 * max lenQ (hgcdMatLen M hM i0) - 1) guard guardLen
  disjointScratch : U64SlicesDisjoint scratch
    (8 * max lenQ (hgcdMatLen M hM i0)) guard guardLen
  disjointAddOutput : t.region ≠ guard.region

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

/-- Descriptor-length effect of the same real row update.  The active branch
uses the actual `_mul` heap transition and the actual `_poly_add`
normalization bound; the inactive branch is the source's pure swap. -/
theorem matRowUpdate_length_bound_of_workspace (this : DenseUPolyZp)
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
      hgcdMatLen result.matrix hResult i0 ≤
        max (hgcdMatLen M hM i1)
          (lenQ + hgcdMatLen M hM i0 - 1) ∧
      hgcdMatLen result.matrix hResult i1 = hgcdMatLen M hM i0 ∧
      ∀ j : Fin 4, j ≠ i0 → j ≠ i1 →
        hgcdMatLen result.matrix hResult j = hgcdMatLen M hM j := by
  by_cases hzero : lenQ = 0 ∨ hgcdMatLen M hM i0 = 0
  · rcases matRowUpdate_zero_exec this M i0 i1 Q lenQ T lenT t scratch
      heap hM hne hzero with
    ⟨actual, hactual, _, _, _, _, hActual, _, hlen0, _, hlen1, hother⟩
    have heq : actual = result := Except.ok.inj (hactual.symm.trans hrun)
    subst actual
    exact ⟨hActual, hlen0.le.trans (le_max_left _ _), hlen1,
      fun j hj0 hj1 => (hother j hj0 hj1).2⟩
  · have hQPos : 0 < lenQ := by omega
    have hEntryPos : 0 < hgcdMatLen M hM i0 := by omega
    rcases matRowUpdate_nonzero_success_shape this M i0 i1 Q lenQ T lenT
      t scratch heap result hM hne (by omega) (by omega) hrun with
      ⟨heap1, heap2, sumLen, hmul, hadd, _, _, _, _, hResult,
        _, hlen0, _, hlen1, hother⟩
    have hProduct := matRowUpdate_mul_result this M hM i0 Q T scratch lenQ
      heap heap1 quotient entry0 hcfg hp hQPos hEntryPos workspace.lenWord
      workspace.validT workspace.validScratch workspace.disjointTQ
      (workspace.disjointTMatrix i0) workspace.disjointTScratch
      workspace.disjointScratchQ (workspace.disjointScratchMatrix i0)
      hQRep hEntry0Rep hmul
    have hOutput1 := (hProduct.1 _ _).mp workspace.validAddOutput
    have hEntry1Valid1 :=
      (hProduct.1 (hgcdMatPtr M hM i1) (hgcdMatLen M hM i1)).mp
        hEntry1Rep.1
    rcases polyAdd_ok this t (hgcdMatPtr M hM i1)
        (hgcdMatLen M hM i1) T
        (lenQ + hgcdMatLen M hM i0 - 1) heap1 hOutput1
        hEntry1Valid1 hProduct.2.1 with
      ⟨heapAdd, lenAdd, haddAdd, _, hsum⟩
    have heqAdd : (heapAdd, lenAdd) = (heap2, sumLen) :=
      Except.ok.inj (haddAdd.symm.trans hadd)
    cases heqAdd
    exact ⟨hResult, hlen0.le.trans hsum, hlen1,
      fun j hj0 hj1 => (hother j hj0 hj1).2⟩

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

/-- Apply the general frame theorem using the packaged row workspace and the
three additional separation facts for an external live polynomial. -/
theorem matRowUpdate_preserves_guard_of_workspaces
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch guard : RawPtr UInt64)
    (guardLen : Nat) (heap : RawHeap) (result : MatRowUpdateResult)
    (hM : M.Valid)
    (quotient entry0 entry1 guardPoly : Polynomial (ZMod this._p.toNat))
    (hne : i0 ≠ i1) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (workspace : MatRowUpdateWorkspace M i0 i1 Q lenQ T t scratch heap hM)
    (guardWorkspace : MatRowUpdateGuardWorkspace M i0 Q lenQ T t scratch
      guard guardLen hM)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1)
    (hGuardRep : RawDensePolyRep this heap guard guardLen guardPoly)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result) :
    RawDensePolyRep this result.heap guard guardLen guardPoly := by
  exact matRowUpdate_preserves_guard this M i0 i1 Q lenQ T lenT t scratch
    guard guardLen heap result hM quotient entry0 entry1 guardPoly hne hcfg
    hp workspace.lenWord workspace.validT workspace.validScratch
    workspace.validAddOutput workspace.disjointTQ
    (workspace.disjointTMatrix i0) (workspace.disjointTMatrix i1)
    workspace.disjointTScratch workspace.disjointScratchQ
    (workspace.disjointScratchMatrix i0)
    (workspace.disjointScratchMatrix i1) guardWorkspace.disjointT
    guardWorkspace.disjointScratch guardWorkspace.disjointAddOutput hQRep
    hEntry0Rep hEntry1Rep hGuardRep hrun

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

/-- An arbitrary external live raw polynomial survives both actual generated
row updates.  The proof follows the first call, reconstructs the unchanged
`(0,1)` descriptors and representations, and then follows the second call. -/
theorem hgcdTwoRowUpdates_preserve_guard (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch guard : RawPtr UInt64)
    (guardLen : Nat) (heap : RawHeap) (row23 row01 : MatRowUpdateResult)
    (hM : M.Valid) (h23 : row23.matrix.Valid)
    (quotient guardPoly : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (workspace23 : MatRowUpdateWorkspace M (2 : Fin 4) (3 : Fin 4)
      Q lenQ T t scratch heap hM)
    (guard23 : MatRowUpdateGuardWorkspace M (2 : Fin 4) Q lenQ T t
      scratch guard guardLen hM)
    (workspace01 : MatRowUpdateWorkspace row23.matrix (0 : Fin 4)
      (1 : Fin 4) Q lenQ row23.T row23.t scratch row23.heap h23)
    (guard01 : MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4) Q
      lenQ row23.T row23.t scratch guard guardLen h23)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hGuardRep : RawDensePolyRep this heap guard guardLen guardPoly)
    (hrow23 : dense_upoly_zp__mat_row_update_ir this M
      (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch heap = .ok row23)
    (hrow01 : dense_upoly_zp__mat_row_update_ir this row23.matrix
      (0 : Fin 4) (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t
      scratch row23.heap = .ok row01) :
    RawDensePolyRep this row01.heap guard guardLen guardPoly := by
  have hE0 := hMatrix (0 : Fin 4)
  have hE1 := hMatrix (1 : Fin 4)
  have hE2 := hMatrix (2 : Fin 4)
  have hE3 := hMatrix (3 : Fin 4)
  have hGuard23 := matRowUpdate_preserves_guard_of_workspaces this M
    (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch guard guardLen heap
    row23 hM quotient (entries 2) (entries 3) guardPoly (by decide) hcfg hp
    workspace23 guard23 hQRep hE2 hE3 hGuardRep hrow23
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
    ⟨h23Frame, _, hframe23, _⟩
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
  exact matRowUpdate_preserves_guard_of_workspaces this row23.matrix
    (0 : Fin 4) (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t scratch
    guard guardLen row23.heap row01 h23 quotient (entries 0) (entries 1)
    guardPoly (by decide) hcfg hp workspace01 guard01 hQ23 hE0_23 hE1_23
    hGuard23 hrow01

/-- Descriptor-length effect of the two source-ordered matrix row updates.
This theorem follows both real generated calls; in particular the second
bound uses the quotient and the untouched `(0,1)` entries framed through the
first heap transition. -/
theorem hgcdTwoRowUpdates_length_bounds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (row23 row01 : MatRowUpdateResult)
    (hM : M.Valid) (h23 : row23.matrix.Valid)
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
      hgcdMatLen row01.matrix h01 (0 : Fin 4) ≤
        max (hgcdMatLen M hM (1 : Fin 4))
          (lenQ + hgcdMatLen M hM (0 : Fin 4) - 1) ∧
      hgcdMatLen row01.matrix h01 (1 : Fin 4) =
        hgcdMatLen M hM (0 : Fin 4) ∧
      hgcdMatLen row01.matrix h01 (2 : Fin 4) ≤
        max (hgcdMatLen M hM (3 : Fin 4))
          (lenQ + hgcdMatLen M hM (2 : Fin 4) - 1) ∧
      hgcdMatLen row01.matrix h01 (3 : Fin 4) =
        hgcdMatLen M hM (2 : Fin 4) := by
  have hE0 := hMatrix (0 : Fin 4)
  have hE1 := hMatrix (1 : Fin 4)
  have hE2 := hMatrix (2 : Fin 4)
  have hE3 := hMatrix (3 : Fin 4)
  rcases matRowUpdate_length_bound_of_workspace this M (2 : Fin 4)
      (3 : Fin 4) Q lenQ T lenT t scratch heap row23 hM quotient
      (entries 2) (entries 3) (by decide) hcfg hp workspace23 hQRep hE2 hE3
      hrow23 with ⟨h23', hlen2, hlen3, hother23⟩
  have hh23 : h23' = h23 := Subsingleton.elim _ _
  subst h23'
  have hQ23 := matRowUpdate_preserves_quotient_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) Q lenQ T lenT t scratch heap row23 hM
    quotient (entries 2) (entries 3) (by decide) hcfg hp workspace23
    hQRep hE2 hE3 hrow23
  have hE0_23 := matRowUpdate_preserves_matrix_entry_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) (0 : Fin 4) Q lenQ T lenT t scratch heap
    row23 hM quotient (entries 2) (entries 3) (entries 0) (by decide)
    (by decide) hcfg hp workspace23 hQRep hE2 hE3 hE0 hrow23
  have hE1_23 := matRowUpdate_preserves_matrix_entry_of_workspace this M
    (2 : Fin 4) (3 : Fin 4) (1 : Fin 4) Q lenQ T lenT t scratch heap
    row23 hM quotient (entries 2) (entries 3) (entries 1) (by decide)
    (by decide) hcfg hp workspace23 hQRep hE2 hE3 hE1 hrow23
  rcases hgcdTwoRowUpdates_descriptor_frame this M Q lenQ T lenT t scratch
      heap row23 row01 hM hrow23 hrow01 with
    ⟨h23Frame, _, hframe23, _⟩
  have hh23Frame : h23Frame = h23 := Subsingleton.elim _ _
  subst h23Frame
  have hE0_23' : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (0 : Fin 4))
      (hgcdMatLen row23.matrix h23 (0 : Fin 4)) (entries 0) := by
    have hptr : hgcdMatPtr row23.matrix h23 (0 : Fin 4) =
        hgcdMatPtr M hM (0 : Fin 4) := by
      simpa only using (hframe23 (0 : Fin 2)).1
    have hlen : hgcdMatLen row23.matrix h23 (0 : Fin 4) =
        hgcdMatLen M hM (0 : Fin 4) := by
      simpa only using (hframe23 (0 : Fin 2)).2
    rw [hptr, hlen]
    exact hE0_23
  have hE1_23' : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (1 : Fin 4))
      (hgcdMatLen row23.matrix h23 (1 : Fin 4)) (entries 1) := by
    have hptr : hgcdMatPtr row23.matrix h23 (1 : Fin 4) =
        hgcdMatPtr M hM (1 : Fin 4) := by
      simpa only using (hframe23 (1 : Fin 2)).1
    have hlen : hgcdMatLen row23.matrix h23 (1 : Fin 4) =
        hgcdMatLen M hM (1 : Fin 4) := by
      simpa only using (hframe23 (1 : Fin 2)).2
    rw [hptr, hlen]
    exact hE1_23
  rcases matRowUpdate_length_bound_of_workspace this row23.matrix
      (0 : Fin 4) (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t scratch
      row23.heap row01 h23 quotient (entries 0) (entries 1) (by decide)
      hcfg hp workspace01 hQ23 hE0_23' hE1_23' hrow01 with
    ⟨h01, hlen0, hlen1, hother01⟩
  have h23len0 : hgcdMatLen row23.matrix h23 (0 : Fin 4) =
      hgcdMatLen M hM (0 : Fin 4) := hother23 (0 : Fin 4) (by decide)
        (by decide)
  have h23len1 : hgcdMatLen row23.matrix h23 (1 : Fin 4) =
      hgcdMatLen M hM (1 : Fin 4) := hother23 (1 : Fin 4) (by decide)
        (by decide)
  have h01len2 := hother01 (2 : Fin 4) (by decide) (by decide)
  have h01len3 := hother01 (3 : Fin 4) (by decide) (by decide)
  have hfinal0 : hgcdMatLen row01.matrix h01 (0 : Fin 4) ≤
      max (hgcdMatLen M hM (1 : Fin 4))
        (lenQ + hgcdMatLen M hM (0 : Fin 4) - 1) := by
    simpa only [h23len0, h23len1] using hlen0
  have hfinal1 : hgcdMatLen row01.matrix h01 (1 : Fin 4) =
      hgcdMatLen M hM (0 : Fin 4) := by
    simpa only [h23len0] using hlen1
  have hfinal2 : hgcdMatLen row01.matrix h01 (2 : Fin 4) ≤
      max (hgcdMatLen M hM (3 : Fin 4))
        (lenQ + hgcdMatLen M hM (2 : Fin 4) - 1) := h01len2 ▸ hlen2
  have hfinal3 : hgcdMatLen row01.matrix h01 (3 : Fin 4) =
      hgcdMatLen M hM (2 : Fin 4) := h01len3.trans hlen3
  exact ⟨h01, hfinal0, hfinal1, hfinal2, hfinal3⟩

/-- One complete nonterminal HGCD Euclidean iteration, from the actual
generated divrem call through both actual row updates, preserves the raw
state, matrix transform, signed determinant, normalized gcd, and the strict
well-founded measure. -/
theorem hgcdIterationCalls_refine (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch A B R : RawPtr UInt64) (lenA lenB : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t : RawPtr UInt64)
    (heap heap1 : RawHeap) (lenQ lenR : Nat)
    (row23 row01 : MatRowUpdateResult) (hM : M.Valid)
    (h23 : row23.matrix.Valid)
    (left right dividend divisor : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (sgn : Int) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat) (hlenB : 0 < lenB)
    (hARep : RawDensePolyRep this heap A lenA dividend)
    (hBRep : RawDensePolyRep this heap B lenB divisor)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (htransform : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right
      dividend divisor (entries 0) (entries 1) (entries 2) (entries 3))
    (hdet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet sgn
      (entries 0) (entries 1) (entries 2) (entries 3))
    (hQ : heap.ValidU64Slice Q (lenA - (lenB - 1)))
    (hR : heap.ValidU64Slice R (Nat.min lenA (lenB - 1)))
    (hW3 : heap.ValidWord3Slice W3 lenA)
    (hqCapacity : lenA - (lenB - 1) < limbBase)
    (hRA : R.region ≠ A.region) (hWA : W3.region ≠ A.region)
    (hWB : W3.region ≠ B.region) (hQB : Q.region ≠ B.region)
    (hQW : Q.region ≠ W3.region) (hRW : R.region ≠ W3.region)
    (hRQ : R.region ≠ Q.region) (hRB : R.region ≠ B.region)
    (hQMatrix : ∀ i : Fin 4, Q.region ≠ (hgcdMatPtr M hM i).region)
    (hRMatrix : ∀ i : Fin 4, R.region ≠ (hgcdMatPtr M hM i).region)
    (hW3Matrix : ∀ i : Fin 4, W3.region ≠ (hgcdMatPtr M hM i).region)
    (workspace23 : MatRowUpdateWorkspace M (2 : Fin 4) (3 : Fin 4)
      Q lenQ T t scratch heap1 hM)
    (workspace01 : MatRowUpdateWorkspace row23.matrix (0 : Fin 4)
      (1 : Fin 4) Q lenQ row23.T row23.t scratch row23.heap h23)
    (divisorGuard23 : MatRowUpdateGuardWorkspace M (2 : Fin 4) Q lenQ
      T t scratch B lenB hM)
    (divisorGuard01 : MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4)
      Q lenQ row23.T row23.t scratch B lenB h23)
    (remainderGuard23 : MatRowUpdateGuardWorkspace M (2 : Fin 4) Q lenQ
      T t scratch R lenR hM)
    (remainderGuard01 : MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4)
      Q lenQ row23.T row23.t scratch R lenR h23)
    (hdiv : Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q R
      A lenA B lenB W3 heap = .ok (heap1, lenQ, lenR))
    (hrow23 : dense_upoly_zp__mat_row_update_ir this M (2 : Fin 4)
      (3 : Fin 4) Q lenQ T lenT t scratch heap1 = .ok row23)
    (hrow01 : dense_upoly_zp__mat_row_update_ir this row23.matrix
      (0 : Fin 4) (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t scratch
      row23.heap = .ok row01) :
    ∃ quotient remainder h01,
      RawDensePolyRep this heap1 Q lenQ quotient ∧
      HgcdMatRawDenseRep this heap1 M entries hM ∧
      HgcdMatRawDenseRep this row01.heap row01.matrix
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
        h01 ∧
      RawDensePolyRep this row01.heap B lenB divisor ∧
      RawDensePolyRep this row01.heap R lenR remainder ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right divisor
        remainder
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 0)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 1)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 2)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 3) ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet (-sgn)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 0)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 1)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 2)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 3) ∧
      normalize (EuclideanDomain.gcd dividend divisor) =
        normalize (EuclideanDomain.gcd divisor remainder) ∧
      lenQ ≤ lenA - (lenB - 1) ∧ lenR < lenB := by
  rcases polyDivrem_next_state this Q R A B lenA lenB W3 heap dividend
      divisor hlenB hARep hBRep hQ hR hW3 hqCapacity hRA hWA hWB hQB hQW
      hRW hRQ hRB hcfg with
    ⟨semanticHeap, semanticLenQ, semanticLenR, quotient, remainder,
      hsemantic, hQRep, hBRep1, hRRep1, hdivision, hgcd, _, hlenQ, _, hlt⟩
  have heq : (semanticHeap, semanticLenQ, semanticLenR) =
      (heap1, lenQ, lenR) := Except.ok.inj (hsemantic.symm.trans hdiv)
  cases heq
  have hMatrix1 := polyDivrem_preserves_hgcdMatRawDenseRep this M Q R A B
    lenA lenB W3 heap heap1 lenQ lenR entries hM hARep.1 hBRep.1 hQ hR
    hW3 hQMatrix hRMatrix hW3Matrix hMatrix hdiv
  rcases hgcdTwoRowUpdates_refine_matrix this M Q lenQ T lenT t scratch
      heap1 row23 row01 hM h23 quotient entries hcfg hp workspace23
      workspace01 hQRep hMatrix1 hrow23 hrow01 with ⟨h01, hMatrix01⟩
  have hDivisor01 := hgcdTwoRowUpdates_preserve_guard this M Q lenQ T
    lenT t scratch B lenB heap1 row23 row01 hM h23 quotient divisor entries
    hcfg hp workspace23 divisorGuard23 workspace01 divisorGuard01 hQRep
    hMatrix1 hBRep1 hrow23 hrow01
  have hRemainder01 := hgcdTwoRowUpdates_preserve_guard this M Q lenQ T
    lenT t scratch R lenR heap1 row23 row01 hM h23 quotient remainder
    entries hcfg hp workspace23 remainderGuard23 workspace01
    remainderGuard01 hQRep hMatrix1 hRRep1 hrow23 hrow01
  have htransform' :=
    CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries_preserves_transform
      left right dividend divisor remainder quotient entries htransform
      hdivision
  have hdet' :=
    CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries_preserves_signedDet
      sgn quotient entries hdet
  exact ⟨quotient, remainder, h01, hQRep, hMatrix1, hMatrix01, hDivisor01,
    hRemainder01,
    htransform', hdet', hgcd, hlenQ, hlt⟩

/-- Semantic invariant carried by the well-founded generated HGCD loop.
The original pair is related to the current raw pair by the actual matrix,
whose determinant sign is synchronized with the source `sgn` field. -/
structure HgcdIterRawInvariant (this : DenseUPolyZp)
    (left right currentA currentB : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (state : HgcdIterState) (hM : state.matrix.Valid) : Prop where
  matrixRep : HgcdMatRawDenseRep this state.heap state.matrix entries hM
  aRep : RawDensePolyRep this state.heap state.A state.lenA currentA
  bRep : RawDensePolyRep this state.heap state.B state.lenB currentB
  transform : CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right
    currentA currentB (entries 0) (entries 1) (entries 2) (entries 3)
  signedDet : CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet state.sgn
    (entries 0) (entries 1) (entries 2) (entries 3)

/-- Descriptor-length invariant needed by the enclosing recursive HGCD.
Each matrix coefficient is paired with the current Euclidean operand whose
degree it complements. -/
structure HgcdMatrixLengthInvariant (inputLength : Nat)
    (state : HgcdIterState) (hM : state.matrix.Valid) : Prop where
  row0A : hgcdMatLen state.matrix hM (0 : Fin 4) + state.lenA ≤
    inputLength + 1
  row1B : hgcdMatLen state.matrix hM (1 : Fin 4) + state.lenB ≤
    inputLength + 1
  row2A : hgcdMatLen state.matrix hM (2 : Fin 4) + state.lenA ≤
    inputLength + 1
  row3B : hgcdMatLen state.matrix hM (3 : Fin 4) + state.lenB ≤
    inputLength + 1

/-- Arithmetic closure for one Euclidean matrix step.  The quotient bound
comes from the real generated divrem call and the four descriptor bounds
come from the two real row updates. -/
theorem hgcdMatrixLengthInvariant_step
    (inputLength lenA lenB lenQ lenR : Nat)
    (M nextM : HgcdMat) (hM : M.Valid) (hNextM : nextM.Valid)
    (hinvariant :
      hgcdMatLen M hM (0 : Fin 4) + lenA ≤ inputLength + 1 ∧
      hgcdMatLen M hM (1 : Fin 4) + lenB ≤ inputLength + 1 ∧
      hgcdMatLen M hM (2 : Fin 4) + lenA ≤ inputLength + 1 ∧
      hgcdMatLen M hM (3 : Fin 4) + lenB ≤ inputLength + 1)
    (horder : lenB ≤ lenA)
    (hlenQ : lenQ ≤ lenA - (lenB - 1)) (hlenR : lenR < lenB)
    (hlen0 : hgcdMatLen nextM hNextM (0 : Fin 4) ≤
      max (hgcdMatLen M hM (1 : Fin 4))
        (lenQ + hgcdMatLen M hM (0 : Fin 4) - 1))
    (hlen1 : hgcdMatLen nextM hNextM (1 : Fin 4) =
      hgcdMatLen M hM (0 : Fin 4))
    (hlen2 : hgcdMatLen nextM hNextM (2 : Fin 4) ≤
      max (hgcdMatLen M hM (3 : Fin 4))
        (lenQ + hgcdMatLen M hM (2 : Fin 4) - 1))
    (hlen3 : hgcdMatLen nextM hNextM (3 : Fin 4) =
      hgcdMatLen M hM (2 : Fin 4)) :
    hgcdMatLen nextM hNextM (0 : Fin 4) + lenB ≤ inputLength + 1 ∧
    hgcdMatLen nextM hNextM (1 : Fin 4) + lenR ≤ inputLength + 1 ∧
    hgcdMatLen nextM hNextM (2 : Fin 4) + lenB ≤ inputLength + 1 ∧
    hgcdMatLen nextM hNextM (3 : Fin 4) + lenR ≤ inputLength + 1 := by
  rcases hinvariant with ⟨h0A, h1B, h2A, h3B⟩
  have hq0 : lenQ + hgcdMatLen M hM (0 : Fin 4) - 1 + lenB ≤
      inputLength + 1 := by omega
  have hq2 : lenQ + hgcdMatLen M hM (2 : Fin 4) - 1 + lenB ≤
      inputLength + 1 := by omega
  have hmax0 : max (hgcdMatLen M hM (1 : Fin 4))
        (lenQ + hgcdMatLen M hM (0 : Fin 4) - 1) + lenB ≤
      inputLength + 1 := by
    rcases le_total (hgcdMatLen M hM (1 : Fin 4))
        (lenQ + hgcdMatLen M hM (0 : Fin 4) - 1) with hle | hle
    · rw [max_eq_right hle]
      exact hq0
    · rw [max_eq_left hle]
      omega
  have hmax2 : max (hgcdMatLen M hM (3 : Fin 4))
        (lenQ + hgcdMatLen M hM (2 : Fin 4) - 1) + lenB ≤
      inputLength + 1 := by
    rcases le_total (hgcdMatLen M hM (3 : Fin 4))
        (lenQ + hgcdMatLen M hM (2 : Fin 4) - 1) with hle | hle
    · rw [max_eq_right hle]
      exact hq2
    · rw [max_eq_left hle]
      omega
  constructor
  · omega
  constructor
  · omega
  constructor <;> omega

/-- The real `_mat_one` descriptor and the two source-ordered input copies
establish the matrix-length invariant. -/
theorem hgcdIterInit_matrixLengthInvariant
    (M : HgcdMat) (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (initial : HgcdIterState) (hM : initial.matrix.Valid)
    (horder : lenB ≤ lenA)
    (hrun : hgcdIterInit M A B T t lenT a lenA b lenB heap = .ok initial) :
    HgcdMatrixLengthInvariant lenA initial hM := by
  have hlens := hgcdIterInit_lengths M A B T t lenT a lenA b lenB heap
    initial hrun
  rcases hlens with ⟨hA, hB, hMatrix⟩
  constructor <;>
    simp [hgcdMatLen, hA, hB, hMatrix] <;> omega

/-- Purely physical obligations for one concrete successful nonterminal
iteration.  Every field concerns validity, capacity, or allocation
separation of the exact buffers used by the generated calls. -/
structure HgcdIterationWorkspace (this : DenseUPolyZp)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (state : HgcdIterState) (heap1 : RawHeap) (lenQ lenR : Nat)
    (row23 row01 : MatRowUpdateResult) (hM : state.matrix.Valid) : Prop where
  validQ : state.heap.ValidU64Slice Q
    (state.lenA - (state.lenB - 1))
  validR : state.heap.ValidU64Slice state.T
    (Nat.min state.lenA (state.lenB - 1))
  validW3 : state.heap.ValidWord3Slice W3 state.lenA
  quotientCapacity : state.lenA - (state.lenB - 1) < limbBase
  rA : state.T.region ≠ state.A.region
  wA : W3.region ≠ state.A.region
  wB : W3.region ≠ state.B.region
  qB : Q.region ≠ state.B.region
  qW : Q.region ≠ W3.region
  rW : state.T.region ≠ W3.region
  rQ : state.T.region ≠ Q.region
  rB : state.T.region ≠ state.B.region
  qMatrix : ∀ i : Fin 4,
    Q.region ≠ (hgcdMatPtr state.matrix hM i).region
  rMatrix : ∀ i : Fin 4,
    state.T.region ≠ (hgcdMatPtr state.matrix hM i).region
  wMatrix : ∀ i : Fin 4,
    W3.region ≠ (hgcdMatPtr state.matrix hM i).region
  matrix23Valid : row23.matrix.Valid
  row23Workspace : MatRowUpdateWorkspace state.matrix (2 : Fin 4)
    (3 : Fin 4) Q lenQ state.A state.t scratch heap1 hM
  row01Workspace : MatRowUpdateWorkspace row23.matrix (0 : Fin 4)
    (1 : Fin 4) Q lenQ row23.T row23.t scratch row23.heap matrix23Valid
  divisorGuard23 : MatRowUpdateGuardWorkspace state.matrix (2 : Fin 4)
    Q lenQ state.A state.t scratch state.B state.lenB hM
  divisorGuard01 : MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4)
    Q lenQ row23.T row23.t scratch state.B state.lenB matrix23Valid
  remainderGuard23 : MatRowUpdateGuardWorkspace state.matrix (2 : Fin 4)
    Q lenQ state.A state.t scratch state.T lenR hM
  remainderGuard01 : MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4)
    Q lenQ row23.T row23.t scratch state.T lenR matrix23Valid

/-- Separation-logic style precondition for the whole loop: for every
concrete successful source step, the caller can discharge the purely
physical workspace obligations above.  It contains no polynomial or L2
result. -/
def HgcdLoopWorkspaceProvider (this : DenseUPolyZp) (m : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch : RawPtr UInt64) : Prop :=
  ∀ (state : HgcdIterState) (hM : state.matrix.Valid)
    (heap1 : RawHeap) (lenQ lenR : Nat)
    (row23 row01 : MatRowUpdateResult),
    state.lenB ≥ m + 1 →
    Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
      state.A state.lenA state.B state.lenB W3 state.heap =
        .ok (heap1, lenQ, lenR) →
    dense_upoly_zp__mat_row_update_ir this state.matrix (2 : Fin 4)
      (3 : Fin 4) Q lenQ state.A state.lenT state.t scratch heap1 =
        .ok row23 →
    dense_upoly_zp__mat_row_update_ir this row23.matrix (0 : Fin 4)
      (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t scratch row23.heap =
        .ok row01 →
    HgcdIterationWorkspace this Q W3 scratch state heap1 lenQ lenR row23
      row01 hM

/-- End-to-end semantic refinement of the actual well-founded generated
`hgcdIterLoop`.  Recursion is on the source measure `state.lenB`; each
recursive call is justified by the exact divrem result's `lenR < lenB`; no
bounded execution counter or alternate L2 execution is introduced. -/
theorem hgcdIterLoop_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (m : Nat) (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch : RawPtr UInt64)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdLoopWorkspaceProvider this m Q W3 scratch) :
    ∀ (state final : HgcdIterState)
      (currentA currentB : Polynomial (ZMod this._p.toNat))
      (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
      (hM : state.matrix.Valid),
      HgcdIterRawInvariant this left right currentA currentB entries state hM →
      hgcdIterLoop this m Q W3 scratch state = .ok final →
      ∃ finalA finalB finalEntries hFinalM,
        HgcdIterRawInvariant this left right finalA finalB finalEntries final
          hFinalM ∧
        normalize (EuclideanDomain.gcd currentA currentB) =
          normalize (EuclideanDomain.gcd finalA finalB) ∧
        final.lenB < m + 1
  | state, final, currentA, currentB, entries, hM, hinvariant, hrun => by
    by_cases hguard : state.lenB ≥ m + 1
    · rcases hgcdIterLoop_step_shape this m Q W3 scratch state final hguard
        hrun with
      ⟨heap1, lenQ, lenR, row23, row01, hdiv, hrow23, hrow01, htail, hlt⟩
      have hworkspace := physical state hM heap1 lenQ lenR row23 row01
        hguard hdiv hrow23 hrow01
      rcases hgcdIterationCalls_refine this state.matrix Q W3 scratch state.A
          state.B state.T state.lenA state.lenB state.A state.lenT state.t
          state.heap heap1 lenQ lenR row23 row01 hM
          hworkspace.matrix23Valid left right currentA currentB entries
          state.sgn hcfg hp (by omega) hinvariant.aRep hinvariant.bRep
          hinvariant.matrixRep hinvariant.transform hinvariant.signedDet
          hworkspace.validQ hworkspace.validR hworkspace.validW3
          hworkspace.quotientCapacity hworkspace.rA hworkspace.wA
          hworkspace.wB hworkspace.qB hworkspace.qW hworkspace.rW
          hworkspace.rQ hworkspace.rB hworkspace.qMatrix
          hworkspace.rMatrix hworkspace.wMatrix hworkspace.row23Workspace
          hworkspace.row01Workspace hworkspace.divisorGuard23
          hworkspace.divisorGuard01 hworkspace.remainderGuard23
          hworkspace.remainderGuard01 hdiv hrow23 hrow01 with
        ⟨quotient, remainder, h01, _, _, hMatrix01, hDivisor01, hRemainder01,
          htransform, hdet, hgcdStep, _, hlt'⟩
      let next : HgcdIterState := {
        heap := row01.heap
        matrix := row01.matrix
        A := state.B
        lenA := state.lenB
        B := state.T
        lenB := lenR
        T := row01.T
        lenT := row01.lenT
        t := row01.t
        sgn := -state.sgn }
      have hnextInvariant : HgcdIterRawInvariant this left right currentB
          remainder
          (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
          next h01 := by
        exact ⟨hMatrix01, hDivisor01, hRemainder01, htransform, hdet⟩
      have htail' : hgcdIterLoop this m Q W3 scratch next = .ok final := by
        simpa [next] using htail
      rcases hgcdIterLoop_refines this m Q W3 scratch left right hcfg hp
          physical next final currentB remainder
          (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
          h01 hnextInvariant htail' with
        ⟨finalA, finalB, finalEntries, hFinalM, hfinalInvariant,
          hgcdRest, hstop⟩
      exact ⟨finalA, finalB, finalEntries, hFinalM, hfinalInvariant,
        hgcdStep.trans hgcdRest, hstop⟩
    · have hstop : state.lenB < m + 1 := by omega
      have hsame := hgcdIterLoop_stop this m Q W3 scratch state hstop
      have hfinal : state = final := Except.ok.inj (hsame.symm.trans hrun)
      subst final
      exact ⟨currentA, currentB, entries, hM, hinvariant, rfl, hstop⟩
termination_by state => state.lenB
decreasing_by exact hlt'

/-- The same well-founded source loop carries the descriptor-length
invariant needed by the enclosing recursive call.  This is deliberately a
separate theorem from semantic refinement so its result cannot be supplied
without following the generated heap execution. -/
theorem hgcdIterLoop_preserves_matrixLength (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (inputLength m : Nat) (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch : RawPtr UInt64)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdLoopWorkspaceProvider this m Q W3 scratch) :
    ∀ (state final : HgcdIterState)
      (currentA currentB : Polynomial (ZMod this._p.toNat))
      (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
      (hM : state.matrix.Valid),
      HgcdIterRawInvariant this left right currentA currentB entries state hM →
      HgcdMatrixLengthInvariant inputLength state hM →
      state.lenB ≤ state.lenA →
      state.lenA ≤ inputLength →
      0 < state.lenA →
      hgcdIterLoop this m Q W3 scratch state = .ok final →
      ∃ hFinalM : final.matrix.Valid,
        HgcdMatrixLengthInvariant inputLength final hFinalM ∧
        final.lenB ≤ final.lenA ∧ final.lenA ≤ inputLength ∧ 0 < final.lenA
  | state, final, currentA, currentB, entries, hM, hraw, hlength, horder,
      hinputBound, hpositive, hrun => by
    by_cases hguard : state.lenB ≥ m + 1
    · rcases hgcdIterLoop_step_shape this m Q W3 scratch state final hguard
        hrun with
      ⟨heap1, lenQ, lenR, row23, row01, hdiv, hrow23, hrow01, htail, hlt⟩
      have hworkspace := physical state hM heap1 lenQ lenR row23 row01
        hguard hdiv hrow23 hrow01
      rcases hgcdIterationCalls_refine this state.matrix Q W3 scratch state.A
          state.B state.T state.lenA state.lenB state.A state.lenT state.t
          state.heap heap1 lenQ lenR row23 row01 hM
          hworkspace.matrix23Valid left right currentA currentB entries
          state.sgn hcfg hp (by omega) hraw.aRep hraw.bRep hraw.matrixRep
          hraw.transform hraw.signedDet hworkspace.validQ hworkspace.validR
          hworkspace.validW3 hworkspace.quotientCapacity hworkspace.rA
          hworkspace.wA hworkspace.wB hworkspace.qB hworkspace.qW
          hworkspace.rW hworkspace.rQ hworkspace.rB hworkspace.qMatrix
          hworkspace.rMatrix hworkspace.wMatrix hworkspace.row23Workspace
          hworkspace.row01Workspace hworkspace.divisorGuard23
          hworkspace.divisorGuard01 hworkspace.remainderGuard23
          hworkspace.remainderGuard01 hdiv hrow23 hrow01 with
        ⟨quotient, remainder, h01, hQRep, hMatrix1, hMatrix01, hDivisor01,
          hRemainder01, htransform, hdet, _, hlenQ, hlenR⟩
      rcases hgcdTwoRowUpdates_length_bounds this state.matrix Q lenQ
          state.A state.lenT state.t scratch heap1 row23 row01 hM
          hworkspace.matrix23Valid quotient entries hcfg hp
          hworkspace.row23Workspace hworkspace.row01Workspace hQRep
          hMatrix1 hrow23 hrow01 with
        ⟨h01', hlen0, hlen1, hlen2, hlen3⟩
      have hh01 : h01' = h01 := Subsingleton.elim _ _
      subst h01'
      let next : HgcdIterState := {
        heap := row01.heap
        matrix := row01.matrix
        A := state.B
        lenA := state.lenB
        B := state.T
        lenB := lenR
        T := row01.T
        lenT := row01.lenT
        t := row01.t
        sgn := -state.sgn }
      have hnextRaw : HgcdIterRawInvariant this left right currentB remainder
          (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
          next h01 := ⟨hMatrix01, hDivisor01, hRemainder01, htransform, hdet⟩
      have hnextBounds := hgcdMatrixLengthInvariant_step inputLength
        state.lenA state.lenB lenQ lenR state.matrix row01.matrix hM h01
        ⟨hlength.row0A, hlength.row1B, hlength.row2A, hlength.row3B⟩
        horder hlenQ hlenR hlen0 hlen1 hlen2 hlen3
      have hnextLength : HgcdMatrixLengthInvariant inputLength next h01 := by
        exact ⟨hnextBounds.1, hnextBounds.2.1, hnextBounds.2.2.1,
          hnextBounds.2.2.2⟩
      have htail' : hgcdIterLoop this m Q W3 scratch next = .ok final := by
        simpa [next] using htail
      exact hgcdIterLoop_preserves_matrixLength this inputLength m Q W3
        scratch left right hcfg hp physical next final currentB remainder
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
        h01 hnextRaw hnextLength (by simpa [next] using hlenR.le)
        (by simpa [next] using horder.trans hinputBound)
        (by simpa [next] using (show 0 < state.lenB by omega)) htail'
    · have hstop : state.lenB < m + 1 := by omega
      have hsame := hgcdIterLoop_stop this m Q W3 scratch state hstop
      have hfinal : state = final := Except.ok.inj (hsame.symm.trans hrun)
      subst final
      exact ⟨hM, hlength, horder, hinputBound, hpositive⟩
termination_by state => state.lenB
decreasing_by exact hlt

/-- End-to-end refinement of generated C++ `_hgcd_iter`: the exact identity
initialization and ordered copies feed the exact well-founded loop theorem. -/
theorem hgcdIter_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (heap : RawHeap) (final : HgcdIterState)
    (left right : Polynomial (ZMod this._p.toNat)) (hM : M.Valid)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenB ≤ lenA)
    (hlenAPos : 0 < lenA)
    (physical : HgcdLoopWorkspaceProvider this (lenA / 2) Q W3 scratch)
    (h0 : heap.ValidU64Slice (hgcdMatPtr M hM (0 : Fin 4)) 1)
    (h3 : heap.ValidU64Slice (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (h03 : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1
      (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (hA : heap.ValidU64Slice A lenA)
    (hB : heap.ValidU64Slice B lenB)
    (hAa : U64SlicesDisjoint A lenA a lenA)
    (hBb : U64SlicesDisjoint B lenB b lenB)
    (hAb : U64SlicesDisjoint A lenA b lenB)
    (hBA : U64SlicesDisjoint B lenB A lenA)
    (h0a : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 a lenA)
    (h3a : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 a lenA)
    (h0b : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 b lenB)
    (h3b : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 b lenB)
    (hAMatrix : ∀ i : Fin 4, U64SlicesDisjoint A lenA
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hBMatrix : ∀ i : Fin 4, U64SlicesDisjoint B lenB
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hMatrixValid : ∀ i : Fin 4, heap.ValidU64Slice
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hLeft : RawDensePolyRep this heap a lenA left)
    (hRight : RawDensePolyRep this heap b lenB right)
    (hrun : dense_upoly_zp__hgcd_iter_ir this M A B T t lenT a lenA b
      lenB Q W3 scratch heap = .ok final) :
    ∃ finalA finalB finalEntries hFinalM,
      HgcdIterRawInvariant this left right finalA finalB finalEntries final
        hFinalM ∧
      normalize (EuclideanDomain.gcd left right) =
        normalize (EuclideanDomain.gcd finalA finalB) ∧
      final.lenB < lenA / 2 + 1 ∧
      HgcdMatrixLengthInvariant lenA final hFinalM ∧
      final.lenB ≤ final.lenA ∧ final.lenA ≤ lenA ∧ 0 < final.lenA := by
  rcases hgcdIterInit_refines this M A B T t lenT a lenA b lenB heap left
      right hM hp h0 h3 h03 hA hB hAa hBb hAb hBA h0a h3a h0b h3b
      hAMatrix hBMatrix hMatrixValid hLeft hRight with
    ⟨initial, hinit, _, _, _, _, _, _, _, _, hInitialM,
      hInitialMatrix, hInitialA, hInitialB, hInitialTransform,
      hInitialDet⟩
  have hloop : hgcdIterLoop this (lenA / 2) Q W3 scratch initial =
      .ok final := by
    simpa [dense_upoly_zp__hgcd_iter_ir, hinit] using hrun
  have hInitialInvariant : HgcdIterRawInvariant this left right left right
      (identityEntries this._p.toNat) initial hInitialM :=
    ⟨hInitialMatrix, hInitialA, hInitialB, hInitialTransform, hInitialDet⟩
  have hInitialLength := hgcdIterInit_matrixLengthInvariant M A B T t lenT
    a lenA b lenB heap initial hInitialM horder hinit
  have hInitialOrder : initial.lenB ≤ initial.lenA := by
    have hlens := hgcdIterInit_lengths M A B T t lenT a lenA b lenB heap
      initial hinit
    omega
  rcases hgcdIterLoop_refines this (lenA / 2) Q W3 scratch left right hcfg hp
      physical initial final left right (identityEntries this._p.toNat)
      hInitialM hInitialInvariant hloop with
    ⟨finalA, finalB, finalEntries, hFinalM, hFinalRaw, hGcd, hStop⟩
  rcases hgcdIterLoop_preserves_matrixLength this lenA (lenA / 2) Q W3
      scratch left right hcfg hp physical initial final left right
      (identityEntries this._p.toNat) hInitialM hInitialInvariant
      hInitialLength hInitialOrder (by
        have hlens := hgcdIterInit_lengths M A B T t lenT a lenA b lenB heap
          initial hinit
        omega) (by
        have hlens := hgcdIterInit_lengths M A B T t lenT a lenA b lenB heap
          initial hinit
        omega) hloop with
    ⟨hFinalM', hFinalLength, hFinalOrder, hFinalInputBound,
      hFinalPositive⟩
  have hhFinal : hFinalM' = hFinalM := Subsingleton.elim _ _
  subst hFinalM'
  exact ⟨finalA, finalB, finalEntries, hFinalM, hFinalRaw, hGcd, hStop,
    hFinalLength, hFinalOrder, hFinalInputBound, hFinalPositive⟩

/-- Purely physical obligations needed after one concrete iterator result to
stabilize its matrix and normalize its two output pointers.  No polynomial
or expected L2 result occurs in this contract. -/
structure HgcdRecursiveIterFinalizeWorkspace
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 stage : RawPtr UInt64)
    (iter : HgcdIterState) (hIter : iter.matrix.Valid)
    (stable : HgcdMatRestoreResult) (hStable : stable.matrix.Valid) : Prop where
  stabilize : HgcdMatStabilizeWorkspace iter.heap original iter.matrix
    hOriginal hIter stage
  stageA : stage.region ≠ iter.A.region
  stageB : stage.region ≠ iter.B.region
  originalA : ∀ j : Fin 4,
    (hgcdMatPtr original hOriginal j).region ≠ iter.A.region
  originalB : ∀ j : Fin 4,
    (hgcdMatPtr original hOriginal j).region ≠ iter.B.region
  validA3 : iter.heap.ValidU64Slice a3 iter.lenA
  validB3 : iter.heap.ValidU64Slice b3 iter.lenB
  pAEq : (iter.A == a3) = true → iter.A = a3
  pBEq : (iter.B == b3) = true → iter.B = b3
  b3PB : U64SlicesDisjoint b3 iter.lenB iter.B iter.lenB
  b3PA : U64SlicesDisjoint b3 iter.lenB iter.A iter.lenA
  a3PA : U64SlicesDisjoint a3 iter.lenA iter.A iter.lenA
  a3PB : U64SlicesDisjoint a3 iter.lenA iter.B iter.lenB
  a3B3 : U64SlicesDisjoint a3 iter.lenA b3 iter.lenB
  a3Matrix : ∀ j : Fin 4, U64SlicesDisjoint a3 iter.lenA
    (hgcdMatPtr stable.matrix hStable j) (hgcdMatLen stable.matrix hStable j)
  b3Matrix : ∀ j : Fin 4, U64SlicesDisjoint b3 iter.lenB
    (hgcdMatPtr stable.matrix hStable j) (hgcdMatLen stable.matrix hStable j)

/-- Separation-logic provider for the exact successful iterator and
stabilization executions exposed by `hgcdRecursiveIterBranch_exec`. -/
def HgcdRecursiveIterFinalizeWorkspaceProvider (this : DenseUPolyZp)
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage : RawPtr UInt64) (heap : RawHeap) : Prop :=
  ∀ (iter : HgcdIterState) (hIter : iter.matrix.Valid)
    (stable : HgcdMatRestoreResult) (hStable : stable.matrix.Valid),
    dense_upoly_zp__hgcd_iter_ir this original a3 b3 T0 T1 0 inputA
        lenInputA inputB lenInputB Q W3 scratch heap = .ok iter →
    hgcdMatStabilize original iter.matrix hOriginal hIter stage iter.heap =
        .ok stable →
    HgcdRecursiveIterFinalizeWorkspace original hOriginal a3 b3 stage iter
      hIter stable hStable

/-- End-to-end refinement of an actual recursive-HGCD iterator arm.  The
proof consumes exactly the generated iterator, stabilization, and
alias-sensitive store executions and returns their raw polynomial meaning,
matrix transform, signed determinant, GCD invariant, and stopping bound. -/
theorem hgcdRecursiveIterBranch_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage : RawPtr UInt64) (heap : RawHeap)
    (result : HgcdRecursiveIterBranchResult)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hInputOrder : lenInputB ≤ lenInputA)
    (hInputAPos : 0 < lenInputA)
    (loopPhysical : HgcdLoopWorkspaceProvider this (lenInputA / 2) Q W3
      scratch)
    (finalizePhysical : HgcdRecursiveIterFinalizeWorkspaceProvider this
      original hOriginal a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1
      scratch stage heap)
    (h0 : heap.ValidU64Slice
      (hgcdMatPtr original hOriginal (0 : Fin 4)) 1)
    (h3 : heap.ValidU64Slice
      (hgcdMatPtr original hOriginal (3 : Fin 4)) 1)
    (h03 : U64SlicesDisjoint
      (hgcdMatPtr original hOriginal (0 : Fin 4)) 1
      (hgcdMatPtr original hOriginal (3 : Fin 4)) 1)
    (hA3 : heap.ValidU64Slice a3 lenInputA)
    (hB3 : heap.ValidU64Slice b3 lenInputB)
    (hAInput : U64SlicesDisjoint a3 lenInputA inputA lenInputA)
    (hBInput : U64SlicesDisjoint b3 lenInputB inputB lenInputB)
    (hAInputB : U64SlicesDisjoint a3 lenInputA inputB lenInputB)
    (hB3A3 : U64SlicesDisjoint b3 lenInputB a3 lenInputA)
    (h0A : U64SlicesDisjoint
      (hgcdMatPtr original hOriginal (0 : Fin 4)) 1 inputA lenInputA)
    (h3A : U64SlicesDisjoint
      (hgcdMatPtr original hOriginal (3 : Fin 4)) 1 inputA lenInputA)
    (h0B : U64SlicesDisjoint
      (hgcdMatPtr original hOriginal (0 : Fin 4)) 1 inputB lenInputB)
    (h3B : U64SlicesDisjoint
      (hgcdMatPtr original hOriginal (3 : Fin 4)) 1 inputB lenInputB)
    (hA3Matrix : ∀ i : Fin 4, U64SlicesDisjoint a3 lenInputA
      (hgcdMatPtr original hOriginal i) (identityEntryLen i))
    (hB3Matrix : ∀ i : Fin 4, U64SlicesDisjoint b3 lenInputB
      (hgcdMatPtr original hOriginal i) (identityEntryLen i))
    (hMatrixValid : ∀ i : Fin 4, heap.ValidU64Slice
      (hgcdMatPtr original hOriginal i) (identityEntryLen i))
    (hLeft : RawDensePolyRep this heap inputA lenInputA left)
    (hRight : RawDensePolyRep this heap inputB lenInputB right)
    (hrun : hgcdRecursiveIterBranch this original hOriginal a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap = .ok result) :
    ∃ finalA finalB finalEntries hResultM,
      HgcdMatRawDenseRep this result.heap result.matrix finalEntries hResultM ∧
      RawDensePolyRep this result.heap a3 result.lenA finalA ∧
      RawDensePolyRep this result.heap b3 result.lenB finalB ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right finalA finalB
        (finalEntries 0) (finalEntries 1) (finalEntries 2) (finalEntries 3) ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet result.sgn
        (finalEntries 0) (finalEntries 1) (finalEntries 2) (finalEntries 3) ∧
      normalize (EuclideanDomain.gcd left right) =
        normalize (EuclideanDomain.gcd finalA finalB) ∧
      result.lenB < lenInputA / 2 + 1 ∧
      hgcdMatLen result.matrix hResultM (0 : Fin 4) + result.lenA ≤
        lenInputA + 1 ∧
      hgcdMatLen result.matrix hResultM (1 : Fin 4) + result.lenB ≤
        lenInputA + 1 ∧
      hgcdMatLen result.matrix hResultM (2 : Fin 4) + result.lenA ≤
        lenInputA + 1 ∧
      hgcdMatLen result.matrix hResultM (3 : Fin 4) + result.lenB ≤
        lenInputA + 1 ∧
      result.lenB ≤ result.lenA ∧ result.lenA ≤ lenInputA ∧
      0 < result.lenA := by
  rcases hgcdRecursiveIterBranch_exec this original hOriginal a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap result hrun with
    ⟨iter, hIter, stable, hiter, hstable, hstore, hResultMatrix,
      hResultLenA, hResultLenB, hResultSgn⟩
  rcases hgcdIter_refines this original a3 b3 T0 T1 0 inputA lenInputA
      inputB lenInputB Q W3 scratch heap iter left right hOriginal hcfg hp
      hInputOrder hInputAPos loopPhysical h0 h3 h03 hA3 hB3 hAInput hBInput hAInputB hB3A3 h0A
      h3A h0B h3B hA3Matrix hB3Matrix hMatrixValid hLeft hRight hiter with
    ⟨finalA, finalB, finalEntries, hFinalIter,
      hInvariant, hGcd, hStop, hMatrixLength, hFinalOrder,
      hFinalInputBound, hFinalPositive⟩
  have hStableDescriptors := hgcdMatStabilize_preserves_descriptors original
    iter.matrix hOriginal hIter stage iter.heap stable hstable
  have hStableValid : stable.matrix.Valid := hStableDescriptors.1
  have hFinalize := finalizePhysical iter hIter stable hStableValid hiter hstable
  rcases hgcdMatStabilize_refines this original iter.matrix hOriginal hIter
      stage finalEntries iter.heap hFinalize.stabilize
      hInvariant.matrixRep with
    ⟨semanticStable, hSemanticStable, hSemanticValid, hStableMatrixRep⟩
  have hStableEq : semanticStable = stable :=
    Except.ok.inj (hSemanticStable.symm.trans hstable)
  subst semanticStable
  have hAFrame := hgcdMatStabilize_preserves_rawDenseRep this original
    iter.matrix hOriginal hIter stage finalEntries iter.heap
    hFinalize.stabilize hInvariant.matrixRep iter.A iter.lenA finalA
    hFinalize.stageA hFinalize.originalA hInvariant.aRep stable hstable
  have hBFrame := hgcdMatStabilize_preserves_rawDenseRep this original
    iter.matrix hOriginal hIter stage finalEntries iter.heap
    hFinalize.stabilize hInvariant.matrixRep iter.B iter.lenB finalB
    hFinalize.stageB hFinalize.originalB hInvariant.bRep stable hstable
  have hA3Stable : stable.heap.ValidU64Slice a3 iter.lenA :=
    (hAFrame.1 a3 iter.lenA).mp hFinalize.validA3
  have hB3Stable : stable.heap.ValidU64Slice b3 iter.lenB :=
    (hBFrame.1 b3 iter.lenB).mp hFinalize.validB3
  rcases hgcdRecursiveStoreIterOutputs_refines this a3 b3 iter.A iter.B
      iter.lenA iter.lenB finalA finalB stable.heap hFinalize.pAEq
      hFinalize.pBEq hA3Stable hB3Stable hAFrame.2 hBFrame.2
      hFinalize.b3PB hFinalize.b3PA hFinalize.a3PA hFinalize.a3PB
      hFinalize.a3B3 with
    ⟨semanticHeap, hSemanticStore, hFinalARep, hFinalBRep⟩
  have hHeapEq : semanticHeap = result.heap :=
    Except.ok.inj (hSemanticStore.symm.trans hstore)
  subst semanticHeap
  have hMatrixAfter := hgcdRecursiveStoreIterOutputs_preserves_matrix this
    a3 b3 iter.A iter.B iter.lenA iter.lenB stable.heap result.heap
    stable.matrix finalEntries hSemanticValid hA3Stable hB3Stable
    hAFrame.2.1 hBFrame.2.1 hFinalize.a3Matrix hFinalize.b3Matrix
    hStableMatrixRep hstore
  have hResultValid : result.matrix.Valid := by
    rw [hResultMatrix]
    exact hSemanticValid
  have hResultLength0 : hgcdMatLen result.matrix hResultValid (0 : Fin 4) +
      result.lenA ≤ lenInputA + 1 := by
    simp only [hgcdMatLen, hResultMatrix, hResultLenA,
      hStableDescriptors.2.2]
    exact hMatrixLength.row0A
  have hResultLength1 : hgcdMatLen result.matrix hResultValid (1 : Fin 4) +
      result.lenB ≤ lenInputA + 1 := by
    simp only [hgcdMatLen, hResultMatrix, hResultLenB,
      hStableDescriptors.2.2]
    exact hMatrixLength.row1B
  have hResultLength2 : hgcdMatLen result.matrix hResultValid (2 : Fin 4) +
      result.lenA ≤ lenInputA + 1 := by
    simp only [hgcdMatLen, hResultMatrix, hResultLenA,
      hStableDescriptors.2.2]
    exact hMatrixLength.row2A
  have hResultLength3 : hgcdMatLen result.matrix hResultValid (3 : Fin 4) +
      result.lenB ≤ lenInputA + 1 := by
    simp only [hgcdMatLen, hResultMatrix, hResultLenB,
      hStableDescriptors.2.2]
    exact hMatrixLength.row3B
  have hMatrixResult : HgcdMatRawDenseRep this result.heap result.matrix
      finalEntries hResultValid := by
    simpa only [hResultMatrix] using hMatrixAfter.2
  have hAResult : RawDensePolyRep this result.heap a3 result.lenA finalA := by
    simpa only [hResultLenA] using hFinalARep
  have hBResult : RawDensePolyRep this result.heap b3 result.lenB finalB := by
    simpa only [hResultLenB] using hFinalBRep
  refine ⟨finalA, finalB, finalEntries, hResultValid, ?_⟩
  refine ⟨hMatrixResult, hAResult, hBResult, ?_, ?_, hGcd, ?_,
    hResultLength0, hResultLength1, hResultLength2, hResultLength3, ?_, ?_,
    ?_⟩
  · simpa only [hResultLenA, hResultLenB] using hInvariant.transform
  · simpa only [hResultSgn] using hInvariant.signedDet
  · simpa only [hResultLenB] using hStop
  · simpa only [hResultLenA, hResultLenB] using hFinalOrder
  · simpa only [hResultLenA] using hFinalInputBound
  · simpa only [hResultLenA] using hFinalPositive

/-- Physical obligations for the two concrete guarded products and the
source-selected tail of one `_mat_mul_entry`. -/
structure HgcdMatMulEntryWorkspace
    (heap : RawHeap) (C P Q R S T scratch : RawPtr UInt64)
    (lenP lenQ lenR lenS : Nat)
    (productPQ productRS : HgcdMulTermResult) : Prop where
  first : HgcdMulTermWorkspace heap C P lenP Q lenQ scratch
  second : HgcdMulTermWorkspace productPQ.heap T R lenR S lenS scratch
  firstDstR : U64SlicesDisjoint C (hgcdMulCapacity lenP lenQ) R lenR
  firstScratchR : U64SlicesDisjoint scratch (8 * max lenP lenQ) R lenR
  firstDstS : U64SlicesDisjoint C (hgcdMulCapacity lenP lenQ) S lenS
  firstScratchS : U64SlicesDisjoint scratch (8 * max lenP lenQ) S lenS
  secondDstC : U64SlicesDisjoint T (hgcdMulCapacity lenR lenS)
    C productPQ.length
  secondScratchC : U64SlicesDisjoint scratch (8 * max lenR lenS)
    C productPQ.length
  finalCValid : productRS.heap.ValidU64Slice C
    (max productPQ.length productRS.length)
  addAliasT : ExactOrDisjoint C T
  copyCT : U64SlicesDisjoint C productRS.length T productRS.length

def HgcdMatMulEntryWorkspaceProvider
    (this : DenseUPolyZp) (heap : RawHeap)
    (C P Q R S T scratch : RawPtr UInt64)
    (lenP lenQ lenR lenS : Nat) : Prop :=
  ∀ productPQ productRS,
    hgcdRecursiveMulTerm this C P lenP Q lenQ scratch heap = .ok productPQ →
    hgcdRecursiveMulTerm this T R lenR S lenS scratch productPQ.heap =
      .ok productRS →
    HgcdMatMulEntryWorkspace heap C P Q R S T scratch lenP lenQ lenR lenS
      productPQ productRS

/-- Complete raw semantic refinement of the actual `_mat_mul_entry` control
flow.  All four source tails produce `P*Q + R*S`; zero products are derived
from their real zero-length raw representations. -/
theorem hgcdMatMulEntry_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (C P Q R S T scratch : RawPtr UInt64)
    (lenP lenQ lenR lenS : Nat) (heap : RawHeap)
    (result : HgcdMatMulEntryResult)
    (polyP polyQ polyR polyS : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatMulEntryWorkspaceProvider this heap C P Q R S T
      scratch lenP lenQ lenR lenS)
    (hP : RawDensePolyRep this heap P lenP polyP)
    (hQ : RawDensePolyRep this heap Q lenQ polyQ)
    (hR : RawDensePolyRep this heap R lenR polyR)
    (hS : RawDensePolyRep this heap S lenS polyS)
    (hrun : hgcdMatMulEntry this C P Q R S T scratch lenP lenQ lenR lenS
      heap = .ok result) :
    RawDensePolyRep this result.heap C result.length
      (polyP * polyQ + polyR * polyS) := by
  rcases hgcdMatMulEntry_exec this C P Q R S T scratch lenP lenQ lenR lenS
      heap result hrun with ⟨productPQ, productRS, hPQ, hRS, htail⟩
  have hwork := physical productPQ productRS hPQ hRS
  rcases hgcdRecursiveMulTerm_refines this C P lenP Q lenQ scratch heap
      polyP polyQ hcfg hp hwork.first hP hQ with
    ⟨semanticPQ, hSemanticPQ, hLayoutPQ, hPQRep⟩
  have hPQEq : semanticPQ = productPQ :=
    Except.ok.inj (hSemanticPQ.symm.trans hPQ)
  subst semanticPQ
  have hRPrefix := hgcdRecursiveMulTerm_preserves_guard this C P lenP Q lenQ
    scratch R lenR heap productPQ hwork.first hP.1 hQ.1 hwork.firstDstR
    hwork.firstScratchR hPQ
  have hSPrefix := hgcdRecursiveMulTerm_preserves_guard this C P lenP Q lenQ
    scratch S lenS heap productPQ hwork.first hP.1 hQ.1 hwork.firstDstS
    hwork.firstScratchS hPQ
  have hR1 := rawDensePolyRep_of_same_prefix this heap productPQ.heap R lenR
    polyR hLayoutPQ hRPrefix hR
  have hS1 := rawDensePolyRep_of_same_prefix this heap productPQ.heap S lenS
    polyS hLayoutPQ hSPrefix hS
  rcases hgcdRecursiveMulTerm_refines this T R lenR S lenS scratch
      productPQ.heap polyR polyS hcfg hp hwork.second hR1 hS1 with
    ⟨semanticRS, hSemanticRS, hLayoutRS, hRSRep⟩
  have hRSEq : semanticRS = productRS :=
    Except.ok.inj (hSemanticRS.symm.trans hRS)
  subst semanticRS
  have hCPrefix := hgcdRecursiveMulTerm_preserves_guard this T R lenR S lenS
    scratch C productPQ.length productPQ.heap productRS hwork.second hR1.1
    hS1.1 hwork.secondDstC hwork.secondScratchC hRS
  have hPQRep2 := rawDensePolyRep_of_same_prefix this productPQ.heap
    productRS.heap C productPQ.length (polyP * polyQ) hLayoutRS hCPrefix hPQRep
  split at htail
  next hboth =>
    rcases htail with ⟨length, hadd, hlength⟩
    subst length
    have hpWord : this._p ≠ 0 := by
      intro hzero
      have hzeroNat := congrArg UInt64.toNat hzero
      simp at hzeroNat
      omega
    exact polyAdd_refines this C C productPQ.length T productRS.length
      productRS.heap result.heap result.length (polyP * polyQ)
      (polyR * polyS) hpWord hwork.finalCValid hPQRep2 hRSRep
      (Or.inl rfl) hwork.addAliasT hadd
  next hnotBoth =>
    split at htail
    next hPQPos =>
      rcases htail with ⟨hheap, hlength⟩
      have hRSLength : productRS.length = 0 := by
        simp at hnotBoth
        omega
      have hzeroRS : polyR * polyS = 0 :=
        slicePolyRep_zero_length productRS.heap T this._p.toNat
          (polyR * polyS) (by simpa [hRSLength] using hRSRep.2.2.1)
      simpa [hheap, hlength, hzeroRS] using hPQRep2
    next hPQZero =>
      have hPQLength : productPQ.length = 0 := by omega
      split at htail
      next hRSPos =>
        rcases htail with ⟨hcopy, hlength⟩
        have hzeroPQ : polyP * polyQ = 0 :=
          slicePolyRep_zero_length productRS.heap C this._p.toNat
            (polyP * polyQ) (by simpa [hPQLength] using hPQRep2.2.2.1)
        have hCValid : productRS.heap.ValidU64Slice C productRS.length :=
          productRS.heap.validU64Slice_mono C
            (max productPQ.length productRS.length) productRS.length
            hwork.finalCValid (by omega)
        rcases copyU64_refines_rawDense this productRS.heap C T
            productRS.length (polyR * polyS) hCValid hwork.copyCT hRSRep with
          ⟨heap', hcopy', _, hrep⟩
        have heq : heap' = result.heap :=
          Except.ok.inj (hcopy'.symm.trans hcopy)
        subst heap'
        simpa [hlength, hzeroPQ] using hrep
      next hRSZero =>
        rcases htail with ⟨hheap, hlength⟩
        have hRSLength : productRS.length = 0 := by omega
        have hzeroPQ : polyP * polyQ = 0 :=
          slicePolyRep_zero_length productRS.heap C this._p.toNat
            (polyP * polyQ) (by simpa [hPQLength] using hPQRep2.2.2.1)
        have hzeroRS : polyR * polyS = 0 :=
          slicePolyRep_zero_length productRS.heap T this._p.toNat
            (polyR * polyS) (by simpa [hRSLength] using hRSRep.2.2.1)
        simpa [hheap, hlength, hzeroPQ, hzeroRS] using
          rawDensePolyRep_zero_length this productRS.heap C
            (productRS.heap.validU64Slice_mono C
              (max productPQ.length productRS.length) 0 hwork.finalCValid
              (by omega))

/-- Purely spatial obligations saying that one concrete `_mat_mul_entry`
does not overwrite an additional live polynomial.  In particular, this
contract contains no polynomial value and no expected matrix product. -/
structure HgcdMatMulEntryGuardWorkspace
    (C P Q R S T scratch guard : RawPtr UInt64)
    (lenP lenQ lenR lenS guardLen : Nat) : Prop where
  firstDst : U64SlicesDisjoint C (hgcdMulCapacity lenP lenQ) guard guardLen
  firstScratch : U64SlicesDisjoint scratch (8 * max lenP lenQ) guard guardLen
  secondDst : U64SlicesDisjoint T (hgcdMulCapacity lenR lenS) guard guardLen
  secondScratch : U64SlicesDisjoint scratch (8 * max lenR lenS) guard guardLen
  addDstRegion : C.region ≠ guard.region

/-- Frame rule for the complete, source-selected `_mat_mul_entry` execution.
The proof follows both guarded products and then the actual add/copy/skip
tail.  It is the bridge needed to retain both input matrices and previously
computed output entries while `_mat_mul` advances. -/
theorem hgcdMatMulEntry_preserves_rawDenseRep (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (C P Q R S T scratch guard : RawPtr UInt64)
    (lenP lenQ lenR lenS guardLen : Nat) (heap : RawHeap)
    (result : HgcdMatMulEntryResult)
    (polyP polyQ polyR polyS guardPoly : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatMulEntryWorkspaceProvider this heap C P Q R S T
      scratch lenP lenQ lenR lenS)
    (frame : HgcdMatMulEntryGuardWorkspace C P Q R S T scratch guard
      lenP lenQ lenR lenS guardLen)
    (hP : RawDensePolyRep this heap P lenP polyP)
    (hQ : RawDensePolyRep this heap Q lenQ polyQ)
    (hR : RawDensePolyRep this heap R lenR polyR)
    (hS : RawDensePolyRep this heap S lenS polyS)
    (hGuard : RawDensePolyRep this heap guard guardLen guardPoly)
    (hrun : hgcdMatMulEntry this C P Q R S T scratch lenP lenQ lenR lenS
      heap = .ok result) :
    RawHeap.SameLayout heap result.heap ∧
      RawDensePolyRep this result.heap guard guardLen guardPoly := by
  rcases hgcdMatMulEntry_exec this C P Q R S T scratch lenP lenQ lenR lenS
      heap result hrun with ⟨productPQ, productRS, hPQ, hRS, htail⟩
  have hwork := physical productPQ productRS hPQ hRS
  rcases hgcdRecursiveMulTerm_refines this C P lenP Q lenQ scratch heap
      polyP polyQ hcfg hp hwork.first hP hQ with
    ⟨semanticPQ, hSemanticPQ, hLayoutPQ, hPQRep⟩
  have hPQEq : semanticPQ = productPQ :=
    Except.ok.inj (hSemanticPQ.symm.trans hPQ)
  subst semanticPQ
  have hGuardPrefix1 := hgcdRecursiveMulTerm_preserves_guard this C P lenP
    Q lenQ scratch guard guardLen heap productPQ hwork.first hP.1 hQ.1
    frame.firstDst frame.firstScratch hPQ
  have hGuard1 := rawDensePolyRep_of_same_prefix this heap productPQ.heap
    guard guardLen guardPoly hLayoutPQ hGuardPrefix1 hGuard
  have hRPrefix := hgcdRecursiveMulTerm_preserves_guard this C P lenP Q lenQ
    scratch R lenR heap productPQ hwork.first hP.1 hQ.1 hwork.firstDstR
    hwork.firstScratchR hPQ
  have hSPrefix := hgcdRecursiveMulTerm_preserves_guard this C P lenP Q lenQ
    scratch S lenS heap productPQ hwork.first hP.1 hQ.1 hwork.firstDstS
    hwork.firstScratchS hPQ
  have hR1 := rawDensePolyRep_of_same_prefix this heap productPQ.heap R lenR
    polyR hLayoutPQ hRPrefix hR
  have hS1 := rawDensePolyRep_of_same_prefix this heap productPQ.heap S lenS
    polyS hLayoutPQ hSPrefix hS
  rcases hgcdRecursiveMulTerm_refines this T R lenR S lenS scratch
      productPQ.heap polyR polyS hcfg hp hwork.second hR1 hS1 with
    ⟨semanticRS, hSemanticRS, hLayoutRS, hRSRep⟩
  have hRSEq : semanticRS = productRS :=
    Except.ok.inj (hSemanticRS.symm.trans hRS)
  subst semanticRS
  have hGuardPrefix2 := hgcdRecursiveMulTerm_preserves_guard this T R lenR
    S lenS scratch guard guardLen productPQ.heap productRS hwork.second
    hR1.1 hS1.1 frame.secondDst frame.secondScratch hRS
  have hGuard2 := rawDensePolyRep_of_same_prefix this productPQ.heap
    productRS.heap guard guardLen guardPoly hLayoutRS hGuardPrefix2 hGuard1
  have hCPrefix := hgcdRecursiveMulTerm_preserves_guard this T R lenR S lenS
    scratch C productPQ.length productPQ.heap productRS hwork.second hR1.1
    hS1.1 hwork.secondDstC hwork.secondScratchC hRS
  have hPQRep2 := rawDensePolyRep_of_same_prefix this productPQ.heap
    productRS.heap C productPQ.length (polyP * polyQ) hLayoutRS hCPrefix
    hPQRep
  have hLayout12 : RawHeap.SameLayout heap productRS.heap :=
    fun ptr count => (hLayoutPQ ptr count).trans (hLayoutRS ptr count)
  split at htail
  next hboth =>
    rcases htail with ⟨length, hadd, hlength⟩
    subst length
    rcases polyAdd_ok this C C productPQ.length T productRS.length
        productRS.heap hwork.finalCValid hPQRep2.1 hRSRep.1 with
      ⟨heap', length, hadd', hLayoutTail, _⟩
    have heq : (heap', length) = (result.heap, result.length) :=
      Except.ok.inj (hadd'.symm.trans hadd)
    have hHeapEq : heap' = result.heap := congrArg Prod.fst heq
    subst heap'
    have hPrefixTail := polyAdd_preserves_prefix_region_ne this C C
      productPQ.length T guard productRS.length guardLen productRS.heap
      result.heap result.length hwork.finalCValid hPQRep2.1 hRSRep.1
      frame.addDstRegion hadd
    exact ⟨fun ptr count =>
        (hLayout12 ptr count).trans (hLayoutTail ptr count),
      rawDensePolyRep_of_same_prefix this productRS.heap result.heap guard
        guardLen guardPoly hLayoutTail hPrefixTail hGuard2⟩
  next hnotBoth =>
    split at htail
    next hPQPos =>
      rcases htail with ⟨hheap, _⟩
      simpa [hheap] using And.intro hLayout12 hGuard2
    next hPQZero =>
      split at htail
      next hRSPos =>
        rcases htail with ⟨hcopy, _⟩
        have hCValid : productRS.heap.ValidU64Slice C productRS.length :=
          productRS.heap.validU64Slice_mono C
            (max productPQ.length productRS.length) productRS.length
            hwork.finalCValid (by omega)
        have hCopyFrame := copyU64_preserves_rawDenseRep this productRS.heap
          result.heap C T productRS.length guard guardLen guardPoly hCValid
          hRSRep.1
          (by intro _ _ _ _; exact Or.inl frame.addDstRegion) hcopy hGuard2
        exact ⟨fun ptr count =>
            (hLayout12 ptr count).trans (hCopyFrame.1 ptr count),
          hCopyFrame.2⟩
      next hRSZero =>
        rcases htail with ⟨hheap, _⟩
        simpa [hheap] using And.intro hLayout12 hGuard2

/-- The L2 entry selected by source iteration `i` of `_mat_mul`. -/
noncomputable def hgcdMatProductEntry {p : Nat}
    (left right : Fin 4 → Polynomial (ZMod p)) (i : Fin 4) :
    Polynomial (ZMod p) :=
  left ⟨2 * (i.val / 2), by omega⟩ * right ⟨i.val % 2, by omega⟩ +
    left ⟨2 * (i.val / 2) + 1, by omega⟩ *
      right ⟨2 + i.val % 2, by omega⟩

/-- Physical obligations for one actual source iteration of `_mat_mul`.
The frame fields mention only addresses and lengths; the expected L2 matrix
product is deliberately absent. -/
structure HgcdMatMulStepWorkspace (this : DenseUPolyZp)
    (A B C : HgcdMat) (hA : A.Valid) (hB : B.Valid) (hC : C.Valid)
    (T scratch : RawPtr UInt64) (i : Nat) (hi : i < 4)
    (heap : RawHeap) : Prop where
  entry : HgcdMatMulEntryWorkspaceProvider this heap
    (hgcdMatPtr C hC ⟨i, hi⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatPtr B hB ⟨i % 2, by omega⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatPtr B hB ⟨2 + i % 2, by omega⟩) T scratch
    (hgcdMatLen A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatLen B hB ⟨i % 2, by omega⟩)
    (hgcdMatLen A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatLen B hB ⟨2 + i % 2, by omega⟩)
  frameA : ∀ j : Fin 4, HgcdMatMulEntryGuardWorkspace
    (hgcdMatPtr C hC ⟨i, hi⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatPtr B hB ⟨i % 2, by omega⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatPtr B hB ⟨2 + i % 2, by omega⟩) T scratch
    (hgcdMatPtr A hA j)
    (hgcdMatLen A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatLen B hB ⟨i % 2, by omega⟩)
    (hgcdMatLen A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatLen B hB ⟨2 + i % 2, by omega⟩)
    (hgcdMatLen A hA j)
  frameB : ∀ j : Fin 4, HgcdMatMulEntryGuardWorkspace
    (hgcdMatPtr C hC ⟨i, hi⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatPtr B hB ⟨i % 2, by omega⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatPtr B hB ⟨2 + i % 2, by omega⟩) T scratch
    (hgcdMatPtr B hB j)
    (hgcdMatLen A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatLen B hB ⟨i % 2, by omega⟩)
    (hgcdMatLen A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatLen B hB ⟨2 + i % 2, by omega⟩)
    (hgcdMatLen B hB j)
  frameDone : ∀ j : Fin 4, j.val < i → HgcdMatMulEntryGuardWorkspace
    (hgcdMatPtr C hC ⟨i, hi⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatPtr B hB ⟨i % 2, by omega⟩)
    (hgcdMatPtr A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatPtr B hB ⟨2 + i % 2, by omega⟩) T scratch
    (hgcdMatPtr C hC j)
    (hgcdMatLen A hA ⟨2 * (i / 2), by omega⟩)
    (hgcdMatLen B hB ⟨i % 2, by omega⟩)
    (hgcdMatLen A hA ⟨2 * (i / 2) + 1, by omega⟩)
    (hgcdMatLen B hB ⟨2 + i % 2, by omega⟩)
    (hgcdMatLen C hC j)

/-- A provider for every concrete state reached by the four-entry loop. -/
def HgcdMatMulLoopWorkspaceProvider (this : DenseUPolyZp)
    (A B : HgcdMat) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) : Prop :=
  ∀ (C : HgcdMat) (hC : C.Valid) (i : Nat) (hi : i < 4)
    (heap : RawHeap),
    HgcdMatMulStepWorkspace this A B C hA hB hC T scratch i hi heap

/-- Semantic refinement of the actual four-step `_mat_mul` loop.  The
induction follows `hgcdMatMulLoop`; each new result comes from executing
`hgcdMatMulEntry`, while frame proofs retain both inputs and earlier outputs. -/
theorem hgcdMatMulLoop_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (A B : HgcdMat) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) (C : HgcdMat) (hC : C.Valid)
    (i : Nat) (heap : RawHeap) (result : HgcdMatMulResult)
    (left right : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatMulLoopWorkspaceProvider this A B hA hB T scratch)
    (hLeft : HgcdMatRawDenseRep this heap A left hA)
    (hRight : HgcdMatRawDenseRep this heap B right hB)
    (hDone : ∀ j : Fin 4, j.val < i → RawDensePolyRep this heap
      (hgcdMatPtr C hC j) (hgcdMatLen C hC j)
      (hgcdMatProductEntry left right j))
    (hrun : hgcdMatMulLoop this A B hA hB T scratch C hC i heap =
      .ok result) :
    HgcdMatRawDenseRep this result.heap A left hA ∧
      HgcdMatRawDenseRep this result.heap B right hB ∧
      ∃ hResult : result.matrix.Valid,
        HgcdMatRawDenseRep this result.heap result.matrix
          (hgcdMatProductEntry left right) hResult := by
  rw [hgcdMatMulLoop] at hrun
  split at hrun
  next hi =>
    dsimp only at hrun
    split at hrun
    next fault hentry => simp at hrun
    next entry hentry =>
      let out : Fin 4 := ⟨i, hi⟩
      let rowBase : Fin 4 := ⟨2 * (i / 2), by omega⟩
      let rowNext : Fin 4 := ⟨2 * (i / 2) + 1, by omega⟩
      let col : Fin 4 := ⟨i % 2, by omega⟩
      let lowerCol : Fin 4 := ⟨2 + i % 2, by omega⟩
      have hstep := physical C hC i hi heap
      have hCurrent := hgcdMatMulEntry_refines this
        (hgcdMatPtr C hC out) (hgcdMatPtr A hA rowBase)
        (hgcdMatPtr B hB col) (hgcdMatPtr A hA rowNext)
        (hgcdMatPtr B hB lowerCol) T scratch
        (hgcdMatLen A hA rowBase) (hgcdMatLen B hB col)
        (hgcdMatLen A hA rowNext) (hgcdMatLen B hB lowerCol) heap entry
        (left rowBase) (right col) (left rowNext) (right lowerCol)
        hcfg hp hstep.entry (hLeft rowBase) (hRight col) (hLeft rowNext)
        (hRight lowerCol) hentry
      have hLeft1 : HgcdMatRawDenseRep this entry.heap A left hA := by
        intro j
        exact (hgcdMatMulEntry_preserves_rawDenseRep this
          (hgcdMatPtr C hC out) (hgcdMatPtr A hA rowBase)
          (hgcdMatPtr B hB col) (hgcdMatPtr A hA rowNext)
          (hgcdMatPtr B hB lowerCol) T scratch (hgcdMatPtr A hA j)
          (hgcdMatLen A hA rowBase) (hgcdMatLen B hB col)
          (hgcdMatLen A hA rowNext) (hgcdMatLen B hB lowerCol)
          (hgcdMatLen A hA j) heap entry (left rowBase) (right col)
          (left rowNext) (right lowerCol) (left j) hcfg hp hstep.entry
          (hstep.frameA j) (hLeft rowBase) (hRight col) (hLeft rowNext)
          (hRight lowerCol) (hLeft j) hentry).2
      have hRight1 : HgcdMatRawDenseRep this entry.heap B right hB := by
        intro j
        exact (hgcdMatMulEntry_preserves_rawDenseRep this
          (hgcdMatPtr C hC out) (hgcdMatPtr A hA rowBase)
          (hgcdMatPtr B hB col) (hgcdMatPtr A hA rowNext)
          (hgcdMatPtr B hB lowerCol) T scratch (hgcdMatPtr B hB j)
          (hgcdMatLen A hA rowBase) (hgcdMatLen B hB col)
          (hgcdMatLen A hA rowNext) (hgcdMatLen B hB lowerCol)
          (hgcdMatLen B hB j) heap entry (left rowBase) (right col)
          (left rowNext) (right lowerCol) (right j) hcfg hp hstep.entry
          (hstep.frameB j) (hLeft rowBase) (hRight col) (hLeft rowNext)
          (hRight lowerCol) (hRight j) hentry).2
      let nextLen := C.len.set i entry.length
        (by rw [hC.2]; exact hi)
      let next : HgcdMat := { C with len := nextLen }
      have hNext : next.Valid := by
        exact ⟨hC.1, by simp [next, nextLen, hC.2]⟩
      have hDone1 : ∀ j : Fin 4, j.val < i + 1 →
          RawDensePolyRep this entry.heap (hgcdMatPtr next hNext j)
            (hgcdMatLen next hNext j) (hgcdMatProductEntry left right j) := by
        intro j hj
        by_cases hji : j = out
        · subst j
          simpa [out, rowBase, rowNext, col, lowerCol, next, nextLen,
            hgcdMatProductEntry, hgcdMatPtr, hgcdMatLen] using hCurrent
        · have hjne : j.val ≠ i := by
            intro hval
            exact hji (Fin.ext hval)
          have hjlt : j.val < i := by omega
          have hpreserved := (hgcdMatMulEntry_preserves_rawDenseRep this
            (hgcdMatPtr C hC out) (hgcdMatPtr A hA rowBase)
            (hgcdMatPtr B hB col) (hgcdMatPtr A hA rowNext)
            (hgcdMatPtr B hB lowerCol) T scratch (hgcdMatPtr C hC j)
            (hgcdMatLen A hA rowBase) (hgcdMatLen B hB col)
            (hgcdMatLen A hA rowNext) (hgcdMatLen B hB lowerCol)
            (hgcdMatLen C hC j) heap entry (left rowBase) (right col)
            (left rowNext) (right lowerCol) (hgcdMatProductEntry left right j)
            hcfg hp hstep.entry (hstep.frameDone j hjlt) (hLeft rowBase)
            (hRight col) (hLeft rowNext) (hRight lowerCol) (hDone j hjlt)
            hentry).2
          simp only [next, nextLen, hgcdMatPtr, hgcdMatLen]
          rw [Array.getElem_set_ne
            (by simpa [hC.2] using hi)
            (by simpa [hC.2] using j.isLt) (Ne.symm hjne)]
          exact hpreserved
      exact hgcdMatMulLoop_refines this A B hA hB T scratch next hNext
        (i + 1) entry.heap result left right hcfg hp physical hLeft1 hRight1
        hDone1 hrun
  next hi =>
    have heq : result = HgcdMatMulResult.mk heap C :=
      (Except.ok.inj hrun).symm
    subst result
    exact ⟨hLeft, hRight, hC, fun j => hDone j (by omega)⟩
termination_by 4 - i
decreasing_by omega

/-- Entry-point theorem for the complete generated `_mat_mul`. -/
theorem hgcdMatMul_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (C A B : HgcdMat) (hC : C.Valid) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) (heap : RawHeap) (result : HgcdMatMulResult)
    (left right : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatMulLoopWorkspaceProvider this A B hA hB T scratch)
    (hLeft : HgcdMatRawDenseRep this heap A left hA)
    (hRight : HgcdMatRawDenseRep this heap B right hB)
    (hrun : hgcdMatMul this C A B hC hA hB T scratch heap = .ok result) :
    HgcdMatRawDenseRep this result.heap A left hA ∧
      HgcdMatRawDenseRep this result.heap B right hB ∧
      ∃ hResult : result.matrix.Valid,
        HgcdMatRawDenseRep this result.heap result.matrix
          (hgcdMatProductEntry left right) hResult := by
  exact hgcdMatMulLoop_refines this A B hA hB T scratch C hC 0 heap result
    left right hcfg hp physical hLeft hRight (by intro j hj; omega) hrun

/-- L2 effect of one quotient-matrix column update. -/
noncomputable def hgcdMatQuotientUpdateEntries {p : Nat}
    (entries : Fin 4 → Polynomial (ZMod p))
    (quotient : Polynomial (ZMod p)) (top bottom : Fin 4) :
    Fin 4 → Polynomial (ZMod p) :=
  Function.update entries top (entries top + quotient * entries bottom)

/-- The quotient-entry lowering has exactly the same physical multiplication
and addition requirements as a row update with `bottom` as the multiplied
entry and the existing `top` buffer as the in-place sum destination. -/
abbrev HgcdMatQuotientEntryWorkspace (S : HgcdMat) (hS : S.Valid)
    (top bottom : Fin 4) (q : RawPtr UInt64) (lenQ : Nat)
    (T scratch : RawPtr UInt64) (heap : RawHeap) : Prop :=
  MatRowUpdateWorkspace S bottom top q lenQ T (hgcdMatPtr S hS top)
    scratch heap hS

/-- Complete raw refinement of one actual guarded quotient-column update.
The inactive branch derives the vanishing product from a real zero-length
representation.  The active branch consumes the generated ordered `_mul`
and in-place `_poly_add`, then frames every non-target matrix entry. -/
theorem hgcdMatQuotientEntry_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (S : HgcdMat) (hS : S.Valid) (top bottom : Fin 4)
    (q : RawPtr UInt64) (lenQ : Nat) (T scratch : RawPtr UInt64)
    (heap : RawHeap) (result : HgcdMatQuotientEntryResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (workspace : HgcdMatQuotientEntryWorkspace S hS top bottom q lenQ T
      scratch heap)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hMatrix : HgcdMatRawDenseRep this heap S entries hS)
    (hrun : hgcdMatQuotientEntry this S hS top bottom q lenQ T scratch
      heap = .ok result) :
    HgcdMatRawDenseRep this result.heap result.matrix
      (hgcdMatQuotientUpdateEntries entries quotient top bottom)
      result.valid ∧
    RawDensePolyRep this result.heap q lenQ quotient := by
  simp only [hgcdMatQuotientEntry] at hrun
  split at hrun
  next hactive =>
    have hparts : decide (lenQ > 0) = true ∧
        decide (hgcdMatLenRaw S hS bottom > 0) = true := by
      simpa using hactive
    have hQPos : 0 < lenQ := by
      simpa only [decide_eq_true_eq] using hparts.1
    have hBottomPos : 0 < hgcdMatLen S hS bottom := by
      simpa only [decide_eq_true_eq, hgcdMatLenRaw, hgcdMatLen] using
        hparts.2
    generalize hmul : Generated.StrictMul.dense_upoly_zp__mul_ir this T
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then q
          else hgcdMatPtrRaw S hS bottom)
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then lenQ
          else hgcdMatLenRaw S hS bottom)
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then
          hgcdMatPtrRaw S hS bottom else q)
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then
          hgcdMatLenRaw S hS bottom else lenQ)
        scratch heap = mulRun at hrun
    cases mulRun with
    | error fault => simp at hrun
    | ok heap1 =>
      simp only [hmul] at hrun
      generalize hadd : Generated.StrictPolyAddSub.dense_upoly_zp__poly_add_ir this
          (hgcdMatPtrRaw S hS top) (hgcdMatPtrRaw S hS top)
          (hgcdMatLenRaw S hS top) T
          (lenQ + hgcdMatLenRaw S hS bottom - 1) heap1 = addRun at hrun
      cases addRun with
      | error fault => simp [hadd] at hrun
      | ok pair =>
        rcases pair with ⟨heap2, sumLen⟩
        simp only [hadd] at hrun
        have heq := Except.ok.inj hrun
        cases heq
        have hProduct := matRowUpdate_mul_result this S hS bottom q T scratch
          lenQ heap heap1 quotient (entries bottom) hcfg hp hQPos hBottomPos
          workspace.lenWord workspace.validT workspace.validScratch
          workspace.disjointTQ (workspace.disjointTMatrix bottom)
          workspace.disjointTScratch workspace.disjointScratchQ
          (workspace.disjointScratchMatrix bottom) hQ (hMatrix bottom)
          (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
              hgcdMatLen] using hmul)
        have hMatrix1 : HgcdMatRawDenseRep this heap1 S entries hS := by
          intro j
          exact matRowUpdate_mul_preserves_guard this S hS bottom q T scratch
            (hgcdMatPtr S hS j) lenQ (hgcdMatLen S hS j) heap heap1
            (entries j) hQPos hBottomPos workspace.validT
            workspace.validScratch hQ.1 (hMatrix bottom).1
            workspace.disjointScratchQ
            (workspace.disjointScratchMatrix bottom)
            (workspace.disjointTMatrix j)
            (workspace.disjointScratchMatrix j) (hMatrix j) hProduct.1
            (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
                hgcdMatLen] using hmul)
        have hQ1 := matRowUpdate_mul_preserves_guard this S hS bottom q T
          scratch q lenQ lenQ heap heap1 quotient hQPos hBottomPos
          workspace.validT workspace.validScratch hQ.1 (hMatrix bottom).1
          workspace.disjointScratchQ
          (workspace.disjointScratchMatrix bottom) workspace.disjointTQ
          workspace.disjointScratchQ hQ hProduct.1
          (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
              hgcdMatLen] using hmul)
        have hOutput1 := (hProduct.1 _ _).mp workspace.validAddOutput
        have hpWord : this._p ≠ 0 := by
          intro hzero
          have hzeroNat := congrArg UInt64.toNat hzero
          simp at hzeroNat
          omega
        have hSum : RawDensePolyRep this heap2 (hgcdMatPtr S hS top)
            sumLen (entries top + quotient * entries bottom) :=
          polyAdd_refines this (hgcdMatPtr S hS top)
            (hgcdMatPtr S hS top) (hgcdMatLen S hS top) T
            (lenQ + hgcdMatLen S hS bottom - 1) heap1 heap2 sumLen
            (entries top) (quotient * entries bottom) hpWord hOutput1
            (hMatrix1 top) hProduct.2 workspace.aliasEntry1
            workspace.aliasProduct
            (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
                hgcdMatLen] using hadd)
        rcases polyAdd_ok this (hgcdMatPtr S hS top)
            (hgcdMatPtr S hS top) (hgcdMatLen S hS top) T
            (lenQ + hgcdMatLen S hS bottom - 1) heap1 hOutput1
            (hMatrix1 top).1 hProduct.2.1 with
          ⟨heapAdd, lenAdd, haddAdd, hLayoutAdd, _⟩
        have heqAdd : (heapAdd, lenAdd) = (heap2, sumLen) :=
          Except.ok.inj (haddAdd.symm.trans (by
            simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
              hgcdMatLen] using hadd))
        have hHeapAdd : heapAdd = heap2 := congrArg Prod.fst heqAdd
        have hLayoutAdd2 : RawHeap.SameLayout heap1 heap2 := by
          intro ptr count
          rw [← hHeapAdd]
          exact hLayoutAdd ptr count
        have hQPrefix := polyAdd_preserves_prefix_region_ne this
          (hgcdMatPtr S hS top) (hgcdMatPtr S hS top)
          (hgcdMatLen S hS top) T q
          (lenQ + hgcdMatLen S hS bottom - 1) lenQ heap1 heap2 sumLen
          hOutput1 (hMatrix1 top).1 hProduct.2.1 workspace.tDisjointQ
          (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
              hgcdMatLen] using hadd)
        have hQ2 := rawDensePolyRep_of_same_prefix this heap1 heap2 q lenQ
          quotient hLayoutAdd2 hQPrefix hQ1
        constructor
        · intro j
          by_cases hj : j = top
          · subst j
            simpa [hgcdMatQuotientUpdateEntries, hgcdMatPtr, hgcdMatLen,
              hS.2] using hSum
          · have hPrefix := polyAdd_preserves_prefix_region_ne this
              (hgcdMatPtr S hS top) (hgcdMatPtr S hS top)
              (hgcdMatLen S hS top) T (hgcdMatPtr S hS j)
              (lenQ + hgcdMatLen S hS bottom - 1) (hgcdMatLen S hS j)
              heap1 heap2 sumLen hOutput1 (hMatrix1 top).1 hProduct.2.1
              (workspace.tDisjointOther j hj)
              (by simpa [hgcdMatPtrRaw, hgcdMatLenRaw, hgcdMatPtr,
                  hgcdMatLen] using hadd)
            have hFrame := rawDensePolyRep_of_same_prefix this heap1 heap2
              (hgcdMatPtr S hS j) (hgcdMatLen S hS j) (entries j)
              hLayoutAdd2 hPrefix (hMatrix1 j)
            simp only [hgcdMatPtr, hgcdMatLen]
            rw [Array.getElem_set_ne
              (by simpa [hS.2] using top.isLt)
              (by simpa [hS.2] using j.isLt) (by
                intro hval
                exact hj (Fin.ext hval.symm))]
            simpa [hgcdMatQuotientUpdateEntries, hj] using hFrame
        · simpa using hQ2
  next hinactive =>
    have heq := Except.ok.inj hrun
    subst result
    have hzero : quotient = 0 ∨ entries bottom = 0 := by
      have hnot : ¬(0 < lenQ ∧ 0 < hgcdMatLen S hS bottom) := by
        simpa only [Bool.and_eq_true, decide_eq_true_eq, hgcdMatLenRaw,
          hgcdMatLen] using hinactive
      rcases not_and_or.mp hnot with hq | hb
      · left
        have hlen : lenQ = 0 := by omega
        exact slicePolyRep_zero_length heap q this._p.toNat quotient
          (by simpa [hlen] using hQ.2.2.1)
      · right
        have hlen : hgcdMatLen S hS bottom = 0 := by omega
        exact slicePolyRep_zero_length heap (hgcdMatPtr S hS bottom)
          this._p.toNat (entries bottom)
          (by simpa [hlen] using (hMatrix bottom).2.2.1)
    constructor
    · intro j
      by_cases hj : j = top
      · subst j
        rcases hzero with hq | hb
        · simpa [hgcdMatQuotientUpdateEntries, hq] using hMatrix top
        · simpa [hgcdMatQuotientUpdateEntries, hb] using hMatrix top
      · simpa [hgcdMatQuotientUpdateEntries, hj] using hMatrix j
    · exact hQ

/-- L2 descriptor permutation induced by the four source row swaps. -/
noncomputable def hgcdMatSwapEntries {p : Nat}
    (entries : Fin 4 → Polynomial (ZMod p)) :
    Fin 4 → Polynomial (ZMod p) :=
  fun j => entries ⟨(j.val + 2) % 4, by omega⟩

/-- The descriptor-only row swaps preserve the same four raw buffers and
only permute their L2 indexing. -/
theorem hgcdMatSwapRows_refines (this : DenseUPolyZp)
    (S : HgcdMat) (hS : S.Valid)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat)) (heap : RawHeap)
    (hMatrix : HgcdMatRawDenseRep this heap S entries hS) :
    HgcdMatRawDenseRep this heap (hgcdMatSwapRows S hS)
      (hgcdMatSwapEntries entries) (hgcdMatSwapRows_valid S hS) := by
  intro j
  fin_cases j <;>
    simpa [hgcdMatSwapRows, hgcdMatSwapEntries, hgcdMatPtr, hgcdMatLen,
      hgcdMatPtrRaw, hgcdMatLenRaw] using
      hMatrix (⟨_, by omega⟩ : Fin 4)

/-- Final L2 entries after both source-order quotient column updates. -/
noncomputable def hgcdMatApplyQuotientEntries {p : Nat}
    (entries : Fin 4 → Polynomial (ZMod p))
    (quotient : Polynomial (ZMod p)) : Fin 4 → Polynomial (ZMod p) :=
  let swapped := hgcdMatSwapEntries entries
  let first := hgcdMatQuotientUpdateEntries swapped quotient
    (0 : Fin 4) (2 : Fin 4)
  hgcdMatQuotientUpdateEntries first quotient (1 : Fin 4) (3 : Fin 4)

/-- Physical provider for the two concrete column executions.  It is
quantified over the actual first result and contains no L2 polynomial. -/
def HgcdMatApplyQuotientWorkspaceProvider (this : DenseUPolyZp)
    (S : HgcdMat) (hS : S.Valid) (q : RawPtr UInt64) (lenQ : Nat)
    (T scratch : RawPtr UInt64) (heap : RawHeap) : Prop :=
  ∀ first,
    hgcdMatQuotientEntry this (hgcdMatSwapRows S hS)
      (hgcdMatSwapRows_valid S hS) (0 : Fin 4) (2 : Fin 4) q lenQ T
      scratch heap = .ok first →
    HgcdMatQuotientEntryWorkspace (hgcdMatSwapRows S hS)
        (hgcdMatSwapRows_valid S hS) (0 : Fin 4) (2 : Fin 4) q lenQ T
        scratch heap ∧
      HgcdMatQuotientEntryWorkspace first.matrix first.valid
        (1 : Fin 4) (3 : Fin 4) q lenQ T scratch first.heap

/-- End-to-end refinement of the exact source block
`S := [[q,1],[1,0]] * S`: descriptor swaps and both real guarded
multiplication/addition executions are all consumed in source order. -/
theorem hgcdMatApplyQuotient_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (S : HgcdMat) (hS : S.Valid) (q : RawPtr UInt64) (lenQ : Nat)
    (T scratch : RawPtr UInt64) (heap : RawHeap)
    (result : HgcdMatQuotientResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatApplyQuotientWorkspaceProvider this S hS q lenQ T
      scratch heap)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hMatrix : HgcdMatRawDenseRep this heap S entries hS)
    (hrun : hgcdMatApplyQuotient this S hS q lenQ T scratch heap =
      .ok result) :
    HgcdMatRawDenseRep this result.heap result.matrix
      (hgcdMatApplyQuotientEntries entries quotient) result.valid ∧
    RawDensePolyRep this result.heap q lenQ quotient := by
  rcases hgcdMatApplyQuotient_exec this S hS q lenQ T scratch heap result
      hrun with ⟨first, hfirst, hsecond⟩
  have hwork := physical first hfirst
  have hSwapped := hgcdMatSwapRows_refines this S hS entries heap hMatrix
  have hFirst := hgcdMatQuotientEntry_refines this (hgcdMatSwapRows S hS)
    (hgcdMatSwapRows_valid S hS) (0 : Fin 4) (2 : Fin 4) q lenQ T
    scratch heap first (hgcdMatSwapEntries entries) quotient hcfg hp hwork.1
    hQ hSwapped hfirst
  have hSecond := hgcdMatQuotientEntry_refines this first.matrix first.valid
    (1 : Fin 4) (3 : Fin 4) q lenQ T scratch first.heap
    (HgcdMatQuotientEntryResult.mk result.heap result.matrix result.valid)
    (hgcdMatQuotientUpdateEntries (hgcdMatSwapEntries entries) quotient
      (0 : Fin 4) (2 : Fin 4)) quotient hcfg hp hwork.2 hFirst.2 hFirst.1
    hsecond
  simpa [hgcdMatApplyQuotientEntries] using hSecond

/-- L2 low parts computed by the two sign-selected reconstruction blocks. -/
noncomputable def hgcdReconstructedLowB {p : Nat}
    (entries : Fin 4 → Polynomial (ZMod p))
    (lowA lowB : Polynomial (ZMod p)) (sgn : Int) :=
  if sgn < 0 then entries 2 * lowA - entries 0 * lowB
  else entries 0 * lowB - entries 2 * lowA

noncomputable def hgcdReconstructedLowA {p : Nat}
    (entries : Fin 4 → Polynomial (ZMod p))
    (lowA lowB : Polynomial (ZMod p)) (sgn : Int) :=
  if sgn < 0 then entries 1 * lowB - entries 3 * lowA
  else entries 3 * lowA - entries 1 * lowB

/-- Physical obligations for the exact four-call reconstruction pair.  All
frame fields mention only heap cells; no expected L2 polynomial occurs in
this contract. -/
structure HgcdRecursiveReconstructPairWorkspace (this : DenseUPolyZp)
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift : Nat)
    (M : HgcdMat) (hM : M.Valid) (sgn : Int) (heap : RawHeap)
    (heap1 : RawHeap) (lowLenB : Nat) (liftedB : HgcdLiftHighResult)
    (heap3 : RawHeap) (lowLenA : Nat) (liftedA : HgcdLiftHighResult) : Prop where
  reconstructB : HgcdReconstructWorkspace heap B T0
    (hgcdMatPtr M hM (2 : Fin 4)) (hgcdMatLen M hM (2 : Fin 4)) lowA lenLowA
    (hgcdMatPtr M hM (0 : Fin 4)) (hgcdMatLen M hM (0 : Fin 4)) lowB lenLowB
    scratch
  highBLayout : RawHeap.SameLayout heap heap1
  highBPrefix : SameU64Prefix heap heap1 highB lenHighB
  liftB : HgcdLiftHighWorkspace heap1 B highB lowLenB shift lenHighB
  aInputLayout : RawHeap.SameLayout heap liftedB.heap
  aMatrixPrefix : ∀ i : Fin 4, SameU64Prefix heap liftedB.heap
    (hgcdMatPtr M hM i) (hgcdMatLen M hM i)
  aLowAPrefix : SameU64Prefix heap liftedB.heap lowA lenLowA
  aLowBPrefix : SameU64Prefix heap liftedB.heap lowB lenLowB
  reconstructA : HgcdReconstructWorkspace liftedB.heap A T0
    (hgcdMatPtr M hM (3 : Fin 4)) (hgcdMatLen M hM (3 : Fin 4)) lowA lenLowA
    (hgcdMatPtr M hM (1 : Fin 4)) (hgcdMatLen M hM (1 : Fin 4)) lowB lenLowB
    scratch
  highALayout : RawHeap.SameLayout heap heap3
  highAPrefix : SameU64Prefix heap heap3 highA lenHighA
  liftA : HgcdLiftHighWorkspace heap3 A highA lowLenA shift lenHighA
  finalBLayout : RawHeap.SameLayout liftedB.heap liftedA.heap
  finalBPrefix : SameU64Prefix liftedB.heap liftedA.heap B liftedB.length

def HgcdRecursiveReconstructPairWorkspaceProvider (this : DenseUPolyZp)
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift : Nat)
    (M : HgcdMat) (hM : M.Valid) (sgn : Int) (heap : RawHeap) : Prop :=
  ∀ heap1 lowLenB liftedB heap3 lowLenA liftedA,
    hgcdRecursiveReconstructB this B T0 (hgcdMatPtr M hM (2 : Fin 4))
        (hgcdMatPtr M hM (0 : Fin 4)) lowA lowB scratch
        (hgcdMatLen M hM (2 : Fin 4)) (hgcdMatLen M hM (0 : Fin 4))
        lenLowA lenLowB sgn heap = .ok (heap1, lowLenB) →
    hgcdRecursiveLiftHigh this B highB lowLenB shift lenHighB heap1 =
        .ok liftedB →
    hgcdRecursiveReconstructA this A T0 (hgcdMatPtr M hM (3 : Fin 4))
        (hgcdMatPtr M hM (1 : Fin 4)) lowA lowB scratch
        (hgcdMatLen M hM (3 : Fin 4)) (hgcdMatLen M hM (1 : Fin 4))
        lenLowA lenLowB sgn liftedB.heap = .ok (heap3, lowLenA) →
    hgcdRecursiveLiftHigh this A highA lowLenA shift lenHighA heap3 =
        .ok liftedA →
    HgcdRecursiveReconstructPairWorkspace this A B T0 lowA lowB highA highB
      scratch lenLowA lenLowB lenHighA lenHighB shift M hM sgn heap heap1
      lowLenB liftedB heap3 lowLenA liftedA

/-- Complete raw semantic refinement of the actual paired reconstruction.
The returned polynomials arise from the four generated calls in source
order, including both final normalizations. -/
theorem hgcdRecursiveReconstructPair_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift : Nat)
    (M : HgcdMat) (hM : M.Valid) (sgn : Int) (heap : RawHeap)
    (result : HgcdRecursiveReconstructPairResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
      lowA lowB highA highB scratch lenLowA lenLowB lenHighA lenHighB shift
      M hM sgn heap)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hLowA : RawDensePolyRep this heap lowA lenLowA polyLowA)
    (hLowB : RawDensePolyRep this heap lowB lenLowB polyLowB)
    (hHighA : RawDensePolyRep this heap highA lenHighA polyHighA)
    (hHighB : RawDensePolyRep this heap highB lenHighB polyHighB)
    (hrun : hgcdRecursiveReconstructPair this A B T0 lowA lowB highA highB
      scratch lenLowA lenLowB lenHighA lenHighB shift M hM sgn heap =
      .ok result) :
    RawDensePolyRep this result.heap A result.lenA
        (hgcdReconstructedLowA entries polyLowA polyLowB sgn +
          Polynomial.X ^ shift * polyHighA) ∧
      RawDensePolyRep this result.heap B result.lenB
        (hgcdReconstructedLowB entries polyLowA polyLowB sgn +
          Polynomial.X ^ shift * polyHighB) ∧
      result.lenB ≤ max (shift + lenHighB)
        (max
          (hgcdMatLen M hM (2 : Fin 4) + lenLowA)
          (hgcdMatLen M hM (0 : Fin 4) + lenLowB)) := by
  rcases hgcdRecursiveReconstructPair_exec this A B T0 lowA lowB highA
      highB scratch lenLowA lenLowB lenHighA lenHighB shift M hM sgn heap
      result hrun with
    ⟨heap1, lowLenB, liftedB, heap3, lowLenA, liftedA,
      hBRun, hLiftBRun, hARun, hLiftARun, hHeap, hLenA, hLenB⟩
  have hwork := physical heap1 lowLenB liftedB heap3 lowLenA liftedA hBRun
    hLiftBRun hARun hLiftARun
  rcases hgcdRecursiveReconstructB_refines this B T0
      (hgcdMatPtr M hM (2 : Fin 4)) (hgcdMatPtr M hM (0 : Fin 4)) lowA
      lowB scratch (hgcdMatLen M hM (2 : Fin 4))
      (hgcdMatLen M hM (0 : Fin 4)) lenLowA lenLowB sgn heap
      (entries 2) (entries 0) polyLowA polyLowB hcfg hp
      hwork.reconstructB (hMatrix 2) (hMatrix 0) hLowA hLowB with
    ⟨heapB, lenB0, hBRun', _, hB0, hLowLenB⟩
  have hEqB : (heapB, lenB0) = (heap1, lowLenB) :=
    Except.ok.inj (hBRun'.symm.trans hBRun)
  cases hEqB
  have hHighB1 := rawDensePolyRep_of_same_prefix this heap heap1 highB
    lenHighB polyHighB hwork.highBLayout hwork.highBPrefix hHighB
  have hpWord : this._p ≠ 0 := by
    intro hzero
    have := congrArg UInt64.toNat hzero
    simp at this
    omega
  rcases hgcdRecursiveLiftHigh_terminates this B highB lowLenB shift lenHighB
      heap1 (hgcdReconstructedLowB entries polyLowA polyLowB sgn) hpWord
      hwork.liftB (by simpa [hgcdReconstructedLowB] using hB0) with
    ⟨boundedB, hBoundedBRun, _, hBoundB⟩
  have hEqBoundedB : boundedB = liftedB :=
    Except.ok.inj (hBoundedBRun.symm.trans hLiftBRun)
  subst boundedB
  rcases hgcdRecursiveLiftHigh_refines this B highB lowLenB shift lenHighB
      heap1 (hgcdReconstructedLowB entries polyLowA polyLowB sgn) polyHighB
      hpWord hwork.liftB (by simpa [hgcdReconstructedLowB] using hB0)
      hHighB1 with ⟨liftedB', hLiftBRun', _, hBFinal⟩
  have hEqLiftB : liftedB' = liftedB :=
    Except.ok.inj (hLiftBRun'.symm.trans hLiftBRun)
  subst liftedB'
  have matrixAtB (i : Fin 4) : RawDensePolyRep this liftedB.heap
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i) (entries i) :=
    rawDensePolyRep_of_same_prefix this heap liftedB.heap
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i) (entries i)
      hwork.aInputLayout (hwork.aMatrixPrefix i) (hMatrix i)
  have lowAAtB := rawDensePolyRep_of_same_prefix this heap liftedB.heap lowA
    lenLowA polyLowA hwork.aInputLayout hwork.aLowAPrefix hLowA
  have lowBAtB := rawDensePolyRep_of_same_prefix this heap liftedB.heap lowB
    lenLowB polyLowB hwork.aInputLayout hwork.aLowBPrefix hLowB
  rcases hgcdRecursiveReconstructA_refines this A T0
      (hgcdMatPtr M hM (3 : Fin 4)) (hgcdMatPtr M hM (1 : Fin 4)) lowA
      lowB scratch (hgcdMatLen M hM (3 : Fin 4))
      (hgcdMatLen M hM (1 : Fin 4)) lenLowA lenLowB sgn liftedB.heap
      (entries 3) (entries 1) polyLowA polyLowB hcfg hp hwork.reconstructA
      (matrixAtB 3) (matrixAtB 1) lowAAtB lowBAtB with
    ⟨heapA, lenA0, hARun', _, hA0, _⟩
  have hEqA : (heapA, lenA0) = (heap3, lowLenA) :=
    Except.ok.inj (hARun'.symm.trans hARun)
  cases hEqA
  have hHighA3 := rawDensePolyRep_of_same_prefix this heap heap3 highA
    lenHighA polyHighA hwork.highALayout hwork.highAPrefix hHighA
  rcases hgcdRecursiveLiftHigh_refines this A highA lowLenA shift lenHighA
      heap3 (hgcdReconstructedLowA entries polyLowA polyLowB sgn) polyHighA
      hpWord hwork.liftA (by simpa [hgcdReconstructedLowA] using hA0)
      hHighA3 with ⟨liftedA', hLiftARun', _, hAFinal⟩
  have hEqLiftA : liftedA' = liftedA :=
    Except.ok.inj (hLiftARun'.symm.trans hLiftARun)
  subst liftedA'
  have hBAtFinal := rawDensePolyRep_of_same_prefix this liftedB.heap
    liftedA.heap B liftedB.length
    (hgcdReconstructedLowB entries polyLowA polyLowB sgn +
      Polynomial.X ^ shift * polyHighB) hwork.finalBLayout
      hwork.finalBPrefix hBFinal
  rw [hHeap, hLenA, hLenB]
  refine ⟨hAFinal, hBAtFinal, ?_⟩
  exact hBoundB.trans (max_le_max (Nat.le_refl _) hLowLenB)

/-- The physical paired-reconstruction bound closes against an enclosing
input length once its shifted high part and both real low products fit that
input.  This is the arithmetic form consumed by the well-founded recursive
call, not a runtime decrease test. -/
theorem hgcdRecursiveReconstructPair_lenB_le_input
    (resultLen inputLength shift lenHighB lenR2 lenLowA lenR0 lenLowB : Nat)
    (hresult : resultLen ≤ max (shift + lenHighB)
      (max (lenR2 + lenLowA) (lenR0 + lenLowB)))
    (hhigh : shift + lenHighB ≤ inputLength)
    (hsecond : lenR2 + lenLowA ≤ inputLength)
    (hzero : lenR0 + lenLowB ≤ inputLength) :
    resultLen ≤ inputLength := by
  exact hresult.trans (max_le hhigh (max_le hsecond hzero))

/-- Close the first reconstruction bound from the returned first-HGCD
matrix invariant.  Here `highLength = inputLength - inputLength / 2` and
the low slices have exactly the lengths selected by the source. -/
theorem hgcdRecursiveFirstReconstruct_lenB_le_input
    (resultLen inputLength inputLengthB returnedLenA returnedLenB lenR2 lenR0 : Nat)
    (hresult : resultLen ≤
      max (inputLength / 2 + returnedLenB)
        (max
          (lenR2 + Nat.min inputLength (inputLength / 2))
          (lenR0 + Nat.min inputLengthB (inputLength / 2))))
    (hreturnedOrder : returnedLenB ≤ returnedLenA)
    (hreturnedBound : returnedLenA ≤ inputLength - inputLength / 2)
    (hreturnedPos : 0 < returnedLenA)
    (hrow2 : lenR2 + returnedLenA ≤
      (inputLength - inputLength / 2) + 1)
    (hrow0 : lenR0 + returnedLenA ≤
      (inputLength - inputLength / 2) + 1) :
    resultLen ≤ inputLength := by
  have hsplit : inputLength / 2 +
      (inputLength - inputLength / 2) = inputLength := by omega
  have hhigh : inputLength / 2 + returnedLenB ≤ inputLength := by omega
  have hsecond : lenR2 + Nat.min inputLength (inputLength / 2) ≤
      inputLength := by
    have hmin : Nat.min inputLength (inputLength / 2) ≤
        inputLength / 2 := Nat.min_le_right _ _
    omega
  have hzero : lenR0 + Nat.min inputLengthB (inputLength / 2) ≤
      inputLength := by
    have hmin : Nat.min inputLengthB (inputLength / 2) ≤
        inputLength / 2 := Nat.min_le_right _ _
    omega
  exact hgcdRecursiveReconstructPair_lenB_le_input resultLen inputLength
    (inputLength / 2) returnedLenB lenR2
    (Nat.min inputLength (inputLength / 2)) lenR0
    (Nat.min inputLengthB (inputLength / 2)) hresult hhigh hsecond hzero

/-- Concrete discharge of the proof-only first-reconstruction bound used by
the strictly decreasing recursive body.  Every premise is either a raw
representation, a physical workspace fact, or the length invariant returned
by the real first HGCD call. -/
theorem hgcdRecursiveFirstReconstruct_bound_of_invariant
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (A B T0 a b highA highB scratch : RawPtr UInt64)
    (lenA lenB : Nat) (first : HgcdRecursiveResult)
    (result : HgcdRecursiveReconstructPairResult)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hinvariant : HgcdRecursiveLengthInvariant (lenA - lenA / 2) first)
    (physical : HgcdRecursiveReconstructPairWorkspaceProvider this A B T0
      a b highA highB scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap)
    (hMatrix : HgcdMatRawDenseRep this first.heap first.matrix entries
      first.valid)
    (hLowA : RawDensePolyRep this first.heap a
      (Nat.min lenA (lenA / 2)) polyLowA)
    (hLowB : RawDensePolyRep this first.heap b
      (Nat.min lenB (lenA / 2)) polyLowB)
    (hHighA : RawDensePolyRep this first.heap highA first.lenA polyHighA)
    (hHighB : RawDensePolyRep this first.heap highB first.lenB polyHighB)
    (hrun : hgcdRecursiveReconstructPair this A B T0 a b highA highB scratch
      (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2)) first.lenA
      first.lenB (lenA / 2) first.matrix first.valid first.sgn first.heap =
        .ok result) :
    result.lenB ≤ lenA := by
  have hrefines := hgcdRecursiveReconstructPair_refines this A B T0 a b
    highA highB scratch (Nat.min lenA (lenA / 2))
    (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
    first.matrix first.valid first.sgn first.heap result entries polyLowA
    polyLowB polyHighA polyHighB hcfg hp physical hMatrix hLowA hLowB
    hHighA hHighB hrun
  exact hgcdRecursiveFirstReconstruct_lenB_le_input result.lenB lenA lenB
    first.lenA first.lenB
    (hgcdMatLen first.matrix first.valid (2 : Fin 4))
    (hgcdMatLen first.matrix first.valid (0 : Fin 4)) hrefines.2.2
    hinvariant.order hinvariant.inputBound hinvariant.positive
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row2A)
    (by simpa [hgcdMatLen, hgcdMatLenRaw] using hinvariant.row0A)

/-- Purely physical obligations for the exact final matrix block.  Besides
the two existing generated-call workspaces, the frame fields state that the
quotient update does not alter any buffer of the left matrix `R`. -/
structure HgcdRecursiveCombineMatrixWorkspace (this : DenseUPolyZp)
    (R S : HgcdMat) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) (modified : HgcdMatQuotientResult) : Prop where
  quotient : HgcdMatApplyQuotientWorkspaceProvider this S hS q lenQ T
    scratch heap
  rightLayout : RawHeap.SameLayout heap modified.heap
  rightPrefix : ∀ i : Fin 4, SameU64Prefix heap modified.heap
    (hgcdMatPtr R hR i) (hgcdMatLen R hR i)
  multiply : HgcdMatMulLoopWorkspaceProvider this R modified.matrix hR
    modified.valid a2 scratch

def HgcdRecursiveCombineMatrixWorkspaceProvider (this : DenseUPolyZp)
    (R S : HgcdMat) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) : Prop :=
  ∀ modified,
    hgcdMatApplyQuotient this S hS q lenQ T scratch heap = .ok modified →
    HgcdRecursiveCombineMatrixWorkspace this R S hR hS q lenQ T a2 scratch
      heap modified

/-- End-to-end semantic refinement of the exact final C++ matrix block.
The proof consumes the real quotient update and the complete four-entry
matrix multiplication; the L2 result is obtained only after both generated
executions have succeeded. -/
theorem hgcdRecursiveCombineMatrix_refines (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) (result : HgcdMatMulResult)
    (right entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveCombineMatrixWorkspaceProvider this R S hR hS
      q lenQ T a2 scratch heap)
    (hRight : HgcdMatRawDenseRep this heap R right hR)
    (hSRep : HgcdMatRawDenseRep this heap S entries hS)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hrun : hgcdRecursiveCombineMatrix this M R S hM hR hS q lenQ T a2
      scratch heap = .ok result) :
    HgcdMatRawDenseRep this result.heap R right hR ∧
      ∃ hResult : result.matrix.Valid,
        HgcdMatRawDenseRep this result.heap result.matrix
          (hgcdMatProductEntry right
            (hgcdMatApplyQuotientEntries entries quotient)) hResult := by
  rcases hgcdRecursiveCombineMatrix_exec this M R S hM hR hS q lenQ T a2
      scratch heap result hrun with ⟨modified, hmodified, hmul⟩
  have hwork := physical modified hmodified
  have hModified := hgcdMatApplyQuotient_refines this S hS q lenQ T scratch
    heap modified entries quotient hcfg hp hwork.quotient hQ hSRep hmodified
  have hRightModified :
      HgcdMatRawDenseRep this modified.heap R right hR := by
    intro i
    exact rawDensePolyRep_of_same_prefix this heap modified.heap
      (hgcdMatPtr R hR i) (hgcdMatLen R hR i) (right i)
      hwork.rightLayout (hwork.rightPrefix i) (hRight i)
  have hProduct := hgcdMatMul_refines this M R modified.matrix hM hR
    modified.valid a2 scratch modified.heap result right
    (hgcdMatApplyQuotientEntries entries quotient) hcfg hp hwork.multiply
    hRightModified hModified.1 hmul
  exact ⟨hProduct.1, hProduct.2.2⟩

/-- The generated recursive-HGCD base helper with matrix computation enabled
is exactly the already-refined iterator initialization prefix, modulo its
smaller return record. -/
theorem hgcdRecursiveBase_true_eq_init (M : HgcdMat)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap) :
    hgcdRecursiveBase M true A B a b lenA lenB heap =
      (hgcdIterInit M A B A A 0 a lenA b lenB heap).map
        HgcdIterState.toRecursiveBaseResult := by
  unfold hgcdRecursiveBase hgcdIterInit
  generalize hone : dense_upoly_zp__mat_one_ir M heap = one
  cases one with
  | error fault =>
    simp only [hone, ↓reduceIte, Except.map]
  | ok pair =>
    rcases pair with ⟨heap1, matrix⟩
    generalize hcopyA : heap1.copyU64 A a lenA = copyA
    cases copyA with
    | error fault =>
      simp only [hone, hcopyA, ↓reduceIte, Except.map]
    | ok heap2 =>
      generalize hcopyB : heap2.copyU64 B b lenB = copyB
      cases copyB with
      | error fault =>
        simp only [hone, hcopyA, hcopyB, ↓reduceIte, Except.map]
      | ok heap3 =>
        simp only [hone, hcopyA, hcopyB, ↓reduceIte, Except.map,
          HgcdIterState.toRecursiveBaseResult]

/-- Semantic refinement of the exact `_hgcd_recursive` base branch when the
source requests its matrix.  The two raw copies retain their C++ order. -/
theorem hgcdRecursiveBase_true_refines (this : DenseUPolyZp)
    (M : HgcdMat) (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (heap : RawHeap) (result : HgcdRecursiveBaseResult)
    (left right : Polynomial (ZMod this._p.toNat)) (hM : M.Valid)
    (hp : 1 < this._p.toNat)
    (h0 : heap.ValidU64Slice (hgcdMatPtr M hM (0 : Fin 4)) 1)
    (h3 : heap.ValidU64Slice (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (h03 : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1
      (hgcdMatPtr M hM (3 : Fin 4)) 1)
    (hA : heap.ValidU64Slice A lenA) (hB : heap.ValidU64Slice B lenB)
    (hAa : U64SlicesDisjoint A lenA a lenA)
    (hBb : U64SlicesDisjoint B lenB b lenB)
    (hAb : U64SlicesDisjoint A lenA b lenB)
    (hBA : U64SlicesDisjoint B lenB A lenA)
    (h0a : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 a lenA)
    (h3a : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 a lenA)
    (h0b : U64SlicesDisjoint (hgcdMatPtr M hM (0 : Fin 4)) 1 b lenB)
    (h3b : U64SlicesDisjoint (hgcdMatPtr M hM (3 : Fin 4)) 1 b lenB)
    (hAMatrix : ∀ i : Fin 4, U64SlicesDisjoint A lenA
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hBMatrix : ∀ i : Fin 4, U64SlicesDisjoint B lenB
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hMatrixValid : ∀ i : Fin 4, heap.ValidU64Slice
      (hgcdMatPtr M hM i) (identityEntryLen i))
    (hLeft : RawDensePolyRep this heap a lenA left)
    (hRight : RawDensePolyRep this heap b lenB right)
    (hrun : hgcdRecursiveBase M true A B a b lenA lenB heap = .ok result) :
    result.lenA = lenA ∧ result.lenB = lenB ∧ result.sgn = 1 ∧
      ∃ hResultM : result.matrix.Valid,
        HgcdMatRawDenseRep this result.heap result.matrix
          (identityEntries this._p.toNat) hResultM ∧
        RawDensePolyRep this result.heap A result.lenA left ∧
        RawDensePolyRep this result.heap B result.lenB right := by
  rcases hgcdIterInit_refines this M A B A A 0 a lenA b lenB heap left
      right hM hp h0 h3 h03 hA hB hAa hBb hAb hBA h0a h3a h0b h3b
      hAMatrix hBMatrix hMatrixValid hLeft hRight with
    ⟨initial, hinit, hIA, hlenIA, hIB, hlenIB, _, _, _, hsgn,
      hInitialM, hMatrix,
      hARep, hBRep, _, _⟩
  have hbase : hgcdRecursiveBase M true A B a b lenA lenB heap =
      .ok initial.toRecursiveBaseResult := by
    rw [hgcdRecursiveBase_true_eq_init, hinit]
    rfl
  have heq := Except.ok.inj (hbase.symm.trans hrun)
  subst result
  refine ⟨hlenIA, hlenIB, hsgn, hInitialM, hMatrix, ?_, ?_⟩
  · simpa [HgcdIterState.toRecursiveBaseResult, hIA, hlenIA] using hARep
  · simpa [HgcdIterState.toRecursiveBaseResult, hIB, hlenIB] using hBRep

/-- Semantic refinement of the exact `_hgcd_recursive` base branch used by
GCD when matrix output is disabled.  No matrix call or matrix specification
is introduced; only the two source-ordered raw copies are executed. -/
theorem hgcdRecursiveBase_false_refines (this : DenseUPolyZp)
    (M : HgcdMat) (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (heap : RawHeap) (result : HgcdRecursiveBaseResult)
    (left right : Polynomial (ZMod this._p.toNat))
    (hA : heap.ValidU64Slice A lenA) (hB : heap.ValidU64Slice B lenB)
    (hAa : U64SlicesDisjoint A lenA a lenA)
    (hBb : U64SlicesDisjoint B lenB b lenB)
    (hAb : U64SlicesDisjoint A lenA b lenB)
    (hBA : U64SlicesDisjoint B lenB A lenA)
    (hLeft : RawDensePolyRep this heap a lenA left)
    (hRight : RawDensePolyRep this heap b lenB right)
    (hrun : hgcdRecursiveBase M false A B a b lenA lenB heap = .ok result) :
    result.matrix = M ∧ result.lenA = lenA ∧ result.lenB = lenB ∧
      result.sgn = 1 ∧
      RawDensePolyRep this result.heap A result.lenA left ∧
      RawDensePolyRep this result.heap B result.lenB right := by
  rcases copyU64_refines_rawDense this heap A a lenA left hA hAa hLeft with
    ⟨heap1, hcopyA, hlayout1, hA1⟩
  have hsameB1 := copyU64_preserves_prefix heap heap1 A a b lenA lenB hA
    hLeft.1 hAb hcopyA
  have hRight1 := rawDensePolyRep_of_same_prefix this heap heap1 b lenB
    right hlayout1 hsameB1 hRight
  have hB1 : heap1.ValidU64Slice B lenB := (hlayout1 B lenB).mp hB
  rcases copyU64_refines_rawDense this heap1 B b lenB right hB1 hBb
      hRight1 with ⟨heap2, hcopyB, hlayout2, hB2⟩
  have hsameA2 := copyU64_preserves_prefix heap1 heap2 B b A lenB lenA
    hB1 hRight1.1 hBA hcopyB
  have hA2 := rawDensePolyRep_of_same_prefix this heap1 heap2 A lenA left
    hlayout2 hsameA2 hA1
  have hactual : hgcdRecursiveBase M false A B a b lenA lenB heap =
      .ok (.mk heap2 M lenA lenB 1) := by
    simp [hgcdRecursiveBase, hcopyA, hcopyB]
  have heq : HgcdRecursiveBaseResult.mk heap2 M lenA lenB 1 = result :=
    Except.ok.inj (hactual.symm.trans hrun)
  subst result
  exact ⟨rfl, rfl, rfl, rfl, hA2, hB2⟩

end CLPoly.Impl.StrictHGCDRawRefinement
