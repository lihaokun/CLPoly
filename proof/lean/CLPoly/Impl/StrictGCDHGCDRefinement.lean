import CLPoly.Generated.StrictGCDHGCD
import CLPoly.Impl.StrictHGCDCheckedRefinement
import CLPoly.Impl.StrictEuclidRefinement

set_option autoImplicit false

namespace CLPoly.Impl.StrictGCDHGCDRefinement

open Generated.StrictHGCD
open Generated.StrictGCDHGCD
open Generated.StrictPolyAddSub
open CLPoly.Impl.RawPolynomialRep
open CLPoly.Impl.StrictDivremRefinement
open CLPoly.Impl.StrictEuclidRefinement
open CLPoly.Impl.StrictHGCDRawRefinement
open CLPoly.Impl.StrictMulRefinement
open CLPoly.Impl.StrictPolyAddSubRefinement
open CLPoly.Impl.StrictWordArithmetic

/-- Total raw execution of generated `_poly_add` from allocation and aliasing
safety.  The proof follows its common loop, selected tail copy, and final
normalization; it does not assume a successful run. -/
theorem polyAdd_succeeds (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (hC : heap.ValidU64Slice C (max lenA lenB))
    (hA : heap.ValidU64Slice A lenA) (hB : heap.ValidU64Slice B lenB)
    (hAliasA : ExactOrDisjoint C A) (hAliasB : ExactOrDisjoint C B) :
    ∃ heap' outLen,
      dense_upoly_zp__poly_add_ir this C A lenA B lenB heap =
        .ok (heap', outLen) ∧ RawHeap.SameLayout heap heap' := by
  rcases Nat.lt_trichotomy lenA lenB with hlt | heq | hgt
  · have hmin : min lenA lenB = lenA := Nat.min_eq_left (Nat.le_of_lt hlt)
    have hmax : max lenA lenB = lenB := Nat.max_eq_right (Nat.le_of_lt hlt)
    have hCFull : heap.ValidU64Slice C lenB := by simpa [hmax] using hC
    rcases addCommonLoop_ok this C A B lenA 0 heap
        (heap.validU64Slice_mono C lenB lenA hCFull (Nat.le_of_lt hlt)) hA
        (heap.validU64Slice_mono B lenB lenA hB (Nat.le_of_lt hlt))
        (by omega) with ⟨heap1, hloop, hlayout1⟩
    rcases addRightLongTail this C A B lenA lenB heap heap1 hlt
        hCFull hA hB hAliasB hloop hlayout1 with
      ⟨heap2, htail, hlayout2, _, _⟩
    have hC2 : heap2.ValidU64Slice C lenB :=
      (hlayout2 C lenB).mp ((hlayout1 C lenB).mp hCFull)
    rcases normaliseU64_ok heap2 C lenB hC2 with ⟨outLen, hnorm, _⟩
    refine ⟨heap2, outLen, ?_, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
    simp [dense_upoly_zp__poly_add_ir, hmin, hmax, hloop, hlt, htail]
    simp [show ¬ lenB < lenA by omega, hnorm]
  · subst lenB
    rcases addCommonLoop_ok this C A B lenA 0 heap (by simpa using hC) hA hB
        (by omega) with ⟨heap1, hloop, hlayout1⟩
    have hC1 : heap1.ValidU64Slice C lenA :=
      (hlayout1 C lenA).mp (by simpa using hC)
    rcases normaliseU64_ok heap1 C lenA hC1 with ⟨outLen, hnorm, _⟩
    refine ⟨heap1, outLen, ?_, hlayout1⟩
    simp [dense_upoly_zp__poly_add_ir, hloop, hnorm]
  · have hmin : min lenA lenB = lenB := Nat.min_eq_right (Nat.le_of_lt hgt)
    have hmax : max lenA lenB = lenA := Nat.max_eq_left (Nat.le_of_lt hgt)
    have hCFull : heap.ValidU64Slice C lenA := by simpa [hmax] using hC
    rcases addCommonLoop_ok this C A B lenB 0 heap
        (heap.validU64Slice_mono C lenA lenB hCFull (Nat.le_of_lt hgt))
        (heap.validU64Slice_mono A lenA lenB hA (Nat.le_of_lt hgt)) hB
        (by omega) with ⟨heap1, hloop, hlayout1⟩
    rcases addLeftLongTail this C A B lenA lenB heap heap1 hgt
        hCFull hA hB hAliasA hloop hlayout1 with
      ⟨heap2, htail, hlayout2, _, _⟩
    have hC2 : heap2.ValidU64Slice C lenA :=
      (hlayout2 C lenA).mp ((hlayout1 C lenA).mp hCFull)
    rcases normaliseU64_ok heap2 C lenA hC2 with ⟨outLen, hnorm, _⟩
    refine ⟨heap2, outLen, ?_, fun ptr length =>
      (hlayout1 ptr length).trans (hlayout2 ptr length)⟩
    simp [dense_upoly_zp__poly_add_ir, hmin, hmax, hloop, hgt, htail,
      hnorm]

/-- Total raw execution of generated `_mat_row_update`.  The inactive branch
is the descriptor swap; the active branch executes the real raw multiplication
and then the total raw addition proved above. -/
theorem matRowUpdate_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (hM : M.Valid)
    (quotient entry0 entry1 : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (workspace : MatRowUpdateWorkspace M i0 i1 Q lenQ T t scratch heap hM)
    (hQRep : RawDensePolyRep this heap Q lenQ quotient)
    (hEntry0Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i0)
      (hgcdMatLen M hM i0) entry0)
    (hEntry1Rep : RawDensePolyRep this heap (hgcdMatPtr M hM i1)
      (hgcdMatLen M hM i1) entry1) :
    ∃ result, dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T
      lenT t scratch heap = .ok result := by
  let p0 := hgcdMatPtr M hM i0
  let p1 := hgcdMatPtr M hM i1
  let l0 := hgcdMatLen M hM i0
  let l1 := hgcdMatLen M hM i1
  by_cases hactive : lenQ ≠ 0 ∧ l0 ≠ 0
  · have hQPos : 0 < lenQ := Nat.pos_of_ne_zero hactive.1
    have h0Pos : 0 < l0 := Nat.pos_of_ne_zero hactive.2
    have hmul : ∃ heap1,
        Generated.StrictMul.dense_upoly_zp__mul_ir this T
          (if lenQ ≥ hgcdMatLen M hM i0 then Q else hgcdMatPtr M hM i0)
          (if lenQ ≥ hgcdMatLen M hM i0 then lenQ else
            hgcdMatLen M hM i0)
          (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatPtr M hM i0 else Q)
          (if lenQ ≥ hgcdMatLen M hM i0 then hgcdMatLen M hM i0 else
            lenQ) scratch heap =
            .ok heap1 ∧ RawHeap.SameLayout heap heap1 := by
      by_cases hge : lenQ ≥ l0
      · have hge' : lenQ ≥ hgcdMatLen M hM i0 := by simpa [l0] using hge
        rcases CLPoly.Impl.StrictMulRefinement.mul_refines_rawDense
            this T Q lenQ p0 l0 scratch heap quotient entry0 hcfg hp hQPos
            h0Pos hge (lt_of_le_of_lt (Nat.le_max_left _ _)
              workspace.lenWord)
            (by simpa [l0, Nat.max_eq_left hge'] using workspace.validT)
            (by simpa [l0, Nat.max_eq_left hge'] using workspace.validScratch)
            (by simpa [l0, Nat.max_eq_left hge'] using workspace.disjointTQ)
            (by simpa [p0, l0, Nat.max_eq_left hge'] using
              workspace.disjointTMatrix i0)
            (by simpa [l0, Nat.max_eq_left hge'] using
              workspace.disjointTScratch)
            (by simpa [l0, Nat.max_eq_left hge'] using
              workspace.disjointScratchQ)
            (by simpa [p0, l0, Nat.max_eq_left hge'] using
              workspace.disjointScratchMatrix i0)
            hQRep (by simpa [p0, l0] using hEntry0Rep) with
          ⟨heap1, hrun, hlayout, _⟩
        exact ⟨heap1, by simpa [hge', p0, l0] using hrun,
          hlayout⟩
      · have hle : lenQ ≤ l0 := by omega
        have hge' : ¬ lenQ ≥ hgcdMatLen M hM i0 := by simpa [l0] using hge
        have hle' : lenQ ≤ hgcdMatLen M hM i0 := by simpa [l0] using hle
        rcases CLPoly.Impl.StrictMulRefinement.mul_refines_rawDense
            this T p0 l0 Q lenQ scratch heap entry0 quotient hcfg hp h0Pos
            hQPos hle (lt_of_le_of_lt (Nat.le_max_right _ _)
              workspace.lenWord)
            (by simpa [l0, Nat.max_eq_right hle'] using workspace.validT)
            (by simpa [l0, Nat.max_eq_right hle'] using workspace.validScratch)
            (by simpa [p0, l0, Nat.max_eq_right hle'] using
              workspace.disjointTMatrix i0)
            (by simpa [l0, Nat.max_eq_right hle'] using workspace.disjointTQ)
            (by simpa [l0, Nat.max_eq_right hle'] using
              workspace.disjointTScratch)
            (by simpa [p0, l0, Nat.max_eq_right hle'] using
              workspace.disjointScratchMatrix i0)
            (by simpa [l0, Nat.max_eq_right hle'] using
              workspace.disjointScratchQ)
            (by simpa [p0, l0] using hEntry0Rep) hQRep with
          ⟨heap1, hrun, hlayout, _⟩
        exact ⟨heap1, by simpa [hge', p0, l0] using hrun, hlayout⟩
    rcases hmul with ⟨heap1, hmul, hlayoutMul⟩
    have hEntry1After := matRowUpdate_mul_preserves_entry1 this M hM i0 i1 Q
      T scratch lenQ heap heap1 entry1 hQPos h0Pos workspace.validT
      workspace.validScratch hQRep.1 hEntry0Rep.1 workspace.disjointScratchQ
      (workspace.disjointScratchMatrix i0) (workspace.disjointTMatrix i1)
      (workspace.disjointScratchMatrix i1) hEntry1Rep hlayoutMul (by
        simpa using hmul)
    have hProduct := matRowUpdate_mul_result this M hM i0 Q T scratch lenQ
      heap heap1 quotient entry0 hcfg hp hQPos h0Pos workspace.lenWord
      workspace.validT workspace.validScratch workspace.disjointTQ
      (workspace.disjointTMatrix i0) workspace.disjointTScratch
      workspace.disjointScratchQ (workspace.disjointScratchMatrix i0) hQRep
      hEntry0Rep (by simpa using hmul)
    have hAddOut := (hlayoutMul t
      (max l1 (lenQ + l0 - 1))).mp (by
        simpa [l0, l1] using workspace.validAddOutput)
    rcases polyAdd_succeeds this t p1 T l1 (lenQ + l0 - 1) heap1 hAddOut
        (by simpa [p1, l1] using hEntry1After.1) hProduct.2.1
        workspace.aliasEntry1 workspace.aliasProduct with
      ⟨heap2, sumLen, hadd, _⟩
    refine ⟨⟨heap2,
      { poly := (M.poly.set i1.val p0 (by rw [hM.1]; exact i1.isLt)).set
          i0.val t (by simp [hM.1]),
        len := (M.len.set i1.val l0 (by rw [hM.2]; exact i1.isLt)).set
          i0.val sumLen (by simp [hM.2]) },
      T, lenQ + l0 - 1, p1⟩, ?_⟩
    have hactiveRaw : lenQ ≠ 0 ∧
        M.len[i0.val]'(by rw [hM.2]; exact i0.isLt) ≠ 0 := by
      simpa [l0, hgcdMatLen] using hactive
    have hmulRaw := hmul
    simp only [hgcdMatPtr, hgcdMatLen] at hmulRaw
    have haddRaw := hadd
    simp only [p1, l1, l0, hgcdMatPtr, hgcdMatLen] at haddRaw
    have hvalid : M.poly.size = 4 ∧ M.len.size = 4 := hM
    rw [dense_upoly_zp__mat_row_update_ir]
    rw [dif_pos hvalid]
    dsimp only
    rw [if_pos hactiveRaw, hmulRaw]
    simp only
    rw [haddRaw]
    rfl
  · refine ⟨⟨heap,
      { poly := (M.poly.set i1.val p0 (by rw [hM.1]; exact i1.isLt)).set
          i0.val p1 (by simp [hM.1]),
        len := (M.len.set i1.val l0 (by rw [hM.2]; exact i1.isLt)).set
          i0.val l1 (by simp [hM.2]) },
      T, lenT, t⟩, ?_⟩
    have hactiveRaw : ¬ (lenQ ≠ 0 ∧
        M.len[i0.val]'(by rw [hM.2]; exact i0.isLt) ≠ 0) := by
      simpa [l0, hgcdMatLen] using hactive
    have hvalid : M.poly.size = 4 ∧ M.len.size = 4 := hM
    rw [dense_upoly_zp__mat_row_update_ir]
    rw [dif_pos hvalid]
    dsimp only
    rw [if_neg hactiveRaw]
    rfl

/-- Non-circular physical precondition for one complete HGCD iterator step.
The row workspaces are exposed only after the preceding *actual* generated
call has returned.  No field assumes either row update succeeds and no field
contains an L2 polynomial or result record. -/
structure HgcdIterationTotalWorkspace (this : DenseUPolyZp)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (state : HgcdIterState) (hM : state.matrix.Valid) : Type where
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
  row23Workspace : ∀ (heap1 : RawHeap) (lenQ lenR : Nat),
    Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
      state.A state.lenA state.B state.lenB W3 state.heap =
        .ok (heap1, lenQ, lenR) →
    MatRowUpdateWorkspace state.matrix (2 : Fin 4) (3 : Fin 4) Q lenQ
      state.A state.t scratch heap1 hM
  row01Workspace : ∀ (heap1 : RawHeap) (lenQ lenR : Nat)
    (row23 : MatRowUpdateResult) (h23 : row23.matrix.Valid),
    Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
      state.A state.lenA state.B state.lenB W3 state.heap =
        .ok (heap1, lenQ, lenR) →
    dense_upoly_zp__mat_row_update_ir this state.matrix (2 : Fin 4)
      (3 : Fin 4) Q lenQ state.A state.lenT state.t scratch heap1 =
        .ok row23 →
    MatRowUpdateWorkspace row23.matrix (0 : Fin 4) (1 : Fin 4) Q lenQ
      row23.T row23.t scratch row23.heap h23
  divisorGuard23 : ∀ (heap1 : RawHeap) (lenQ lenR : Nat)
    (row23 : MatRowUpdateResult) (_h23 : row23.matrix.Valid),
    Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
      state.A state.lenA state.B state.lenB W3 state.heap =
        .ok (heap1, lenQ, lenR) →
    dense_upoly_zp__mat_row_update_ir this state.matrix (2 : Fin 4)
      (3 : Fin 4) Q lenQ state.A state.lenT state.t scratch heap1 =
        .ok row23 →
    MatRowUpdateGuardWorkspace state.matrix (2 : Fin 4) Q lenQ state.A
      state.t scratch state.B state.lenB hM
  divisorGuard01 : ∀ (heap1 : RawHeap) (lenQ lenR : Nat)
    (row23 : MatRowUpdateResult) (h23 : row23.matrix.Valid),
    Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
      state.A state.lenA state.B state.lenB W3 state.heap =
        .ok (heap1, lenQ, lenR) →
    dense_upoly_zp__mat_row_update_ir this state.matrix (2 : Fin 4)
      (3 : Fin 4) Q lenQ state.A state.lenT state.t scratch heap1 =
        .ok row23 →
    MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4) Q lenQ row23.T
      row23.t scratch state.B state.lenB h23
  remainderGuard23 : ∀ (heap1 : RawHeap) (lenQ lenR : Nat)
    (row23 : MatRowUpdateResult) (_h23 : row23.matrix.Valid),
    Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
      state.A state.lenA state.B state.lenB W3 state.heap =
        .ok (heap1, lenQ, lenR) →
    dense_upoly_zp__mat_row_update_ir this state.matrix (2 : Fin 4)
      (3 : Fin 4) Q lenQ state.A state.lenT state.t scratch heap1 =
        .ok row23 →
    MatRowUpdateGuardWorkspace state.matrix (2 : Fin 4) Q lenQ state.A
      state.t scratch state.T lenR hM
  remainderGuard01 : ∀ (heap1 : RawHeap) (lenQ lenR : Nat)
    (row23 : MatRowUpdateResult) (h23 : row23.matrix.Valid),
    Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
      state.A state.lenA state.B state.lenB W3 state.heap =
        .ok (heap1, lenQ, lenR) →
    dense_upoly_zp__mat_row_update_ir this state.matrix (2 : Fin 4)
      (3 : Fin 4) Q lenQ state.A state.lenT state.t scratch heap1 =
        .ok row23 →
    MatRowUpdateGuardWorkspace row23.matrix (0 : Fin 4) Q lenQ row23.T
      row23.t scratch state.T lenR h23

/-- One nonterminal generated HGCD iterator step is total from staged physical
safety.  The returned semantic state is extracted from the same divrem and row
update executions exhibited by the theorem. -/
theorem hgcdIteration_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (left right currentA currentB : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (state : HgcdIterState) (hM : state.matrix.Valid)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hinvariant : HgcdIterRawInvariant this left right currentA currentB
      entries state hM)
    (hlenB : 0 < state.lenB)
    (physical : HgcdIterationTotalWorkspace this Q W3 scratch state hM) :
    ∃ heap1 lenQ lenR row23 row01 quotient remainder h01,
      Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q state.T
        state.A state.lenA state.B state.lenB W3 state.heap =
          .ok (heap1, lenQ, lenR) ∧
      dense_upoly_zp__mat_row_update_ir this state.matrix (2 : Fin 4)
        (3 : Fin 4) Q lenQ state.A state.lenT state.t scratch heap1 =
          .ok row23 ∧
      dense_upoly_zp__mat_row_update_ir this row23.matrix (0 : Fin 4)
        (1 : Fin 4) Q lenQ row23.T row23.lenT row23.t scratch row23.heap =
          .ok row01 ∧
      HgcdMatRawDenseRep this row01.heap row01.matrix
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
        h01 ∧
      RawDensePolyRep this row01.heap state.B state.lenB currentB ∧
      RawDensePolyRep this row01.heap state.T lenR remainder ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdTransform left right currentB
        remainder
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 0)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 1)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 2)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 3) ∧
      CLPoly.Impl.StrictHGCDRefinement.HgcdSignedDet (-state.sgn)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 0)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 1)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 2)
        (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries 3) ∧
      normalize (EuclideanDomain.gcd currentA currentB) =
        normalize (EuclideanDomain.gcd currentB remainder) ∧
      lenR < state.lenB := by
  rcases polyDivrem_next_state this Q state.T state.A state.B state.lenA
      state.lenB W3 state.heap currentA currentB hlenB hinvariant.aRep
      hinvariant.bRep physical.validQ physical.validR physical.validW3
      physical.quotientCapacity physical.rA physical.wA physical.wB
      physical.qB physical.qW physical.rW physical.rQ physical.rB hcfg with
    ⟨heap1, lenQ, lenR, quotient, remainder, hdiv, hQRep, hBRep1, hRRep1,
      _, _, _, _, _, hlt⟩
  have hMatrix1 := polyDivrem_preserves_hgcdMatRawDenseRep this state.matrix
    Q state.T state.A state.B state.lenA state.lenB W3 state.heap heap1 lenQ
    lenR entries hM hinvariant.aRep.1 hinvariant.bRep.1 physical.validQ
    physical.validR physical.validW3 physical.qMatrix physical.rMatrix
    physical.wMatrix hinvariant.matrixRep hdiv
  have hws23 := physical.row23Workspace heap1 lenQ lenR hdiv
  rcases matRowUpdate_succeeds this state.matrix (2 : Fin 4) (3 : Fin 4) Q
      lenQ state.A state.lenT state.t scratch heap1 hM quotient (entries 2)
      (entries 3) hcfg hp hws23 hQRep (hMatrix1 2) (hMatrix1 3) with
    ⟨row23, hrow23⟩
  have h23 := matRowUpdate_result_valid this state.matrix (2 : Fin 4)
    (3 : Fin 4) Q lenQ state.A state.lenT state.t scratch heap1 row23 hrow23
  have hQ23 := matRowUpdate_preserves_quotient_of_workspace this state.matrix
    (2 : Fin 4) (3 : Fin 4) Q lenQ state.A state.lenT state.t scratch heap1
    row23 hM quotient (entries 2) (entries 3) (by decide) hcfg hp hws23
    hQRep (hMatrix1 2) (hMatrix1 3) hrow23
  have hE0Old23 := matRowUpdate_preserves_matrix_entry_of_workspace this
    state.matrix (2 : Fin 4) (3 : Fin 4) (0 : Fin 4) Q lenQ state.A
    state.lenT state.t scratch heap1 row23 hM quotient (entries 2)
    (entries 3) (entries 0) (by decide) (by decide) hcfg hp hws23 hQRep
    (hMatrix1 2) (hMatrix1 3) (hMatrix1 0) hrow23
  have hE1Old23 := matRowUpdate_preserves_matrix_entry_of_workspace this
    state.matrix (2 : Fin 4) (3 : Fin 4) (1 : Fin 4) Q lenQ state.A
    state.lenT state.t scratch heap1 row23 hM quotient (entries 2)
    (entries 3) (entries 1) (by decide) (by decide) hcfg hp hws23 hQRep
    (hMatrix1 2) (hMatrix1 3) (hMatrix1 1) hrow23
  rcases matRowUpdate_preserves_other_descriptor this state.matrix
      (2 : Fin 4) (3 : Fin 4) (0 : Fin 4) Q lenQ state.A state.lenT state.t
      scratch heap1 row23 hM (by decide) (by decide) (by decide) hrow23 with
    ⟨h23zero, hptr0, hlen0⟩
  have hh23zero : h23zero = h23 := Subsingleton.elim _ _
  subst h23zero
  rcases matRowUpdate_preserves_other_descriptor this state.matrix
      (2 : Fin 4) (3 : Fin 4) (1 : Fin 4) Q lenQ state.A state.lenT state.t
      scratch heap1 row23 hM (by decide) (by decide) (by decide) hrow23 with
    ⟨h23one, hptr1, hlen1⟩
  have hh23one : h23one = h23 := Subsingleton.elim _ _
  subst h23one
  have hE0_23 : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (0 : Fin 4))
      (hgcdMatLen row23.matrix h23 (0 : Fin 4)) (entries 0) := by
    rw [hptr0, hlen0]
    exact hE0Old23
  have hE1_23 : RawDensePolyRep this row23.heap
      (hgcdMatPtr row23.matrix h23 (1 : Fin 4))
      (hgcdMatLen row23.matrix h23 (1 : Fin 4)) (entries 1) := by
    rw [hptr1, hlen1]
    exact hE1Old23
  have hws01 := physical.row01Workspace heap1 lenQ lenR row23 h23 hdiv
    hrow23
  rcases matRowUpdate_succeeds this row23.matrix (0 : Fin 4) (1 : Fin 4) Q
      lenQ row23.T row23.lenT row23.t scratch row23.heap h23 quotient
      (entries 0) (entries 1) hcfg hp hws01 hQ23 hE0_23 hE1_23 with
    ⟨row01, hrow01⟩
  let hworkspace : HgcdIterationWorkspace this Q W3 scratch state heap1 lenQ
      lenR row23 row01 hM := {
    validQ := physical.validQ
    validR := physical.validR
    validW3 := physical.validW3
    quotientCapacity := physical.quotientCapacity
    rA := physical.rA
    wA := physical.wA
    wB := physical.wB
    qB := physical.qB
    qW := physical.qW
    rW := physical.rW
    rQ := physical.rQ
    rB := physical.rB
    qMatrix := physical.qMatrix
    rMatrix := physical.rMatrix
    wMatrix := physical.wMatrix
    matrix23Valid := h23
    row23Workspace := hws23
    row01Workspace := hws01
    divisorGuard23 := physical.divisorGuard23 heap1 lenQ lenR row23 h23
      hdiv hrow23
    divisorGuard01 := physical.divisorGuard01 heap1 lenQ lenR row23 h23
      hdiv hrow23
    remainderGuard23 := physical.remainderGuard23 heap1 lenQ lenR row23 h23
      hdiv hrow23
    remainderGuard01 := physical.remainderGuard01 heap1 lenQ lenR row23 h23
      hdiv hrow23 }
  rcases hgcdIterationCalls_refine this state.matrix Q W3 scratch state.A
      state.B state.T state.lenA state.lenB state.A state.lenT state.t
      state.heap heap1 lenQ lenR row23 row01 hM h23 left right currentA
      currentB entries state.sgn hcfg hp hlenB hinvariant.aRep
      hinvariant.bRep hinvariant.matrixRep hinvariant.transform
      hinvariant.signedDet hworkspace.validQ hworkspace.validR
      hworkspace.validW3 hworkspace.quotientCapacity hworkspace.rA
      hworkspace.wA hworkspace.wB hworkspace.qB hworkspace.qW hworkspace.rW
      hworkspace.rQ hworkspace.rB hworkspace.qMatrix hworkspace.rMatrix
      hworkspace.wMatrix hworkspace.row23Workspace hworkspace.row01Workspace
      hworkspace.divisorGuard23 hworkspace.divisorGuard01
      hworkspace.remainderGuard23 hworkspace.remainderGuard01 hdiv hrow23
      hrow01 with
    ⟨quotient', remainder', h01, _, _, hMatrix01, hDivisor01,
      hRemainder01, htransform, hdet, hgcd, _, hlt'⟩
  exact ⟨heap1, lenQ, lenR, row23, row01, quotient', remainder', h01, hdiv,
    hrow23, hrow01, hMatrix01, hDivisor01, hRemainder01, htransform, hdet,
    hgcd, hlt'⟩

/-- Staged allocation safety at every semantically represented state reached
by the generated iterator.  The provider returns only the non-circular
physical record above. -/
def HgcdIterTotalWorkspaceProvider (this : DenseUPolyZp) (m : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch : RawPtr UInt64) : Type :=
  ∀ (left right currentA currentB : Polynomial (ZMod this._p.toNat))
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (state : HgcdIterState) (hM : state.matrix.Valid),
    state.lenB ≥ m + 1 →
    HgcdIterRawInvariant this left right currentA currentB entries state hM →
    HgcdIterationTotalWorkspace this Q W3 scratch state hM

/-- Total semantic execution of the actual well-founded generated HGCD
iterator.  Recursion follows the source state and is justified by the
remainder length returned by that same raw divrem call. -/
theorem hgcdIterLoop_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (m : Nat) (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (scratch : RawPtr UInt64)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdIterTotalWorkspaceProvider this m Q W3 scratch) :
    ∀ (left right currentA currentB : Polynomial (ZMod this._p.toNat))
      (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
      (state : HgcdIterState) (hM : state.matrix.Valid),
      HgcdIterRawInvariant this left right currentA currentB entries state hM →
      ∃ final finalA finalB finalEntries hFinalM,
        hgcdIterLoop this m Q W3 scratch state = .ok final ∧
        HgcdIterRawInvariant this left right finalA finalB finalEntries final
          hFinalM ∧
        normalize (EuclideanDomain.gcd currentA currentB) =
          normalize (EuclideanDomain.gcd finalA finalB) ∧
        final.lenB < m + 1
  | left, right, currentA, currentB, entries, state, hM, hinvariant => by
      by_cases hguard : state.lenB ≥ m + 1
      · have hlenB : 0 < state.lenB := by omega
        have hphysical := physical left right currentA currentB entries state
          hM hguard hinvariant
        rcases hgcdIteration_succeeds this Q W3 scratch left right currentA
            currentB entries state hM hcfg hp hinvariant hlenB hphysical with
          ⟨heap1, lenQ, lenR, row23, row01, quotient, remainder, h01, hdiv,
            hrow23, hrow01, hMatrix01, hDivisor01, hRemainder01, htransform,
            hdet, hgcdStep, hlt⟩
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
        have hnext : HgcdIterRawInvariant this left right currentB remainder
            (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
            next h01 := ⟨hMatrix01, hDivisor01, hRemainder01, htransform, hdet⟩
        rcases hgcdIterLoop_succeeds this m Q W3 scratch hcfg hp physical left
            right currentB remainder
            (CLPoly.Impl.StrictHGCDRefinement.hgcdStepEntries quotient entries)
            next h01 hnext with
          ⟨final, finalA, finalB, finalEntries, hFinalM, htail,
            hfinalInvariant, hgcdRest, hstop⟩
        have hrow23Any : ∀ (i j : Fin 4), i.val = 2 → j.val = 3 →
            dense_upoly_zp__mat_row_update_ir this state.matrix i j Q lenQ
              state.A state.lenT state.t scratch heap1 = .ok row23 := by
          intro i j hi hj
          have hieq : i = (2 : Fin 4) := Fin.ext hi
          have hjeq : j = (3 : Fin 4) := Fin.ext hj
          subst i
          subst j
          exact hrow23
        have hrow01Any : ∀ (i j : Fin 4), i.val = 0 → j.val = 1 →
            dense_upoly_zp__mat_row_update_ir this row23.matrix i j Q lenQ
              row23.T row23.lenT row23.t scratch row23.heap = .ok row01 := by
          intro i j hi hj
          have hieq : i = (0 : Fin 4) := Fin.ext hi
          have hjeq : j = (1 : Fin 4) := Fin.ext hj
          subst i
          subst j
          exact hrow01
        refine ⟨final, finalA, finalB, finalEntries, hFinalM, ?_,
          hfinalInvariant, hgcdStep.trans hgcdRest, hstop⟩
        rw [hgcdIterLoop]
        rw [if_pos hguard, hdiv]
        simp only
        rw [hrow23Any _ _ rfl rfl]
        simp only
        rw [hrow01Any _ _ rfl rfl]
        simpa only [next] using htail
      · have hstop : state.lenB < m + 1 := by omega
        refine ⟨state, currentA, currentB, entries, hM, ?_, hinvariant, rfl,
          hstop⟩
        rw [hgcdIterLoop, if_neg hguard]
termination_by left right currentA currentB entries state hM _hinvariant =>
  state.lenB
decreasing_by exact hlt

/-- Total execution of the complete generated `_hgcd_iter`: the exact
identity initialization and ordered copies feed the total well-founded loop
above. -/
theorem hgcdIter_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M : HgcdMat) (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (heap : RawHeap) (left right : Polynomial (ZMod this._p.toNat))
    (hM : M.Valid) (hcfg : DensePreinvConfigured this)
    (hp : 1 < this._p.toNat)
    (physical : HgcdIterTotalWorkspaceProvider this (lenA / 2) Q W3 scratch)
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
    (hRight : RawDensePolyRep this heap b lenB right) :
    ∃ final finalA finalB finalEntries hFinalM,
      dense_upoly_zp__hgcd_iter_ir this M A B T t lenT a lenA b lenB Q W3
        scratch heap = .ok final ∧
      HgcdIterRawInvariant this left right finalA finalB finalEntries final
        hFinalM ∧
      normalize (EuclideanDomain.gcd left right) =
        normalize (EuclideanDomain.gcd finalA finalB) ∧
      final.lenB < lenA / 2 + 1 := by
  rcases hgcdIterInit_refines this M A B T t lenT a lenA b lenB heap left
      right hM hp h0 h3 h03 hA hB hAa hBb hAb hBA h0a h3a h0b h3b
      hAMatrix hBMatrix hMatrixValid hLeft hRight with
    ⟨initial, hinit, _, _, _, _, _, _, _, _, hInitialM, hInitialMatrix,
      hInitialA, hInitialB, hInitialTransform, hInitialDet⟩
  have hInitialInvariant : HgcdIterRawInvariant this left right left right
      (identityEntries this._p.toNat) initial hInitialM :=
    ⟨hInitialMatrix, hInitialA, hInitialB, hInitialTransform, hInitialDet⟩
  rcases hgcdIterLoop_succeeds this (lenA / 2) Q W3 scratch hcfg hp physical
      left right left right (identityEntries this._p.toNat) initial hInitialM
      hInitialInvariant with
    ⟨final, finalA, finalB, finalEntries, hFinalM, hloop, hFinalInvariant,
      hgcd, hstop⟩
  refine ⟨final, finalA, finalB, finalEntries, hFinalM, ?_, hFinalInvariant,
    hgcd, hstop⟩
  simp [dense_upoly_zp__hgcd_iter_ir, hinit, hloop]

/-- Non-circular staged safety for the four source-ordered calls in one
recursive reconstruction pair.  Every later workspace is indexed only by
executions that precede it in the generated C++ control flow. -/
structure HgcdRecursiveReconstructPairTotalWorkspace (this : DenseUPolyZp)
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift : Nat)
    (M : HgcdMat) (hM : M.Valid) (sgn : Int) (heap : RawHeap) : Type where
  reconstructB : HgcdReconstructWorkspace heap B T0
    (hgcdMatPtr M hM (2 : Fin 4)) (hgcdMatLen M hM (2 : Fin 4)) lowA lenLowA
    (hgcdMatPtr M hM (0 : Fin 4)) (hgcdMatLen M hM (0 : Fin 4)) lowB lenLowB
    scratch
  afterB : ∀ (heap1 : RawHeap) (lowLenB : Nat),
    hgcdRecursiveReconstructB this B T0 (hgcdMatPtr M hM (2 : Fin 4))
      (hgcdMatPtr M hM (0 : Fin 4)) lowA lowB scratch
      (hgcdMatLen M hM (2 : Fin 4)) (hgcdMatLen M hM (0 : Fin 4))
      lenLowA lenLowB sgn heap = .ok (heap1, lowLenB) →
    RawHeap.SameLayout heap heap1 ∧
      SameU64Prefix heap heap1 highB lenHighB ∧
      HgcdLiftHighWorkspace heap1 B highB lowLenB shift lenHighB
  afterLiftB : ∀ (heap1 : RawHeap) (lowLenB : Nat)
    (liftedB : HgcdLiftHighResult),
    hgcdRecursiveReconstructB this B T0 (hgcdMatPtr M hM (2 : Fin 4))
      (hgcdMatPtr M hM (0 : Fin 4)) lowA lowB scratch
      (hgcdMatLen M hM (2 : Fin 4)) (hgcdMatLen M hM (0 : Fin 4))
      lenLowA lenLowB sgn heap = .ok (heap1, lowLenB) →
    hgcdRecursiveLiftHigh this B highB lowLenB shift lenHighB heap1 =
      .ok liftedB →
    RawHeap.SameLayout heap liftedB.heap ∧
      (∀ i : Fin 4, SameU64Prefix heap liftedB.heap
        (hgcdMatPtr M hM i) (hgcdMatLen M hM i)) ∧
      SameU64Prefix heap liftedB.heap lowA lenLowA ∧
      SameU64Prefix heap liftedB.heap lowB lenLowB ∧
      HgcdReconstructWorkspace liftedB.heap A T0
        (hgcdMatPtr M hM (3 : Fin 4)) (hgcdMatLen M hM (3 : Fin 4))
        lowA lenLowA (hgcdMatPtr M hM (1 : Fin 4))
        (hgcdMatLen M hM (1 : Fin 4)) lowB lenLowB scratch
  afterA : ∀ (heap1 : RawHeap) (lowLenB : Nat)
    (liftedB : HgcdLiftHighResult) (heap3 : RawHeap) (lowLenA : Nat),
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
    RawHeap.SameLayout heap heap3 ∧
      SameU64Prefix heap heap3 highA lenHighA ∧
      HgcdLiftHighWorkspace heap3 A highA lowLenA shift lenHighA
  finalFrame : ∀ (heap1 : RawHeap) (lowLenB : Nat)
    (liftedB : HgcdLiftHighResult) (heap3 : RawHeap) (lowLenA : Nat)
    (liftedA : HgcdLiftHighResult),
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
    RawHeap.SameLayout liftedB.heap liftedA.heap ∧
      SameU64Prefix liftedB.heap liftedA.heap B liftedB.length ∧
      RawHeap.SameLayout heap liftedA.heap ∧
      (∀ i : Fin 4, SameU64Prefix heap liftedA.heap
        (hgcdMatPtr M hM i) (hgcdMatLen M hM i))

/-- Total semantic execution of the exact four-call reconstruction pair from
the staged physical contract. -/
theorem hgcdRecursiveReconstructPair_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift : Nat)
    (M : HgcdMat) (hM : M.Valid) (sgn : Int) (heap : RawHeap)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveReconstructPairTotalWorkspace this A B T0 lowA
      lowB highA highB scratch lenLowA lenLowB lenHighA lenHighB shift M hM
      sgn heap)
    (hMatrix : HgcdMatRawDenseRep this heap M entries hM)
    (hLowA : RawCanonicalPolySlice this heap lowA lenLowA polyLowA)
    (hLowB : RawCanonicalPolySlice this heap lowB lenLowB polyLowB)
    (hHighA : RawDensePolyRep this heap highA lenHighA polyHighA)
    (hHighB : RawDensePolyRep this heap highB lenHighB polyHighB) :
    ∃ result,
      hgcdRecursiveReconstructPair this A B T0 lowA lowB highA highB scratch
        lenLowA lenLowB lenHighA lenHighB shift M hM sgn heap = .ok result ∧
      RawDensePolyRep this result.heap A result.lenA
        (hgcdReconstructedLowA entries polyLowA polyLowB sgn +
          Polynomial.X ^ shift * polyHighA) ∧
      RawDensePolyRep this result.heap B result.lenB
        (hgcdReconstructedLowB entries polyLowA polyLowB sgn +
          Polynomial.X ^ shift * polyHighB) ∧
      HgcdMatRawDenseRep this result.heap M entries hM := by
  rcases hgcdRecursiveReconstructB_refines this B T0
      (hgcdMatPtr M hM (2 : Fin 4)) (hgcdMatPtr M hM (0 : Fin 4)) lowA
      lowB scratch (hgcdMatLen M hM (2 : Fin 4))
      (hgcdMatLen M hM (0 : Fin 4)) lenLowA lenLowB sgn heap (entries 2)
      (entries 0) polyLowA polyLowB hcfg hp physical.reconstructB
      (hMatrix 2) (hMatrix 0) hLowA hLowB with
    ⟨heap1, lowLenB, hBRun, _, hBReconstructed, _⟩
  rcases physical.afterB heap1 lowLenB hBRun with
    ⟨hLayoutB, hHighBPrefix, hLiftBWork⟩
  have hHighB1 : RawDensePolyRep this heap1 highB lenHighB polyHighB :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
      heap heap1 highB lenHighB polyHighB hLayoutB hHighBPrefix hHighB
  have hpWord : this._p ≠ 0 := by
    intro hzero
    have hzeroNat := congrArg UInt64.toNat hzero
    simp at hzeroNat
    omega
  rcases hgcdRecursiveLiftHigh_refines this B highB lowLenB shift lenHighB
      heap1
      (if sgn < 0 then entries 2 * polyLowA - entries 0 * polyLowB
       else entries 0 * polyLowB - entries 2 * polyLowA)
      polyHighB hpWord hLiftBWork hBReconstructed hHighB1 with
    ⟨liftedB, hLiftBRun, _, hBFinal⟩
  rcases physical.afterLiftB heap1 lowLenB liftedB hBRun hLiftBRun with
    ⟨hLayoutLiftB, hMatrixPrefix, hLowAPrefix, hLowBPrefix, hReconstructA⟩
  have hEntry3 : RawDensePolyRep this liftedB.heap
      (hgcdMatPtr M hM (3 : Fin 4)) (hgcdMatLen M hM (3 : Fin 4))
      (entries 3) :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
      heap liftedB.heap
      (hgcdMatPtr M hM (3 : Fin 4)) (hgcdMatLen M hM (3 : Fin 4))
      (entries 3) hLayoutLiftB (hMatrixPrefix 3) (hMatrix 3)
  have hEntry1 : RawDensePolyRep this liftedB.heap
      (hgcdMatPtr M hM (1 : Fin 4)) (hgcdMatLen M hM (1 : Fin 4))
      (entries 1) :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
      heap liftedB.heap
      (hgcdMatPtr M hM (1 : Fin 4)) (hgcdMatLen M hM (1 : Fin 4))
      (entries 1) hLayoutLiftB (hMatrixPrefix 1) (hMatrix 1)
  have hLowA2 := rawCanonicalPolySlice_of_same_prefix this heap liftedB.heap
    lowA lenLowA polyLowA hLayoutLiftB hLowAPrefix hLowA
  have hLowB2 := rawCanonicalPolySlice_of_same_prefix this heap liftedB.heap
    lowB lenLowB polyLowB hLayoutLiftB hLowBPrefix hLowB
  rcases hgcdRecursiveReconstructA_refines this A T0
      (hgcdMatPtr M hM (3 : Fin 4)) (hgcdMatPtr M hM (1 : Fin 4)) lowA
      lowB scratch (hgcdMatLen M hM (3 : Fin 4))
      (hgcdMatLen M hM (1 : Fin 4)) lenLowA lenLowB sgn liftedB.heap
      (entries 3) (entries 1) polyLowA polyLowB hcfg hp hReconstructA
      hEntry3 hEntry1 hLowA2 hLowB2 with
    ⟨heap3, lowLenA, hARun, _, hAReconstructed, _⟩
  rcases physical.afterA heap1 lowLenB liftedB heap3 lowLenA hBRun hLiftBRun
      hARun with ⟨hLayoutA, hHighAPrefix, hLiftAWork⟩
  have hHighA3 : RawDensePolyRep this heap3 highA lenHighA polyHighA :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
      heap heap3 highA lenHighA polyHighA hLayoutA hHighAPrefix hHighA
  rcases hgcdRecursiveLiftHigh_refines this A highA lowLenA shift lenHighA
      heap3
      (if sgn < 0 then entries 1 * polyLowB - entries 3 * polyLowA
       else entries 3 * polyLowA - entries 1 * polyLowB)
      polyHighA hpWord hLiftAWork hAReconstructed hHighA3 with
    ⟨liftedA, hLiftARun, _, hAFinal⟩
  rcases physical.finalFrame heap1 lowLenB liftedB heap3 lowLenA liftedA
      hBRun hLiftBRun hARun hLiftARun with
    ⟨hFinalBLayout, hFinalBPrefix, hFinalMatrixLayout, hFinalMatrixPrefix⟩
  let result : HgcdRecursiveReconstructPairResult :=
    ⟨liftedA.heap, liftedA.length, liftedB.length⟩
  have hBFinal' : RawDensePolyRep this liftedA.heap B liftedB.length
      (hgcdReconstructedLowB entries polyLowA polyLowB sgn +
        Polynomial.X ^ shift * polyHighB) :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
      liftedB.heap liftedA.heap B liftedB.length _ hFinalBLayout
      hFinalBPrefix (by
        simpa [hgcdReconstructedLowB] using hBFinal)
  have hMatrixFinal : HgcdMatRawDenseRep this liftedA.heap M entries hM := by
    intro i
    exact CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix
      this heap liftedA.heap
      (hgcdMatPtr M hM i) (hgcdMatLen M hM i) (entries i)
      hFinalMatrixLayout (hFinalMatrixPrefix i) (hMatrix i)
  have hBRunRaw : hgcdRecursiveReconstructB this B T0
      (hgcdMatPtrRaw M hM (2 : Fin 4)) (hgcdMatPtrRaw M hM (0 : Fin 4))
      lowA lowB scratch (hgcdMatLenRaw M hM (2 : Fin 4))
      (hgcdMatLenRaw M hM (0 : Fin 4)) lenLowA lenLowB sgn heap =
        .ok (heap1, lowLenB) := by
    simpa only [hgcdMatPtrRaw, hgcdMatPtr, hgcdMatLenRaw, hgcdMatLen] using
      hBRun
  have hARunRaw : hgcdRecursiveReconstructA this A T0
      (hgcdMatPtrRaw M hM (3 : Fin 4)) (hgcdMatPtrRaw M hM (1 : Fin 4))
      lowA lowB scratch (hgcdMatLenRaw M hM (3 : Fin 4))
      (hgcdMatLenRaw M hM (1 : Fin 4)) lenLowA lenLowB sgn liftedB.heap =
        .ok (heap3, lowLenA) := by
    simpa only [hgcdMatPtrRaw, hgcdMatPtr, hgcdMatLenRaw, hgcdMatLen] using
      hARun
  refine ⟨result, ?_, ?_, ?_, hMatrixFinal⟩
  · simp [result, hgcdRecursiveReconstructPair, hBRunRaw, hLiftBRun, hARunRaw,
      hLiftARun]
  · simpa [result, hgcdReconstructedLowA] using hAFinal
  · simpa [result] using hBFinal'

/-- Total execution of the exact middle divrem and source pointer arithmetic
between the two recursive HGCD children. -/
theorem hgcdRecursiveMiddle_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (W : RawPtr UInt64) (lenA m : Nat)
    (reconstructed : HgcdRecursiveReconstructPairResult)
    (polyA polyB : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this)
    (physical : HgcdRecursiveMiddleWorkspace W lenA reconstructed)
    (hlenB : 0 < reconstructed.lenB)
    (hA : RawDensePolyRep this reconstructed.heap
      (hgcdRecursiveWorkspace W lenA).a2 reconstructed.lenA polyA)
    (hB : RawDensePolyRep this reconstructed.heap
      (hgcdRecursiveWorkspace W lenA).b2 reconstructed.lenB polyB) :
    ∃ result quotient remainder,
      hgcdRecursiveMiddle this (hgcdRecursiveWorkspace W lenA).q
        (hgcdRecursiveWorkspace W lenA).d
        (hgcdRecursiveWorkspace W lenA).a2
        (hgcdRecursiveWorkspace W lenA).b2 reconstructed.lenA
        reconstructed.lenB m (hgcdRecursiveWorkspace W lenA).W3
        reconstructed.heap = .ok result ∧
      RawDensePolyRep this result.heap (hgcdRecursiveWorkspace W lenA).q
        result.lenQ quotient ∧
      RawDensePolyRep this result.heap (hgcdRecursiveWorkspace W lenA).b2
        reconstructed.lenB polyB ∧
      RawDensePolyRep this result.heap (hgcdRecursiveWorkspace W lenA).d
        result.lenD remainder ∧
      normalize (EuclideanDomain.gcd polyA polyB) =
        normalize (EuclideanDomain.gcd polyB remainder) ∧
      RawHeap.SameLayout reconstructed.heap result.heap ∧
      result.lenD < reconstructed.lenB := by
  let ws := hgcdRecursiveWorkspace W lenA
  rcases polyDivrem_next_state this ws.q ws.d ws.a2 ws.b2
      reconstructed.lenA reconstructed.lenB ws.W3 reconstructed.heap polyA
      polyB hlenB hA hB physical.validQ physical.validD physical.validW3
      physical.quotientCapacity physical.dA physical.wA physical.wB
      physical.qB physical.qW physical.dW physical.dQ physical.dB hcfg with
    ⟨heap1, lenQ, lenD, quotient, remainder, hdiv, hQRep, hBRep, hDRep, _,
      hgcd, hlayout, _, _, hlt⟩
  let result : HgcdRecursiveMiddleResult := {
    heap := heap1
    lenQ := lenQ
    lenD := lenD
    k := 2 * m - reconstructed.lenB + 1
    c0 := ws.b2.add (2 * m - reconstructed.lenB + 1)
    lenC0 := if reconstructed.lenB ≥ 2 * m - reconstructed.lenB + 1 then
      reconstructed.lenB - (2 * m - reconstructed.lenB + 1) else 0
    d0 := ws.d.add (2 * m - reconstructed.lenB + 1)
    lenD0 := if lenD ≥ 2 * m - reconstructed.lenB + 1 then
      lenD - (2 * m - reconstructed.lenB + 1) else 0 }
  refine ⟨result, quotient, remainder, ?_, ?_, ?_, ?_, hgcd, ?_, ?_⟩
  · simp [result, hgcdRecursiveMiddle, ws, hdiv]
  · simpa [result, ws] using hQRep
  · simpa [result, ws] using hBRep
  · simpa [result, ws] using hDRep
  · simpa [result] using hlayout
  · simpa [result] using hlt

/-- Total execution of one generated quotient-matrix column update.  The
active arithmetic is the same real raw multiplication/addition sequence used
by `_mat_row_update`; only the final descriptor construction differs. -/
theorem hgcdMatQuotientEntry_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (S : HgcdMat) (hS : S.Valid) (top bottom : Fin 4)
    (q : RawPtr UInt64) (lenQ : Nat) (T scratch : RawPtr UInt64)
    (heap : RawHeap) (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hne : top ≠ bottom)
    (workspace : HgcdMatQuotientEntryWorkspace S hS top bottom q lenQ T
      scratch heap)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hMatrix : HgcdMatRawDenseRep this heap S entries hS) :
    ∃ result,
      hgcdMatQuotientEntry this S hS top bottom q lenQ T scratch heap =
        .ok result := by
  by_cases hactive : 0 < lenQ ∧ 0 < hgcdMatLen S hS bottom
  · rcases matRowUpdate_succeeds this S bottom top q lenQ T 0
      (hgcdMatPtr S hS top) scratch heap hS quotient (entries bottom)
      (entries top) hcfg hp workspace hQ (hMatrix bottom) (hMatrix top) with
      ⟨row, hrow⟩
    rcases matRowUpdate_nonzero_success_shape this S bottom top q lenQ T 0
        (hgcdMatPtr S hS top) scratch heap row hS hne.symm
        (Nat.ne_of_gt hactive.1) (Nat.ne_of_gt hactive.2) hrow with
      ⟨heap1, heap2, sumLen, hmul, hadd, _, _, _, _, _⟩
    let len' := S.len.set top.val sumLen (by rw [hS.2]; exact top.isLt)
    let next : HgcdMat := { S with len := len' }
    have hNext : next.Valid := ⟨hS.1, by simp [next, len', hS.2]⟩
    let result : HgcdMatQuotientEntryResult := ⟨heap2, next, hNext⟩
    refine ⟨result, ?_⟩
    have hguardRaw : 0 < lenQ ∧ 0 < hgcdMatLenRaw S hS bottom := by
      simpa [hgcdMatLenRaw, hgcdMatLen] using hactive
    have hguardBool : (decide (lenQ > 0) &&
        decide (hgcdMatLenRaw S hS bottom > 0)) = true := by
      simp [hguardRaw]
    have hmulRaw : Generated.StrictMul.dense_upoly_zp__mul_ir this T
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then q else
          hgcdMatPtrRaw S hS bottom)
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then lenQ else
          hgcdMatLenRaw S hS bottom)
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then
          hgcdMatPtrRaw S hS bottom else q)
        (if lenQ ≥ hgcdMatLenRaw S hS bottom then
          hgcdMatLenRaw S hS bottom else lenQ) scratch heap = .ok heap1 := by
      simpa only [hgcdMatPtrRaw, hgcdMatPtr, hgcdMatLenRaw, hgcdMatLen] using
        hmul
    have haddRaw : dense_upoly_zp__poly_add_ir this
        (hgcdMatPtrRaw S hS top) (hgcdMatPtrRaw S hS top)
        (hgcdMatLenRaw S hS top) T
        (lenQ + hgcdMatLenRaw S hS bottom - 1) heap1 =
          .ok (heap2, sumLen) := by
      simpa only [hgcdMatPtrRaw, hgcdMatPtr, hgcdMatLenRaw, hgcdMatLen] using
        hadd
    rw [hgcdMatQuotientEntry]
    rw [if_pos hguardBool]
    dsimp only
    rw [hmulRaw]
    dsimp only
    rw [haddRaw]
  · let result : HgcdMatQuotientEntryResult := ⟨heap, S, hS⟩
    refine ⟨result, ?_⟩
    have hguardRaw : ¬(0 < lenQ ∧ 0 < hgcdMatLenRaw S hS bottom) := by
      simpa [hgcdMatLenRaw, hgcdMatLen] using hactive
    have hguardBool : (decide (lenQ > 0) &&
        decide (hgcdMatLenRaw S hS bottom > 0)) ≠ true := by
      simpa [Bool.and_eq_true] using hguardRaw
    rw [hgcdMatQuotientEntry, if_neg hguardBool]

/-- Staged, non-circular workspace for the two concrete quotient columns. -/
structure HgcdMatApplyQuotientTotalWorkspace (this : DenseUPolyZp)
    (S : HgcdMat) (hS : S.Valid) (q : RawPtr UInt64) (lenQ : Nat)
    (T scratch : RawPtr UInt64) (heap : RawHeap) : Type where
  first : HgcdMatQuotientEntryWorkspace (hgcdMatSwapRows S hS)
    (hgcdMatSwapRows_valid S hS) (0 : Fin 4) (2 : Fin 4) q lenQ T scratch
    heap
  second : ∀ first,
    hgcdMatQuotientEntry this (hgcdMatSwapRows S hS)
      (hgcdMatSwapRows_valid S hS) (0 : Fin 4) (2 : Fin 4) q lenQ T scratch
      heap = .ok first →
    HgcdMatQuotientEntryWorkspace first.matrix first.valid (1 : Fin 4)
      (3 : Fin 4) q lenQ T scratch first.heap

/-- Total semantic execution of the exact row swap and two generated quotient
column updates. -/
theorem hgcdMatApplyQuotient_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (S : HgcdMat) (hS : S.Valid) (q : RawPtr UInt64) (lenQ : Nat)
    (T scratch : RawPtr UInt64) (heap : RawHeap)
    (entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatApplyQuotientTotalWorkspace this S hS q lenQ T
      scratch heap)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hMatrix : HgcdMatRawDenseRep this heap S entries hS) :
    ∃ result,
      hgcdMatApplyQuotient this S hS q lenQ T scratch heap = .ok result ∧
      HgcdMatRawDenseRep this result.heap result.matrix
        (hgcdMatApplyQuotientEntries entries quotient) result.valid ∧
      RawDensePolyRep this result.heap q lenQ quotient := by
  let swapped := hgcdMatSwapRows S hS
  let hSwapped := hgcdMatSwapRows_valid S hS
  have hSwappedRep := hgcdMatSwapRows_refines this S hS entries heap hMatrix
  rcases hgcdMatQuotientEntry_succeeds this swapped hSwapped (0 : Fin 4)
      (2 : Fin 4) q lenQ T scratch heap (hgcdMatSwapEntries entries)
      quotient hcfg hp (by decide) physical.first hQ hSwappedRep with
    ⟨first, hfirst⟩
  have hFirstSemantic := hgcdMatQuotientEntry_refines this swapped hSwapped
    (0 : Fin 4) (2 : Fin 4) q lenQ T scratch heap first
    (hgcdMatSwapEntries entries) quotient hcfg hp physical.first hQ
    hSwappedRep hfirst
  have hSecondWork := physical.second first hfirst
  rcases hgcdMatQuotientEntry_succeeds this first.matrix first.valid
      (1 : Fin 4) (3 : Fin 4) q lenQ T scratch first.heap
      (hgcdMatQuotientUpdateEntries (hgcdMatSwapEntries entries) quotient
        (0 : Fin 4) (2 : Fin 4))
      quotient hcfg hp (by decide) hSecondWork hFirstSemantic.2
      hFirstSemantic.1 with ⟨second, hsecond⟩
  have hSecondSemantic := hgcdMatQuotientEntry_refines this first.matrix
    first.valid (1 : Fin 4) (3 : Fin 4) q lenQ T scratch first.heap second
    (hgcdMatQuotientUpdateEntries (hgcdMatSwapEntries entries) quotient
      (0 : Fin 4) (2 : Fin 4))
    quotient hcfg hp hSecondWork hFirstSemantic.2 hFirstSemantic.1 hsecond
  let result : HgcdMatQuotientResult :=
    ⟨second.heap, second.matrix, second.valid⟩
  refine ⟨result, ?_, ?_, ?_⟩
  · simp [result, hgcdMatApplyQuotient, swapped, hfirst, hsecond]
  · simpa [result, hgcdMatApplyQuotientEntries] using hSecondSemantic.1
  · simpa [result] using hSecondSemantic.2

/-- Non-circular physical stages for the two products and selected tail of
one generated `_mat_mul_entry`. -/
structure HgcdMatMulEntryTotalWorkspace (this : DenseUPolyZp)
    (heap : RawHeap) (C P Q R S T scratch : RawPtr UInt64)
    (lenP lenQ lenR lenS : Nat) : Type where
  first : HgcdMulTermWorkspace heap C P lenP Q lenQ scratch
  afterFirst : ∀ productPQ,
    hgcdRecursiveMulTerm this C P lenP Q lenQ scratch heap = .ok productPQ →
    HgcdMulTermWorkspace productPQ.heap T R lenR S lenS scratch ∧
      U64SlicesDisjoint C (hgcdMulCapacity lenP lenQ) R lenR ∧
      U64SlicesDisjoint scratch (8 * max lenP lenQ) R lenR ∧
      U64SlicesDisjoint C (hgcdMulCapacity lenP lenQ) S lenS ∧
      U64SlicesDisjoint scratch (8 * max lenP lenQ) S lenS
  afterSecond : ∀ productPQ productRS,
    hgcdRecursiveMulTerm this C P lenP Q lenQ scratch heap = .ok productPQ →
    hgcdRecursiveMulTerm this T R lenR S lenS scratch productPQ.heap =
      .ok productRS →
    U64SlicesDisjoint T (hgcdMulCapacity lenR lenS) C productPQ.length ∧
      U64SlicesDisjoint scratch (8 * max lenR lenS) C productPQ.length ∧
      productRS.heap.ValidU64Slice C
        (max productPQ.length productRS.length) ∧
      ExactOrDisjoint C T ∧
      U64SlicesDisjoint C productRS.length T productRS.length

/-- Total semantic execution of both guarded products and the exact one of
four source tails selected by their returned lengths. -/
theorem hgcdMatMulEntry_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (C P Q R S T scratch : RawPtr UInt64)
    (lenP lenQ lenR lenS : Nat) (heap : RawHeap)
    (polyP polyQ polyR polyS : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatMulEntryTotalWorkspace this heap C P Q R S T scratch
      lenP lenQ lenR lenS)
    (hP : RawDensePolyRep this heap P lenP polyP)
    (hQ : RawDensePolyRep this heap Q lenQ polyQ)
    (hR : RawDensePolyRep this heap R lenR polyR)
    (hS : RawDensePolyRep this heap S lenS polyS) :
    ∃ result,
      hgcdMatMulEntry this C P Q R S T scratch lenP lenQ lenR lenS heap =
        .ok result ∧
      RawDensePolyRep this result.heap C result.length
        (polyP * polyQ + polyR * polyS) := by
  rcases hgcdRecursiveMulTerm_refines this C P lenP Q lenQ scratch heap polyP
      polyQ hcfg hp physical.first hP hQ with
    ⟨productPQ, hPQ, hLayoutPQ, hPQRep⟩
  rcases physical.afterFirst productPQ hPQ with
    ⟨hSecondWork, hFirstDstR, hFirstScratchR, hFirstDstS, hFirstScratchS⟩
  have hRPrefix := hgcdRecursiveMulTerm_preserves_guard this C P lenP Q lenQ
    scratch R lenR heap productPQ physical.first hP.1 hQ.1 hFirstDstR
    hFirstScratchR hPQ
  have hSPrefix := hgcdRecursiveMulTerm_preserves_guard this C P lenP Q lenQ
    scratch S lenS heap productPQ physical.first hP.1 hQ.1 hFirstDstS
    hFirstScratchS hPQ
  have hR1 :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
      heap productPQ.heap R lenR polyR hLayoutPQ hRPrefix hR
  have hS1 :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
      heap productPQ.heap S lenS polyS hLayoutPQ hSPrefix hS
  rcases hgcdRecursiveMulTerm_refines this T R lenR S lenS scratch
      productPQ.heap polyR polyS hcfg hp hSecondWork hR1 hS1 with
    ⟨productRS, hRS, hLayoutRS, hRSRep⟩
  rcases physical.afterSecond productPQ productRS hPQ hRS with
    ⟨hSecondDstC, hSecondScratchC, hFinalCValid, hAddAliasT, hCopyCT⟩
  have hCPrefix := hgcdRecursiveMulTerm_preserves_guard this T R lenR S lenS
    scratch C productPQ.length productPQ.heap productRS hSecondWork hR1.1
    hS1.1 hSecondDstC hSecondScratchC hRS
  have hPQRep2 :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
      productPQ.heap productRS.heap C productPQ.length (polyP * polyQ)
      hLayoutRS hCPrefix hPQRep
  have oldPhysical : HgcdMatMulEntryWorkspaceProvider this heap C P Q R S T
      scratch lenP lenQ lenR lenS := by
    intro actualPQ actualRS hActualPQ hActualRS
    have hpqEq : actualPQ = productPQ :=
      Except.ok.inj (hActualPQ.symm.trans hPQ)
    subst actualPQ
    have hrsEq : actualRS = productRS :=
      Except.ok.inj (hActualRS.symm.trans hRS)
    subst actualRS
    exact {
      first := physical.first
      second := hSecondWork
      firstDstR := hFirstDstR
      firstScratchR := hFirstScratchR
      firstDstS := hFirstDstS
      firstScratchS := hFirstScratchS
      secondDstC := hSecondDstC
      secondScratchC := hSecondScratchC
      finalCValid := hFinalCValid
      addAliasT := hAddAliasT
      copyCT := hCopyCT }
  by_cases hPQPos : 0 < productPQ.length
  · by_cases hRSPos : 0 < productRS.length
    · rcases polyAdd_succeeds this C C T productPQ.length productRS.length
          productRS.heap hFinalCValid hPQRep2.1 hRSRep.1 (Or.inl rfl)
          hAddAliasT with ⟨heap3, length, hadd, _⟩
      let result : HgcdMatMulEntryResult := ⟨heap3, length⟩
      have hrun : hgcdMatMulEntry this C P Q R S T scratch lenP lenQ lenR
          lenS heap = .ok result := by
        simp [result, hgcdMatMulEntry, hPQ, hRS, hPQPos, hRSPos, hadd]
      exact ⟨result, hrun, hgcdMatMulEntry_refines this C P Q R S T scratch
        lenP lenQ lenR lenS heap result polyP polyQ polyR polyS hcfg hp
        oldPhysical hP hQ hR hS hrun⟩
    · let result : HgcdMatMulEntryResult :=
        ⟨productRS.heap, productPQ.length⟩
      have hrun : hgcdMatMulEntry this C P Q R S T scratch lenP lenQ lenR
          lenS heap = .ok result := by
        simp [result, hgcdMatMulEntry, hPQ, hRS, hPQPos, hRSPos]
      exact ⟨result, hrun, hgcdMatMulEntry_refines this C P Q R S T scratch
        lenP lenQ lenR lenS heap result polyP polyQ polyR polyS hcfg hp
        oldPhysical hP hQ hR hS hrun⟩
  · by_cases hRSPos : 0 < productRS.length
    · have hCValid := productRS.heap.validU64Slice_mono C
          (max productPQ.length productRS.length) productRS.length
          hFinalCValid (Nat.le_max_right _ _)
      rcases copyU64_refines_rawDense this productRS.heap C T productRS.length
          (polyR * polyS) hCValid hCopyCT hRSRep with
        ⟨heap3, hcopy, _, _⟩
      let result : HgcdMatMulEntryResult := ⟨heap3, productRS.length⟩
      have hrun : hgcdMatMulEntry this C P Q R S T scratch lenP lenQ lenR
          lenS heap = .ok result := by
        simp [result, hgcdMatMulEntry, hPQ, hRS, hPQPos, hRSPos, hcopy]
      exact ⟨result, hrun, hgcdMatMulEntry_refines this C P Q R S T scratch
        lenP lenQ lenR lenS heap result polyP polyQ polyR polyS hcfg hp
        oldPhysical hP hQ hR hS hrun⟩
    · let result : HgcdMatMulEntryResult := ⟨productRS.heap, 0⟩
      have hrun : hgcdMatMulEntry this C P Q R S T scratch lenP lenQ lenR
          lenS heap = .ok result := by
        simp [result, hgcdMatMulEntry, hPQ, hRS, hPQPos, hRSPos]
      exact ⟨result, hrun, hgcdMatMulEntry_refines this C P Q R S T scratch
        lenP lenQ lenR lenS heap result polyP polyQ polyR polyS hcfg hp
        oldPhysical hP hQ hR hS hrun⟩

/-- Total entry workspaces plus the purely spatial frame provider for every
concrete state of the four-entry matrix multiplication loop. -/
structure HgcdMatMulLoopTotalWorkspace (this : DenseUPolyZp)
    (A B : HgcdMat) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) : Type where
  conditional : HgcdMatMulLoopWorkspaceProvider this A B hA hB T scratch
  entry : ∀ (C : HgcdMat) (hC : C.Valid) (i : Nat) (hi : i < 4)
    (heap : RawHeap),
    HgcdMatMulEntryTotalWorkspace this heap
      (hgcdMatPtr C hC ⟨i, hi⟩)
      (hgcdMatPtr A hA ⟨2 * (i / 2), by omega⟩)
      (hgcdMatPtr B hB ⟨i % 2, by omega⟩)
      (hgcdMatPtr A hA ⟨2 * (i / 2) + 1, by omega⟩)
      (hgcdMatPtr B hB ⟨2 + i % 2, by omega⟩) T scratch
      (hgcdMatLen A hA ⟨2 * (i / 2), by omega⟩)
      (hgcdMatLen B hB ⟨i % 2, by omega⟩)
      (hgcdMatLen A hA ⟨2 * (i / 2) + 1, by omega⟩)
      (hgcdMatLen B hB ⟨2 + i % 2, by omega⟩)

/-- Total execution of the generated four-entry matrix loop. -/
theorem hgcdMatMulLoop_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (A B : HgcdMat) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64)
    (left right : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatMulLoopTotalWorkspace this A B hA hB T scratch) :
    ∀ (C : HgcdMat) (hC : C.Valid) (i : Nat) (heap : RawHeap),
      HgcdMatRawDenseRep this heap A left hA →
      HgcdMatRawDenseRep this heap B right hB →
      ∃ result,
        hgcdMatMulLoop this A B hA hB T scratch C hC i heap = .ok result
  | C, hC, i, heap, hLeft, hRight => by
      by_cases hi : i < 4
      · let out : Fin 4 := ⟨i, hi⟩
        let rowBase : Fin 4 := ⟨2 * (i / 2), by omega⟩
        let rowNext : Fin 4 := ⟨2 * (i / 2) + 1, by omega⟩
        let col : Fin 4 := ⟨i % 2, by omega⟩
        let lowerCol : Fin 4 := ⟨2 + i % 2, by omega⟩
        have htotal := physical.entry C hC i hi heap
        rcases hgcdMatMulEntry_succeeds this (hgcdMatPtr C hC out)
            (hgcdMatPtr A hA rowBase) (hgcdMatPtr B hB col)
            (hgcdMatPtr A hA rowNext) (hgcdMatPtr B hB lowerCol) T scratch
            (hgcdMatLen A hA rowBase) (hgcdMatLen B hB col)
            (hgcdMatLen A hA rowNext) (hgcdMatLen B hB lowerCol) heap
            (left rowBase) (right col) (left rowNext) (right lowerCol) hcfg hp
            htotal (hLeft rowBase) (hRight col) (hLeft rowNext)
            (hRight lowerCol) with ⟨entry, hentry, _⟩
        have hstep := physical.conditional C hC i hi heap
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
        let nextLen := C.len.set i entry.length (by rw [hC.2]; exact hi)
        let next : HgcdMat := { C with len := nextLen }
        have hNext : next.Valid :=
          ⟨hC.1, by simp [next, nextLen, hC.2]⟩
        rcases hgcdMatMulLoop_succeeds this A B hA hB T scratch left right
            hcfg hp physical next hNext (i + 1) entry.heap hLeft1 hRight1 with
          ⟨result, htail⟩
        refine ⟨result, ?_⟩
        rw [hgcdMatMulLoop]
        rw [dif_pos hi]
        dsimp only
        have hentryRaw : hgcdMatMulEntry this
            (hgcdMatPtrRaw C hC out) (hgcdMatPtrRaw A hA rowBase)
            (hgcdMatPtrRaw B hB col) (hgcdMatPtrRaw A hA rowNext)
            (hgcdMatPtrRaw B hB lowerCol) T scratch
            (hgcdMatLenRaw A hA rowBase) (hgcdMatLenRaw B hB col)
            (hgcdMatLenRaw A hA rowNext) (hgcdMatLenRaw B hB lowerCol) heap =
              .ok entry := by
          simpa only [hgcdMatPtrRaw, hgcdMatPtr, hgcdMatLenRaw, hgcdMatLen]
            using hentry
        rw [hentryRaw]
        simpa [next, nextLen, out, rowBase, rowNext, col, lowerCol] using htail
      · refine ⟨⟨heap, C⟩, ?_⟩
        rw [hgcdMatMulLoop, dif_neg hi]
termination_by C _hC i heap _hLeft _hRight => 4 - i
decreasing_by omega

/-- Total semantic execution of the complete generated `_mat_mul`. -/
theorem hgcdMatMul_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (C A B : HgcdMat) (hC : C.Valid) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) (heap : RawHeap)
    (left right : Fin 4 → Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdMatMulLoopTotalWorkspace this A B hA hB T scratch)
    (hLeft : HgcdMatRawDenseRep this heap A left hA)
    (hRight : HgcdMatRawDenseRep this heap B right hB) :
    ∃ result hResult,
      hgcdMatMul this C A B hC hA hB T scratch heap = .ok result ∧
      HgcdMatRawDenseRep this result.heap result.matrix
        (hgcdMatProductEntry left right) hResult := by
  rcases hgcdMatMulLoop_succeeds this A B hA hB T scratch left right hcfg hp
      physical C hC 0 heap hLeft hRight with ⟨result, hrun⟩
  have hsemantic := hgcdMatMul_refines this C A B hC hA hB T scratch heap
    result left right hcfg hp physical.conditional hLeft hRight (by
      simpa [hgcdMatMul] using hrun)
  rcases hsemantic.2.2 with ⟨hResult, hProduct⟩
  exact ⟨result, hResult, by simpa [hgcdMatMul] using hrun, hProduct⟩

/-- Physical facts and the remaining workspace after the quotient update in
the final recursive combine block. -/
structure HgcdRecursiveCombineMatrixAfterQuotient (this : DenseUPolyZp)
    (R : HgcdMat) (hR : R.Valid) (a2 scratch : RawPtr UInt64)
    (heap : RawHeap) (modified : HgcdMatQuotientResult) : Type where
  layout : RawHeap.SameLayout heap modified.heap
  rightPrefix : ∀ i : Fin 4, SameU64Prefix heap modified.heap
    (hgcdMatPtr R hR i) (hgcdMatLen R hR i)
  multiply : HgcdMatMulLoopTotalWorkspace this R modified.matrix hR
    modified.valid a2 scratch

/-- Staged physical safety for the quotient update followed by the complete
matrix multiplication in the final recursive combine block. -/
structure HgcdRecursiveCombineMatrixTotalWorkspace (this : DenseUPolyZp)
    (R S : HgcdMat) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) : Type where
  quotient : HgcdMatApplyQuotientTotalWorkspace this S hS q lenQ T scratch
    heap
  afterQuotient : ∀ modified,
    hgcdMatApplyQuotient this S hS q lenQ T scratch heap = .ok modified →
    HgcdRecursiveCombineMatrixAfterQuotient this R hR a2 scratch heap modified

/-- Total semantic execution of the exact quotient update and full matrix
product selected by `_hgcd_recursive`. -/
theorem hgcdRecursiveCombineMatrix_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) (right entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveCombineMatrixTotalWorkspace this R S hR hS q
      lenQ T a2 scratch heap)
    (hRight : HgcdMatRawDenseRep this heap R right hR)
    (hSRep : HgcdMatRawDenseRep this heap S entries hS)
    (hQ : RawDensePolyRep this heap q lenQ quotient) :
    ∃ result hResult,
      hgcdRecursiveCombineMatrix this M R S hM hR hS q lenQ T a2 scratch
        heap = .ok result ∧
      HgcdMatRawDenseRep this result.heap result.matrix
        (hgcdMatProductEntry right
          (hgcdMatApplyQuotientEntries entries quotient)) hResult := by
  rcases hgcdMatApplyQuotient_succeeds this S hS q lenQ T scratch heap entries
      quotient hcfg hp physical.quotient hQ hSRep with
    ⟨modified, hmodified, hModifiedRep, _⟩
  rcases physical.afterQuotient modified hmodified with
    ⟨hLayout, hRightPrefix, hMultiply⟩
  have hRightModified : HgcdMatRawDenseRep this modified.heap R right hR := by
    intro i
    exact CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix
      this heap modified.heap (hgcdMatPtr R hR i) (hgcdMatLen R hR i)
      (right i) hLayout (hRightPrefix i) (hRight i)
  rcases hgcdMatMul_succeeds this M R modified.matrix hM hR modified.valid a2
      scratch modified.heap right (hgcdMatApplyQuotientEntries entries quotient)
      hcfg hp hMultiply hRightModified hModifiedRep with
    ⟨result, hResult, hmul, hProduct⟩
  refine ⟨result, hResult, ?_, hProduct⟩
  simp [hgcdRecursiveCombineMatrix, hmodified, hmul]

/-- Concrete frames and the optional remaining workspace after the actual
four-call reconstruction in the recursive finish block. -/
structure HgcdRecursiveFinishAfterReconstruct (this : DenseUPolyZp)
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (computeM : Bool) (A B q T0 a2 scratch : RawPtr UInt64) (lenQ : Nat)
    (heap : RawHeap) (reconstructed : HgcdRecursiveReconstructPairResult) :
    Type where
  layout : RawHeap.SameLayout heap reconstructed.heap
  rightPrefix : ∀ i : Fin 4, SameU64Prefix heap reconstructed.heap
    (hgcdMatPtr R hR i) (hgcdMatLen R hR i)
  secondPrefix : ∀ i : Fin 4, SameU64Prefix heap reconstructed.heap
    (hgcdMatPtr S hS i) (hgcdMatLen S hS i)
  quotientPrefix : SameU64Prefix heap reconstructed.heap q lenQ
  combine : computeM = true →
    HgcdRecursiveCombineMatrixTotalWorkspace this R S hR hS q lenQ T0 a2
      scratch reconstructed.heap
  afterCombine : ∀ (_hcompute : computeM = true) combined,
    hgcdRecursiveCombineMatrix this M R S hM hR hS q lenQ T0 a2 scratch
        reconstructed.heap = .ok combined →
    RawHeap.SameLayout reconstructed.heap combined.heap ∧
      SameU64Prefix reconstructed.heap combined.heap A reconstructed.lenA ∧
      SameU64Prefix reconstructed.heap combined.heap B reconstructed.lenB

/-- Non-circular total workspace for the exact generated recursive finish:
the later combine obligations are indexed by the preceding reconstruction. -/
structure HgcdRecursiveFinishTotalWorkspace (this : DenseUPolyZp)
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (computeM : Bool) (A B T0 lowA lowB highA highB q : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift lenQ : Nat)
    (a2 scratch : RawPtr UInt64) (sgnS : Int) (heap : RawHeap) : Type where
  reconstruct : HgcdRecursiveReconstructPairTotalWorkspace this A B T0 lowA
    lowB highA highB scratch lenLowA lenLowB lenHighA lenHighB shift S hS
    sgnS heap
  afterReconstruct : ∀ reconstructed,
    hgcdRecursiveReconstructPair this A B T0 lowA lowB highA highB scratch
        lenLowA lenLowB lenHighA lenHighB shift S hS sgnS heap =
      .ok reconstructed →
    HgcdRecursiveFinishAfterReconstruct this M R S hM hR hS computeM A B q
      T0 a2 scratch lenQ heap reconstructed

/-- Total semantic execution of the exact generated recursive finish block.
Neither reconstruction nor the optional matrix block is assumed to succeed. -/
theorem hgcdRecursiveFinish_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (computeM : Bool) (A B T0 lowA lowB highA highB q : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift lenQ : Nat)
    (a2 scratch : RawPtr UInt64) (sgnR sgnS : Int) (heap : RawHeap)
    (right entries : Fin 4 → Polynomial (ZMod this._p.toNat))
    (quotient polyLowA polyLowB polyHighA polyHighB :
      Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (physical : HgcdRecursiveFinishTotalWorkspace this M R S hM hR hS
      computeM A B T0 lowA lowB highA highB q lenLowA lenLowB lenHighA
      lenHighB shift lenQ a2 scratch sgnS heap)
    (hRRep : HgcdMatRawDenseRep this heap R right hR)
    (hSRep : HgcdMatRawDenseRep this heap S entries hS)
    (hQ : RawDensePolyRep this heap q lenQ quotient)
    (hLowA : RawCanonicalPolySlice this heap lowA lenLowA polyLowA)
    (hLowB : RawCanonicalPolySlice this heap lowB lenLowB polyLowB)
    (hHighA : RawDensePolyRep this heap highA lenHighA polyHighA)
    (hHighB : RawDensePolyRep this heap highB lenHighB polyHighB) :
    ∃ result,
      hgcdRecursiveFinish this M R S hM hR hS computeM A B T0 lowA lowB
          highA highB q lenLowA lenLowB lenHighA lenHighB shift lenQ a2
          scratch sgnR sgnS heap = .ok result ∧
      RawDensePolyRep this result.heap A result.lenA
        (hgcdReconstructedLowA entries polyLowA polyLowB sgnS +
          Polynomial.X ^ shift * polyHighA) ∧
      RawDensePolyRep this result.heap B result.lenB
        (hgcdReconstructedLowB entries polyLowA polyLowB sgnS +
          Polynomial.X ^ shift * polyHighB) ∧
      result.sgn = -(sgnR * sgnS) ∧
      (computeM = true →
        HgcdMatRawDenseRep this result.heap result.matrix
          (hgcdMatProductEntry right
            (hgcdMatApplyQuotientEntries entries quotient)) result.valid) := by
  rcases hgcdRecursiveReconstructPair_succeeds this A B T0 lowA lowB highA
      highB scratch lenLowA lenLowB lenHighA lenHighB shift S hS sgnS heap
      entries polyLowA polyLowB polyHighA polyHighB hcfg hp
      physical.reconstruct hSRep hLowA hLowB hHighA hHighB with
    ⟨reconstructed, hreconstruct, hAReconstructed, hBReconstructed, _⟩
  have hafter := physical.afterReconstruct reconstructed hreconstruct
  have hRReconstructed : HgcdMatRawDenseRep this reconstructed.heap R right
      hR := by
    intro i
    exact CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix
      this heap reconstructed.heap (hgcdMatPtr R hR i) (hgcdMatLen R hR i)
      (right i) hafter.layout (hafter.rightPrefix i) (hRRep i)
  have hSReconstructed : HgcdMatRawDenseRep this reconstructed.heap S entries
      hS := by
    intro i
    exact CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix
      this heap reconstructed.heap (hgcdMatPtr S hS i) (hgcdMatLen S hS i)
      (entries i) hafter.layout (hafter.secondPrefix i) (hSRep i)
  have hQReconstructed : RawDensePolyRep this reconstructed.heap q lenQ
      quotient :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
      heap reconstructed.heap q lenQ quotient hafter.layout
      hafter.quotientPrefix hQ
  cases computeM with
  | false =>
      refine ⟨{
        heap := reconstructed.heap
        matrix := M
        valid := hM
        lenA := reconstructed.lenA
        lenB := reconstructed.lenB
        sgn := -(sgnR * sgnS) }, ?_, hAReconstructed, hBReconstructed,
        rfl, ?_⟩
      · simp [hgcdRecursiveFinish, hreconstruct]
      · intro htrue
        simp at htrue
  | true =>
      rcases hgcdRecursiveCombineMatrix_succeeds this M R S hM hR hS q lenQ
          T0 a2 scratch reconstructed.heap right entries quotient hcfg hp
          (hafter.combine rfl) hRReconstructed hSReconstructed
          hQReconstructed with
        ⟨combined, hCombinedValid, hcombine, hCombinedRep⟩
      rcases hafter.afterCombine rfl combined hcombine with
        ⟨hLayoutCombined, hAPrefix, hBPrefix⟩
      have hAFinal :=
        CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
          reconstructed.heap combined.heap A reconstructed.lenA
          (hgcdReconstructedLowA entries polyLowA polyLowB sgnS +
            Polynomial.X ^ shift * polyHighA)
          hLayoutCombined hAPrefix hAReconstructed
      have hBFinal :=
        CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
          reconstructed.heap combined.heap B reconstructed.lenB
          (hgcdReconstructedLowB entries polyLowA polyLowB sgnS +
            Polynomial.X ^ shift * polyHighB)
          hLayoutCombined hBPrefix hBReconstructed
      have hfinishExists : ∃ result,
          hgcdRecursiveFinish this M R S hM hR hS true A B T0 lowA lowB
              highA highB q lenLowA lenLowB lenHighA lenHighB shift lenQ a2
              scratch sgnR sgnS heap = .ok result := by
        simp only [hgcdRecursiveFinish, hreconstruct, ↓reduceIte]
        split
        next fault hfault => simp [hcombine] at hfault
        next combined' hsuccess => exact ⟨_, rfl⟩
      rcases hfinishExists with ⟨result, hfinish⟩
      rcases hgcdRecursiveFinish_exec this M R S hM hR hS true A B T0
          lowA lowB highA highB q lenLowA lenLowB lenHighA lenHighB shift lenQ
          a2 scratch sgnR sgnS heap result hfinish with
        ⟨reconstructed', hreconstruct', hlenA, hlenB, hsgn, htail⟩
      have hreconstructedEq : reconstructed' = reconstructed :=
        Except.ok.inj (hreconstruct'.symm.trans hreconstruct)
      subst reconstructed'
      simp at htail
      rcases htail with ⟨combined', hcombine', hheap, hmatrix⟩
      have hcombinedEq : combined' = combined :=
        Except.ok.inj (hcombine'.symm.trans hcombine)
      subst combined'
      refine ⟨result, hfinish, ?_, ?_, hsgn, ?_⟩
      · simpa [hheap, hlenA] using hAFinal
      · simpa [hheap, hlenB] using hBFinal
      · intro _
        simpa [hheap, hmatrix] using hCombinedRep

/-- The stabilization workspace is available after the actual iterator;
the remaining finalize workspace is selected only after stabilization. -/
structure HgcdRecursiveIterBranchAfterIter (this : DenseUPolyZp)
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage : RawPtr UInt64) (heap : RawHeap)
    (iter : HgcdIterState) (hIter : iter.matrix.Valid) : Type where
  stabilize : HgcdMatStabilizeWorkspace iter.heap original iter.matrix
    hOriginal hIter stage
  afterStable : ∀ (stable : HgcdMatRestoreResult)
    (hStable : stable.matrix.Valid),
    hgcdMatStabilize original iter.matrix hOriginal hIter stage iter.heap =
        .ok stable →
    HgcdRecursiveIterFinalizeWorkspace original hOriginal a3 b3 stage iter
      hIter stable hStable

/-- Non-circular total physical contract for the exact iterator cutoff arm. -/
structure HgcdRecursiveIterBranchTotalWorkspace (this : DenseUPolyZp)
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat)) : Type where
  inputAPos : 0 < lenInputA
  iterTotal : HgcdIterTotalWorkspaceProvider this (lenInputA / 2) Q W3
    scratch
  loopPhysical : HgcdLoopWorkspaceProvider this (lenInputA / 2) Q W3 scratch
  valid0 : heap.ValidU64Slice
    (hgcdMatPtr original hOriginal (0 : Fin 4)) 1
  valid3 : heap.ValidU64Slice
    (hgcdMatPtr original hOriginal (3 : Fin 4)) 1
  disjoint03 : U64SlicesDisjoint
    (hgcdMatPtr original hOriginal (0 : Fin 4)) 1
    (hgcdMatPtr original hOriginal (3 : Fin 4)) 1
  validA3 : heap.ValidU64Slice a3 lenInputA
  validB3 : heap.ValidU64Slice b3 lenInputB
  a3InputA : U64SlicesDisjoint a3 lenInputA inputA lenInputA
  b3InputB : U64SlicesDisjoint b3 lenInputB inputB lenInputB
  a3InputB : U64SlicesDisjoint a3 lenInputA inputB lenInputB
  b3A3 : U64SlicesDisjoint b3 lenInputB a3 lenInputA
  row0InputA : U64SlicesDisjoint
    (hgcdMatPtr original hOriginal (0 : Fin 4)) 1 inputA lenInputA
  row3InputA : U64SlicesDisjoint
    (hgcdMatPtr original hOriginal (3 : Fin 4)) 1 inputA lenInputA
  row0InputB : U64SlicesDisjoint
    (hgcdMatPtr original hOriginal (0 : Fin 4)) 1 inputB lenInputB
  row3InputB : U64SlicesDisjoint
    (hgcdMatPtr original hOriginal (3 : Fin 4)) 1 inputB lenInputB
  a3Matrix : ∀ i : Fin 4, U64SlicesDisjoint a3 lenInputA
    (hgcdMatPtr original hOriginal i) (identityEntryLen i)
  b3Matrix : ∀ i : Fin 4, U64SlicesDisjoint b3 lenInputB
    (hgcdMatPtr original hOriginal i) (identityEntryLen i)
  matrixValid : ∀ i : Fin 4, heap.ValidU64Slice
    (hgcdMatPtr original hOriginal i) (identityEntryLen i)
  leftRep : RawDensePolyRep this heap inputA lenInputA left
  rightRep : RawDensePolyRep this heap inputB lenInputB right
  afterIter : ∀ (iter : HgcdIterState) (hIter : iter.matrix.Valid),
    dense_upoly_zp__hgcd_iter_ir this original a3 b3 T0 T1 0 inputA
        lenInputA inputB lenInputB Q W3 scratch heap = .ok iter →
    HgcdRecursiveIterBranchAfterIter this original hOriginal a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap iter hIter

/-- Total execution and L2 invariant for the exact generated cutoff iterator
arm, including stabilization and alias-sensitive output normalization. -/
theorem hgcdRecursiveIterBranch_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenInputB < lenInputA)
    (physical : HgcdRecursiveIterBranchTotalWorkspace this original hOriginal
      a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap
      left right) :
    ∃ result finalA finalB entries hResultM,
      hgcdRecursiveIterBranch this original hOriginal a3 b3 inputA inputB
          lenInputA lenInputB Q W3 T0 T1 scratch stage heap = .ok result ∧
      HgcdRecursiveRawInvariant this left right finalA finalB entries true a3
        b3 lenInputA (result.toResult hResultM) := by
  rcases hgcdIter_succeeds this original a3 b3 T0 T1 0 inputA lenInputA
      inputB lenInputB Q W3 scratch heap left right hOriginal hcfg hp
      physical.iterTotal physical.valid0 physical.valid3 physical.disjoint03
      physical.validA3 physical.validB3 physical.a3InputA physical.b3InputB
      physical.a3InputB physical.b3A3 physical.row0InputA physical.row3InputA
      physical.row0InputB physical.row3InputB physical.a3Matrix
      physical.b3Matrix physical.matrixValid physical.leftRep physical.rightRep
      with
    ⟨iter, finalA, finalB, entries, hIter, hiter, hIterInvariant, _, _⟩
  have hafter := physical.afterIter iter hIter hiter
  rcases hgcdMatStabilize_refines this original iter.matrix hOriginal hIter
      stage entries iter.heap hafter.stabilize hIterInvariant.matrixRep with
    ⟨stable, hstable, hStable, hStableRep⟩
  have hfinalize := hafter.afterStable stable hStable hstable
  have hAStable := hgcdMatStabilize_preserves_rawDenseRep this original
    iter.matrix hOriginal hIter stage entries iter.heap hafter.stabilize
    hIterInvariant.matrixRep iter.A iter.lenA finalA hfinalize.stageA
    hfinalize.originalA hIterInvariant.aRep stable hstable
  have hBStable := hgcdMatStabilize_preserves_rawDenseRep this original
    iter.matrix hOriginal hIter stage entries iter.heap hafter.stabilize
    hIterInvariant.matrixRep iter.B iter.lenB finalB hfinalize.stageB
    hfinalize.originalB hIterInvariant.bRep stable hstable
  have hA3Stable : stable.heap.ValidU64Slice a3 iter.lenA :=
    (hAStable.1 a3 iter.lenA).mp hfinalize.validA3
  have hB3Stable : stable.heap.ValidU64Slice b3 iter.lenB :=
    (hBStable.1 b3 iter.lenB).mp hfinalize.validB3
  rcases hgcdRecursiveStoreIterOutputs_refines this a3 b3 iter.A iter.B
      iter.lenA iter.lenB finalA finalB stable.heap hfinalize.pAEq
      hfinalize.pBEq hA3Stable hB3Stable hAStable.2 hBStable.2
      hfinalize.b3PB hfinalize.b3PA hfinalize.a3PA hfinalize.a3PB
      hfinalize.a3B3 with
    ⟨heap1, hstore, _, _⟩
  let result : HgcdRecursiveIterBranchResult := {
    heap := heap1
    matrix := stable.matrix
    lenA := iter.lenA
    lenB := iter.lenB
    sgn := iter.sgn }
  have hbranch : hgcdRecursiveIterBranch this original hOriginal a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap = .ok result := by
    rw [hgcdRecursiveIterBranch]
    split
    next fault hiter' => simp [hiter] at hiter'
    next iter' hiter' =>
      have hiterEq : iter' = iter := Except.ok.inj (hiter'.symm.trans hiter)
      subst iter'
      dsimp only
      split
      next fault hstable' =>
        have hstableFalse :
            hgcdMatStabilize original iter.matrix hOriginal hIter stage
                iter.heap = .error fault := by
          simpa only [Subsingleton.elim (hgcdIter_result_valid this original
            a3 b3 T0 T1 0 inputA lenInputA inputB lenInputB Q W3 scratch heap
            iter hiter') hIter] using hstable'
        rw [hstable] at hstableFalse
        simp at hstableFalse
      next stable' hstable' =>
        have hstableRun :
            hgcdMatStabilize original iter.matrix hOriginal hIter stage
                iter.heap = .ok stable' := by
          simpa only [Subsingleton.elim (hgcdIter_result_valid this original
            a3 b3 T0 T1 0 inputA lenInputA inputB lenInputB Q W3 scratch heap
            iter hiter') hIter] using hstable'
        have hstableEq : stable' = stable :=
          Except.ok.inj (hstableRun.symm.trans hstable)
        subst stable'
        simp [hstore, result]
  have hFinalizeProvider : HgcdRecursiveIterFinalizeWorkspaceProvider this
      original hOriginal a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1
      scratch stage heap := by
    intro iter' hIter' stable' hStable' hiter' hstable'
    have hiterEq : iter' = iter := Except.ok.inj (hiter'.symm.trans hiter)
    subst iter'
    have hhIter : hIter' = hIter := Subsingleton.elim _ _
    subst hIter'
    have hstableEq : stable' = stable :=
      Except.ok.inj (hstable'.symm.trans hstable)
    subst stable'
    have hhStable : hStable' = hStable := Subsingleton.elim _ _
    subst hStable'
    exact hfinalize
  let oldPhysical : HgcdRecursiveDispatchIterWorkspace this original hOriginal
      a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap
      left right := {
    inputAPos := physical.inputAPos
    loopPhysical := physical.loopPhysical
    finalizePhysical := hFinalizeProvider
    valid0 := physical.valid0
    valid3 := physical.valid3
    disjoint03 := physical.disjoint03
    validA3 := physical.validA3
    validB3 := physical.validB3
    a3InputA := physical.a3InputA
    b3InputB := physical.b3InputB
    a3InputB := physical.a3InputB
    b3A3 := physical.b3A3
    row0InputA := physical.row0InputA
    row3InputA := physical.row3InputA
    row0InputB := physical.row0InputB
    row3InputB := physical.row3InputB
    a3Matrix := physical.a3Matrix
    b3Matrix := physical.b3Matrix
    matrixValid := physical.matrixValid
    leftRep := physical.leftRep
    rightRep := physical.rightRep }
  rcases hgcdRecursiveIterBranch_refines this original hOriginal a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap result left
      right hcfg hp (Nat.le_of_lt horder) oldPhysical.inputAPos
      oldPhysical.loopPhysical oldPhysical.finalizePhysical oldPhysical.valid0
      oldPhysical.valid3 oldPhysical.disjoint03 oldPhysical.validA3
      oldPhysical.validB3 oldPhysical.a3InputA oldPhysical.b3InputB
      oldPhysical.a3InputB oldPhysical.b3A3 oldPhysical.row0InputA
      oldPhysical.row3InputA oldPhysical.row0InputB oldPhysical.row3InputB
      oldPhysical.a3Matrix oldPhysical.b3Matrix oldPhysical.matrixValid
      oldPhysical.leftRep oldPhysical.rightRep hbranch with
    ⟨finalA', finalB', entries', hResultM, hInvariant⟩
  exact ⟨result, finalA', finalB', entries', hResultM, hbranch, hInvariant⟩

/-- Total well-founded induction contract at one strictly smaller recursive
call.  Unlike the older conditional refinement contract, it constructs the
actual child execution as well as its semantic invariant. -/
def HgcdRecursiveCallbackSucceedsAt (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)] (bound : Nat)
    (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (WNext scratch : RawPtr UInt64) (heap : RawHeap)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound)
    (left right : Polynomial (ZMod this._p.toNat)) : Prop :=
  ∃ result finalA finalB entries,
    recurse matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
        WNext scratch heap horder hdecrease = .ok result ∧
    HgcdRecursiveRawInvariant this left right finalA finalB entries true a3
      b3 lenInputA result

/-- Total execution of the exact cutoff dispatch.  The small source arm is
the generated iterator block; the large arm is the strictly smaller
well-founded recursive call. -/
theorem hgcdRecursiveDispatchBelow_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound)
    (iterPhysical : HgcdRecursiveIterBranchTotalWorkspace this matrix hMatrix
      a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap
      left right)
    (recursiveSucceeds : HgcdRecursiveCallbackSucceedsAt this bound recurse
      matrix hMatrix a3 b3 inputA inputB lenInputA lenInputB WNext scratch
      heap horder hdecrease left right) :
    ∃ result finalA finalB entries,
      hgcdRecursiveDispatchBelow this bound recurse matrix hMatrix a3 b3
          inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage WNext
          heap horder hdecrease = .ok result ∧
      HgcdRecursiveRawInvariant this left right finalA finalB entries true a3
        b3 lenInputA result := by
  by_cases hsmall : lenInputA < hgcdRecursiveCutoff
  · rcases hgcdRecursiveIterBranch_succeeds this matrix hMatrix a3 b3 inputA
        inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap left right
        hcfg hp horder iterPhysical with
      ⟨iter, finalA, finalB, entries, hIterValid, hiter, hInvariant⟩
    let hGeneratedValid := hgcdRecursiveIterBranch_result_valid this matrix
      hMatrix a3 b3 inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch
      stage heap iter hiter
    let result := iter.toResult hGeneratedValid
    refine ⟨result, finalA, finalB, entries, ?_, ?_⟩
    · rw [hgcdRecursiveDispatchBelow, if_pos hsmall]
      split
      next fault hiter' => simp [hiter] at hiter'
      next iter' hiter' =>
        have hiterEq : iter' = iter := Except.ok.inj (hiter'.symm.trans hiter)
        subst iter'
        apply congrArg Except.ok
        apply HgcdRecursiveResult.ext_value
        rfl
    · simpa [result] using hInvariant
  · rcases recursiveSucceeds with
      ⟨result, finalA, finalB, entries, hrun, hInvariant⟩
    exact ⟨result, finalA, finalB, entries, by
      simpa [hgcdRecursiveDispatchBelow, hsmall] using hrun, hInvariant⟩

/-- Staged total safety for the first child dispatch and the immediately
following four-call paired reconstruction of a non-base invocation. -/
structure HgcdRecursiveFirstCallTotalWorkspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (package : HgcdRecursiveNonBasePackage this bound recurse a b W scratch
      lenA lenB heap) : Type where
  iter :
    let ws := hgcdRecursiveWorkspace W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    HgcdRecursiveIterBranchTotalWorkspace this ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 heap
      package.inputHighA package.inputHighB
  reconstruct : ∀ (first : HgcdRecursiveResult)
    (hchildOrder :
      (hgcdRecursiveHighInput a b lenA lenB).lenB0 <
        (hgcdRecursiveHighInput a b lenA lenB).lenA0)
    (hchildDecrease :
      (hgcdRecursiveHighInput a b lenA lenB).lenA0 < bound),
    let ws := hgcdRecursiveWorkspace W lenA
    hgcdRecursiveDispatchBelow this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3
        (hgcdRecursiveHighInput a b lenA lenB).a0
        (hgcdRecursiveHighInput a b lenA lenB).b0
        (hgcdRecursiveHighInput a b lenA lenB).lenA0
        (hgcdRecursiveHighInput a b lenA lenB).lenB0 ws.q ws.W3 ws.T0 ws.T1
        scratch ws.a2 ws.next heap hchildOrder hchildDecrease = .ok first →
    HgcdRecursiveReconstructPairTotalWorkspace this ws.a2 ws.b2 ws.T0 a b
      ws.a3 ws.b3 scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap

/-- Total first-child dispatch and exact first paired reconstruction. -/
theorem hgcdRecursiveFirstCall_succeeds (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (bound : Nat) (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (hbase : ¬ lenB < lenA / 2 + 1)
    (package : HgcdRecursiveNonBasePackage this bound recurse a b W scratch
      lenA lenB heap)
    (physical : HgcdRecursiveFirstCallTotalWorkspace this bound recurse a b W
      scratch lenA lenB heap package)
    (recursiveSucceeds :
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      HgcdRecursiveCallbackSucceedsAt this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.next scratch heap hchildOrder hchildDecrease
        package.inputHighA package.inputHighB) :
    ∃ first reconstructed outputHighA outputHighB entries,
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∃ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      hgcdRecursiveDispatchBelow this bound recurse ws.R
          (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
          high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
          heap hchildOrder hchildDecrease = .ok first ∧
      hgcdRecursiveReconstructPair this ws.a2 ws.b2 ws.T0 a b ws.a3 ws.b3
          scratch (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2))
          first.lenA first.lenB (lenA / 2) first.matrix first.valid first.sgn
          first.heap = .ok reconstructed ∧
      HgcdRecursiveRawInvariant this package.inputHighA package.inputHighB
        outputHighA outputHighB entries true ws.a3 ws.b3 high.lenA0 first ∧
      RawDensePolyRep this reconstructed.heap ws.a2 reconstructed.lenA
        (hgcdReconstructedLowA entries package.lowPolyA package.lowPolyB
            first.sgn + Polynomial.X ^ (lenA / 2) * outputHighA) ∧
      RawDensePolyRep this reconstructed.heap ws.b2 reconstructed.lenB
        (hgcdReconstructedLowB entries package.lowPolyA package.lowPolyB
            first.sgn + Polynomial.X ^ (lenA / 2) * outputHighB) ∧
      HgcdMatRawDenseRep this reconstructed.heap first.matrix entries
        first.valid ∧
      HgcdFirstReconstructionInvariant lenA first reconstructed := by
  let ws := hgcdRecursiveWorkspace W lenA
  let high := hgcdRecursiveHighInput a b lenA lenB
  have hchildOrder : high.lenB0 < high.lenA0 :=
    hgcdRecursiveHighInput_order a b lenA lenB horder
  have hchildDecrease : high.lenA0 < bound := by
    rw [← hbound]
    exact hgcdRecursiveHighInput_len_lt a b lenA lenB horder (by omega)
  rcases hgcdRecursiveDispatchBelow_succeeds this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap package.inputHighA package.inputHighB hcfg hp hchildOrder
      hchildDecrease physical.iter (recursiveSucceeds hchildOrder
        hchildDecrease) with
    ⟨first, outputHighA, outputHighB, entries, hfirst, hFirstInvariant⟩
  have hframe := package.workspace.frame hchildOrder hchildDecrease first
    hfirst
  have hLowAFirst : RawDensePolyRep this first.heap a
      (Nat.min lenA (lenA / 2)) package.lowPolyA :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
      heap first.heap a (Nat.min lenA (lenA / 2)) package.lowPolyA
      hframe.layout hframe.lowAPrefix package.workspace.lowA
  have hLowBFirst : RawDensePolyRep this first.heap b
      (Nat.min lenB (lenA / 2)) package.lowPolyB :=
    CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix this
      heap first.heap b (Nat.min lenB (lenA / 2)) package.lowPolyB
      hframe.layout hframe.lowBPrefix package.workspace.lowB
  rcases hgcdRecursiveReconstructPair_succeeds this ws.a2 ws.b2 ws.T0 a b
      ws.a3 ws.b3 scratch (Nat.min lenA (lenA / 2))
      (Nat.min lenB (lenA / 2)) first.lenA first.lenB (lenA / 2)
      first.matrix first.valid first.sgn first.heap entries package.lowPolyA
      package.lowPolyB outputHighA outputHighB hcfg hp
      (physical.reconstruct first hchildOrder hchildDecrease hfirst)
      (hFirstInvariant.matrixSemantics rfl).1
      ⟨hLowAFirst.1, hLowAFirst.2.1, hLowAFirst.2.2.1⟩
      ⟨hLowBFirst.1, hLowBFirst.2.1, hLowBFirst.2.2.1⟩
      hFirstInvariant.aRep hFirstInvariant.bRep with
    ⟨reconstructed, hreconstruct, hAReconstructed, hBReconstructed,
      hMatrixReconstructed⟩
  have hactual : HgcdFirstDispatchResult this bound recurse a b W scratch
      lenA lenB heap first := by
    exact ⟨hchildOrder, hchildDecrease, by simpa [ws, high] using hfirst⟩
  have hOldReconstruct := package.workspace.reconstruct first hactual
  have hReconstructionInvariant :=
    hgcdRecursiveFirstReconstruct_invariant_of_execution this ws.a2 ws.b2
      ws.T0 a b ws.a3 ws.b3 scratch lenA lenB first reconstructed entries
      package.lowPolyA package.lowPolyB outputHighA outputHighB hcfg hp horder
      (by simpa [high, hgcdRecursiveHighInput] using
        hFirstInvariant.lengths rfl)
      hOldReconstruct (hFirstInvariant.matrixSemantics rfl).1 hLowAFirst
      hLowBFirst hFirstInvariant.aRep hFirstInvariant.bRep hreconstruct
  exact ⟨first, reconstructed, outputHighA, outputHighB, entries,
    hchildOrder, hchildDecrease, by simpa [ws, high] using hfirst,
    by simpa [ws] using hreconstruct, by simpa [ws, high] using hFirstInvariant,
    by simpa [ws] using hAReconstructed, by simpa [ws] using hBReconstructed,
    by simpa [ws] using hMatrixReconstructed, hReconstructionInvariant⟩

/-- Physical divrem storage available at every represented state reached by
the source HGCD-GCD loop.  The provider supplies only allocation and aliasing
facts; quotient and remainder semantics still come from the actual raw call. -/
def GcdHgcdLoopDivremWorkspaceProvider (this : DenseUPolyZp)
    (G J Q R : RawPtr UInt64) (W3 : RawPtr Word3) (capacity : Nat) : Prop :=
  ∀ (heap : RawHeap) (lenG lenJ : Nat)
    (left right : Polynomial (ZMod this._p.toNat)),
    lenG ≤ capacity → lenJ ≤ capacity →
    RawDensePolyRep this heap G lenG left →
    RawDensePolyRep this heap J lenJ right →
    EuclidWorkspace heap Q W3 G J R capacity

/-- Physical local-vector storage for the exact `_gcd_euclid` call selected
by the small branch.  It is conditional on the represented operands and does
not state the call's result. -/
def GcdHgcdLoopEuclidWorkspaceProvider (this : DenseUPolyZp)
    (G J R euclidA euclidB euclidQ euclidR : RawPtr UInt64)
    (euclidW3 : RawPtr Word3) (capacity : Nat) : Prop :=
  ∀ (heap : RawHeap) (lenJ lenR : Nat)
    (divisor remainder : Polynomial (ZMod this._p.toNat)),
    lenJ ≤ capacity → lenR ≤ capacity →
    RawDensePolyRep this heap J lenJ divisor →
    RawDensePolyRep this heap R lenR remainder →
    GcdEuclidRawWorkspace heap G J R euclidA euclidB euclidQ euclidR
      euclidW3 capacity

/-- Exact-size physical contract for one dynamically reached loop-head
division.  Unlike the fixed-capacity convenience provider above, this record
remains usable after an HGCD call returns new descriptor lengths. -/
structure GcdHgcdDivremWorkspace (heap : RawHeap)
    (G J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenG lenJ : Nat) : Prop where
  validQ : heap.ValidU64Slice Q (lenG - (lenJ - 1))
  validR : heap.ValidU64Slice R (Nat.min lenG (lenJ - 1))
  validW3 : heap.ValidWord3Slice W3 lenG
  validGForDivisor : heap.ValidU64Slice G lenJ
  quotientCapacity : lenG - (lenJ - 1) < limbBase
  r_g : R.region ≠ G.region
  w3_g : W3.region ≠ G.region
  w3_j : W3.region ≠ J.region
  q_j : Q.region ≠ J.region
  q_w3 : Q.region ≠ W3.region
  r_w3 : R.region ≠ W3.region
  r_q : R.region ≠ Q.region
  r_j : R.region ≠ J.region
  g_j : G.region ≠ J.region

/-- Dynamic physical safety for every represented loop state.  It carries no
semantic output and is therefore stable under the lengths returned by the
actual checked HGCD execution. -/
def GcdHgcdDynamicDivremWorkspaceProvider (this : DenseUPolyZp)
    (G J Q R : RawPtr UInt64) (W3 : RawPtr Word3) : Prop :=
  ∀ (heap : RawHeap) (lenG lenJ : Nat)
    (left right : Polynomial (ZMod this._p.toNat)),
    RawDensePolyRep this heap G lenG left →
    RawDensePolyRep this heap J lenJ right →
    GcdHgcdDivremWorkspace heap G J Q R W3 lenG lenJ

/-- Dynamic local-vector safety for the small Euclid branch. -/
def GcdHgcdDynamicEuclidWorkspaceProvider (this : DenseUPolyZp)
    (G J R euclidA euclidB euclidQ euclidR : RawPtr UInt64)
    (euclidW3 : RawPtr Word3) : Prop :=
  ∀ (heap : RawHeap) (lenJ lenR : Nat)
    (divisor remainder : Polynomial (ZMod this._p.toNat)),
    RawDensePolyRep this heap J lenJ divisor →
    RawDensePolyRep this heap R lenR remainder →
    ∃ capacity,
      GcdEuclidRawWorkspace heap G J R euclidA euclidB euclidQ euclidR
          euclidW3 capacity ∧
      lenJ ≤ capacity ∧ lenR ≤ capacity ∧ capacity < limbBase

/-- Concrete recursive safety package for the exact large-branch HGCD call
at one dynamically reached divisor/remainder state. -/
structure GcdHgcdDynamicHgcdWorkspace (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid) (G J R : RawPtr UInt64)
    (W scratch : RawPtr UInt64) (heap : RawHeap) (lenJ lenR : Nat)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (horder : lenR < lenJ) : Type where
  invocation : HgcdRecursiveInvocationWorkspace this
    (hgcdRecursiveCallChecked this) lenJ M hM false G J J R lenJ lenR W
    scratch heap divisor remainder hcfg hp rfl horder
  descendants : HgcdRecursiveInvocationWorkspaceProviderBelow this
    (hgcdRecursiveCallChecked this) hcfg hp lenJ

/-- Dynamic physical provider for every represented large-branch state.  It
contains recursive allocation safety but neither an HGCD result nor a length
claim. -/
def GcdHgcdDynamicHgcdWorkspaceProvider (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid) (G J R : RawPtr UInt64)
    (W scratch : RawPtr UInt64) : Type :=
  ∀ (heap : RawHeap) (lenJ lenR : Nat)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (horder : lenR < lenJ),
    RawDensePolyRep this heap J lenJ divisor →
    RawDensePolyRep this heap R lenR remainder →
    GcdHgcdDynamicHgcdWorkspace this hcfg hp M hM G J R W scratch heap lenJ
      lenR divisor remainder horder

/-- The exact generated base arm cannot fault when its physical workspace is
valid.  This is the first total-execution component needed to upgrade checked
HGCD from conditional semantic correctness to total raw refinement. -/
theorem hgcdRecursiveBase_succeeds (this : DenseUPolyZp)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (hp : 1 < this._p.toNat)
    (workspace : HgcdRecursiveBaseCallWorkspace this M hM A B a b lenA lenB
      heap left right) :
    ∃ result, hgcdRecursiveBase M computeM A B a b lenA lenB heap =
      .ok result := by
  cases computeM with
  | false =>
      rcases copyU64_refines_rawDense this heap A a lenA left
          workspace.validA workspace.Aa workspace.leftRep with
        ⟨heap1, hcopyA, hlayout1, hARep⟩
      have hright1 :=
        CLPoly.Impl.StrictHGCDRawRefinement.rawDensePolyRep_of_same_prefix
          this heap heap1 b lenB right hlayout1
        (copyU64_preserves_prefix heap heap1 A a b lenA lenB
          workspace.validA workspace.leftRep.1 workspace.Ab hcopyA)
        workspace.rightRep
      have hB1 : heap1.ValidU64Slice B lenB :=
        (hlayout1 B lenB).mp workspace.validB
      rcases copyU64_refines_rawDense this heap1 B b lenB right hB1
          workspace.Bb hright1 with ⟨heap2, hcopyB, _, _⟩
      exact ⟨⟨heap2, M, lenA, lenB, 1⟩, by
        simp [hgcdRecursiveBase, hcopyA, hcopyB]⟩
  | true =>
      rcases hgcdIterInit_refines this M A B A A 0 a lenA b lenB heap left
          right hM hp workspace.valid0 workspace.valid3 workspace.disjoint03
          workspace.validA workspace.validB workspace.Aa workspace.Bb
          workspace.Ab workspace.BA workspace.row0a workspace.row3a
          workspace.row0b workspace.row3b workspace.aMatrix workspace.bMatrix
          workspace.matrixValid workspace.leftRep workspace.rightRep with
        ⟨initial, hinit, _⟩
      refine ⟨initial.toRecursiveBaseResult, ?_⟩
      rw [hgcdRecursiveBase_true_eq_init, hinit]
      rfl

/-- Dynamic form of the loop-head divrem refinement. -/
theorem gcdHgcdLoop_dynamic_divrem_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (G J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenG lenJ : Nat) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (physical : GcdHgcdDynamicDivremWorkspaceProvider this G J Q R W3)
    (hlenJPos : 0 < lenJ)
    (hleft : RawDensePolyRep this heap G lenG left)
    (hright : RawDensePolyRep this heap J lenJ right)
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenQ lenR quotient remainder,
      Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q R G lenG J
          lenJ W3 heap = .ok (heap', lenQ, lenR) ∧
      RawDensePolyRep this heap' Q lenQ quotient ∧
      RawDensePolyRep this heap' J lenJ right ∧
      RawDensePolyRep this heap' R lenR remainder ∧
      normalize (EuclideanDomain.gcd left right) =
        normalize (EuclideanDomain.gcd right remainder) ∧
      RawHeap.SameLayout heap heap' ∧ lenR < lenJ := by
  have ws := physical heap lenG lenJ left right hleft hright
  rcases polyDivrem_next_state this Q R G J lenG lenJ W3 heap left right
      hlenJPos hleft hright ws.validQ ws.validR ws.validW3
      ws.quotientCapacity ws.r_g ws.w3_g ws.w3_j ws.q_j ws.q_w3 ws.r_w3
      ws.r_q ws.r_j hcfg with
    ⟨heap', lenQ, lenR, quotient, remainder, hrun, hQRep, hJRep, hRRep,
      _, hgcd, hlayout, _, _, hlt⟩
  exact ⟨heap', lenQ, lenR, quotient, remainder, hrun, hQRep, hJRep,
    hRRep, hgcd, hlayout, hlt⟩

/-- Dynamic form of the exact small raw-Euclid branch refinement. -/
theorem gcdHgcdLoop_dynamic_euclid_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (G J R euclidA euclidB euclidQ euclidR : RawPtr UInt64)
    (euclidW3 : RawPtr Word3) (lenJ lenR : Nat) (heap : RawHeap)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (physical : GcdHgcdDynamicEuclidWorkspaceProvider this G J R euclidA
      euclidB euclidQ euclidR euclidW3)
    (hdivisor : RawDensePolyRep this heap J lenJ divisor)
    (hremainder : RawDensePolyRep this heap R lenR remainder)
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenG result,
      strictGcdEuclidRaw this G J R euclidA euclidB euclidQ euclidR
          euclidW3 lenJ lenR heap = .ok ⟨heap', lenG⟩ ∧
      RawDensePolyRep this heap' G lenG result ∧
      normalize (EuclideanDomain.gcd divisor remainder) = normalize result ∧
      RawHeap.SameLayout heap heap' := by
  rcases physical heap lenJ lenR divisor remainder hdivisor hremainder with
    ⟨capacity, hworkspace, hlenJ, hlenR, hcapacity⟩
  exact strictGcdEuclidRaw_refines this G J R euclidA euclidB euclidQ euclidR
    euclidW3 lenJ lenR capacity heap divisor remainder hworkspace hlenJ hlenR
    hcapacity hdivisor hremainder hcfg

/-- The actual loop-head divrem call produces the next represented pair and
preserves its normalized gcd. -/
theorem gcdHgcdLoop_divrem_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (G J Q R : RawPtr UInt64) (W3 : RawPtr Word3)
    (capacity lenG lenJ : Nat) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (physical : GcdHgcdLoopDivremWorkspaceProvider this G J Q R W3 capacity)
    (hlenG : lenG ≤ capacity) (hlenJ : lenJ ≤ capacity)
    (hcapacity : capacity < limbBase) (hlenJPos : 0 < lenJ)
    (hleft : RawDensePolyRep this heap G lenG left)
    (hright : RawDensePolyRep this heap J lenJ right)
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenQ lenR quotient remainder,
      Generated.StrictDivrem.dense_upoly_zp__poly_divrem_ir this Q R G lenG J
          lenJ W3 heap = .ok (heap', lenQ, lenR) ∧
      RawDensePolyRep this heap' Q lenQ quotient ∧
      RawDensePolyRep this heap' J lenJ right ∧
      RawDensePolyRep this heap' R lenR remainder ∧
      normalize (EuclideanDomain.gcd left right) =
        normalize (EuclideanDomain.gcd right remainder) ∧
      RawHeap.SameLayout heap heap' ∧ lenR < lenJ := by
  have ws := physical heap lenG lenJ left right hlenG hlenJ hleft hright
  have hqBound : lenG - (lenJ - 1) ≤ capacity :=
    (Nat.sub_le lenG (lenJ - 1)).trans hlenG
  have hrBound : Nat.min lenG (lenJ - 1) ≤ capacity :=
    (Nat.min_le_left lenG (lenJ - 1)).trans hlenG
  rcases polyDivrem_next_state this Q R G J lenG lenJ W3 heap left right
      hlenJPos hleft hright
      (heap.validU64Slice_mono Q capacity _ ws.validQ hqBound)
      (heap.validU64Slice_mono R capacity _ ws.validR hrBound)
      (heap.validU64Slice_mono (RawPtr.reinterpret W3) (3 * capacity)
        (3 * lenG) ws.validW3 (by omega))
      (hqBound.trans_lt hcapacity) ws.regions.a_r.symm ws.regions.w3_a
      ws.regions.w3_b ws.regions.q_b ws.regions.q_w3
      ws.regions.w3_r.symm ws.regions.q_r.symm ws.regions.b_r.symm hcfg with
    ⟨heap', lenQ, lenR, quotient, remainder, hrun, hQRep, hJRep, hRRep,
      _, hgcd, hlayout, _, _, hlt⟩
  exact ⟨heap', lenQ, lenR, quotient, remainder, hrun, hQRep, hJRep,
    hRRep, hgcd, hlayout, hlt⟩

/-- The small source branch is the already-refined exact raw Euclid helper on
the divisor/remainder pair produced by the loop-head division. -/
theorem gcdHgcdLoop_euclid_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (G J R euclidA euclidB euclidQ euclidR : RawPtr UInt64)
    (euclidW3 : RawPtr Word3) (capacity lenJ lenR : Nat)
    (heap : RawHeap)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (physical : GcdHgcdLoopEuclidWorkspaceProvider this G J R euclidA
      euclidB euclidQ euclidR euclidW3 capacity)
    (hlenJ : lenJ ≤ capacity) (hlenR : lenR ≤ capacity)
    (hcapacity : capacity < limbBase)
    (hdivisor : RawDensePolyRep this heap J lenJ divisor)
    (hremainder : RawDensePolyRep this heap R lenR remainder)
    (hcfg : DensePreinvConfigured this) :
    ∃ heap' lenG result,
      strictGcdEuclidRaw this G J R euclidA euclidB euclidQ euclidR
          euclidW3 lenJ lenR heap = .ok ⟨heap', lenG⟩ ∧
      RawDensePolyRep this heap' G lenG result ∧
      normalize (EuclideanDomain.gcd divisor remainder) = normalize result ∧
      RawHeap.SameLayout heap heap' := by
  exact strictGcdEuclidRaw_refines this G J R euclidA euclidB euclidQ euclidR
    euclidW3 lenJ lenR capacity heap divisor remainder
    (physical heap lenJ lenR divisor remainder hlenJ hlenR hdivisor
      hremainder)
    hlenJ hlenR hcapacity hdivisor hremainder hcfg

/-- The zero-remainder source branch copies the current divisor into `G`; a
zero-length represented remainder is the zero polynomial, so this divisor is
the normalized gcd of the current pair. -/
theorem gcdHgcdLoop_zero_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (G J R : RawPtr UInt64) (lenJ : Nat) (heap : RawHeap)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (hG : heap.ValidU64Slice G lenJ) (hGJ : G.region ≠ J.region)
    (hdivisor : RawDensePolyRep this heap J lenJ divisor)
    (hremainder : RawDensePolyRep this heap R 0 remainder) :
    ∃ heap', heap.copyU64 G J lenJ = .ok heap' ∧
      RawDensePolyRep this heap' G lenJ divisor ∧
      normalize (EuclideanDomain.gcd divisor remainder) = normalize divisor ∧
      RawHeap.SameLayout heap heap' := by
  have hremainderZero : remainder = 0 :=
    slicePolyRep_zero_length heap R this._p.toNat remainder
      hremainder.2.2.1
  subst remainder
  rcases copyU64_refines_rawDense_of_region_ne this heap G J lenJ divisor hG
      hGJ hdivisor with ⟨heap', hcopy, hlayout, hresult⟩
  exact ⟨heap', hcopy, hresult, by rw [EuclideanDomain.gcd_zero_right],
    hlayout⟩

/-- The actual result of a successful large-branch checked HGCD call strictly
decreases `_gcd_hgcd`'s source loop measure.  The bound is obtained from the
concrete well-founded HGCD invariant and the source cutoff, rather than from a
separate execution. -/
theorem hgcdRecursiveCallChecked_lenB_lt_of_large
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (left right : Polynomial (ZMod this._p.toNat))
    (horder : lenB < lenA)
    (hlarge : ¬ lenA < hgcdRecursiveCutoff)
    (workspace : HgcdRecursiveInvocationWorkspace this
      (hgcdRecursiveCallChecked this) lenA M hM false A B a b lenA lenB W
      scratch heap left right hcfg hp rfl horder)
    (provider : HgcdRecursiveInvocationWorkspaceProviderBelow this
      (hgcdRecursiveCallChecked this) hcfg hp lenA)
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveCallChecked this M hM false A B a b lenA lenB W
      scratch heap = .ok result) :
    result.lenB < lenA := by
  rcases hgcdRecursiveCallChecked_rawInvariant_wf this hcfg hp M hM false A B
      a b lenA lenB W scratch heap left right horder workspace provider result
      hrun with ⟨finalA, finalB, entries, hinvariant⟩
  have hhalf : lenA / 2 + 1 < lenA := by
    simp only [hgcdRecursiveCutoff] at hlarge
    omega
  exact hinvariant.stopped.trans hhalf

/-- One successful large source branch simultaneously preserves the raw
polynomial pair's normalized gcd and decreases the exact loop length. -/
theorem hgcdRecursiveCallChecked_gcd_step_of_large
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid)
    (G J R : RawPtr UInt64) (lenJ lenR : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (horder : lenR < lenJ)
    (hlarge : ¬ lenJ < hgcdRecursiveCutoff)
    (workspace : HgcdRecursiveInvocationWorkspace this
      (hgcdRecursiveCallChecked this) lenJ M hM false G J J R lenJ lenR W
      scratch heap divisor remainder hcfg hp rfl horder)
    (provider : HgcdRecursiveInvocationWorkspaceProviderBelow this
      (hgcdRecursiveCallChecked this) hcfg hp lenJ)
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveCallChecked this M hM false G J J R lenJ lenR W
      scratch heap = .ok result) :
    ∃ (finalG finalJ : Polynomial (ZMod this._p.toNat)),
      RawDensePolyRep this result.heap G result.lenA finalG ∧
      RawDensePolyRep this result.heap J result.lenB finalJ ∧
      normalize (EuclideanDomain.gcd divisor remainder) =
        normalize (EuclideanDomain.gcd finalG finalJ) ∧
      result.lenB < lenJ := by
  rcases hgcdRecursiveCallChecked_rawInvariant_wf this hcfg hp M hM false G J
      J R lenJ lenR W scratch heap divisor remainder horder workspace provider
      result hrun with ⟨finalG, finalJ, entries, hinvariant⟩
  refine ⟨finalG, finalJ, hinvariant.aRep, hinvariant.bRep,
    hinvariant.gcdPreserved, ?_⟩
  have hhalf : lenJ / 2 + 1 < lenJ := by
    simp only [hgcdRecursiveCutoff] at hlarge
    omega
  exact hinvariant.stopped.trans hhalf

/-- Reachable large-branch termination and semantics are derived from the
dynamic physical provider and the one actual checked HGCD result. -/
theorem gcdHgcdLoop_dynamic_hgcd_step (this : DenseUPolyZp)
    [Fact (Nat.Prime this._p.toNat)]
    (hcfg : DensePreinvConfigured this) (hp : 1 < this._p.toNat)
    (M : HgcdMat) (hM : M.Valid) (G J R : RawPtr UInt64)
    (W scratch : RawPtr UInt64)
    (physical : GcdHgcdDynamicHgcdWorkspaceProvider this hcfg hp M hM G J R
      W scratch)
    (heap : RawHeap) (lenJ lenR : Nat)
    (divisor remainder : Polynomial (ZMod this._p.toNat))
    (horder : lenR < lenJ) (hlarge : ¬ lenJ < hgcdRecursiveCutoff)
    (hdivisor : RawDensePolyRep this heap J lenJ divisor)
    (hremainder : RawDensePolyRep this heap R lenR remainder)
    (result : HgcdRecursiveResult)
    (hrun : hgcdRecursiveCallChecked this M hM false G J J R lenJ lenR W
      scratch heap = .ok result) :
    ∃ (finalG finalJ : Polynomial (ZMod this._p.toNat)),
      RawDensePolyRep this result.heap G result.lenA finalG ∧
      RawDensePolyRep this result.heap J result.lenB finalJ ∧
      normalize (EuclideanDomain.gcd divisor remainder) =
        normalize (EuclideanDomain.gcd finalG finalJ) ∧
      result.lenB < lenJ := by
  let node := physical heap lenJ lenR divisor remainder horder hdivisor
    hremainder
  exact hgcdRecursiveCallChecked_gcd_step_of_large this hcfg hp M hM G J R
    lenJ lenR W scratch heap divisor remainder horder hlarge node.invocation
    node.descendants result hrun

end CLPoly.Impl.StrictGCDHGCDRefinement
