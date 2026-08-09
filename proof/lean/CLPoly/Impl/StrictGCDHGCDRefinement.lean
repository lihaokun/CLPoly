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
