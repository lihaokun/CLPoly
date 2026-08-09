import CLPoly.Generated.StrictMul
import CLPoly.Generated.StrictPolyAddSub
import CLPoly.Generated.StrictDivrem

set_option autoImplicit false

namespace Generated.StrictHGCD

open Generated.StrictMul
open Generated.StrictPolyAddSub
open Generated.StrictDivrem

/-- Raw lowering of the four field updates in C++ `_mat_one`.  The matrix
descriptor itself is a value, while its polynomial entries live in RawHeap. -/
def dense_upoly_zp__mat_one_ir (M : HgcdMat) (heap : RawHeap) :
    RawExec (RawHeap × HgcdMat) :=
  if hvalid : M.poly.size = 4 ∧ M.len.size = 4 then
    let p0 := M.poly[0]'(by omega)
    let p3 := M.poly[3]'(by omega)
    match heap.writeU64 p0 0 1 with
    | .error fault => .error fault
    | .ok heap1 =>
      match heap1.writeU64 p3 0 1 with
      | .error fault => .error fault
      | .ok heap2 =>
        .ok (heap2, { M with len := #[1, 0, 0, 1] })
  else
    .error .assertionFailure

/-- Reference-eliminated state returned by `_mat_row_update`.  `T` and `t`
are the values of the two C++ pointer variables after its swaps. -/
structure MatRowUpdateResult where
  heap : RawHeap
  matrix : HgcdMat
  T : RawPtr UInt64
  lenT : Nat
  t : RawPtr UInt64

/-- Raw lowering of C++ `_mat_row_update`.  The `Fin 4` indices are the
checked form of its two matrix indices; all multiplication, addition and
pointer swaps remain explicit. -/
def dense_upoly_zp__mat_row_update_ir (this : DenseUPolyZp)
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) : RawExec MatRowUpdateResult :=
  if hvalid : M.poly.size = 4 ∧ M.len.size = 4 then
    let p0 := M.poly[i0.val]'(by rw [hvalid.1]; exact i0.isLt)
    let p1 := M.poly[i1.val]'(by rw [hvalid.1]; exact i1.isLt)
    let l0 := M.len[i0.val]'(by rw [hvalid.2]; exact i0.isLt)
    let l1 := M.len[i1.val]'(by rw [hvalid.2]; exact i1.isLt)
    if lenQ ≠ 0 ∧ l0 ≠ 0 then
      let left := if lenQ ≥ l0 then Q else p0
      let leftLen := if lenQ ≥ l0 then lenQ else l0
      let right := if lenQ ≥ l0 then p0 else Q
      let rightLen := if lenQ ≥ l0 then l0 else lenQ
      match dense_upoly_zp__mul_ir this T left leftLen right rightLen scratch heap with
      | .error fault => .error fault
      | .ok heap1 =>
        let productLen := lenQ + l0 - 1
        match dense_upoly_zp__poly_add_ir this t p1 l1 T productLen heap1 with
        | .error fault => .error fault
        | .ok (heap2, sumLen) =>
          let poly' := (M.poly.set i1.val p0 (by omega)).set i0.val t (by simp; omega)
          let len' := (M.len.set i1.val l0 (by omega)).set i0.val sumLen (by simp; omega)
          .ok (MatRowUpdateResult.mk heap2
            ({ poly := poly', len := len' } : HgcdMat) T productLen p1)
    else
      let poly' := (M.poly.set i1.val p0 (by omega)).set i0.val p1 (by simp; omega)
      let len' := (M.len.set i1.val l0 (by omega)).set i0.val l1 (by simp; omega)
      .ok (MatRowUpdateResult.mk heap
        ({ poly := poly', len := len' } : HgcdMat) T lenT t)
  else
    .error .assertionFailure

/-- Reference-eliminated state of C++ `_hgcd_iter`.  The pointer fields are
the current values of `*pA`, `*pB`, `*pT`, and `*pt`. -/
structure HgcdIterState where
  heap : RawHeap
  matrix : HgcdMat
  A : RawPtr UInt64
  lenA : Nat
  B : RawPtr UInt64
  lenB : Nat
  T : RawPtr UInt64
  lenT : Nat
  t : RawPtr UInt64
  sgn : Int

/-- Exact well-founded lowering of the C++ `_hgcd_iter` while-loop.  One
iteration performs the real divrem, pointer rotation, and both matrix-row
updates before recurring on the strictly shorter remainder. -/
def hgcdIterLoop (this : DenseUPolyZp) (m : Nat) (Q : RawPtr UInt64)
    (W3 : RawPtr Word3) (scratch : RawPtr UInt64) :
    (state : HgcdIterState) → RawExec HgcdIterState
  | state =>
    if state.lenB ≥ m + 1 then
      match hdiv : dense_upoly_zp__poly_divrem_ir this Q state.T state.A
          state.lenA state.B state.lenB W3 state.heap with
      | .error fault => .error fault
      | .ok (heap1, lenQ, lenR) =>
        let rotatedA := state.B
        let rotatedB := state.T
        let rotatedT := state.A
        match dense_upoly_zp__mat_row_update_ir this state.matrix
            ⟨2, by omega⟩ ⟨3, by omega⟩ Q lenQ rotatedT state.lenT
            state.t scratch heap1 with
        | .error fault => .error fault
        | .ok row23 =>
          match dense_upoly_zp__mat_row_update_ir this row23.matrix
              ⟨0, by omega⟩ ⟨1, by omega⟩ Q lenQ row23.T row23.lenT
              row23.t scratch row23.heap with
          | .error fault => .error fault
          | .ok row01 =>
            hgcdIterLoop this m Q W3 scratch {
              heap := row01.heap
              matrix := row01.matrix
              A := rotatedA
              lenA := state.lenB
              B := rotatedB
              lenB := lenR
              T := row01.T
              lenT := row01.lenT
              t := row01.t
              sgn := -state.sgn
            }
    else
      .ok state
termination_by state => state.lenB
decreasing_by
  exact polyDivrem_remainder_lt this Q state.T state.A state.lenA state.B
    state.lenB W3 state.heap heap1 lenQ lenR hdiv

/-- Raw lowering of the initialization prefix of C++ `_hgcd_iter`.  The
source order of the two copies is kept because it is required for the
documented `{B,a}` aliasing case. -/
def hgcdIterInit (M : HgcdMat)
    (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) : RawExec HgcdIterState :=
  match dense_upoly_zp__mat_one_ir M heap with
  | .error fault => .error fault
  | .ok (heap1, matrix) =>
    match heap1.copyU64 A a lenA with
    | .error fault => .error fault
    | .ok heap2 =>
      match heap2.copyU64 B b lenB with
      | .error fault => .error fault
      | .ok heap3 =>
        .ok {
          heap := heap3
          matrix := matrix
          A := A
          lenA := lenA
          B := B
          lenB := lenB
          T := T
          lenT := lenT
          t := t
          sgn := 1
        }

/-- Raw lowering of the complete C++ `_hgcd_iter` initialization followed by
its well-founded Euclidean loop. -/
def dense_upoly_zp__hgcd_iter_ir (this : DenseUPolyZp) (M : HgcdMat)
    (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (heap : RawHeap) : RawExec HgcdIterState :=
  match hgcdIterInit M A B T t lenT a lenA b lenB heap with
  | .error fault => .error fault
  | .ok initial => hgcdIterLoop this (lenA / 2) Q W3 scratch initial

theorem matOne_result_valid (M : HgcdMat) (heap heap' : RawHeap)
    (matrix : HgcdMat)
    (hrun : dense_upoly_zp__mat_one_ir M heap = .ok (heap', matrix)) :
    matrix.Valid := by
  simp only [dense_upoly_zp__mat_one_ir] at hrun
  split at hrun
  next hvalid =>
    split at hrun
    next fault hwrite0 => simp at hrun
    next heap1 hwrite0 =>
      split at hrun
      next fault hwrite3 => simp at hrun
      next heap2 hwrite3 =>
        have heq := Except.ok.inj hrun
        cases heq
        exact ⟨hvalid.1, by simp⟩
  next hvalid => simp at hrun

theorem matOne_result_len (M : HgcdMat) (heap heap' : RawHeap)
    (matrix : HgcdMat)
    (hrun : dense_upoly_zp__mat_one_ir M heap = .ok (heap', matrix)) :
    matrix.len = #[1, 0, 0, 1] := by
  simp only [dense_upoly_zp__mat_one_ir] at hrun
  split at hrun
  next hvalid =>
    split at hrun
    next fault hwrite0 => simp at hrun
    next heap1 hwrite0 =>
      split at hrun
      next fault hwrite3 => simp at hrun
      next heap2 hwrite3 =>
        have heq := Except.ok.inj hrun
        cases heq
        rfl
  next hvalid => simp at hrun

/-- Descriptor and operand lengths installed by the real iterator prefix. -/
theorem hgcdIterInit_lengths (M : HgcdMat)
    (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (result : HgcdIterState)
    (hrun : hgcdIterInit M A B T t lenT a lenA b lenB heap = .ok result) :
    result.lenA = lenA ∧ result.lenB = lenB ∧
      result.matrix.len = #[1, 0, 0, 1] := by
  simp only [hgcdIterInit] at hrun
  split at hrun
  next fault hone => simp at hrun
  next heap1 matrix hone =>
    split at hrun
    next fault hcopyA => simp at hrun
    next heap2 hcopyA =>
      split at hrun
      next fault hcopyB => simp at hrun
      next heap3 hcopyB =>
        have heq := Except.ok.inj hrun
        subst result
        exact ⟨rfl, rfl, matOne_result_len M heap heap1 matrix hone⟩

theorem matRowUpdate_result_valid (this : DenseUPolyZp)
    (M : HgcdMat) (i0 i1 : Fin 4) (Q : RawPtr UInt64) (lenQ : Nat)
    (T : RawPtr UInt64) (lenT : Nat) (t scratch : RawPtr UInt64)
    (heap : RawHeap) (result : MatRowUpdateResult)
    (hrun : dense_upoly_zp__mat_row_update_ir this M i0 i1 Q lenQ T lenT
      t scratch heap = .ok result) :
    result.matrix.Valid := by
  simp only [dense_upoly_zp__mat_row_update_ir] at hrun
  split at hrun
  next hvalid =>
    split at hrun
    next hnonzero =>
      split at hrun
      next fault hmul => simp at hrun
      next heap1 hmul =>
        split at hrun
        next fault hadd => simp at hrun
        next heap2 sumLen hadd =>
          have heq := Except.ok.inj hrun
          cases heq
          exact ⟨by simp [hvalid.1], by simp [hvalid.2]⟩
    next hzero =>
      have heq := Except.ok.inj hrun
      cases heq
      exact ⟨by simp [hvalid.1], by simp [hvalid.2]⟩
  next hvalid => simp at hrun

theorem hgcdIterLoop_result_valid (this : DenseUPolyZp) (m : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (state result : HgcdIterState) (hState : state.matrix.Valid)
    (hrun : hgcdIterLoop this m Q W3 scratch state = .ok result) :
    result.matrix.Valid := by
  rw [hgcdIterLoop] at hrun
  split at hrun
  next hcontinue =>
    split at hrun
    next fault hdiv => simp at hrun
    next heap1 lenQ lenR hdiv =>
      dsimp only at hrun
      split at hrun
      next fault hrow23 => simp at hrun
      next row23 hrow23 =>
        split at hrun
        next fault hrow01 => simp at hrun
        next row01 hrow01 =>
          have hvalid01 := matRowUpdate_result_valid this row23.matrix
            ⟨0, by omega⟩ ⟨1, by omega⟩ Q lenQ row23.T row23.lenT
            row23.t scratch row23.heap row01 hrow01
          have hlt : lenR < state.lenB :=
            polyDivrem_remainder_lt this Q state.T state.A state.lenA state.B
              state.lenB W3 state.heap heap1 lenQ lenR hdiv
          exact hgcdIterLoop_result_valid this m Q W3 scratch _ result
            hvalid01 hrun
  next hstop =>
    have heq : state = result := Except.ok.inj hrun
    subst result
    exact hState
termination_by state.lenB
decreasing_by assumption

theorem hgcdIterInit_result_valid (M : HgcdMat)
    (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) (result : HgcdIterState)
    (hrun : hgcdIterInit M A B T t lenT a lenA b lenB heap = .ok result) :
    result.matrix.Valid := by
  simp only [hgcdIterInit] at hrun
  split at hrun
  next fault hone => simp at hrun
  next heap1 matrix hone =>
    split at hrun
    next fault hcopyA => simp at hrun
    next heap2 hcopyA =>
      split at hrun
      next fault hcopyB => simp at hrun
      next heap3 hcopyB =>
        have heq := Except.ok.inj hrun
        subst result
        exact matOne_result_valid M heap heap1 matrix hone

theorem hgcdIter_result_valid (this : DenseUPolyZp) (M : HgcdMat)
    (A B T t : RawPtr UInt64) (lenT : Nat)
    (a : RawPtr UInt64) (lenA : Nat) (b : RawPtr UInt64) (lenB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3) (scratch : RawPtr UInt64)
    (heap : RawHeap) (result : HgcdIterState)
    (hrun : dense_upoly_zp__hgcd_iter_ir this M A B T t lenT a lenA b lenB
      Q W3 scratch heap = .ok result) :
    result.matrix.Valid := by
  simp only [dense_upoly_zp__hgcd_iter_ir] at hrun
  split at hrun
  next fault hinit => simp at hrun
  next initial hinit =>
    exact hgcdIterLoop_result_valid this (lenA / 2) Q W3 scratch initial
      result (hgcdIterInit_result_valid M A B T t lenT a lenA b lenB heap
        initial hinit) hrun

/-- Reference-eliminated return state of the source `_hgcd_recursive` base
branch. -/
structure HgcdRecursiveBaseResult where
  heap : RawHeap
  matrix : HgcdMat
  lenA : Nat
  lenB : Nat
  sgn : Int

def HgcdIterState.toRecursiveBaseResult (state : HgcdIterState) :
    HgcdRecursiveBaseResult :=
  .mk state.heap state.matrix state.lenA state.lenB state.sgn

/-- Exact lowering of `_hgcd_recursive`'s `len_b < len_a / 2 + 1` branch.
The optional identity initialization and the alias-sensitive order of the two
source `memcpy` calls remain explicit. -/
def hgcdRecursiveBase (M : HgcdMat) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap) :
    RawExec HgcdRecursiveBaseResult :=
  let continueWith (heap1 : RawHeap) (matrix : HgcdMat) :=
    match heap1.copyU64 A a lenA with
    | .error fault => .error fault
    | .ok heap2 =>
      match heap2.copyU64 B b lenB with
      | .error fault => .error fault
      | .ok heap3 => .ok {
          heap := heap3
          matrix := matrix
          lenA := lenA
          lenB := lenB
          sgn := 1 }
  if computeM then
    match dense_upoly_zp__mat_one_ir M heap with
    | .error fault => .error fault
    | .ok (heap1, matrix) => continueWith heap1 matrix
  else
    continueWith heap M

theorem hgcdRecursiveBase_result_valid (M : HgcdMat) (hM : M.Valid)
    (computeM : Bool) (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (heap : RawHeap) (result : HgcdRecursiveBaseResult)
    (hrun : hgcdRecursiveBase M computeM A B a b lenA lenB heap =
      .ok result) :
    result.matrix.Valid := by
  simp only [hgcdRecursiveBase] at hrun
  split at hrun
  next hcompute =>
    split at hrun
    next fault hone => simp at hrun
    next heap1 matrix hone =>
      split at hrun
      next fault hcopyA => simp at hrun
      next heap2 hcopyA =>
        split at hrun
        next fault hcopyB => simp at hrun
        next heap3 hcopyB =>
          have heq := Except.ok.inj hrun
          subst result
          exact matOne_result_valid M heap heap1 matrix hone
  next hcompute =>
    split at hrun
    next fault hcopyA => simp at hrun
    next heap1 hcopyA =>
      split at hrun
      next fault hcopyB => simp at hrun
      next heap2 hcopyB =>
        have heq := Except.ok.inj hrun
        subst result
        exact hM

/-- Exact pointer slices created at the start of the non-base
`_hgcd_recursive` branch. -/
structure HgcdRecursiveWorkspace where
  n : Nat
  half : Nat
  a2 : RawPtr UInt64
  b2 : RawPtr UInt64
  a3 : RawPtr UInt64
  b3 : RawPtr UInt64
  q : RawPtr UInt64
  d : RawPtr UInt64
  T0 : RawPtr UInt64
  T1 : RawPtr UInt64
  R : HgcdMat
  S : HgcdMat
  W3 : RawPtr Word3
  next : RawPtr UInt64

/-- Faithful lowering of the source's pointer-arithmetic workspace block.
The zero length arrays are the total Lean value for C++ length slots that are
written by the recursive/iterator call before any source read. -/
def hgcdRecursiveWorkspace (W : RawPtr UInt64) (lenA : Nat) :
    HgcdRecursiveWorkspace :=
  let n := lenA
  let half := (n + 1) / 2
  let a2 := W
  let b2 := a2.add n
  let a3 := b2.add n
  let b3 := a3.add n
  let q := b3.add n
  let d := q.add half
  let T0 := d.add n
  let T1 := T0.add n
  let r0 := T1.add half
  let r1 := r0.add half
  let r2 := r1.add half
  let r3 := r2.add half
  let R : HgcdMat := { poly := #[r0, r1, r2, r3], len := #[0, 0, 0, 0] }
  let s0 := r3.add half
  let s1 := s0.add half
  let s2 := s1.add half
  let s3 := s2.add half
  let S : HgcdMat := { poly := #[s0, s1, s2, s3], len := #[0, 0, 0, 0] }
  let W3 : RawPtr Word3 := RawPtr.reinterpret (s3.add half)
  let next := W.add (6 * n + 10 * half + 3 * n)
  { n, half, a2, b2, a3, b3, q, d, T0, T1, R, S, W3, next }

/-- Auditable equations for every source pointer slice and both four-entry
matrix descriptors. -/
theorem hgcdRecursiveWorkspace_layout (W : RawPtr UInt64) (lenA : Nat) :
    let ws := hgcdRecursiveWorkspace W lenA
    let r3 := (((ws.T1.add ws.half).add ws.half).add ws.half).add ws.half
    ws.n = lenA ∧ ws.half = (lenA + 1) / 2 ∧
    ws.a2 = W ∧ ws.b2 = ws.a2.add ws.n ∧
    ws.a3 = ws.b2.add ws.n ∧ ws.b3 = ws.a3.add ws.n ∧
    ws.q = ws.b3.add ws.n ∧ ws.d = ws.q.add ws.half ∧
    ws.T0 = ws.d.add ws.n ∧ ws.T1 = ws.T0.add ws.n ∧
    ws.R.poly = #[ws.T1.add ws.half,
      (ws.T1.add ws.half).add ws.half,
      ((ws.T1.add ws.half).add ws.half).add ws.half,
      (((ws.T1.add ws.half).add ws.half).add ws.half).add ws.half] ∧
    ws.S.poly = #[r3.add ws.half,
        (r3.add ws.half).add ws.half,
        ((r3.add ws.half).add ws.half).add ws.half,
        (((r3.add ws.half).add ws.half).add ws.half).add ws.half]
      ∧ ws.R.Valid ∧ ws.S.Valid ∧
    ws.next = W.add (6 * lenA + 10 * ((lenA + 1) / 2) + 3 * lenA) := by
  simp [hgcdRecursiveWorkspace, HgcdMat.Valid]

/-- High-half input passed to the first recursive/iterator call in the
non-base source branch. -/
structure HgcdRecursiveHighInput where
  m : Nat
  a0 : RawPtr UInt64
  lenA0 : Nat
  b0 : RawPtr UInt64
  lenB0 : Nat

def hgcdRecursiveHighInput (a b : RawPtr UInt64) (lenA lenB : Nat) :
    HgcdRecursiveHighInput :=
  let m := lenA / 2
  { m := m
    a0 := a.add m
    lenA0 := lenA - m
    b0 := b.add m
    lenB0 := if lenB ≥ m then lenB - m else 0 }

theorem hgcdRecursiveHighInput_layout (a b : RawPtr UInt64)
    (lenA lenB : Nat) :
    let input := hgcdRecursiveHighInput a b lenA lenB
    input.m = lenA / 2 ∧ input.a0 = a.add input.m ∧
      input.lenA0 = lenA - input.m ∧ input.b0 = b.add input.m ∧
      input.lenB0 = if lenB ≥ input.m then lenB - input.m else 0 := by
  simp [hgcdRecursiveHighInput]

/-- The first true recursive call decreases the source `len_a` measure. -/
theorem hgcdRecursiveHighInput_len_lt (a b : RawPtr UInt64)
    (lenA lenB : Nat) (horder : lenB < lenA)
    (hnonbase : lenA / 2 + 1 ≤ lenB) :
    (hgcdRecursiveHighInput a b lenA lenB).lenA0 < lenA := by
  simp only [hgcdRecursiveHighInput]
  omega

/-- Removing the same source high-half offset preserves the strict operand
length order required by the first recursive HGCD call. -/
theorem hgcdRecursiveHighInput_order (a b : RawPtr UInt64)
    (lenA lenB : Nat) (horder : lenB < lenA) :
    (hgcdRecursiveHighInput a b lenA lenB).lenB0 <
      (hgcdRecursiveHighInput a b lenA lenB).lenA0 := by
  simp only [hgcdRecursiveHighInput]
  split <;> omega

theorem hgcdRecursiveHighInput_lenA_pos (a b : RawPtr UInt64)
    (lenA lenB : Nat) (horder : lenB < lenA) :
    0 < (hgcdRecursiveHighInput a b lenA lenB).lenA0 := by
  simp only [hgcdRecursiveHighInput]
  omega

/-- Heap and accumulated offset after the first matrix-stabilization loop. -/
structure HgcdMatStageResult where
  heap : RawHeap
  off : Nat

def hgcdMatPtrRaw (M : HgcdMat) (hM : M.Valid) (i : Fin 4) :
    RawPtr UInt64 := M.poly[i.val]'(by rw [hM.1]; exact i.isLt)

def hgcdMatLenRaw (M : HgcdMat) (hM : M.Valid) (i : Fin 4) : Nat :=
  M.len[i.val]'(by rw [hM.2]; exact i.isLt)

/-- Source-order offset of an entry in the contiguous stabilization buffer.
Indices at or beyond four denote the end of the four-entry block. -/
def hgcdMatStageOffset (M : HgcdMat) (hM : M.Valid) : Nat → Nat
  | 0 => 0
  | 1 => hgcdMatLenRaw M hM ⟨0, by omega⟩
  | 2 => hgcdMatLenRaw M hM ⟨0, by omega⟩ +
      hgcdMatLenRaw M hM ⟨1, by omega⟩
  | 3 => hgcdMatLenRaw M hM ⟨0, by omega⟩ +
      hgcdMatLenRaw M hM ⟨1, by omega⟩ +
      hgcdMatLenRaw M hM ⟨2, by omega⟩
  | _ => hgcdMatLenRaw M hM ⟨0, by omega⟩ +
      hgcdMatLenRaw M hM ⟨1, by omega⟩ +
      hgcdMatLenRaw M hM ⟨2, by omega⟩ +
      hgcdMatLenRaw M hM ⟨3, by omega⟩

def hgcdMatStageSize (M : HgcdMat) (hM : M.Valid) : Nat :=
  hgcdMatStageOffset M hM 4

theorem hgcdMatStageOffset_step (M : HgcdMat) (hM : M.Valid)
    (i : Nat) (hi : i < 4) :
    hgcdMatStageOffset M hM (i + 1) = hgcdMatStageOffset M hM i +
      hgcdMatLenRaw M hM ⟨i, hi⟩ := by
  have hcases : i = 0 ∨ i = 1 ∨ i = 2 ∨ i = 3 := by omega
  rcases hcases with rfl | rfl | rfl | rfl <;>
    simp [hgcdMatStageOffset] <;> omega

theorem hgcdMatStageOffset_entry_le_size (M : HgcdMat) (hM : M.Valid)
    (i : Nat) (hi : i < 4) :
    hgcdMatStageOffset M hM i + hgcdMatLenRaw M hM ⟨i, hi⟩ ≤
      hgcdMatStageSize M hM := by
  have hcases : i = 0 ∨ i = 1 ∨ i = 2 ∨ i = 3 := by omega
  rcases hcases with rfl | rfl | rfl | rfl <;>
    simp [hgcdMatStageSize, hgcdMatStageOffset] <;> omega

theorem hgcdMatStageOffset_entry_le_later (M : HgcdMat) (hM : M.Valid)
    (j i : Nat) (hj : j < i) (hi : i < 4) (hj4 : j < 4) :
    hgcdMatStageOffset M hM j + hgcdMatLenRaw M hM ⟨j, hj4⟩ ≤
      hgcdMatStageOffset M hM i := by
  have hicases : i = 1 ∨ i = 2 ∨ i = 3 := by omega
  have hjcases : j = 0 ∨ j = 1 ∨ j = 2 := by omega
  rcases hicases with rfl | rfl | rfl <;>
    rcases hjcases with rfl | rfl | rfl <;>
    simp [hgcdMatStageOffset] <;> omega

/-- Exact first `for (i=0; i<4; ++i)` stabilization loop: copy every live
matrix entry to consecutive cells of the source `stage` buffer. -/
def hgcdMatStageLoop (M : HgcdMat) (hM : M.Valid)
    (stage : RawPtr UInt64) :
    (i off : Nat) → (heap : RawHeap) → RawExec HgcdMatStageResult
  | i, off, heap =>
    if hi : i < 4 then
      let index : Fin 4 := ⟨i, hi⟩
      match heap.copyU64 (stage.add off) (hgcdMatPtrRaw M hM index)
          (hgcdMatLenRaw M hM index) with
      | .error fault => .error fault
      | .ok heap1 => hgcdMatStageLoop M hM stage (i + 1)
          (off + hgcdMatLenRaw M hM index) heap1
    else
      .ok { heap := heap, off := off }
termination_by i off heap => 4 - i
decreasing_by omega

/-- The actual first generated loop advances its offset by exactly the sum of
the four live descriptor lengths. -/
theorem hgcdMatStageLoop_final_off (M : HgcdMat) (hM : M.Valid)
    (stage : RawPtr UInt64) (i off : Nat) (heap : RawHeap)
    (result : HgcdMatStageResult) (hi : i ≤ 4)
    (hrun : hgcdMatStageLoop M hM stage i off heap = .ok result) :
    result.off + hgcdMatStageOffset M hM i =
      off + hgcdMatStageSize M hM := by
  rw [hgcdMatStageLoop] at hrun
  split at hrun
  next hlt =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      have hrec := hgcdMatStageLoop_final_off M hM stage (i + 1)
        (off + hgcdMatLenRaw M hM ⟨i, hlt⟩) heap1 result
        (by omega) hrun
      rw [hgcdMatStageOffset_step M hM i hlt] at hrec
      omega
  next hstop =>
    have hi4 : i = 4 := by omega
    subst i
    have heq : ({ heap := heap, off := off } : HgcdMatStageResult) = result :=
      Except.ok.inj hrun
    subst result
    simp [hgcdMatStageSize]
termination_by 4 - i
decreasing_by omega

/-- Heap, restored descriptor, and accumulated offset after the second
matrix-stabilization loop. -/
structure HgcdMatRestoreResult where
  heap : RawHeap
  matrix : HgcdMat
  off : Nat

/-- Exact second stabilization loop: copy staged entries to the four saved
pre-iterator pointers and restore each descriptor pointer after its copy. -/
def hgcdMatRestoreLoop (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64) :
    (i off : Nat) → (heap : RawHeap) → RawExec HgcdMatRestoreResult
  | i, off, heap =>
    if hi : i < 4 then
      let index : Fin 4 := ⟨i, hi⟩
      match heap.copyU64 (hgcdMatPtrRaw original hOriginal index)
          (stage.add off) (hgcdMatLenRaw current hCurrent index) with
      | .error fault => .error fault
      | .ok heap1 =>
        let poly' := current.poly.set i (hgcdMatPtrRaw original hOriginal index)
          (by rw [hCurrent.1]; omega)
        let next : HgcdMat := { current with poly := poly' }
        have hNext : next.Valid := by
          exact ⟨by simp [next, poly', hCurrent.1], hCurrent.2⟩
        hgcdMatRestoreLoop original next hOriginal hNext stage (i + 1)
          (off + hgcdMatLenRaw current hCurrent index) heap1
    else
      .ok { heap := heap, matrix := current, off := off }
termination_by i off heap => 4 - i
decreasing_by omega

/-- The complete two-loop stabilization block used after either iterator
call in `_hgcd_recursive`. -/
def hgcdMatStabilize (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64) (heap : RawHeap) : RawExec HgcdMatRestoreResult :=
  match hgcdMatStageLoop current hCurrent stage 0 0 heap with
  | .error fault => .error fault
  | .ok staged => hgcdMatRestoreLoop original current hOriginal hCurrent
      stage 0 0 staged.heap

/-- Exact post-iterator pointer normalization in `_hgcd_recursive`.  The
cross-protection branch must save `pB = a3` into `b3` before overwriting `a3`
from `pA`; all other branches retain the source's two conditional copies. -/
def hgcdRecursiveStoreIterOutputs (a3 b3 pA pB : RawPtr UInt64)
    (lenA3 lenB3 : Nat) (heap : RawHeap) : RawExec RawHeap :=
  if !(pA == a3) && (pB == a3) then
    match heap.copyU64 b3 pB lenB3 with
    | .error fault => .error fault
    | .ok heap1 => heap1.copyU64 a3 pA lenA3
  else
    match if !(pA == a3) then heap.copyU64 a3 pA lenA3 else .ok heap with
    | .error fault => .error fault
    | .ok heap1 =>
      if !(pB == b3) then heap1.copyU64 b3 pB lenB3 else .ok heap1

/-- Return state of either iterator arm inside `_hgcd_recursive`. -/
structure HgcdRecursiveIterBranchResult where
  heap : RawHeap
  matrix : HgcdMat
  lenA : Nat
  lenB : Nat
  sgn : Int

/-- Exact source composition for an iterator arm of `_hgcd_recursive`:
run `_hgcd_iter`, stabilize all four matrix entries at their saved pointers,
then perform the alias-sensitive `pA/pB` copies. -/
def hgcdRecursiveIterBranch (this : DenseUPolyZp)
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage : RawPtr UInt64) (heap : RawHeap) :
    RawExec HgcdRecursiveIterBranchResult :=
  match hiter : dense_upoly_zp__hgcd_iter_ir this original a3 b3 T0 T1 0
      inputA lenInputA inputB lenInputB Q W3 scratch heap with
  | .error fault => .error fault
  | .ok iterResult =>
    have hIterValid : iterResult.matrix.Valid :=
      hgcdIter_result_valid this original a3 b3 T0 T1 0 inputA lenInputA
        inputB lenInputB Q W3 scratch heap iterResult hiter
    match hgcdMatStabilize original iterResult.matrix hOriginal hIterValid
        stage iterResult.heap with
    | .error fault => .error fault
    | .ok stable =>
      match hgcdRecursiveStoreIterOutputs a3 b3 iterResult.A iterResult.B
          iterResult.lenA iterResult.lenB stable.heap with
      | .error fault => .error fault
      | .ok heap1 => .ok {
          heap := heap1
          matrix := stable.matrix
          lenA := iterResult.lenA
          lenB := iterResult.lenB
          sgn := iterResult.sgn }

/-- A successful iterator arm necessarily contains those same three source
executions, in order; this rules out later replacing stabilization or output
copies with a specification-only result. -/
theorem hgcdRecursiveIterBranch_exec (this : DenseUPolyZp)
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage : RawPtr UInt64) (heap : RawHeap)
    (result : HgcdRecursiveIterBranchResult)
    (hrun : hgcdRecursiveIterBranch this original hOriginal a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap = .ok result) :
    ∃ iterResult, ∃ hIterValid : iterResult.matrix.Valid, ∃ stable,
      dense_upoly_zp__hgcd_iter_ir this original a3 b3 T0 T1 0 inputA
          lenInputA inputB lenInputB Q W3 scratch heap = .ok iterResult ∧
      hgcdMatStabilize original iterResult.matrix hOriginal
          hIterValid
          stage iterResult.heap = .ok stable ∧
      hgcdRecursiveStoreIterOutputs a3 b3 iterResult.A iterResult.B
          iterResult.lenA iterResult.lenB stable.heap = .ok result.heap ∧
      result.matrix = stable.matrix ∧ result.lenA = iterResult.lenA ∧
      result.lenB = iterResult.lenB ∧ result.sgn = iterResult.sgn := by
  simp only [hgcdRecursiveIterBranch] at hrun
  split at hrun
  next fault hiter => simp at hrun
  next iterResult hiter =>
    split at hrun
    next fault hstable => simp at hrun
    next stable hstable =>
      split at hrun
      next fault hstore => simp at hrun
      next heap1 hstore =>
        have heq := Except.ok.inj hrun
        cases heq
        let hIterValid := hgcdIter_result_valid this original a3 b3 T0 T1 0
          inputA lenInputA inputB lenInputB Q W3 scratch heap iterResult hiter
        exact ⟨iterResult, hIterValid, stable, hiter, hstable, hstore,
          rfl, rfl, rfl, rfl⟩

theorem hgcdRecursiveStoreIterOutputs_cross_exec
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (heap heap1 heap2 : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = true)
    (hcopyB : heap.copyU64 b3 pB lenB3 = .ok heap1)
    (hcopyA : heap1.copyU64 a3 pA lenA3 = .ok heap2) :
    hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
      .ok heap2 := by
  simp [hgcdRecursiveStoreIterOutputs, hcross, hcopyB, hcopyA]

theorem hgcdRecursiveStoreIterOutputs_regular_exec
    (a3 b3 pA pB : RawPtr UInt64) (lenA3 lenB3 : Nat)
    (heap heap1 heap2 : RawHeap)
    (hcross : (!(pA == a3) && (pB == a3)) = false)
    (hfirst : (if (pA == a3) = false then heap.copyU64 a3 pA lenA3
      else .ok heap) = .ok heap1)
    (hsecond : (if (pB == b3) = false then heap1.copyU64 b3 pB lenB3
      else .ok heap1) = .ok heap2) :
    hgcdRecursiveStoreIterOutputs a3 b3 pA pB lenA3 lenB3 heap =
      .ok heap2 := by
  simp [hgcdRecursiveStoreIterOutputs, hcross, hfirst, hsecond]

structure HgcdMulTermResult where
  heap : RawHeap
  length : Nat

/-- Exact guarded multiplication idiom used four times while reconstructing
the low halves after the first HGCD call. -/
def hgcdRecursiveMulTerm (this : DenseUPolyZp) (dst left : RawPtr UInt64)
    (lenLeft : Nat) (right : RawPtr UInt64) (lenRight : Nat)
    (scratch : RawPtr UInt64) (heap : RawHeap) : RawExec HgcdMulTermResult :=
  if lenLeft > 0 && lenRight > 0 then
    let run := if lenLeft ≥ lenRight then
      dense_upoly_zp__mul_ir this dst left lenLeft right lenRight scratch heap
    else
      dense_upoly_zp__mul_ir this dst right lenRight left lenLeft scratch heap
    match run with
    | .error fault => .error fault
    | .ok heap1 => .ok { heap := heap1, length := lenLeft + lenRight - 1 }
  else
    .ok { heap := heap, length := 0 }

/-- The physical output length selected by the guarded source block fits the
full product buffer for the longer operand. -/
theorem hgcdRecursiveMulTerm_length_le (this : DenseUPolyZp)
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat)
    (scratch : RawPtr UInt64) (heap : RawHeap) (result : HgcdMulTermResult)
    (hrun : hgcdRecursiveMulTerm this dst left lenLeft right lenRight
      scratch heap = .ok result) :
    result.length ≤ 2 * max lenLeft lenRight - 1 := by
  simp only [hgcdRecursiveMulTerm] at hrun
  split at hrun
  next hnonzero =>
    split at hrun
    next fault hmul => simp at hrun
    next heap1 hmul =>
      have heq : result =
          HgcdMulTermResult.mk heap1 (lenLeft + lenRight - 1) :=
        (Except.ok.inj hrun).symm
      subst result
      simp at hnonzero
      rcases Nat.le_total lenLeft lenRight with hle | hle
      · simp [Nat.max_eq_right hle]
        omega
      · simp [Nat.max_eq_left hle]
        omega
  next hzero =>
    have heq : result = HgcdMulTermResult.mk heap 0 :=
      (Except.ok.inj hrun).symm
    subst result
    exact Nat.zero_le _

/-- Sharper source-length bound used by the enclosing recursive measure. -/
theorem hgcdRecursiveMulTerm_length_le_sum (this : DenseUPolyZp)
    (dst left : RawPtr UInt64) (lenLeft : Nat)
    (right : RawPtr UInt64) (lenRight : Nat)
    (scratch : RawPtr UInt64) (heap : RawHeap) (result : HgcdMulTermResult)
    (hrun : hgcdRecursiveMulTerm this dst left lenLeft right lenRight
      scratch heap = .ok result) :
    result.length ≤ lenLeft + lenRight := by
  simp only [hgcdRecursiveMulTerm] at hrun
  split at hrun
  next hnonzero =>
    split at hrun
    next fault hmul => simp at hrun
    next heap1 hmul =>
      have heq : result =
          HgcdMulTermResult.mk heap1 (lenLeft + lenRight - 1) :=
        (Except.ok.inj hrun).symm
      subst result
      change lenLeft + lenRight - 1 ≤ lenLeft + lenRight
      omega
  next hzero =>
    have heq : result = HgcdMulTermResult.mk heap 0 :=
      (Except.ok.inj hrun).symm
    subst result
    change 0 ≤ lenLeft + lenRight
    omega

/-- Exact `b2` low-half reconstruction block:
`R[2]*a_lo` and `R[0]*b_lo`, followed by the sign-selected subtraction. -/
def hgcdRecursiveReconstructB (this : DenseUPolyZp)
    (b2 T0 r2 r0 aLo bLo scratch : RawPtr UInt64)
    (lenR2 lenR0 lenALo lenBLo : Nat) (sgn : Int)
    (heap : RawHeap) : RawExec (RawHeap × Nat) :=
  match hgcdRecursiveMulTerm this b2 r2 lenR2 aLo lenALo scratch heap with
  | .error fault => .error fault
  | .ok term1 =>
    match hgcdRecursiveMulTerm this T0 r0 lenR0 bLo lenBLo scratch term1.heap with
    | .error fault => .error fault
    | .ok term2 =>
      if sgn < 0 then
        dense_upoly_zp__poly_sub_ir this b2 b2 term1.length T0 term2.length
          term2.heap
      else
        dense_upoly_zp__poly_sub_ir this b2 T0 term2.length b2 term1.length
          term2.heap

/-- Exact `a2` low-half reconstruction block:
`R[3]*a_lo` and `R[1]*b_lo`, with the source's opposite sign orientation. -/
def hgcdRecursiveReconstructA (this : DenseUPolyZp)
    (a2 T0 r3 r1 aLo bLo scratch : RawPtr UInt64)
    (lenR3 lenR1 lenALo lenBLo : Nat) (sgn : Int)
    (heap : RawHeap) : RawExec (RawHeap × Nat) :=
  match hgcdRecursiveMulTerm this a2 r3 lenR3 aLo lenALo scratch heap with
  | .error fault => .error fault
  | .ok term1 =>
    match hgcdRecursiveMulTerm this T0 r1 lenR1 bLo lenBLo scratch term1.heap with
    | .error fault => .error fault
    | .ok term2 =>
      if sgn < 0 then
        dense_upoly_zp__poly_sub_ir this a2 T0 term2.length a2 term1.length
          term2.heap
      else
        dense_upoly_zp__poly_sub_ir this a2 a2 term1.length T0 term2.length
          term2.heap

structure HgcdLiftHighResult where
  heap : RawHeap
  length : Nat

/-- Exact lowering of the two identical “zero-fill, add the shifted high
half, normalize the whole buffer” blocks used for `b2` and `a2`. -/
def hgcdRecursiveLiftHigh (this : DenseUPolyZp)
    (out high : RawPtr UInt64) (lowLength m highLength : Nat)
    (heap : RawHeap) : RawExec HgcdLiftHighResult :=
  let required := m + highLength
  let fullLength := max required lowLength
  let zeroRun := if lowLength < required then
    mulZeroPadLoop out lowLength (required - lowLength) 0 heap
  else .ok heap
  match zeroRun with
  | .error fault => .error fault
  | .ok heap1 =>
    let outHigh := out.add m
    let oldHighLength := if m ≤ lowLength then lowLength - m else 0
    match dense_upoly_zp__poly_add_ir this outHigh outHigh oldHighLength
        high highLength heap1 with
    | .error fault => .error fault
    | .ok (heap2, _) =>
      match heap2.normaliseU64 out fullLength with
      | .error fault => .error fault
      | .ok length => .ok { heap := heap2, length := length }

/-- Return state of either complete A/B reconstruction pair in
`_hgcd_recursive`. -/
structure HgcdRecursiveReconstructPairResult where
  heap : RawHeap
  lenA : Nat
  lenB : Nat

/-- Exact source-order composition shared by both reconstruction sites:
reconstruct and lift B first, then reconstruct and lift A.  Instantiating
`shift=m` with the first matrix gives the first site; instantiating `shift=k`
with the second matrix gives the final output site. -/
def hgcdRecursiveReconstructPair (this : DenseUPolyZp)
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift : Nat)
    (M : HgcdMat) (hM : M.Valid) (sgn : Int) (heap : RawHeap) :
    RawExec HgcdRecursiveReconstructPairResult :=
  match hgcdRecursiveReconstructB this B T0
      (hgcdMatPtrRaw M hM (2 : Fin 4))
      (hgcdMatPtrRaw M hM (0 : Fin 4)) lowA lowB scratch
      (hgcdMatLenRaw M hM (2 : Fin 4))
      (hgcdMatLenRaw M hM (0 : Fin 4)) lenLowA lenLowB sgn heap with
  | .error fault => .error fault
  | .ok (heap1, lowLenB) =>
    match hgcdRecursiveLiftHigh this B highB lowLenB shift lenHighB heap1 with
    | .error fault => .error fault
    | .ok liftedB =>
      match hgcdRecursiveReconstructA this A T0
          (hgcdMatPtrRaw M hM (3 : Fin 4))
          (hgcdMatPtrRaw M hM (1 : Fin 4)) lowA lowB scratch
          (hgcdMatLenRaw M hM (3 : Fin 4))
          (hgcdMatLenRaw M hM (1 : Fin 4)) lenLowA lenLowB sgn
          liftedB.heap with
      | .error fault => .error fault
      | .ok (heap3, lowLenA) =>
        match hgcdRecursiveLiftHigh this A highA lowLenA shift lenHighA
            heap3 with
        | .error fault => .error fault
        | .ok liftedA => .ok {
            heap := liftedA.heap
            lenA := liftedA.length
            lenB := liftedB.length }

/-- Successful pair reconstruction exposes the exact four source blocks in
their C++ order and pins both returned lengths to the two normalization
results. -/
theorem hgcdRecursiveReconstructPair_exec (this : DenseUPolyZp)
    (A B T0 lowA lowB highA highB scratch : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift : Nat)
    (M : HgcdMat) (hM : M.Valid) (sgn : Int) (heap : RawHeap)
    (result : HgcdRecursiveReconstructPairResult)
    (hrun : hgcdRecursiveReconstructPair this A B T0 lowA lowB highA highB
      scratch lenLowA lenLowB lenHighA lenHighB shift M hM sgn heap =
        .ok result) :
    ∃ heap1 lowLenB liftedB heap3 lowLenA liftedA,
      hgcdRecursiveReconstructB this B T0
          (hgcdMatPtrRaw M hM (2 : Fin 4))
          (hgcdMatPtrRaw M hM (0 : Fin 4)) lowA lowB scratch
          (hgcdMatLenRaw M hM (2 : Fin 4))
          (hgcdMatLenRaw M hM (0 : Fin 4)) lenLowA lenLowB sgn heap =
        .ok (heap1, lowLenB) ∧
      hgcdRecursiveLiftHigh this B highB lowLenB shift lenHighB heap1 =
        .ok liftedB ∧
      hgcdRecursiveReconstructA this A T0
          (hgcdMatPtrRaw M hM (3 : Fin 4))
          (hgcdMatPtrRaw M hM (1 : Fin 4)) lowA lowB scratch
          (hgcdMatLenRaw M hM (3 : Fin 4))
          (hgcdMatLenRaw M hM (1 : Fin 4)) lenLowA lenLowB sgn
          liftedB.heap = .ok (heap3, lowLenA) ∧
      hgcdRecursiveLiftHigh this A highA lowLenA shift lenHighA heap3 =
        .ok liftedA ∧
      result.heap = liftedA.heap ∧ result.lenA = liftedA.length ∧
      result.lenB = liftedB.length := by
  simp only [hgcdRecursiveReconstructPair] at hrun
  generalize hB : hgcdRecursiveReconstructB this B T0
      (hgcdMatPtrRaw M hM (2 : Fin 4))
      (hgcdMatPtrRaw M hM (0 : Fin 4)) lowA lowB scratch
      (hgcdMatLenRaw M hM (2 : Fin 4))
      (hgcdMatLenRaw M hM (0 : Fin 4)) lenLowA lenLowB sgn heap = runB
    at hrun
  cases runB with
  | error fault => simp [hB] at hrun
  | ok pairB =>
    rcases pairB with ⟨heap1, lowLenB⟩
    simp only [hB] at hrun
    generalize hLiftB : hgcdRecursiveLiftHigh this B highB lowLenB shift
      lenHighB heap1 = runLiftB at hrun
    cases runLiftB with
    | error fault => simp [hLiftB] at hrun
    | ok liftedB =>
      simp only [hLiftB] at hrun
      generalize hA : hgcdRecursiveReconstructA this A T0
          (hgcdMatPtrRaw M hM (3 : Fin 4))
          (hgcdMatPtrRaw M hM (1 : Fin 4)) lowA lowB scratch
          (hgcdMatLenRaw M hM (3 : Fin 4))
          (hgcdMatLenRaw M hM (1 : Fin 4)) lenLowA lenLowB sgn
          liftedB.heap = runA at hrun
      cases runA with
      | error fault => simp [hA] at hrun
      | ok pairA =>
        rcases pairA with ⟨heap3, lowLenA⟩
        simp only [hA] at hrun
        generalize hLiftA : hgcdRecursiveLiftHigh this A highA lowLenA shift
          lenHighA heap3 = runLiftA at hrun
        cases runLiftA with
        | error fault => simp [hLiftA] at hrun
        | ok liftedA =>
          have hactual :
              (Except.ok {
                heap := liftedA.heap
                lenA := liftedA.length
                lenB := liftedB.length } :
                RawExec HgcdRecursiveReconstructPairResult) = .ok result := by
            simpa [hLiftA] using hrun
          have heq := Except.ok.inj hactual
          cases heq
          exact ⟨heap1, lowLenB, liftedB, heap3, lowLenA, liftedA,
            rfl, hLiftB, hA, hLiftA, rfl, rfl, rfl⟩

/-- Raw result of the C++ `_mat_mul_entry` helper. -/
structure HgcdMatMulEntryResult where
  heap : RawHeap
  length : Nat

/-- Exact lowering of `_mat_mul_entry`: compute the two guarded products in
source order, then select add, first-only, copied-second, or zero. -/
def hgcdMatMulEntry (this : DenseUPolyZp)
    (C P Q R S T scratch : RawPtr UInt64)
    (lenP lenQ lenR lenS : Nat) (heap : RawHeap) :
    RawExec HgcdMatMulEntryResult :=
  match hgcdRecursiveMulTerm this C P lenP Q lenQ scratch heap with
  | .error fault => .error fault
  | .ok productPQ =>
    match hgcdRecursiveMulTerm this T R lenR S lenS scratch productPQ.heap with
    | .error fault => .error fault
    | .ok productRS =>
      if productPQ.length > 0 && productRS.length > 0 then
        match dense_upoly_zp__poly_add_ir this C C productPQ.length T
            productRS.length productRS.heap with
        | .error fault => .error fault
        | .ok (heap3, length) => .ok { heap := heap3, length := length }
      else if productPQ.length > 0 then
        .ok { heap := productRS.heap, length := productPQ.length }
      else if productRS.length > 0 then
        match productRS.heap.copyU64 C T productRS.length with
        | .error fault => .error fault
        | .ok heap3 => .ok { heap := heap3, length := productRS.length }
      else
        .ok { heap := productRS.heap, length := 0 }

/-- Successful `_mat_mul_entry` exposes both actual product executions and
the exact source post-product branch. -/
theorem hgcdMatMulEntry_exec (this : DenseUPolyZp)
    (C P Q R S T scratch : RawPtr UInt64)
    (lenP lenQ lenR lenS : Nat) (heap : RawHeap)
    (result : HgcdMatMulEntryResult)
    (hrun : hgcdMatMulEntry this C P Q R S T scratch lenP lenQ lenR lenS
      heap = .ok result) :
    ∃ productPQ productRS,
      hgcdRecursiveMulTerm this C P lenP Q lenQ scratch heap = .ok productPQ ∧
      hgcdRecursiveMulTerm this T R lenR S lenS scratch productPQ.heap =
        .ok productRS ∧
      (if productPQ.length > 0 && productRS.length > 0 then
        ∃ length, dense_upoly_zp__poly_add_ir this C C productPQ.length T
            productRS.length productRS.heap = .ok (result.heap, length) ∧
          result.length = length
       else if productPQ.length > 0 then
        result.heap = productRS.heap ∧ result.length = productPQ.length
       else if productRS.length > 0 then
        productRS.heap.copyU64 C T productRS.length = .ok result.heap ∧
          result.length = productRS.length
       else
        result.heap = productRS.heap ∧ result.length = 0) := by
  simp only [hgcdMatMulEntry] at hrun
  generalize hPQ : hgcdRecursiveMulTerm this C P lenP Q lenQ scratch heap =
    runPQ at hrun ⊢
  cases runPQ with
  | error fault => simp at hrun
  | ok productPQ =>
    generalize hRS : hgcdRecursiveMulTerm this T R lenR S lenS scratch
      productPQ.heap = runRS at hrun ⊢
    cases runRS with
    | error fault => simp [hRS] at hrun
    | ok productRS =>
      simp only [hRS] at hrun
      refine ⟨productPQ, productRS, rfl, hRS, ?_⟩
      split at hrun
      next hboth =>
        simp only [hboth, ↓reduceIte]
        generalize hadd : dense_upoly_zp__poly_add_ir this C C
          productPQ.length T productRS.length productRS.heap = addRun at hrun ⊢
        cases addRun with
        | error fault => simp at hrun
        | ok pair =>
          rcases pair with ⟨heap3, length⟩
          have heq := Except.ok.inj hrun
          cases heq
          exact ⟨length, rfl, rfl⟩
      next hboth =>
        simp only [hboth, ↓reduceIte]
        split at hrun
        next hPQPos =>
          simp only [hPQPos, ↓reduceIte]
          have heq := Except.ok.inj hrun
          cases heq
          exact ⟨rfl, rfl⟩
        next hPQPos =>
          simp only [hPQPos, ↓reduceIte]
          split at hrun
          next hRSPos =>
            simp only [hRSPos, ↓reduceIte]
            generalize hcopy : productRS.heap.copyU64 C T productRS.length =
              copyRun at hrun ⊢
            cases copyRun with
            | error fault => simp at hrun
            | ok heap3 =>
              have heq := Except.ok.inj hrun
              cases heq
              exact ⟨rfl, rfl⟩
          next hRSPos =>
            simp only [hRSPos, ↓reduceIte]
            have heq := Except.ok.inj hrun
            cases heq
            exact ⟨rfl, rfl⟩

/-- Raw result of the four-entry C++ `_mat_mul`. -/
structure HgcdMatMulResult where
  heap : RawHeap
  matrix : HgcdMat

/-- Exact source-order lowering of `_mat_mul`.  Entry `i` uses row
`i/2` and column `i%2`, hence the two products
`A[2r]*B[c] + A[2r+1]*B[2+c]`. -/
def hgcdMatMulLoop (this : DenseUPolyZp)
    (A B : HgcdMat) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) :
    (C : HgcdMat) → (hC : C.Valid) → (i : Nat) →
      (heap : RawHeap) → RawExec HgcdMatMulResult
  | C, hC, i, heap =>
    if hi : i < 4 then
      let out : Fin 4 := ⟨i, hi⟩
      let rowBase : Fin 4 := ⟨2 * (i / 2), by omega⟩
      let rowNext : Fin 4 := ⟨2 * (i / 2) + 1, by omega⟩
      let col : Fin 4 := ⟨i % 2, by omega⟩
      let lowerCol : Fin 4 := ⟨2 + i % 2, by omega⟩
      match hgcdMatMulEntry this
          (hgcdMatPtrRaw C hC out)
          (hgcdMatPtrRaw A hA rowBase) (hgcdMatPtrRaw B hB col)
          (hgcdMatPtrRaw A hA rowNext) (hgcdMatPtrRaw B hB lowerCol)
          T scratch
          (hgcdMatLenRaw A hA rowBase) (hgcdMatLenRaw B hB col)
          (hgcdMatLenRaw A hA rowNext) (hgcdMatLenRaw B hB lowerCol)
          heap with
      | .error fault => .error fault
      | .ok entry =>
        let nextLen := C.len.set i entry.length
          (by rw [hC.2]; exact hi)
        let next : HgcdMat := { C with len := nextLen }
        have hNext : next.Valid := by
          exact ⟨hC.1, by simp [next, nextLen, hC.2]⟩
        hgcdMatMulLoop this A B hA hB T scratch next hNext (i + 1)
          entry.heap
    else
      .ok { heap := heap, matrix := C }
termination_by C hC i heap => 4 - i
decreasing_by omega

def hgcdMatMul (this : DenseUPolyZp)
    (C A B : HgcdMat) (hC : C.Valid) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) (heap : RawHeap) : RawExec HgcdMatMulResult :=
  hgcdMatMulLoop this A B hA hB T scratch C hC 0 heap

theorem hgcdMatMulLoop_result_valid (this : DenseUPolyZp)
    (A B : HgcdMat) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) (C : HgcdMat) (hC : C.Valid)
    (i : Nat) (heap : RawHeap) (result : HgcdMatMulResult)
    (hrun : hgcdMatMulLoop this A B hA hB T scratch C hC i heap =
      .ok result) : result.matrix.Valid := by
  rw [hgcdMatMulLoop] at hrun
  split at hrun
  next hi =>
    dsimp only at hrun
    split at hrun
    next fault hentry => simp at hrun
    next entry hentry =>
      let nextLen := C.len.set i entry.length
        (by rw [hC.2]; exact hi)
      let next : HgcdMat := { C with len := nextLen }
      have hNext : next.Valid := by
        exact ⟨hC.1, by simp [next, nextLen, hC.2]⟩
      exact hgcdMatMulLoop_result_valid this A B hA hB T scratch next hNext
        (i + 1) entry.heap result hrun
  next hstop =>
    have heq := Except.ok.inj hrun
    subst result
    exact hC
termination_by 4 - i
decreasing_by omega

theorem hgcdMatMul_result_valid (this : DenseUPolyZp)
    (C A B : HgcdMat) (hC : C.Valid) (hA : A.Valid) (hB : B.Valid)
    (T scratch : RawPtr UInt64) (heap : RawHeap) (result : HgcdMatMulResult)
    (hrun : hgcdMatMul this C A B hC hA hB T scratch heap = .ok result) :
    result.matrix.Valid := by
  exact hgcdMatMulLoop_result_valid this A B hA hB T scratch C hC 0 heap
    result hrun

/-- Result of one source column update in
`[[q,1],[1,0]] * S`: `top := top + q * bottom`. -/
structure HgcdMatQuotientEntryResult where
  heap : RawHeap
  matrix : HgcdMat
  valid : matrix.Valid

/-- Exact lowering of one of the two guarded quotient updates after the
source has swapped the rows of `S`.  The multiplication operand order and
the in-place `_poly_add` are retained verbatim. -/
def hgcdMatQuotientEntry (this : DenseUPolyZp)
    (S : HgcdMat) (hS : S.Valid) (top bottom : Fin 4)
    (q : RawPtr UInt64) (lenQ : Nat) (T scratch : RawPtr UInt64)
    (heap : RawHeap) : RawExec HgcdMatQuotientEntryResult :=
  let bottomLen := hgcdMatLenRaw S hS bottom
  if lenQ > 0 && bottomLen > 0 then
    let left := if lenQ ≥ bottomLen then q else hgcdMatPtrRaw S hS bottom
    let leftLen := if lenQ ≥ bottomLen then lenQ else bottomLen
    let right := if lenQ ≥ bottomLen then hgcdMatPtrRaw S hS bottom else q
    let rightLen := if lenQ ≥ bottomLen then bottomLen else lenQ
    match dense_upoly_zp__mul_ir this T left leftLen right rightLen scratch
        heap with
    | .error fault => .error fault
    | .ok heap1 =>
      let lenT := lenQ + bottomLen - 1
      match dense_upoly_zp__poly_add_ir this (hgcdMatPtrRaw S hS top)
          (hgcdMatPtrRaw S hS top) (hgcdMatLenRaw S hS top) T lenT
          heap1 with
      | .error fault => .error fault
      | .ok (heap2, sumLen) =>
        let len' := S.len.set top.val sumLen (by rw [hS.2]; exact top.isLt)
        let next : HgcdMat := { S with len := len' }
        have hNext : next.Valid := ⟨hS.1, by simp [next, len', hS.2]⟩
        .ok { heap := heap2, matrix := next, valid := hNext }
  else
    .ok { heap := heap, matrix := S, valid := hS }

theorem hgcdMatQuotientEntry_result_valid (this : DenseUPolyZp)
    (S : HgcdMat) (hS : S.Valid) (top bottom : Fin 4)
    (q : RawPtr UInt64) (lenQ : Nat) (T scratch : RawPtr UInt64)
    (heap : RawHeap) (result : HgcdMatQuotientEntryResult)
    (hrun : hgcdMatQuotientEntry this S hS top bottom q lenQ T scratch
      heap = .ok result) : result.matrix.Valid := by
  exact result.valid

/-- Descriptor-only effect of the four source `std::swap` calls. -/
def hgcdMatSwapRows (S : HgcdMat) (hS : S.Valid) : HgcdMat :=
  { poly := #[hgcdMatPtrRaw S hS (2 : Fin 4),
      hgcdMatPtrRaw S hS (3 : Fin 4),
      hgcdMatPtrRaw S hS (0 : Fin 4),
      hgcdMatPtrRaw S hS (1 : Fin 4)]
    len := #[hgcdMatLenRaw S hS (2 : Fin 4),
      hgcdMatLenRaw S hS (3 : Fin 4),
      hgcdMatLenRaw S hS (0 : Fin 4),
      hgcdMatLenRaw S hS (1 : Fin 4)] }

theorem hgcdMatSwapRows_valid (S : HgcdMat) (hS : S.Valid) :
    (hgcdMatSwapRows S hS).Valid := by
  simp [hgcdMatSwapRows, HgcdMat.Valid]

/-- Result of the exact source block constructing
`S_modified = [[q,1],[1,0]] * S`. -/
structure HgcdMatQuotientResult where
  heap : RawHeap
  matrix : HgcdMat
  valid : matrix.Valid

/-- Exact source-order composition: swap rows, update column zero, then
update column one, reusing `T` and `scratch` exactly as C++ does. -/
def hgcdMatApplyQuotient (this : DenseUPolyZp)
    (S : HgcdMat) (hS : S.Valid) (q : RawPtr UInt64) (lenQ : Nat)
    (T scratch : RawPtr UInt64) (heap : RawHeap) :
    RawExec HgcdMatQuotientResult :=
  let swapped := hgcdMatSwapRows S hS
  have hSwapped : swapped.Valid := hgcdMatSwapRows_valid S hS
  match hgcdMatQuotientEntry this swapped hSwapped (0 : Fin 4) (2 : Fin 4)
      q lenQ T scratch heap with
  | .error fault => .error fault
  | .ok first =>
    match hgcdMatQuotientEntry this first.matrix first.valid (1 : Fin 4)
        (3 : Fin 4) q lenQ T scratch first.heap with
    | .error fault => .error fault
    | .ok second => .ok {
        heap := second.heap, matrix := second.matrix, valid := second.valid }

/-- Successful quotient application exposes the two real update executions
and their generated validity witnesses. -/
theorem hgcdMatApplyQuotient_exec (this : DenseUPolyZp)
    (S : HgcdMat) (hS : S.Valid) (q : RawPtr UInt64) (lenQ : Nat)
    (T scratch : RawPtr UInt64) (heap : RawHeap)
    (result : HgcdMatQuotientResult)
    (hrun : hgcdMatApplyQuotient this S hS q lenQ T scratch heap =
      .ok result) :
    let swapped := hgcdMatSwapRows S hS
    let hSwapped := hgcdMatSwapRows_valid S hS
    ∃ first,
      hgcdMatQuotientEntry this swapped hSwapped (0 : Fin 4) (2 : Fin 4)
          q lenQ T scratch heap = .ok first ∧
      hgcdMatQuotientEntry this first.matrix first.valid (1 : Fin 4) (3 : Fin 4)
          q lenQ T scratch first.heap =
        .ok (HgcdMatQuotientEntryResult.mk result.heap result.matrix
          result.valid) := by
  simp only [hgcdMatApplyQuotient] at hrun
  split at hrun
  next fault hfirst => simp at hrun
  next first hfirst =>
    split at hrun
    next fault hsecond => simp at hrun
    next second hsecond =>
      have heq := Except.ok.inj hrun
      subst result
      exact ⟨first, hfirst, hsecond⟩

theorem hgcdMatApplyQuotient_result_valid (this : DenseUPolyZp)
    (S : HgcdMat) (hS : S.Valid) (q : RawPtr UInt64) (lenQ : Nat)
    (T scratch : RawPtr UInt64) (heap : RawHeap)
    (result : HgcdMatQuotientResult)
    (hrun : hgcdMatApplyQuotient this S hS q lenQ T scratch heap =
      .ok result) : result.matrix.Valid := by
  exact result.valid

/-- Exact lowering of the final matrix block of `_hgcd_recursive`:
first form `S := [[q,1],[1,0]] * S` using the two real guarded column
updates, then execute the complete generated `_mat_mul(M, R, S, a2,
scratch)`. -/
def hgcdRecursiveCombineMatrix (this : DenseUPolyZp)
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) : RawExec HgcdMatMulResult :=
  match hgcdMatApplyQuotient this S hS q lenQ T scratch heap with
  | .error fault => .error fault
  | .ok modified =>
    hgcdMatMul this M R modified.matrix hM hR modified.valid a2 scratch
      modified.heap

/-- A successful final matrix block exposes both actual generated calls in
their C++ order. -/
theorem hgcdRecursiveCombineMatrix_exec (this : DenseUPolyZp)
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) (result : HgcdMatMulResult)
    (hrun : hgcdRecursiveCombineMatrix this M R S hM hR hS q lenQ T a2
      scratch heap = .ok result) :
    ∃ modified,
      hgcdMatApplyQuotient this S hS q lenQ T scratch heap = .ok modified ∧
      hgcdMatMul this M R modified.matrix hM hR modified.valid a2 scratch
        modified.heap = .ok result := by
  simp only [hgcdRecursiveCombineMatrix] at hrun
  split at hrun
  next fault hmodified => simp at hrun
  next modified hmodified => exact ⟨modified, hmodified, hrun⟩

theorem hgcdRecursiveCombineMatrix_result_valid (this : DenseUPolyZp)
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (q : RawPtr UInt64) (lenQ : Nat) (T a2 scratch : RawPtr UInt64)
    (heap : RawHeap) (result : HgcdMatMulResult)
    (hrun : hgcdRecursiveCombineMatrix this M R S hM hR hS q lenQ T a2
      scratch heap = .ok result) : result.matrix.Valid := by
  rcases hgcdRecursiveCombineMatrix_exec this M R S hM hR hS q lenQ T a2
      scratch heap result hrun with ⟨modified, _, hmul⟩
  exact hgcdMatMul_result_valid this M R modified.matrix hM hR modified.valid
    a2 scratch modified.heap result hmul

/-- Return state of the complete tail after the second high-half HGCD call. -/
structure HgcdRecursiveFinishResult where
  heap : RawHeap
  matrix : HgcdMat
  valid : matrix.Valid
  lenA : Nat
  lenB : Nat
  sgn : Int

/-- Exact source-order tail of `_hgcd_recursive`: reconstruct `B`, reconstruct
`A`, optionally execute the quotient update and full matrix product, then
return `-(sgnR * sgnS)`. -/
def hgcdRecursiveFinish (this : DenseUPolyZp)
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (computeM : Bool) (A B T0 lowA lowB highA highB q : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift lenQ : Nat)
    (a2 scratch : RawPtr UInt64) (sgnR sgnS : Int) (heap : RawHeap) :
    RawExec HgcdRecursiveFinishResult :=
  match hgcdRecursiveReconstructPair this A B T0 lowA lowB highA highB
      scratch lenLowA lenLowB lenHighA lenHighB shift S hS sgnS heap with
  | .error fault => .error fault
  | .ok reconstructed =>
    if computeM then
      match hcombine : hgcdRecursiveCombineMatrix this M R S hM hR hS q
          lenQ T0 a2 scratch reconstructed.heap with
      | .error fault => .error fault
      | .ok combined => .ok {
          heap := combined.heap
          matrix := combined.matrix
          valid := hgcdRecursiveCombineMatrix_result_valid this M R S hM hR
            hS q lenQ T0 a2 scratch reconstructed.heap combined hcombine
          lenA := reconstructed.lenA
          lenB := reconstructed.lenB
          sgn := -(sgnR * sgnS) }
    else .ok {
      heap := reconstructed.heap
      matrix := M
      valid := hM
      lenA := reconstructed.lenA
      lenB := reconstructed.lenB
      sgn := -(sgnR * sgnS) }

/-- Successful execution of the recursive tail exposes the actual pair
reconstruction and, exactly when requested by C++, the actual matrix block. -/
theorem hgcdRecursiveFinish_exec (this : DenseUPolyZp)
    (M R S : HgcdMat) (hM : M.Valid) (hR : R.Valid) (hS : S.Valid)
    (computeM : Bool) (A B T0 lowA lowB highA highB q : RawPtr UInt64)
    (lenLowA lenLowB lenHighA lenHighB shift lenQ : Nat)
    (a2 scratch : RawPtr UInt64) (sgnR sgnS : Int) (heap : RawHeap)
    (result : HgcdRecursiveFinishResult)
    (hrun : hgcdRecursiveFinish this M R S hM hR hS computeM A B T0 lowA
      lowB highA highB q lenLowA lenLowB lenHighA lenHighB shift lenQ a2
      scratch sgnR sgnS heap = .ok result) :
    ∃ reconstructed,
      hgcdRecursiveReconstructPair this A B T0 lowA lowB highA highB scratch
          lenLowA lenLowB lenHighA lenHighB shift S hS sgnS heap =
        .ok reconstructed ∧
      result.lenA = reconstructed.lenA ∧
      result.lenB = reconstructed.lenB ∧
      result.sgn = -(sgnR * sgnS) ∧
      (if computeM then
        ∃ combined,
          hgcdRecursiveCombineMatrix this M R S hM hR hS q lenQ T0 a2
              scratch reconstructed.heap = .ok combined ∧
          result.heap = combined.heap ∧ result.matrix = combined.matrix
       else result.heap = reconstructed.heap ∧ result.matrix = M) := by
  simp only [hgcdRecursiveFinish] at hrun
  split at hrun
  next fault hreconstruct => simp at hrun
  next reconstructed hreconstruct =>
    split at hrun
    next hcompute =>
      split at hrun
      next fault hcombine => simp at hrun
      next combined hcombine =>
        have heq := Except.ok.inj hrun
        subst result
        exact ⟨reconstructed, hreconstruct, rfl, rfl, rfl, by
          simp [hcompute, hcombine]⟩
    next hcompute =>
      have heq := Except.ok.inj hrun
      subst result
      exact ⟨reconstructed, hreconstruct, rfl, rfl, rfl, by
        simp [hcompute]⟩

structure HgcdEarlyMatrixResult where
  heap : RawHeap
  matrix : HgcdMat

/-- Exact source-order loop used by the early-return branch to copy `R` into
the caller's matrix buffers and install the four returned lengths. -/
def hgcdEarlyMatrixLoop (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid)
    (i : Nat) (heap : RawHeap) : RawExec HgcdEarlyMatrixResult :=
  if hi : i < 4 then
    let index : Fin 4 := ⟨i, hi⟩
    match heap.copyU64 (hgcdMatPtrRaw M hM index)
        (hgcdMatPtrRaw R hR index) (hgcdMatLenRaw R hR index) with
    | .error fault => .error fault
    | .ok heap1 =>
      let nextLen := M.len.set i (hgcdMatLenRaw R hR index)
        (by rw [hM.2]; exact hi)
      let next : HgcdMat := { M with len := nextLen }
      have hNext : next.Valid := by
        exact ⟨hM.1, by simp [next, nextLen, hM.2]⟩
      hgcdEarlyMatrixLoop next R hNext hR (i + 1) heap1
  else
    .ok { heap := heap, matrix := M }
termination_by 4 - i
decreasing_by omega

theorem hgcdEarlyMatrixLoop_result_valid (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid) (i : Nat) (heap : RawHeap)
    (result : HgcdEarlyMatrixResult)
    (hrun : hgcdEarlyMatrixLoop M R hM hR i heap = .ok result) :
    result.matrix.Valid := by
  rw [hgcdEarlyMatrixLoop] at hrun
  split at hrun
  next hi =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      exact hgcdEarlyMatrixLoop_result_valid _ R _ hR (i + 1)
        heap1 result hrun
  next hi =>
    have heq : result = HgcdEarlyMatrixResult.mk heap M :=
      (Except.ok.inj hrun).symm
    subst result
    exact hM
termination_by 4 - i
decreasing_by omega

/-- Descriptor invariant of the real early matrix-copy loop.  At entry `i`,
all earlier slots already carry `R`'s lengths; a successful suffix finishes
all four slots. -/
theorem hgcdEarlyMatrixLoop_lengths (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid) (i : Nat) (heap : RawHeap)
    (result : HgcdEarlyMatrixResult)
    (hdone : ∀ j : Fin 4, j.val < i →
      hgcdMatLenRaw M hM j = hgcdMatLenRaw R hR j)
    (hrun : hgcdEarlyMatrixLoop M R hM hR i heap = .ok result) :
    result.matrix.len = R.len := by
  rw [hgcdEarlyMatrixLoop] at hrun
  split at hrun
  next hi =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hi⟩
      let nextLen := M.len.set i (hgcdMatLenRaw R hR index)
        (by rw [hM.2]; exact hi)
      let next : HgcdMat := { M with len := nextLen }
      have hNext : next.Valid := by
        exact ⟨hM.1, by simp [next, nextLen, hM.2]⟩
      have hdoneNext : ∀ j : Fin 4, j.val < i + 1 →
          hgcdMatLenRaw next hNext j = hgcdMatLenRaw R hR j := by
        intro j hj
        by_cases hji : j.val = i
        · subst i
          simp [hgcdMatLenRaw, next, nextLen, index, hM.2]
        · have hold := hdone j (by omega)
          simp only [hgcdMatLenRaw, next, nextLen]
          rw [Array.getElem_set_ne
            (by simpa [hM.2] using hi)
            (by simpa [hM.2] using j.isLt)
            (by exact Ne.symm hji)]
          exact hold
      exact hgcdEarlyMatrixLoop_lengths next R hNext hR (i + 1) heap1
        result hdoneNext hrun
  next hi =>
    have heq : result = HgcdEarlyMatrixResult.mk heap M :=
      (Except.ok.inj hrun).symm
    subst result
    apply Array.ext
    · rw [hM.2, hR.2]
    · intro k hkM hkR
      have hk4 : k < 4 := by simpa [hR.2] using hkR
      let j : Fin 4 := ⟨k, hk4⟩
      have hj := hdone j (by omega)
      simpa [hgcdMatLenRaw, j] using hj
termination_by 4 - i
decreasing_by omega

structure HgcdRecursiveEarlyResult where
  heap : RawHeap
  matrix : HgcdMat
  lenA : Nat
  lenB : Nat
  sgn : Int

/-- Exact lowering of the `_hgcd_recursive` early return after both shifted
reconstructions.  Output copies precede the optional matrix loop. -/
def hgcdRecursiveEarlyReturn (M R : HgcdMat)
    (hM : M.Valid) (hR : R.Valid) (computeM : Bool)
    (A B a2 b2 : RawPtr UInt64) (lenA2 lenB2 : Nat) (sgn : Int)
    (heap : RawHeap) : RawExec HgcdRecursiveEarlyResult :=
  match heap.copyU64 A a2 lenA2 with
  | .error fault => .error fault
  | .ok heap1 =>
    match heap1.copyU64 B b2 lenB2 with
    | .error fault => .error fault
    | .ok heap2 =>
      if computeM then
        match hgcdEarlyMatrixLoop M R hM hR 0 heap2 with
        | .error fault => .error fault
        | .ok matrixResult => .ok {
            heap := matrixResult.heap
            matrix := matrixResult.matrix
            lenA := lenA2
            lenB := lenB2
            sgn := sgn }
      else .ok {
        heap := heap2
        matrix := M
        lenA := lenA2
        lenB := lenB2
        sgn := sgn }

/-- State produced by the real middle divrem and the source pointer arithmetic
for the second high-half HGCD call. -/
structure HgcdRecursiveMiddleResult where
  heap : RawHeap
  lenQ : Nat
  lenD : Nat
  k : Nat
  c0 : RawPtr UInt64
  lenC0 : Nat
  d0 : RawPtr UInt64
  lenD0 : Nat

/-- Exact lowering of `_hgcd_recursive` from its failed early-stop test through
the middle `_poly_divrem` and construction of `c0`/`d0`. -/
def hgcdRecursiveMiddle (this : DenseUPolyZp)
    (q d a2 b2 : RawPtr UInt64) (lenA2 lenB2 m : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap) :
    RawExec HgcdRecursiveMiddleResult :=
  match dense_upoly_zp__poly_divrem_ir this q d a2 lenA2 b2 lenB2 W3
      heap with
  | .error fault => .error fault
  | .ok (heap1, lenQ, lenD) =>
    let k := 2 * m - lenB2 + 1
    .ok {
      heap := heap1
      lenQ := lenQ
      lenD := lenD
      k := k
      c0 := b2.add k
      lenC0 := if lenB2 ≥ k then lenB2 - k else 0
      d0 := d.add k
      lenD0 := if lenD ≥ k then lenD - k else 0 }

theorem hgcdRecursiveMiddle_layout (this : DenseUPolyZp)
    (q d a2 b2 : RawPtr UInt64) (lenA2 lenB2 m : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (result : HgcdRecursiveMiddleResult)
    (hrun : hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3 heap =
      .ok result) :
    dense_upoly_zp__poly_divrem_ir this q d a2 lenA2 b2 lenB2 W3 heap =
        .ok (result.heap, result.lenQ, result.lenD) ∧
      result.k = 2 * m - lenB2 + 1 ∧
      result.c0 = b2.add result.k ∧
      result.lenC0 = (if lenB2 ≥ result.k then lenB2 - result.k else 0) ∧
      result.d0 = d.add result.k ∧
      result.lenD0 = (if result.lenD ≥ result.k then
        result.lenD - result.k else 0) := by
  simp only [hgcdRecursiveMiddle] at hrun
  generalize hdiv : dense_upoly_zp__poly_divrem_ir this q d a2 lenA2 b2
    lenB2 W3 heap = divResult at hrun
  cases divResult with
  | error fault => simp at hrun
  | ok value =>
    rcases value with ⟨heap1, lenQ, lenD⟩
    have heq := Except.ok.inj hrun
    subst result
    simp [hdiv]

/-- The second recursive call also decreases the enclosing `len_a` measure.
The source always chooses `k = 2*m-len_b2+1`, hence a nonempty `c0` is a
strict suffix of `b2`; the reconstruction bound then places it below `len_a`. -/
theorem hgcdRecursiveMiddle_lenC0_lt (this : DenseUPolyZp)
    (q d a2 b2 : RawPtr UInt64) (lenA2 lenB2 m lenA : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (result : HgcdRecursiveMiddleResult)
    (hlenA : 0 < lenA) (hlenB2 : lenB2 ≤ lenA)
    (hrun : hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3 heap =
      .ok result) :
    result.lenC0 < lenA := by
  have hlayout := hgcdRecursiveMiddle_layout this q d a2 b2 lenA2 lenB2 m
    W3 heap result hrun
  rw [hlayout.2.2.2.1]
  split
  next hk =>
    have hkPos : 0 < result.k := by
      rw [hlayout.2.1]
      omega
    omega
  next hk =>
    omega

/-- Failure of the source's early-return guard makes the second-call
divisor suffix nonempty.  This follows from the exact source choice
`k = 2*m-lenB2+1`, including truncated natural subtraction. -/
theorem hgcdRecursiveMiddle_lenC0_pos (this : DenseUPolyZp)
    (q d a2 b2 : RawPtr UInt64) (lenA2 lenB2 m : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (result : HgcdRecursiveMiddleResult)
    (hm : 0 < m)
    (hnonearly : m + 1 ≤ lenB2)
    (hrun : hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3 heap =
      .ok result) :
    0 < result.lenC0 := by
  have hlayout := hgcdRecursiveMiddle_layout this q d a2 b2 lenA2 lenB2 m
    W3 heap result hrun
  rw [hlayout.2.2.2.1, hlayout.2.1]
  split
  next hk => omega
  next hk => omega

/-- Shifting divisor and remainder by the same source offset preserves their
strict length order whenever the divisor suffix is nonempty.  This is the
ordering required by the second recursive HGCD branch. -/
theorem hgcdRecursiveMiddle_lenD0_lt_lenC0 (this : DenseUPolyZp)
    (q d a2 b2 : RawPtr UInt64) (lenA2 lenB2 m : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (result : HgcdRecursiveMiddleResult)
    (hlenD : result.lenD < lenB2) (hlenC0 : 0 < result.lenC0)
    (hrun : hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3 heap =
      .ok result) :
    result.lenD0 < result.lenC0 := by
  have hlayout := hgcdRecursiveMiddle_layout this q d a2 b2 lenA2 lenB2 m
    W3 heap result hrun
  have hk := hlayout.2.1
  have hc := hlayout.2.2.2.1
  have hd := hlayout.2.2.2.2.2
  rw [hk] at hc hd
  rw [hc] at hlenC0
  rw [hd, hc]
  split at hlenC0
  next hkB =>
    split
    next hkD => omega
    next hkD => omega
  next hkB =>
    simp at hlenC0

/-- All measure facts required at the second recursive call, tied to the
same concrete middle divrem execution.  The only incoming algorithmic facts
are the first reconstruction bound and the real remainder strict decrease. -/
theorem hgcdRecursiveMiddle_second_call_bounds (this : DenseUPolyZp)
    (q d a2 b2 : RawPtr UInt64) (lenA2 lenB2 m outerLenA : Nat)
    (W3 : RawPtr Word3) (heap : RawHeap)
    (result : HgcdRecursiveMiddleResult)
    (houterPos : 0 < outerLenA) (hm : 0 < m)
    (hreconstructed : lenB2 ≤ outerLenA)
    (hnonearly : m + 1 ≤ lenB2)
    (hremainder : result.lenD < lenB2)
    (hrun : hgcdRecursiveMiddle this q d a2 b2 lenA2 lenB2 m W3 heap =
      .ok result) :
    0 < result.lenC0 ∧ result.lenD0 < result.lenC0 ∧
      result.lenC0 < outerLenA := by
  have hcpos := hgcdRecursiveMiddle_lenC0_pos this q d a2 b2 lenA2 lenB2
    m W3 heap result hm hnonearly hrun
  exact ⟨hcpos,
    hgcdRecursiveMiddle_lenD0_lt_lenC0 this q d a2 b2 lenA2 lenB2 m W3
      heap result hremainder hcpos hrun,
    hgcdRecursiveMiddle_lenC0_lt this q d a2 b2 lenA2 lenB2 m outerLenA
      W3 heap result houterPos hreconstructed hrun⟩

/-- Every successful suffix of the real restore loop returns a valid matrix
and leaves the complete length descriptor byte-for-byte unchanged. -/
theorem hgcdMatRestoreLoop_preserves_valid_len
    (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64) (i off : Nat) (heap : RawHeap)
    (result : HgcdMatRestoreResult)
    (hrun : hgcdMatRestoreLoop original current hOriginal hCurrent stage
      i off heap = .ok result) :
    result.matrix.Valid ∧ result.matrix.len = current.len := by
  rw [hgcdMatRestoreLoop] at hrun
  split at hrun
  next hi =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hi⟩
      let poly' := current.poly.set i (hgcdMatPtrRaw original hOriginal index)
        (by rw [hCurrent.1]; omega)
      let next : HgcdMat := { current with poly := poly' }
      have hNext : next.Valid := by
        exact ⟨by simp [next, poly', hCurrent.1], hCurrent.2⟩
      have hrec := hgcdMatRestoreLoop_preserves_valid_len original next
        hOriginal hNext stage (i + 1)
        (off + hgcdMatLenRaw current hCurrent index) heap1 result hrun
      exact ⟨hrec.1, hrec.2.trans (by rfl)⟩
  next hstop =>
    have heq : ({ heap := heap, matrix := current, off := off } :
        HgcdMatRestoreResult) = result := Except.ok.inj hrun
    subst result
    exact ⟨hCurrent, rfl⟩
termination_by 4 - i
decreasing_by omega

/-- Pure descriptor effect of the generated restore loop.  Keeping it separate
from the heap operation makes the exact four pointer writes explicit without
assigning any polynomial meaning to the copied bytes. -/
def hgcdMatRestorePointers (original : HgcdMat) (hOriginal : original.Valid) :
    (poly : Array (RawPtr UInt64)) → (hPoly : poly.size = 4) →
      (i : Nat) → Array (RawPtr UInt64)
  | poly, hPoly, i =>
    if hi : i < 4 then
      hgcdMatRestorePointers original hOriginal
        (poly.set i (hgcdMatPtrRaw original hOriginal ⟨i, hi⟩)
          (by rw [hPoly]; exact hi))
        (by simp [hPoly]) (i + 1)
    else
      poly
termination_by poly hPoly i => 4 - i
decreasing_by omega

theorem hgcdMatRestorePointers_size (original : HgcdMat)
    (hOriginal : original.Valid) (poly : Array (RawPtr UInt64))
    (hPoly : poly.size = 4) (i : Nat) :
    (hgcdMatRestorePointers original hOriginal poly hPoly i).size = poly.size := by
  rw [hgcdMatRestorePointers]
  split
  next hi =>
    rw [hgcdMatRestorePointers_size]
    simp
  next hstop => rfl
termination_by 4 - i
decreasing_by omega

/-- A successful generated restore loop has exactly the descriptor-pointer
effect described by `hgcdMatRestorePointers`. -/
theorem hgcdMatRestoreLoop_poly_eq
    (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64) (i off : Nat) (heap : RawHeap)
    (result : HgcdMatRestoreResult)
    (hrun : hgcdMatRestoreLoop original current hOriginal hCurrent stage
      i off heap = .ok result) :
    result.matrix.poly = hgcdMatRestorePointers original hOriginal current.poly
      hCurrent.1 i := by
  rw [hgcdMatRestoreLoop] at hrun
  split at hrun
  next hi =>
    dsimp only at hrun
    split at hrun
    next fault hcopy => simp at hrun
    next heap1 hcopy =>
      let index : Fin 4 := ⟨i, hi⟩
      let poly' := current.poly.set i (hgcdMatPtrRaw original hOriginal index)
        (by rw [hCurrent.1]; omega)
      let next : HgcdMat := { current with poly := poly' }
      have hNext : next.Valid := by
        exact ⟨by simp [next, poly', hCurrent.1], hCurrent.2⟩
      have hrec := hgcdMatRestoreLoop_poly_eq original next hOriginal hNext
        stage (i + 1) (off + hgcdMatLenRaw current hCurrent index)
        heap1 result hrun
      rw [hgcdMatRestorePointers]
      simp only [hi, ↓reduceDIte]
      exact hrec
  next hstop =>
    have heq : ({ heap := heap, matrix := current, off := off } :
        HgcdMatRestoreResult) = result := Except.ok.inj hrun
    subst result
    rw [hgcdMatRestorePointers]
    simp only [hstop, ↓reduceDIte]
termination_by 4 - i
decreasing_by omega

/-- At the generated entry index, the pure descriptor effect is the four
source-order pointer writes. -/
theorem hgcdMatRestorePointers_zero_effect (original : HgcdMat)
    (hOriginal : original.Valid) (poly : Array (RawPtr UInt64))
    (hPoly : poly.size = 4) :
    hgcdMatRestorePointers original hOriginal poly hPoly 0 =
      (((poly.set 0 (hgcdMatPtrRaw original hOriginal ⟨0, by omega⟩)
          (by omega)).set 1
          (hgcdMatPtrRaw original hOriginal ⟨1, by omega⟩) (by simp [hPoly])).set 2
          (hgcdMatPtrRaw original hOriginal ⟨2, by omega⟩) (by simp [hPoly])).set 3
          (hgcdMatPtrRaw original hOriginal ⟨3, by omega⟩) (by simp [hPoly]) := by
  simp [hgcdMatRestorePointers]

/-- Consequently all four returned descriptors use the exact pointers saved
before the iterator call. -/
theorem hgcdMatRestorePointers_zero (original : HgcdMat)
    (hOriginal : original.Valid) (poly : Array (RawPtr UInt64))
    (hPoly : poly.size = 4) :
    hgcdMatRestorePointers original hOriginal poly hPoly 0 = original.poly := by
  rw [hgcdMatRestorePointers_zero_effect]
  apply Array.ext
  · simp [hPoly, hOriginal.1]
  · intro i hleft hright
    have hi : i < 4 := by simpa [hOriginal.1] using hright
    have hcases : i = 0 ∨ i = 1 ∨ i = 2 ∨ i = 3 := by omega
    rcases hcases with rfl | rfl | rfl | rfl <;>
      simp only [Array.getElem_set, hgcdMatPtrRaw] <;> simp

/-- The complete second generated loop restores all saved pointers while
retaining the iterator-produced lengths. -/
theorem hgcdMatRestoreLoop_zero_descriptors
    (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64) (heap : RawHeap)
    (result : HgcdMatRestoreResult)
    (hrun : hgcdMatRestoreLoop original current hOriginal hCurrent stage
      0 0 heap = .ok result) :
    result.matrix.Valid ∧ result.matrix.poly = original.poly ∧
      result.matrix.len = current.len := by
  have hvalidLen := hgcdMatRestoreLoop_preserves_valid_len original current
    hOriginal hCurrent stage 0 0 heap result hrun
  have hpoly := hgcdMatRestoreLoop_poly_eq original current hOriginal hCurrent
    stage 0 0 heap result hrun
  exact ⟨hvalidLen.1,
    hpoly.trans (hgcdMatRestorePointers_zero original hOriginal current.poly hCurrent.1),
    hvalidLen.2⟩

/-- Thus the exact two-loop source stabilization block restores the pre-call
pointers and preserves the post-call lengths on every successful execution. -/
theorem hgcdMatStabilize_preserves_descriptors
    (original current : HgcdMat)
    (hOriginal : original.Valid) (hCurrent : current.Valid)
    (stage : RawPtr UInt64) (heap : RawHeap)
    (result : HgcdMatRestoreResult)
    (hrun : hgcdMatStabilize original current hOriginal hCurrent stage heap =
      .ok result) :
    result.matrix.Valid ∧ result.matrix.poly = original.poly ∧
      result.matrix.len = current.len := by
  rw [hgcdMatStabilize] at hrun
  split at hrun
  next fault hstage => simp at hrun
  next staged hstage =>
    exact hgcdMatRestoreLoop_zero_descriptors original current hOriginal hCurrent
      stage staged.heap result hrun

theorem hgcdRecursiveIterBranch_result_valid (this : DenseUPolyZp)
    (original : HgcdMat) (hOriginal : original.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage : RawPtr UInt64) (heap : RawHeap)
    (result : HgcdRecursiveIterBranchResult)
    (hrun : hgcdRecursiveIterBranch this original hOriginal a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap =
        .ok result) :
    result.matrix.Valid := by
  rcases hgcdRecursiveIterBranch_exec this original hOriginal a3 b3 inputA
      inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap result hrun with
    ⟨iter, hIter, stable, _, hstable, _, hmatrix, _, _, _⟩
  rw [hmatrix]
  exact (hgcdMatStabilize_preserves_descriptors original iter.matrix
    hOriginal hIter stage iter.heap stable hstable).1

theorem hgcdRecursiveEarlyReturn_result_valid
    (M R : HgcdMat) (hM : M.Valid) (hR : R.Valid) (computeM : Bool)
    (A B a2 b2 : RawPtr UInt64) (lenA2 lenB2 : Nat) (sgn : Int)
    (heap : RawHeap) (result : HgcdRecursiveEarlyResult)
    (hrun : hgcdRecursiveEarlyReturn M R hM hR computeM A B a2 b2 lenA2
      lenB2 sgn heap = .ok result) :
    result.matrix.Valid := by
  simp only [hgcdRecursiveEarlyReturn] at hrun
  split at hrun
  next fault hcopyA => simp at hrun
  next heap1 hcopyA =>
    split at hrun
    next fault hcopyB => simp at hrun
    next heap2 hcopyB =>
      split at hrun
      next hcompute =>
        split at hrun
        next fault hmatrix => simp at hrun
        next matrixResult hmatrix =>
          have heq := Except.ok.inj hrun
          subst result
          exact hgcdEarlyMatrixLoop_result_valid M R hM hR 0 heap2
            matrixResult hmatrix
      next hcompute =>
        have heq := Except.ok.inj hrun
        subst result
        exact hM

/-- Common return record for the complete recursive lowering.  The validity
field is proof-only; all computational fields are copied byte-for-byte from
the concrete source branch that produced them. -/
structure HgcdRecursiveResult where
  heap : RawHeap
  matrix : HgcdMat
  valid : matrix.Valid
  lenA : Nat
  lenB : Nat
  sgn : Int

/-- Proof-erased observable return state.  Refinement and unfolding theorems
compare this record, while `HgcdRecursiveResult.valid` remains available to
feed the next concrete matrix call. -/
structure HgcdRecursiveValue where
  heap : RawHeap
  matrix : HgcdMat
  lenA : Nat
  lenB : Nat
  sgn : Int

def HgcdRecursiveResult.value (result : HgcdRecursiveResult) :
    HgcdRecursiveValue :=
  ⟨result.heap, result.matrix, result.lenA, result.lenB, result.sgn⟩

@[ext] theorem HgcdRecursiveResult.ext_value
    (left right : HgcdRecursiveResult) (hvalue : left.value = right.value) :
    left = right := by
  cases left with
  | mk leftHeap leftMatrix leftValid leftLenA leftLenB leftSgn =>
    cases right with
    | mk rightHeap rightMatrix rightValid rightLenA rightLenB rightSgn =>
      simp only [HgcdRecursiveResult.value, HgcdRecursiveValue.mk.injEq] at hvalue
      rcases hvalue with ⟨rfl, rfl, rfl, rfl, rfl⟩
      rfl

def HgcdRecursiveBaseResult.toResult (result : HgcdRecursiveBaseResult)
    (hvalid : result.matrix.Valid) : HgcdRecursiveResult :=
  ⟨result.heap, result.matrix, hvalid, result.lenA, result.lenB, result.sgn⟩

def HgcdRecursiveIterBranchResult.toResult
    (result : HgcdRecursiveIterBranchResult)
    (hvalid : result.matrix.Valid) : HgcdRecursiveResult :=
  ⟨result.heap, result.matrix, hvalid, result.lenA, result.lenB, result.sgn⟩

def HgcdRecursiveEarlyResult.toResult (result : HgcdRecursiveEarlyResult)
    (hvalid : result.matrix.Valid) : HgcdRecursiveResult :=
  ⟨result.heap, result.matrix, hvalid, result.lenA, result.lenB, result.sgn⟩

def HgcdRecursiveFinishResult.toResult (result : HgcdRecursiveFinishResult) :
    HgcdRecursiveResult :=
  ⟨result.heap, result.matrix, result.valid, result.lenA, result.lenB,
    result.sgn⟩

/-- Computational signature of one recursive self-call.  It is used only to
factor the source's identical cutoff dispatch at the two call sites; the
well-founded main definition supplies itself as this argument. -/
abbrev HgcdRecursiveCall :=
  (M : HgcdMat) → M.Valid → (computeM : Bool) →
  (A B a b : RawPtr UInt64) → (lenA lenB : Nat) →
  (W scratch : RawPtr UInt64) → RawHeap → RawExec HgcdRecursiveResult

/-- Recursive callback available while defining an invocation of measure
`bound`.  Its type makes every recursive call carry a strict `lenA` decrease;
there is no counter and no executable guard for this proof-only argument. -/
abbrev HgcdRecursiveCallBelow (bound : Nat) :=
  (M : HgcdMat) → M.Valid → (computeM : Bool) →
  (A B a b : RawPtr UInt64) → (lenA lenB : Nat) →
  (W scratch : RawPtr UInt64) → RawHeap → lenB < lenA → lenA < bound →
  RawExec HgcdRecursiveResult

def hgcdRecursiveCutoff : Nat := 100

theorem hgcdRecursiveWorkspace_R_valid (W : RawPtr UInt64) (lenA : Nat) :
    (hgcdRecursiveWorkspace W lenA).R.Valid := by
  simp [hgcdRecursiveWorkspace, HgcdMat.Valid]

theorem hgcdRecursiveWorkspace_S_valid (W : RawPtr UInt64) (lenA : Nat) :
    (hgcdRecursiveWorkspace W lenA).S.Valid := by
  simp [hgcdRecursiveWorkspace, HgcdMat.Valid]

/-- Exact lowering of either `len < HGCD_CUTOFF` dispatch.  The small arm
executes the generated iterator/stabilization/store block; the large arm
uses the supplied recursive call with `computeM=true`. -/
def hgcdRecursiveDispatch (this : DenseUPolyZp)
    (recurse : HgcdRecursiveCall)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap) :
    RawExec HgcdRecursiveResult :=
  if lenInputA < hgcdRecursiveCutoff then
    match hrun : hgcdRecursiveIterBranch this matrix hMatrix a3 b3 inputA
        inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap with
    | .error fault => .error fault
    | .ok result => .ok (result.toResult
        (hgcdRecursiveIterBranch_result_valid this matrix hMatrix a3 b3
          inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap
          result hrun))
  else
    recurse matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
      WNext scratch heap

/-- Strictly decreasing form of the same cutoff dispatch.  The small arm is
definitionally the identical iterator block; the large arm can invoke the
callback only with the supplied well-founded decrease proof. -/
def hgcdRecursiveDispatchBelow (this : DenseUPolyZp) (bound : Nat)
    (recurse : HgcdRecursiveCallBelow bound)
    (matrix : HgcdMat) (hMatrix : matrix.Valid)
    (a3 b3 inputA inputB : RawPtr UInt64) (lenInputA lenInputB : Nat)
    (Q : RawPtr UInt64) (W3 : RawPtr Word3)
    (T0 T1 scratch stage WNext : RawPtr UInt64) (heap : RawHeap)
    (horder : lenInputB < lenInputA) (hdecrease : lenInputA < bound) :
    RawExec HgcdRecursiveResult :=
  if lenInputA < hgcdRecursiveCutoff then
    match hrun : hgcdRecursiveIterBranch this matrix hMatrix a3 b3 inputA
        inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap with
    | .error fault => .error fault
    | .ok result => .ok (result.toResult
        (hgcdRecursiveIterBranch_result_valid this matrix hMatrix a3 b3
          inputA inputB lenInputA lenInputB Q W3 T0 T1 scratch stage heap
          result hrun))
  else
    recurse matrix hMatrix true a3 b3 inputA inputB lenInputA lenInputB
      WNext scratch heap horder hdecrease

/-- Length facts exported by either real iterator arm or a recursively
refined child.  They are exactly the facts consumed by reconstruction and
the next well-founded call. -/
structure HgcdRecursiveLengthInvariant (inputLength : Nat)
    (result : HgcdRecursiveResult) : Prop where
  row0A : hgcdMatLenRaw result.matrix result.valid (0 : Fin 4) + result.lenA ≤
    inputLength + 1
  row1B : hgcdMatLenRaw result.matrix result.valid (1 : Fin 4) + result.lenB ≤
    inputLength + 1
  row2A : hgcdMatLenRaw result.matrix result.valid (2 : Fin 4) + result.lenA ≤
    inputLength + 1
  row3B : hgcdMatLenRaw result.matrix result.valid (3 : Fin 4) + result.lenB ≤
    inputLength + 1
  row1A : hgcdMatLenRaw result.matrix result.valid (1 : Fin 4) + result.lenA ≤
    inputLength + 1
  row3A : hgcdMatLenRaw result.matrix result.valid (3 : Fin 4) + result.lenA ≤
    inputLength + 1
  inputBound : result.lenA ≤ inputLength
  stopped : result.lenB < inputLength / 2 + 1
  positive : 0 < result.lenA
  aboveHalf : inputLength / 2 < result.lenA
  coeffBound : ∀ i : Fin 4,
    hgcdMatLenRaw result.matrix result.valid i ≤
      inputLength - inputLength / 2

theorem hgcdRecursiveLengthInvariant_toResult_proof_irrel
    (inputLength : Nat) (result : HgcdRecursiveIterBranchResult)
    (hleft hright : result.matrix.Valid)
    (hinvariant : HgcdRecursiveLengthInvariant inputLength
      (result.toResult hleft)) :
    HgcdRecursiveLengthInvariant inputLength (result.toResult hright) := by
  have heq : hleft = hright := Subsingleton.elim _ _
  subst hright
  exact hinvariant

/-- Facts exported by the actual first paired reconstruction.  Besides the
strict outer decrease used by well-founded recursion, the exact leading-A
length and operand order are needed by the following generated divrem. -/
structure HgcdFirstReconstructionInvariant (outerLength : Nat)
    (first : HgcdRecursiveResult)
    (reconstructed : HgcdRecursiveReconstructPairResult) : Prop where
  leadingA : reconstructed.lenA = outerLength / 2 + first.lenA
  positiveA : 0 < reconstructed.lenA
  ordered : reconstructed.lenB ≤ reconstructed.lenA
  decreases : reconstructed.lenB < outerLength

/-- The successful result relation of the actual first cutoff dispatch. -/
def HgcdFirstDispatchResult (this : DenseUPolyZp) (bound : Nat)
    (recurse : HgcdRecursiveCallBelow bound)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat) (heap : RawHeap)
    (first : HgcdRecursiveResult) : Prop :=
  let ws := hgcdRecursiveWorkspace W lenA
  let high := hgcdRecursiveHighInput a b lenA lenB
  ∃ (hchildOrder : high.lenB0 < high.lenA0)
    (hchildDecrease : high.lenA0 < bound),
    hgcdRecursiveDispatchBelow this bound recurse ws.R
      (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
      high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
      heap hchildOrder hchildDecrease = .ok first

/-- Proof-only algorithmic invariant needed between the two recursive call
sites.  The result must come from the actual first dispatch; this prevents a
caller from supplying a bound for an unrelated, preselected record. -/
def HgcdFirstReconstructionBoundProvider (this : DenseUPolyZp)
    (a b W scratch : RawPtr UInt64) (lenA lenB : Nat)
    (actualFirst : HgcdRecursiveResult → Prop) : Prop :=
  let ws := hgcdRecursiveWorkspace W lenA
  ∀ (first : HgcdRecursiveResult)
    (reconstructed : HgcdRecursiveReconstructPairResult),
    actualFirst first →
    HgcdRecursiveLengthInvariant (lenA - lenA / 2) first →
    hgcdRecursiveReconstructPair this ws.a2 ws.b2 ws.T0 a b ws.a3 ws.b3
      scratch (Nat.min lenA (lenA / 2)) (Nat.min lenB (lenA / 2))
      first.lenA first.lenB (lenA / 2) first.matrix first.valid first.sgn
      first.heap = .ok reconstructed →
    HgcdFirstReconstructionInvariant lenA first reconstructed

/-- Strictly decreasing version of the complete recursive body.  Both
recursive dispatches receive a proof that their source `lenA` is smaller
than the current one; the proofs are erased and introduce no executable
counter or alternate branch. -/
def hgcdRecursiveBodyBelow (this : DenseUPolyZp) (bound : Nat)
    (recurse : HgcdRecursiveCallBelow bound)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap)
    (hbound : lenA = bound) (horder : lenB < lenA)
    (firstLength : ¬ lenB < lenA / 2 + 1 → ∀ first,
      let ws := hgcdRecursiveWorkspace W lenA
      let high := hgcdRecursiveHighInput a b lenA lenB
      ∀ (hchildOrder : high.lenB0 < high.lenA0)
        (hchildDecrease : high.lenA0 < bound),
      hgcdRecursiveDispatchBelow this bound recurse ws.R
        (hgcdRecursiveWorkspace_R_valid W lenA) ws.a3 ws.b3 high.a0 high.b0
        high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2 ws.next
        heap hchildOrder hchildDecrease = .ok first →
      HgcdRecursiveLengthInvariant high.lenA0 first)
    (reconstructionBound : ¬ lenB < lenA / 2 + 1 →
      HgcdFirstReconstructionBoundProvider this a b W
      scratch lenA lenB (HgcdFirstDispatchResult this bound recurse a b W
        scratch lenA lenB heap)) :
    RawExec HgcdRecursiveResult :=
  let m := lenA / 2
  if hbaseGuard : lenB < m + 1 then
    match hbase : hgcdRecursiveBase M computeM A B a b lenA lenB heap with
    | .error fault => .error fault
    | .ok result => .ok (result.toResult
        (hgcdRecursiveBase_result_valid M hM computeM A B a b lenA lenB heap
          result hbase))
  else
    let ws := hgcdRecursiveWorkspace W lenA
    have hR : ws.R.Valid := hgcdRecursiveWorkspace_R_valid W lenA
    have hS : ws.S.Valid := hgcdRecursiveWorkspace_S_valid W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    have hfirstDecrease : high.lenA0 < bound := by
      rw [← hbound]
      exact hgcdRecursiveHighInput_len_lt a b lenA lenB horder (by omega)
    match hfirst : hgcdRecursiveDispatchBelow this bound recurse ws.R hR ws.a3
        ws.b3 high.a0 high.b0 high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1
        scratch ws.a2 ws.next heap
        (hgcdRecursiveHighInput_order a b lenA lenB horder) hfirstDecrease with
    | .error fault => .error fault
    | .ok first =>
      have hfirstLength : HgcdRecursiveLengthInvariant high.lenA0 first :=
        firstLength hbaseGuard first
          (hgcdRecursiveHighInput_order a b lenA lenB horder)
          hfirstDecrease hfirst
      match hreconstruct : hgcdRecursiveReconstructPair this ws.a2 ws.b2 ws.T0
          a b ws.a3 ws.b3 scratch (Nat.min lenA m) (Nat.min lenB m)
          first.lenA first.lenB m first.matrix first.valid first.sgn first.heap with
      | .error fault => .error fault
      | .ok reconstructed =>
        have hreconstructedInvariant :
            HgcdFirstReconstructionInvariant lenA first reconstructed :=
          reconstructionBound hbaseGuard first reconstructed ⟨_, _, hfirst⟩ (by
            simpa [high, hgcdRecursiveHighInput] using hfirstLength)
            hreconstruct
        if hearlyGuard : reconstructed.lenB < m + 1 then
          match hearly : hgcdRecursiveEarlyReturn M first.matrix hM first.valid
              computeM A B ws.a2 ws.b2 reconstructed.lenA reconstructed.lenB
              first.sgn reconstructed.heap with
          | .error fault => .error fault
          | .ok result => .ok (result.toResult
              (hgcdRecursiveEarlyReturn_result_valid M first.matrix hM
                first.valid computeM A B ws.a2 ws.b2 reconstructed.lenA
                reconstructed.lenB first.sgn reconstructed.heap result hearly))
        else
          match hmiddle : hgcdRecursiveMiddle this ws.q ws.d ws.a2 ws.b2
              reconstructed.lenA reconstructed.lenB m ws.W3 reconstructed.heap with
          | .error fault => .error fault
          | .ok middle =>
            have hreconstructedStrict : reconstructed.lenB < lenA :=
              hreconstructedInvariant.decreases
            have hreconstructed : reconstructed.lenB ≤ lenA :=
              Nat.le_of_lt hreconstructedStrict
            have hremainder : middle.lenD < reconstructed.lenB :=
              polyDivrem_remainder_lt this ws.q ws.d ws.a2 reconstructed.lenA
                ws.b2 reconstructed.lenB ws.W3 reconstructed.heap middle.heap
                middle.lenQ middle.lenD
                (hgcdRecursiveMiddle_layout this ws.q ws.d ws.a2 ws.b2
                  reconstructed.lenA reconstructed.lenB m ws.W3
                  reconstructed.heap middle hmiddle).1
            have hsecondBounds := hgcdRecursiveMiddle_second_call_bounds this
              ws.q ws.d ws.a2 ws.b2 reconstructed.lenA reconstructed.lenB m
              lenA ws.W3 reconstructed.heap middle (by omega) (by omega)
              hreconstructed (by omega) hremainder hmiddle
            have hsecondDecrease : middle.lenC0 < bound := by
              rw [← hbound]
              exact hsecondBounds.2.2
            match hgcdRecursiveDispatchBelow this bound recurse ws.S hS ws.a3
                ws.b3 middle.c0 middle.d0 middle.lenC0 middle.lenD0 ws.a2 ws.W3
                ws.T0 ws.T1 scratch ws.a2 ws.next middle.heap
                hsecondBounds.2.1 hsecondDecrease with
            | .error fault => .error fault
            | .ok second =>
              match hgcdRecursiveFinish this M first.matrix second.matrix hM
                  first.valid second.valid computeM A B ws.T0 ws.b2 ws.d ws.a3
                  ws.b3 ws.q (Nat.min reconstructed.lenB middle.k)
                  (Nat.min middle.lenD middle.k) second.lenA second.lenB middle.k
                  middle.lenQ ws.a2 scratch first.sgn second.sgn second.heap with
              | .error fault => .error fault
              | .ok result => .ok result.toResult

/-- The complete source body of `_hgcd_recursive`, parameterized only at its
two genuine recursive call sites.  Every other operation is one of the exact
generated helpers above and occurs in C++ order. -/
def hgcdRecursiveBody (this : DenseUPolyZp) (recurse : HgcdRecursiveCall)
    (M : HgcdMat) (hM : M.Valid) (computeM : Bool)
    (A B a b : RawPtr UInt64) (lenA lenB : Nat)
    (W scratch : RawPtr UInt64) (heap : RawHeap) :
    RawExec HgcdRecursiveResult :=
  let m := lenA / 2
  if lenB < m + 1 then
    match hbase : hgcdRecursiveBase M computeM A B a b lenA lenB heap with
    | .error fault => .error fault
    | .ok result => .ok (result.toResult
        (hgcdRecursiveBase_result_valid M hM computeM A B a b lenA lenB heap
          result hbase))
  else
    let ws := hgcdRecursiveWorkspace W lenA
    have hR : ws.R.Valid := hgcdRecursiveWorkspace_R_valid W lenA
    have hS : ws.S.Valid := hgcdRecursiveWorkspace_S_valid W lenA
    let high := hgcdRecursiveHighInput a b lenA lenB
    match hgcdRecursiveDispatch this recurse ws.R hR ws.a3 ws.b3 high.a0
        high.b0 high.lenA0 high.lenB0 ws.q ws.W3 ws.T0 ws.T1 scratch ws.a2
        ws.next heap with
    | .error fault => .error fault
    | .ok first =>
      let lenLowA := Nat.min lenA m
      let lenLowB := Nat.min lenB m
      match hgcdRecursiveReconstructPair this ws.a2 ws.b2 ws.T0 a b ws.a3
          ws.b3 scratch lenLowA lenLowB first.lenA first.lenB m first.matrix
          first.valid first.sgn first.heap with
      | .error fault => .error fault
      | .ok reconstructed =>
        if reconstructed.lenB < m + 1 then
          match hearly : hgcdRecursiveEarlyReturn M first.matrix hM first.valid
              computeM A B ws.a2 ws.b2 reconstructed.lenA reconstructed.lenB
              first.sgn reconstructed.heap with
          | .error fault => .error fault
          | .ok result => .ok (result.toResult
              (hgcdRecursiveEarlyReturn_result_valid M first.matrix hM
                first.valid computeM A B ws.a2 ws.b2 reconstructed.lenA
                reconstructed.lenB first.sgn reconstructed.heap result hearly))
        else
          match hgcdRecursiveMiddle this ws.q ws.d ws.a2 ws.b2
              reconstructed.lenA reconstructed.lenB m ws.W3 reconstructed.heap with
          | .error fault => .error fault
          | .ok middle =>
            match hgcdRecursiveDispatch this recurse ws.S hS ws.a3 ws.b3
                middle.c0 middle.d0 middle.lenC0 middle.lenD0 ws.a2 ws.W3
                ws.T0 ws.T1 scratch ws.a2 ws.next middle.heap with
            | .error fault => .error fault
            | .ok second =>
              match hgcdRecursiveFinish this M first.matrix second.matrix hM
                  first.valid second.valid computeM A B ws.T0 ws.b2 ws.d ws.a3
                  ws.b3 ws.q (Nat.min reconstructed.lenB middle.k)
                  (Nat.min middle.lenD middle.k) second.lenA second.lenB middle.k
                  middle.lenQ ws.a2 scratch first.sgn second.sgn second.heap with
              | .error fault => .error fault
              | .ok result => .ok result.toResult

end Generated.StrictHGCD
