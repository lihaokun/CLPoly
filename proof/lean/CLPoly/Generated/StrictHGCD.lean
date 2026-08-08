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

/-- Heap and accumulated offset after the first matrix-stabilization loop. -/
structure HgcdMatStageResult where
  heap : RawHeap
  off : Nat

def hgcdMatPtrRaw (M : HgcdMat) (hM : M.Valid) (i : Fin 4) :
    RawPtr UInt64 := M.poly[i.val]'(by rw [hM.1]; exact i.isLt)

def hgcdMatLenRaw (M : HgcdMat) (hM : M.Valid) (i : Fin 4) : Nat :=
  M.len[i.val]'(by rw [hM.2]; exact i.isLt)

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

end Generated.StrictHGCD
