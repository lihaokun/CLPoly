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

end Generated.StrictHGCD
