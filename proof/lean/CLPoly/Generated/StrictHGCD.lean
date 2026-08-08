import CLPoly.Generated.StrictMul
import CLPoly.Generated.StrictPolyAddSub

set_option autoImplicit false

namespace Generated.StrictHGCD

open Generated.StrictMul
open Generated.StrictPolyAddSub

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

end Generated.StrictHGCD
