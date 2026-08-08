import CLPoly.Model

set_option autoImplicit false

namespace Generated.StrictHGCD

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

end Generated.StrictHGCD
