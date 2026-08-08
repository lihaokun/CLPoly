-- AST-validated raw lowering of dense_upoly_zp multiplication.
import CLPoly.Generated.StrictGCD

namespace Generated.StrictMul

open Generated.StrictGCD

def classicalDotLoop (heap : RawHeap) (A B : RawPtr UInt64)
    (k stop j : Nat) (acc : Word3) : RawExec Word3 :=
  if h : j ≤ stop then
    match heap.readU64 A j with
    | .error fault => .error fault
    | .ok a =>
      match heap.readU64 B (k - j) with
      | .error fault => .error fault
      | .ok b =>
        let product := dense_upoly_zp__umul128_ir 0 0 a b
        let acc' := dense_upoly_zp__add_carry3_ir acc product.fst product.snd
        classicalDotLoop heap A B k stop (j + 1) acc'
  else
    .ok acc
termination_by stop + 1 - j
decreasing_by omega

def classicalOuterLoop (this : DenseUPolyZp) (C A B : RawPtr UInt64)
    (lenA lenB lenC k : Nat) (heap : RawHeap) : RawExec RawHeap :=
  if h : k < lenC then
    let jMin := if k ≥ lenB then k - lenB + 1 else 0
    let jMax := if k < lenA then k else lenA - 1
    match classicalDotLoop heap A B k jMax jMin { lo := 0, mid := 0, hi := 0 } with
    | .error fault => .error fault
    | .ok acc =>
      let value := dense_upoly_zp__lll_mod_preinv_ir
        acc.hi acc.mid acc.lo this._p this._ninv this._norm
      match heap.writeU64 C k value with
      | .error fault => .error fault
      | .ok heap' => classicalOuterLoop this C A B lenA lenB lenC (k + 1) heap'
  else
    .ok heap
termination_by lenC - k
decreasing_by omega

def dense_upoly_zp__classical_mul_ir (this : DenseUPolyZp)
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) : RawExec RawHeap :=
  if lenA = 0 ∨ lenB = 0 then .error .assertionFailure
  else classicalOuterLoop this C A B lenA lenB (lenA + lenB - 1) 0 heap

end Generated.StrictMul
