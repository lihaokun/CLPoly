-- AST-validated raw lowering of dense_upoly_zp::_poly_add and _poly_sub.
import CLPoly.Model

namespace Generated.StrictPolyAddSub

def dense_upoly_zp_nmod_add_ir (this : DenseUPolyZp)
    (a b : UInt64) : UInt64 :=
  let neg := this._p - a
  if neg > b then a + b else b - neg

def dense_upoly_zp_nmod_sub_ir (this : DenseUPolyZp)
    (a b : UInt64) : UInt64 :=
  if a ≥ b then a - b else this._p - b + a

def addCommonLoop (this : DenseUPolyZp) (C A B : RawPtr UInt64)
    (limit i : Nat) (heap : RawHeap) : RawExec RawHeap :=
  if h : i < limit then
    match heap.readU64 A i with
    | .error fault => .error fault
    | .ok a =>
      match heap.readU64 B i with
      | .error fault => .error fault
      | .ok b =>
        match heap.writeU64 C i (dense_upoly_zp_nmod_add_ir this a b) with
        | .error fault => .error fault
        | .ok heap' => addCommonLoop this C A B limit (i + 1) heap'
  else
    .ok heap
termination_by limit - i
decreasing_by omega

def subCommonLoop (this : DenseUPolyZp) (C A B : RawPtr UInt64)
    (limit i : Nat) (heap : RawHeap) : RawExec RawHeap :=
  if h : i < limit then
    match heap.readU64 A i with
    | .error fault => .error fault
    | .ok a =>
      match heap.readU64 B i with
      | .error fault => .error fault
      | .ok b =>
        match heap.writeU64 C i (dense_upoly_zp_nmod_sub_ir this a b) with
        | .error fault => .error fault
        | .ok heap' => subCommonLoop this C A B limit (i + 1) heap'
  else
    .ok heap
termination_by limit - i
decreasing_by omega

def subNegTailLoop (this : DenseUPolyZp) (C B : RawPtr UInt64)
    (limit i : Nat) (heap : RawHeap) : RawExec RawHeap :=
  if h : i < limit then
    match heap.readU64 B i with
    | .error fault => .error fault
    | .ok b =>
      let value := if b = 0 then 0 else this._p - b
      match heap.writeU64 C i value with
      | .error fault => .error fault
      | .ok heap' => subNegTailLoop this C B limit (i + 1) heap'
  else
    .ok heap
termination_by limit - i
decreasing_by omega

def dense_upoly_zp__poly_add_ir (this : DenseUPolyZp)
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) : RawExec (RawHeap × Nat) :=
  let minLen := min lenA lenB
  let maxLen := max lenA lenB
  match addCommonLoop this C A B minLen 0 heap with
  | .error fault => .error fault
  | .ok heap1 =>
    let tail : RawExec RawHeap :=
      if lenA > lenB then
        if C.sameAddress A then .ok heap1
        else heap1.copyU64 (C.add minLen) (A.add minLen) (lenA - minLen)
      else if lenB > lenA then
        if C.sameAddress B then .ok heap1
        else heap1.copyU64 (C.add minLen) (B.add minLen) (lenB - minLen)
      else .ok heap1
    match tail with
    | .error fault => .error fault
    | .ok heap2 =>
      match heap2.normaliseU64 C maxLen with
      | .error fault => .error fault
      | .ok length => .ok (heap2, length)

def dense_upoly_zp__poly_sub_ir (this : DenseUPolyZp)
    (C A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (heap : RawHeap) : RawExec (RawHeap × Nat) :=
  let minLen := min lenA lenB
  let maxLen := max lenA lenB
  match subCommonLoop this C A B minLen 0 heap with
  | .error fault => .error fault
  | .ok heap1 =>
    let tail : RawExec RawHeap :=
      if lenA > lenB then
        if C.sameAddress A then .ok heap1
        else heap1.copyU64 (C.add minLen) (A.add minLen) (lenA - minLen)
      else if lenB > lenA then
        subNegTailLoop this C B lenB minLen heap1
      else .ok heap1
    match tail with
    | .error fault => .error fault
    | .ok heap2 =>
      match heap2.normaliseU64 C maxLen with
      | .error fault => .error fault
      | .ok length => .ok (heap2, length)

end Generated.StrictPolyAddSub
