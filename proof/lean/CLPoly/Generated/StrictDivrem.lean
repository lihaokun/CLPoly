-- AST-validated specialized lowering of dense_upoly_zp::_poly_divrem.
import CLPoly.Generated.StrictGCD

namespace Generated.StrictDivrem

open Generated.StrictGCD

def initW3Loop (heap : RawHeap) (A : RawPtr UInt64) (W3 : RawPtr Word3)
    (lenA i : Nat) : RawExec RawHeap :=
  if h : i < lenA then
    match heap.readU64 A i with
    | .error fault => .error fault
    | .ok value =>
      match heap.writeWord3 W3 i { lo := value, mid := 0, hi := 0 } with
      | .error fault => .error fault
      | .ok heap' => initW3Loop heap' A W3 lenA (i + 1)
  else
    .ok heap
termination_by lenA - i
decreasing_by omega

def addMulLoop (heap : RawHeap) (B : RawPtr UInt64) (W3 : RawPtr Word3)
    (i d j : Nat) (c : UInt64) : RawExec RawHeap :=
  if h : j ≤ d then
    match heap.readU64 B j with
    | .error fault => .error fault
    | .ok bj =>
      let product := dense_upoly_zp__umul128_ir 0 0 c bj
      match heap.readWord3 W3 (i + j) with
      | .error fault => .error fault
      | .ok accum =>
        let accum' := dense_upoly_zp__add_carry3_ir accum product.fst product.snd
        match heap.writeWord3 W3 (i + j) accum' with
        | .error fault => .error fault
        | .ok heap' => addMulLoop heap' B W3 i d (j + 1) c
  else
    .ok heap
termination_by d + 1 - j
decreasing_by omega

def quotientLoop (this : DenseUPolyZp) (Q : RawPtr UInt64)
    (B : RawPtr UInt64) (W3 : RawPtr Word3) (d : Nat) (invLc : UInt64) :
    (heap : RawHeap) → (ii : Nat) → RawExec RawHeap
  | heap, 0 => .ok heap
  | heap, ii + 1 =>
    let i := ii
    match heap.readWord3 W3 (i + d) with
    | .error fault => .error fault
    | .ok accum =>
      let r := dense_upoly_zp__lll_mod_preinv_ir
        accum.hi accum.mid accum.lo this._p this._ninv this._norm
      let qi := dense_upoly_zp_nmod_mul_ir this r invLc
      match heap.writeU64 Q i qi with
      | .error fault => .error fault
      | .ok heap' =>
        if qi != 0 then
          match addMulLoop heap' B W3 i d 0 (this._p - qi) with
          | .error fault => .error fault
          | .ok heap'' => quotientLoop this Q B W3 d invLc heap'' ii
        else
          quotientLoop this Q B W3 d invLc heap' ii

def remainderLoop (this : DenseUPolyZp) (R : RawPtr UInt64)
    (W3 : RawPtr Word3) (d i : Nat) (heap : RawHeap) : RawExec RawHeap :=
  if h : i < d then
    match heap.readWord3 W3 i with
    | .error fault => .error fault
    | .ok accum =>
      let value := dense_upoly_zp__lll_mod_preinv_ir
        accum.hi accum.mid accum.lo this._p this._ninv this._norm
      match heap.writeU64 R i value with
      | .error fault => .error fault
      | .ok heap' => remainderLoop this R W3 d (i + 1) heap'
  else
    .ok heap
termination_by d - i
decreasing_by omega

def dense_upoly_zp__poly_divrem_ir (this : DenseUPolyZp)
    (Q R : RawPtr UInt64) (A : RawPtr UInt64) (lenA : Nat)
    (B : RawPtr UInt64) (lenB : Nat) (W3 : RawPtr Word3)
    (heap : RawHeap) : RawExec (RawHeap × Nat × Nat) :=
  match lenB with
  | 0 => .error .assertionFailure
  | d + 1 =>
    if lenA < d + 1 then
      match heap.copyU64 R A lenA with
      | .error fault => .error fault
      | .ok heap' => .ok (heap', 0, lenA)
    else
      let qLen := lenA - d
      match heap.readU64 B d with
      | .error fault => .error fault
      | .ok lead =>
        let invLc := dense_upoly_zp_nmod_inv_ir this lead
        match initW3Loop heap A W3 lenA 0 with
        | .error fault => .error fault
        | .ok heap1 =>
          match quotientLoop this Q B W3 d invLc heap1 qLen with
          | .error fault => .error fault
          | .ok heap2 =>
            match remainderLoop this R W3 d 0 heap2 with
            | .error fault => .error fault
            | .ok heap3 =>
              match heap3.normaliseU64 Q qLen with
              | .error fault => .error fault
              | .ok lenQ =>
                match heap3.normaliseU64 R d with
                | .error fault => .error fault
                | .ok lenR => .ok (heap3, lenQ, lenR)

end Generated.StrictDivrem
