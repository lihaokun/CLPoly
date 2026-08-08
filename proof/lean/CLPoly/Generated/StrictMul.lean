-- AST-validated raw lowering of dense_upoly_zp multiplication.
import CLPoly.Generated.StrictGCD
import CLPoly.Generated.StrictPolyAddSub

namespace Generated.StrictMul

open Generated.StrictGCD
open Generated.StrictPolyAddSub

def karAddHalvesLoop (this : DenseUPolyZp)
    (A B t1 t2 : RawPtr UInt64) (m i : Nat) (heap : RawHeap) :
    RawExec RawHeap :=
  if h : i < m then
    match heap.readU64 A i with
    | .error fault => .error fault
    | .ok alo =>
      match heap.readU64 A (m + i) with
      | .error fault => .error fault
      | .ok ahi =>
        match heap.readU64 B i with
        | .error fault => .error fault
        | .ok blo =>
          match heap.readU64 B (m + i) with
          | .error fault => .error fault
          | .ok bhi =>
            match heap.writeU64 t1 i (dense_upoly_zp_nmod_add_ir this alo ahi) with
            | .error fault => .error fault
            | .ok heap1 =>
              match heap1.writeU64 t2 i
                  (dense_upoly_zp_nmod_add_ir this blo bhi) with
              | .error fault => .error fault
              | .ok heap2 => karAddHalvesLoop this A B t1 t2 m (i + 1) heap2
  else
    .ok heap
termination_by m - i
decreasing_by omega

def karSubLoop (this : DenseUPolyZp) (dst sub : RawPtr UInt64)
    (count i : Nat) (heap : RawHeap) : RawExec RawHeap :=
  if h : i < count then
    match heap.readU64 dst i with
    | .error fault => .error fault
    | .ok a =>
      match heap.readU64 sub i with
      | .error fault => .error fault
      | .ok b =>
        match heap.writeU64 dst i (dense_upoly_zp_nmod_sub_ir this a b) with
        | .error fault => .error fault
        | .ok heap1 => karSubLoop this dst sub count (i + 1) heap1
  else
    .ok heap
termination_by count - i
decreasing_by omega

def karAssembleLoop (this : DenseUPolyZp) (C sP1 : RawPtr UInt64)
    (m count i : Nat) (heap : RawHeap) : RawExec RawHeap :=
  if h : i < count then
    match heap.readU64 C (m + i) with
    | .error fault => .error fault
    | .ok base =>
      match heap.readU64 sP1 i with
      | .error fault => .error fault
      | .ok cross =>
        match heap.writeU64 C (m + i)
            (dense_upoly_zp_nmod_add_ir this base cross) with
        | .error fault => .error fault
        | .ok heap1 => karAssembleLoop this C sP1 m count (i + 1) heap1
  else
    .ok heap
termination_by count - i
decreasing_by omega

def karOddTail (A B t1 t2 : RawPtr UInt64) (m h : Nat)
    (heap : RawHeap) : RawExec RawHeap :=
  if h > m then
    match heap.readU64 A (m + m) with
    | .error fault => .error fault
    | .ok aTail =>
      match heap.readU64 B (m + m) with
      | .error fault => .error fault
      | .ok bTail =>
        match heap.writeU64 t1 m aTail with
        | .error fault => .error fault
        | .ok heap1 => heap1.writeU64 t2 m bTail
  else
    .ok heap

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

def dense_upoly_zp__kar_mul_ir (this : DenseUPolyZp)
    (C A B : RawPtr UInt64) (n : Nat) (scratch : RawPtr UInt64)
    (heap : RawHeap) : RawExec RawHeap :=
  if n < 16 then
    dense_upoly_zp__classical_mul_ir this C A n B n heap
  else
    let m := n / 2
    let h := n - m
    let t1 := scratch
    let t2 := t1.add h
    let sP0 := t2.add h
    let sP1 := sP0.add (2 * m - 1)
    let recScratch := sP1.add (2 * h - 1)
    match karAddHalvesLoop this A B t1 t2 m 0 heap with
    | .error fault => .error fault
    | .ok heap1 =>
      match karOddTail A B t1 t2 m h heap1 with
      | .error fault => .error fault
      | .ok heap2 =>
        match dense_upoly_zp__kar_mul_ir this sP0 A B m recScratch heap2 with
        | .error fault => .error fault
        | .ok heap3 =>
          match dense_upoly_zp__kar_mul_ir this sP1 t1 t2 h recScratch heap3 with
          | .error fault => .error fault
          | .ok heap4 =>
            match dense_upoly_zp__kar_mul_ir this (C.add (2 * m))
                (A.add m) (B.add m) h recScratch heap4 with
            | .error fault => .error fault
            | .ok heap5 =>
              match karSubLoop this sP1 sP0 (2 * m - 1) 0 heap5 with
              | .error fault => .error fault
              | .ok heap6 =>
                match karSubLoop this sP1 (C.add (2 * m))
                    (2 * h - 1) 0 heap6 with
                | .error fault => .error fault
                | .ok heap7 =>
                  match heap7.copyU64 C sP0 (2 * m - 1) with
                  | .error fault => .error fault
                  | .ok heap8 =>
                    match heap8.writeU64 C (2 * m - 1) 0 with
                    | .error fault => .error fault
                    | .ok heap9 =>
                      karAssembleLoop this C sP1 m (2 * h - 1) 0 heap9
termination_by n
decreasing_by
  all_goals
    have hn : 0 < n := by omega
    have hm : 0 < n / 2 := Nat.div_pos (by omega) (by omega)
    omega

end Generated.StrictMul
