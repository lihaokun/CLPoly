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

/-- Any successful raw normalization returns a prefix no longer than the
declared input.  This is a control-flow fact and needs no heap validity
assumption. -/
theorem normaliseU64_result_le_raw (heap : RawHeap) (ptr : RawPtr UInt64)
    (len result : Nat) (hrun : heap.normaliseU64 ptr len = .ok result) :
    result ≤ len := by
  cases len with
  | zero =>
      simp only [RawHeap.normaliseU64] at hrun
      exact Nat.le_of_eq (Except.ok.inj hrun).symm
  | succ n =>
      simp only [RawHeap.normaliseU64] at hrun
      generalize hread : heap.readU64 ptr n = readResult at hrun
      cases readResult with
      | error fault => simp [hread] at hrun
      | ok value =>
          simp only [hread] at hrun
          split at hrun
          next =>
            exact Nat.le_trans
              (normaliseU64_result_le_raw heap ptr n result hrun) (by omega)
          next =>
            have heq : result = n + 1 := Except.ok.inj hrun.symm
            omega

/-- The remainder length returned by the actual generated `_poly_divrem` is
strictly smaller than the nonempty divisor length.  This supplies the
well-founded measure used by HGCD's source while-loop. -/
theorem polyDivrem_remainder_lt (this : DenseUPolyZp)
    (Q R A : RawPtr UInt64) (lenA : Nat) (B : RawPtr UInt64) (lenB : Nat)
    (W3 : RawPtr Word3) (heap heap' : RawHeap) (lenQ lenR : Nat)
    (hrun : dense_upoly_zp__poly_divrem_ir this Q R A lenA B lenB W3 heap =
      .ok (heap', lenQ, lenR)) :
    lenR < lenB := by
  cases lenB with
  | zero => simp [dense_upoly_zp__poly_divrem_ir] at hrun
  | succ d =>
      by_cases hsmall : lenA < d + 1
      · simp only [dense_upoly_zp__poly_divrem_ir, hsmall, ↓reduceIte] at hrun
        generalize hcopy : heap.copyU64 R A lenA = copyResult at hrun
        cases copyResult with
        | error fault => simp [hcopy] at hrun
        | ok copied =>
            have hrun' : (.ok (copied, 0, lenA) :
                RawExec (RawHeap × Nat × Nat)) = .ok (heap', lenQ, lenR) := by
              simpa [hcopy] using hrun
            have heq := Except.ok.inj hrun'
            have hlenR := congrArg (fun x : RawHeap × Nat × Nat => x.2.2) heq
            have hlenR' : lenA = lenR := by simpa using hlenR
            omega
      · simp only [dense_upoly_zp__poly_divrem_ir, hsmall, ↓reduceIte] at hrun
        generalize hlead : heap.readU64 B d = leadResult at hrun
        cases leadResult with
        | error fault => simp [hlead] at hrun
        | ok lead =>
          simp only [hlead] at hrun
          generalize hinit : initW3Loop heap A W3 lenA 0 = initResult at hrun
          cases initResult with
          | error fault => simp [hinit] at hrun
          | ok heap1 =>
            simp only [hinit] at hrun
            generalize hquot : quotientLoop this Q B W3 d
              (dense_upoly_zp_nmod_inv_ir this lead) heap1 (lenA - d) =
                quotResult at hrun
            cases quotResult with
            | error fault => simp [hquot] at hrun
            | ok heap2 =>
              simp only [hquot] at hrun
              generalize hrem : remainderLoop this R W3 d 0 heap2 =
                remResult at hrun
              cases remResult with
              | error fault => simp [hrem] at hrun
              | ok heap3 =>
                simp only [hrem] at hrun
                generalize hnormQ : heap3.normaliseU64 Q (lenA - d) =
                  normQResult at hrun
                cases normQResult with
                | error fault => simp [hnormQ] at hrun
                | ok observedQ =>
                  simp only [hnormQ] at hrun
                  generalize hnormR : heap3.normaliseU64 R d =
                    normRResult at hrun
                  cases normRResult with
                  | error fault => simp [hnormR] at hrun
                  | ok observedR =>
                    have hrun' : (.ok (heap3, observedQ, observedR) :
                        RawExec (RawHeap × Nat × Nat)) =
                        .ok (heap', lenQ, lenR) := by
                      simpa [hnormR] using hrun
                    have heq := Except.ok.inj hrun'
                    have hlenR := congrArg
                      (fun x : RawHeap × Nat × Nat => x.2.2) heq
                    have hlenR' : observedR = lenR := by simpa using hlenR
                    have hle := normaliseU64_result_le_raw heap3 R d
                      observedR hnormR
                    omega

end Generated.StrictDivrem
