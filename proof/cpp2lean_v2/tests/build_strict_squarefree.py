"""Generate the strict raw control flow for the C++ ``__squarefree_Zp``.

The generated shell preserves the source branches and recursion while making
allocation-owning algebraic calls explicit ``RawExec`` operations.  Its two
dependent state packages carry only representation/termination invariants;
they cannot select an algebraic result or invoke the L2 SQF model.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


OUT = (
    Path(__file__).resolve().parents[2] / "lean" / "CLPoly" / "Generated" /
    "StrictSquarefreeZp.lean"
)


def generate_strict_squarefree() -> str:
    return r'''-- Auto-generated strict raw C++ control flow for `__squarefree_Zp`.
import CLPoly.Model

set_option autoImplicit false
set_option maxErrors 2000

namespace Generated.StrictSquarefreeZp

def squarefreeMeasure (f : SparsePolyZp) : Nat :=
  if f.isEmpty then 0 else f[0]!.fst.deg + 1

def scaleMultiplicityLoop (index : Nat)
    (source output : Array (SparsePolyZp × UInt64)) (prime : UInt64) :
    Array (SparsePolyZp × UInt64) :=
  if hmore : index < source.size then
    let item := source[index]
    scaleMultiplicityLoop (index + 1) source
      (output.push (item.1, item.2 * prime)) prime
  else
    output
termination_by source.size - index
decreasing_by omega

def certifyRawExec {α : Type} (run : RawExec α) :
    RawExec { output : α // run = .ok output } :=
  match hrun : run with
  | .error fault => .error fault
  | .ok output => .ok ⟨output, by simpa using hrun⟩

def certifyBool (value : Bool) : { result : Bool // value = result } :=
  ⟨value, rfl⟩

structure SQFRawOps where
  derivative : SparsePolyZp → SparsePolyZp
  extractPthRoot : SparsePolyZp → RawExec SparsePolyZp
  makeMonic : SparsePolyZp → RawExec SparsePolyZp
  gcd : SparsePolyZp → SparsePolyZp → RawExec SparsePolyZp
  exactDiv : SparsePolyZp → SparsePolyZp → RawExec SparsePolyZp
  EntryInvariant : SparsePolyZp → Prop
  YunInvariant : SparsePolyZp → UInt64 → SparsePolyZp → SparsePolyZp →
    Array (SparsePolyZp × UInt64) → Prop
  derivativeZeroStep : ∀ (f root rootMonic : SparsePolyZp),
    EntryInvariant f →
    (derivative f).isEmpty = true →
    extractPthRoot f = .ok root →
    makeMonic root = .ok rootMonic →
    EntryInvariant rootMonic ∧
      squarefreeMeasure rootMonic < squarefreeMeasure f
  yunInitial : ∀ (f derivativeOut c wRaw : SparsePolyZp),
    EntryInvariant f →
    (derivative f).isEmpty = false →
    derivativeOut = derivative f →
    gcd f derivativeOut = .ok c →
    exactDiv f c = .ok wRaw →
    YunInvariant f 1 (SparsePolyZp.normalization wRaw) c #[]
  yunFactorStep : ∀ (source : SparsePolyZp) (multiplicity : UInt64) (w c : SparsePolyZp)
      (result : Array (SparsePolyZp × UInt64))
      (y zRaw zMonic cRaw : SparsePolyZp),
    YunInvariant source multiplicity w c result →
    (!w.isEmpty && get_deg w > 0) = true →
    gcd w c = .ok y →
    exactDiv w y = .ok zRaw →
    (!zRaw.normalization.isEmpty &&
      get_deg zRaw.normalization > 0) = true →
    makeMonic zRaw.normalization = .ok zMonic →
    exactDiv c y = .ok cRaw →
    YunInvariant source (multiplicity + 1) y cRaw.normalization
      (result.push (zMonic, multiplicity)) ∧
    squarefreeMeasure y + squarefreeMeasure cRaw.normalization <
      squarefreeMeasure w + squarefreeMeasure c
  yunNoFactorStep : ∀ (source : SparsePolyZp) (multiplicity : UInt64) (w c : SparsePolyZp)
      (result : Array (SparsePolyZp × UInt64))
      (y zRaw cRaw : SparsePolyZp),
    YunInvariant source multiplicity w c result →
    (!w.isEmpty && get_deg w > 0) = true →
    gcd w c = .ok y →
    exactDiv w y = .ok zRaw →
    (!zRaw.normalization.isEmpty &&
      get_deg zRaw.normalization > 0) = false →
    exactDiv c y = .ok cRaw →
    YunInvariant source (multiplicity + 1) y cRaw.normalization result ∧
    squarefreeMeasure y + squarefreeMeasure cRaw.normalization <
      squarefreeMeasure w + squarefreeMeasure c
  remainderRootStep : ∀ (f derivativeOut : SparsePolyZp)
      (multiplicity : UInt64) (w c : SparsePolyZp)
      (result : Array (SparsePolyZp × UInt64))
      (root rootMonic : SparsePolyZp),
    EntryInvariant f →
    (derivative f).isEmpty = false →
    YunInvariant f multiplicity w c result →
    (!w.isEmpty && get_deg w > 0) = false →
    (!c.isEmpty && get_deg c > 0) = true →
    extractPthRoot c = .ok root →
    makeMonic root = .ok rootMonic →
    EntryInvariant rootMonic ∧
      squarefreeMeasure rootMonic < squarefreeMeasure f

structure YunRawState (ops : SQFRawOps) (source : SparsePolyZp) where
  multiplicity : UInt64
  w : SparsePolyZp
  c : SparsePolyZp
  result : Array (SparsePolyZp × UInt64)
  valid : ops.YunInvariant source multiplicity w c result

def yunStateMeasure {ops : SQFRawOps} {source : SparsePolyZp}
    (state : YunRawState ops source) : Nat :=
  squarefreeMeasure state.w + squarefreeMeasure state.c

structure YunRawFinalState (ops : SQFRawOps) (source : SparsePolyZp) where
  state : YunRawState ops source
  done : (!state.w.isEmpty && get_deg state.w > 0) = false

def _loop___squarefree_Zp_1_raw_ir_state (ops : SQFRawOps)
    {source : SparsePolyZp} (state : YunRawState ops source) :
    RawExec (YunRawFinalState ops source) :=
  match certifyBool (!state.w.isEmpty && get_deg state.w > 0) with
  | ⟨true, hguard⟩ =>
    match certifyRawExec (ops.gcd state.w state.c) with
    | .error fault => .error fault
    | .ok hyRun =>
      let y := hyRun.val
      have hy := hyRun.property
      match certifyRawExec (ops.exactDiv state.w y) with
      | .error fault => .error fault
      | .ok hzRun =>
        let zRaw := hzRun.val
        have hz := hzRun.property
        let z := zRaw.normalization
        match certifyBool (!z.isEmpty && get_deg z > 0) with
        | ⟨true, hzGuard⟩ =>
          match certifyRawExec (ops.makeMonic z) with
          | .error fault => .error fault
          | .ok hmonicRun =>
            let zMonic := hmonicRun.val
            have hmonic := hmonicRun.property
            match certifyRawExec (ops.exactDiv state.c y) with
            | .error fault => .error fault
            | .ok hcRun =>
              let cRaw := hcRun.val
              have hc := hcRun.property
              have hstep := ops.yunFactorStep source state.multiplicity state.w
                state.c state.result y zRaw zMonic cRaw state.valid hguard
                hy hz hzGuard hmonic hc
              _loop___squarefree_Zp_1_raw_ir_state ops
                ⟨state.multiplicity + 1, y, cRaw.normalization,
                  state.result.push (zMonic, state.multiplicity), hstep.1⟩
        | ⟨false, hzGuardFalse⟩ =>
          match certifyRawExec (ops.exactDiv state.c y) with
          | .error fault => .error fault
          | .ok hcRun =>
            let cRaw := hcRun.val
            have hc := hcRun.property
            have hstep := ops.yunNoFactorStep source state.multiplicity state.w
              state.c state.result y zRaw cRaw state.valid hguard hy hz
              hzGuardFalse hc
            _loop___squarefree_Zp_1_raw_ir_state ops
              ⟨state.multiplicity + 1, y, cRaw.normalization,
                state.result, hstep.1⟩
  | ⟨false, hguardFalse⟩ =>
    .ok ⟨state, hguardFalse⟩
termination_by yunStateMeasure state
decreasing_by
  all_goals exact hstep.2

def _loop___squarefree_Zp_1_raw_ir (ops : SQFRawOps)
    (source : SparsePolyZp) (multiplicity : UInt64) (w c : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64))
    (hvalid : ops.YunInvariant source multiplicity w c result) :
    RawExec (SparsePolyZp × Array (SparsePolyZp × UInt64)) :=
  match _loop___squarefree_Zp_1_raw_ir_state ops
      ⟨multiplicity, w, c, result, hvalid⟩ with
  | .error fault => .error fault
  | .ok finalState => .ok (finalState.state.c, finalState.state.result)

structure SQFRawState (ops : SQFRawOps) where
  f : SparsePolyZp
  valid : ops.EntryInvariant f

def sqfStateMeasure {ops : SQFRawOps} (state : SQFRawState ops) : Nat :=
  squarefreeMeasure state.f

def __squarefree_Zp_raw_ir_state (ops : SQFRawOps)
    (state : SQFRawState ops) :
    RawExec (Array (SparsePolyZp × UInt64)) :=
  let f := state.f
  let prime := f[0]!.2.prime
  let derivativeOut := ops.derivative f
  match certifyBool derivativeOut.isEmpty with
  | ⟨true, hderivative⟩ =>
    match certifyRawExec (ops.extractPthRoot f) with
    | .error fault => .error fault
    | .ok hrootRun =>
      let root := hrootRun.val
      have hroot := hrootRun.property
      match certifyRawExec (ops.makeMonic root) with
      | .error fault => .error fault
      | .ok hmonicRun =>
        let rootMonic := hmonicRun.val
        have hmonic := hmonicRun.property
        have hstep := ops.derivativeZeroStep f root rootMonic state.valid
          hderivative hroot hmonic
        match __squarefree_Zp_raw_ir_state ops ⟨rootMonic, hstep.1⟩ with
        | .error fault => .error fault
        | .ok factors => .ok (scaleMultiplicityLoop 0 factors #[] prime)
  | ⟨false, hderivativeFalse⟩ =>
    match certifyRawExec (ops.gcd f derivativeOut) with
    | .error fault => .error fault
    | .ok hcRun =>
      let c := hcRun.val
      have hc := hcRun.property
      match certifyRawExec (ops.exactDiv f c) with
      | .error fault => .error fault
      | .ok hwRun =>
        let wRaw := hwRun.val
        have hw := hwRun.property
        let w := wRaw.normalization
        have hyunInitial := ops.yunInitial f derivativeOut c wRaw state.valid
          hderivativeFalse rfl hc hw
        match _loop___squarefree_Zp_1_raw_ir_state ops
            ⟨1, w, c, #[], hyunInitial⟩ with
        | .error fault => .error fault
        | .ok finalState =>
          match certifyBool (!finalState.state.c.isEmpty &&
              get_deg finalState.state.c > 0) with
          | ⟨true, hguard⟩ =>
            match certifyRawExec (ops.extractPthRoot finalState.state.c) with
            | .error fault => .error fault
            | .ok hrootRun =>
              let root := hrootRun.val
              have hroot := hrootRun.property
              match certifyRawExec (ops.makeMonic root) with
              | .error fault => .error fault
              | .ok hmonicRun =>
                let rootMonic := hmonicRun.val
                have hmonic := hmonicRun.property
                have hstep := ops.remainderRootStep f derivativeOut
                  finalState.state.multiplicity finalState.state.w
                  finalState.state.c finalState.state.result root rootMonic
                  state.valid hderivativeFalse finalState.state.valid
                  finalState.done hguard hroot hmonic
                match __squarefree_Zp_raw_ir_state ops
                    ⟨rootMonic, hstep.1⟩ with
                | .error fault => .error fault
                | .ok factors =>
                  .ok (scaleMultiplicityLoop 0 factors
                    finalState.state.result prime)
          | ⟨false, _⟩ =>
            .ok finalState.state.result
termination_by sqfStateMeasure state
decreasing_by
  all_goals exact hstep.2

def __squarefree_Zp_raw_ir (ops : SQFRawOps) (f : SparsePolyZp)
    (hinitial : ¬f.isEmpty → ops.EntryInvariant f) :
    RawExec (Array (SparsePolyZp × UInt64)) :=
  if hf : f.isEmpty then
    .error .assertionFailure
  else
    __squarefree_Zp_raw_ir_state ops ⟨f, hinitial hf⟩

end Generated.StrictSquarefreeZp
'''


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source = generate_strict_squarefree()
    forbidden = ("sorry", "admit", "axiom", "partial def", "sqfZp")
    if any(token in source for token in forbidden):
        raise RuntimeError("strict SQF output contains a placeholder or L2 fallback")
    if args.check:
        if not args.output.exists() or args.output.read_text() != source:
            print(f"FAIL: {args.output} is not reproducible", file=sys.stderr)
            return 1
        print(f"PASS: {args.output} is reproducible and placeholder-free")
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(source)
    print(f"generated {args.output} ({source.count(chr(10)) + 1} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
