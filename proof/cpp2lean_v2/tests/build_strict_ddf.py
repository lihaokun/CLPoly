"""Generate the strict DDF C++ L1 entry and its proved L1→L2 contract.

Allocation-owning calls remain explicit ``RawExec`` operations and are bound
by the refinement layer to verified raw implementations.  The generated entry
therefore never dispatches to the hand-written L2 polynomial algorithms.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))
sys.path.insert(0, str(V2_ROOT / "tests"))

from build_pass8_corpus import generate_corpus
from pass8_codegen import codegen_corpus


OUT = V2_ROOT.parent / "lean" / "CLPoly" / "Generated" / "StrictDDF.lean"
CONTRACT_OUT = (
    V2_ROOT.parent / "lean" / "CLPoly" / "Generated" /
    "StrictDDFRefinement.lean"
)
CPP_ENTRY = "__ddf_Zp_raw_ir"
STRICT_DDF_ROOTS = {
    "__make_zp",
    "__upoly_subtract_x",
}


def _raw_ddf_adapter() -> str:
    """Generated effectful shell for the source DDF control flow.

    Pass 8 still emits pure calls for allocation-owning C++ helpers.  The
    strict artifact instead exposes those calls as ``RawExec`` operations;
    refinement instantiates them only with the proved raw implementations.
    The two decrease fields are execution invariants, not result or semantic
    oracles.
    """
    return r'''
def ddfRawMeasure (fStar : SparsePolyZp) (d : UInt64) : Nat :=
  if fStar.isEmpty then 0 else (get_deg fStar).toNat + 1 - 2 * d.toNat

structure DDFRawOps where
  powmod : SparsePolyZp → Nat → SparsePolyZp → RawExec SparsePolyZp
  gcd : SparsePolyZp → SparsePolyZp → RawExec SparsePolyZp
  makeMonic : SparsePolyZp → RawExec SparsePolyZp
  exactDiv : SparsePolyZp → SparsePolyZp → RawExec SparsePolyZp
  mod : SparsePolyZp → SparsePolyZp → RawExec SparsePolyZp
  Invariant : UInt64 → SparsePolyZp → SparsePolyZp →
    Array (SparsePolyZp × UInt64) → UInt64 → Prop
  splitStep : ∀ d fStar h result p hPow gdRaw gd quotient hNext,
    Invariant d fStar h result p →
    ¬get_deg fStar < (2 * d).toInt64 →
    powmod h p.toNat fStar = .ok hPow →
    gcd (__upoly_subtract_x_ir hPow p) fStar = .ok gdRaw →
    (!gdRaw.isEmpty && get_deg gdRaw > 0) = true →
    makeMonic gdRaw = .ok gd →
    exactDiv fStar gd = .ok quotient →
    mod hPow (SparsePolyZp.normalization quotient) = .ok hNext →
    Invariant (d + 1) (SparsePolyZp.normalization quotient) hNext
        (result.push (gd, d)) p ∧
      ddfRawMeasure (SparsePolyZp.normalization quotient) (d + 1) <
        ddfRawMeasure fStar d
  noSplitStep : ∀ d fStar h result p hPow gdRaw,
    Invariant d fStar h result p →
    ¬get_deg fStar < (2 * d).toInt64 →
    powmod h p.toNat fStar = .ok hPow →
    gcd (__upoly_subtract_x_ir hPow p) fStar = .ok gdRaw →
    (!gdRaw.isEmpty && get_deg gdRaw > 0) = false →
    Invariant (d + 1) fStar hPow result p ∧
      ddfRawMeasure fStar (d + 1) < ddfRawMeasure fStar d

structure DDFRawState (ops : DDFRawOps) (p : UInt64) where
  d : UInt64
  fStar : SparsePolyZp
  h : SparsePolyZp
  result : Array (SparsePolyZp × UInt64)
  valid : ops.Invariant d fStar h result p

def ddfRawStateMeasure {ops : DDFRawOps} {p : UInt64}
    (state : DDFRawState ops p) : Nat :=
  ddfRawMeasure state.fStar state.d

/-- Turn a raw execution into the same execution with its successful result
certified by the actual run equation.  The proof is erased; the runtime
operation and error behavior are unchanged. -/
def certifyRawExec {α : Type} (run : RawExec α) :
    RawExec { output : α // run = .ok output } :=
  match hrun : run with
  | .error fault => .error fault
  | .ok output => .ok ⟨output, by simpa using hrun⟩

/-- Carry the evaluated C++ boolean together with its equality. -/
def certifyBool (value : Bool) : { result : Bool // value = result } :=
  ⟨value, rfl⟩

/-- Pack the dependent invariant with its runtime state.  This keeps match
equations out of the recursive function's type while preserving the same
well-founded proof and exactly the same C++ control flow. -/
def _loop___ddf_Zp_raw_ir_state (ops : DDFRawOps) (p : UInt64)
    (state : DDFRawState ops p) :
    RawExec (SparsePolyZp × Array (SparsePolyZp × UInt64)) :=
  if hterm : get_deg state.fStar < (2 * state.d).toInt64 then
    .ok (state.fStar, state.result)
  else
    match certifyRawExec (ops.powmod state.h p.toNat state.fStar) with
    | .error fault => .error fault
    | .ok hPowRun =>
      let hPow := hPowRun.val
      have hpow := hPowRun.property
      let hMinusX := __upoly_subtract_x_ir hPow p
      match certifyRawExec (ops.gcd hMinusX state.fStar) with
      | .error fault => .error fault
      | .ok hgcdRun =>
        let gdRaw := hgcdRun.val
        have hgcd := hgcdRun.property
        match certifyBool (!gdRaw.isEmpty && get_deg gdRaw > 0) with
        | ⟨true, hsplit⟩ =>
          match certifyRawExec (ops.makeMonic gdRaw) with
          | .error fault => .error fault
          | .ok hmonicRun =>
            let gd := hmonicRun.val
            have hmonic := hmonicRun.property
            match certifyRawExec (ops.exactDiv state.fStar gd) with
            | .error fault => .error fault
            | .ok hdivRun =>
              let quotient := hdivRun.val
              have hdiv := hdivRun.property
              let fNext := SparsePolyZp.normalization quotient
              match certifyRawExec (ops.mod hPow fNext) with
              | .error fault => .error fault
              | .ok hmodRun =>
                let hNext := hmodRun.val
                have hmod := hmodRun.property
                have hstep := ops.splitStep state.d state.fStar state.h
                  state.result p hPow gdRaw gd quotient hNext state.valid hterm
                  hpow hgcd hsplit hmonic hdiv hmod
                _loop___ddf_Zp_raw_ir_state ops p
                  ⟨state.d + 1, fNext, hNext,
                    state.result.push (gd, state.d), hstep.1⟩
        | ⟨false, hsplitFalse⟩ =>
          have hstep := ops.noSplitStep state.d state.fStar state.h
            state.result p hPow gdRaw state.valid hterm hpow hgcd hsplitFalse
          _loop___ddf_Zp_raw_ir_state ops p
            ⟨state.d + 1, state.fStar, hPow, state.result, hstep.1⟩
termination_by ddfRawStateMeasure state
decreasing_by
  all_goals exact hstep.2

def _loop___ddf_Zp_raw_ir (ops : DDFRawOps) (d : UInt64)
    (fStar h : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64)) (p : UInt64)
    (hvalid : ops.Invariant d fStar h result p) :
    RawExec (SparsePolyZp × Array (SparsePolyZp × UInt64)) :=
  _loop___ddf_Zp_raw_ir_state ops p ⟨d, fStar, h, result, hvalid⟩

def __ddf_Zp_raw_ir (ops : DDFRawOps) (f : SparsePolyZp)
    (hinitial : ¬f.isEmpty →
      ops.Invariant 1 f
        #[(UMonomial.mk (1 : Int64),
          __make_zp_ir (1 : Int64) f[0]!.2.prime)] #[] f[0]!.2.prime) :
    RawExec (Array (SparsePolyZp × UInt64)) :=
  if hf : f.isEmpty then
    .error .assertionFailure
  else
    let p := f[0]!.2.prime
    let h : SparsePolyZp :=
      #[(UMonomial.mk (1 : Int64), __make_zp_ir (1 : Int64) p)]
    match _loop___ddf_Zp_raw_ir ops 1 f h #[] p (hinitial hf) with
    | .error fault => .error fault
    | .ok (fStar, result) =>
      if !fStar.isEmpty && get_deg fStar > 0 then
        match ops.makeMonic fStar with
        | .error fault => .error fault
        | .ok monic =>
          .ok (result.push (monic, get_deg monic |>.toUInt64))
      else
        .ok result
'''


def generate_strict_ddf() -> str:
    _, skipped, roots = generate_corpus()
    if skipped:
        details = ", ".join(f"{name}:{reason}" for name, reason in skipped)
        raise RuntimeError(f"full MIR generation skipped functions: {details}")
    selected = [f for f in roots if f.base_name in STRICT_DDF_ROOTS]
    found = {f.base_name for f in selected}
    missing = STRICT_DDF_ROOTS - found
    if missing:
        raise RuntimeError(f"missing strict DDF roots: {sorted(missing)}")
    source = codegen_corpus(selected, namespace="Generated.StrictDDF")
    marker = "end Generated.StrictDDF"
    if not source.endswith(marker):
        raise RuntimeError("strict DDF namespace footer changed")
    source = source.removesuffix(marker) + _raw_ddf_adapter() + "\n" + marker
    if "sorry" in source or "partial def" in source or "      default" in source:
        raise RuntimeError("strict DDF output contains an opaque placeholder")
    if "polynomial_GCD" in source:
        raise RuntimeError("strict DDF output contains the untranslated C++ GCD boundary")
    forbidden_dispatch = (
        "pair_vec_div", "HasPolyDivmod", "SparsePolyZp.divmod",
        "Array.get!", "Array.set!", "front!", "back!",
    )
    if any(token in source for token in forbidden_dispatch):
        raise RuntimeError("strict DDF output contains an L2 polynomial-division dispatch")
    return source


def generate_strict_ddf_contract() -> str:
    """Emit the public theorem under the exact original C++ entry name."""
    theorem_name = f"{CPP_ENTRY}_refines_ddf"
    return f'''-- Auto-generated strict refinement contract for `{CPP_ENTRY}`.
import CLPoly.Refinement.DDF

set_option autoImplicit false

namespace Refinement.StrictDDF

open CLPoly.Math

/-- The cpp2lean-generated C++ DDF entry terminates and decodes to L2 `ddf`.
The double underscores are retained verbatim from the original C++ symbol. -/
theorem {theorem_name}
    (this : DenseUPolyZp) [Fact (Nat.Prime this._p.toNat)]
    (providers : DDFRawProviders this) (f : SparsePolyZp)
    (hfPrime : f[0]!.2.prime = this._p)
    (hfCanonical : SparsePolyZp.Canonical this._p.toNat f)
    (hfDegree : ∀ term ∈ f.toList, term.1.deg < 2 ^ 62)
    (hfMonic : (SparsePolyZp.toPoly this._p.toNat f).Monic)
    (hfSquarefree : Squarefree (SparsePolyZp.toPoly this._p.toNat f)) :
    ∃ output,
      Generated.StrictDDF.{CPP_ENTRY}
          (strictDDFRawOps this providers
            (SparsePolyZp.toPoly this._p.toNat f)) f
          (fun _ => DDFLoopInvariant.initial this f f[0]!.2.prime hfPrime
            hfCanonical hfDegree hfMonic hfSquarefree) = .ok output ∧
      ddfResultToL2 this._p.toNat output =
        ddf (SparsePolyZp.toPoly this._p.toNat f) := by
  exact strictDDFEntryIR_refines_ddf this providers f hfPrime hfCanonical
    hfDegree hfMonic hfSquarefree

end Refinement.StrictDDF
'''


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--contract-output", type=Path, default=CONTRACT_OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source = generate_strict_ddf()
    contract = generate_strict_ddf_contract()
    if args.check:
        if not args.output.exists() or args.output.read_text() != source:
            print(f"FAIL: {args.output} is not reproducible", file=sys.stderr)
            return 1
        if (not args.contract_output.exists() or
                args.contract_output.read_text() != contract):
            print(
                f"FAIL: {args.contract_output} is not reproducible",
                file=sys.stderr,
            )
            return 1
        print(
            f"PASS: {args.output} and {args.contract_output} are "
            "reproducible and placeholder-free"
        )
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(source)
    args.contract_output.parent.mkdir(parents=True, exist_ok=True)
    args.contract_output.write_text(contract)
    print(f"generated {args.output} ({source.count(chr(10)) + 1} lines)")
    print(
        f"generated {args.contract_output} "
        f"({contract.count(chr(10)) + 1} lines)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
