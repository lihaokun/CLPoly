"""Generate only the currently honest strict C++ L1 dependencies for DDF.

The DDF entry itself must not be emitted until its C++ ``polynomial_GCD``
dependency is translated as L1.  Mapping that call to ``CLPoly.Model`` would
silently turn the purported L1 proof into an L2 fallback.

The same restriction applies to ``__upoly_mod`` and ``__upoly_powmod``: their
current translation reaches ``pair_vec_div5`` and hence the hand-written L2
``SparsePolyZp.divmod`` instance.  They are deliberately excluded until the
raw four-loop division implementation is connected to their source call.
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
  Invariant : UInt64 → SparsePolyZp → SparsePolyZp → UInt64 → Prop
  splitStep : ∀ d fStar h p hPow gdRaw gd quotient hNext,
    Invariant d fStar h p →
    powmod h p.toNat fStar = .ok hPow →
    gcd (__upoly_subtract_x_ir hPow p) fStar = .ok gdRaw →
    (!gdRaw.isEmpty && get_deg gdRaw > 0) = true →
    makeMonic gdRaw = .ok gd →
    exactDiv fStar gd = .ok quotient →
    mod hPow (SparsePolyZp.normalization quotient) = .ok hNext →
    Invariant (d + 1) (SparsePolyZp.normalization quotient) hNext p ∧
      ddfRawMeasure (SparsePolyZp.normalization quotient) (d + 1) <
        ddfRawMeasure fStar d
  noSplitStep : ∀ d fStar h p hPow gdRaw,
    Invariant d fStar h p →
    powmod h p.toNat fStar = .ok hPow →
    gcd (__upoly_subtract_x_ir hPow p) fStar = .ok gdRaw →
    (!gdRaw.isEmpty && get_deg gdRaw > 0) = false →
    Invariant (d + 1) fStar hPow p ∧
      ddfRawMeasure fStar (d + 1) < ddfRawMeasure fStar d

def _loop___ddf_Zp_raw_ir (ops : DDFRawOps) (d : UInt64)
    (fStar h : SparsePolyZp)
    (result : Array (SparsePolyZp × UInt64)) (p : UInt64)
    (hvalid : ops.Invariant d fStar h p) :
    RawExec (SparsePolyZp × Array (SparsePolyZp × UInt64)) :=
  if hterm : get_deg fStar < (2 * d).toInt64 then
    .ok (fStar, result)
  else
    match hpow : ops.powmod h p.toNat fStar with
    | .error fault => .error fault
    | .ok hPow =>
      let hMinusX := __upoly_subtract_x_ir hPow p
      match hgcd : ops.gcd hMinusX fStar with
      | .error fault => .error fault
      | .ok gdRaw =>
        if hsplit : !gdRaw.isEmpty && get_deg gdRaw > 0 then
          match hmonic : ops.makeMonic gdRaw with
          | .error fault => .error fault
          | .ok gd =>
            match hdiv : ops.exactDiv fStar gd with
            | .error fault => .error fault
            | .ok quotient =>
              let fNext := SparsePolyZp.normalization quotient
              match hmod : ops.mod hPow fNext with
              | .error fault => .error fault
              | .ok hNext =>
                have hstep := ops.splitStep d fStar h p hPow gdRaw gd
                  quotient hNext hvalid hpow hgcd hsplit hmonic hdiv hmod
                _loop___ddf_Zp_raw_ir ops (d + 1) fNext hNext
                  (result.push (gd, d)) p hstep.1
        else
          have hsplitFalse :
              (!gdRaw.isEmpty && get_deg gdRaw > 0) = false := by
            cases hb : (!gdRaw.isEmpty && get_deg gdRaw > 0) <;> simp_all
          have hstep := ops.noSplitStep d fStar h p hPow gdRaw hvalid
            hpow hgcd hsplitFalse
          _loop___ddf_Zp_raw_ir ops (d + 1) fStar hPow result p hstep.1
termination_by ddfRawMeasure fStar d
decreasing_by
  all_goals exact hstep.2

def __ddf_Zp_raw_ir (ops : DDFRawOps) (f : SparsePolyZp)
    (hinitial : ¬f.isEmpty →
      ops.Invariant 1 f
        #[(UMonomial.mk (1 : Int64),
          __make_zp_ir (1 : Int64) f[0]!.2.prime)] f[0]!.2.prime) :
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    source = generate_strict_ddf()
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
