"""Generate strict source-shaped L1 entries for C++ ZZ factorization."""

from __future__ import annotations

import argparse
from pathlib import Path


OUT = (
    Path(__file__).resolve().parents[2]
    / "lean" / "CLPoly" / "Generated" / "StrictFactorZZ.lean"
)


def generate_strict_factor_zz() -> str:
    return r'''-- Auto-generated strict raw C++ control flow for
-- `__lll_factorize` and `__factor_squarefree_primitive_ZZ`.
import CLPoly.Model

set_option autoImplicit false

namespace Generated.StrictFactorZZ

/-- Concrete effectful callees of the two C++ entries.  These fields expose
only raw executions and cannot manufacture an L2 factorization witness. -/
structure FactorZZRawOps where
  selectPrime : SparsePolyZZ → Bool → RawExec PrimeSelectionResult
  heuristicStartingPrecision : SparsePolyZZ → Int32 → UInt64 →
    RawExec (Int32 × Int32)
  henselLift : SparsePolyZZ → Array SparsePolyZp → UInt64 → Int32 →
    RawExec (Array SparsePolyZZ × ZZ)
  vanHoeijRecombine : SparsePolyZZ → Array SparsePolyZZ → ZZ →
    RawExec (Array SparsePolyZZ)

/-- Exact source-shaped lowering of C++ `__lll_factorize`: heuristic
precision, first Hensel/recombine pass, and the conditional full-precision
second pass. -/
def __lll_factorize_raw_ir (ops : FactorZZRawOps) (f : SparsePolyZZ)
    (factors : Array SparsePolyZp) (p : UInt64) :
    RawExec (Array SparsePolyZZ) :=
  let r : Int32 := factors.size.toUInt32.toInt32
  match ops.heuristicStartingPrecision f r p with
  | .error fault => .error fault
  | .ok (aH, aMig) =>
    match ops.henselLift f factors p aH with
    | .error fault => .error fault
    | .ok (liftedH, mH) =>
      match ops.vanHoeijRecombine f liftedH mH with
      | .error fault => .error fault
      | .ok result =>
        if result.size.toUInt32.toInt32 < r && aH < aMig then
          match ops.henselLift f factors p 0 with
          | .error fault => .error fault
          | .ok (liftedMig, mMig) =>
            ops.vanHoeijRecombine f liftedMig mMig
        else
          .ok result

/-- Exact source-shaped lowering of C++
`__factor_squarefree_primitive_ZZ`: select a prime, return the input on the
irreducible/single-factor branch, otherwise execute `__lll_factorize`. -/
def __factor_squarefree_primitive_ZZ_raw_ir (ops : FactorZZRawOps)
    (useLargePrime : Bool) (f : SparsePolyZZ) : RawExec (Array SparsePolyZZ) :=
  if f.isEmpty || get_deg f < 2 then
    .error .assertionFailure
  else
    match ops.selectPrime f useLargePrime with
    | .error fault => .error fault
    | .ok selection =>
      if selection.irreducible || selection.factors.size ≤ 1 then
        .ok #[f]
      else
        __lll_factorize_raw_ir ops f selection.factors selection.prime

end Generated.StrictFactorZZ
'''


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    generated = generate_strict_factor_zz()
    if args.check:
        if not OUT.exists() or OUT.read_text() != generated:
            raise SystemExit(f"generated output is stale: {OUT}")
        return
    OUT.write_text(generated)


if __name__ == "__main__":
    main()
