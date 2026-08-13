"""Generate the strict, well-founded L1 control flow for C++ ``__factor_Zp``.

The generated entry preserves the source make-monic/SQF/DDF/EDF call order
and all three nested range-for loops.  Component implementations are supplied
as ``RawExec`` operations, so this artifact cannot dispatch to an L2 model.
"""

from __future__ import annotations

import argparse
from pathlib import Path


OUT = (
    Path(__file__).resolve().parents[2]
    / "lean"
    / "CLPoly"
    / "Generated"
    / "StrictFactorZp.lean"
)


def generate_strict_factor_zp() -> str:
    return r'''-- Auto-generated strict raw C++ control flow for `__factor_Zp`.
import CLPoly.Generated.StrictSquarefreeZp
import CLPoly.Generated.StrictDDF
import CLPoly.Generated.StrictEDF

set_option autoImplicit false

namespace Generated.StrictFactorZp

/-- The concrete callees of the C++ `__factor_Zp` entry.  Every field is an
effectful L1 execution; no L2 polynomial algorithm occurs in this interface. -/
structure FactorZpRawOps (State : Type) where
  makeMonic : SparsePolyZp → RawExec (Zp × SparsePolyZp)
  squarefree : SparsePolyZp → RawExec (Array (SparsePolyZp × UInt64))
  ddf : SparsePolyZp → RawExec (Array (SparsePolyZp × UInt64))
  edf : Array SparsePolyZp → SparsePolyZp → UInt64 → State →
    RawExec (Array SparsePolyZp × State)
  sortByDegree : Array (SparsePolyZp × UInt64) →
    RawExec (Array (SparsePolyZp × UInt64))

/-- Exact lowering of the innermost C++ range-for which attaches the SQF
multiplicity to every EDF factor. -/
def _loop___factor_Zp_0_raw_ir (factors : Array SparsePolyZp)
    (multiplicity : UInt64) (index : Nat)
    (result : Array (SparsePolyZp × UInt64)) :
    Array (SparsePolyZp × UInt64) :=
  if hindex : index < factors.size then
    _loop___factor_Zp_0_raw_ir factors multiplicity (index + 1)
      (result.push (factors[index], multiplicity))
  else
    result
termination_by factors.size - index
decreasing_by omega

/-- Exact lowering of the DDF-component loop. -/
def _loop___factor_Zp_1_raw_ir {State : Type} (ops : FactorZpRawOps State)
    (components : Array (SparsePolyZp × UInt64)) (index : Nat) (rng : State)
    (result : Array (SparsePolyZp × UInt64)) :
    RawExec (Array (SparsePolyZp × UInt64) × State) :=
  if hindex : index < components.size then
    let component := components[index]
    match ops.edf #[] component.1 component.2 rng with
    | .error fault => .error fault
    | .ok (factors, rng') =>
      _loop___factor_Zp_1_raw_ir ops components (index + 1) rng'
        (_loop___factor_Zp_0_raw_ir factors component.2 0 result)
  else
    .ok (result, rng)
termination_by components.size - index
decreasing_by omega

/-- Exact lowering of the SQF-component loop. -/
def _loop___factor_Zp_2_raw_ir {State : Type} (ops : FactorZpRawOps State)
    (squarefreeFactors : Array (SparsePolyZp × UInt64)) (index : Nat)
    (rng : State) (result : Array (SparsePolyZp × UInt64)) :
    RawExec (Array (SparsePolyZp × UInt64)) :=
  if hindex : index < squarefreeFactors.size then
    let factor := squarefreeFactors[index]
    match ops.ddf factor.1 with
    | .error fault => .error fault
    | .ok components =>
      match _loop___factor_Zp_1_raw_ir ops components 0 rng result with
      | .error fault => .error fault
      | .ok (result', rng') =>
        _loop___factor_Zp_2_raw_ir ops squarefreeFactors (index + 1) rng'
          result'
  else
    .ok result
termination_by squarefreeFactors.size - index
decreasing_by omega

/-- Strict generated L1 for the original C++ `__factor_Zp` entry.  The
control flow is source-shaped: constant branch, make-monic, SQF, nested
DDF/EDF loops, multiplicity attachment, then sorting. -/
def __factor_Zp_raw_ir {State : Type} (ops : FactorZpRawOps State)
    (initialRng : State) (f : SparsePolyZp) :
    RawExec (Zp × Array (SparsePolyZp × UInt64)) :=
  if f.isEmpty then
    .error .assertionFailure
  else if get_deg f ≤ 0 then
    .ok (f[0]!.2, #[])
  else
    match ops.makeMonic f with
    | .error fault => .error fault
    | .ok (leadingCoefficient, monic) =>
      match ops.squarefree monic with
      | .error fault => .error fault
      | .ok squarefreeFactors =>
        match _loop___factor_Zp_2_raw_ir ops squarefreeFactors 0 initialRng #[] with
        | .error fault => .error fault
        | .ok result =>
          match ops.sortByDegree result with
          | .error fault => .error fault
          | .ok sorted => .ok (leadingCoefficient, sorted)

end Generated.StrictFactorZp
'''


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    generated = generate_strict_factor_zp()
    if args.check:
        if not OUT.exists() or OUT.read_text() != generated:
            raise SystemExit(f"generated output is stale: {OUT}")
        return
    OUT.write_text(generated)


if __name__ == "__main__":
    main()
