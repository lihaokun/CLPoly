"""Generate strict raw control flow for the C++ ``__hensel_step`` entry.

The legacy corpus is intentionally unsuitable for refinement: it is partial
and its Hensel aggregate conversion discards every field.  This generator
keeps the source operation order, all six coefficient loops, and every raw
division call explicit.  Mathematical operation laws are proved separately.
"""

from __future__ import annotations

import argparse
from pathlib import Path


OUT = (
    Path(__file__).resolve().parents[2] / "lean" / "CLPoly" / "Generated" /
    "StrictHensel.lean"
)


def generate_strict_hensel() -> str:
    return r'''-- Auto-generated strict raw C++ control flow for `__hensel_step`.
import CLPoly.Model

set_option autoImplicit false
set_option maxHeartbeats 1600000

namespace Generated.StrictHensel

def certifyRawExec {α : Type} (run : RawExec α) :
    RawExec { output : α // run = .ok output } :=
  match hrun : run with
  | .error fault => .error fault
  | .ok output => .ok ⟨output, by simpa using hrun⟩

/-- Exact lowering of the source coefficient loop implementing
`coefficient := fdiv_q(coefficient,m); coefficient := fdiv_r(coefficient,m)`
followed by the in-place zero compaction. -/
def divideCoeffs (f : SparsePolyZZ) (m : ZZ) : SparsePolyZZ :=
  f.map fun term => (term.fst, ZZ.fdiv_q term.snd term.snd m)

def divideThenReduceCoeffs (f : SparsePolyZZ) (m : ZZ) : SparsePolyZZ :=
  f.filterMap fun term =>
    let quotient := ZZ.fdiv_q term.snd term.snd m
    let coefficient := ZZ.fdiv_r quotient quotient m
    if coefficient != 0 then some (term.fst, coefficient) else none

/-- Exact total lowering shared by the four source range-for loops which
multiply every coefficient by `m`. -/
def scaleCoeffs (f : SparsePolyZZ) (m : ZZ) : SparsePolyZZ :=
  f.map fun term => (term.fst, term.snd * m)

structure HenselStepRawOps where
  mul : SparsePolyZZ → SparsePolyZZ → RawExec SparsePolyZZ
  add : SparsePolyZZ → SparsePolyZZ → RawExec SparsePolyZZ
  sub : SparsePolyZZ → SparsePolyZZ → RawExec SparsePolyZZ
  modCoeff : SparsePolyZZ → ZZ → RawExec SparsePolyZZ
  divmodMod : SparsePolyZZ → SparsePolyZZ → ZZ →
    RawExec (SparsePolyZZ × SparsePolyZZ)

/-- Strict sequential translation of C++ `__hensel_step`.  It neither calls
the L2 Hensel model nor manufactures an output when a raw operation fails. -/
def __hensel_step_raw_ir (ops : HenselStepRawOps) (node : HenselNode)
    (f : SparsePolyZZ) (m : ZZ) : RawExec HenselNode := do
  let m2 := m * m
  let gh ← ops.mul node.g node.h
  let difference ← ops.sub f gh
  let e := divideThenReduceCoeffs difference m
  let se ← ops.mul node.s e
  let qr ← ops.divmodMod se node.h m
  let te ← ops.mul node.t e
  let qg ← ops.mul qr.1 node.g
  let tauRaw ← ops.add te qg
  let tau ← ops.modCoeff tauRaw m
  let gRaw ← ops.add node.g (scaleCoeffs tau m)
  let gNew ← ops.modCoeff gRaw m2
  let hRaw ← ops.add node.h (scaleCoeffs qr.2 m)
  let hNew ← ops.modCoeff hRaw m2
  let factorNode := { node with g := gNew, h := hNew }
  let sg ← ops.mul factorNode.s factorNode.g
  let th ← ops.mul factorNode.t factorNode.h
  let one : SparsePolyZZ := #[(UMonomial.mk 0, 1)]
  let oneMinusSg ← ops.sub one sg
  let bezoutDifference ← ops.sub oneMinusSg th
  let ep := divideThenReduceCoeffs bezoutDifference m
  let sep ← ops.mul factorNode.s ep
  let qrBezout ← ops.divmodMod sep factorNode.h m
  let sRaw ← ops.add factorNode.s (scaleCoeffs qrBezout.2 m)
  let sNew ← ops.modCoeff sRaw m2
  let tep ← ops.mul factorNode.t ep
  let qpg ← ops.mul qrBezout.1 factorNode.g
  let tau2Raw ← ops.add tep qpg
  let tau2 ← ops.modCoeff tau2Raw m
  let tRaw ← ops.add factorNode.t (scaleCoeffs tau2 m)
  let tNew ← ops.modCoeff tRaw m2
  return { factorNode with s := sNew, t := tNew }

end Generated.StrictHensel
'''


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    generated = generate_strict_hensel()
    if args.check:
        if not OUT.exists() or OUT.read_text() != generated:
            raise SystemExit(f"stale generated file: {OUT}")
        return
    OUT.write_text(generated)


if __name__ == "__main__":
    main()
