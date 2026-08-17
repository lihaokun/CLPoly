# Generated FactorZZ end-to-end contract

## What changed

The Pass 9 verified-contract registry now emits the public theorem
`Refinement.__factor_squarefree_primitive_ZZ_raw_ir_refines_FactorZZCorrect`
into `CLPoly/Refinement/Generated.lean`.  Its executable side fixes the actual
strict select-prime, modular factorization, Hensel-lift, van-Hoeij, and
Zassenhaus implementations.  `CLPoly/Pipeline/L1.lean` directly re-exports
that contract as `strictFactorSquarefreePrimitiveZZStage`.

The Hensel entry no longer receives factor-count or topology machine bounds
through its semantic readiness structure.  A successful literal raw run
proves both the source `2 <= factors.size` guard and the checked `INT_MAX`
range guard.  The tree topology bound follows from the new well-founded proof
that the number of allocated non-singleton interval nodes is strictly less
than the factor count.

The generated SelectPrime, FactorZp, and EDF wrappers were synchronized with
their strengthened physical/canonical conclusions.  Obsolete `h2p` and `hp2`
arguments were removed from Hensel EEA, tree-build, and lift wrappers.

## Why

The final C++-named theorem previously existed only inside the FactorZZ proof
module and was not discoverable in the centralized generated API.  Machine
control-flow facts also remained mixed into semantic readiness even though
the strict C++ L1 execution itself checks them.  Central generation and
execution-derived guards make the refinement boundary mechanically traceable
and prevent callers from supplying facts that the source did not establish.

## Verification

- `lake build CLPoly.Refinement.Hensel`: passed.
- `lake build CLPoly.Refinement.FactorZZ CLPoly.Refinement.Generated`: passed.
- `lake build CLPoly.Pipeline.L1`: passed.
- `python3 proof/cpp2lean_v2/tests/check_verified_refinements.py`: passed.
- Python syntax compilation for `class_map.py` and `pass9_refine_gen.py`:
  passed with the bytecode cache redirected to `/tmp`.
- Relevant banned-construct search found no actual `sorry`, `axiom`,
  `partial`, `fuel`, or `native_decide`; matches were negative comments.
- No C++ production file changed in this step, so there is no new native or
  C++/Lean B2B change surface.

## Remaining work

The public FactorZZ contract still accepts `hhenselInvariant`, quantified over
the actual selected prime and actual Hensel intermediate executions.  It
cannot choose an output, but the next stage must derive all fields available
from selection/input facts and reduce the parameter to genuinely necessary
arithmetic runtime readiness before final full-build and B2B closure.

## Files

- `proof/cpp2lean_v2/class_map.py`
- `proof/cpp2lean_v2/passes/pass9_refine_gen.py`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `proof/lean/CLPoly/Refinement/FactorZZ.lean`
- `proof/lean/CLPoly/Refinement/Generated.lean`
- `proof/lean/CLPoly/Pipeline/L1.lean`
- `proof/lean/CLPoly/Refinement/AxiomAudit.lean`
