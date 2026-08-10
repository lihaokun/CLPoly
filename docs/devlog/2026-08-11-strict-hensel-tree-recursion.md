# Strict Hensel tree-recursion refinement

## What changed

- Added the strict generated L1 program
  `Generated.StrictHensel.__hensel_lift_recursive_raw_ir`.
- Replaced the legacy partial recursion with structural recursion over a finite
  `HenselLiftTree` topology certificate.  The certificate contains only node
  indices and branch shape; it contains no polynomial result or L2 output.
- Preserved the exact C++ order: lift and store the current node, traverse the
  left child using the updated `g`, re-read the parent after the left traversal,
  then traverse the right child using that concrete `h`.
- Proved a raw-to-safe execution bridge from source safety conditions, with no
  fuel parameter.
- Added `HenselLiftRecursiveCorrect`, a semantic execution trace tying every
  visited generated raw step to `HenselStepCorrect`.
- Proved the complete recursive L1-to-L2 theorem and published it through
  cpp2lean Pass 9 as
  `Refinement.__hensel_lift_recursive_raw_ir_refines` in the centralized
  `Refinement/Generated.lean` module.

## Why

The full `__hensel_lift` entry calls this top-down tree traversal inside its
quadratic precision loop.  Closing the traversal is therefore required before
the outer `m := m * m` loop can be refined.  A fuel-based iterator or a model
fallback would not establish refinement of the actual C++ recursion.

## Key decisions

- Termination evidence is restricted to finite topology.  All polynomial
  values are produced only by executing the strict generated program.
- Descendant preconditions are universally quantified over exact raw outputs;
  they cannot choose a preferred output.
- The public theorem exposes both the exact generated run equation and the
  recursive semantic trace, making skipped or substituted node updates
  impossible.

## Files

- `proof/lean/CLPoly/Generated/StrictHensel.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `proof/lean/CLPoly/Refinement/Generated.lean`
- `proof/cpp2lean_v2/class_map.py`
- `proof/cpp2lean_v2/passes/pass9_refine_gen.py`
- `docs/devlog/2026-08-11-strict-hensel-tree-recursion.md`
