# Strict generated `__factor_Zp` entry

Date: 2026-08-13

## What changed

- Added a reproducible generator for the strict L1 translation of the C++
  `__factor_Zp` entry.
- Generated the source-shaped make-monic/SQF/DDF/EDF control flow, including
  all three nested range-for loops and an explicit effectful boundary for the
  final `std::sort` call.
- Replaced translated `partial def` loops with structural well-founded
  recursion measured by the number of unvisited array elements.

## Why

The existing `Generated.Corpus.__factor_Zp_ir` is partial and invokes legacy
non-strict component translations.  It therefore cannot be the executable
side of a genuine end-to-end C++ L1-to-L2 refinement theorem.

## Decisions

- Component calls remain explicit `RawExec` operations.  The later refinement
  layer must instantiate them with the already verified strict SQF, DDF and
  EDF executions; the generated L1 contains no L2 factorization oracle.
- `std::sort` is likewise an explicit L1 operation.  Its refinement provider
  must prove that the returned array is a permutation ordered by source degree;
  Lean's executable `Array.qsort` has no public permutation theorem and is not
  silently assumed correct.
- No fuel parameter is introduced.  Each finite C++ range-for loop uses
  `array.size - index` as its termination measure.
- The generator supports `--check`, making the checked-in Lean artifact
  reproducible and detecting stale manual edits.

## Verification

- `python3 cpp2lean_v2/tests/build_strict_factor_zp.py --check`
- `lake env lean CLPoly/Generated/StrictFactorZp.lean`

## Files

- `proof/cpp2lean_v2/tests/build_strict_factor_zp.py`
- `proof/lean/CLPoly/Generated/StrictFactorZp.lean`
