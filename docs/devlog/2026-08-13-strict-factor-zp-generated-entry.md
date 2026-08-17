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

## Follow-up refinement infrastructure

`CLPoly.Refinement.FactorZp` now proves directly over the generated loops:

- the innermost loop appends exactly the concrete EDF factors with the SQF
  multiplicity attached;
- the DDF-component loop can advance only from an actual successful EDF run;
- the SQF-component loop can advance only from actual successful DDF and
  nested EDF-loop runs.

The latter two theorems are parametric in an accumulator predicate, but that
predicate cannot produce or alter execution results.  This permits the next
step to instantiate the predicate with product and irreducibility facts while
keeping every executable output tied to a strict L1 run equation.

`factor_Zp_correct_of_concrete_components` provides the matching L2 algebraic
composition boundary.  Unlike the older function-parametric theorem, it takes
the concrete SQF list and a dependent list of concrete DDF/EDF expansions;
its conclusion refers to exactly their flattened output.

## End-to-end composition finding

While proving the middle DDF-to-EDF loop, the proof detected that the first
strict generator revision attached the DDF degree `gk_dk.second` to final
factors.  The C++ source attaches the outer SQF multiplicity `sj_ej.second`.
The generated middle loop now carries that outer multiplicity explicitly and
uses it in the innermost append loop.  This is a semantic correction, not a
proof-only rewrite.

The refinement layer now also contains proof-erased DDF/EDF call adapters and
a proved `factorLoop1_refines` theorem.  Every factor in its result originates
from a successful strict EDF execution for an element of the actual strict DDF
output, and the theorem proves both exact decoding and product/irreducibility
properties.

`factorLoop2_refines` closes the generated outer SQF loop.  For every concrete
SQF array element it executes the strict DDF adapter and the proved middle EDF
loop, preserves the original SQF multiplicity, and proves that the exact
flattened output has the same powered product as the unvisited SQF suffix with
irreducible monic factors.

## Files

- `proof/cpp2lean_v2/tests/build_strict_factor_zp.py`
- `proof/lean/CLPoly/Generated/StrictFactorZp.lean`
- `proof/lean/CLPoly/Refinement/FactorZp.lean`
- `proof/lean/CLPoly/Pipeline/FactorZp.lean`
