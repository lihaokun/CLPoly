# Strict C++ L1 → Lean L2 factorization refinement

This directory contains the kernel-checked refinement layer for CLPoly's two
public univariate factorization entries.  The shortest route to the final
result is:

| C++ entry | Generated strict L1 | Public final theorem | L2 conclusion |
|---|---|---|---|
| `__factor_Zp` | [`Generated/StrictFactorZp.lean`](../Generated/StrictFactorZp.lean) | [`__factor_Zp_raw_ir_refines_FactorZpCorrect`](Generated.lean) | `FactorZpCorrect` |
| `__factor_squarefree_primitive_ZZ` | [`Generated/StrictVanHoeijFactorZZ.lean`](../Generated/StrictVanHoeijFactorZZ.lean) | [`__factor_squarefree_primitive_ZZ_raw_ir_refines_FactorZZCorrect`](Generated.lean) | `FactorZZCorrect` |

Both public theorems are generated into [`Generated.lean`](Generated.lean).
[`Pipeline/L1.lean`](../Pipeline/L1.lean) re-exports the same strict programs;
it does not reconstruct the result from an L2 existence theorem.

## Executed stages

The `Zp` entry executes make-monic, squarefree factorization, DDF, EDF,
multiplicity attachment, RNG threading, and the final degree sort.

The `ZZ` entry executes prime enumeration and candidate validation, the full
`Zp` factorization path, multifactor Hensel lifting, and certificate-only
van Hoeij retry with a checked heuristic precision prefix. Recursive source loops are represented with structural or
well-founded recursion.  No fuel counter, `partial def`, semantic factor
oracle, or L2 fallback is part of the strict proof boundary.

The main implementation files are organized by responsibility:

- [`SquarefreeZp.lean`](SquarefreeZp.lean): verified raw SQF primitives;
- [`SquarefreeZpEntry.lean`](SquarefreeZpEntry.lean): the generated SQF entry simulation;
- [`DDF.lean`](DDF.lean) and [`EDF.lean`](EDF.lean): modular factor stages;
- [`SelectPrime.lean`](SelectPrime.lean): concrete candidate pipeline;
- [`Hensel.lean`](Hensel.lean): tree construction and quadratic lifting;
- [`Recombine.lean`](Recombine.lean): shared arithmetic, LLL, and legacy recombination contracts;
- [`VanHoeijEntryTotal.lean`](VanHoeijEntryTotal.lean) and [`VanHoeijRetryTotal.lean`](VanHoeijRetryTotal.lean): actual certified entry and retry success;
- [`FactorZp.lean`](FactorZp.lean) and [`VanHoeijFactorZZ.lean`](VanHoeijFactorZZ.lean): current final composition proofs.

## Trust and acceptance gates

The formal boundary is checked by:

```text
python3 proof/cpp2lean_v2/tests/check_verified_refinements.py
python3 proof/cpp2lean_v2/tests/check_strict_refinement_boundary.py
python3 proof/cpp2lean_v2/tests/check_b2b_strict_runtime_isolation.py
cd proof/lean && lake build
cd proof/lean && lake env lean CLPoly/Refinement/FinalAxiomAudit.lean
python3 proof/b2b/runner/runner.py proof/b2b/vectors
```

The axiom report for both final theorems contains only Lean/Mathlib's standard
`propext`, `Classical.choice`, and `Quot.sound`.  Executable B2B support lives
under `proof/lean/B2B`; an isolation check prevents its test-only erased proof
fields from entering `Generated`, `Refinement`, or `Pipeline`.

## C++ implementation issues exposed by refinement

The formalization was not merely a proof of an unchanged model.  Following
source-level obligations through the strict L1 semantics exposed the following
production C++ defects, all of which were corrected in C++, mirrored in the
generated Lean program, and covered by native or C++/Lean B2B regression tests:

- **Large-prime signed narrowing in DDF/EDF.**  Constructing `-1 mod p` as
  `(int64_t)(p - 1)` overflowed for primes above `2^63` and could make EDF fail
  to terminate.  The implementation now constructs the residue through the
  unsigned `Zp(p - 1, p)` path; large-prime factorization cases exercise the
  corrected branch.
- **Insufficient Karatsuba scratch storage.**  The unequal-length dense
  multiplication path allocated `6*n` words for a recursive routine whose
  verified workspace bound is `7*n`.  The production allocation and its HGCD
  caller accounting were increased to the proved bound.
- **Van-Hoeij lattice and completion certificate (2026-10-16).** The new
  controller uses modulus rows, safely retains the short-vector space, and
  accepts a complete validated partition regardless of its factor count.
  [Retry correctness](VanHoeijRetryCorrectness.lean) connects actual Hensel,
  lattice, and validation executions to product recovery and irreducibility,
  with source and recursive execution-readiness conditions explicit. The
  previous cardinality safety helper and its B2B fixture have been retired.
  Legacy `StrictFactorZZ` controller theorems remain historical dependencies
  of shared mathematical helpers; they are not evidence for the new path.
  [Actual retry totality](VanHoeijRetryTotal.lean) now constructs the checked
  prefix and pure retry success under explicit input precision bounds and the
  existing Hensel readiness contract. [Index width](VanHoeijIndexWidth.lean)
  connects the actual matrix dimension to the widened C++ container indices.
  Final B2B passed 3107/3107, native regression passed all 27 files and 406
  profiling assertions, and the proof audits passed. Three-version performance
  measurements retain baseline timeouts, small-input regressions, and host-load
  limitations; see the [acceptance record](../../../../docs/test-reports/vanhoeij-final-acceptance.md)
  and [source correspondence record](../../../../docs/test-reports/vanhoeij-source-correspondence.md).
- **Unchecked Hensel factor-count narrowing.**  A `size_t` factor count was
  converted to the signed indices used by the Hensel tree without first proving
  it was at most `INT_MAX`.  Both production entry points now reject an
  out-of-range count before conversion.  Boundary tests cover `INT_MAX`,
  `INT_MAX + 1`, and `UINT64_MAX`, and the Lean refinement derives the bound
  from successful execution rather than assuming it at prime selection.

The detailed records are in
[`formal-proof-ddf-edf.md`](../../../../docs/design/factorization/formal-proof-ddf-edf.md),
[`2026-08-08-strict-mul-lowering.md`](../../../../docs/devlog/2026-08-08-strict-mul-lowering.md),
[`2026-08-17-vanhoeij-coordinate-and-safety-fix-plan.md`](../../../../docs/devlog/2026-08-17-vanhoeij-coordinate-and-safety-fix-plan.md),
and
[`2026-08-17-hensel-factor-count-guard.md`](../../../docs/devlog/2026-08-17-hensel-factor-count-guard.md).
