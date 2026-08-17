# Strict C++ L1 → Lean L2 factorization refinement

This directory contains the kernel-checked refinement layer for CLPoly's two
public univariate factorization entries.  The shortest route to the final
result is:

| C++ entry | Generated strict L1 | Public final theorem | L2 conclusion |
|---|---|---|---|
| `__factor_Zp` | [`Generated/StrictFactorZp.lean`](../Generated/StrictFactorZp.lean) | [`__factor_Zp_raw_ir_refines_FactorZpCorrect`](Generated.lean) | `FactorZpCorrect` |
| `__factor_squarefree_primitive_ZZ` | [`Generated/StrictFactorZZ.lean`](../Generated/StrictFactorZZ.lean) | [`__factor_squarefree_primitive_ZZ_raw_ir_refines_FactorZZCorrect`](Generated.lean) | `FactorZZCorrect` |

Both public theorems are generated into [`Generated.lean`](Generated.lean).
[`Pipeline/L1.lean`](../Pipeline/L1.lean) re-exports the same strict programs;
it does not reconstruct the result from an L2 existence theorem.

## Executed stages

The `Zp` entry executes make-monic, squarefree factorization, DDF, EDF,
multiplicity attachment, RNG threading, and the final degree sort.

The `ZZ` entry executes prime enumeration and candidate validation, the full
`Zp` factorization path, multifactor Hensel lifting, and van-Hoeij/Zassenhaus
recombination.  Recursive source loops are represented with structural or
well-founded recursion.  No fuel counter, `partial def`, semantic factor
oracle, or L2 fallback is part of the strict proof boundary.

The main implementation files are organized by responsibility:

- [`SquarefreeZp.lean`](SquarefreeZp.lean): verified raw SQF primitives;
- [`SquarefreeZpEntry.lean`](SquarefreeZpEntry.lean): the generated SQF entry simulation;
- [`DDF.lean`](DDF.lean) and [`EDF.lean`](EDF.lean): modular factor stages;
- [`SelectPrime.lean`](SelectPrime.lean): concrete candidate pipeline;
- [`Hensel.lean`](Hensel.lean): tree construction and quadratic lifting;
- [`Recombine.lean`](Recombine.lean): Zassenhaus and van-Hoeij execution;
- [`FactorZp.lean`](FactorZp.lean) and [`FactorZZ.lean`](FactorZZ.lean): final composition proofs.

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
