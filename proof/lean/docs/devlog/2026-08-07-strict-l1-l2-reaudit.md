# Strict C++ L1 → Lean L2 re-audit

Date: 2026-08-07

## Acceptance rule

A module is complete only if a theorem starts from the cpp2lean-generated entry
point (or a total definition proved equation-by-equation equivalent to it) and
derives the L2 value/specification under explicit representation and machine
preconditions.

The following do **not** count as refinement:

- `if Correct candidate then candidate else certifiedL2`;
- assuming `Correct candidate` in order to prove that the candidate is adopted;
- replacing a generated function by a separately written safe function without
  a raw-to-safe trace/equation theorem;
- using `Classical.choose` or an L2 algorithm as the implementation branch;
- proving only the final Pipeline contract while ignoring which branch produced
  the value.

## Current strict status

| Module | Strict status | Missing evidence |
|---|---|---|
| SQF | substantially closed | final re-audit of the exact generated entry theorem |
| DDF | incomplete | `Generated.__ddf_Zp_ir` / raw loop → `ddfZpSafe` equation/trace |
| EDF | incomplete | remove `if EDFCorrect`; prove random/trace/powmod split and recursive output invariants |
| Hensel | incomplete | remove `if HenselCorrect`; prove tree construction, step, recursion, extraction, and `p^k` projection |

## Mechanical rejection scan

Before completion, the strict scope must contain no implementation-path matches
for:

```text
if EDFCorrect
if HenselCorrect
edfCertified
else hensel_lift
```

Any fallback must itself be a proved C++/L1 execution branch, not an L2
existence result.
