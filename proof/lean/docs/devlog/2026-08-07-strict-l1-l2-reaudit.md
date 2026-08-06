# Strict C++ L1 → Lean L2 re-audit

Date: 2026-08-07

## Acceptance rule

A module is complete only if a theorem starts from the cpp2lean-generated entry
point (or a total definition proved equation-by-equation equivalent to it) and
derives the L2 value/specification under explicit representation and machine
preconditions.

“Generated” means reproducibly emitted from the checked-in C++ source and the
checked-in cpp2lean passes.  Merely inserting a handwritten definition into
`Generated/Corpus.lean` does not make it L1.  A completion claim therefore also
requires regenerating the corpus in a temporary location and obtaining a
zero diff against the checked-in generated artifacts.

The following do **not** count as refinement:

- `if Correct candidate then candidate else certifiedL2`;
- assuming `Correct candidate` in order to prove that the candidate is adopted;
- replacing a generated function by a separately written safe function without
  a raw-to-safe trace/equation theorem;
- using `Classical.choose` or an L2 algorithm as the implementation branch;
- proving only the final Pipeline contract while ignoring which branch produced
  the value.
- editing `Generated/Corpus.lean` directly without making the same definition a
  deterministic result of the cpp2lean pipeline.

## Current strict status

| Module | Strict status | Missing evidence |
|---|---|---|
| SQF | incomplete | move the handwritten totalization into cpp2lean lowering, regenerate, then re-audit the exact generated entry theorem |
| DDF | incomplete | generate a total/trace form from the DDF MIR and prove its terminating execution refines `ddf` |
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

The final mechanical gate must also establish all of the following:

```text
regenerate C++ AST/MIR/Lean artifacts in a temporary directory
diff regenerated artifacts against the checked-in artifacts (zero diff)
scan the strict implementation path for sorry/admit/oracle/L2 fallback
build the complete strict Pipeline dependency closure
print axioms for every exported L1→L2 theorem
```
