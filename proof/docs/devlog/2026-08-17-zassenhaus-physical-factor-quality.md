# Zassenhaus physical factor quality

## What changed

The strict recombination refinement now proves a single execution-level
quality predicate for every factor physically returned by the generated
Zassenhaus implementation.  Each returned integer polynomial is:

- a nonunit over `Int`;
- still a nonunit after reduction modulo the selected prime; and
- headed by a coefficient that remains nonzero modulo that prime.

The proof covers the literal fixed-size candidate scan, successful removal,
the well-founded outer recursion, the zero/one-lifted-factor fast path,
fallback appending, final-factor appending, and degree sorting.

## Why

The final `FactorZZCorrect` composition must turn equal physical/modular
factor cardinality into pointwise modular irreducibility.  That argument
requires nonunit and nonzero-leading facts about the actual generated output;
assuming them at the recombination boundary would weaken the requested real
C++ L1-to-L2 refinement.

## Key decisions

- `PhysicalFactorQuality` is stated directly over the generated array, so
  sorting and appending can be proved by membership/permutation arguments.
- Candidate quality is derived from the exact returned candidate indices,
  their product of live irreducible Hensel factors, and the actual primitive
  exact-division equation.
- The loop theorem follows `zassenhausLoop` with the same lexicographic
  well-founded measure.  No fuel, partial definition, existence witness, or
  semantic result provider is introduced.

## Verification

- `lake env lean CLPoly/Refinement/Recombine.lean`: passed.
- The new public lemmas were added to `Refinement/AxiomAudit.lean`.
- Production C++ changes: none.  Therefore this proof-only step creates no
  new native or C++/Lean B2B change surface; the previously verified C++
  execution remains unchanged.

## Files

- `proof/lean/CLPoly/Refinement/Recombine.lean`
- `proof/lean/CLPoly/Refinement/AxiomAudit.lean`
- `proof/docs/devlog/2026-08-17-zassenhaus-physical-factor-quality.md`

## Progress and estimate

This closes the Zassenhaus half of the physical output-quality invariant.
The next step is to propagate the same predicate through candidate validation,
van-Hoeij retry/extraction, and its literal Zassenhaus fallback, then consume
it in the equal-cardinality final `FactorZZCorrect` theorem.  Overall goal
progress is approximately 99.68%; estimated remaining time is 3--5 full
working days, dominated by final FactorZZ composition, concrete pipeline
wiring, and complete regression/audit passes.
