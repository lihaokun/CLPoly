# Hensel algebraic and runtime readiness split

## What changed

The complete strict Hensel entry invariant is now factored into two explicit
structures:

- `HenselLiftAlgebraicInvariant` contains source, target, and adjusted-factor
  facts that the outer FactorZZ refinement must derive from the concrete
  selected candidate and source input;
- `HenselLiftRuntimeReadiness` contains only the universally quantified
  low-level lift-loop and normalization execution contracts.

`HenselLiftEntryInvariant` extends both structures, so the existing genuine
raw-execution refinement theorem retains exactly the same assumptions and
conclusion while making their provenance independently auditable.

## Why

The generated public FactorZZ theorem still exposes a single Hensel readiness
argument.  Several of its fields are consequences of the actual
`SelectionCorrect`/`SelectionPhysical` result or of successful source control
flow, while the lift and normalize fields describe genuine provider-level
execution readiness.  Separating the two classes is the prerequisite for
removing the derivable fields from the public theorem without weakening the
literal C++ L1 execution proof.

No output, tree, modulus, or factorization witness was introduced.  There is
no fuel, partial definition, semantic execution oracle, or new axiom.

## Verification

- Direct `CLPoly.Refinement.Hensel` compilation passed.
- `lake build CLPoly.Refinement.FactorZZ CLPoly.Refinement.Generated
  CLPoly.Pipeline.L1` passed: 3551 jobs.
- Production C++ changes: none, so this step has no new native regression or
  C++/Lean B2B change surface.

## Files

- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `proof/docs/devlog/2026-08-17-hensel-readiness-split.md`

## Next step

Construct the algebraic half from the actual selected prime result and source
execution, then change the generated public FactorZZ wrapper to request only
the irreducible low-level runtime readiness that cannot be derived there.
