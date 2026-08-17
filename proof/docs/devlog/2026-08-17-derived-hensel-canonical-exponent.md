# Derive Hensel adjusted canonicality and exponent legality

## What changed

The strict Hensel refinement now proves two remaining entry facts directly
from concrete execution and representation invariants.

First, the literal `scaleZpCoeffs` map preserves the degree order and runtime
prime/reduced-residue conditions of a canonical sparse modular polynomial.
The following source `normalization` filter removes exactly zero coefficients,
so its result is physically canonical.  Combined with the existing frame
theorem for positive array slots, this proves that the complete adjusted
factor array is canonical.

Second, a successful generated `__hensel_lift_upoly_raw_ir` execution proves
`aTarget = 0 ∨ aTarget > 0`.  A negative `Int32` target exponent executes the
literal assertion-failure branch of `__hensel_explicit_target_raw_ir`, so it
cannot be hidden behind an external source-contract premise.

## Why

Both properties were previously bundled into `HenselLiftEntryInvariant` and
therefore exposed through the generated FactorZZ theorem's readiness
argument.  They are now available for construction from actual upstream
physical factors and successful source control flow.

The canonicality proof is representation-level: it reasons about the actual
array map/filter and modular multiplication.  It does not replace the C++
operation with a semantic polynomial normalization oracle.

## Verification

- Isolated exponent-control-flow check passed and the temporary check file
  was removed.
- `lake env lean CLPoly/Refinement/Hensel.lean`: passed in full.
- No production C++ change; there is no new native/B2B change surface.

## Files

- `proof/lean/CLPoly/Refinement/Hensel.lean`
- `proof/docs/devlog/2026-08-17-derived-hensel-canonical-exponent.md`

## Next step

Assemble the complete internal algebraic invariant from selection/input/run
facts and change the generated public FactorZZ theorem to accept only the two
genuine runtime readiness fields for lift and normalize execution.
