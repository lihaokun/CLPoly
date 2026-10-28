# Hensel factor-correction algebra

## Outcome

Added the modulus-projection algebra needed to connect the concrete generated
factor-correction phase to the quadratic Hensel invariant:

- a coefficient in `ZMod (m^2)` that vanishes modulo `m` is annihilated by
  multiplication by `m`;
- the corresponding statement holds coefficientwise for polynomials;
- the exact error quotient, concrete modular divmod identity, concrete `tau`
  identity, and input Bezout identity imply that the generated updates
  `g + m*tau` and `h + m*r` multiply to the target modulo `m^2`.

The theorem accepts the semantic equations produced by the raw arithmetic
bridges.  It does not choose correction polynomials, invoke an L2 Hensel
implementation, or hide termination behind fuel.

## Verification

```text
lake env lean CLPoly/Refinement/Hensel.lean
```

The module compiles successfully without `sorry` or additional axioms.

## Progress

Strict Hensel refinement is approximately 88% complete.  The next step is to
thread the already-proved raw multiplication/addition/subtraction/divmod
equations through this algebra theorem for the generated factor phase, then
repeat the composition for the Bezout phase.  Estimated remaining focused
time: 0.5-1 day.
