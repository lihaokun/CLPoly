# Hensel phase execution bridges

## Outcome

The generated strict `__hensel_step` body is now split at the exact boundary
between its contiguous factor-correction and Bezout-correction statement
ranges.  The public raw entry still binds those two ranges in the original
C++ operation order; no operation was reordered, replaced, or omitted.

Both generated phases and the complete entry now have raw-to-safe execution
theorems.  The complete bridge consumes a safety-only invariant recording the
two concrete modular-division assertions.  It contains no semantic output,
oracle, fuel counter, or L2 fallback.

The refinement layer also gained:

- modular decoding of the generated coefficient-scaling loop at an arbitrary
  target modulus;
- projection of a decoded sparse polynomial from a modulus to any divisor;
- proof that coefficient reduction modulo a stronger modulus preserves the
  polynomial at every divisor modulus.

These are required to prove that the `m^2` writes made by the real C++ step
both establish the quadratic invariant and reduce to the original fields
modulo `m`.

## Verification

```text
python3 proof/cpp2lean_v2/tests/build_strict_hensel.py --check
lake env lean CLPoly/Refinement/Hensel.lean
```

Both checks pass.  No `sorry`, `admit`, oracle, fuel implementation, or L2
fallback was added.

## Progress

Strict Hensel refinement is approximately 86% complete.  Remaining work is
the semantic factor-phase proof, the semantic Bezout-phase composition, the
generated public `__hensel_step` contract, Pipeline wiring, and the final
repository-wide build/placeholder/axiom audit.  Estimated remaining focused
time: 0.5-1.25 days.
