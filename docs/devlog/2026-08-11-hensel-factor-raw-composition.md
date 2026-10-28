# Hensel factor raw composition

## Outcome

Completed the semantic composition of every raw operation in the first
contiguous generated `__hensel_step` phase:

- heap multiplication and sparse subtraction producing the exact error;
- generated coefficient division/reduction producing `e`;
- multiplication `s*e` and the exact well-founded modular divmod producing
  the concrete `q,r` arrays;
- the concrete `t*e + q*g` add path and coefficient reduction producing
  `tau`;
- concrete `g + m*tau` and `h + m*r` adds;
- coefficient stores modulo `m^2` producing `gNew,hNew`.

The resulting theorem proves directly that `gNew*hNew = f` modulo `m^2` and
that both new factors reduce to their input arrays modulo `m`.  Every named
intermediate is constrained by its generated raw-run equality, so none can be
selected by the proof as a specification oracle.

Explicit-output coefficient-reduction bridges were added to preserve the
unique raw values while moving between modulus `m^2` and modulus `m`.

## Verification

```text
lake env lean CLPoly/Refinement/Hensel.lean
```

The module compiles successfully without `sorry`, oracle, L2 fallback, or
fuel recursion.

## Progress

Strict Hensel refinement is approximately 93% complete.  Remaining work is
to extract these raw-run equations from the generated factor-phase execution,
compose the analogous Bezout correction, publish the generated final theorem,
wire Pipeline, and run global audits.  Estimated remaining focused time:
0.3-0.7 days.
