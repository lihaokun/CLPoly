# Hensel Bezout raw composition

## Outcome

Completed the semantic composition of every concrete raw operation in the
second contiguous generated `__hensel_step` phase:

- `s*g`, `t*h`, and both sparse subtractions producing the Bezout error;
- exact coefficient division/reduction producing `ep`;
- `s*ep` and the second finite well-founded modular divmod;
- the concrete `s + m*r` update and modulus-`m^2` store;
- concrete `t*ep + q*g`, its coefficient reduction to `tau2`, and the
  `t + m*tau2` modulus-`m^2` store.

The theorem proves that the concrete `sNew,tNew` establish
`sNew*g+tNew*h=1` modulo `m^2` and both reduce to their input values modulo
`m`.  Every intermediate is tied to a generated raw-run equation; the only
input semantic assumption is the existing Bezout component of the node
invariant modulo `m`.

## Verification

```text
lake env lean CLPoly/Refinement/Hensel.lean
```

The module compiles successfully without `sorry`, oracle, L2 fallback, or
fuel recursion.

## Progress

Strict Hensel refinement is approximately 97.5% complete.  Remaining work is
to quantify the deterministic locals inside the complete generated Bezout
phase, compose the two generated phases into the final public entry theorem,
and complete Pipeline/global audits.  Estimated remaining focused time:
0.15-0.4 days.
