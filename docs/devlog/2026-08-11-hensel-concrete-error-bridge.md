# Hensel concrete error bridge

## Outcome

The generated `__hensel_step` error path now has a direct semantic bridge:

1. the actual heap multiplication output decodes to `g * h`;
2. the actual sparse subtraction output decodes to `f - g * h`;
3. coefficientwise exact divisibility by `m` justifies the generated GMP
   quotient loop;
4. reducing and compacting that quotient modulo `m` preserves the equation
   after multiplication by `m` in `ZMod (m^2)`.

Consequently the exact generated intermediate `e` satisfies
`f = g*h + m*e` modulo `m^2`.  No specification-chosen error polynomial is
introduced.  Explicit-output variants of the generated arithmetic and
modular-divmod refinement theorems were also added for threading the unique
values appearing in a successful raw execution.

## Verification

```text
lake env lean CLPoly/Refinement/Hensel.lean
```

The module compiles successfully.  The new bridge contains no `sorry`,
oracle, fallback, or fuel recursion.

## Progress

Strict Hensel refinement is approximately 90% complete.  Remaining work is
to thread the first concrete divmod and `tau/g/h` stores into the factor
algebra theorem, close the analogous Bezout correction, generate the public
entry theorem, and complete Pipeline/global audits.  Estimated remaining
focused time: 0.4-0.9 days.
