# Define the concrete Hensel-step L2 contract

## Problem

`Algorithm.hensel_step` proves only that some lifted factors exist, while
`Spec.HenselCorrect` describes a whole factor list and its length.  Neither is
precise enough to state that the node returned by the generated C++
`__hensel_step` is correct.

## Change

`Refinement.StrictHensel.HenselNodeInvariant` now decodes the four concrete
node arrays and requires both:

- `g * h = f` modulo the current modulus;
- `s * g + t * h = 1` modulo the current modulus.

`Refinement.StrictHensel.HenselStepCorrect` requires the concrete output to
satisfy that invariant modulo `m^2`, while each of `g`, `h`, `s`, and `t`
reduces to the corresponding concrete input polynomial modulo `m`.

This is the target of the strict raw refinement theorem.  It does not choose
an existential L2 result and does not execute an L2 fallback.

## Verification

```text
lake build CLPoly.Refinement.Hensel
Build completed successfully (1923 jobs).
```
