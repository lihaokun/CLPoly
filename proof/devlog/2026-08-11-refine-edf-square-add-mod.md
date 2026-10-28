# Refine the characteristic-two EDF square/add/mod step

`strictSquareAddModIR` implements one source trace-map iteration in its exact
order: verified raw sparse multiplication for `g*g`, sparse polynomial
addition with `r`, and verified raw division/remainder reduction by `f`.

`strictSquareAddModIR_refines` proves that the concrete returned array is
canonical and denotes `(g^2 + r) %ₘ f`. It uses the physical multiplication
and modular-reduction workspaces from `DDFRawProviders`; no L2 result is used
as an executable replacement.

Verification:

```text
lake build CLPoly.Refinement.EDF
```
