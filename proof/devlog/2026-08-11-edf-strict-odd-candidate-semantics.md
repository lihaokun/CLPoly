# Refine the concrete EDF odd candidate pipeline

`strictOddCandidateIR` now composes only verified raw C++ boundaries:

1. DDF's strict generated multiply/reduce powmod implementation;
2. EDF's generated `__upoly_subtract_one_raw_ir`;
3. DDF's strict generated public GCD implementation.

Its refinement theorem proves successful execution, canonical output, and the
exact L2 candidate formula
`normalize (gcd ((r^((p^d-1)/2) mod f) - 1) f)`.

The physical providers select multiplication, reduction, and GCD workspaces;
they do not contain candidate values or semantic factorization results.

Files changed:

- `proof/lean/CLPoly/Refinement/EDF.lean`
