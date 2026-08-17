# Prove the strict odd EDF candidate is a factor

The exact odd-characteristic pipeline already returned the normalized GCD of
the generated powmod result minus one and the current polynomial.
`strictOddCandidateIR_factor` now derives the recursive facts required by
`EDFSplitLaw` for that same raw output: canonical representation, monicity,
and divisibility of the current polynomial.

The proof uses `gcd_dvd_right` and the previously established execution
equality. It does not replace the returned candidate with a mathematical GCD
witness.

Verification:

```text
lake build CLPoly.Refinement.EDF
```
