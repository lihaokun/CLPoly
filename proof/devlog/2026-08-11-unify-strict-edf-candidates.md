# Unify strict EDF candidate branches

The generated odd-characteristic candidate is now connected to the certified
strict powmod/subtract-one/GCD pipeline. The source `Int` exponent is proved to
convert exactly to the nonnegative `Nat` exponent used by the strict powmod
boundary; no exponent is approximated or truncated.

`strictCandidateRun_factor` follows the generated characteristic test and
combines this proof with the characteristic-two mod/trace/GCD proof. For every
canonical random polynomial it produces an actual successful generated
candidate execution whose concrete output is canonical, monic, and divides the
EDF input.

Verification:

```text
lake env lean CLPoly/Refinement/EDF.lean
lake env lean /private/tmp/CheckEDFCandidates.lean
```
