# Define the strict recursive EDF invariant

`EDFEntryInvariant` now records the obligations required at every concrete
recursive `__edf_Zp` call: canonical sparse representation, matching modulus,
signed-degree safety, monicity, positive degree, positive target degree,
squarefreeness, and the equal-degree condition on every irreducible divisor.

Derived lemmas prove that an invariant input is nonempty and that its L2
natural degree fits the signed 64-bit source guard. These facts will discharge
the generated base guard and provide the premises for the exact division,
make-monic, and recursive split proofs.

Verification:

```text
lake build CLPoly.Refinement.EDF
```
