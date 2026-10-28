# Certify the strict EDF characteristic-two candidate

`strictCandidateRun_charTwo_factor` now follows the generated characteristic-
two branch through its actual raw modular reduction, well-founded trace loop,
and raw GCD call. The concrete returned GCD output is canonical, monic, and
divides the recursive EDF input.

The theorem uses `strictEDFRawOps`; it neither manufactures an L2 factor nor
replaces a failed raw call. The trace loop is the generated loop already proved
terminating by the exact source distance `d-i`.

Verification:

```text
lake env lean CLPoly/Refinement/EDF.lean
lake env lean /private/tmp/CheckEDFCharTwo.lean
```
