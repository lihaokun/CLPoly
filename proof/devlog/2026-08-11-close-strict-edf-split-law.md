# Close the strict EDF split law

`strictEDFSplitLaw` now certifies the two recursive children created by a
successful generated EDF retry. The law retains and checks the actual random
engine transition, candidate execution, exact division, and both make-monic
executions.

The concrete candidate is proved canonical, monic, and a proper divisor. Its
exact quotient is canonical and monic; both children inherit squarefreeness and
the equal-degree irreducible-factor invariant. Their generated `edfMeasure`
values are proved strictly smaller via exact signed C++ degree-to-L2 degree
bridges.

The generated split-law interface now carries the successful random execution
instead of discarding it. This is necessary to derive canonicality of the
specific random polynomial passed to `candidateRun`. The strict EDF generator
was updated and its reproducibility/placeholder check passes.

Verification:

```text
python3 proof/cpp2lean_v2/tests/build_strict_edf.py --check
lake build CLPoly.Refinement.EDF
lake env lean /private/tmp/CheckEDFSplitLaw.lean
```
