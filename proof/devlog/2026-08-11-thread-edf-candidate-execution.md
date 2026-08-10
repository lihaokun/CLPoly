# Thread the concrete EDF candidate execution into split preservation

The proper-degree guard alone does not prove that the candidate divides the
current polynomial. The previous operation-record layout also hid the exact
`candidateRun` execution from its recursive split proof.

The generator now separates computational operations (`EDFRawOps`) from the
proof layer (`EDFSplitLaw`). `splitStep` receives the successful
`candidateRun = .ok g` evidence stored in `SplitState`, as well as the proper
guard and the exact division/make-monic executions. The strict implementation
must therefore derive divisibility from the real candidate pipeline before it
can establish recursive invariants.

Verification:

```text
python3 ../cpp2lean_v2/tests/build_strict_edf.py --check
lake build CLPoly.Refinement.EDF
```
