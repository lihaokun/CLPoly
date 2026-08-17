# Replace the EDF retry rank with an exact execution trace

The previous generated EDF shell accepted `EDFRetryLaw`, including an
arbitrary natural-valued retry rank supplied by its caller. Although each
decrease obligation mentioned concrete executions, that interface was too
abstract for the final C++ refinement boundary.

The generator now emits `RetryTrace`. Its constructors record every concrete
`ops.random` execution and, for nonempty random polynomials, every concrete
`candidateRun` execution. Empty and unsuccessful attempts contain the trace
from the advanced RNG state; a successful attempt contains the exact proper
factor test and terminal RNG state. `retryLoop` is structurally recursive on
this finite trace.

`EDFTermination.retryTrace` supplies such a trace at each recursive EDF node.
This expresses the correct conditional termination boundary of the source
`while (true)` loop without a bounded counter and without claiming that every
possible pseudo-random trajectory terminates. The refinement theorem
`retryLoop_terminates` is the corresponding raw-to-safe execution bridge.

Verification:

```text
python3 ../cpp2lean_v2/tests/build_strict_edf.py --check
lake build CLPoly.Refinement.EDF
```
