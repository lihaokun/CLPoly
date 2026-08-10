# Generated SQF entry refinement

The nonzero-derivative outer branch now executes the concrete raw GCD and
exact division prefix, the already refined generated Yun loop, and both
post-loop outcomes.  A positive remainder executes the real p-th-root and
make-monic calls, recurses through the same generated entry at a strictly
smaller measure, and appends the generated multiplicity-scaled result.  A
zero-degree remainder returns the generated accumulator directly.

`generatedSQFState_refines` combines both outer branches by strong induction
on the generated source measure.  The public theorem
`__squarefree_Zp_raw_ir_refines_sqfZp` now proves successful execution of
`Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir` and exact equality of
its decoded result with L2 `sqfZp`.  The theorem name retains the original C++
double underscores.

The public preconditions include positive source degree.  This matches the
actual recursive C++ entry's domain: calling `__squarefree_Zp` on a nonzero
constant would recurse on the same constant and is not a terminating source
execution.

Verification:

```text
lake env lean CLPoly/Refinement/StrictSquarefreeGenerated.lean
```

The file checks without placeholders.  The next step is to make Pass 9 emit
this generated raw contract, then run reproducibility, build, placeholder,
fallback, and axiom audits for SQF.
