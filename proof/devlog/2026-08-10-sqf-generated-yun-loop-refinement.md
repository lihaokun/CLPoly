# Generated SQF Yun loop refinement

The proof now closes the complete cpp2lean-generated Yun state loop by strong
induction on its generated well-founded measure.  Each active iteration
executes the concrete raw object GCD, both `pairVecDivIR` calls, and the source
make-monic call when its C++ guard is true.  Recursive use is justified only
by the strict transition measure already carried by `SQFRawOps`.

The invariant target is no longer existential.  `sourceYunTarget` defines the
unique proof-level Yun result associated with the top-level SQF source, and
every generated transition preserves equality to that target.  This function
does not occur in any raw operation field and therefore cannot serve as an L2
execution fallback.

`generatedYunLoopState_refines` proves both successful raw execution and exact
decoding of the final result array and remainder to `sourceYunTarget`.  It
uses the final generated false guard to unfold the L2 Yun loop at degree zero.

Verification:

```text
lake env lean CLPoly/Refinement/StrictSquarefreeGenerated.lean
git diff --check
```

Both checks pass.  The next step is the outer well-founded
`__squarefree_Zp_raw_ir_state` refinement, including multiplicity scaling and
recursive p-th-root calls.
