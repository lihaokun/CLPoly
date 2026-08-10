# Generated SQF Yun termination and post-root bridge

The generated well-founded Yun loop now returns a dependent final state whose
`done` field proves that the original C++ loop guard is false.  Its source
polynomial is an index of the loop state, so every recursive transition
preserves the exact top-level `__squarefree_Zp` input rather than an unrelated
existential witness.

`StrictSquarefreeGenerated.lean` strengthens the Yun invariant with the full
remaining `yunLoop` semantics.  The factor and no-factor transitions preserve
that continuation.  At termination, `remainderRootStep` uses the generated
false guard to prove that `w` has degree zero, identifies the concrete `c`
with the L2 remainder, and closes the real `extractPthRootIR` and
`upolyMakeMonicIR` executions with strict source-measure descent.

The concrete `strictSQFRawOps` bundle now contains only checked L1 operations:
`derivativeIR`, `extractPthRootIR`, `upolyMakeMonicIR`, the raw object GCD, and
`pairVecDivIR`.  It contains no call to `sqfZp`, no L2 execution fallback, and
no fuel parameter.

Verification:

```text
python3 ../cpp2lean_v2/tests/build_strict_squarefree.py
lake build CLPoly.Generated.StrictSquarefreeZp
lake env lean CLPoly/Refinement/StrictSquarefreeGenerated.lean
```

All commands pass.  The remaining SQF work is to prove total execution and
the result relation for the generated Yun loop and top-level well-founded
entry, then switch Pass 9 to the public generated contract named
`__squarefree_Zp_raw_ir_refines_sqfZp` and run the full audits.
