# Generated SQF derivative-zero branch

The outer generated `__squarefree_Zp_raw_ir_state` derivative-zero branch is
now semantically composed.  Its execution follows the actual
`extractPthRootIR` and projected `upolyMakeMonicIR` calls, invokes the same
generated state entry recursively at the strictly smaller generated measure,
and runs the generated multiplicity-copy loop.

The result proof identifies the raw p-th root with polynomial contraction,
derives the multiplication no-wrap bound from source degree, and reuses the
proved concrete scaling loop theorem to match the characteristic-p branch of
`sqfZp`.  Neither the old handwritten `strictSquarefreeZpIR` entry nor an L2
execution oracle occurs on the executable side.

Verification:

```text
lake env lean CLPoly/Refinement/StrictSquarefreeGenerated.lean
```

The check passes.  The remaining outer branch is the nonzero-derivative path:
raw prefix, generated Yun loop, optional remainder contraction recursion, and
accumulator scaling.
