# Generated SQF contract and audit

Pass 9 now emits `__squarefree_Zp_raw_ir_refines_sqfZp`.  The generated public
contract imports `StrictSquarefreeGenerated` and its executable side is
`Generated.StrictSquarefreeZp.__squarefree_Zp_raw_ir` instantiated with the
concrete raw operation bundle.  The previous public wrapper around the
handwritten `strictSquarefreeZpIR` is no longer generated.

The contract includes the positive-degree domain precondition of the actual
recursive C++ entry.  The original function would recurse without descent on
a nonzero constant, while all pipeline call sites supply positive-degree
factors.

Two finite word-size facts in the reused raw prefix theorem and all such facts
in the new generated refinement were changed from `native_decide` to kernel
checked `norm_num`.  The final axiom report contains only the standard
mathlib foundations `propext`, `Classical.choice`, and `Quot.sound`; there are
no generated native-decision axioms.

Audit evidence:

```text
python3 ../cpp2lean_v2/tests/build_strict_squarefree.py --check
lake build CLPoly.Refinement.Generated.SquarefreeZp
lake env lean /private/tmp/CLPolySQFAxiomAudit.lean
```

Pass 9 output also compares byte-for-byte with `_emit_verified_contract`.
Searches find no `sorry`, `admit`, declared `axiom`, `partial def`, `fuel`, old
`strictSquarefreeZpIR`, or old public theorem name in the generated SQF
contract path.
