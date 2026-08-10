# Centralize final refinement contracts

All closed, public cpp2lean C++ L1 → Lean L2 contracts now appear directly in
`CLPoly/Refinement/Generated.lean`. The former per-stage modules below
`CLPoly/Refinement/Generated/` were removed so readers no longer need to search
an import-only index and a second directory layer.

Pass 9 now generates this single module deterministically. Both completed
contracts are generator-owned and retain the original double-underscore C++
entry names:

- `Refinement.__squarefree_Zp_raw_ir_refines_sqfZp`
- `Refinement.__ddf_Zp_raw_ir_refines_ddf`

Detailed proofs remain in their stage-specific refinement modules. A stage is
added to the centralized file only after its strict top-level theorem closes.

Verification:

```text
python3 -m py_compile proof/cpp2lean_v2/passes/pass9_refine_gen.py \
  proof/cpp2lean_v2/class_map.py
generated-output equality check: passed
lake build CLPoly.Refinement.Generated
lake env lean /private/tmp/CheckCentralGenerated.lean
```

Both public contracts use only `propext`, `Classical.choice`, and `Quot.sound`.
