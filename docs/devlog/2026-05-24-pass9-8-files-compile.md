# Pass 9: 8 refinement skeleton files compile

Date: 2026-05-24

## What

Rewrote `pass9_refine_gen.py` to generate correct L1→L2 refinement theorem skeletons for 8 functions across 4 files. All 3313 `lake build` jobs pass with 0 errors.

## Changes

### `proof/cpp2lean_v2/passes/pass9_refine_gen.py`
- **`_PARAM_TABLE`**: Replaced single-param inference (`_FIRST_PARAM`) with explicit per-function `[(name, type_str)]` dict, matching actual Corpus signatures. Fixed duplicate keys and removed unused entries.
- **`l1_call_params`**: Changed from comma-joined (Python style) to space-joined (Lean function application syntax).
- **`_emit_theorem_body`**: Added `l2_call` lookup in `REFINEMENT_MAP` — when present, used directly instead of auto-inferred params; avoids brittle `_first_sparse_poly_param` logic. For `pair_eq`/`map_eq_pair`, RHS is now wrapped in `SparsePolyZp.toPoly p` for type matching.
- Removed unused `_first_sparse_poly_param` function.

### `proof/cpp2lean_v2/class_map.py`
- `REFINEMENT_MAP` expanded from 2 → 9 entries with full metadata (`l1_name`, `l2_name`, `l2_call`, `l2_import`, `refinement_file`, `result_kind`, `cpp_source`, `doc`).
- `l2_call` explicitly specifies the L2 function invocation with correct argument names, bypassing fragile param inference.
- `result_kind` values: `direct_eq` (6 uses), `map_eq` (2), `map_eq_pair` (1), `pair_eq` (1).

### `proof/lean/CLPoly/Model.lean`
Added three `ZZ` namespace stubs:
- `ZZ.symmetricMod (a m : ZZ) : ZZ` — symmetric modulo
- `ZZ.binomial (n k : ZZ) : ZZ` — integer binomial (stub returns 1)
- `ZZ.isqrtCeil (n : ZZ) : ZZ` — integer square root ceiling (stub returns n)

### `proof/lean/CLPoly.lean`
Imports now include `Refinement.ZZArith` and `Refinement.ZpArith`.

## Files generated
| File | Theorems | Status |
|------|----------|--------|
| `Refinement/DDF.lean` | `__ddf_Zp_ir_refines` | ✓ compiles |
| `Refinement/SquarefreeZp.lean` | `__squarefree_Zp_ir_refines` | ✓ compiles |
| `Refinement/ZZArith.lean` | `__symmetric_mod_ir_refines`, `__binomial_ir_refines`, `__isqrt_ceil_ir_refines` | ✓ compiles |
| `Refinement/ZpArith.lean` | `__make_zp_ir_refines`, `__upoly_make_monic_ir_refines`, `__upoly_divmod_ir_refines` | ✓ compiles |

## Key decisions
1. **`l2_call` over inference**: Each refinement entry now specifies the exact L2 call string, avoiding name collision bugs (e.g., `p` param shadowing the prime `p`) and type mismatches.
2. **Stubs instead of real impl**: `ZZ.binomial` and `ZZ.isqrtCeil` are minimal stubs (`Nat.choose`/`Nat.sqrt` not available in this Lean version); refinement skeletons use `sorry` pending real L1 body translation.
3. **Not included**: `EDF`, `FactorZZ`, `FactorZp` — need fundamentally different refinement forms (result buffers, theorem-with-args pattern), deferred to follow-up.

## Next
- Fill `sorry` proofs for each skeleton (verify L1 body matches L2 spec).
- Add `EDF`, `FactorZZ`, `FactorZp` entries with appropriate `result_kind`.
