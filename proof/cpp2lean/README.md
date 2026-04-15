# cpp2lean — CLPoly C++ to Lean 4 translator

Translates CLPoly's polynomial factorization C++ code into Lean 4 IR (intermediate representation), for L1 implementation-level verification.

## Prerequisites

- **Clang** (with `-ast-dump=json` support, tested with clang 14+)
- **Python 3.10+**
- **libclang Python bindings**: `pip install libclang`
- CLPoly source tree (the translator runs from `proof/cpp2lean/` and expects `../../clpoly/` to exist)

## Quick start

```bash
cd proof/cpp2lean

# Full translation: all 66 functions (Zp + Univar + Wang modules)
python3 gen_full.py > output.lean 2>log.txt

# Zp module only (fast, 13 functions)
python3 gen_full.py --zp-only > output_zp.lean

# Check translation log
cat log.txt
```

## How it works

The translator is a 5-stage pipeline:

```
C++ source
    │
    ▼  [1] clang_hybrid.py — Clang AST JSON extraction
instantiate.cc  ──► clang++ -ast-dump=json -ast-dump-filter=FUNC
    │
    ▼  [2] clang_ast.py — AST → IR parsing
FuncIR (statements + expressions)
    │
    ▼  [3] ssa_transform.py — SSA transformation
SSAFunc (let-chain, loops extracted as standalone functions)
    │
    ▼  [4] ub_collector.py — UB obligation injection
SSAFunc + Require nodes (division-by-zero guards etc.)
    │
    ▼  [5] lean_codegen.py — IR → Lean 4 code generation
Lean 4 source (partial def functions)
```

### Key design decisions

- **Loop-as-function**: All C++ loops (for, while, range-for) are extracted as standalone `partial def` functions with explicit state passing, not inline `let rec`.
- **SSA transform**: C++ mutable variables become Lean let-chains with version-bumped names (`x_1`, `x_2`, ...).
- **Output parameters**: C++ `void f(T& out)` is transformed to `def f(...) : T` returning the output value.
- **Type mapping**: All C++ types are resolved from Clang AST, not inferred. See `class_map.py` for the complete mapping table.
- **Trusted base**: `clpoly_model.lean` defines the Lean-side type models (Zp, SparsePolyZp, MvPolyZZ, etc.) that the generated code imports.

## File overview

| File | Role |
|------|------|
| `gen_full.py` | Main entry point, orchestrates the pipeline |
| `instantiate.cc` | Forces C++ template instantiation for concrete types |
| `clang_hybrid.py` | Clang AST extraction (libclang for discovery, JSON for parsing) |
| `clang_ast.py` | AST node → IR node translation with dispatch tables |
| `ssa_transform.py` | SSA transformation, loop extraction, closure capture |
| `lean_codegen.py` | IR → Lean 4 code generation |
| `class_map.py` | Type/method/function mapping tables (`CLASS_MAP`, `FUNC_MAP`, `TRANSLATION_SCOPE`) |
| `ir_types.py` | IR node dataclass definitions |
| `clpoly_model.lean` | Lean trusted base (type definitions + primitive operations) |
| `ub_collector.py` | UB obligation detection (div-by-zero, array OOB) |

## Translation scope

66 functions across 3 C++ headers:

- `polynomial_factorize_zp.hh` — 13 functions (Zp field operations, DDF, EDF)
- `polynomial_factorize_univar.hh` — 34 functions (Hensel lifting, LLL, recombination)
- `polynomial_factorize_wang.hh` — 19 functions (MTSHL, Wang leading coeff, eval point selection)

5 dead-code functions (classical Wang Hensel path, replaced by MTSHL) are excluded from translation scope. See `TRANSLATION_SCOPE` in `class_map.py`.

## Configuring translation scope

Edit `class_map.py`:

- **`TRANSLATION_SCOPE`**: Set of function names to translate
- **`TRANSLATION_SCOPE_OUTPUT_PARAMS`**: Output parameter indices per function
- **`CLASS_MAP`**: C++ class → Lean type + method mappings
- **`FUNC_MAP`**: Standalone function mappings

## Design documents

- `docs/design/l1-translation-validation/loop-extraction-design.md` — Loop-as-function architecture
- `docs/design/l1-translation-validation/ast-dispatch-design.md` — AST dispatch mechanism
- `docs/design/l1-translation-validation/mutation-model-design.md` — Mutation model (G1-G6)
- `docs/design/l1-translation-validation/class-model-design.md` — Class/type model design
