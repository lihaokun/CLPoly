# REFINEMENT_MAP 扩展：__upoly_mod 及精化证明填写

## 做了什么

1. 在 `REFINEMENT_MAP` 中添加 `__upoly_mod`（class_map.py）
   - L1 `__upoly_mod_ir` → L2 `(SparsePolyZp.divmod f g).snd`
   - `result_kind: "direct_eq"` — 两侧都是 SparsePolyZp，无需 `toPoly` bridge

2. 对应 `_PARAM_TABLE` 更新（pass9_refine_gen.py）

3. 填写了 3 个 trivial 精化证明（ZpArith.lean）：
   - `__make_zp_ir_refines` — L1 直接调用 `Zp.ofInt`，`rfl`
   - `__upoly_divmod_ir_refines` — L1 直接调用 `pair_vec_div5` → `SparsePolyZp.divmod`
   - `__upoly_mod_ir_refines` — 同上，取 `.snd`

4. 更新了 AGENTS.md（工作规范）

## 为什么做

为 L1→L2 精化框架建立第一个已证明的原子精化定理。`__upoly_mod_ir_refines` 是最简单的底层函数，可作为其余精化定理的基础构件。

## 关键决策

- 仅添加有明确 L2 等价的入口（`__upoly_mod`）。其余 Zp 侧函数（`__upoly_powmod`、`__extract_pth_root`、`__upoly_subtract_one` 等）无直接 L2 等价，不由 REFINEMENT_MAP 覆盖，改由调用方精化定理内联处理。
- `__upoly_mod` 使用 `direct_eq` 而非 `map_eq`，因此**不加入** `_SPARSE_POLY_FUNCS` 表，避免生成多余的 `SparsePolyZp.toPoly` 包裹。

## 涉及文件

- `proof/cpp2lean_v2/class_map.py` — REFINEMENT_MAP 新增 `__upoly_mod`
- `proof/cpp2lean_v2/passes/pass9_refine_gen.py` — _PARAM_TABLE 新增 `__upoly_mod`
- `proof/lean/CLPoly/Refinement/ZpArith.lean` — 3 个精化定理已证明，1 个仍为 sorry
- `proof/lean/CLPoly/Refinement/DDF.lean` — 骨架（未变）
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean` — 骨架（未变）
- `proof/lean/CLPoly/Refinement/ZZArith.lean` — 骨架（未变）
- `proof/lean/CLPoly/Refinement/Basic.lean` — bridge 引理（已有）
