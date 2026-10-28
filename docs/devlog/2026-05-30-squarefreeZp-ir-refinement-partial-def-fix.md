# 2026-05-30: SquarefreeZp `partial def` 问题修复与精化引理框架

## 问题

`__squarefree_Zp_ir_refines` 主定理在 `SquarefreeZp.lean` 中无法编译，原因：

1. **`partial def` 透明性问题**：`__squarefree_Zp_ir` 位于 `Generated/Corpus.lean` 的 6908 行 `mutual` 块内，定义为 `partial def`。递归 `partial def` 在 Lean 4 中完全透明（`unfold`/`dsimp`/`.eq_1` 均失败）。

2. **引理签名错误**：`normalization_toPoly_eq` 将 `SparsePolyZp.normalization`（零系数过滤）等同于 `Polynomial.normalize`（monic 化），实际两者语义不同，旧签名 false。

## 修复

### Corpus.lean 结构改动

- 在 `mutual` 块之前（line 100 前）插入 `_def` 变体（`def` + `termination_by`），共 7 个函数
- 将 `mutual` 块内对应 `partial def` 改为非递归 wrapper（`__squarefree_Zp_ir f := __squarefree_Zp_ir_def f`）
- 非递归 `partial def` 可以被 `unfold`

### 修复的引理

| 引理 | 状态 | 说明 |
|------|------|------|
| `toPolyList_append` | ✅ 完成 | `simp [toPolyList, Array.map_append]` |
| `toPolyList_scaleExp` | ✅ 完成 | 需 `h_no_overflow` 条件避免 UInt64 溢出 |
| `normalization_toPoly_eq` | ✅ 完成 | 改为 `toPoly p (normalization f) = toPoly p f` |
| `upoly_make_monic_snd_toPoly_eq` | ✅ 完成 | 委托 `__upoly_make_monic_ir_refines` |
| `derivative_toPoly_eq` | ❌ sorry | 需 `SparsePolyZp.derivative` 与 `Polynomial.derivative` 的 toPoly 对应 |
| `gcd_toPoly_eq` | ❌ sorry | 需 `polynomial_GCD` 与 `EuclideanDomain.gcd` 的对应 |
| `pair_vec_div_toPoly_eq` | ❌ sorry | 需 `pair_vec_div` 与 `/ₘ` 的对应 |
| `extract_pth_root_toPoly_eq` | ❌ sorry | 需 `__extract_pth_root_ir` 与 `contract p` 的对应 |
| `yunLoop_toPoly_list_eq` | ❌ sorry | 复杂，依赖 `_loop___squarefree_Zp_1_ir_def` (partial def) |
| `__squarefree_Zp_ir_refines` | ❌ sorry | 主定理 body（强归纳骨架已就位） |
| `get_deg_toPoly` | ❌ sorry | 简单但主定理不需要 |

### 关键决策

- **UInt64 溢出处理**：`toPolyList_scaleExp` 添加 `h_no_overflow` 条件。在 `__squarefree_Zp_ir` 的实际使用中，被 scale 的指数来自 Yun 循环的计数值（小整数），`p * exponent < 2^64` 在合理范围内成立。
- **Yun 循环**：`_loop___squarefree_Zp_1_ir_def` 仍为 `partial def`（终止性无法形式化），对应的 `yunLoop_toPolyList_eq` 暂留 sorry。

## 涉及文件

- `CLPoly/Generated/Corpus.lean` — `_def` 变体插入、wrapper 替换
- `CLPoly/Refinement/SquarefreeZp.lean` — 主定理骨架、辅助引理修复
- `CLPoly/Refinement/ZpArith.lean` — 适配 `unfold` + `rw [eq_1]` 模式
