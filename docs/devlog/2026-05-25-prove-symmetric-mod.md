# 2026-05-25: 证明 __symmetric_mod_ir_refines

## 做了什么

1. **修复 `ZZ.symmetricMod` L2 模型**：将 `a % m` 改为 `a.fmod m`（`Int.fmod`），与 C++ 翻译（使用 `Int.fmod`）一致。原实现用 `a % m` 在负模数时与 `Int.fmod` 产生不同结果。
2. **证明 `__symmetric_mod_ir_refines`**：使用 `Int.fdiv_eq_ediv_of_nonneg` 将 `Int.fdiv m 2` 转为 `m / 2`，然后 `omega` 完成条件等价性证明。
3. **尝试证明 `__upoly_make_monic_ir_refines`**：
   - 成功证明了循环不变量 `loop_toList_invariant`（使用 `measure` well-founded induction）
   - 成功证明了 `loop_toPoly`（从循环不变量导出）
   - 受阻于 `toPoly_scalarMul`：`scalarMul` 使用 `filterMap`（去除零系数项），与 C++ 循环的 `Array.set!` 行为不同。需要证明 `listSum(filterMap ...) = listSum(map ...)`，这需要 `x.prime = c.prime`。从 `WellFormed` 可推导 `x.prime = UInt64.ofNat p`，但证明 `c.prime = UInt64.ofNat p`（从 `hc_prime : c.prime.toNat = p`）需要 `UInt64.eq_of_toNat_eq` 等不存在的旧版 API。

## 当前状态

| 定理 | 状态 |
|------|------|
| `__symmetric_mod_ir_refines` | ✅ 已证明 |
| `__binomial_ir_refines` | ❌ L2 模型不完整 (`ZZ.binomial` return 1) |
| `__isqrt_ceil_ir_refines` | ❌ L2 模型不完整 (`ZZ.isqrtCeil` return n) |
| `__upoly_make_monic_ir_refines` | ❌ 受阻于 UInt64 API |
| `__squarefree_Zp_ir_refines` | ❌ 依赖诸多下层 refinement |
| `__ddf_Zp_ir_refines` | ❌ 依赖诸多下层 refinement |
| Pipeline 定理 | ❌ 依赖 sqf/ddf refinement |

## 涉及文件
- proof/lean/CLPoly/Model.lean
- proof/lean/CLPoly/Refinement/ZZArith.lean
