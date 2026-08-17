# 2026-05-25 ZpArith 重构：listSum_map_mul_prime_only

## 做了什么

- 将 `listSum_map_mul` 重命名为 `listSum_map_mul_prime_only`
- 将其参数从 `hc_red : Zp.Reduced p c` 简化为 `hc_prime : c.prime.toNat = p`
- 填充 `__upoly_make_monic_ir_refines` 的证明（后因 well-founded induction 的 `Array` API 问题回退为 `sorry`）

## 为何回退

`loop_toPoly` 需要 well-founded induction 证明 `_loop___upoly_make_monic_0_ir` 的语义等价于 `scalarMul`。该 loop 使用 `Array.set!` 原地修改，与 `scalarMul` 的 `filterMap` 实现不同。well-founded induction 的 `measure` 方案因 `arr` binder 在 `revert` 后不可用而失败。此外 `Array.mem_iff_get`、`Array.size_pos_of_nonempty` 等 API 在所用 Lean 版本中不存在。

## 涉及文件

- proof/lean/CLPoly/Refinement/ZpArith.lean
