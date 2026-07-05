# 2026-06-03: __squarefree_Zp_ir_safe 重构 + 关键引理

## 做了什么

1. **safe wrapper 重构**：`h_main` 改用 `__squarefree_Zp_ir_safe` 作为强归纳结论
2. **恢复关键引理**：`upoly_make_monic_allReduced` 完整证明（含 hp_size + UInt64 roundtrip）
3. **新增辅助引理**：`drop_eq_get_cons'`, `toPoly_push`, `listSum_singleton`

## 当前状态

- 构建通过（3464 jobs）
- Branch A 和 Branch B 仍为 `sorry`
- 原始文件中的其他 `admit` 保持不变

## 教训

git checkout 撤销了所有未提交的修改。完成的重要工作应及时提交。
