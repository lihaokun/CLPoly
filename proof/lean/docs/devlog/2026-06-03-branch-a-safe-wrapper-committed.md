# 2026-06-03: Branch A safe wrapper + 关键引理已提交

## 完成的工作

1. safe wrapper 重构（`h_main` 使用 `__squarefree_Zp_ir_safe`）
2. `upoly_make_monic_allReduced` 完整证明
3. `drop_eq_get_cons'`、`toPoly_push`、`listSum_singleton`
4. Branch A 证明骨架（`h_eq_nd`、`h_contract_deg_lt`、`h_result` calc）
5. L1 pipeline 更新（使用 `__squarefree_Zp_ir_safe`、增加 `hp2` 参数）
6. Corpus.lean 补充 `_def` variants

## 当前状态

- 构建通过，所有修改已提交
- Branch A：`admit` 剩待填（`h_extract_eq`、`h_nd_g2`、IH 子条件、`h_loop_eq` 等）
- Branch B：`sorry`
- 原始文件 `admit`：`get_deg_toPoly`、`Zp_toZMod_inv_mul_self`、`listSum_derivative_eq`、`divmod'_wf_deg_lt`

## 关键教训

- `partial def` 不能被 `simp`/`dsimp` 展开，导致 `simpa`/`rw`/`change` 无法匹配展开前后的表达式
- 使用 `unfold` + `rw [eq_1]` 展开后，需用 `admit` 处理返回目标
