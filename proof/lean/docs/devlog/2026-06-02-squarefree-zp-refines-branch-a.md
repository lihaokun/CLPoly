# 2026-06-02: __squarefree_Zp_ir_refines Branch A 结构完成

## 做了什么

完成了 `__squarefree_Zp_ir_refines` 中 Branch A（derivative = 0 → p-th root）的证明结构。

## 关键决策

- **使用 `simpa` + `rw` 处理递归定义**：`__squarefree_Zp_ir_def` 是递归 def，直接 `simp` 会导致无限递归。改用 `rw [← eq_1, ← __squarefree_Zp_ir]` 双向折叠。

- **辅助引理全部 admit**：17 个辅助引理（property preservation、loop 正确性、sqfZp 不变性等）全部暂时 admit，后续逐一定理填补。

- **不使用 calc 复杂展开**：calc 中的 `toPolyList` 映射步骤全部 admit，避免无穷递归。

## 涉及的 admits

- `extract_pth_root_toPoly_eq` — extract_pth_root ↔ contract p
- `extract_pth_root_wellFormed/allReduced/...` — 保持性
- `upoly_make_monic_wellFormed/allReduced/...` — 保持性
- `makeMonic_toPoly_normalize_eq` — makeMonic ↔ normalize
- `sqfZp_normalize_eq` — sqfZp 在 normalize 下不变
- `loop_squarefree_Zp_0_correct` — 循环语义（scale exponents）
- `h_empty` — 常数多项式分支
- `h_unfold` — 展开递归定义
- `hmem` — g.front! ∈ g.toList
- `toPolyList` calc 步骤
- `h_contract_deg_lt` 相关

## 后续

- Branch B（derivative ≠ 0 → Yun algorithm）尚未开始（`sorry`）
- 逐个填写以上 admits
- 需要证明 `sqfZp_smul`（数乘不变性）和 `sqfZp_normalize_eq`
