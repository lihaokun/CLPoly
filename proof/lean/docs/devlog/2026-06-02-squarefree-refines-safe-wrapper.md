# 2026-06-02: __squarefree_Zp_ir_safe 重构 + 常数 case

## 做了什么

1. 将 `h_main` 从使用 `__squarefree_Zp_ir` 改为使用 `__squarefree_Zp_ir_safe`，避免 `partial def` 在常数多项式上不停机
2. 常数 case 直接用 `simp` 证明（`toPoly g = 0` 时两边均为 `[]`）
3. 通过 `simpa`/`calc`/`rw` 处理了所有结构性证明错误
4. 移除了不可证的 `h_empty` admit（`__squarefree_Zp_ir g = #[]`）

## 关键决策

- **使用 `__squarefree_Zp_ir_safe`** 作为强归纳的结论，而非 `__squarefree_Zp_ir`。这样常数 case 自然成立，无需证明 `partial def` 在非终止输入上的行为

## 剩余工作

21 个 admits 按重要性排列：
1. `loop_squarefree_Zp_0_correct` — 循环语义（中等难度）
2. `extract_pth_root_toPoly_eq` — 核心代数引理（高难度）
3. `extract_pth_root_wellFormed/allReduced/deg_bound/no_overflow` — 4 个保持性
4. `makeMonic_toPoly_normalize_eq` — makeMonic ↔ normalize
5. `sqfZp_normalize_eq` — sqfZp 不变性
6. `h_unfold` + calc_step + `hmem` — 辅助
7. Branch B — Yun 算法
