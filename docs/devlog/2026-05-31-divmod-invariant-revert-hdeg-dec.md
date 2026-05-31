# divmodAux_invariant 用 Nat 强归纳重构 + h_deg_dec 清理

## 日期
2026-05-31

## 做了什么
1. **保留 `h_deg_dec` 为 admit**：经过对 `divmodAux'_wf.induct` 的实验，确定 `h_deg_dec`（首项抵消导致的度递减）需要 `a * a⁻¹ = 1` 在 Zp 中的正确性证明（`modInv` / `extGcdAux` 的 Bezout 等式），这是一个独立的、非平凡的证明工作。
2. **尝试 `divmodAux'_wf.induct`**：`divmodAux'_wf` 内置的递归归纳原理本可完全避免 `h_deg_dec`，但因 `let` 绑定作用域在 `induct` case3 中与用户定义的 `let` 冲突，导致 `hred_r'` 被误识别为 `ℕ`。该路径暂不可行。
3. **回滚到 `Nat.strong_induction_on`**：保留了现有的强归纳结构，`h_deg_dec` 保持为干净的 `admit`。
4. **改进了 case 证明**：case 1 和 case 2 使用 `rw [divmodAux'_wf.eq_1]` 代替 `unfold` 展开 `divmodAux'_wf`。

## 为什么做
`h_deg_dec` 是 divmod 循环不变量证明中唯一剩余的子目标。填它需要：
- `Zp` 中 `a * a⁻¹ = 1` 的引理（`modInv` 正确性）
- 以及 `mergeAdd` 的 zero-coefficient 删除性质

## 关键决策
- `h_deg_dec` 暂时 admit，不阻塞其他工作
- `divmodAux'_wf.induct` 因 `let` 作用域问题暂时放弃

## 涉及的文件
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`

## 当前状态
- 构建零错误
- `h_deg_dec` 保留 admit（单个 admit）
- 其余 17 个 admit/sorry 分布在文件各处
