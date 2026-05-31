# h_deg_dec 分析 + extGcdAux Bezout 恒等式骨架

## 日期
2026-05-31

## 做了什么
1. **尝试多种 h_deg_dec 策略**：
   - `divmodAux'_wf.induct`：内置递归归纳可避免 `h_deg_dec`，但 `let` 绑定在 `induct` case3 中与用户 `let` 作用域冲突，导致 `hred_r'` 被误认为 `ℕ`
   - `Nat.strong_induction_on`：最干净的方案，但需要 `h_deg_dec`（首项抵消）
2. **分析 h_deg_dec 依赖关系**：需要 `a * a⁻¹ ≡ 1 (mod p)`，即 `modInv(a.val, p)` 的正确性
3. **添加 `extGcdAux_bezout` 引理骨架**：使用 `Zp.extGcdAux.induct` 做归纳。已在文件中添加声明和 base case（`r=0`），recursive case 需要追踪原始 `(A,B)` 的归纳不变量
4. **发现关键性质**：`Zp.extGcdAux` 的不变量是 `old_s * A + s * B = r`（`A,B` 为原始输入），`motive` 需要携带 `A,B` 参数，使归纳复杂化

## 不变量说明
`Zp.extGcdAux old_r r old_s s` 中（设 `A = old_r, B = r` 为原始值）：
- 每一步有 `old_s * A + s * B = r`（当前余数）
- 终止时 `r = 0`，返回 `(old_r, old_s)`
- `old_s`（返回值 .2）是上一步的 `s`，它满足 `s * B + old_s_prev * A = gcd`
- 当 `gcd = 1`：`(old_s).toNat * B.toNat ≡ 1 (mod A.toNat)`，即 `modInv` 正确

## 涉及的文件
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`：添加 `extGcdAux_bezout` 引理

## 当前状态
- 构建零错误
- `h_deg_dec` 仍为 admit
- `extGcdAux_bezout` 已作为 admit 添加（归纳证明非平凡，需追踪原始 `A,B`）
- 其余 admit 不变
