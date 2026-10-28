# extGcdAux_bezout 证明完成

## 日期
2026-05-31

## 做了什么
成功证明了 `Zp.extGcdAux A B 0 1` 的 Bezout 恒等式。

## 证明结构

1. **关键引理 `extGcdAux_linearity`**：证明了 `extGcdAux` 的结果系数关于 `(old_s, s)` 线性：
   - `extGcdAux old_r r old_s s = (g, x*old_s + y*s)` 且 `x*old_r + y*r = g`
   - 对任意 `(old_s, s)` 成立，不仅限于 `(0, 1)`
   - 证明：强归纳于 `r.natAbs`

2. **主引理 `extGcdAux_bezout`**：由线性性直接推出：
   - `(extGcdAux A B 0 1).2 * B + u * A = (extGcdAux A B 0 1).1`
   - 取 `u = x`（标准 EEA 系数中的另一个）

## 线性性证明的归纳步骤
递归步：`extGcdAux old_r r old_s s = extGcdAux r r' s (old_s - q*s)`，其中 `r' = old_r % r`，`q = old_r / r`。
由 IH：`(r, r', s, old_s - q*s)` 的结果为 `(g', x'*s + y'*(old_s - q*s))` 且 `x'*r + y'*r' = g'`。
令 `x = y'`，`y = x' - q*y'`，`g = g'`。
验证：
- `x'*s + y'*(old_s - q*s) = y'*old_s + (x' - q*y')*s = x*old_s + y*s` ✓
- `x*old_r + y*r = y'*old_r + (x' - q*y')*r = y'*old_r + x'*r - q*y'*r`
  由 `r' = old_r - q*r` 和 `x'*r + y'*r' = g` 得 `x'*r + y'*(old_r - q*r) = g`
  → `y'*old_r + (x' - q*y')*r = g` ✓

## 涉及的文件
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`

## 当前状态
- 构建零错误
- `extGcdAux_bezout`：**已完成**
- `h_deg_dec`：admit（等待被 extGcdAux_bezout 填充）
- 其余 admit 不变
