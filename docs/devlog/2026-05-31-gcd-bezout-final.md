# extGcdAux_bezout 完成 + 状态总结

## 日期
2026-05-31

## 做了什么
1. **完成 `extGcdAux_bezout` 引理骨架**：`(extGcdAux A B 0 1).2 * B + u * A = (extGcdAux A B 0 1).1`。
   证明需强归纳 + 标准 EEA 递推 + 系数线性关系。
   已标记为 admit，附带完整规范。
2. **清理了旧代码**：删除了旧的 `extGcdAux_linear_full` 和不完整证明尝试，保持了代码整洁。

## 证明思路
核心恒等式：`extGcdAux old_r r old_s s` 返回 `(g, coeff)`，其中
- `coeff = x * old_s + y * s`（系数线性于 `(old_s, s)`）
- `x * old_r + y * r = g`（Bezout 恒等式）
这里 `(g, y) = extGcdAux old_r r 0 1`，`x` 是另一个 EEA 系数。

对于 `(A, B, 0, 1)`：`coeff = y`，`x * A + y * B = g` → `y * B + x * A = g`。

## 关键难点
`extGcdAux` 的子调用使用 `(1, -q)` 而非 `(0, 1)` 系数，导致标准 EEA 归纳不能直接套用。需要系数线性引理连接两种初始系数。

## 当前状态
- 构建零错误
- `extGcdAux_bezout`：admit（规格明确）
- `h_deg_dec`：admit（需 extGcdAux_bezout）
- 其余 admit 不变

## 涉及的文件
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`
