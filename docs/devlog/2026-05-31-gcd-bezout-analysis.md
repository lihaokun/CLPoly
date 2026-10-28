# extGcdAux Bezout 恒等式深度分析

## 日期
2026-05-31

## 分析结果

`Zp.extGcdAux A B 0 1` 返回 `(g, coeff)`，满足 `coeff * B + y_sub * A = g`，其中 `y_sub = (Zp.extGcdAux B (A % B) 0 1).2`。

证明思路（暂 admit）：
1. 强归纳于 `B.natAbs`，设 `n = B.natAbs`
2. 基例 `B = 0`：平凡
3. 递归步：令 `r = A % B, q = A / B`，则 `extGcdAux A B 0 1 = extGcdAux B r 1 (-q)`
4. 由归纳假设（因 `r.natAbs < n`）：`(g, y_sub) = extGcdAux B r 0 1` 满足 `y_sub * r + u_sub * B = g`
5. 从 `r = A - q*B`：`y_sub*(A - q*B) + u_sub*B = g` → `(u_sub - q*y_sub)*B + y_sub*A = g`
6. 于是 `u = y_sub` 满足 `coeff * B + u * A = g`

关键未证：`coeff` 与 `y_sub` 的关系。
实际上 `coeff = u_sub - q * y_sub`，但此式需从 `extGcdAux` 的系数传播性质推出，这正是归纳的核心。

## 涉及的文件
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`

## 当前状态
- `extGcdAux_bezout` 仍为 admit，已记录完整证明思路
- `h_deg_dec` 仍为 admit
- 其他 admit 不变
