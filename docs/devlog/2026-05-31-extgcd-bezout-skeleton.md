# extGcdAux Bezout 恒等式骨架 + 状态总结

## 日期
2026-05-31

## 做了什么
1. **添加 `extGcdAux_bezout` 引理**：证明 `Zp.extGcdAux A B 0 1` 的返回值 `(g, s')` 满足 `s'*B + t*A = g`（Bezout 恒等式）。保留了 `admit`，附带了详细的规范。
2. **分析 extGcdAux 的不变量**：追踪发现不变量为 `old_s * A + s * B = r`（`A,B` 为原始输入）。证明需要用 `extGcdAux.induct` 做强归纳，并处理 base case（`r=0`）下 `old_s*A + s*B = 0` 到 `old_s*B + s*A = old_r` 的转换。

## 关键技术细节
- `extGcdAux` 是 `Zp` 命名空间内的 `Int` 扩展欧几里得算法，用 `r.natAbs` 良基递归
- `extGcdAux.induct` 提供了归纳原理，但 `motive` 只包含 `old_r r old_s s`，不包含原始 `A,B`
- 递归步骤的不变量保持：`s*A + (old_s - q*s)*B = old_r - q*r` 由 `old_s*A + s*B = r` 和 `q = old_r / r` 推出
- 需要添加 `extGcdAux_bezout` 到 `h_deg_dec` 的桥梁：`(a * a.inv).toZMod p = 1`

## 当前 admit/sorry 分布

| 类别 | 数量 | 关键阻塞 |
|------|------|----------|
| divmod/gcd | ~10 | h_deg_dec (extGcdAux 不变量) |
| squarefree | ~9 | extract_pth_root_toPoly_eq, Branch A/B |
| 辅助 | 2 | get_deg_toPoly, extGcdAux_bezout |
| 未使用 | 1 | toPolyList_scaleExp |

## 涉及的文件
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`
