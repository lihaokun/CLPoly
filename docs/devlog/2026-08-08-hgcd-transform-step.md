# HGCD 变换不变式单步

日期：2026-08-08

## 做了什么

- 定义 C++ HGCD 矩阵的精确方向：原始多项式对等于矩阵乘以当前 Euclid 多项式对。
- 证明一次除法恒等式和源码两次行更新保持该变换关系。
- 证明两次源码行更新使 2×2 行列式符号翻转。
- 加强零分支执行定理，给出交换后两个矩阵描述符的精确指针和长度。
- 闭合零分支语义：实际描述符交换表示 `[entry1 + quotient*entry0, entry0]`。

## 为什么做

HGCD 循环的 raw 执行只有在每轮都保持“原始对—矩阵—当前对”的关系时，才能最终连接到 GCD 及 SQF。零分支也必须从真实 L1 表示推导，不能当作特殊规格捷径跳过。

## 关键决策

- 采用源码实际的右乘方向：每一行 `[u,v]` 更新为 `[v+Q*u,u]`，吸收 `A=Q*B+R`。
- 零分支中的消失乘积由零长度 `SlicePolyRep` 推出 quotient 或 entry 为零，而不是增加 L2 前提。
- 行列式翻转单独证明，后续将与状态中的 `sgn := -sgn` 对齐。

## 涉及文件

- `proof/lean/CLPoly/Impl/StrictHGCDRefinement.lean`
- `proof/lean/CLPoly/Impl/StrictHGCDRawRefinement.lean`
- `proof/cpp2lean_v2/tests/check_strict_hgcd_source.py`
- `docs/devlog/2026-08-08-hgcd-transform-step.md`
