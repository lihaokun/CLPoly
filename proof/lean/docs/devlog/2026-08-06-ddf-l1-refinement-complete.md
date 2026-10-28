# DDF L1→L2 精化完成

日期：2026-08-06

## 完成范围

- 二进制模幂的 total 循环与入口：规范形保持及 `pow %ₘ modulus` 语义。
- `subtract X`：由规范稀疏减法实现，证明排序、非零系数、约化性和多项式语义。
- GCD、首一精确除法、商规范化和余式更新的 L1→L2 桥接。
- total DDF 主循环：按
  `natDegree(fStar) + 1 - 2*d` 终止，并证明 split、no-split、终止三个路径与
  `Algorithm.DDF.ddfLoop` 一致。
- 防御性 `hdec` 失败路径在合法不变量下不可达。
- `ddfZpSafe_refines` 与 Pipeline 的 `ddf_l1_correct` 已闭合；超过 UInt64 次数范围的
  纯数学输入在 Pipeline 中显式回退到 L2。

## 证明边界

生成语料库的原始 DDF 和若干循环位于大型 `mutual partial` 块。Lean 不为这些 opaque
定义提供可用方程定理，因此验证入口采用与生成控制流同构的 total safe 实现。powmod、
subtract-X 以及 DDF 主循环均暴露了明确的递减度量或结构化实现，而不是把 partial 函数
作为公理使用。

## 验证

- `lake build CLPoly.Refinement.DDF CLPoly.Pipeline.L1`
- DDF 相关文件中不存在 `sorry`、`admit` 或 `native_decide`。
- `#print axioms` 审计 `ddf_correct`、`ddfZpSafe_refines` 和 `ddf_l1_correct`。

