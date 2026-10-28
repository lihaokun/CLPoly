# 严格 Hensel EEA 控制流

日期：2026-08-11

## 做了什么

- 用良基递归重建 Hensel 树构造调用的 Zp 扩展欧几里得主循环。
- 精确保留 quotient/remainder、`q*s`、`q*t`、减法、normalization、六个状态赋值以及最终三次首一缩放。
- 建立实际 raw 执行可达前缀、安全不变量和逐步语义 trace，并证明 raw→safe。

## 为什么做

旧 `Corpus` 的 `polynomial_GCD_eea` 最终落到 `partial def`，不能用于真正的 `__hensel_tree_build` 精化。该提交先移除 EEA 循环自身的部分递归。

## 关键决策

- 终止性由实际 EEA 状态度量给出；下降条件对 raw divmod 的每个成功结果量化，不是 fuel。
- 安全不变量对所有实际可达状态量化，不能提供预期 GCD、Bézout 系数或最终输出。
- quotient/remainder 仍保留为一个可执行 raw 边界；在公开树构造合同生成前，必须实例化为严格 C++ `pair_vec_div` 的双输出实现。当前没有把抽象边界冒充为最终闭合证明。

## 涉及文件

- `proof/lean/CLPoly/Generated/StrictHensel.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
