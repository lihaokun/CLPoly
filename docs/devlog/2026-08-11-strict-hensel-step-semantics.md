# 严格 Hensel 单步生成语义

## 日期

2026-08-11

## 做了什么

- 审计旧 Hensel generated corpus，确认其不能作为精化证明的 L1。
- 新增 `build_strict_hensel.py`，生成 C++ `__hensel_step` 的严格 raw
  控制流 `Generated.StrictHensel.__hensel_step_raw_ir`。
- 显式保留 C++ 中两次“精确除法后模约化”、四次系数乘 `m` 循环、两次
  modular divmod，以及因子和 Bézout 系数更新的原始顺序。
- 证明首个 raw→L2 桥 `scaleCoeffs_toPoly`。

## 为什么做

旧 `Corpus.lean` 中的 Hensel 函数是 `partial def`；更严重的是，C++
`__hensel_node` 聚合初始化经由 lossy coercion 直接返回默认节点，丢弃全部
字段。继续基于它证明会退化为错误语义上的“精化”。严格版本必须先重建
可信、可终止且逐操作失败透明的 L1。

## 关键决策

- 所有尚未闭合的底层多项式运算通过 `HenselStepRawOps` 暴露为
  `RawExec`；任一调用失败就原样失败，不允许 L2 fallback。
- 系数遍历直接由 `Array.map`/`filterMap` 总函数表达，不使用 fuel。
- 暂不发布最终 Hensel 定理；只有底层运算律、单步不变式、树递归和入口
  全部闭合后才进入集中 `Refinement/Generated.lean`。

## 验证

- `build_strict_hensel.py --check` 能验证生成文件同步。
- `CLPoly/Generated/StrictHensel.lean` 单文件编译成功。
- `CLPoly/Refinement/Hensel.lean` 单文件编译成功。

## 涉及文件

- `proof/cpp2lean_v2/tests/build_strict_hensel.py`
- `proof/lean/CLPoly/Generated/StrictHensel.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
