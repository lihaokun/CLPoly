# 严格 Hensel 显式精度目标

日期：2026-08-11

## 做了什么

- 严格翻译 `__hensel_lift` 中 `a_target > 0` 的循环，计算 `p^a_target - 1`。
- 证明循环结果精确等于整数幂，从而证明 C++ 停止条件 `m <= target` 与目标精度一致。

## 为什么做

完整入口必须覆盖显式精度分支，不能把调用者给定的 `a_target` 替换成 L2 直接幂运算。

## 关键决策

- 结构递归参数是源 `for` 循环尚需执行的真实次数，不是人为 fuel。
- 非正输入保留为 raw assertion fault；完整入口只会在 `a_target > 0` 分支调用它。
- 未复用旧 `Corpus` 中含 `partial def` 的树构造/EEA 路径。

## 涉及文件

- `proof/lean/CLPoly/Generated/StrictHensel.lean`
- `proof/lean/CLPoly/Refinement/Hensel.lean`
