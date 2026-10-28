# SQF 生成 Yun raw 转移

日期：2026-08-10

## 做了什么

- 从实际 raw GCD 和首次精确除法返回值构造生成 Yun 初始状态。
- 完成 Yun 公共 raw 单步准备定理，证明实际 GCD 与两次精确除法输出的
  canonical、monic、normalization、长度界和严格 measure 下降。
- 分别完成生成 Yun 的有因子与无因子转移，包括 make-monic 输出对齐、结果
  push 的 canonical 保持、multiplicity 加一无回绕及预算保持。

## 为什么做

这些转移是实例化 `Generated.StrictSquarefreeZp.SQFRawOps` 的执行核心；没有
它们，生成 Yun 循环仍无法绑定到真实 C++ raw 操作。

## 关键决策

- 所有中间多项式通过 `Except.ok` 单射与 raw 执行结果对齐。
- measure 下降来自 GCD 的整除关系和两次真实精确除法，不是额外的 fuel 或
  调用者提供的结果预言。

## 验证

- `lake env lean CLPoly/Refinement/StrictSquarefreeGenerated.lean`：通过。

## 涉及文件

- `proof/lean/CLPoly/Refinement/StrictSquarefreeGenerated.lean`
- `docs/devlog/2026-08-10-sqf-generated-yun-raw-transitions.md`
