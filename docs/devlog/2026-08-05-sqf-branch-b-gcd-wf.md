# SQF Branch B：GCD 良基递归前置

日期：2026-08-05

## 做了什么

- 将 `SparsePolyZp.gcdAux` 从 `partial def` 改为以除数首项次数为度量的良基递归定义。
- 在 Euclid 递归处显式检查余式为空或首项次数严格下降；合法规范形输入保持原有递归路径，非规范输入保守终止。
- 删除 `Refinement/SquarefreeZp.lean` 中未使用且带 `admit` 的 `get_deg_toPoly`。

## 为什么做

SQF 导数非零分支需要证明 C++ Yun 循环中的 `polynomial_GCD` 与 L2 多项式 GCD 对应。原 `partial def gcdAux` 没有可用归纳原理，无法建立 Euclid 循环的精化证明。当前 `divmodAux_snd_deg_lt` 已经给出合法除法余式的严格次数界，因此可以把 GCD 总化为真正的良基递归。

## 关键决策

- 度量使用：空数组为 `0`，非空数组为首项次数加一。
- 余式为空时直接返回当前除数，这与原定义下一次递归命中 base case 的结果一致。
- 若表示不变量不成立、`divmod` 回退且次数没有下降，则返回当前除数，避免原 partial 模型在异常输入上潜在不终止；合法输入不经过该分支。
- 未使用的 `get_deg_toPoly` 不属于主证明依赖，直接删除而不是为死代码引入数组排序前提。

## 验证

- `lake env lean CLPoly/Model.lean`：通过，GCD 示例恢复可执行求值。
- 全量 `lake build`：见本次提交前验证结果。

## 涉及文件

- `proof/lean/CLPoly/Model.lean`
- `proof/lean/CLPoly/Refinement/SquarefreeZp.lean`

## 度量

- 耗时：约 0.5 小时（现状核对、设计复核、形式化与调试）
- 迭代：3 轮 `lake env lean` 编译-修复
- Lean 新增/修改行数：约 35 行
- 对应 C++ 行数：约 5 行 Euclid GCD 控制流
- 放弃的方案：恢复旧的影子 `gcdAux'_wf`；该方案仍需证明影子定义等于实际 partial 定义，不能直接闭合真实 L1 精化链。
