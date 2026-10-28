# 实际单元素分组构造的内容与合法索引

日期：2026-10-15。

## 变更

- [VanHoeijSingletonConstruction.lean](../../lean/CLPoly/Refinement/VanHoeijSingletonConstruction.lean)：
  对实际 Int32 循环证明前缀内容不变量。每步把当前位置的空组变为
  `[index]`，机器整数往返保证与 `[Int32.ofNat position]` 一致；其余
  位置不变。终止时所有合法位置都属于已完成前缀。
- [VanHoeijSingletonConditions.lean](../../lean/CLPoly/Refinement/VanHoeijSingletonConditions.lean)：
  由内容定理导出每组大小为 1；由外层大小和 Int32 范围导出索引非负、
  不越界及转换往返，提供快速验证入口所需条件。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：先记录不变量与推导。

## 问题及处理

展开源循环时，已有非负前提使 Lean 自动消去非负分支，归纳分支必须
相应调整；数组更新需在转换索引前使用当前空组等式。修正过程没有
更改模型、生产代码或测试向量。

## 已确认验证与边界

上一轮 B2B 已返回终态，九组 750/750 全部通过。原生回归仍为上一轮
重跑的 382/382；本轮仅改证明，没有重复原生测试。内容、单元素形状及
合法索引不变量的定向构建已通过（3535 jobs），`git diff --check` 通过。
完整不可约性依赖重建已沿上一轮进程完成，3575 jobs 全部通过。
四条新增内容/形状/索引定理的公理审计已完成，仅含 `propext`、
`Classical.choice`、`Quot.sound`，无新增公理或 sorryAx。日志为
`/tmp/clpoly-vh-singleton-contents-audit.log`；本步构建结果见
`/tmp/clpoly-vh-singleton-conditions-build.log`。

后续控制器模型应独立于稳定的辅助模型模块，减少全链条无必要重建。
还需把实际单元素构造条件接入验证入口，以及完成列循环、参数选择和
精度重试模型。充分精度终止、参考对照及性能验收仍待完成，未合并。
