# 真实模原子分离与验证结果物理支撑

日期：2026-10-15。

## 变更与理由

- [VanHoeijHenselAtoms.lean](../../lean/CLPoly/Refinement/VanHoeijHenselAtoms.lean)：
  由选素数/Hensel 来源链的模不可约性推出模原子素性；实际调整因子
  互素性传递到提升数组，排除不同物理索引相伴。
- [VanHoeijSupportMultiplicity.lean](../../lean/CLPoly/Refinement/VanHoeijSupportMultiplicity.lean)：
  从非单位原子乘积无平方推出支撑索引无重复；由完整验证乘积等式
  推出返回因子整除输入，继而保持好素数下的无平方性。
- [VanHoeijPhysicalSupport.lean](../../lean/CLPoly/Refinement/VanHoeijPhysicalSupport.lean)：
  实际分组按原位置转换为 `Fin` 列表，证明与实际模乘积相等。
  将上述事实组合为成功验证因子的实际组来源、合法索引、物理支撑
  无重复及模相伴关系。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：补充上述论证。

## 决策与边界

定向构建最终通过（3568 jobs）。七条主定理执行公理审计，仅使用
`propext`、`Classical.choice`、`Quot.sound` 的子集，无新增公理或
sorryAx。日志为本机 `/tmp/clpoly-vh-physical-support-build.log` 和
`/tmp/clpoly-vh-physical-support-audit.log`。`git diff --check` 通过。

物理列表乘积证明中，函数复合规约后的数组索引需要用 `change` 显式
展开，再应用已知合法索引的 `getElem!` 等式；未修改模型规避类型问题。

无重复性通过验证输出与好素数证明，不另外增加候选收集的状态不变量。
有界索引转换保留顺序和重复；不能把提前去重后的列表冒充运行时组。
仍需将真因子支撑的成员一致性转换到这些物理列表，并完成非单位补
因子的矛盾，才能声称整数输出不可约。

本步没有修改 C++ 或生成模型，未重复原生/B2B 测试。严格顶层控制器、
纯 van Hoeij 充分精度终止、参考对照及性能验收仍待完成，未合并。
