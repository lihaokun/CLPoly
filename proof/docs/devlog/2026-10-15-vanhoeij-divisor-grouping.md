# 同一真因子支撑的模投影、物理转换与格组合

日期：2026-10-15。

## 完成内容

- [VanHoeijDivisorSupport.lean](../../lean/CLPoly/Refinement/VanHoeijDivisorSupport.lean)：
  合法自然数数组逐位置构造 Fin 列表；证明 val、成员资格及乘积保持。
  将原模数上的真因子等式投影到所选素数，消去非零首项单位，得到
  同一物理支撑的模相伴关系。
- [VanHoeijDivisorGrouping.lean](../../lean/CLPoly/Refinement/VanHoeijDivisorGrouping.lean)：
  接入真实 Hensel、高度、CLD、截断、初始格、LLL 及候选提取执行，
  对每个正次数整数因子产出同一个物理支撑的模相伴和所有类内的一致性。
- [VanHoeijValidationNonunit.lean](../../lean/CLPoly/Refinement/VanHoeijValidationNonunit.lean)：
  完整验证循环的非单位输出不变量。空组累计次数零，与成功分支矛盾；
  接受的非空组给出非单位模乘积，因此恢复的整数因子非单位。另补
  本原非单位多项式正次数、整数除数保持好素数首项条件两个接口引理。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：先补充相应论证。

## 关键决策与边界

定向构建通过：真因子支撑组合 3572 jobs，验证非单位性及两条辅助
定理 3563 jobs。八条主定理公理审计仅包含 `propext`、`Classical.choice`、
`Quot.sound`；无新增公理或 sorryAx。`git diff --check` 通过。
构建/审计日志分别保存在本机 `/tmp/clpoly-vh-divisor-grouping-build.log`、
`/tmp/clpoly-vh-validation-nonunit-build.log`、
`/tmp/clpoly-vh-divisor-grouping-audit.log`。

模数投影和物理转换保留原始真因子支撑见证，不能重新选择模 p 支撑，
否则格证明的成员一致性无法直接适用。输出非单位性来自实际次数检查
及恢复运行，而非调用方假设。

实际执行支撑定理仍保留明确机器次数/索引范围及 LLL 输入合法性等
入口条件。尚待将其与完整验证输出、本原性和不可约性判据组合。
严格顶层控制器、充分精度终止、参考对照及性能验收仍未完成。

本步没有修改 C++ 或生成模型，未重复运行原生/B2B 测试。未合并，
目标继续进行。
