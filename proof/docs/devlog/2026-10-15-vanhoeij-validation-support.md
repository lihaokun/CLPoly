# 实际验证循环的模支撑与物理组来源

日期：2026-10-15。

## 变更与动机

验证乘积正确性不足以推出输出不可约；还必须把输出与 LLL 产生的
实际分组连接。本步骤补齐这种对应，而非假设重构因子具有预期支撑。

- [VanHoeijValidationSupport.lean](../../lean/CLPoly/Refinement/VanHoeijValidationSupport.lean)：
  从成功验证的乘积等式推出输出本原；将新增次数累积的分组乘积循环
  投影到已有乘积执行；证明任意整除实际模数的基模上的缩放支撑等式；
  证明所选素数上 content 可逆的条件及模支撑相伴；精确除法直接保持
  恢复因子与未再本原化余商的首项系数模 p 非零。
- [VanHoeijValidationProvenance.lean](../../lean/CLPoly/Refinement/VanHoeijValidationProvenance.lean)：
  对实际循环归纳，跟踪余商规范性、首项系数不变量与已输出因子来源。
  最后通过排序置换保持来源，推出成功入口的每个返回因子都与某个
  输入物理组的选中模乘积相伴。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：先补充相应论证。

## 决策与问题

复用旧模块仅限乘积及恢复子程序的执行语义；不使用旧顶层控制器
来声称新控制器正确。新余商没有再次 primitive，因此直接从零余数
精确除法的乘积等式证明首项系数不变量，未插入不存在的运行步骤。

Int32 比较需要使用其 toInt 比较等价式；一般线性序推理不能直接
应用。Int.cast 的函数应用也需先规约后才能重写乘法。

## 验证与边界

定向 `lake build CLPoly.Refinement.VanHoeijValidationProvenance` 通过
（3563 jobs）。`git diff --check` 通过。没有 C++ 或生成模型改动，
故本步不重复原生/B2B 回归；此前结果不等同于本步新增测试。

对八条主要新定理执行 `#print axioms`，结果仅包含 `propext`、
`Classical.choice`、`Quot.sound`（执行投影不依赖 choice）；无新增
公理或 sorryAx。审计输出保存在本机临时日志
`/tmp/clpoly-vh-validation-support-audit.log`。

入口仍明确要求各组合法索引、提升因子模 p 不可约、输入规范及好素数
首项系数条件。这里不可约性针对模因子，不是假设整数输出不可约。
尚需从实际候选提取/Hensel 输入导出这些条件，并将支撑相伴与真因子
同组支撑一致性组合为最终不可约性。严格控制器、纯 van Hoeij 充分
精度终止、参考对照及性能门槛仍待完成。未合并，目标继续进行。
