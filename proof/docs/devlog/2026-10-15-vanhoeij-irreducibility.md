# 实际格路径与单元素分组的输出不可约性

日期：2026-10-15。

## 变更

- [VanHoeijIrreducibility.lean](../../lean/CLPoly/Refinement/VanHoeijIrreducibility.lean)：
  `hensel_lattice_validation_irreducible` 从真实 Hensel、CLD、高度、截断、
  构格、LLL、候选提取和完整验证执行推出乘积正确及逐输出不可约性。
  输出本原、非单位、好素数首项及模无平方性质均由验证和输入推出。
  对每个非单位除数，通过传递整除和本原正次数接入实际真因子支撑，
  实例化同组支撑不可分解判据，消去之前显式的除数支撑假设。
- [VanHoeijSingleton.lean](../../lean/CLPoly/Refinement/VanHoeijSingleton.lean)：
  单元素组的选中模乘积就是一个已知不可约模因子；验证后的模相伴、
  本原性和首项不消失推出整数不可约。未使用因子数量相等来证明一般组。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：先补充组合论证。

## 范围与后续

格路径组合定理构建通过（3574 jobs），单元素路径构建通过（3575 jobs）。
四条主定理的公理审计仅含 `propext`、`Classical.choice`、`Quot.sound`，
没有新增公理或 sorryAx。`git diff --check` 通过。
本机日志：`/tmp/clpoly-vh-irreducibility-build.log`、
`/tmp/clpoly-vh-singleton-build.log`、`/tmp/clpoly-vh-irreducibility-audit.log`。

格路径不以输出不可约、除数支撑存在或同组支撑一致性作为独立入口
假设；这些性质均从列出的实际执行前提导出。保留机器索引/次数范围、
LLL 输入合法性、截断参数/短向量界公式等明确前提，必须由严格顶层
控制器提供。单元素路径仍需控制器模型建立分组构造与参数条件。

本步骤不证明手写生成模板与 Clang AST 的语义等价，不证明尚未建模
的顶层重试控制器，也不把有限回退当成纯 van Hoeij 充分精度终止。
原生/参考对照/性能门槛仍待最终验收。没有 C++ 或生成模型改动，未
重复原生/B2B 测试。目标保持进行中，未合并。
