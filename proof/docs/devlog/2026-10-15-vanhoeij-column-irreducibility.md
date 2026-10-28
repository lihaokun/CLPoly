# 实际列循环到不可约性

日期：2026-10-15。

## 变更

- [列数范围证明](../../lean/CLPoly/Refinement/VanHoeijColumnBounds.lean)：构格循环保持 added≤target，实际返回列数不超过 target；Int32 请求值对应的返回列数小于 2^31，转换回 Int32 的整数值不变。
- [列循环组合](../../lean/CLPoly/Refinement/VanHoeijColumnIrreducibility.lean)：从初始未完成、最终完成的实际列循环取出成功轮次，再追踪该轮初始化、构格、LLL、提取和完整验证。实际构格证明提供 LLL 合法性，数组大小及非负计数消去计数转换，列数界消去 added 转换，然后调用已有格验证不可约性定理。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：先记录范围不变量、转换无溢出及组合推理。

## 边界

组合定理保留真实选素数/Hensel 运行、CLD/截断运行、源表示与次数范围、
参数条件及源计数/数组大小对应。它不是全入口控制器证明，参数选择、
宽度扫描、精度重试及 C++ 对应仍需补齐。充分精度分离与性能验收不因
此定理而完成。没有新增输出不可约性假设或独立 LLL 输入合法性假设。

## 实现及验证

构格列数证明已构建通过。Lean 的 `Int32.toInt_ofNat` 处理数值字面量，
不能用于这里的运行时 `Int32.ofNat`；改用已有范围引理
`Int32.toInt_ofNat_of_lt`。未扩大递归限制，也未改变机器转换模型。
列循环组合入口定向构建通过（3581 jobs），日志为
`/tmp/clpoly-vh-column-irred-build.log`；`git diff --check` 通过。
四条新增定理公理审计通过：列数界只依赖 `propext`、`Quot.sound`；
机器转换及组合定理另依赖 `Classical.choice`。无新增公理或 sorryAx，
日志为 `/tmp/clpoly-vh-column-irred-audit.log`。控制器模板 `--check` 通过。
本轮未修改生产代码，未重复原生/B2B 测试。
