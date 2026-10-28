# 控制器单元素构造和列数增长的严格辅助模型

日期：2026-10-15。

## 实现

- [生产头文件](../../../clpoly/polynomial_factorize_univar.hh)：提取
  `__vanhoeij_singleton_groups` 和 `__vanhoeij_next_target`，生产路径调用
  这两个函数。保留原循环/表达式；单元素辅助入口明确断言非负计数。
- [生成模板](../../cpp2lean_v2/tests/build_strict_vanhoeij.py) 和
  [生成模型](../../lean/CLPoly/Generated/StrictVanHoeij.lean)：新增相应
  Int32 控制流。单元素循环检查数组访问并保留递增，证明 index<count
  时递增不溢出，用整数差的自然数截断作为终止量。
- [控制器算术证明](../../lean/CLPoly/Refinement/VanHoeijControllerArithmetic.lean)：
  证明 0<target<width 时下一轮列数严格增长、不超宽度，以及剩余
  宽度递减。明确处理有符号除法和乘法的 bmod 语义，不假设无溢出。
- [单元素构造证明](../../lean/CLPoly/Refinement/VanHoeijSingletonConstruction.lean)：
  成功的源循环只 set 外层数组条目，所以保持数组大小；入口成功后
  外层大小等于 count。这不替代仍待证明的每组内容与合法索引不变量。
- [原生注册表](../../b2b/registry/cpp_registry.cc)、
  [Lean 注册表](../../lean/B2B/Registry.lean) 和两个新增向量文件：
  单元素构造 51 例，列数增长 223 例，后者覆盖 INT_MAX 附近值。
  构造测试数量上限 4096 只是测试资源约束，不是形式模型中的限制。

## 边界

已完成验证：原生重组回归 382/382；生成模型构建 598 jobs；控制器算术
构建 599 jobs；数组大小不变量构建 600 jobs；模板 `--check` 及
`git diff --check` 通过。递增、构造循环定义及列数进展/递减量公理审计
仅含 `propext`、`Classical.choice`、`Quot.sound`。
数组大小不变量两条定理的公理审计也通过，仅依赖同一标准公理集合。

九组 B2B 已完成，750/750 通过。受模型改动影响的完整不可约性依赖重建
已沿原进程完成，3575 jobs 通过。对应本机日志：
`/tmp/clpoly-vh-controller-helpers-b2b.log`、
`/tmp/clpoly-vh-controller-helpers-refinement.log`。
原生结果见 `/tmp/clpoly-vh-controller-helpers-native.log`。

手写模板仍不构成 AST 翻译或 C++ 对应证明。当前只补齐控制器两个
实际辅助步骤；单元素数组内容及合法索引已在后续步骤证明，尚需与验证
入口组合。完整列循环、截断参数、精度
重试控制器及其入口前提仍须连接。纯 van Hoeij 充分精度终止未完成，
不能由这些局部递减量或三次尝试后的回退冒充。

生成模型变化需重建既有不可约性证明链；原生辅助函数提取后重跑
重组回归。目标保持进行中，未合并。
