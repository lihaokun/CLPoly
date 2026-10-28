# 实际单元素快速路径：成功构造与不可约性

日期：2026-10-15。

## 变更与理由

- [构造证明](../../lean/CLPoly/Refinement/VanHoeijSingletonConstruction.lean)：新增循环和入口成功性定理；非负计数下，数组大小与有界递增不变量排除断言失败及越界，得到实际成功返回。
- [组合入口](../../lean/CLPoly/Refinement/VanHoeijSingletonExecution.lean)：从真实选素数、Hensel 输出、计数对应、单元素构造和完整验证成功，推出乘积等于输入且所有返回因子不可约。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：先记录成功性归纳和组合入口推理。

组形状、合法索引由实际构造推出；模数正性、素数整除、模因子不可约
由 Hensel 与选素数正确性推出。仍保留输入本原、规范表示、机器次数
范围及真实运行等式等明确前提。没有新增输出不可约性假设。

## 验证

构造总性构建通过（600 jobs），组合入口构建通过（3579 jobs）。
日志分别为 `/tmp/clpoly-vh-singleton-total-build.log` 和
`/tmp/clpoly-vh-singleton-execution-build.log`。`git diff --check` 通过。
三条新增定理公理审计通过，仅依赖 `propext`、`Classical.choice`、
`Quot.sound`，无新增公理或 sorryAx；审计日志为
`/tmp/clpoly-vh-singleton-execution-audit.log`。
本轮仅改证明与文档，未重复原生测试；最近生产辅助函数修改后的结果
为原生 382/382、九组 B2B 750/750。

## 剩余义务

接下来在独立控制器模型中连接参数选择、列循环与精度重试，避免继续
修改稳定辅助模块引发全链条重建。手写控制流模板加 B2B 不等同于 AST
翻译或 C++ 等价性证明。充分精度下最终分离、参考实现对照及 r=24/32
性能验收仍未完成；三次重试后回退不能代替前一项，尚未合并。
