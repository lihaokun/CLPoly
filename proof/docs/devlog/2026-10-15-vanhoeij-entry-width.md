# 宽度对应与完整尝试入口模型

日期：2026-10-15。

## 变更

- [宽度证明](../../lean/CLPoly/Refinement/VanHoeijWidth.lean)：证明 Int32 的 min/max 投影、首项次数转换及加一、实际宽度扫描与构格自然数宽度一致，以及初始目标列数范围。
- [入口模板](../../cpp2lean_v2/tests/build_strict_vanhoeij_entry.py) 与 [生成入口](../../lean/CLPoly/Generated/StrictVanHoeijEntry.lean)：按源顺序连接断言、高度计算、单元素构造与验证、参数选择、CLD、截断、宽度和列循环。
- [具体入口执行](../../lean/CLPoly/Refinement/VanHoeijEntryExecution.lean)：从已有实际运行定理提供 CLD/截断大小、参数模数非零与初始目标范围。擦除的证明字段不选择数据或返回结果。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：先记录扫描不变量和完整入口组合次序。

## 条件与剩余工作

宽度对应仍显式要求每个非空 CLD 首项次数加一小于 2^31，需从实际
CLD/截断运行与输入次数导出。入口保留 lifted.size<2^31 的擦除范围
前提；它不是额外运行时检查，也尚未从 Hensel 入口导出。

这是完整 try-recombine 尝试入口的执行模型，不是其完整正确性定理。
下一步需从实际入口成功运行追踪单元素或列循环成功分支，再接入两条
已有不可约性证明。外层精度重试、全入口 B2B、C++ 对应、充分精度
最终分离、参考对照与性能验收仍未完成。模板不是 AST 翻译。
后续入口分支追踪与条件不可约性连接已完成，见
[入口不可约性日志](2026-10-15-vanhoeij-entry-irreducibility.md)。

## 验证

宽度证明构建通过（601 jobs）；入口模型与具体执行实例构建通过
（3544 jobs）。日志为 `/tmp/clpoly-vh-width-build.log` 和
`/tmp/clpoly-vh-entry-execution-build.log`。入口模板 `--check`、
`git diff --check` 通过。本轮未改生产代码，未重复上一轮原生/B2B。
七条宽度/范围定理及具体入口执行实例公理审计通过，仅依赖
`propext`、`Classical.choice`、`Quot.sound` 的子集，无新增公理或
sorryAx；日志为 `/tmp/clpoly-vh-entry-width-audit.log`。
