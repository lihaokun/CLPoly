# 回退指数来源和实际 Hensel 恢复界

日期：2026-10-16。

## 证明链

先更新[草稿](../../nl-proof/vanhoeij-certificate-repair.md)，随后证明：

- [精度来源](../../lean/CLPoly/Refinement/VanHoeijRetryPrecision.lean)：成功的
  活动重试轮必具有非负请求指数；维持 attempt≥1 时 precision≥full。
  前两次失败通过已证 nextPrecision 算术建立不变量；第三次保持。
  初始 attempt=0 的前态条件真空成立，且成功顶层保留实际启发式运行、
  起始指数非负和加强后的返回来源。
- [恢复界](../../lean/CLPoly/Refinement/VanHoeijRetryRecovery.lean)：假设实际
  Hensel 语义证书、非负 full 以及 target≥full。target=0 使用真实
  Mignotte 目标分支；target>0 由输出模数>p^target−1 得到模数≥p^target。
  实际启发式给出 p^full 超过恢复目标，指数单调性完成连接。

这里复用旧模块的启发式幂循环、Hensel 模数界和整数恢复引理；不使用
旧控制器或其因子数分支的正确性定理。起始指数非负结合启发式 min
性质可推出 full 非负，仍需在新的端到端定理中完成最终组装。

## 验证

- 精度来源：3340 jobs 构建通过。
- 恢复界及[扩展审计](../../lean/CLPoly/Refinement/VanHoeijRetryAudit.lean)：
  3567 jobs 构建通过；13 项仅出现 `propext`、`Classical.choice`、
  `Quot.sound` 的子集，无 sorryAx。
- 日志 `/tmp/clpoly-vh-retry-precision-build.log`、
  `/tmp/clpoly-vh-retry-recovery-build.log`。
- 本步未改变任何运行函数，沿用此前 27/27 完整重试 B2B。
- 尚未声明整个重试入口无额外前提的不可约性，也未解决纯 van Hoeij
  的充分精度最终分离或机器范围的全部上游推导。

## 度量

- 耗时：约 0.2 小时（草稿、形式化与编译检查）。
- 迭代：精度来源 2 轮构建，恢复界及审计 1 轮构建。
- Lean：新增两个证明模块，修改公理审计；生产 C++ 修改 0 行。
- 问题与处理：Int32 比较的整数投影需先重写比较定理再规范化字面量；
  cases 已重写的运行等式用 rfl 提供，避免重复使用原式。
- 放弃方案：无。
