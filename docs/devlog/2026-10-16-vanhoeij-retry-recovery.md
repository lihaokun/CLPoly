# 重试回退的指数与恢复界

日期：2026-10-16。

## 变更

继续补齐[重试来源](2026-10-16-vanhoeij-hensel-count-retry.md)中尚未证明
的回退指数条件。实际 Hensel 成功排除负请求指数；前两次失败后的精度
更新至少达到 full，第三次失败保持该不变量。因此完整重试的回退来源
来自某次请求指数不小于 full 的实际提升。

新增恢复界引理允许请求指数大于启发式 full，不要求二者相等：启发式
幂循环越过 Mignotte 恢复目标，实际 Hensel 模数至少为请求的素数幂，
指数单调性连接二者，得到 LiveRecoveryPrecision。零请求指数单独使用
Hensel 的 Mignotte 目标分支。

## 验证及限制

精度来源构建通过（3340 jobs），恢复界和扩展公理审计构建通过
（3567 jobs）。13 项审计仅包含标准公理子集，无新增公理或 sorryAx。
日志 `/tmp/clpoly-vh-retry-recovery-build.log`。本步未改生产代码或测试
运行路径，未重复运行刚通过的 27/27 完整重试 B2B。

两部分已分别证明，尚需在新的端到端重试定理中组装：从实际 Hensel
执行导出语义证书与次数范围，分别连接完成验证和 Zassenhaus 分支。
纯 van Hoeij 充分精度分离、参考对照及性能验收仍未完成。

## 文件

- [指数来源](../../proof/lean/CLPoly/Refinement/VanHoeijRetryPrecision.lean)
- [恢复界](../../proof/lean/CLPoly/Refinement/VanHoeijRetryRecovery.lean)
- [审计](../../proof/lean/CLPoly/Refinement/VanHoeijRetryAudit.lean)
- [自然语言草稿](../../proof/nl-proof/vanhoeij-certificate-repair.md)
