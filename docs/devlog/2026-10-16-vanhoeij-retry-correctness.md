# 新重试入口的条件式端到端正确性

日期：2026-10-16。

## 完成内容

证明实际 Hensel 输出因子的次数不超过输入多项式次数。输出为首一，
模好素数后次数保持；其模多项式整除输入的模像，由此得到次数界。
规范稀疏首项的次数等于数学次数，因此从输入 degree<2^31 推出输出
物理首项次数<2^63，删除单次尝试定理中的独立输出次数前提。

新增真实 Hensel 连接和重试入口正确性定理：从成功运行取得计数检查、
指数检查与输出数组，结合实际素数选择及底层执行就绪，建立同一数组
的 Hensel 语义证书。完成验证分支使用新 van Hoeij 证书，回退分支
使用已推导的指数下界与恢复界，最终返回数组满足乘积还原及逐因子
不可约性。没有使用旧控制器的等基数判据，也没有添加输出不可约性假设。

## 验证与边界

新增次数模块构建通过（3547 jobs）；新的完整重试正确性构建通过
（3602 jobs），日志 `/tmp/clpoly-vh-retry-correctness-build.log`。
扩展公理审计见[证明日志](../../proof/docs/devlog/2026-10-16-vanhoeij-retry-correctness.md)。
本轮只修改证明与文档，沿用上一轮完整重试 27/27 B2B，不重复运行未变代码。

这是带明确输入与执行就绪条件的正确性定理。它仍要求素数选择正确性、
稠密运算配置和递归 Hensel 执行就绪；尚未迁移外层 squarefree FactorZZ
调用模型，也未证明纯 van Hoeij 充分精度最终分离。参考对照、性能验收
及完整原生构建问题仍待完成，研究分支未合并。

## 文件

- [Hensel 次数界](../../proof/lean/CLPoly/Refinement/VanHoeijHenselDegree.lean)
- [实际 Hensel 连接](../../proof/lean/CLPoly/Refinement/VanHoeijRetryHensel.lean)
- [新重试正确性](../../proof/lean/CLPoly/Refinement/VanHoeijRetryCorrectness.lean)
- [加强后的单次尝试](../../proof/lean/CLPoly/Refinement/VanHoeijEntryIrreducibility.lean)
- [公理审计](../../proof/lean/CLPoly/Refinement/VanHoeijRetryAudit.lean)
- [自然语言草稿](../../proof/nl-proof/vanhoeij-certificate-repair.md)
