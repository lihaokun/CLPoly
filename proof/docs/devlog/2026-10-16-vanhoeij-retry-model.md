# 外层三次精度重试模型

日期：2026-10-16。

## 变更

- [重试模板](../../cpp2lean_v2/tests/build_strict_vanhoeij_retry.py) 与 [生成模型](../../lean/CLPoly/Generated/StrictVanHoeijRetry.lean)：复用原启发式精度计算，每轮执行新的严格 Hensel 提升与新 try-recombine；成功立即返回，不比较因子数。前两次失败经溢出检查计算 max(full,2*precision)，第三次失败直接结束循环，再对最后一次提升执行严格 Zassenhaus。
- [精度算术证明](../../lean/CLPoly/Refinement/VanHoeijRetryArithmetic.lean)：非负精度下成功更新的 Int32 投影精确等于整数 max，并至少为 full 与旧 precision；超过 INT_MAX/2 时返回明确溢出 fault。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：先记录控制流顺序、attempt 的严格递减量与精度转换条件。

## 决策及边界

重组调用固定为新尝试入口和实际 Zassenhaus 模型，不通过可替换因子
输出的回调。重试执行参数仍包含既有 Hensel 步骤/建树算术和证明义务：
实际提升输出计数符合 Int32 范围。具体实例及整体正确性尚未完成，不能
将参数结构的存在作为该义务已解决的证据。

循环的 `(3-attempt.toInt).toNat` 严格递减，只证明有限重试控制器的
终止。未证明纯 van Hoeij 在充分精度最终分离，也不以单次 Zassenhaus
回退替代该证明。重试 B2B、上游范围推导、参考对照及性能验收仍待完成。

## 验证

生成模型构建通过（603 jobs），精度算术构建通过（606 jobs）。日志
为 `/tmp/clpoly-vh-retry-build.log` 和
`/tmp/clpoly-vh-retry-arithmetic-build.log`。模板 `--check` 与
`git diff --check` 通过。本步未修改生产代码。
重试循环与两条算术定理公理审计通过，仅依赖 `propext`、
`Classical.choice`、`Quot.sound` 的子集，无新增公理或 sorryAx；
日志为 `/tmp/clpoly-vh-retry-arithmetic-audit.log`。
