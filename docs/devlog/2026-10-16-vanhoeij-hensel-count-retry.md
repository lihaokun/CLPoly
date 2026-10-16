# Hensel 计数与完整重试入口

日期：2026-10-16。

## 变更与动机

证明实际 Hensel 成功输出的因子数等于输入模因子数，将入口已有的
Int32 计数检查传递到 van Hoeij 调用。用该证明构造具体重试执行实例，
消去原先悬置的输出计数条件。随后证明完整重试成功结果必来自某次
实际完成的 van Hoeij 验证，或某次实际 Hensel 输出上的 Zassenhaus。

新增完整重试 B2B 注册和 27 个独立预期用例，覆盖模因子数大于真因子数、
非首一多项式、单元素快速路径、不同模因子顺序和 16/48/96 位高度扰动。
正二次多项式及第八分圆多项式提供已知不可约因子；生成器独立检查模根
简单性和模乘积重构。两端结果均与已知因子比较，不只互相比对。

## 验证与边界

- 具体重试实例构建通过（3548 jobs）；来源证明和 B2B 注册构建通过（3566 jobs）。
- 完整重试 B2B：27/27 通过，日志 `/tmp/clpoly-vh-retry-b2b.log`。
- 九个新证明/实例的公理审计仅出现标准 `propext`、`Classical.choice`、
  `Quot.sound` 的子集，无新增公理或 `sorryAx`。
- 本步未修改生产算法；沿用此前生产版本原生 382/382 通过的结果。
- 来源证明尚未推出最后提升达到恢复界，也未连接为整个重试入口的不可约性定理。
- 27 个用例不声称穷尽三次失败回退分支。纯 van Hoeij 充分精度最终分离、
  上游范围条件、参考对照与性能验收仍未完成，不合并分支。
- 严格模型仍为手写控制流模板，B2B 通过不等于 AST 翻译或普遍等价性证明。

## 文件

- [计数证明](../../proof/lean/CLPoly/Refinement/VanHoeijHenselCount.lean)
- [具体重试实例](../../proof/lean/CLPoly/Refinement/VanHoeijRetryExecution.lean)
- [成功来源](../../proof/lean/CLPoly/Refinement/VanHoeijRetryTrace.lean)
- [公理审计](../../proof/lean/CLPoly/Refinement/VanHoeijRetryAudit.lean)
- [自然语言草稿](../../proof/nl-proof/vanhoeij-certificate-repair.md)
- [测试生成器](../../proof/b2b/generate_vanhoeij_retry_vectors.py)
- [测试向量](../../proof/b2b/vectors/__lll_factorize.json)
- [C++ 注册](../../proof/b2b/registry/cpp_registry.cc)、[Lean 注册](../../proof/lean/B2B/Registry.lean)
