# 实际 Hensel 计数保持与重试成功来源

日期：2026-10-16。

## 证明

先在[自然语言草稿](../../nl-proof/vanhoeij-certificate-repair.md)记录推理，
再形式化以下步骤。

1. 非空区间拓扑的内部节点数加一等于区间长度。中点严格位于区间内，
   两个子区间归纳等式相加得到父区间等式。
2. 任意实际叶因子提取成功时，输出大小等于输入大小加节点数加一。
   四种子树组合分别处理 push 和递归运行。
3. 最终规范化返回原数组或 set! 第零项，保持外层大小。
4. 原始 Hensel 完整入口成功时，从空数组提取再规范化，输出因子数恰好
   等于输入模因子数；原入口计数检查传递到输出，构造具体重试实例。
5. 重试变式强归纳：每次递归记录刚完成的 Hensel 运行。成功分支保留
   实际 try-recombine 的 complete 和因子数组等式；退出分支保留实际
   Zassenhaus 运行。初始 attempt=0 无需假设已有提升来源。

## 实现问题与决策

只使用已有内部节点数上界已经足够满足 Int32 范围，但进一步证明精确
计数以匹配草稿且防止掩盖丢失叶因子的错误。局部 let 需显式展开才可
使用 omega；拓扑归纳使用现有左右子树严格递减引理，避免将子节点计数
无谓展开成递归表达式。两个命名空间同时存在 nextPrecision，引用时
限定新重试命名空间。

## 验证及范围

- [计数](../../lean/CLPoly/Refinement/VanHoeijHenselCount.lean)和
  [实例](../../lean/CLPoly/Refinement/VanHoeijRetryExecution.lean)：3548 jobs 构建通过。
- [来源](../../lean/CLPoly/Refinement/VanHoeijRetryTrace.lean)和新 B2B 注册：3566 jobs 构建通过。
- [审计](../../lean/CLPoly/Refinement/VanHoeijRetryAudit.lean)：九项均仅依赖标准公理子集，
  日志 `/tmp/clpoly-vh-retry-audit.log`。
- [生成器](../../b2b/generate_vanhoeij_retry_vectors.py)生成完整重试输入，
  27/27 B2B 通过，日志 `/tmp/clpoly-vh-retry-b2b.log`。
- 测试运行时沿用既有物理堆的擦除证明设施；新输出计数条件来自已证定理。
  测试运行时不被形式化证明模块导入。
- 本步未证明回退精度足够、完整重试输出不可约或纯 van Hoeij 最终分离。
  也未主张 27 个用例覆盖全部回退路径、模板等于 AST 自动翻译或普遍实现等价。

## 度量

- 耗时：约 0.3 小时（本轮续接，含计数修正、来源证明和 B2B）。
- 迭代：本轮 6 次相关构建，修复计数展开、存在量词及命名空间歧义。
- Lean 新增：276 行（计数、实例、来源、审计），另修改 B2B 注册。
- 对应 C++：本步生产代码修改 0 行；验证既有 Hensel 提取/规范化及重试控制流。
- 放弃方案：仅用计数上界满足实例，改为精确等式；未改变运行算法。
