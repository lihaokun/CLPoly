# 非空保留前缀、精确支持与实际组乘积

日期：2026-10-16。

## 变更与理由

在上一轮充分精度分组证明上，先补充
[自然语言证明](../../nl-proof/vanhoeij-eventual-separation.md) §14–16，
再新增五个模块、17 项定理，并加强既有实际支持分组入口。

- [VanHoeijRetainedNonempty](../../lean/CLPoly/Refinement/VanHoeijRetainedNonempty.lean)：
  将非零初始短向量经实际幺模 transform 换到约化基，取最后非零坐标。
  其 GS 范数上界与“全部行都被删除”的实际扫描条件矛盾。另从模 p
  的首项非零及次数保持证明真实正次数因子的物理支持非空。最后直接
  对整个 f 使用已有 Hensel 短向量证明，导出实际保留前缀非空，无需
  预先找到整数不可约因子。
- [VanHoeijExactGrouping](../../lean/CLPoly/Refinement/VanHoeijExactGrouping.lean)：
  核对来源索引经 Int32 转换的有界往返，证明实际输出组恰是某个列
  的完整相等列类。非空支持的短指示向量提供反向排除，结式分离提供
  正向包含，从而一个实际输出组恰好等于完整支持。
- [VanHoeijGroupMultiplicity](../../lean/CLPoly/Refinement/VanHoeijGroupMultiplicity.lean)：
  不仅证明成员相同，还按实际 collect 递归证明每行等于按来源顺序
  过滤出的索引列表。解码后列表无重复；该结论在任何恢复验证之前
  成立，避免借已验证因子的无平方性反推待验证组无重复的循环论证。
- [VanHoeijGroupSupportProduct](../../lean/CLPoly/Refinement/VanHoeijGroupSupportProduct.lean)：
  成员双向对应加两边 Nodup，给出实际解码组与合法 Hensel 支持的
  排列等价，随后证明 SelectedProductMod 与该支持乘积完全相同。
- [VanHoeijActualExactGrouping](../../lean/CLPoly/Refinement/VanHoeijActualExactGrouping.lean)：
  同一次 hensel_factor_truncated_short 同时提供支持、模乘积等式和
  短向量。将该支持直接代入充分精度分离，产出真实 groups 数组、
  真实提取运行等式、精确组成员和组的模乘积等式。没有分别选两个
  支持后假定它们相同，也不要求调用者假设分组成功或保留行非空。

既有 [VanHoeijActualSupportGrouping](../../lean/CLPoly/Refinement/VanHoeijActualSupportGrouping.lean)
已删除保留行非空前提，改用上述实际执行证明。必要的 CLD 高度值、
`m>height` 和 cutoff 下界显式保留；此前结式证明只给 cutoff 上界，
不能由此排除 cutoff 为零的情况。实际参数选择同时满足这两个方向
的界，后续控制器组合仍需显式传入相关已证明性质。

## 验证与边界

所有新模块及修改入口均已构建通过；包含实际组乘积的最终组合目标
为 3600 jobs，含依赖缓存。严格边界、refinement 注册和补丁格式检查
通过。合并公理审计 114/114 通过，3602 jobs；对完成的日志进行机器
解析，全部仅依赖 `propext`、`Classical.choice`、`Quot.sound`。
日志：`/tmp/clpoly-vh-exact-grouping-audit.log`。

本轮没有生产 C++、生成模板、B2B 注册或测试变更，未重复原生回归。

**剩余义务**：

1. 当前结论是“对每个被分析的真实不可约因子，能找到对应的精确
   实际输出组”。还需从完整输入分解与物理组覆盖推出所有组的对应，
   并在逐组提取过程中保持剩余因子的整除关系。
2. 实际 partitionProductLoop 的次数累加、对称恢复、primitiveRaw、
   精确除法和末尾剩余 1 的检查仍需执行完备性证明。组乘积等式是
   它们的输入证据，不能直接当作整个验证循环已经成功。
3. 当前 D 随 m 增长的截断策略与三次尝试后回退仍未改变；参数策略
   的阈值可达性、生成对应最终审计和性能验收仍未完成，尚未合并。

## 度量

- 耗时：未单独记录净证明时间；精确分组模块一次成功编译约 25 秒，
  组支持乘积模块约 2 秒，最终组合模块约 15 秒。
- 迭代：共 8 次目标构建调用，含最终合并审计。非空性模块首次通过；
  精确分组修正一次 lambda/依赖数组索引化简；来源列表模块修正两轮
  List API 隐式参数及 map 恒等式用法。一次检查脚本的 Python 路径
  拼写错误已修正后重跑通过，不算证明构建迭代。
- Lean：新增五个模块共 608 行、17 项定理、一份来源后缀定义；修改
  既有实际支持入口以消去非空前提，审计新增 17 项输出及一条 import。
- 对应 C++：已有实际 LLL 保留扫描、分组收集与模乘积，无原生逻辑变更。
- 关键决策：使用整个 f 的短向量证明非空；同一个支持同时服务于短向量
  与结式分离；从实际 collect 顺序证明无重复性。
- 放弃的方案：不从“已有验证成功”倒推组无重复，不以集合相等代替
  带重复次数的乘积相等，不引入输出不可约性或分组成功 oracle。
