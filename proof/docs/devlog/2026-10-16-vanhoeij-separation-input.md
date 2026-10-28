# 从实际 Hensel/CLD 执行导出分离前提

日期：2026-10-16。

## 变更与理由

此前完整残差证明仍需调用方给出 CLD 次数、局部因子整除和互素条件。
本轮先补[自然语言证明](../../nl-proof/vanhoeij-eventual-separation.md) §9–10，
再用现有真实执行证书导出这些条件，新增十项定理及一份仅含 Prop 的
分离输入证书。

- [VanHoeijCldDegree](../../lean/CLPoly/Refinement/VanHoeijCldDegree.lean)：
  首一除法的商次数与导数次数给出模 CLD 的严格次数界。实际对称
  约化使高次整数系数 c 满足 `2|c|≤m`；其模像为零，又有 m|c，
  因而 c=0。由此得到原整数 CLD、实际截断宽度和完整残差的次数界。
- [VanHoeijHenselLocal](../../lean/CLPoly/Refinement/VanHoeijHenselLocal.lean)：
  从实际规范化乘积证书消去单位，证明每个物理因子在返回模数下整除
  源多项式。将已有模 p 两两互素性提升到真实返回的 p^a。好素数下
  的无平方性在有限域上给出与导数互素，再提升到 p^a。物理首一性
  保持次数，结合已有模 p 因子不可约性推出局部因子正次数。
- [VanHoeijSeparationInput](../../lean/CLPoly/Refinement/VanHoeijSeparationInput.lean)：
  `CldSeparationData` 记录同一个实际 Hensel 因子数组和 CLD 返回数组
  的性质；构造定理接受真实 Hensel 证书、真实 cldPolys 运行及源输入
  条件，导出所有字段。再与真实 truncateCld 运行组合，得到完整残差
  的次数界、局部整除及导数互素性。

证书不包含最终输出不可约、保留行已分离或分组成功等假设；没有引入
另一套数学上任选的因子数组。用到的模 p 不可约性来自原始有限域因式
分解的既有正确性证明，并用于分析这些实际提升因子。

## 边界与剩余义务

上述源输入条件包括规范表示、非空、正次数、次数小于 2^31、正确
选素数与完整 Hensel 执行证书。相关机器界用于实际 CLD 底层执行。
这还不是“足够精度下当前控制器一定成功”的定理。

下一步需组合实际维数与 cutoff 的统一输入界，给结式使用的残差差值
设统一系数阈值，并把真实因子 Hensel 支持代入，推出保留行在同一支持
上恒定。还需实际分组验证完备性，以及保留充分精度信息的参数/重试策略。
当前 D 可能随 m 增大，三次尝试后回退仍未改变；性能验收仍未完成。

## 验证

- 三个模块构建通过；统一执行证书目标共 3568 jobs，含依赖缓存。
- [合并公理审计](../../lean/CLPoly/Refinement/VanHoeijSeparationAudit.lean)
  构建通过，共 3574 jobs；64 项定理均仅含标准公理
  `propext`、`Classical.choice`、`Quot.sound`，已对日志机器解析核对。
  日志：`/tmp/clpoly-vh-separation-input-audit.log`。
- strict refinement boundary、verified refinements 和 `git diff --check` 通过。
- 未修改生产 C++ 或生成模板，本轮不重复无关的原生/B2B 回归。

## 度量

- 耗时：未单独记录净证明时间。
- 迭代：共 6 次目标构建，含两次 CLD 次数证明修正和最终合并审计。
- Lean：新增三个模块共 329 行、十项定理、一份 Prop 证书；审计增加 11 行。
- 对应 C++：已有实际 Hensel 返回数组、规范化乘积和 CLD 运算，无新增原生逻辑。
- 问题与解决：natDegree_mul_le 的参数是隐式参数；投影 coeff 时先
  化简 lambda，再用 coeff_map 连接模像与整数系数，避免错误地把
  非定义等式交给 change。
- 放弃的方案：不新增抽象互素 oracle，复用已有素数幂互素提升定理。
