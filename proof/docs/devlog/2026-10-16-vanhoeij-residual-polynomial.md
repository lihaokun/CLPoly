# 完整残差多项式与局部 CLD 整除

日期：2026-10-16。

## 变更

先补充[自然语言证明](../../nl-proof/vanhoeij-eventual-separation.md) §7–8，
再新增四个证明模块：

- [VanHoeijColumnCoverage](../../lean/CLPoly/Refinement/VanHoeijColumnCoverage.lean)：
  实际交替低/高次列遍历的逆索引；成功追加循环的精确列数；实际
  宽度扫描的上下界；宽度外的截断系数为零；floor 为零时原系数
  是小于 D 的非负余数；从原 CLD 次数界推出实际截断宽度界。
- [VanHoeijResidualPolynomial](../../lean/CLPoly/Refinement/VanHoeijResidualPolynomial.lean)：
  以原 CLD 组合 P 和实际模数行系数构造 `H=P+C(m)Z`，同时证明
  同一个 H 的模同余和全系数高度界。接到实际 LLL 返回的 transform、
  保留扫描和已有 CLD 证书，不再把 CLD 规范性/系数界分别假设。
  原 CLD 的次数界成立时，也推出 H 的次数界。
- [VanHoeijLocalCld](../../lean/CLPoly/Algorithm/VanHoeijLocalCld.lean)：
  任意交换系数环上的对角、非对角 CLD 整除引理和加权有限和。
  适用于 Z/p^aZ，没有要求它是域。
- [VanHoeijLocalResidual](../../lean/CLPoly/Refinement/VanHoeijLocalResidual.lean)：
  使用实际 CLD semantic 字段及 H 的模同余，推出局部因子整除
  `H−e_i f'`；首一、整除、互素前提目前仍须从上游 Hensel 输出组合。

## 关键决定

H 使用原 CLD 的完整多项式，不能直接把所有高于截断宽度的项删掉。
这些项可能原本非零、仅在 floor 后消失。它们以原系数参与 H，
由 floor=0 的余数界控制，保持模同余。

有限和的模数修正只覆盖 d<W。W=0 时修正多项式为零，所有次数
走余数界；空 CLD 单独处理。全系数高度界不需要先假设原 CLD 的
次数上界，次数上界保留为独立待连接性质。

## 剩余义务

1. 将原 CLD 次数界从实际除法、导数和对称约化导出。
2. 从实际 Hensel/好素数导出局部整除与互素条件，包括 q_i 与 f' 的互素。
   已找到项目现有 `Algorithm/Hensel.isCoprime_lift_pk`，可提升模 p 互素性。
3. 组合维数、cutoff、系数界，接结式阈值和实际分组验证完备性。
4. 解决 D 随模数增大的策略问题；当前三次重试后回退尚未改变。
5. 完整生成对应关系与性能验收仍未完成，不合并研究分支。

## 验证

- 四个模块逐项构建通过；完整残差模块 3543 jobs，局部残差连接
  3545 jobs，均含依赖缓存。
- [合并公理审计](../../lean/CLPoly/Refinement/VanHoeijSeparationAudit.lean)
  构建通过，3551 jobs；54 项定理均仅依赖标准公理
  `propext`、`Classical.choice`、`Quot.sound`，已对日志机器解析核对。
  日志：`/tmp/clpoly-vh-residual-polynomial-audit.log`。
- strict refinement boundary、verified refinements、`git diff --check` 通过。
- 本轮没有生产 C++ 改动，不重复原生/B2B 回归。

## 度量

- 耗时：未单独记录净证明时间。
- 迭代：共 8 次目标构建，包含逐步新增功能后的构建和最终合并审计。
- Lean：四模块合计 539 行，新增 23 项定理；审计增加 24 行。
- 对应 C++：既有列遍历、共同截断、实际 LLL 和 CLD 输出，没有新增原生逻辑。
- 问题与解决：嵌套 get! 索引需括号；Bool false=true 分支用完整 simp
  约化；规范稀疏系数桥接先显式 rw，避免 simp 改写数组索引后失配。
- 放弃的方案：不使用先截掉原 CLD 高次项的构造；改用原组合加模数
  倍的有限修正多项式，保持每个次数的模同余。
