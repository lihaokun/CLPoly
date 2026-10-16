# van Hoeij 验证器乘积证明

日期：2026-10-15。目标 active，未合并。

## 本轮变化

沿实际 validatePartitionLoop 运行，对剩余组数归纳，证明返回因子
数组的乘积等于进入该循环时的已接受因子乘积乘 remaining。

每轮使用实际 exactDivmodRaw 输出的商与余数；空余数给出当前余商
等于新因子乘商。商的规范性由实际除法证明获得，商整除非零余商
给出次数不增，因而把初始次数范围传递到终点。

终点按实际物理检查证明余商为数学常数 1：数组长度为 1、首项系数
为 1、get_deg=0，结合规范表示及 Int64 合法次数范围还原次数。
排序的排列性质保持乘积。入口定理证明 validatePartition 返回成功
时，实际 after.factors 的乘积等于输入多项式。

## 边界

该结论是乘积正确性，不是不可约性。还需验证因子本原性、物理组与
模支持对应，以及组合已有 LLL 分组性质。初始规范表示及次数小于
2^63 是明确输入前提；后续入口已有更强的 32 位次数范围可提供它。

## 验证

- 单文件检查与完整目标构建通过（3561 jobs，含已有依赖），
  git diff --check 通过。新增证明模块 154 行。
- 三项关键定理公理审计仅含 propext、Classical.choice、Quot.sound，
  无 sorry 或新公理。
- 调试修正构造器等式、Option.none 与 some 的不可能分支，以及
  toPoly 展开后的改写匹配；未弱化数学目标。
- 只修改证明和文档，未重复上一轮原生 382/382 与 B2B 476/476。
- 总耗时未精确记录，不估造。

## 文件

- [乘积证明](../../lean/CLPoly/Refinement/VanHoeijValidationProduct.lean)
- [严格模型](../../lean/CLPoly/Generated/StrictVanHoeij.lean)
- [证明草稿](../../nl-proof/vanhoeij-certificate-repair.md)
