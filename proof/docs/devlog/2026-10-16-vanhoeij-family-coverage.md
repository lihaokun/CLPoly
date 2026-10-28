# van Hoeij 实际类覆盖与支持唯一性

日期：2026-10-16。

## 目的与结果

补齐“每个真因子都有精确组”到“每个实际组都有因子来源”的反向链。
先在 [自然语言证明](../../nl-proof/vanhoeij-eventual-separation.md) §20–22
写出循环不变量、边界分支和素元整除论证，再实现三个 Lean 模块，
共 6 项定理。

- [VanHoeijGroupCoverage](../../lean/CLPoly/Refinement/VanHoeijGroupCoverage.lean)：
  实际 partition 循环的每个已分配类都有物理列见证。seed 创建新见证，
  assign 保持已有槽；实际 collector 将见证放入返回组。空 rows 返回
  空数组，非空性结论不需要额外排除此分支。机器索引界保证同一列
  不能属于两个不同实际组。
- [VanHoeijSupportUnique](../../lean/CLPoly/Refinement/VanHoeijSupportUnique.lean)：
  两个正首项本原整除因子若有相同模支持，在同一结式阈值以上必相等。
  先交叉乘首项，对称区间唯一性得到整数等式；取 content 得到首项
  相等，再消去非零常数。没有不可约性或成功恢复前提。
- [VanHoeijFactorFamilyCoverage](../../lean/CLPoly/Refinement/VanHoeijFactorFamilyCoverage.lean)：
  从实际组的来源列出发，其原素数下的模因子是素元且整除输入。
  对数学因子列表乘积的素元整除证明找到对应因子；投影 Hensel 支持
  得到列在该因子精确组内，由实际类唯一性识别原组。

## 条件范围与剩余义务

覆盖定理以数学因子列表乘积等于输入、原素数下各模因子的素性和
不相伴性、以及逐因子的精确分组结论为条件。它不假设任意输出组
已有因子，也不假设验证成功；但仍需从实际 selection/Hensel 数据和
整数因子族实例化这些条件，证明完整有序目标数组的排列和乘积，
接入上一轮的验证完备性定理。实际精度控制器达到阈值的证明与修复
仍未完成，不能据本轮宣布纯 van Hoeij 终止或最终验收。

## 验证

- 三个新增模块首次编译均通过，最终目标 3601 jobs。
- 严格 refinement 边界及 verified-refinements 注册检查通过。
- 合并公理审计 139/139 通过，3611 jobs；已解析全部输出，仅使用
  propext、Classical.choice、Quot.sound。
- `git diff --check` 通过。
- 本轮无生产 C++、生成模板或 B2B 输入修改，未重复原生性能回归。

## 度量

- 耗时：本轮包含终态编译报告的约 15 分钟 family 模块构建；无单独
  可靠的活动编码计时，不将等待或上下文恢复计成精确工时。
- 迭代：3 次模块构建均首次通过，另 1 次合并审计构建。
- Lean 新增：3 模块 250 行、6 项定理；审计增加 1 个 import 和 6 项。
- 对应 C++：0 行修改；验证既有 partition/assign/collect 的实际模型。
- 放弃方案：无。没有把输出不可约或成功验证当作覆盖性前提。
