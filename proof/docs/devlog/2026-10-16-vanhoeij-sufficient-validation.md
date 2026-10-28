# van Hoeij 充分精度下实际分组与验证成功

日期：2026-10-16。

## 结果与范围

新增 8 项定理，完成从原输入到完整有序目标、再到实际验证成功的组合。
最终 [hensel_lll_validation_complete](../../lean/CLPoly/Refinement/VanHoeijSufficientPrecisionValidation.lean)
同时返回实际 extractCandidates 的组数组和实际 validatePartition 的
complete=true 结果，不再假设完整目标数组、输出不可约性或验证成功。

该结论以同一实际 Hensel 输出、CLD 数据和截断执行、全宽规范格、
实际 LLL 执行及显式精度界为前提。它证明单次充分精度执行成功，
尚未证明实际列循环和精度控制器必能到达这些条件。

## 证明链与文件

先补充 [自然语言证明](../../nl-proof/vanhoeij-eventual-separation.md) §23–25。

- [VanHoeijFactorGroupBijection](../../lean/CLPoly/Refinement/VanHoeijFactorGroupBijection.lean)：
  保留原因子数组下标证明 assigned 满射；同一实际组的完整模支持与
  足够精度使两个因子相等，原数组 Nodup 给出下标相等，从而单射。
- [VanHoeijTargetPermutation](../../lean/CLPoly/Refinement/VanHoeijTargetPermutation.lean)：
  沿双射逆映射构造 targets，有限乘积重排证明乘积保持；逐项继承
  本原性、正首项、正次数和对应实际组的支持。
- [VanHoeijTargetsFromGroups](../../lean/CLPoly/Refinement/VanHoeijTargetsFromGroups.lean)：
  从逐因子精确分组选择 assigned 和支持，实例化双射并构造完整
  ExactPartitionTargets，消去原先整份目标数组的假设。
- [VanHoeijInputFactorFamily](../../lean/CLPoly/Refinement/VanHoeijInputFactorFamily.lean)：
  本原输入的好素数无平方因子性推出整数输入无平方因子；正首项使
  normalize f=f；规范因子列表无重复，乘积恰为 f，各项本原、正首项、
  正次数且在 Q 上不可约。这是用于分析真实执行的数学因子族，未用
  UFD 存在性替代实际重组算法。
- [VanHoeijSufficientPrecisionValidation](../../lean/CLPoly/Refinement/VanHoeijSufficientPrecisionValidation.lean)：
  实际 LLL 形状保证组提取成功；逐因子 hensel_exact_support_group 返回
  的数组由函数确定性识别为同一 groups。实际 Hensel 素元、支持和
  物理头性质消去模因子条件与表示条件，完整验证定理返回成功。

## 验证及问题处理

- 最终组合构建通过：3614 jobs。
- 严格 refinement 边界与 verified-refinements 注册检查通过。
- 合并公理审计 147/147 通过，3616 jobs；逐项解析确认仅使用
  propext、Classical.choice、Quot.sound。
- `git diff --check` 通过。
- 未改动生产 C++、严格生成模板或 B2B 输入，本轮没有重复原生回归。
- 编译修复仅涉及 Lean API：Array.size_ofFn/toList_ofFn 使用隐式参数；
  列表转数组再转列表是 List.toList_toArray；Multiset 的 Nodup 通过
  coe_nodup 和 coe_toList 连接。常数多项式整除式改写后补系数化简。

## 剩余关键义务

1. 把单次全宽成功接入真实 CLD/列循环入口，并明确所有机器界。
2. 当前参数除数 D 随 m 增大，现有阈值含 D，不能据 m 增长断言最终
   超过阈值；实际三次重试之后的 Zassenhaus 回退也不是纯 van Hoeij
   终止证明。必须修复参数/重试策略并证明阈值可达。
3. 更新生产和严格模型对应后，重新执行相应 B2B、公开证明构建、
   公理审计、原生回归、固定参考对照与 r=24/32 性能验收。
4. 全部验收满足前不合并。

## 度量

- 耗时：未单独记录可靠活动工时；五个最终模块的编译报告分别为
  34、16、3、19、11 秒，不能把它们当作总工作时间。
- 迭代：5 次目标构建调用，另 1 次合并审计构建；3 次目标构建因 API
  使用失败后修复，最后组合目标通过。
- Lean 新增：5 模块 437 行、8 项定理；审计增加 import 和 8 项输出。
- 对应 C++：0 行修改，组合既有真实 Hensel/LLL/分组/恢复验证模型。
- 放弃方案：无；没有把整份目标数组或输出不可约性保留为新入口前提。
