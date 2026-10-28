# 实际恢复步骤与完整验证循环的条件完备性

日期：2026-10-16。

## 变更与理由

上一轮已把实际输出组连接到真实支持乘积，但这不等于验证循环必定
成功。本轮先写[自然语言证明](../../nl-proof/vanhoeij-eventual-separation.md)
§17–19，再新增六个模块、19 项定理，覆盖机器次数、非首一恢复和
完整验证循环的执行。

- [VanHoeijPartitionProductExecution](../../lean/CLPoly/Refinement/VanHoeijPartitionProductExecution.lean)：
  有效物理索引和非零模数使实际 partitionProductLoop 返回；次数值
  精确等于机器加法 foldl，再证明它等于自然数次数和的一次 Int64 转换。
  此阶段不把截断后的机器值直接当作无限精度整数。
- [VanHoeijPartitionDegree](../../lean/CLPoly/Refinement/VanHoeijPartitionDegree.lean)：
  投影真实支持等式到原素数，利用首项非零和物理首一性保持次数，
  得到该和恰为真实因子次数；次数小于 2^63 时证明实际正次数检查。
- [VanHoeijPrimitiveRecovery](../../lean/CLPoly/Refinement/VanHoeijPrimitiveRecovery.lean)：
  对正首项的缩放本原因子，证明实际 contentLoop、符号分支及
  primitiveRaw 恰好消去缩放常数。非首一情形恢复的是 lc(Q)F，不能
  假设本原化没有作用。
- [VanHoeijRecoveryHeight](../../lean/CLPoly/Refinement/VanHoeijRecoveryHeight.lean)：
  F、Q 均整除原输入时，统一 B 控制它们的系数，故缩放目标满足
  `2|lc(Q)F_k|≤2B²≤(2n)!B^(2n)`。已有结式阈值同时覆盖恢复所需的
  严格半模数界，不新增“恢复成功”假设。
- [VanHoeijPartitionRecoveryStep](../../lean/CLPoly/Refinement/VanHoeijPartitionRecoveryStep.lean)：
  实际构造 product、symmetric、primitive 和精确除法结果，证明因子
  是 F、商是 Q，空因子、次数不符、非正次数与非空余数分支都不可达。
  结论是实际 validatePartitionLoop 与下一次调用的等式，未假设已有
  成功运行轨迹。
- [VanHoeijPartitionValidationComplete](../../lean/CLPoly/Refinement/VanHoeijPartitionValidationComplete.lean)：
  记录与实际 groups 对齐的数学 targets：等长、乘积为原输入、本原、
  正首项、正次数及每项的精确模支持。以 remaining=targets 后缀乘积
  为不变量，逐组应用真实步骤，证明整个循环及 validatePartition
  入口成功。末尾从规范表示与整数多项式等于 1，逐项证明实际数组
  size=1、get_deg=0、头系数=1，未跳过机器检查。

## 条件范围与剩余工作

完整验证定理的 targets 条件没有不可约性字段，也没有任何成功执行
字段。它是本轮明确的条件完备性定理，**尚不是从实际分组无条件导出
完整 targets 的最终组合**。

下一步需要证明每个实际输出组非空、同一物理列不会出现在不同组，
并把此前“每个真实不可约因子有一个精确组”的方向扩展为全部组与
完整因子分解的对应，从而构造 targets 记录。随后仍需连接实际参数、
完整列控制器和可达的精度策略。当前随 m 增大的 D 与三次失败后
Zassenhaus 回退尚未修改；有限回退仍不代替纯 van Hoeij 终止证明。
生成对应最终审计及性能验收仍未完成，研究分支未合并。

## 验证

- 六个模块均构建通过。完整验证入口目标共 3597 jobs，含依赖缓存；
  最终单模块编译约 24 秒。
- strict refinement boundary、verified refinements 和 `git diff --check`
  均通过。
- 合并分离/分组/恢复审计 133/133 通过，3608 jobs；仅含 propext、Classical.choice、Quot.sound，已逐项解析核对。
- 本轮无生产 C++、生成模板、B2B 注册或原生测试变更，未重复原生回归。

## 度量

- 耗时：未单独记录净证明时间。
- 迭代：共 13 次目标构建调用，含合并审计。组乘积 2 次；次数和本原
  恢复各首次构建后合并修正 1 次；恢复高度 1 次；单步验证 4 次
  （首次因依赖失败而未进入该模块）；完整循环 2 次；审计 1 次。
- Lean：新增六个模块共 623 行、19 项定理及一个仅含数学性质的 Prop
  记录；审计增加 19 项输出及一条 import。
- 对应 C++：已有 partitionProductLoop 与 validatePartitionLoop，
  没有新增原生逻辑。
- 问题与解决：对整数首项、绝对值转 natAbs 和有符号机器比较使用
  显式类型，避免隐式 toNat 转换改变中间命题；模常数用显式目标
  统一 castRingHom 应用与整数 cast；末尾物理检查显式搬运 getElem!。
- 放弃的方案：不假设本原化是恒等操作，不忽略 Int64 截断，不用
  整数多项式等于 1 直接代替实际末尾检查，不把条件 targets 记录
  宣称为已经从所有实际组构造完毕。
