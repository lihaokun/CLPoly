# van Hoeij 真实入口的充分精度不可约性

日期：2026-10-16。

## 证明与边界

此前正目标列循环和零列宽分支分别已有完成性，但调用者还须给出
实际 CLD、截断和列循环执行。本次先补充
[自然语言证明 §32](../../nl-proof/vanhoeij-eventual-separation.md)，再直接
展开真实 tryRecombine，从其正常返回执行中提取全部中间事实。

- [VanHoeijEntryWidth](../../lean/CLPoly/Refinement/VanHoeijEntryWidth.lean)：
  由实际 CLD 次数、截断宽度界和输入次数范围，证明实际 Int32 扫描
  与自然完整列宽相等，无需调用者另给各截断多项式首项界。
- [VanHoeijEntrySufficient](../../lean/CLPoly/Refinement/VanHoeijEntrySufficient.lean)：
  识别实际高度，排除选参拒绝；若单元素验证成功直接得到 complete；
  否则从真实 CLD 执行构造分离数据，零列宽与已记录验证拒绝矛盾，
  正列宽调用列循环完成性。再连接已有实际输出不可约性证明，得到
  complete=true、输出乘积为输入、每个实际输出因子在 Z[x] 不可约。

最终定理没有 complete、输出不可约性、CLD/截断/LLL 执行或完整
目标数组假设。它保留实际 Hensel 提升、源输入合法性、机器大小界、
统一充分精度以及**入口正常返回**前提。内部操作无错执行、外层
指数可达性尚须证明；本次不是无条件终止结论。

## 验证

宽度目标 3571 jobs、入口目标 3634 jobs 通过；加入不可约性组合后
合并审计 168/168，3643 jobs，逐项确认只依赖 propext、Classical.choice、
Quot.sound。没有修改 C++，原生/B2B 验证沿用此前相同运行代码结果。

## 度量

- 耗时：未单独记录可靠活动工时。
- 迭代：2 次目标构建和 1 次合并审计，均首次通过。
- Lean：新增 2 模块、185 行、3 项定理；审计追加 3 项。
- C++：未修改，连接已有真实入口模型。
- 放弃方案：无。明确保留正常返回前提，没有把错误分支或机器溢出
  忽略为已证终止。
