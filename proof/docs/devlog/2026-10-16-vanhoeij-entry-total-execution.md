# van Hoeij 重组入口的完整执行与充分精度成功

日期：2026-10-16。

## 缺口与解决

上一阶段的充分精度入口定理仍以实际入口正常返回为前提。本阶段先
补写[自然语言证明 §37–42](../../nl-proof/vanhoeij-eventual-separation.md)，
逐个排除其内部算术、断言和数组访问错误，再组合成返回存在性。

- [VanHoeijBuildExecution](../../lean/CLPoly/Refinement/VanHoeijBuildExecution.lean)：
  实际填写 CLD 列、追加模数行与规范格构造正常返回，覆盖零目标等边界。
- [VanHoeijDivmodReduced](../../lean/CLPoly/Refinement/VanHoeijDivmodReduced.lean) 与
  [VanHoeijDivmodComplete](../../lean/CLPoly/Refinement/VanHoeijDivmodComplete.lean)：
  实际模除法的余式系数位于 [0,m)。首一除数的模像整除源时，实际
  余式数组为空。模环可有零因子，次数论证显式使用 monic 乘积次数
  定理，不错误地假设模环是整环。
- [VanHoeijCldExecution](../../lean/CLPoly/Refinement/VanHoeijCldExecution.lean)：
  模除法、空余式断言、导数、模乘与对称恢复依次正常返回，全部条件
  从同一次实际 Hensel 输出导出，使用 concreteDivmodTermination。
- [VanHoeijTruncateExecution](../../lean/CLPoly/Refinement/VanHoeijTruncateExecution.lean)：
  正除数保证真实系数截断与逐多项式规范化循环正常返回。
- [VanHoeijValidationExecution](../../lean/CLPoly/Refinement/VanHoeijValidationExecution.lean)：
  任意合法索引分组的验证正常返回，包括错误分组被正常拒绝，不预设
  候选就是因子或整除当前余项。
- [VanHoeijControllerTotal](../../lean/CLPoly/Refinement/VanHoeijControllerTotal.lean)：
  连接实际构格、LLL、候选提取与验证，证明任意列数回合及实际列数
  增长循环正常返回；沿用真实列数剩余量严格下降。
- [VanHoeijEntryTotal](../../lean/CLPoly/Refinement/VanHoeijEntryTotal.lean)：
  构造实际高度界、单元素分组检查、参数选择、CLD/截断及列循环返回。
  最终定理从输入及实际 Hensel 条件推出入口返回、complete、精确乘积
  和实际输出不可约，不再假设入口或任何内部步骤已正常返回。

## 验证与限制

入口组合目标已通过编译，合并公理审计 206/206、3656 jobs 通过，
逐项只依赖 propext、Classical.choice、Quot.sound。严格证明边界、
verified-refinements、补丁格式和新增文件无 sorry/axiom/admit 检查
已通过。本阶段只新增证明与文档，未修改 C++ 或严格运行模型，因此
沿用上一轮 383/383 原生重组、27 文件完整原生回归、1738/1738 B2B
结果，不重复未受影响的运行测试。

本结果以实际 Hensel 输出正确性为输入，并明确保留源次数和因子数
机器范围条件。它完成一次重组入口的执行和充分精度成功证明；外层
仍是三次尝试后 Zassenhaus 回退，不能视为纯 van Hoeij 终止。
真实精度增长、机器上限处理、最终生成对应审计与性能验收仍须完成，
尚未合并。手写严格模板的核对和 B2B 不构成普遍 AST 翻译等价证明。

## 度量

- 耗时：未单独记录可靠活动工时。
- 迭代：16 次目标/组合构建调用，含最终审计；构格 2 次，模除法
  5 次，CLD/验证 3 次，控制器/入口 5 次，审计 1 次。
- Lean：新增 8 模块、755 行、24 项定理；审计新增 24 项。
- C++：本阶段修改 0 行，验证已有重组入口、CLD、截断、构格及验证循环。
- 工程问题：限定根命名空间以消除 Generated 重名；通过公开等式使用
  私有数组辅助函数；展开 dependent match 并识别实际返回值；统一局部
  cutoff 别名以衔接同一实际 LLL 输出。
- 放弃方案：不使用整环的一般整除次数界处理 Z/mZ，改用首一乘积
  次数等式；已先修正自然语言草稿，再完成相应形式化。
