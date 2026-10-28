# van Hoeij 实际 LLL 接口的完整执行证明

日期：2026-10-16。

## 缺口与解决

已有 LLL 的行列式/索引严格下降证明与返回值正确性，并不自动排除
实际数组访问的 RawExec 错误。要消去重组入口的正常返回前提，必须
证明它的内部操作确实能返回。本阶段完成 LLL 部分，先编写
[自然语言证明 §33–36](../../nl-proof/vanhoeij-eventual-separation.md)。

- [VanHoeijInitialExecution](../../lean/CLPoly/Refinement/VanHoeijInitialExecution.lean)：
  方阵上对角赋值、实际初始矩阵构造必返回；任意非空方阵的实际
  LLL 初始化不会因形状访问出错，复用已有 Gram–Schmidt 初始化执行。
- [VanHoeijRowExecution](../../lean/CLPoly/Refinement/VanHoeijRowExecution.lean)：
  四行同步消元、矩阵/变换同步修改、实际 sizeReduceAt 均正常返回。
  舍入使用实际 roundQQ_eq_round，mu 更新使用实际循环的数组等式。
- [VanHoeijStepExecution](../../lean/CLPoly/Refinement/VanHoeijStepExecution.lean)：
  下降索引的额外消元循环以及整数/有理矩阵行交换正常返回；证明
  QQ 行交换保持方阵尺寸，为后续写入提供物理边界。
- [VanHoeijLLLStepTotal](../../lean/CLPoly/Refinement/VanHoeijLLLStepTotal.lean)：
  Lovász 推进与交换两分支均正常执行且返回 k>0。复用已有真实变换
  关系、数学不变量保持和严格 concreteLLLRank 下降，构造主循环返回。
- [VanHoeijLLLExecutionComplete](../../lean/CLPoly/Refinement/VanHoeijLLLExecutionComplete.lean)：
  正常完成短行插入与收集，再连接实际初始化、主循环及实际最终
  cached norms，得到 lllReduceWithNorms 正常返回存在性。

最终 LLL 定理只要求非空方阵、非零行列式的 ConcreteLLLInputValid
和任意 cutoff，不要求 LLL 成功执行、结果或缓存范数假设。循环终止
使用原有实际度量，不添加 fuel、回退或结果 oracle。

## 验证与限制

完整 LLL 接口目标通过；合并公理审计 182/182，3648 jobs，逐项
确认只依赖 propext、Classical.choice、Quot.sound。新增模块无 sorry
或新公理声明，补丁格式与行尾空白检查通过。严格边界和
verified-refinements 检查已通过。没有修改 C++ 或严格运行模型，
不重复未受影响的原生/B2B 回归。

这只消去 **LLL 内部**的正常执行前提，尚须证明实际 CLD 生成、
构格和整个分区验证不会出错，再接完整重组入口。外层目前仍是三次
尝试后回退，必须继续修改与证明精度增长；机器范围、最终对应审计
和安静环境性能验收未完成，尚未合并。

## 度量

- 耗时：未单独记录可靠活动工时。
- 迭代：共 12 次目标/组合构建调用，含最终合并审计；初始化 2 次，
  行消元至完整接口 10 次。各失败均保留相同源码目标继续修复。
- Lean：新增 5 模块、357 行、14 项定理；审计追加 14 项。
- C++：没有修改，证明既有实际初始化/消元/交换/主循环/收集调用。
- 工程问题：限定实际矩阵及状态名称消除与数学模型重名；依赖
  `match` 使用分支展开与证明无关性衔接，不直接作不合法的替换；
  区分带证明的 `dif` 与普通 `if`；显式识别实际状态和终止实例。
- 放弃方案：无。未将数学不变量或有限递归本身误作错误不可达证明。
