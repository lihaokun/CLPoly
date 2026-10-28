# van Hoeij 实际短向量入口组合

日期：2026-10-15。持续目标 active，尚未合并。

## 本轮变化

CLD 执行桥从实际 Hensel 输出推导非空性、首项系数为 1、首项次数
最大及模首一性。展开真实 invert 实现，显式证明模数大于 1 时系数
1 的求逆成功。由此消费实际 cldPolys 运行，提供整个 CLD 数组性质。

新增组合入口同时消费实际 cldBound、cldPolys、truncateCld、
makeInitialMatrix 和 buildCldMatrix 的运行结果，推导高度公式、
CLD 不变量及规范格，再获得真实因子短向量。结论保留同一个合法
物理索引列表的模支持等式、indicator 系数及范数界。

源输入机器次数界由规范表示、首项次数与 natDegree 对应以及
natDegree<2^31 推出；提升因子的机器次数界仍明确保留为入口前提。
本轮消除了调用者单独提供 CLD 不变量、格布局和高度公式的需求。

## 边界

当前入口是短向量阶段组合，不是完整 factorize 控制器。LLL 安全保留
已有局部证明，但还需与该入口和实际分组连接，再完成完整候选验证及
最终不可约性。重试控制器模型、纯 van Hoeij 分离界及性能验收仍未完成。

## 验证

- 完整目标构建通过（3557 jobs，含已有依赖），git diff --check 通过。
- 求逆特例最初以递归函数作为 simp 规则触发无限展开；改为显式三步
  欧几里得展开，不增加递归上限或引入假设。
- 新增两个 Lean 模块共 140 行；未修改生产 C++，未重复原生测试。
- 三项新增定理公理审计仅含 propext、Classical.choice、Quot.sound，
  无 sorry 或新公理。总耗时未精确记录，不估造。

## 文件

- [CLD 实际执行桥](../../lean/CLPoly/Refinement/VanHoeijExecution.lean)
- [组合入口](../../lean/CLPoly/Refinement/VanHoeijAssembly.lean)
- [证明草稿](../../nl-proof/vanhoeij-certificate-repair.md)
