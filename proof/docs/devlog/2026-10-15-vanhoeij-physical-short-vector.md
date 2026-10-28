# van Hoeij 实际截断格短向量连接

日期：2026-10-15。目标 active，尚未合并。

## 本轮变化

新增实际 CLD 与截断格之间的组合证明：

- raw 坐标直接取原 CLD 数组的 sparseCoeff。
- 从真实截断循环的成功运行推出输出长度及系数 floor 对应。
- 原系数绝对值界与规范表示直接来自 ModularCldValue。
- 从支持集合上的逐次数同余选取进位，构造前部为 indicator、后部
  为负进位的整数系数向量，并证明其实际矩阵像满足所需平方范数界。

另证明无重复有界物理索引列表到 Fin 格行下标集合的求和一致性，
并对合法候选推出支持基数界、下标界和集合形式的进位等式。因此上一轮
Hensel 产生的合法候选可通过该桥消费，无需假设集合与数组选择一致。

## 边界

新短向量定理仍明确消费实际 CLD 循环不变量、规范初始格性质、支持
同余与整数 CLD 系数界；本轮消除了其中独立的 raw 系数界、floor
对应假设，并补齐列表到集合的桥。后续仍需将 Hensel、CLD、构格运行
与高度计算在单一入口组合，然后连接安全保留及候选完整验证。
没有以局部范数界替代端到端不可约性或纯 van Hoeij 分离定理。

## 验证

- 此轮仅新增证明与文档，未改 C++，未重复原生测试。
- 调试修正局部系数函数展开与求和桥显式参数；数学陈述保持不变。
- 完整目标构建通过（3554 jobs，含已有依赖），git diff --check 通过。
- 三个新增定理公理审计仅含 propext、Classical.choice、Quot.sound，
  无 sorry 或新公理。
- 新增 Lean 模块 124 行。总耗时未精确记录，不估造。

## 文件

- [短向量组合证明](../../lean/CLPoly/Refinement/VanHoeijShortVector.lean)
- [Hensel 支持](../../lean/CLPoly/Refinement/VanHoeijHensel.lean)
- [截断执行](../../lean/CLPoly/Refinement/VanHoeijTruncation.lean)
- [证明草稿](../../nl-proof/vanhoeij-certificate-repair.md)
