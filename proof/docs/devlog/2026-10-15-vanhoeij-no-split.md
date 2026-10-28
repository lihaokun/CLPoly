# 同组物理支撑的一致性及不可分解判据

日期：2026-10-15。

## 变更与动机

- [VanHoeijPhysicalGrouping.lean](../../lean/CLPoly/Refinement/VanHoeijPhysicalGrouping.lean)：
  从物理列表成员恢复实际候选位置，将格证明的同组成员等价关系通过
  Int32 往返转换转移到有界物理索引。保留 count≤2^31 的机器范围。
- [VanHoeijNoSplit.lean](../../lean/CLPoly/Refinement/VanHoeijNoSplit.lean)：
  无平方模像与支撑相伴给出无重复列表；非单位真因子有非空支撑。
  包含关系及组内一致性使其占满整组，消去非零真因子得到补因子模像
  为单位。再利用本原性/好素数确保非单位整数因子的模像非单位，得到
  条件式整数不可约性判据。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md)：记录数学链条及边界。

## 证明边界

定向构建 `CLPoly.Refinement.VanHoeijNoSplit` 最终通过（3570 jobs）。
三条主定理公理审计仅含 `propext`、`Classical.choice`、`Quot.sound`；
无新增公理或 sorryAx。`git diff --check` 通过。
日志为本机 `/tmp/clpoly-vh-no-split-build.log` 和
`/tmp/clpoly-vh-no-split-audit.log`。

整数判据中的除数支撑条件仍是显式参数，需要由真实 Hensel 支撑、
真因子短向量、LLL 和候选提取组合实例化。没有把该参数或判据包装成
已完成的端到端证明。还需将同一个模 p^a 支撑降至 p，并连接自然数
支撑与物理 Fin 列表，保持见证一致。

本步无 C++/生成模型变更，未重复原生或 B2B 测试。严格控制器、终止性、
参考对照、性能验收仍未完成；目标继续进行，未合并。
