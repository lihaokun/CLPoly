# 实际候选/Hensel 验证与有限模支撑识别

日期：2026-10-15。

## 变更及理由

- [VanHoeijCandidateValidation.lean](../../lean/CLPoly/Refinement/VanHoeijCandidateValidation.lean)
  从实际 `extractCandidates` 的成员来源与 Int32 转换范围推出合法索引。
  保留行为空时候选为空，与完整验证成功矛盾，无需额外假设非空。
  组合真实选素数/Hensel 输出，导出模因子不可约、正模数、素数整除
  模数与首项非零条件，再接入完整验证的因子物理组来源定理。
- [VanHoeijSupportSeparation.lean](../../lean/CLPoly/Refinement/VanHoeijSupportSeparation.lean)
  对有限列表归纳证明：不同物理索引的素元不相伴时，原子整除选中
  乘积等价于索引属于支撑；带相伴关系的因子整除因此推出支撑包含。
  非空支撑包含于组且组内成员资格一致时，无重复性给出列表置换，
  从而两者模像相伴。
- [自然语言证明](../../nl-proof/vanhoeij-certificate-repair.md) 先记录对应论证。

## 验证及范围

候选/Hensel 组合模块构建通过（3564 jobs）；三条主定理公理审计
仅含 `propext`、`Classical.choice`、`Quot.sound`。未新增公理或 sorry。
有限支撑模块最终构建通过（3565 jobs）；原子整除识别、支撑包含仅
依赖 `propext`，同组支撑相伴仅依赖 `propext` 和 `Quot.sound`。
`git diff --check` 通过。本步没有改动 C++ 或生成模型，不重复原生和
B2B 测试，不将之前测试说成本步重新运行。

有限模支撑引理仅处理给定模原子列表，不使用整数输出不可约性假设，
也未调用整数不可约分解存在性替代算法。

组合入口仍保留真实 Hensel 执行证书、输入规范性、机器索引范围；
有限支撑引理的无重复、不同索引不相伴及组内支撑一致性尚需与实际
候选、好素数和短向量证明实例化。完整不可约性、严格顶层控制器、
纯 van Hoeij 充分精度终止和性能验收仍未完成。目标继续进行，未合并。
