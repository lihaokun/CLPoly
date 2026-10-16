# van Hoeij 完成证书基础引理

日期：2026-10-15。

从[证明草稿](../../nl-proof/vanhoeij-certificate-repair.md)形式化了两个模块：

- [VanHoeijBound](../../lean/CLPoly/Algorithm/VanHoeijBound.lean)：通过系数卷积、加权 binomial 恒等式、Vandermonde、Mahler 乘法性证明非首一 CLD 界，最后得到整数 n·2^(n−1)·L1 版本。
- [VanHoeijCertificate](../../lean/CLPoly/Refinement/VanHoeijCertificate.lean)：尾部扫描规格、截断误差、安全尾部移除、列相等的线性闭包、同组真因子一致性及完整支持分组。

两个模块完整 `lake build` 通过。关键导出定理的 `#print axioms` 仅包含 `propext`、`Classical.choice`、`Quot.sound`。未引入输出不可约性公理或 sorry。

这不是新原生实现的端到端证明：严格生成控制器仍是旧版；格包含、数组投影、实际 candidate 支持与最终 FactorZZ 定理仍须连接。旧模型的编译成功不证明新 C++。

## 度量

- 耗时：未记录精确起止，不估填。
- 迭代：修正过递归 simp 过度展开和 Mahler 非负性显式证明两个问题；独立引理检查及最终模块构建均通过。
- Lean 新增：279 行。
- 对应 C++：CLD 界、尾部保留、分组和控制器约 440 行 diff，尚未完成 1:1 翻译。
- 放弃的方案：直接展开复根索引的 CLD 界证明；等价的卷积/Vandermonde 路径更直接复用当前 Mathlib。
