# 结式分离与实际保留 GS 范数界

日期：2026-10-16。

## 草稿与数学边界

先完成[自然语言草稿](../../nl-proof/vanhoeij-eventual-separation.md)，再
形式化本轮引理。正式发表论文只作证明路线参考；所有新结论均由 Lean
检查，未引入论文定理公理。基础结式定理中不可约性描述待分析的真实
输入除数，绝不假设重组输出不可约。

## 新模块

[VanHoeijSeparation](../../lean/CLPoly/Algorithm/VanHoeijSeparation.lean)，207 行：

- 首一正次多项式整除常数，迫使该常数为零；适用于有零因子的环。
- 在 `ZMod M` 中用结式 Bézout 恒等式推出整个 M 整除整数结式。
- 从 Leibniz 展开给出矩阵及 Sylvester 结式的显式 `d! B^d` 界。
- 充分模数下结式为零；映到 Q 后由真实因子的不可约性得到整除。
- Gauss 下降覆盖零被除式及非本原被除式，不偷加残差本原条件。
- 两个有界整数坐标在同一首一局部因子上产生相同 CLD 残差，且该局部
  因子与导数互素时，模数超过二倍坐标界即可推出坐标相等。

[VanHoeijReducedBounds](../../lean/CLPoly/Refinement/VanHoeijReducedBounds.lean)，155 行：

- 从实际 `lllReduceWithNorms` 运行取出同一最终状态的全部 size-reduced
  和 Lovász 条件；初态 k=1 的两个前缀性质真空成立。
- 结合真实 GS 范数正性，推出实际相邻范数 `norm[i−1]≤2 norm[i]`。
- 实际尾部扫描的非空保留前缀，末项范数不超过 cutoff。
- 迭代相邻界，得到每个实际保留 GS 范数
  `norm[i]≤2^(t−1−i) cutoff`，没有假设范数单调。

[审计模块](../../lean/CLPoly/Refinement/VanHoeijSeparationAudit.lean)，20 行。

## 验证

- 两个模块与十五项导出引理审计构建成功，3539 jobs（含依赖缓存）。
- 仅 `propext`、`Classical.choice`、`Quot.sound`；没有 `sorryAx` 或新公理。
- 六个 van Hoeij 模板、四个向量生成器可重现检查通过；最终严格边界
  与公开注册检查同步执行。没有改动这些检查规则。

## 尚未证明

GS 范数还须经实际 Gram 等式转成基行坐标界；实际截断格的小残差
高度界、Hensel 支持到局部 CLD 条件、分组验证的执行完备性及充分
精度参数策略尚待连接。新结式模块是该连接的基础，不是完整控制器
终止定理，也没有消除手写模板与通用 AST 等价之间的边界。

## 度量

- 耗时：与外层 B2B、原生回归并行推进，未单独记录净证明时间。
- 迭代：结式模块 7 轮构建；实际 GS 模块 2 轮；独立及合并审计各一次。
- Lean 新增：382 行，含 20 行公理审计。
- 对应 C++：复用现有实际 LLL、GS 缓存导出和尾部保留扫描；无算法改动。
- 放弃的方案：直接套用早期预印本或渐进格算法的终止定理；两者都不能
  替代当前具体截断和每轮重建策略的执行证明。采用较松但初等可检验的
  行列式界；不先追求更紧的 Hadamard 常数。
