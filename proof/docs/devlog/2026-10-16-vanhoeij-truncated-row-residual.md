# 实际截断格保留行的残差界

日期：2026-10-16。

## 证明与实现

先扩展[自然语言草稿](../../nl-proof/vanhoeij-eventual-separation.md)
§5–6，再新增 [VanHoeijResidualBounds](../../lean/CLPoly/Refinement/VanHoeijResidualBounds.lean)
八项定理，并在 [VanHoeijParameters](../../lean/CLPoly/Refinement/VanHoeijParameters.lean)
补四项参数定理。

- 对任意有符号组合系数，证明 floor 的绝对值界、带权有限和界、
  模数行系数界与精确的反截断恒等式。没有沿用只适用于 0/1
  indicator 的截断误差证明。
- 若物理数据坐标绝对值≤Y、因子组合系数绝对值≤E，且项数≤r，
  则同一个残差 `H=Σ e_i C_i+z m` 同时满足所需模同余及
  `|H|≤D(2Y+3rE)`。
- 由同一次实际 LLL 返回的 transform 和原始规范格，证明实际物理
  行等于该整数组合。利用前 r 列的 `scale*e_i` 和 scale≥1，
  从物理坐标界导出组合系数界，没有额外假设 transform 系数很小。
- 接上实际 truncateCld 系数定理、实际保留扫描及已证坐标界 K，
  得每个已加入列的残差界 `D(2+3N)K`，其中 `K=N·2^N·B+1`。
- 实际参数选择始终满足 `height≤32D`，所以 residual≤32+2r。
  在 height≥0、r≥0、加入列数 J≤n 下，实际平方 cutoff 满足
  `B≤r·32²+n·(32+2r)²`。

## 边界与剩余义务

上述残差执行定理仍要求规范原 CLD 及其系数界；已有实际 CLD 语义
定理提供这些性质，但本轮尚未把所有上游调用组合进同一个充分精度
入口。残差定理仅覆盖已加入的列，下一步还要处理完整次数覆盖和
截断后消失的高次项，构造一个统一残差多项式。

参数界中的 J≤n 尚须从实际 CLD 宽度及列循环导出。K 的输入统一界
还需结合维数上界。最关键的是 H 的界仍含 D；当前位数限制使 D
可能随模数增长，因此本轮仍没有证明当前重试策略的纯 van Hoeij
充分精度终止。生产策略、三次重试预算均未修改。

## 验证

- [分离审计](../../lean/CLPoly/Refinement/VanHoeijSeparationAudit.lean)
  构建通过，3545 jobs（包含依赖缓存）。
- 31 项分离基础及执行定理全部仅依赖标准公理
  `propext`、`Classical.choice`、`Quot.sound`；对日志作了机器解析核对。
- 运行日志：`/tmp/clpoly-vh-residual-audit.log`。
- 公开 Pipeline、FinalAxiomAudit、VanHoeijRetryAudit 重建通过，3611 jobs；
  当前重试/FactorZZ 链条 21 项定理仍仅含标准公理。
  日志：`/tmp/clpoly-vh-residual-public-build.log`。
- strict refinement boundary 与 verified refinements 两项检查通过。
- `git diff --check` 通过。无原生代码变化，本轮未重复原生/B2B 回归。

## 度量

- 耗时：未单独记录净证明时间。
- 迭代：残差模块 5 轮构建；参数扩展 3 轮（最后一轮并入审计）。
- Lean：新增残差模块 288 行；参数模块增加 63 行；审计增加 14 行。
- 对应 C++：已有共同截断、规范格和 LLL 物理输出；没有新增生产代码。
- 问题与解决：绝对值三角不等式的隐式参数导致一次 elaboration
  超时，改为显式表达式；使用双边 `add_le_add` 避免同名左右平移
  引理的方向混淆；对 ZZ 定义显式展开后再交给整数算术策略。
- 放弃的方案：不把 GS 范数当作物理行范数，不把截断残差与其模
  关系独立假定，也不把本轮含 D 的界宣称为最终统一高度界。
