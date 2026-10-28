# van Hoeij 输入固定的截断上界

日期：2026-10-16。

## 修复与理由

原参数规则随 m 增大而持续扩大 D，使有效格模数约被限制在 r+16 位。
上一轮的分离阈值含 D，因此不能据 m 增长推出最终达到充分精度。

保留位数截断，在 [生产参数选择](../../../clpoly/polynomial_factorize_univar.hh)
中增加仅依赖输入的上界 D₀²；D₀ 为原高度倍增循环的结果。实际
D=min(D₀²,max(D₀,位数截断除数))，未触发位数截断时仍为 D₀。
D₀≥1 保证夹界不破坏原先 32D≥H，残差和短向量证明保持。

同步修改 [严格模板](../../cpp2lean_v2/tests/build_strict_vanhoeij_parameters.py)
并重新生成实际模型。该模板仍是显式模型，不是普遍 AST 翻译或
C++/Lean 等价性证明。

## 证明

先补充[自然语言证明 §26](../../nl-proof/vanhoeij-eventual-separation.md)，再形式化。

- [VanHoeijParameters](../../lean/CLPoly/Refinement/VanHoeijParameters.lean)：
  新增实际 D 的输入固定夹界，重接正性与高度界；原 cutoff 定理通过。
- [VanHoeijPrecisionThreshold](../../lean/CLPoly/Refinement/VanHoeijPrecisionThreshold.lean)：
  参数选择最终成功；结式阈值对非负 D 单调；用 D₀² 给出只依赖输入
  的统一门槛；每个 p>1 都存在正指数 N，使所有 a≥N 的 p^a 让实际
  参数选择成功，且实际 D 对应的结式、坐标界成立。

没有假设位数结果等于数学对数，夹界证明覆盖严格模型的实际位运算。
统一门槛只用于证明，不在生产代码中求值。

## 验证

- 公开 Pipeline.L1、B2B 注册/驱动和合并分离审计重建成功：3654 jobs。
- 公理输出逐项解析：154/154，仅 propext、Classical.choice、Quot.sound。
- 模板、向量生成器的 --check、严格边界、verified-refinements 和补丁
  格式检查通过。
- [独立闭式向量](../../b2b/generate_vanhoeij_parameters_vectors.py) 增加
  固定上界过渡边界，参数 1534 项，宽度 108 项；连同重组 51 项、LLL
  27 项和外层入口 18 项，实际 C++/Lean B2B 共 1738/1738 通过。
- [原生重组测试](../../../test/test_recombine.cc) 增加有效精度增长回归。
  新编译 profiling 模式 383/383 通过；仓库原样完整原生测试脚本
  27 个测试文件通过，0 失败。
- r=24/32 的三族固定 FLINT 对照再次 18/18 通过，新版均证书完成且零回退。
  同配置旧/新二进制六组中位数变化约 −2.1% 至 +2.6%；本轮运行期间有
  后台 B2B 解释器负载，只作诊断，不作最终安静环境性能验收。
  详见[参考报告](../../../docs/test-reports/vanhoeij-flint-reference.md)。

## 剩余义务

此次修复证明了实际参数策略存在充分精度，仍须将单次充分精度成功
接入真实列循环，并修改/证明实际重试控制器能到达该指数，明确机器
指数溢出边界。当前三次重试后 Zassenhaus 回退没有因此成为纯
van Hoeij 终止证明，最终对应审计和性能验收仍未结束。尚未合并。

## 度量

- 耗时：未单独记录可靠活动工时；包含公开证明重建、原生重编译和
  Lean 解释器 B2B 的等待，不把这些误报为编码时间。
- 迭代：1 次原参数目标构建；统一阈值目标 2 次（结构类型标注修复）；
  加入素数幂存在性后组合构建 2 次（幂表达式括号修复）。
- Lean：新增阈值模块 127 行、6 项定理；原参数模块新增 1 项夹界定理
  并重写 2 项旧证明，审计新增 7 项。
- C++：参数选择局部上界修改，新增 5 条精度增长回归断言。
- 放弃方案：无。没有用有限重试或回退替代充分精度存在性。
