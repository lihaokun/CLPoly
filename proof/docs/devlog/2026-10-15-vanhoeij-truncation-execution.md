# van Hoeij 实际截断循环与系数桥接

日期：2026-10-15。持续目标 active，尚未完成整体验收。

## 生产接口与严格模型

将 try-recombine 中的截断循环原样提取为 `__vanhoeij_truncate_cld`，生产
路径直接调用该 helper。这样 B2B 与严格模型验证同一执行入口，而非另写
一份测试专用截断。每个多项式逐项原位更新系数，随后立即规范化；顺序
保持不变，正除数前提显式断言。

严格模板新增内层系数循环、外层多项式循环及入口。保留循环边界、原位
数组赋值和规范化时机。模板 --check 通过，但仍不是 AST 翻译等价证明。

## 证明

- 内循环保留长度与次数，准确给出每个位置修改前/后的值。
- 正除数时模型的 Int.fdiv 等于先前截断误差证明使用的 Int.ediv。
- 外循环精确输出逐多项式系数变换与 normalization，保持原数组位置。
- 规范化保持多项式语义；严格降序次数排除重复次数，因此保持零的任意
  系数变换与取系数兼容。实例化为除以 divisor 后得到真实输出
  sparseCoeff = 原 sparseCoeff / divisor。
- 不对任意含重复次数的项数组假设 floor 与求和交换。

这提供了规范格短向量定理所需截断系数对应的实际执行证明，仍须在顶层
组合调用。真实 Hensel 支持、CLD 子集同余、完整候选验证及顶层控制器
也仍未闭合；整个算法尚不能合并。

## 验证

- 原生重组回归重新编译运行：382/382。
- 原有 190 个 B2B 与首批 61 个截断用例共同运行：251/251。
- 随后补充内部空多项式位置保持用例，截断测试集 62/62。两轮合计覆盖
  252 个不同用例。
- 单文件 Lean 检查通过；最终目标构建通过（3536 jobs，含已有依赖）；
  `git diff --check` 通过。五个关键定理的公理审计均仅含 propext、
  Classical.choice、Quot.sound，无 sorry 或新公理。
- 本轮第一个归纳证明编译失败来自递归假设的参数顺序，按实际签名修正，
  没有改动测试预期或弱化数学陈述。

## 度量

- 耗时：未精确计时，不估造。
- 迭代：4 次证明检查（初次目标构建失败，其后 3 次单文件检查通过），
  另有最终构建和公理审计。
- Lean：新增截断证明模块 199 行；生成模型新增三个函数。
- 对应 C++：提取约 14 行 helper，替换原内联循环为单次调用。
- 放弃方案：不使用独立的 B2B 截断副本；测试生产调用的同一 helper。

## 文件

- [C++ 生产实现](../../../clpoly/polynomial_factorize_univar.hh)
- [严格模板](../../cpp2lean_v2/tests/build_strict_vanhoeij.py)
- [截断证明](../../lean/CLPoly/Refinement/VanHoeijTruncation.lean)
- [证明草稿](../../nl-proof/vanhoeij-certificate-repair.md)
- [C++ 注册表](../../b2b/registry/cpp_registry.cc)
- [Lean 注册表](../../lean/B2B/Registry.lean)
- [截断向量](../../b2b/vectors/__vanhoeij_truncate_cld.json)
