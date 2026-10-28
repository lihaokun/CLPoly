# Vandermonde sorry 消除 + 多变量形式化完成

**日期**：2026-03-24

## 做了什么

消除了 `Wang.lean` 中最后一个 sorry：`vandermonde_solve_unique`。该定理证明 Vandermonde 系统 `V·d = v` 在节点两两不同时有唯一解。

**根因**：`∃!` 宏展开生成 `(fun d => P d) witness` 形式的 beta-redex，`rw [mulVec_mulVec]` 模式匹配失败（找不到 `M *ᵥ (N *ᵥ v)` 子项）。

**修复**：在 `rw` 前加 `show (Matrix.vandermonde θ).mulVec ((Matrix.vandermonde θ)⁻¹.mulVec v) = v` 强制 beta 归约。同时将 `mul_nonsing_inv` 的参数从 `_ (roundabout conversion)` 改为显式 `(Matrix.vandermonde θ) h_det_unit`。

## 当前总体状态

**5745 行 Lean，0 sorry（验证文件）。`lake build` 3071 jobs 全通过。**

### L2 模型与 C++ 算法一致性评估

#### 完全 1:1 的模块

| 模块 | 文件 | 对应 C++ | 一致性 |
|------|------|---------|--------|
| **Hensel 提升** | Hensel.lean (887 行) | polynomial_factorize_univar.hh 307-597 | 完全 1:1：单步 + 度数保持 + LC 保持 + 互素传播 + 多因子归纳 |
| **DDF** | DDF.lean (435 行) | polynomial_factorize_zp.hh 247-292 | 完全 1:1：ddfLoop 递归结构、6 不变量、终止条件 |
| **SQF** | SquarefreeZp.lean (1501 行) | polynomial_factorize_zp.hh 108-180 | 完全 1:1：Yun 算法循环、expand/contract、导数为零检测 |
| **Hensel 唯一性** | Recombine.lean | polynomial_factorize_univar.hh | 完全 1:1：归纳 on k，互素 + LC 控制 |
| **因子恢复** | Recombine.lean | symmetric_recovery + factor_recovery | 完全 1:1：对称约化 + Mignotte 精度判定 |
| **MTSHL Newton 迭代** | Wang.lean 627-786 | polynomial_factorize_wang.hh 866-932 | 完全 1:1：不变量初始化 → 单步推进(Leibniz+MDP+因子定理) → 终止 |
| **MDP 存在性** | Wang.lean 318-342 | Bézout 构造的数学正确性 | 完全 1:1：归纳 Bézout 构造 |
| **Vandermonde 可逆性** | Wang.lean 810-840 | __si_vandermonde_solve 数学正确性 | 完全 1:1：det ≠ 0 → 唯一解 |

#### 规约层面的模块（验证后置条件，不建模内部算法）

| 模块 | 模型方式 | 不建模的 C++ 内容 | 理由 |
|------|---------|------------------|------|
| **EDF** | `edf_correct` 建模 Cantor-Zassenhaus 结构；`unconditional` 用 UFD | 随机数选取 | Lean 无原生随机；UFD 路径证明 EDF 存在性 |
| **LC 分配** | `lcDistribCore` 建模递归提取；后置条件规约 | GCD 提取 LC 公因子 | GCD 算法单独验证 |
| **试除** | `TrialDivResult` 后置条件 | Gosper's hack 子集枚举 + divmod | 控制流优化，不影响数学正确性 |
| **MvHenselOutput** | 输出后置条件 | MTSHL 内部提升过程 | 由 mtshl_step_invariant 间接覆盖 |
| **Recombination** | `recombine_correct` 用 UFD | Zassenhaus 子集穷举 | 控制流；因子恢复已 1:1 验证 |

#### 简化处理的 MDP 求解器变体

| 定理 | 对应 C++ | 简化方式 | 影响 |
|------|---------|---------|------|
| `sparse_int_correct` | `__mtshl_sparse_int` | 归结到 `mdp_exists` | 未建模 θ-array 求值、Vandermonde 恢复 |
| `multi_bdp_correct` | `__mtshl_multi_bdp` | 归结到 `mdp_exists` | 未建模二变量 Taylor 循环 |
| `wmds_correct` | `__mtshl_wmds` | 归结到 `mdp_exists` | 未建模递归 WMDS 结构 |
| `mdp_cascade_correct` | cascade fallback | 归结到 `mdp_exists` | 未建模级联控制流 |

**注**：三个 MDP 求解器（sparse_int / multi_bdp / wmds）在数学上等价——都是对互素因子做偏分式分解，只是计算效率不同（稀疏插值 vs Taylor 循环 vs 递归）。`mdp_exists` 证明了偏分式分解的存在性和正确性，三种变体的区别纯属实现优化。但这意味着 **Lean 未验证这三种求解器各自特有的算法逻辑**（θ-array、Vandermonde Gauss-Jordan、递归 WMDS）。

### 未覆盖的内容（按设计排除）

- **GCD 算法**：Euclid/HGCD 的 L2 验证（独立模块，后续处理）
- **多项式算术**：pair_vec_div / multiplies 的正确性
- **L1 实现模型**：uint64 语义、数组越界、move 语义
- **素数选取**：`__select_prime` 控制流
- **随机性**：EDF 的随机数生成、eval point 选取

## 关键决策

1. **beta-redex 处理**：`∃!` 展开的 lambda 需要 `show` 强制归约，不能直接 `rw`
2. **mul_nonsing_inv 参数**：必须显式传入矩阵和 `IsUnit det` 证明，不能用 `_` 推断

## 度量

- 耗时：~0.5 小时（诊断 beta-redex + 修复）
- 迭代：2 轮编译-修复
- Lean 修改行数：~10 行
- 放弃的方案：无

## 涉及文件

- `CLPoly/Algorithm/Wang.lean`：消除最后 1 sorry
