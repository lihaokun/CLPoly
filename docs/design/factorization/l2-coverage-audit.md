# L2 覆盖审计：C++ 因式分解模块 vs Lean 证明

> 日期：2026-03-25
> 审计范围：5 个 C++ 文件（5,998 行）vs Lean L2 模型（6,276 行）

---

## 1. polynomial_factorize_zp.hh（396 行，14 函数）

| C++ 函数 | 行 | 用途 | L2 状态 | Lean 位置 |
|----------|-----|------|---------|----------|
| `__squarefree_Zp` | 123-180 | Yun SQF 分解 | **已覆盖** | SquarefreeZp.lean: yunLoop + sqf_correct |
| `__ddf_Zp` | 247-292 | 等度分解 DDF | **已覆盖** | DDF.lean: ddfLoop + ddf_correct |
| `__edf_Zp` | 295-354 | Cantor-Zassenhaus EDF | **已覆盖** | EDF.lean: edf + edf_correct + exists_nonQR_poly |
| `__factor_Zp` | 357-391 | SQF→DDF→EDF 管线 | **已覆盖** | FactorZp.lean: factor_Zp_correct |
| `__upoly_make_monic` | 27-36 | 首一化 | Mathlib 代替 | normalize / Monic typeclass |
| `__upoly_mod` | 39-46 | 多项式取模 | Mathlib 代替 | modByMonic |
| `__upoly_divmod` | 49-56 | 多项式除法 | Mathlib 代替 | divByMonic |
| `__upoly_powmod` | 59-85 | 模幂运算 | Mathlib 代替 | DDF 不变量隐含 |
| `__extract_pth_root` | 109-120 | p 次根提取 | Mathlib 代替 | contract p（SQF 中使用） |
| `__make_zp` | 20 | Zp 构造辅助 | 不需要 | 内部工具 |
| `__upoly_random` | 88-102 | 随机多项式 | 不需要 | Classical.choice 替代 |
| `__upoly_subtract_x` | 183-214 | h - X 计算 | 不需要 | 算术辅助 |
| `__upoly_subtract_one` | 217-244 | h - 1 计算 | 不需要 | 算术辅助 |

## 2. polynomial_factorize_univar.hh（1661 行，~30 函数）

### 已覆盖

| C++ 函数 | 行 | 用途 | Lean 位置 |
|----------|-----|------|----------|
| `__hensel_step` | 404-497 | 二次 Hensel 单步 | Hensel.lean: hensel_step + hensel_step_with_degree (H1-H6) |
| `__hensel_lift` | 537-599 | 二次 Hensel 提升到 Mignotte 精度 | Hensel.lean: hensel_two_factor + hensel_multifactor |
| `__mignotte_bound` | 176-184 | Mignotte 系数界 | Recombine.lean: mignotte_bound_l2 + mignotte_bound |
| `__zassenhaus_recombine` | 751-882 | Zassenhaus 子集枚举重组（r≤10） | Recombine.lean: ZassenhausInvariant + init/extract/terminate |
| `__factor_recombine` | 1349-1359 | 路由 Zassenhaus/van Hoeij | Recombine.lean: recombine_correct |
| `__lll_factorize` | 1484-1528 | van Hoeij 管线入口 | Pipeline/FactorZZ.lean |
| `__factor_squarefree_primitive_ZZ` | 1538-1558 | 素数选取→Hensel→重组 | Pipeline/FactorZZ.lean: factor_ZZ_correct |
| `factorize` (upolynomial) | 1564-1657 | 顶层单变量因式分解 | Pipeline/FactorZZInstantiate.lean |

### 部分覆盖

| C++ 函数 | 行 | 用途 | L2 状态 | 差距 |
|----------|-----|------|---------|------|
| `__vanhoeij_recombine` | 1188-1341 | van Hoeij LLL 主循环 | **黑盒** | 仅输入/输出条件，内部 CLD+格+LLL+候选提取未建模 |
| `__hensel_tree_build` | 334-401 | Hensel 二叉树构造 | **隐含** | hensel_multifactor 概念覆盖，树结构未单独建模 |
| `__select_prime` | 1388-1470 | 素数选取枚举 | **隐含** | 存在性假设，枚举逻辑未建模 |

### 未覆盖（van Hoeij 子函数）

| C++ 函数 | 行 | 用途 | 说明 |
|----------|-----|------|------|
| `__cld_polys` | 900-927 | CLD 系数对数导数计算 | van Hoeij 核心：C_i = (f/h_i)·h_i' mod m |
| `__build_cld_matrix` | 937-983 | 格矩阵构造（螺旋 CLD 列插入） | van Hoeij 核心：格基扩展 |
| `__lll_reduce` | 992-1129 | LLL 格基约化（Cohen §2.6.3） | 可作可信原语（DJTY Isabelle 已验证） |
| `__extract_candidates` | 1138-1180 | 从短向量提取因子子集候选 | van Hoeij 核心：列等价分组 |

### 不需要（未使用/辅助/优化）

| C++ 函数 | 行 | 用途 | 说明 |
|----------|-----|------|------|
| `__hensel_step_linear` | 636-679 | 线性 Hensel 步 | **未被调用**：预留给 P1b 完整版 |
| `__hensel_lift_linear_recursive` | 682-694 | 线性 Hensel 递归提升 | **未被调用**：同上 |
| `__symmetric_mod` | 94-103 | 标量对称模约化 | 实现辅助 |
| `__upoly_symmetric_mod` | 106-119 | 多项式对称模约化 | 实现辅助 |
| `__upoly_norm_l2_sq` | 126-132 | L2 范数平方 | Mignotte 内部 |
| `__upoly_norm_l1` | 701-711 | L1 范数 | 内部辅助 |
| `__binomial` | 139-151 | 二项式系数 | Mignotte 内部 |
| `__isqrt_ceil` | 154-173 | 整数平方根上取整 | Mignotte 内部 |
| `__upoly_mod_coeff` | 191-206 | 系数模约化 | Hensel 内部 |
| `__upoly_divmod_mod` | 209-304 | 模长除法 | Hensel 内部 |
| `__upoly_mul_mod` | 323-331 | 模乘法 | Hensel 内部 |
| `__hensel_extract_factors` | 500-515 | 从树提取叶因子 | 树遍历 |
| `__upoly_primitive` | 714-722 | 提取 content | 辅助 |
| `__subset_product_mod` | 725-740 | 子集乘积 mod m | Zassenhaus 内部 |
| `__upoly_const_term` | 743-748 | 常数项提取 | 辅助 |
| `__heuristic_starting_precision` | 607-631 | 启发式精度计算 | 性能优化 |
| `__upoly_to_poly` | 1366-1374 | 类型转换 | 实现细节 |
| `print/reset_factorize_profile` | 45-75 | 性能 profile | 与正确性无关 |

## 3. polynomial_factorize_wang.hh（2560 行，28 函数）

### 已覆盖

| C++ 函数 | 行 | 用途 | Lean 位置 |
|----------|-----|------|----------|
| `__si_vandermonde_solve` | 66-123 | Vandermonde 求解（Gauss-Jordan） | Wang.lean: vandermonde_solve_unique |
| `__si_theta_array_eval` | 133-202 | θ-array 批量求值 | Wang.lean: evalAtBetaPow + evalAtBetaPow_linearTerm |
| `__mtshl_zp_univar_mdp` | 244-301 | 单变量 Bezout MDP | Wang.lean: mdp_exists |
| `__mtshl_step_j` | 763-936 | MTSHL 单步（Taylor 循环 + MDP 级联） | Wang.lean: mtshl_step_invariant |
| `__mtshl_lift` | 1017-1183 | 多变量 Hensel 提升主循环 | Wang.lean: mtshl_invariant_init/terminates |
| `__wang_leading_coeff` | 1325-1528 | LC 分配（valuation 提取） | Wang.lean: lcDistribCore + lc_distrib_prod_correct |
| `__wang_core` | 2228-2439 | Wang 算法主循环 | Wang.lean: wang_correct |
| `__factor_multivar` | 2444-2556 | 多变量因式分解入口 | FactorMv.lean: mv_factor_correct |
| `__multivar_diophantine` | 1641-1830 | 递归多变量丢番图求解 | Wang.lean: mdp_exists + partialEval_linearTerm |
| `__hensel_lc_correct` | 1834-1856 | LC 校正 | Wang.lean: mtshl_step_invariant 内 |
| `__hensel_lift_one_var` | 1864-1984 | 单变量 Hensel 提升 | Wang.lean: mtshl_invariant_terminates |
| `__multivar_hensel_lift` | 1990-2136 | 多变量 Hensel 入口 | Hensel.lean: hensel_multifactor 概念 |

### 部分覆盖（算法模型已建，但不如 C++ 细）

| C++ 函数 | 行 | 用途 | Lean 位置 | 差距 |
|----------|-----|------|----------|------|
| `__mtshl_sparse_int` | 449-605 | θ-array 稀疏插值 MDP | Wang.lean: sparse_int_correct | evalAtBetaPow_linearTerm 建模了数学链，但 θ-array 内部优化和 Vandermonde 恢复步骤未逐步建模 |
| `__mtshl_multi_bdp` | 307-438 | 二变量 Taylor MDP | Wang.lean: multi_bdp_correct | MultiBdpInvariant 建模了循环不变量，但 Taylor 系数提取步骤未逐步建模 |
| `__mtshl_wmds` | 615-758 | 递归 WMDS MDP | Wang.lean: wmds_correct | 复用 MultiBdpInvariant，递归降维结构文档化但未逐步建模 |

### 隐含覆盖（属性已嵌入其他定理，未单独建模）

| C++ 函数 | 行 | 用途 |
|----------|-----|------|
| `__taylor_coeff` / `__taylor_coeff_zp` | 25-53 / 211-239 | Taylor 系数提取 |
| `__mtshl_coeff_bound` | 941-954 | Mignotte 系数界 |
| `__select_eval_point` | 1193-1311 | 求值点选取 |
| `__polynomial_to_zp` / `__assign_partial_zp` | 958-991 | 模约化/偏求值 |
| `__symmetric_mod_poly` | 995-1009 | 对称模约化 |
| `__pseudo_remainder_x1` | 1540-1632 | 伪余式 |
| `__extract_monomial_content` | 2148-2223 | 单项式 content 提取 |
| `poly_convert` | 多处 | 类型转换 |
| `__wang_lc_result` | 1315-1321 | LC 结果结构体 |

## 4. polynomial_factorize.hh（235 行，~6 函数）

| C++ 函数 | 行 | 用途 | L2 状态 | Lean 位置 |
|----------|-----|------|---------|----------|
| `factorize` (ZZ lex) | 22-128 | 顶层入口：content + SQF + 因式分解 | **已覆盖** | FactorZZ.lean + FactorZZInstantiate.lean |
| `factorize` (ZZ generic) | 134-154 | 任意序转 lex 再因式分解 | **已覆盖** | FactorZZInstantiate.lean |
| `factorize` (QQ) | 160-232 | QQ→ZZ 转换再因式分解 | **未覆盖** | LCD 转换逻辑未建模 |
| `__factor_multivar` dispatch | 47 | 多变量 dispatch | **已覆盖** | FactorMv.lean |

## 5. polynomial_gcd.hh（1146 行）— 可信原语

| 状态 | 说明 |
|------|------|
| **整体未覆盖** | Euclidean GCD, HGCD, 模 GCD, CRT 均作为 Mathlib `GCDMonoid.gcd` 可信原语 |

---

## 覆盖缺失总结

### 活跃代码路径中的缺失（必须补充）

| 编号 | 缺失 | C++ 函数 | 行数 | 难度 |
|------|------|---------|------|------|
| **G1** | van Hoeij CLD 多项式 | `__cld_polys` | ~30 | 中 |
| **G2** | van Hoeij 格构造 | `__build_cld_matrix` | ~50 | 中 |
| **G3** | van Hoeij 候选提取 | `__extract_candidates` | ~40 | 中 |
| **G4** | van Hoeij 主循环（回溯逻辑） | `__vanhoeij_recombine` 内部 | ~150 | 高 |
| **G5** | LLL 格基约化 | `__lll_reduce` | ~140 | 高（可作可信原语） |
| **G6** | 素数选取枚举 | `__select_prime` | ~80 | 中 |
| **G7** | QQ[x] 因式分解 | `factorize(QQ)` | ~70 | 低 |

### 按设计排除（可信原语/未使用代码）

| 编号 | 排除项 | 理由 |
|------|--------|------|
| E1 | GCD 算法（1146 行） | 独立模块，Mathlib 可信原语 |
| E2 | 多项式算术辅助 | Mathlib Polynomial 操作 |
| E3 | 线性 Hensel（P1b） | 代码存在但**未被调用** |
| E4 | Profile 函数 | 与正确性无关 |
| E5 | 类型转换函数 | 实现细节 |
