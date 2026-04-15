# L2 覆盖最终审计：C++ 因式分解模块 vs Lean 证明

> 日期：2026-03-26（最终版）
> 审计范围：4 个 C++ 文件（4,852 行）vs Lean L2 模型（6,860 行）

---

## 审计结果总览

| 文件 | C++ 行 | COVERED | MATHLIB | NOT_NEEDED | MISSING |
|------|--------|---------|---------|------------|---------|
| polynomial_factorize_zp.hh | 396 | 4 | 7 | 2 | 0 |
| polynomial_factorize_univar.hh | 1661 | 35 | 1 | 7 | 2* |
| polynomial_factorize_wang.hh | 2560 | 24 | 4 | 3 | 0 |
| polynomial_factorize.hh | 235 | 4 | 0 | 0 | 0 |
| **总计** | **4,852** | **67** | **12** | **12** | **2*** |

*2 个 MISSING 是**线性 Hensel**（`__hensel_step_linear`/`__hensel_lift_linear_recursive`）——C++ 管线中**未被调用**的预留代码。

**活跃代码路径 MISSING = 0。**

## 分类说明

### COVERED（67 个函数）— 有对应 Lean 定义/定理

核心算法全部覆盖：SQF(Yun)、DDF(Shoup)、EDF(Cantor-Zassenhaus)、Hensel(二次+多因子)、Mignotte、Zassenhaus 重组、Van Hoeij LLL(CLD+格+LLL+候选提取+回溯)、MTSHL Newton、稀疏插值、Multi-BDP/WMDS、MDP Bézout、Vandermonde、Wang 组合、素数选取、QQ 因式分解、辅助函数(Taylor系数、系数界、求值点、伪余式、单项式content、LC结果结构体)。

### MATHLIB（12 个函数）— 由 Mathlib 标准 API 覆盖

多项式算术操作：`__upoly_mod`(modByMonic)、`__upoly_divmod`(divByMonic)、`__upoly_make_monic`(normalize)、`__upoly_powmod`(DDF 内隐含)、`__extract_pth_root`(contract p)、`__upoly_subtract_x/one`(多项式减法)、`__polynomial_to_zp`(MvPolynomial.map)、`__assign_partial_zp`(partialEval + map)、`__symmetric_mod_poly`(ZMod.valMinAbs)、`poly_convert`(Polynomial.map)、`__select_prime`(GoodPrime 条件)。

### NOT_NEEDED（12 个函数）— 未调用/优化/配置

线性 Hensel(2, 未调用)、启发式精度(1)、类型转换(3)、Profile(2)、全局配置(1)、阈值常数(1)、结构体字段(2)。

## Lean 代码统计

| 文件 | 行数 | sorry |
|------|------|-------|
| Algorithm/SquarefreeZp.lean | 1,501 | 0 |
| Algorithm/Wang.lean | 1,273 | 0 |
| Algorithm/Hensel.lean | 887 | 0 |
| Algorithm/Recombine.lean | 1,066 | 0 |
| Algorithm/DDF.lean | 435 | 0 |
| Algorithm/EDF.lean | 431 | 0 |
| Pipeline/*.lean | 554 | 0 |
| Math/*.lean | 413 | 0 |
| Spec.lean | 300 | 0 |
| **总计** | **6,860** | **0** |
