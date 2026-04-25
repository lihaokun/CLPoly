# Pass 1 + 2 + 3 + 4 + 5 全量烟测

- 目标函数数：**65**（factorize 3 实例 → 67 HIRs）
- OK: **67** / FAIL: **0**
- Gap 总数: **81**（B 策略残余，Pass 8 codegen 时输出 sorry）

## Gap 分类

| 类别 | 数量 |
|---|---|
| `op_miss` | 56 |
| `method_miss` | 21 |
| `ctor_miss` | 4 |

## Top 20 Gap 详情

| 次数 | 类别 | 详情 |
|---|---|---|
| 21 | op_miss | `('None', '*')` |
| 10 | op_miss | `('None', '+')` |
| 10 | method_miss | `('ZZ', 'operator bool')` |
| 7 | op_miss | `('None', '-')` |
| 5 | op_miss | `('None', '/')` |
| 4 | method_miss | `('BaseType.BOOL', 'operator bool')` |
| 3 | method_miss | `('SparsePolyZZ', 'erase')` |
| 3 | method_miss | `('None', 'sizeinbase')` |
| 3 | op_miss | `('None', '%')` |
| 3 | op_miss | `('None', '<<')` |
| 2 | op_miss | `('None', '<')` |
| 2 | op_miss | `('None', '-=')` |
| 2 | op_miss | `('None', '==')` |
| 2 | ctor_miss | `('default_init_const std::allocator<char>', 0)` |
| 2 | ctor_miss | `('construct_const std::string', 2)` |
| 1 | method_miss | `('None', 'push_back')` |
| 1 | op_miss | `('None', '!=')` |

## Per-function Gap 数

| 函数 | gap |
|---|---|
| `__lll_reduce` | 13 |
| `__vanhoeij_recombine` | 10 |
| `__upoly_divmod_mod` | 9 |
| `__mtshl_lift` | 6 |
| `__si_vandermonde_solve` | 5 |
| `__si_theta_array_eval` | 4 |
| `__wang_leading_coeff` | 4 |
| `factorize_upoly` | 4 |
| `__hensel_step_linear` | 3 |
| `__wang_core` | 3 |
| `__edf_Zp` | 2 |
| `__heuristic_starting_precision` | 2 |
| `__isqrt_ceil` | 2 |
| `__mtshl_sparse_int` | 2 |
| `__mtshl_wmds` | 2 |
| `__upoly_mod_coeff` | 2 |
| `__hensel_lift` | 1 |
| `__hensel_step` | 1 |
| `__select_prime` | 1 |
| `__upoly_powmod` | 1 |
| `__upoly_subtract_one` | 1 |
| `__upoly_subtract_x` | 1 |
| `__upoly_symmetric_mod` | 1 |
| `__zassenhaus_recombine` | 1 |
