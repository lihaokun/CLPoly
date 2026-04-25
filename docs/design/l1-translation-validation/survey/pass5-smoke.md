# Pass 1 + 2 + 3 + 4 + 5 全量烟测

- 目标函数数：**65**（factorize 3 实例 → 67 HIRs）
- OK: **67** / FAIL: **0**
- Gap 总数: **24**（B 策略残余，Pass 8 codegen 时输出 sorry）

## Gap 分类

| 类别 | 数量 |
|---|---|
| `op_miss` | 12 |
| `method_miss` | 8 |
| `ctor_miss` | 4 |

## Top 20 Gap 详情

| 次数 | 类别 | 详情 |
|---|---|---|
| 7 | op_miss | `('None', '+')` |
| 4 | method_miss | `('BaseType.BOOL', 'operator bool')` |
| 3 | method_miss | `('SparsePolyZZ', 'erase')` |
| 2 | op_miss | `('None', '*')` |
| 2 | op_miss | `('None', '-')` |
| 2 | ctor_miss | `('default_init_const std::allocator<char>', 0)` |
| 2 | ctor_miss | `('construct_const std::string', 2)` |
| 1 | op_miss | `('None', '<')` |
| 1 | method_miss | `('None', 'push_back')` |

## Per-function Gap 数

| 函数 | gap |
|---|---|
| `__mtshl_lift` | 4 |
| `factorize_upoly` | 4 |
| `__lll_reduce` | 3 |
| `__vanhoeij_recombine` | 3 |
| `__wang_core` | 3 |
| `__mtshl_wmds` | 2 |
| `__hensel_step_linear` | 1 |
| `__mtshl_sparse_int` | 1 |
| `__upoly_divmod_mod` | 1 |
| `__upoly_mod_coeff` | 1 |
| `__zassenhaus_recombine` | 1 |
