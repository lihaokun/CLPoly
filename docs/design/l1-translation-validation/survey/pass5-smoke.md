# Pass 1 + 2 + 3 + 4 + 5 全量烟测

- 目标函数数：**66**（factorize 3 实例 → 68 HIRs）
- OK: **68** / FAIL: **0**
- Gap 总数: **15**（B 策略残余，Pass 8 codegen 时输出 sorry）

## Gap 分类

| 类别 | 数量 |
|---|---|
| `method_miss` | 12 |
| `ctor_miss` | 3 |

## Top 20 Gap 详情

| 次数 | 类别 | 详情 |
|---|---|---|
| 6 | method_miss | `('polynomial_<clpoly::ZZ, clpoly::lex_<clpoly::less>>', 'front')` |
| 3 | method_miss | `('polynomial_<clpoly::ZZ, clpoly::lex_<clpoly::less>>', 'empty')` |
| 2 | ctor_miss | `('construct_const value_type', 1)` |
| 2 | method_miss | `('SparsePoly_clpoly::Zp', 'empty')` |
| 1 | method_miss | `('SparsePoly_clpoly::ZZ', 'empty')` |
| 1 | ctor_miss | `('construct_LLLMatrix', 2)` |

## Per-function Gap 数

| 函数 | gap |
|---|---|
| `__select_eval_point` | 5 |
| `__wang_leading_coeff` | 5 |
| `__factor_squarefree_primitive_ZZ` | 1 |
| `__select_prime` | 1 |
| `__squarefree_Zp` | 1 |
| `__vanhoeij_recombine` | 1 |
| `__zassenhaus_recombine` | 1 |
