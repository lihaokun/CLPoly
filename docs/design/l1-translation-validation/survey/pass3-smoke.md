# Pass 1 + Pass 2 + Pass 3 全量烟测

- 目标函数数（TRANSLATION_SCOPE）：**65**（factorize 展开 3 实例 → 67 HIRs）
- OK: **67**
- FAIL: **0**

- 有 lambda 的宿主 HIR 数: **15**
- 提升的 lambda 总数: **26**

## 宿主函数 lambda 提升详情

| 宿主函数 | 提升数 | 生成的 aux_lambdas |
|---|---|---|
| `__build_cld_matrix` | 1 | `_lambda___build_cld_matrix_1`(2) |
| `__factor_Zp` | 1 | `_lambda___factor_Zp_1`(2) |
| `__factor_multivar` | 1 | `_lambda___factor_multivar_1`(2) |
| `__lll_reduce` | 5 | `_lambda___lll_reduce_1`(2), `_lambda___lll_reduce_2`(1), `_lambda___lll_reduce_3`(6), `_lambda___lll_reduce_4`(4), `_lambda___lll_reduce_5`(4) |
| `__mtshl_multi_bdp` | 1 | `_lambda___mtshl_multi_bdp_1`(4) |
| `__mtshl_sparse_int` | 2 | `_lambda___mtshl_sparse_int_1`(3), `_lambda___mtshl_sparse_int_2`(4) |
| `__mtshl_step_j` | 2 | `_lambda___mtshl_step_j_1`(5), `_lambda___mtshl_step_j_2`(2) |
| `__mtshl_wmds` | 1 | `_lambda___mtshl_wmds_1`(4) |
| `__select_prime` | 1 | `_lambda___select_prime_1`(2) |
| `__vanhoeij_recombine` | 2 | `_lambda___vanhoeij_recombine_1`(2), `_lambda___vanhoeij_recombine_2`(2) |
| `__wang_core` | 3 | `_lambda___wang_core_1`(1), `_lambda___wang_core_2`(2), `_lambda___wang_core_3`(2) |
| `__wang_leading_coeff` | 2 | `_lambda___wang_leading_coeff_1`(2), `_lambda___wang_leading_coeff_2`(2) |
| `__zassenhaus_recombine` | 2 | `_lambda___zassenhaus_recombine_1`(2), `_lambda___zassenhaus_recombine_2`(2) |
| `factorize_lex` | 1 | `_lambda_factorize_1`(2) |
| `factorize_upoly` | 1 | `_lambda_factorize_1`(2) |

_(params count 包含 captures + lambda 自身参数)_

## 与 `lambdas.md` 预期值核对

- 本轮测得 lambda 总数: **26**
- `lambdas.md` in-scope 预期: **26**

✅ 与预期值完全一致
