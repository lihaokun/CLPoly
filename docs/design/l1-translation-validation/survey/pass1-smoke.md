# Pass 1 parse 全量烟测 (65 函数)

- OK (0 Unknown): **63**
- WARN (有 Unknown 节点): **2**
- FAIL (parse 崩溃): **0**

## 全局 Unknown 节点直方图

### Unknown Stmt kinds
| kind | count |
|---|---|
| (none) | 0 |

### Unknown Expr kinds
| kind | count |
|---|---|
| (none) | 0 |

### Unknown Type qualTypes (top 20)
| qualType | count |
|---|---|
| `std::_Bit_reference` | 6 |

## Ref 参数自动检测 vs TRANSLATION_SCOPE_OUTPUT_PARAMS

⚠️  有 **7** 个函数不匹配：

| 函数 | auto | configured |
|---|---|---|
| `__extract_monomial_content` | [1] | [] |
| `__mtshl_multi_bdp` | [5] | [3, 8] |
| `__mtshl_sparse_int` | [6] | [3, 9] |
| `__mtshl_step_j` | [1] | [3, 5] |
| `__mtshl_wmds` | [5] | [3, 8] |
| `__upoly_divmod` | [0, 1] | [] |
| `__upoly_random` | [2] | [] |

## 函数明细

### FAIL 列表

(无)

### WARN 列表（按 total_unknowns 降序）

| 函数 | total Unknown | stmt kinds | expr kinds | unknown type refs |
|---|---|---|---|---|
| `__vanhoeij_recombine` | 3 | - | - | 3 |
| `__wang_core` | 3 | - | - | 3 |

### OK 列表

- `__assign_partial_zp`
- `__binomial`
- `__build_cld_matrix`
- `__cld_polys`
- `__ddf_Zp`
- `__edf_Zp`
- `__extract_candidates`
- `__extract_monomial_content`
- `__extract_pth_root`
- `__factor_Zp`
- `__factor_multivar`
- `__factor_recombine`
- `__factor_squarefree_primitive_ZZ`
- `__hensel_extract_factors`
- `__hensel_lift`
- `__hensel_lift_linear_recursive`
- `__hensel_lift_recursive`
- `__hensel_step`
- `__hensel_step_linear`
- `__hensel_tree_build`
- `__hensel_tree_build_recursive`
- `__heuristic_starting_precision`
- `__isqrt_ceil`
- `__lll_factorize`
- `__lll_reduce`
- `__make_zp`
- `__mignotte_bound`
- `__mtshl_coeff_bound`
- `__mtshl_lift`
- `__mtshl_multi_bdp`
- `__mtshl_sparse_int`
- `__mtshl_step_j`
- `__mtshl_wmds`
- `__mtshl_zp_univar_mdp`
- `__polynomial_to_zp`
- `__select_eval_point`
- `__select_prime`
- `__si_theta_array_eval`
- `__si_vandermonde_solve`
- `__squarefree_Zp`
- `__subset_product_mod`
- `__symmetric_mod`
- `__symmetric_mod_poly`
- `__taylor_coeff_zp`
- `__upoly_const_term`
- `__upoly_divmod`
- `__upoly_divmod_mod`
- `__upoly_make_monic`
- `__upoly_mod`
- `__upoly_mod_coeff`
- `__upoly_mul_mod`
- `__upoly_norm_l1`
- `__upoly_norm_l2_sq`
- `__upoly_powmod`
- `__upoly_primitive`
- `__upoly_random`
- `__upoly_subtract_one`
- `__upoly_subtract_x`
- `__upoly_symmetric_mod`
- `__upoly_to_poly`
- `__wang_leading_coeff`
- `__zassenhaus_recombine`
- `factorize`