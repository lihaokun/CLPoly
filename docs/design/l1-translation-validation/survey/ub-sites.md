# UB 站点全量枚举

> 参考 `ub-catalog.md` 的 UB-1 到 UB-8 分类定义。
> 本文件是 65 函数的**实测**UB 站点清单。

## 全局统计

| UB 类型 | 站点数 |
|---|---|
| UB-1 Div0 | 31 |
| UB-2 OOB | 468 |
| UB-3 EmptyContainer | 42 |
| UB-6 SignedOverflow | 113 |
| UB-7 UnsignedToSigned | 8 |
| UB-8 Assert | 23 |
| **合计** | **685** |

## 按函数汇总（UB 站点数）

| 函数 | UB-1 Div | UB-2 OOB | UB-3 Empty | UB-4 Shift | UB-6 Signed | UB-7 U→S | UB-8 Assert | 合计 |
|---|---|---|---|---|---|---|---|---|
| `__lll_reduce` | 3 | 108 | 0 | 0 | 21 | 0 | 0 | **132** |
| `__wang_core` | 1 | 35 | 2 | 0 | 15 | 0 | 0 | **53** |
| `__wang_leading_coeff` | 7 | 31 | 8 | 0 | 2 | 3 | 2 | **53** |
| `__si_vandermonde_solve` | 0 | 38 | 0 | 0 | 4 | 0 | 1 | **43** |
| `__mtshl_sparse_int` | 0 | 38 | 0 | 0 | 0 | 0 | 0 | **38** |
| `__mtshl_wmds` | 0 | 31 | 2 | 0 | 0 | 0 | 0 | **33** |
| `__zassenhaus_recombine` | 1 | 14 | 2 | 0 | 14 | 0 | 0 | **31** |
| `__mtshl_lift` | 0 | 23 | 0 | 0 | 7 | 0 | 0 | **30** |
| `__mtshl_step_j` | 0 | 25 | 0 | 0 | 0 | 0 | 0 | **25** |
| `__vanhoeij_recombine` | 1 | 9 | 1 | 0 | 11 | 2 | 0 | **24** |
| `__hensel_tree_build_recursive` | 1 | 14 | 0 | 0 | 7 | 0 | 0 | **22** |
| `__mtshl_multi_bdp` | 0 | 22 | 0 | 0 | 0 | 0 | 0 | **22** |
| `__si_theta_array_eval` | 0 | 19 | 0 | 0 | 0 | 0 | 1 | **20** |
| `__extract_candidates` | 0 | 15 | 0 | 0 | 4 | 0 | 0 | **19** |
| `__mtshl_zp_univar_mdp` | 0 | 14 | 1 | 0 | 0 | 0 | 1 | **16** |
| `__build_cld_matrix` | 3 | 3 | 1 | 0 | 8 | 0 | 0 | **15** |
| `__hensel_lift` | 0 | 6 | 4 | 0 | 0 | 0 | 3 | **13** |
| `__select_eval_point` | 1 | 4 | 3 | 0 | 4 | 0 | 1 | **13** |
| `__upoly_divmod_mod` | 4 | 0 | 2 | 0 | 3 | 0 | 2 | **11** |
| `__hensel_lift_recursive` | 0 | 7 | 0 | 0 | 2 | 0 | 0 | **9** |
| `__hensel_lift_linear_recursive` | 0 | 7 | 0 | 0 | 0 | 0 | 0 | **7** |
| `__extract_pth_root` | 2 | 0 | 1 | 0 | 0 | 1 | 1 | **5** |
| `__binomial` | 0 | 0 | 0 | 0 | 4 | 0 | 0 | **4** |
| `__heuristic_starting_precision` | 0 | 0 | 1 | 0 | 2 | 1 | 0 | **4** |
| `__factor_multivar` | 0 | 0 | 3 | 0 | 0 | 0 | 1 | **4** |
| `__upoly_powmod` | 2 | 0 | 1 | 0 | 0 | 0 | 1 | **4** |
| `__hensel_extract_factors` | 0 | 1 | 0 | 0 | 2 | 0 | 0 | **3** |
| `__select_prime` | 0 | 0 | 1 | 0 | 1 | 0 | 1 | **3** |
| `__ddf_Zp` | 0 | 0 | 1 | 0 | 0 | 1 | 1 | **3** |
| `__factor_Zp` | 0 | 0 | 2 | 0 | 0 | 0 | 1 | **3** |
| `__isqrt_ceil` | 3 | 0 | 0 | 0 | 0 | 0 | 0 | **3** |
| `__assign_partial_zp` | 0 | 2 | 0 | 0 | 0 | 0 | 0 | **2** |
| `__extract_monomial_content` | 0 | 1 | 0 | 0 | 1 | 0 | 0 | **2** |
| `__edf_Zp` | 1 | 0 | 1 | 0 | 0 | 0 | 0 | **2** |
| `__squarefree_Zp` | 0 | 0 | 1 | 0 | 0 | 0 | 1 | **2** |
| `__upoly_const_term` | 0 | 0 | 2 | 0 | 0 | 0 | 0 | **2** |
| `__upoly_make_monic` | 0 | 0 | 1 | 0 | 0 | 0 | 1 | **2** |
| `__mignotte_bound` | 1 | 0 | 0 | 0 | 0 | 0 | 1 | **2** |
| `__subset_product_mod` | 0 | 1 | 0 | 0 | 0 | 0 | 0 | **1** |
| `__upoly_random` | 0 | 0 | 0 | 0 | 1 | 0 | 0 | **1** |
| `__upoly_primitive` | 0 | 0 | 1 | 0 | 0 | 0 | 0 | **1** |
| `__cld_polys` | 0 | 0 | 0 | 0 | 0 | 0 | 1 | **1** |
| `__factor_squarefree_primitive_ZZ` | 0 | 0 | 0 | 0 | 0 | 0 | 1 | **1** |
| `__hensel_tree_build` | 0 | 0 | 0 | 0 | 0 | 0 | 1 | **1** |

## UB-1 Div0 站点样本（前 10）

总计 31，分布于 14 个函数。

| 宿主 | 行 | 细节 |
|---|---|---|
| `__build_cld_matrix` | None | op=`%` / type=`int` |
| `__build_cld_matrix` | None | op=`/` / type=`int` |
| `__build_cld_matrix` | None | op=`/` / type=`int` |
| `__edf_Zp` | None | op=`operator/` / type=`ZZ` |
| `__extract_pth_root` | None | op=`%` / type=`uint64_t` |
| `__extract_pth_root` | None | op=`/` / type=`uint64_t` |
| `__hensel_tree_build_recursive` | None | op=`/` / type=`int` |
| `__isqrt_ceil` | None | op=`/` / type=`size_t` |
| `__isqrt_ceil` | None | op=`operator/` / type=`ZZ` |
| `__isqrt_ceil` | None | op=`operator/` / type=`ZZ` |

## UB-2 OOB 站点样本（前 10）

总计 468，分布于 24 个函数。

| 宿主 | 行 | 细节 |
|---|---|---|
| `__assign_partial_zp` | None | op=`operator[]` / type=`const value_type` |
| `__assign_partial_zp` | None | op=`operator[]` / type=`const value_type` |
| `__build_cld_matrix` | None | op=`operator[]` / type=`value_type` |
| `__build_cld_matrix` | None | op=`operator[]` / type=`const value_type` |
| `__build_cld_matrix` | None | op=`operator[]` / type=`value_type` |
| `__extract_candidates` | None | op=`operator[]` / type=`value_type` |
| `__extract_candidates` | None | op=`operator[]` / type=`value_type` |
| `__extract_candidates` | None | op=`operator[]` / type=`const value_type` |
| `__extract_candidates` | None | op=`operator[]` / type=`const value_type` |
| `__extract_candidates` | None | op=`operator[]` / type=`const value_type` |

## UB-3 EmptyContainer 站点样本（前 10）

总计 42，分布于 22 个函数。

| 宿主 | 行 | 细节 |
|---|---|---|
| `__build_cld_matrix` | None | method=`front` / type=`const std::pair<umonomial, ZZ>` |
| `__ddf_Zp` | None | method=`front` / type=`const std::pair<umonomial, Zp>` |
| `__edf_Zp` | None | method=`front` / type=`const std::pair<umonomial, Zp>` |
| `__extract_pth_root` | None | method=`front` / type=`const std::pair<umonomial, Zp>` |
| `__factor_Zp` | None | method=`front` / type=`std::pair<umonomial, Zp>` |
| `__factor_Zp` | None | method=`front` / type=`std::pair<umonomial, Zp>` |
| `__factor_multivar` | None | method=`front` / type=`std::pair<basic_monomial<lex_<less>>, ZZ>` |
| `__factor_multivar` | None | method=`front` / type=`std::pair<basic_monomial<lex_<less>>, ZZ>` |
| `__factor_multivar` | None | method=`front` / type=`std::pair<basic_monomial<lex_<less>>, ZZ>` |
| `__hensel_lift` | None | method=`front` / type=`const std::pair<umonomial, ZZ>` |

## UB-6 SignedOverflow 站点样本（前 10）

总计 113，分布于 19 个函数。

| 宿主 | 行 | 细节 |
|---|---|---|
| `__binomial` | None | op=`-` / type=`int64_t` |
| `__binomial` | None | op=`-` / type=`int64_t` |
| `__binomial` | None | op=`-` / type=`int64_t` |
| `__binomial` | None | op=`+` / type=`int64_t` |
| `__build_cld_matrix` | None | op=`+` / type=`int` |
| `__build_cld_matrix` | None | op=`-` / type=`int` |
| `__build_cld_matrix` | None | op=`-` / type=`int` |
| `__build_cld_matrix` | None | op=`-` / type=`int` |
| `__build_cld_matrix` | None | op=`+` / type=`int` |
| `__build_cld_matrix` | None | op=`+` / type=`int` |

## UB-7 UnsignedToSigned 站点样本（前 10）

总计 8，分布于 5 个函数。

| 宿主 | 行 | 细节 |
|---|---|---|
| `__ddf_Zp` | None | source=`uint64_t` / target=`int64_t` |
| `__extract_pth_root` | None | source=`uint64_t` / target=`int64_t` |
| `__heuristic_starting_precision` | None | source=`size_t` / target=`int` |
| `__vanhoeij_recombine` | None | source=`size_t` / target=`int` |
| `__vanhoeij_recombine` | None | source=`size_t` / target=`int` |
| `__wang_leading_coeff` | None | source=`size_t` / target=`int` |
| `__wang_leading_coeff` | None | source=`size_t` / target=`int` |
| `__wang_leading_coeff` | None | source=`unsigned long` / target=`int` |

## UB-8 Assert 站点样本（前 10）

总计 23，分布于 19 个函数。

| 宿主 | 行 | 细节 |
|---|---|---|
| `__cld_polys` | None |  |
| `__ddf_Zp` | None |  |
| `__extract_pth_root` | None |  |
| `__factor_Zp` | None |  |
| `__factor_multivar` | None |  |
| `__factor_squarefree_primitive_ZZ` | None |  |
| `__hensel_lift` | None |  |
| `__hensel_lift` | None |  |
| `__hensel_lift` | None |  |
| `__hensel_tree_build` | None |  |
