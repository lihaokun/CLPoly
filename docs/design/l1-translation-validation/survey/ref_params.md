# 参数传递模式扫描

为 HIR `ref_elim` Pass 设计提供输入。

## 全局统计

- 总参数数: **192**
- `value`: 40 (21%)
- `const-ref`: 125 (65%)
- `ref`: 26 (14%)
- `pointer`: 1 (1%)
- `rvalue-ref`: 0 (0%)
- `unknown`: 0 (0%)

- 含 non-const ref 参数的函数（疑似输出参数）: **22**
- 含指针参数的函数: **1**
- 含 rvalue-ref 参数的函数: **0**

## 输出参数对比（AST 扫描 vs TRANSLATION_SCOPE_OUTPUT_PARAMS）

| 函数 | AST 扫出的 non-const ref 位置 | 已配置的 OUTPUT_PARAMS | 是否一致 |
|---|---|---|---|
| `__build_cld_matrix` | [0] | [0] | ✅ |
| `__edf_Zp` | [0, 3] | [0, 3] | ✅ |
| `__extract_monomial_content` | [1] | [] | ❌ |
| `__hensel_extract_factors` | [2] | [2] | ✅ |
| `__hensel_lift_linear_recursive` | [0] | [0] | ✅ |
| `__hensel_lift_recursive` | [0] | [0] | ✅ |
| `__hensel_step` | [0] | [0] | ✅ |
| `__hensel_step_linear` | [0] | [0] | ✅ |
| `__hensel_tree_build_recursive` | [0] | [0] | ✅ |
| `__lll_reduce` | [0, 1] | [0, 1] | ✅ |
| `__mtshl_multi_bdp` | [5] | [3, 8] | ❌ |
| `__mtshl_sparse_int` | [6] | [3, 9] | ❌ |
| `__mtshl_step_j` | [1] | [3, 5] | ❌ |
| `__mtshl_wmds` | [5] | [3, 8] | ❌ |
| `__mtshl_zp_univar_mdp` | [2] | [2] | ✅ |
| `__si_theta_array_eval` | [5] | [6] | ❌ |
| `__si_vandermonde_solve` | [2] | [2] | ✅ |
| `__upoly_divmod` | [0, 1] | [] | ❌ |
| `__upoly_divmod_mod` | [0, 1] | [0, 1] | ✅ |
| `__upoly_make_monic` | [0] | [0] | ✅ |
| `__upoly_mod_coeff` | [0] | [0] | ✅ |
| `__upoly_random` | [2] | [] | ❌ |

## 按函数参数详情

### `__assign_partial_zp` (4 params)
- [0] `f` : `const polynomial_<ZZ, lex_<less>> &` — CREF
- [1] `vars` : `const std::vector<variable> &` — CREF
- [2] `alphas` : `const std::vector<Zp> &` — CREF
- [3] `p` : `uint64_t` — VAL

### `__binomial` (2 params)
- [0] `n` : `int64_t` — VAL
- [1] `k` : `int64_t` — VAL

### `__build_cld_matrix` (5 params)
- [0] `M` : `LLLMatrix &` — **REF (output?)**
- [1] `cld` : `const std::vector<upolynomial_<ZZ>> &` — CREF
- [2] `J_cur` : `int` — VAL
- [3] `J_target` : `int` — VAL
- [4] `<unnamed>` : `const ZZ &` — CREF

### `__cld_polys` (3 params)
- [0] `f_star` : `const upolynomial_<ZZ> &` — CREF
- [1] `active_factors` : `const std::vector<upolynomial_<ZZ>> &` — CREF
- [2] `m` : `const ZZ &` — CREF

### `__ddf_Zp` (1 params)
- [0] `f` : `const upolynomial_<Zp> &` — CREF

### `__edf_Zp` (4 params)
- [0] `result` : `std::vector<upolynomial_<Zp>> &` — **REF (output?)**
- [1] `f` : `const upolynomial_<Zp> &` — CREF
- [2] `d` : `uint64_t` — VAL
- [3] `rng` : `std::mt19937 &` — **REF (output?)**

### `__extract_candidates` (3 params)
- [0] `short_rows` : `const std::vector<int> &` — CREF
- [1] `U` : `const LLLMatrix &` — CREF
- [2] `r` : `int` — VAL

### `__extract_monomial_content` (2 params)
- [0] `f` : `const polynomial_<ZZ, lex_<less>> &` — CREF
- [1] `var_factors` : `std::vector<std::pair<variable, int64_t>> &` — **REF (output?)**

### `__extract_pth_root` (1 params)
- [0] `f` : `const upolynomial_<Zp> &` — CREF

### `__factor_Zp` (1 params)
- [0] `f` : `upolynomial_<Zp>` — VAL

### `__factor_multivar` (1 params)
- [0] `f_input` : `const polynomial_<ZZ, lex_<less>> &` — CREF

### `__factor_recombine` (3 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF
- [1] `lifted` : `const std::vector<upolynomial_<ZZ>> &` — CREF
- [2] `m` : `const ZZ &` — CREF

### `__factor_squarefree_primitive_ZZ` (1 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF

### `__hensel_extract_factors` (3 params)
- [0] `nodes` : `const std::vector<__hensel_node> &` — CREF
- [1] `idx` : `int` — VAL
- [2] `factors` : `std::vector<upolynomial_<ZZ>> &` — **REF (output?)**

### `__hensel_lift` (4 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF
- [1] `factors` : `const std::vector<upolynomial_<Zp>> &` — CREF
- [2] `p` : `uint64_t` — VAL
- [3] `a_target` : `int` — VAL

### `__hensel_lift_linear_recursive` (5 params)
- [0] `nodes` : `std::vector<__hensel_node> &` — **REF (output?)**
- [1] `idx` : `int` — VAL
- [2] `f` : `const upolynomial_<ZZ> &` — CREF
- [3] `m` : `const ZZ &` — CREF
- [4] `p` : `uint64_t` — VAL

### `__hensel_lift_recursive` (4 params)
- [0] `nodes` : `std::vector<__hensel_node> &` — **REF (output?)**
- [1] `idx` : `int` — VAL
- [2] `target` : `const upolynomial_<ZZ> &` — CREF
- [3] `m` : `const ZZ &` — CREF

### `__hensel_step` (3 params)
- [0] `node` : `__hensel_node &` — **REF (output?)**
- [1] `f` : `const upolynomial_<ZZ> &` — CREF
- [2] `m` : `const ZZ &` — CREF

### `__hensel_step_linear` (4 params)
- [0] `node` : `__hensel_node &` — **REF (output?)**
- [1] `f` : `const upolynomial_<ZZ> &` — CREF
- [2] `m` : `const ZZ &` — CREF
- [3] `p` : `uint64_t` — VAL

### `__hensel_tree_build` (2 params)
- [0] `factors` : `const std::vector<upolynomial_<Zp>> &` — CREF
- [1] `p` : `uint64_t` — VAL

### `__hensel_tree_build_recursive` (6 params)
- [0] `nodes` : `std::vector<__hensel_node> &` — **REF (output?)**
- [1] `factors` : `const std::vector<upolynomial_<Zp>> &` — CREF
- [2] `p` : `uint64_t` — VAL
- [3] `start` : `int` — VAL
- [4] `end` : `int` — VAL
- [5] `parent_idx` : `int` — VAL

### `__heuristic_starting_precision` (3 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF
- [1] `r` : `int` — VAL
- [2] `p` : `uint64_t` — VAL

### `__isqrt_ceil` (1 params)
- [0] `n` : `const ZZ &` — CREF

### `__lll_factorize` (3 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF
- [1] `factors` : `const std::vector<upolynomial_<Zp>> &` — CREF
- [2] `p` : `uint64_t` — VAL

### `__lll_reduce` (3 params)
- [0] `M` : `LLLMatrix &` — **REF (output?)**
- [1] `U` : `LLLMatrix &` — **REF (output?)**
- [2] `B` : `const ZZ &` — CREF

### `__make_zp` (2 params)
- [0] `val` : `int64_t` — VAL
- [1] `p` : `uint64_t` — VAL

### `__mignotte_bound` (1 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF

### `__mtshl_coeff_bound` (1 params)
- [0] `f` : `const polynomial_<ZZ, lex_<less>> &` — CREF

### `__mtshl_lift` (6 params)
- [0] `f_scaled` : `const polynomial_<ZZ, lex_<less>> &` — CREF
- [1] `scaled_factors` : `const std::vector<upolynomial_<ZZ>> &` — CREF
- [2] `lc_targets` : `const std::vector<polynomial_<ZZ, lex_<less>>> &` — CREF
- [3] `eval_point` : `const std::map<variable, ZZ> &` — CREF
- [4] `main_var` : `const variable &` — CREF
- [5] `p` : `uint64_t` — VAL

### `__mtshl_multi_bdp` (6 params)
- [0] `F` : `const std::vector<polynomial_<Zp, lex_<less>>> &` — CREF
- [1] `c` : `const polynomial_<Zp, lex_<less>> &` — CREF
- [2] `x1` : `const variable &` — CREF
- [3] `x2` : `const variable &` — CREF
- [4] `alpha2` : `const Zp &` — CREF
- [5] `result` : `std::vector<polynomial_<Zp, lex_<less>>> &` — **REF (output?)**

### `__mtshl_sparse_int` (7 params)
- [0] `F` : `const std::vector<polynomial_<Zp, lex_<less>>> &` — CREF
- [1] `c` : `const polynomial_<Zp, lex_<less>> &` — CREF
- [2] `forms` : `const std::vector<std::vector<basic_monomial<lex_<less>>>> &` — CREF
- [3] `x1` : `const variable &` — CREF
- [4] `aux_vars` : `const std::vector<variable> &` — CREF
- [5] `p` : `uint64_t` — VAL
- [6] `result` : `std::vector<polynomial_<Zp, lex_<less>>> &` — **REF (output?)**

### `__mtshl_step_j` (9 params)
- [0] `aj` : `const polynomial_<Zp, lex_<less>> &` — CREF
- [1] `F` : `std::vector<polynomial_<Zp, lex_<less>>> &` — **REF (output?)**
- [2] `lc_tau` : `const std::vector<polynomial_<Zp, lex_<less>>> &` — CREF
- [3] `xj` : `const variable &` — CREF
- [4] `alpha_j` : `const Zp &` — CREF
- [5] `x1` : `const variable &` — CREF
- [6] `aux_vars` : `const std::vector<variable> &` — CREF
- [7] `ideal_alphas_zp` : `const std::vector<Zp> &` — CREF
- [8] `p` : `uint64_t` — VAL

### `__mtshl_wmds` (6 params)
- [0] `F` : `const std::vector<polynomial_<Zp, lex_<less>>> &` — CREF
- [1] `c` : `const polynomial_<Zp, lex_<less>> &` — CREF
- [2] `x1` : `const variable &` — CREF
- [3] `aux_vars` : `const std::vector<variable> &` — CREF
- [4] `ideal_alphas_zp` : `const std::vector<Zp> &` — CREF
- [5] `result` : `std::vector<polynomial_<Zp, lex_<less>>> &` — **REF (output?)**

### `__mtshl_zp_univar_mdp` (3 params)
- [0] `F` : `const std::vector<upolynomial_<Zp>> &` — CREF
- [1] `c` : `const upolynomial_<Zp> &` — CREF
- [2] `sigma` : `std::vector<upolynomial_<Zp>> &` — **REF (output?)**

### `__polynomial_to_zp` (2 params)
- [0] `f` : `const polynomial_<ZZ, lex_<less>> &` — CREF
- [1] `p` : `uint64_t` — VAL

### `__select_eval_point` (4 params)
- [0] `f` : `const polynomial_<ZZ, lex_<less>> &` — CREF
- [1] `main_var` : `const variable &` — CREF
- [2] `skip` : `int` — VAL
- [3] `lc_coprime_mod` : `uint64_t` — VAL

### `__select_prime` (2 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF
- [1] `use_large_prime` : `bool` — VAL

### `__si_theta_array_eval` (6 params)
- [0] `f` : `const polynomial_<Zp, lex_<less>> &` — CREF
- [1] `x1` : `const variable &` — CREF
- [2] `aux_vars` : `const std::vector<variable> &` — CREF
- [3] `sparse_betas` : `const std::vector<Zp> &` — CREF
- [4] `s` : `int` — VAL
- [5] `images` : `std::vector<upolynomial_<Zp>> &` — **REF (output?)**

### `__si_vandermonde_solve` (3 params)
- [0] `values` : `const std::vector<Zp> &` — CREF
- [1] `thetas` : `const std::vector<Zp> &` — CREF
- [2] `coeffs` : `std::vector<Zp> &` — **REF (output?)**

### `__squarefree_Zp` (1 params)
- [0] `f` : `const upolynomial_<Zp> &` — CREF

### `__subset_product_mod` (4 params)
- [0] `factors` : `const std::vector<upolynomial_<ZZ>> &` — CREF
- [1] `subset` : `const std::vector<size_t> &` — CREF
- [2] `lc_f` : `const ZZ &` — CREF
- [3] `m` : `const ZZ &` — CREF

### `__symmetric_mod` (2 params)
- [0] `a` : `const ZZ &` — CREF
- [1] `m` : `const ZZ &` — CREF

### `__symmetric_mod_poly` (2 params)
- [0] `f` : `const polynomial_<Zp, lex_<less>> &` — CREF
- [1] `p` : `uint64_t` — VAL

### `__taylor_coeff_zp` (4 params)
- [0] `f` : `const polynomial_<Zp, lex_<less>> &` — CREF
- [1] `xk` : `const variable &` — CREF
- [2] `alpha_k` : `const Zp &` — CREF
- [3] `j` : `int` — VAL

### `__upoly_const_term` (1 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF

### `__upoly_divmod` (4 params)
- [0] `q` : `upolynomial_<Zp> &` — **REF (output?)**
- [1] `r` : `upolynomial_<Zp> &` — **REF (output?)**
- [2] `f` : `const upolynomial_<Zp> &` — CREF
- [3] `g` : `const upolynomial_<Zp> &` — CREF

### `__upoly_divmod_mod` (5 params)
- [0] `q` : `upolynomial_<ZZ> &` — **REF (output?)**
- [1] `r` : `upolynomial_<ZZ> &` — **REF (output?)**
- [2] `f` : `const upolynomial_<ZZ> &` — CREF
- [3] `g` : `const upolynomial_<ZZ> &` — CREF
- [4] `m` : `const ZZ &` — CREF

### `__upoly_make_monic` (1 params)
- [0] `f` : `upolynomial_<Zp> &` — **REF (output?)**

### `__upoly_mod` (2 params)
- [0] `f` : `const upolynomial_<Zp> &` — CREF
- [1] `g` : `const upolynomial_<Zp> &` — CREF

### `__upoly_mod_coeff` (2 params)
- [0] `f` : `upolynomial_<ZZ> &` — **REF (output?)**
- [1] `m` : `const ZZ &` — CREF

### `__upoly_mul_mod` (3 params)
- [0] `a` : `const upolynomial_<ZZ> &` — CREF
- [1] `b` : `const upolynomial_<ZZ> &` — CREF
- [2] `m` : `const ZZ &` — CREF

### `__upoly_norm_l1` (1 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF

### `__upoly_norm_l2_sq` (1 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF

### `__upoly_powmod` (3 params)
- [0] `base` : `const upolynomial_<Zp> &` — CREF
- [1] `exp` : `const ZZ &` — CREF
- [2] `modpoly` : `const upolynomial_<Zp> &` — CREF

### `__upoly_primitive` (1 params)
- [0] `f` : `upolynomial_<ZZ>` — VAL

### `__upoly_random` (3 params)
- [0] `max_deg` : `int64_t` — VAL
- [1] `p` : `uint64_t` — VAL
- [2] `rng` : `std::mt19937 &` — **REF (output?)**

### `__upoly_subtract_one` (2 params)
- [0] `h` : `const upolynomial_<Zp> &` — CREF
- [1] `p` : `uint64_t` — VAL

### `__upoly_subtract_x` (2 params)
- [0] `h` : `const upolynomial_<Zp> &` — CREF
- [1] `p` : `uint64_t` — VAL

### `__upoly_symmetric_mod` (2 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF
- [1] `m` : `const ZZ &` — CREF

### `__upoly_to_poly` (3 params)
- [0] `up` : `const upolynomial_<ZZ> &` — CREF
- [1] `var` : `const variable &` — CREF
- [2] `comp_ptr` : `const lex_<less> *` — PTR

### `__vanhoeij_recombine` (3 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF
- [1] `lifted` : `const std::vector<upolynomial_<ZZ>> &` — CREF
- [2] `m` : `const ZZ &` — CREF

### `__wang_core` (1 params)
- [0] `g` : `const polynomial_<ZZ, lex_<less>> &` — CREF

### `__wang_leading_coeff` (5 params)
- [0] `f` : `const polynomial_<ZZ, lex_<less>> &` — CREF
- [1] `univar_factors` : `const std::vector<upolynomial_<ZZ>> &` — CREF
- [2] `eval_point` : `const std::map<variable, ZZ> &` — CREF
- [3] `main_var` : `const variable &` — CREF
- [4] `uni_content` : `const ZZ &` — CREF

### `__zassenhaus_recombine` (3 params)
- [0] `f` : `const upolynomial_<ZZ> &` — CREF
- [1] `lifted` : `const std::vector<upolynomial_<ZZ>> &` — CREF
- [2] `m` : `const ZZ &` — CREF

### `factorize` (1 params)
- [0] `F` : `const polynomial_<ZZ, grlex_<less>> &` — CREF
