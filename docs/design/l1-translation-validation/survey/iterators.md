# 迭代器模式扫描

为 HIR `iter_recognize` Pass 设计输入。

## 全局统计

| 模式 | 数量 | 备注 |
|---|---|---|
| Range-for `for (auto& x : v)` | **92** | AST `CXXForRangeStmt` 统计 |
| Range-for 含结构化绑定 `for (auto& [k,v] : m)` | **24** | 需要 HIR 特殊处理 |
| Classic iterator loop `for (auto it = X.begin();...;++it)` | **1** | 源码正则 |
| Iterator 初始化 `auto X = Y.begin()` | **7** | 源码正则（潜在双指针 compact 信号）|
| One-liner 双指针 (`auto it = X.begin(), out = ...`) | **1** | 源码正则 |
| Compact erase (`X.erase(out, X.end())`) | **4** | 源码正则（严格匹配）|
| 任意 `.erase(...)` 调用 | **10** | 源码正则（含位置删除、range 删除）|
| Structured-binding range-for | **38** | 源码正则 |

## STL 算法调用直方图

| 算法 | 总次数 | 首次出现函数 |
|---|---|---|
| `std::move` | 63 | `factorize`(8), `__factor_multivar`(6), `__upoly_divmod_mod`(6), ... (共 22) |
| `std::sort` | 8 | `__factor_multivar`(2), `__zassenhaus_recombine`(1), `__lll_reduce`(1), ... (共 7) |
| `std::max` | 7 | `__heuristic_starting_precision`(1), `__build_cld_matrix`(1), `__lll_reduce`(1), ... (共 7) |
| `std::iota` | 5 | `__zassenhaus_recombine`(2), `__wang_core`(2), `__vanhoeij_recombine`(1) |
| `std::swap` | 4 | `__lll_reduce`(3), `__si_vandermonde_solve`(1) |
| `std::min` | 3 | `__heuristic_starting_precision`(1), `__vanhoeij_recombine`(1), `__extract_monomial_content`(1) |

## Range-for 分布（按宿主函数）

| 宿主 | range-for | 其中含结构化绑定 |
|---|---|---|
| `__factor_multivar` | 9 | 7 |
| `__mtshl_lift` | 7 | 1 |
| `__extract_monomial_content` | 6 | 6 |
| `__hensel_step` | 6 | 0 |
| `__mtshl_sparse_int` | 6 | 1 |
| `__vanhoeij_recombine` | 6 | 1 |
| `__wang_core` | 6 | 2 |
| `__build_cld_matrix` | 4 | 0 |
| `__select_eval_point` | 4 | 2 |
| `__factor_Zp` | 3 | 0 |
| `__si_theta_array_eval` | 3 | 1 |
| `__hensel_lift` | 2 | 0 |
| `__hensel_step_linear` | 2 | 0 |
| `__mtshl_multi_bdp` | 2 | 0 |
| `__mtshl_step_j` | 2 | 0 |
| `__mtshl_wmds` | 2 | 0 |
| `__select_prime` | 2 | 1 |
| `__squarefree_Zp` | 2 | 0 |
| `__wang_leading_coeff` | 2 | 1 |
| `__zassenhaus_recombine` | 2 | 0 |
| `__cld_polys` | 1 | 0 |
| `__extract_pth_root` | 1 | 0 |
| `__mtshl_coeff_bound` | 1 | 0 |
| `__polynomial_to_zp` | 1 | 0 |
| `__subset_product_mod` | 1 | 0 |
| `__symmetric_mod_poly` | 1 | 0 |
| `__upoly_make_monic` | 1 | 0 |
| `__upoly_norm_l1` | 1 | 0 |
| `__upoly_norm_l2_sq` | 1 | 0 |
| `__upoly_primitive` | 1 | 0 |
| `__upoly_subtract_one` | 1 | 0 |
| `__upoly_subtract_x` | 1 | 0 |
| `__upoly_symmetric_mod` | 1 | 0 |
| `factorize` | 1 | 1 |

## Classic iterator loop 详情

### `__extract_monomial_content` (1)
- polynomial_factorize_wang.hh:2193 — `for (auto it = min_deg.begin()`

## 双指针 compact 详情

### `__hensel_step_linear` (1)
- polynomial_factorize_univar.hh:649 — `auto it = e.data().begin(), out =`

## Compact erase 详情（X.erase(out, X.end()) 严格匹配）

### `__hensel_step` (2)
- polynomial_factorize_univar.hh:429 — `e.data().erase(out, e.data().end())`
- polynomial_factorize_univar.hh:475 — `ep.data().erase(out, ep.data().end())`

### `__hensel_step_linear` (1)
- polynomial_factorize_univar.hh:656 — `e.data().erase(out, e.data().end())`

### `__upoly_mod_coeff` (1)
- polynomial_factorize_univar.hh:205 — `f.data().erase(out, f.data().end())`

## Iterator 初始化详情（auto X = Y.begin()，潜在 compact 信号）

### `__upoly_divmod_mod` (2)
- polynomial_factorize_univar.hh:245 — `auto r_it = r.data().begin()`
- polynomial_factorize_univar.hh:246 — `auto g_it = g.begin()`

### `__extract_monomial_content` (2)
- polynomial_factorize_wang.hh:2181 — `auto it = min_deg.begin()`
- polynomial_factorize_wang.hh:2193 — `auto it = min_deg.begin()`

### `__upoly_mod_coeff` (1)
- polynomial_factorize_univar.hh:193 — `auto it = f.data().begin()`

### `__hensel_step` (1)
- polynomial_factorize_univar.hh:423 — `auto it = e.data().begin()`

### `__hensel_step_linear` (1)
- polynomial_factorize_univar.hh:649 — `auto it = e.data().begin()`

## `.erase(...)` 全部调用

### `__hensel_step` (2)
- polynomial_factorize_univar.hh:429 — `e.data().erase(`
- polynomial_factorize_univar.hh:475 — `ep.data().erase(`

### `__extract_monomial_content` (2)
- polynomial_factorize_wang.hh:2185 — `min_deg.erase(`
- polynomial_factorize_wang.hh:2196 — `min_deg.erase(`

### `__upoly_mod_coeff` (1)
- polynomial_factorize_univar.hh:205 — `f.data().erase(`

### `__upoly_divmod_mod` (1)
- polynomial_factorize_univar.hh:236 — `r.data().erase(`

### `__hensel_step_linear` (1)
- polynomial_factorize_univar.hh:656 — `e.data().erase(`

### `__zassenhaus_recombine` (1)
- polynomial_factorize_univar.hh:856 — `T.erase(`

### `__vanhoeij_recombine` (1)
- polynomial_factorize_univar.hh:1304 — `active.erase(`

### `__wang_core` (1)
- polynomial_factorize_wang.hh:2396 — `mv_T.erase(`

## Structured-binding range-for 详情

### `__extract_monomial_content` (6)
- polynomial_factorize_wang.hh:2162 — `for (const auto& [mono, coeff] :`
- polynomial_factorize_wang.hh:2166 — `for (const auto& [var, deg] :`
- polynomial_factorize_wang.hh:2173 — `for (const auto& [var, deg] :`
- polynomial_factorize_wang.hh:2205 — `for (const auto& [mono, coeff] :`
- polynomial_factorize_wang.hh:2208 — `for (const auto& [var, deg] :`
- polynomial_factorize_wang.hh:2220 — `for (const auto& [var, deg] :`

### `__factor_multivar` (9)
- polynomial_factorize.hh:71 — `for (auto& [sqf_factor, mult] :`
- polynomial_factorize.hh:118 — `for (const auto& [fi, ei] :`
- polynomial_factorize_wang.hh:2453 — `for (auto& [gk, mk] :`
- polynomial_factorize_wang.hh:2469 — `for (auto& [fi, ei] :`
- polynomial_factorize_wang.hh:2479 — `for (auto& [var, vdeg] :`
- polynomial_factorize_wang.hh:2497 — `for (auto& [fi, ei] :`
- polynomial_factorize_wang.hh:2512 — `for (auto& [fi, ei] :`
- polynomial_factorize_wang.hh:2525 — `for (auto& [fac, mult] :`
- polynomial_factorize_wang.hh:2547 — `for (const auto& [fi, ei] :`

### `__mtshl_lift` (1)
- polynomial_factorize_wang.hh:1032 — `for (const auto& [v, a] :`

### `__mtshl_sparse_int` (1)
- polynomial_factorize_wang.hh:555 — `for (auto& [d, indices] :`

### `__select_eval_point` (2)
- polynomial_factorize_wang.hh:1204 — `for (auto& [v, d] :`
- polynomial_factorize_wang.hh:1218 — `for (auto& [lj, ej] :`

### `__si_theta_array_eval` (1)
- polynomial_factorize_wang.hh:198 — `for (auto& [deg, coeff] :`

### `__upoly_make_monic` (1)
- polynomial_factorize_univar.hh:1438 — `for (auto& [gk, dk] :`

### `__wang_core` (2)
- polynomial_factorize_wang.hh:2246 — `for (auto& [v, d] :`
- polynomial_factorize_wang.hh:2304 — `for (auto& [fi, ei] :`

### `__wang_leading_coeff` (10)
- polynomial_factorize_wang.hh:1383 — `for (auto& [lj2, ej2] :`
- polynomial_factorize_wang.hh:1560 — `for (auto& [var, deg] :`
- polynomial_factorize_wang.hh:1572 — `for (auto& [d, p] :`
- polynomial_factorize_wang.hh:1587 — `for (auto& [deg, p] :`
- polynomial_factorize_wang.hh:1603 — `for (auto& [mono, coeff] :`
- polynomial_factorize_wang.hh:1616 — `for (auto& [d, p] :`
- polynomial_factorize_wang.hh:1618 — `for (auto& [mono, coeff] :`
- polynomial_factorize_wang.hh:1624 — `for (auto& [var, deg] :`
- polynomial_factorize_wang.hh:1912 — `for (auto& [v, val] :`
- polynomial_factorize_wang.hh:2092 — `for (auto& [v, d] :`

### `factorize` (5)
- polynomial_factorize.hh:148 — `for (auto& [fac, mult] :`
- polynomial_factorize.hh:199 — `for (auto& [fac_zz, mult] :`
- polynomial_factorize.hh:222 — `for (const auto& [fi, ei] :`
- polynomial_factorize_univar.hh:1598 — `for (auto& [sqf_factor, mult] :`
- polynomial_factorize_univar.hh:1643 — `for (const auto& [fi, ei] :`
