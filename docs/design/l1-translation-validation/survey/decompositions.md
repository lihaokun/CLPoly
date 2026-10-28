# 结构化绑定（DecompositionDecl）扫描

总计 **30** 处，分布于 **13** 个函数。

- 在 range-for 中: **25**
- 独立 DeclStmt: **5**

## 绑定 container 类型

- pair: 30

## 常见绑定名模式（前 20）

- `[fi, ei]`: 5
- `[var, deg]`: 4
- `[mono, coeff]`: 2
- `[fac, mult]`: 2
- `[v, d]`: 2
- `[c_g, pp_g]`: 2
- `[c_q, pp_q]`: 2
- `[gk, mk]`: 1
- `[var, vdeg]`: 1
- `[a_h, a_mig]`: 1
- `[lifted_h, m_h]`: 1
- `[lifted_mig, m_mig]`: 1
- `[v, a]`: 1
- `[d, indices]`: 1
- `[lj, ej]`: 1
- `[gk, dk]`: 1
- `[deg, coeff]`: 1
- `[lj2, ej2]`: 1

## 按函数详情

### `__extract_monomial_content` (6)
- `[mono, coeff]` [range-for] — container: `const std::pair<clpoly::basic_monomial<clpoly::lex_<clpoly::less>>, clpoly::ZZ>`
- `[var, deg]` [range-for] — container: `const std::pair<clpoly::variable, long>`
- `[var, deg]` [range-for] — container: `const std::pair<clpoly::variable, long>`
- `[mono, coeff]` [range-for] — container: `const std::pair<clpoly::basic_monomial<clpoly::lex_<clpoly::less>>, clpoly::ZZ>`
- `[var, deg]` [range-for] — container: `const std::pair<clpoly::variable, long>`
- `[var, deg]` [range-for] — container: `std::pair<const clpoly::variable, long> const`

### `__factor_multivar` (7)
- `[gk, mk]` [range-for] — container: `std::pair<clpoly::basic_polynomial<clpoly::basic_monomial<clpoly::lex_<clpoly::l`
- `[fi, ei]` [range-for] — container: `std::pair<clpoly::basic_polynomial<clpoly::basic_monomial<clpoly::lex_<clpoly::l`
- `[var, vdeg]` [range-for] — container: `std::pair<clpoly::variable, long>`
- `[fi, ei]` [range-for] — container: `std::pair<clpoly::basic_polynomial<clpoly::basic_monomial<clpoly::lex_<clpoly::l`
- `[fi, ei]` [range-for] — container: `std::pair<clpoly::basic_polynomial<clpoly::basic_monomial<clpoly::lex_<clpoly::l`
- `[fac, mult]` [range-for] — container: `std::pair<clpoly::basic_polynomial<clpoly::basic_monomial<clpoly::lex_<clpoly::l`
- `[fi, ei]` [range-for] — container: `std::pair<clpoly::basic_polynomial<clpoly::basic_monomial<clpoly::lex_<clpoly::l`

### `__lll_factorize` (3)
- `[a_h, a_mig]` — container: `std::pair<int, int>`
- `[lifted_h, m_h]` — container: `std::pair<std::vector<upolynomial_<ZZ>>, ZZ>`
- `[lifted_mig, m_mig]` — container: `std::pair<std::vector<upolynomial_<ZZ>>, ZZ>`

### `__mtshl_lift` (1)
- `[v, a]` [range-for] — container: `const std::pair<const clpoly::variable, clpoly::ZZ>`

### `__mtshl_sparse_int` (1)
- `[d, indices]` [range-for] — container: `std::pair<const long, std::vector<int>>`

### `__select_eval_point` (2)
- `[v, d]` [range-for] — container: `std::pair<clpoly::variable, long>`
- `[lj, ej]` [range-for] — container: `std::pair<clpoly::basic_polynomial<clpoly::basic_monomial<clpoly::lex_<clpoly::l`

### `__select_prime` (1)
- `[gk, dk]` [range-for] — container: `std::pair<clpoly::basic_polynomial<clpoly::umonomial, clpoly::Zp, clpoly::uless>`

### `__si_theta_array_eval` (1)
- `[deg, coeff]` [range-for] — container: `std::pair<const long, clpoly::Zp>`

### `__vanhoeij_recombine` (2)
- `[c_g, pp_g]` [range-for] — container: `std::pair<ZZ, upolynomial_<ZZ>>`
- `[c_q, pp_q]` [range-for] — container: `std::pair<ZZ, upolynomial_<ZZ>>`

### `__wang_core` (2)
- `[v, d]` [range-for] — container: `std::pair<clpoly::variable, long>`
- `[fi, ei]` [range-for] — container: `std::pair<clpoly::basic_polynomial<clpoly::umonomial, clpoly::ZZ, clpoly::uless>`

### `__wang_leading_coeff` (1)
- `[lj2, ej2]` [range-for] — container: `std::pair<clpoly::basic_polynomial<clpoly::basic_monomial<clpoly::lex_<clpoly::l`

### `__zassenhaus_recombine` (2)
- `[c_g, pp_g]` — container: `std::pair<ZZ, upolynomial_<ZZ>>`
- `[c_q, pp_q]` — container: `std::pair<ZZ, upolynomial_<ZZ>>`

### `factorize` (1)
- `[fac, mult]` [range-for] — container: `std::pair<clpoly::basic_polynomial<clpoly::basic_monomial<clpoly::lex_<clpoly::l`
