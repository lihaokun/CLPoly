# Lambda 扫描

总计 **28** 个 lambda，其中 **26** 个在 TRANSLATION_SCOPE 内。

> 源码正则扫描（Clang JSON 对 LambdaExpr 的 captures 字段省略，无法用 AST）。

## 全局统计（仅 in-scope）

- Generic（`auto` 参数）: **0** — 若仍有请再修 Clang 单态化问题
- 捕获模式 `[]`（无捕获）: **13**
- 捕获模式 `[&]`（默认引用）: **10**
- 捕获模式 `[=]`（默认值）: **0**
- 捕获模式 `[x, y, ...]`（具名）: **3**
- 捕获模式 `[&, x, ...]` / `[=, &x, ...]`（混合）: **0**

## Body 行数直方图（in-scope）

- 1 行: 2
- 2-3 行: 6
- 4-10 行: 12
- 11-30 行: 6

## 按宿主函数详情

### `__build_cld_matrix` (1 lambda)

**Lambda 1** — polynomial_factorize_univar.hh:947 — `[]` — none
- params: `const upolynomial_<ZZ>& p, int deg`
- body: 5 行
- body preview: `{ |             for (const auto& term : p) |                 if ((int)term.first.deg() == deg) return term.second; |             return ZZ(0); |      `

### `__factor_Zp` (1 lambda)

**Lambda 1** — polynomial_factorize_zp.hh:386 — `[]` — none
- params: `const std::pair<upolynomial_<Zp>, uint64_t>& a,
               const std::pair<upolynomial_<Zp>, uin`
- body: 3 行
- body preview: `{ |                 return get_deg(a.first) < get_deg(b.first); |             }`

### `__factor_multivar` (2 lambdas)

**Lambda 1** — polynomial_factorize.hh:108 — `[]` — none
- params: `const std::pair<Poly, uint64_t>& a,
               const std::pair<Poly, uint64_t>& b`
- body: 3 行
- body preview: `{ |                 return degree(a.first) < degree(b.first); |             }`

**Lambda 2** — polynomial_factorize_wang.hh:2537 — `[]` — none
- params: `const std::pair<Poly, uint64_t>& a,
               const std::pair<Poly, uint64_t>& b`
- body: 3 行
- body preview: `{ |                 return degree(a.first) < degree(b.first); |             }`

### `__lll_reduce` (5 lambdas)

**Lambda 1** — polynomial_factorize_univar.hh:1000 — `[j]` — explicit
- params: `i > j`
- body: 6 行
- 具名捕获: =j
- body preview: `{ |             ZZ s(0); |             for (int k = 0; k < (int)a.size(); ++k) |                 s += a[k] * b[k]; |             return s; |         }`

**Lambda 2** — polynomial_factorize_univar.hh:1019 — `[]` — none
- params: `const QQ& q`
- body: 7 行
- body preview: `{ |             ZZ a = q.get_num() * ZZ(2) + q.get_den(); |             ZZ b = q.get_den() * ZZ(2); |             ZZ result; |             ZZ::fdiv_q(`

**Lambda 3** — polynomial_factorize_univar.hh:1028 — `[&]` — implicit-ref
- params: `int i, int j, const ZZ& c`
- body: 6 行
- body preview: `{ |             for (int k = 0; k < n; ++k) { |                 M[i][k] -= c * M[j][k]; |                 U[i][k] -= c * U[j][k]; |             } |   `

**Lambda 4** — polynomial_factorize_univar.hh:1034 — `[&]` — implicit-ref
- params: `int i, int j`
- body: 4 行
- body preview: `{ |             std::swap(M[i], M[j]); |             std::swap(U[i], U[j]); |         }`

**Lambda 5** — polynomial_factorize_univar.hh:1127 — `[&]` — implicit-ref
- params: `int a, int b`
- body: 1 行
- body preview: `{ return dot(M[a], M[a]) < dot(M[b], M[b]); }`

### `__mtshl_multi_bdp` (1 lambda)

**Lambda 1** — polynomial_factorize_wang.hh:361 — `[&]` — implicit-ref
- params: ``
- body: 12 行
- body preview: `{ |             PolyZp e = c; |             for (int i = 0; i < r; i++) |             { |                 if (result[i].empty()) continue; |          `

### `__mtshl_sparse_int` (2 lambdas)

**Lambda 1** — polynomial_factorize_wang.hh:519 — `[p]` — explicit
- params: `const UPZp& poly, int64_t d`
- body: 6 行
- 具名捕获: =p
- body preview: `{ |             for (const auto& term : poly) |                 if (term.first.deg() == d) |                     return term.second; |             ret`

**Lambda 2** — polynomial_factorize_wang.hh:527 — `[&]` — implicit-ref
- params: `const basic_monomial<lex_<var_order>>& mono`
- body: 12 行
- body preview: `{ |             Zp tm(1, p); |             for (size_t k = 0; k < aux_vars.size(); k++) |             { |                 int64_t ek = 0; |           `

### `__mtshl_step_j` (2 lambdas)

**Lambda 1** — polynomial_factorize_wang.hh:783 — `[&]` — implicit-ref
- params: `PolyZp& Gi, const PolyZp& lc_target`
- body: 16 行
- body preview: `{ |             if (Gi.empty()) return; |             auto lc_cur = leadcoeff(Gi, x1); |             auto diff = lc_target - lc_cur; |             dif`

**Lambda 2** — polynomial_factorize_wang.hh:828 — `[&]` — implicit-ref
- params: ``
- body: 9 行
- body preview: `{ |             PolyZp prod = F[0]; |             for (int i = 1; i < r; i++) |             { |                 prod = prod * F[i]; |                 `

### `__mtshl_wmds` (1 lambda)

**Lambda 1** — polynomial_factorize_wang.hh:691 — `[&]` — implicit-ref
- params: ``
- body: 12 行
- body preview: `{ |             PolyZp e = c; |             for (int i = 0; i < r; i++) |             { |                 if (result[i].empty()) continue; |          `

### `__select_prime` (1 lambda)

**Lambda 1** — polynomial_factorize_univar.hh:1405 — `[use_large_prime]` — explicit
- params: `uint64_t cur`
- body: 6 行
- 具名捕获: =use_large_prime
- body preview: `{ |             if (use_large_prime) |                 return prev_prime_64(cur); |             else |                 return next_prime_64(cur); |   `

### `__vanhoeij_recombine` (2 lambdas)

**Lambda 1** — polynomial_factorize_univar.hh:1212 — `[]` — none
- params: `int rr, int U_exp_`
- body: 7 行
- body preview: `{ |             LLLMatrix M(rr, std::vector<ZZ>(rr, ZZ(0))); |             ZZ scale = ZZ(1) << U_exp_; |             for (int i = 0; i < rr; ++i) |   `

**Lambda 2** — polynomial_factorize_univar.hh:1337 — `[]` — none
- params: `const upolynomial_<ZZ>& a, const upolynomial_<ZZ>& b`
- body: 3 行
- body preview: `{ |                 return get_deg(a) < get_deg(b); |             }`

### `__wang_core` (3 lambdas)

**Lambda 1** — polynomial_factorize_wang.hh:2258 — `[]` — none
- params: `Poly& h`
- body: 6 行
- body preview: `{ |             h = pp(h); |             if (!h.empty() && h.front().second < 0) |                 for (auto& term : h.data()) |                     t`

**Lambda 2** — polynomial_factorize_wang.hh:2283 — `[&]` — implicit-ref
- params: `uint64_t p`
- body: 5 行
- body preview: `{ |                         for (const auto& term : L) |                             if (term.second.fdiv_ui(p) != 0) return false; |                 `

**Lambda 3** — polynomial_factorize_wang.hh:2350 — `[]` — none
- params: `std::vector<int>& idx, int n`
- body: 11 行
- body preview: `{ |                         int sz = (int)idx.size(); |                         int i = sz - 1; |                         while (i >= 0 && idx[i] == n`

### `__wang_leading_coeff` (2 lambdas)

**Lambda 1** — polynomial_factorize_wang.hh:1335 — `[&]` — implicit-ref
- params: `const ZZ& val`
- body: 6 行
- body preview: `{ |             Poly p(comp_ptr); |             if (val != 0) |                 p.push_back({basic_monomial<lex_<var_order>>(comp_ptr), val}); |      `

**Lambda 2** — polynomial_factorize_wang.hh:1372 — `[]` — none
- params: `const std::pair<Poly, uint64_t>& a,
                   const std::pair<Poly, uint64_t>& b`
- body: 1 行
- body preview: `{ return a.second > b.second; }`

### `__zassenhaus_recombine` (2 lambdas)

**Lambda 1** — polynomial_factorize_univar.hh:774 — `[]` — none
- params: `std::vector<int>& idx, int n`
- body: 11 行
- body preview: `{ |             int s = (int)idx.size(); |             int i = s - 1; |             while (i >= 0 && idx[i] == n - s + i) |                 --i; |    `

**Lambda 2** — polynomial_factorize_univar.hh:877 — `[]` — none
- params: `const upolynomial_<ZZ>& a, const upolynomial_<ZZ>& b`
- body: 3 行
- body preview: `{ |                 return get_deg(a) < get_deg(b); |             }`

### `factorize` (1 lambda)

**Lambda 1** — polynomial_factorize_univar.hh:1628 — `[]` — none
- params: `const std::pair<upolynomial_<ZZ>, uint64_t>& a,
               const std::pair<upolynomial_<ZZ>, uin`
- body: 3 行
- body preview: `{ |                 return degree(a.first) < degree(b.first); |             }`

## Out-of-scope lambda（宿主函数不在 TRANSLATION_SCOPE）

- `None` @ polynomial_factorize_univar.hh:55 — `[&]`
- `None` @ polynomial_factorize_univar.hh:56 — `[]`
