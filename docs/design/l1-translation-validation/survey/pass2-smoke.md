# Pass 1 + Pass 2 全量烟测（65 函数）

- 总函数数：**65**
- OK: **65**
  - 有 ref 参数被消除: **22**
  - 无 ref 参数（passthrough）: **43**
- FAIL: **0**

## 转换后函数详情（按 ref 参数数分组）

### 0 ref 参数（43 函数）

- `__assign_partial_zp`
- `__binomial`
- `__cld_polys`
- `__ddf_Zp`
- `__extract_candidates`
- `__extract_pth_root`
- `__factor_Zp`
- `__factor_multivar`
- `__factor_recombine`
- `__factor_squarefree_primitive_ZZ`
- `__hensel_lift`
- `__hensel_tree_build`
- `__heuristic_starting_precision`
- `__isqrt_ceil`
- `__lll_factorize`
- `__make_zp`
- `__mignotte_bound`
- `__mtshl_coeff_bound`
- `__mtshl_lift`
- `__polynomial_to_zp`
- `__select_eval_point`
- `__select_prime`
- `__squarefree_Zp`
- `__subset_product_mod`
- `__symmetric_mod`
- `__symmetric_mod_poly`
- `__taylor_coeff_zp`
- `__upoly_const_term`
- `__upoly_mod`
- `__upoly_mul_mod`
- `__upoly_norm_l1`
- `__upoly_norm_l2_sq`
- `__upoly_powmod`
- `__upoly_primitive`
- `__upoly_subtract_one`
- `__upoly_subtract_x`
- `__upoly_symmetric_mod`
- `__upoly_to_poly`
- `__vanhoeij_recombine`
- `__wang_core`
- `__wang_leading_coeff`
- `__zassenhaus_recombine`
- `factorize`

### 1 ref 参数（18 函数）

| 函数 | ref 参数名 | Pass 1 返回类型 → Pass 2 返回类型 |
|---|---|---|
| `__build_cld_matrix` | M | `BaseType.INT32` → `PairType(fst=<BaseType.INT32: 'Int32'>, snd=NamedType(name='LLLMatrix'))` |
| `__extract_monomial_content` | var_factors | `NamedType(name='MvPolyZZ')` → `PairType(fst=NamedType(name='MvPolyZZ'), snd=ArrayType(elem=PairType(fst=NamedType(name='Variable'), snd=<BaseType.INT64: 'Int64'>)))` |
| `__hensel_extract_factors` | factors | `BaseType.UNIT` → `ArrayType(elem=NamedType(name='SparsePolyZZ'))` |
| `__hensel_lift_linear_recursive` | nodes | `BaseType.UNIT` → `ArrayType(elem=NamedType(name='HenselNode'))` |
| `__hensel_lift_recursive` | nodes | `BaseType.UNIT` → `ArrayType(elem=NamedType(name='HenselNode'))` |
| `__hensel_step` | node | `BaseType.UNIT` → `NamedType(name='HenselNode')` |
| `__hensel_step_linear` | node | `BaseType.UNIT` → `NamedType(name='HenselNode')` |
| `__hensel_tree_build_recursive` | nodes | `BaseType.UNIT` → `ArrayType(elem=NamedType(name='HenselNode'))` |
| `__mtshl_multi_bdp` | result | `BaseType.BOOL` → `PairType(fst=<BaseType.BOOL: 'Bool'>, snd=ArrayType(elem=NamedType(name='MvPolyZp')))` |
| `__mtshl_sparse_int` | result | `BaseType.BOOL` → `PairType(fst=<BaseType.BOOL: 'Bool'>, snd=ArrayType(elem=NamedType(name='MvPolyZp')))` |
| `__mtshl_step_j` | F | `BaseType.BOOL` → `PairType(fst=<BaseType.BOOL: 'Bool'>, snd=ArrayType(elem=NamedType(name='MvPolyZp')))` |
| `__mtshl_wmds` | result | `BaseType.BOOL` → `PairType(fst=<BaseType.BOOL: 'Bool'>, snd=ArrayType(elem=NamedType(name='MvPolyZp')))` |
| `__mtshl_zp_univar_mdp` | sigma | `BaseType.BOOL` → `PairType(fst=<BaseType.BOOL: 'Bool'>, snd=ArrayType(elem=NamedType(name='SparsePolyZp')))` |
| `__si_theta_array_eval` | images | `BaseType.UNIT` → `ArrayType(elem=NamedType(name='SparsePolyZp'))` |
| `__si_vandermonde_solve` | coeffs | `BaseType.BOOL` → `PairType(fst=<BaseType.BOOL: 'Bool'>, snd=ArrayType(elem=NamedType(name='Zp')))` |
| `__upoly_make_monic` | f | `NamedType(name='Zp')` → `PairType(fst=NamedType(name='Zp'), snd=NamedType(name='SparsePolyZp'))` |
| `__upoly_mod_coeff` | f | `BaseType.UNIT` → `NamedType(name='SparsePolyZZ')` |
| `__upoly_random` | rng | `NamedType(name='SparsePolyZp')` → `PairType(fst=NamedType(name='SparsePolyZp'), snd=NamedType(name='Rng'))` |

### 2 ref 参数（4 函数）

| 函数 | ref 参数名 | Pass 1 返回类型 → Pass 2 返回类型 |
|---|---|---|
| `__edf_Zp` | result, rng | `BaseType.UNIT` → `PairType(fst=ArrayType(elem=NamedType(name='SparsePolyZp')), snd=NamedType(name='Rng'))` |
| `__lll_reduce` | M, U | `ArrayType(elem=<BaseType.INT32: 'Int32'>)` → `TupleType(elems=(ArrayType(elem=<BaseType.INT32: 'Int32'>), NamedType(name='LLLMatrix'), NamedType(name='LLLMatrix')))` |
| `__upoly_divmod` | q, r | `BaseType.UNIT` → `PairType(fst=NamedType(name='SparsePolyZp'), snd=NamedType(name='SparsePolyZp'))` |
| `__upoly_divmod_mod` | q, r | `BaseType.UNIT` → `PairType(fst=NamedType(name='SparsePolyZZ'), snd=NamedType(name='SparsePolyZZ'))` |

## 与 `ref_params.md` 预期值核对

- 本轮测得 ref 参数总数: **26**
- `ref_params.md` 预期（Week 2 Day 1）: 26

✅ 与预期值完全一致
