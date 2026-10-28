# STL 依赖调研

为 `clpoly_model.lean` 端 Lean shim 设计提供完整清单。

总计 **20** 个 std::* 符号。

## 容器

| 符号 | 类型出现 | 调用出现 | 合计 | 宿主函数 |
|---|---|---|---|---|
| `std::vector` | 2066 | 0 | 2066 | `__assign_partial_zp`, `__build_cld_matrix`, `__cld_polys` … (共 54) |
| `std::pair` | 685 | 0 | 685 | `__build_cld_matrix`, `__ddf_Zp`, `__edf_Zp` … (共 42) |
| `std::tuple_element` | 235 | 0 | 235 | `__extract_monomial_content`, `__factor_multivar`, `__lll_factorize` … (共 13) |
| `std::map` | 94 | 0 | 94 | `__extract_monomial_content`, `__mtshl_lift`, `__mtshl_sparse_int` … (共 7) |
| `std::list` | 34 | 0 | 34 | `__factor_multivar`, `__select_eval_point`, `__wang_core` |
| `std::set` | 6 | 0 | 6 | `__extract_monomial_content` |

## 算法

| 符号 | 类型出现 | 调用出现 | 合计 | 宿主函数 |
|---|---|---|---|---|
| `std::remove_reference` | 106 | 0 | 106 | `__build_cld_matrix`, `__cld_polys`, `__ddf_Zp` … (共 22) |

## 随机

| 符号 | 类型出现 | 调用出现 | 合计 | 宿主函数 |
|---|---|---|---|---|
| `std::mt19937` | 27 | 0 | 27 | `__edf_Zp`, `__factor_Zp`, `__mtshl_sparse_int` … (共 5) |
| `std::uniform_int_distribution` | 6 | 0 | 6 | `__mtshl_sparse_int`, `__upoly_random` |
| `std::random_device` | 3 | 0 | 3 | `__mtshl_sparse_int` |

## 工具

| 符号 | 类型出现 | 调用出现 | 合计 | 宿主函数 |
|---|---|---|---|---|
| `std::size_t` | 12 | 0 | 12 | `__mtshl_wmds`, `__si_theta_array_eval` |
| `std::initializer_list` | 3 | 0 | 3 | `__ddf_Zp`, `__mtshl_zp_univar_mdp`, `__upoly_powmod` |

## 未分类

| 符号 | 类型出现 | 调用出现 | 合计 | 宿主 |
|---|---|---|---|---|
| `std::_Rb_tree_iterator` | 10 | 0 | 10 | `__extract_monomial_content`, `__mtshl_sparse_int`, `__si_theta_array_eval` |
| `std::_Bit_reference` | 4 | 0 | 4 | `__vanhoeij_recombine`, `__wang_core` |
| `std::mersenne_twister_engine` | 4 | 0 | 4 | `__mtshl_sparse_int`, `__upoly_random` |
| `std::vec` | 2 | 0 | 2 | `__taylor_coeff_zp`, `__wang_core` |
| `std::vect` | 2 | 0 | 2 | `__mtshl_lift`, `__mtshl_step_j` |
| `std::_List_iterator` | 1 | 0 | 1 | `__select_eval_point` |
| `std::v` | 1 | 0 | 1 | `__mtshl_lift` |
| `std::_Rb_tree_const_iterator` | 1 | 0 | 1 | `__mtshl_lift` |

## Lean shim 设计建议

对每类 STL 符号，推荐的 Lean 端实现方式：

| 类别 | Lean 对应 | 备注 |
|---|---|---|
| `std::vector<T>` | `Array T` | 直接等价；`.push_back` → `.push`；`.erase(it, end)` → `.take n` |
| `std::pair<A,B>` | `A × B` | 原生 product type |
| `std::tuple<A,B,C>` | `A × B × C` / 结构体 | 三元 product 或自定义 |
| `std::map<K,V>` | `StdMap K V`（自定义）| 按 key 排序的 Array；CLPoly 用它做线性探针 |
| `std::sort` + comparator | `Array.qsortWith` | lean 4.x 有 `Array.qsort` |
| `std::iota(b, e, v)` | `Array.range` | 生成 [v, v+1, ..., ) |
| `std::move` | `id` | Lean 值语义，move 是 no-op |
| `std::swap(a, b)` | let-rebind | Lean 无 in-place swap，用 `(b, a) := (a, b)` |
| `std::max`/`std::min` | `max` / `min` | Mathlib/stdlib |
| `std::mt19937` + `uniform_int_distribution` | `Rng` 结构体 + `Rng.next` | 公理化（EDF 已用）|
| `std::random_device` | `axiom Rng.seed : Nat` | 仅作 seed 源 |
| `std::log` | `Float.log` / `Nat.log` | 依上下文选具体数值类型 |
| `std::ceil` / `std::floor` | `Nat.ceilDiv` / 截断 | 注意 double 语义 |
