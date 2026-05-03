# cpp2lean v2 翻译器 1 对 1 人工审视主索引

**目标**：对 346 个翻译函数（64 顶层 + 247 loops + 32 lambdas + 3 其他）逐个对照 C++ 源代码做白盒审视，发现 lake build 通不到的语义偏差（silent miss）。

**起因**：Stage F-G 修了大量 latent bug 后，lake build 100% 通过 ≠ 翻译忠实。需要白盒检查每个函数的控制流 / 数据流 / 调用链是否 1-对-1 对应 C++ 源。

**工作流**：参考 user memory `feedback_pass_audit_loop.md`，每个 Pass / 每个函数完成后必须做语义层 1-对-1 审视。

## 工程量与分批

| 类别 | 数量 | 优先级 | 状态 |
|------|------|-------|------|
| 顶层函数（`__xxx_ir`）| 64 | P0 | **64/64 ✅ 全部完成** |
| Loop 辅助（`_loop_*_ir`）| 247 | P1 顺带 | - |
| Lambda 辅助（`_lambda_*_ir`）| 32 | P1 顺带 | - |
| 其他（factorize / partial_def）| 3 | P0 | - |

**策略**：按顶层函数 + 它依赖的 loops/lambdas 一起审视（共 64 批次）。

## 审视顺序（按依赖深度）

按 `back2back-design.md` 已有的 L0-L4 分层 + 依赖图：

### L0 原子（无 stub 依赖）— 1 批

- `__make_zp` ✅ 待审

### L1 单变量基础 — 8 批

- `__upoly_make_monic`
- `__upoly_mod`
- `__upoly_divmod`
- `__upoly_powmod`
- `__upoly_random`
- `__upoly_subtract_x`
- `__upoly_subtract_one`
- `__upoly_const_term`
- `__upoly_norm_l1`
- `__upoly_norm_l2_sq`
- `__upoly_primitive`
- `__upoly_to_poly`
- `__upoly_mod_coeff`
- `__upoly_divmod_mod`
- `__upoly_mul_mod`
- `__upoly_symmetric_mod`

### L2 单变量算法 — 14 批

- `__extract_pth_root`
- `__squarefree_Zp`
- `__symmetric_mod`
- `__binomial`
- `__isqrt_ceil`
- `__mignotte_bound`
- `__heuristic_starting_precision`
- `__select_prime`
- `__hensel_step`
- `__hensel_step_linear`
- `__hensel_tree_build`
- `__hensel_tree_build_recursive`
- `__hensel_lift`
- `__hensel_lift_recursive`
- `__hensel_lift_linear_recursive`
- `__hensel_extract_factors`
- `__cld_polys`
- `__build_cld_matrix`
- `__lll_reduce`
- `__lll_factorize`
- `__extract_candidates`
- `__zassenhaus_recombine`
- `__vanhoeij_recombine`
- `__factor_recombine`
- `__factor_squarefree_primitive_ZZ`
- `__ddf_Zp`
- `__edf_Zp`
- `__factor_Zp`

### L3 单变量顶层 — 1 批

- `factorize` (univar 实例)
- `factorize_grlex` / `factorize_lex` / `factorize_upoly`

### L4 多变量原语 — 12 批

- `__assign_partial_zp`
- `__polynomial_to_zp`
- `__symmetric_mod_poly`
- `__taylor_coeff_zp`
- `__select_eval_point`
- `__si_theta_array_eval`
- `__si_vandermonde_solve`
- `__extract_monomial_content`
- `__mtshl_coeff_bound`
- `__mtshl_lift`
- `__mtshl_multi_bdp`
- `__mtshl_sparse_int`
- `__mtshl_step_j`
- `__mtshl_wmds`
- `__mtshl_zp_univar_mdp`

### L5 Wang 顶层 — 4 批

- `__wang_core`
- `__wang_leading_coeff`
- `__factor_multivar`

### 其他散点 — 4 批

- `factorize`（Lean 顶层 entry）
- 等

## 审视模板

每个函数审视项：

```markdown
## ✅ / ⚠️ / ❌ `__function_name_<inst>_ir`

**C++ 源**：`clpoly/file.hh:LINE_FROM-LINE_TO`
**Lean 翻译**：`Generated/Corpus.lean:LINE_FROM-LINE_TO`
**关联辅助函数**：`_loop_..._N_ir`, `_lambda_..._M_ir` ...

### 签名对照

| 项 | C++ | Lean | 一致 |
|----|-----|------|------|
| 参数 | `(a : T1, b : T2)` | `(a_ir : T1, b_ir : T2)` | ✅/❌ |
| 返回 | `T3` | `T3` | ✅/❌ |

### 控制流对照

| C++ 行 | Lean 行 | 翻译 | 一致 |
|--------|--------|------|------|
| ... | ... | for→loop | ✅ |

### 调用对照

| C++ 调用 | Lean 翻译 | 原语类型 | 备注 |
|---------|----------|---------|------|
| `polynomial_GCD(a, b)` | `polynomial_GCD a b` | stub→default | ⚠️ 待实现真实算法 |

### 偏差 / 已修 / 待修

- ✅ 无偏差
- ⚠️ 待修：xxx
- ❌ 翻译错：xxx，已 commit ABC
```

## 进度追踪

| 批次 | 函数 | 状态 | commit |
|------|------|------|--------|
| - | - | - | - |

（每完成一批更新本表）

## 输出文档

每批审视输出独立 .md：
- `audit/L0-make_zp.md`
- `audit/L1-upoly-basic.md`（合并多个 upoly_xxx）
- `audit/L2-squarefree.md`
- ...

最终审视报告汇总：`audit/SUMMARY.md`。
