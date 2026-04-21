# AST 调研总览

> 生成时间：2026-04-21 10:12

## 覆盖情况

- 目标函数总数（TRANSLATION_SCOPE）: **65**
- 成功 dump AST 的函数数: **65**
- 失败/缺失的函数数: **0**

## 全局统计

- 不同的 AST kind 种数: **58**
- 不同的运算符种数: **45**
- 不同的类型（qualType）种数: **1199**
- AST 节点总数（累计）: **24560**

## 潜在难点提示

### 控制流构造
- ✅ `IfStmt`: 275
- ✅ `WhileStmt`: 20
- ✅ `ForStmt`: 145
- ✅ `CXXForRangeStmt`: 92
- ✅ `DoStmt`: 2
- ✅ `BreakStmt`: 26
- ✅ `ContinueStmt`: 52
- ✅ `ReturnStmt`: 165
- — `SwitchStmt`: 0
- — `CaseStmt`: 0
- — `DefaultStmt`: 0
- — `CXXTryStmt`: 0
- — `CXXThrowExpr`: 0
- — `CXXCatchStmt`: 0
- — `GotoStmt`: 0

### Lambda / 闭包
- ✅ `LambdaExpr`: 24
- ✅ `CXXRecordDecl`: 24
- ✅ `CallExpr`: 507
- ✅ `CXXOperatorCallExpr`: 1332

### 结构化绑定 / 迭代器
- ✅ `DecompositionDecl`: 30
- ✅ `BindingDecl`: 60
- ✅ `CXXBindTemporaryExpr`: 542
- ✅ `MaterializeTemporaryExpr`: 522
- ✅ `CXXConstructExpr`: 627
- — `CXXDependentScopeMemberExpr`: 0

### 模板相关
- — `FunctionTemplateDecl`: 0
- — `ClassTemplateDecl`: 0
- — `ClassTemplateSpecializationDecl`: 0
- — `TemplateTypeParmDecl`: 0
- — `UnresolvedLookupExpr`: 0
- — `CXXDependentScopeMemberExpr`: 0
- — `DependentScopeDeclRefExpr`: 0

### 内存/生命周期
- — `CXXNewExpr`: 0
- — `CXXDeleteExpr`: 0
- ✅ `ExprWithCleanups`: 535
- — `CXXThisExpr`: 0
- — `UnresolvedMemberExpr`: 0

## 每函数使用的 kind 数量

| Function | Distinct Kinds |
|---|---|
| __wang_leading_coeff | 52 |
| __wang_core | 46 |
| __factor_multivar | 45 |
| __mtshl_sparse_int | 43 |
| __mtshl_wmds | 42 |
| __select_eval_point | 41 |
| __vanhoeij_recombine | 41 |
| __mtshl_multi_bdp | 40 |
| __mtshl_step_j | 40 |
| __zassenhaus_recombine | 40 |
| __select_prime | 38 |
| __lll_reduce | 37 |
| __mtshl_lift | 37 |
| __si_theta_array_eval | 37 |
| __mtshl_zp_univar_mdp | 35 |
| __factor_Zp | 34 |
| __si_vandermonde_solve | 34 |
| __extract_monomial_content | 33 |
| __build_cld_matrix | 32 |
| __ddf_Zp | 32 |
| __hensel_lift | 31 |
| __factor_squarefree_primitive_ZZ | 29 |
| __squarefree_Zp | 29 |
| __upoly_divmod_mod | 29 |
| __upoly_powmod | 29 |
| __hensel_tree_build | 28 |
| __taylor_coeff_zp | 28 |
| __edf_Zp | 27 |
| __symmetric_mod_poly | 27 |
| __assign_partial_zp | 26 |
| __cld_polys | 26 |
| __extract_candidates | 26 |
| __extract_pth_root | 26 |
| __heuristic_starting_precision | 25 |
| __mtshl_coeff_bound | 25 |
| __polynomial_to_zp | 25 |
| __upoly_make_monic | 25 |
| factorize | 25 |
| __mignotte_bound | 24 |
| __upoly_random | 24 |
| __upoly_subtract_x | 24 |
| __hensel_tree_build_recursive | 23 |
| __isqrt_ceil | 23 |
| __upoly_subtract_one | 23 |
| __hensel_step | 20 |
| __lll_factorize | 20 |
| __upoly_primitive | 20 |
| __binomial | 19 |
| __hensel_step_linear | 19 |
| __subset_product_mod | 19 |
| __upoly_norm_l1 | 19 |
| __upoly_symmetric_mod | 18 |
| __symmetric_mod | 17 |
| __upoly_norm_l2_sq | 17 |
| __upoly_mod_coeff | 16 |
| __hensel_extract_factors | 15 |
| __upoly_const_term | 15 |
| __factor_recombine | 14 |
| __upoly_mul_mod | 13 |
| __hensel_lift_recursive | 12 |
| __upoly_mod | 12 |
| __upoly_to_poly | 12 |
| __hensel_lift_linear_recursive | 11 |
| __upoly_divmod | 8 |
| __make_zp | 7 |