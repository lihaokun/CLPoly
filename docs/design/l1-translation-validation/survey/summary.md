# AST 调研总览

> 生成时间：2026-05-24 20:33

## 覆盖情况

- 目标函数总数（TRANSLATION_SCOPE）: **66**
- 成功 dump AST 的函数数: **66**
- 失败/缺失的函数数: **0**

## 全局统计

- 不同的 AST kind 种数: **56**
- 不同的运算符种数: **44**
- 不同的类型（qualType）种数: **1199**
- AST 节点总数（累计）: **25059**

## 潜在难点提示

### 控制流构造
- ✅ `IfStmt`: 282
- ✅ `WhileStmt`: 20
- ✅ `ForStmt`: 146
- ✅ `CXXForRangeStmt`: 92
- ✅ `DoStmt`: 2
- ✅ `BreakStmt`: 27
- ✅ `ContinueStmt`: 52
- ✅ `ReturnStmt`: 168
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
- ✅ `CallExpr`: 538
- ✅ `CXXOperatorCallExpr`: 1347

### 结构化绑定 / 迭代器
- ✅ `DecompositionDecl`: 30
- ✅ `BindingDecl`: 60
- ✅ `CXXBindTemporaryExpr`: 558
- ✅ `MaterializeTemporaryExpr`: 536
- ✅ `CXXConstructExpr`: 639
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
- ✅ `ExprWithCleanups`: 530
- — `CXXThisExpr`: 0
- — `UnresolvedMemberExpr`: 0

## 每函数使用的 kind 数量

| Function | Distinct Kinds |
|---|---|
| __wang_leading_coeff | 49 |
| __wang_core | 46 |
| __factor_multivar | 45 |
| __mtshl_sparse_int | 42 |
| __mtshl_wmds | 41 |
| __select_eval_point | 40 |
| __vanhoeij_recombine | 40 |
| __mtshl_multi_bdp | 39 |
| __mtshl_step_j | 39 |
| __zassenhaus_recombine | 39 |
| __select_prime | 37 |
| __mtshl_lift | 36 |
| __lll_reduce | 34 |
| __si_theta_array_eval | 34 |
| __extract_monomial_content | 33 |
| __factor_Zp | 32 |
| __mtshl_zp_univar_mdp | 32 |
| __build_cld_matrix | 31 |
| __ddf_Zp | 30 |
| __si_vandermonde_solve | 30 |
| __hensel_lift | 29 |
| __taylor_coeff_zp | 28 |
| __upoly_powmod | 28 |
| __edf_Zp | 27 |
| __squarefree_Zp | 27 |
| __symmetric_mod_poly | 27 |
| __upoly_divmod_mod | 27 |
| squarefreefactorize | 27 |
| __assign_partial_zp | 26 |
| __factor_squarefree_primitive_ZZ | 26 |
| __extract_candidates | 25 |
| __extract_pth_root | 25 |
| __hensel_tree_build | 25 |
| __heuristic_starting_precision | 25 |
| __mtshl_coeff_bound | 25 |
| __polynomial_to_zp | 25 |
| factorize | 25 |
| __cld_polys | 24 |
| __upoly_random | 24 |
| __upoly_subtract_x | 24 |
| __hensel_tree_build_recursive | 23 |
| __isqrt_ceil | 23 |
| __upoly_make_monic | 23 |
| __upoly_subtract_one | 23 |
| __mignotte_bound | 22 |
| __lll_factorize | 21 |
| __hensel_step | 20 |
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