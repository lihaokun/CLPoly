# 运算符直方图

> 生成时间：2026-05-11 18:39

包括 `CXXOperatorCallExpr`（用户类型重载）、`UnaryOperator`（基本类型 `!` `-` `~` `++` 等）、`BinaryOperator`（`+ - * / == ` 等）、`CompoundAssignOperator`（`+= -= *=` 等）。

| Operator | Count | First Seen In |
|---|---|---|
| CXXOperatorCallExpr::operator[] | 468 | __assign_partial_zp |
| CXXOperatorCallExpr::operator* | 207 | __build_cld_matrix |
| CXXOperatorCallExpr::operator= | 173 | __assign_partial_zp |
| BinaryOperator::< | 146 | __assign_partial_zp |
| UnaryOperator::++ | 145 | __assign_partial_zp |
| CXXOperatorCallExpr::operator!= | 131 | __build_cld_matrix |
| CXXOperatorCallExpr::operator++ | 110 | __build_cld_matrix |
| BinaryOperator::= | 74 | __binomial |
| UnaryOperator::! | 70 | __build_cld_matrix |
| BinaryOperator::- | 64 | __binomial |
| BinaryOperator::== | 40 | __binomial |
| BinaryOperator::+ | 40 | __binomial |
| CXXOperatorCallExpr::operator() | 40 | __build_cld_matrix |
| CXXOperatorCallExpr::operator- | 40 | __edf_Zp |
| BinaryOperator::> | 39 | __binomial |
| BinaryOperator::&& | 36 | __build_cld_matrix |
| CXXOperatorCallExpr::operator+ | 30 | __edf_Zp |
| CXXOperatorCallExpr::operator-> | 29 | __extract_monomial_content |
| BinaryOperator::>= | 24 | __factor_squarefree_primitive_ZZ |
| CXXOperatorCallExpr::operator*= | 23 | __binomial |
| UnaryOperator::__extension__ | 23 | __cld_polys |
| CXXOperatorCallExpr::operator== | 23 | __extract_monomial_content |
| BinaryOperator::<= | 19 | __edf_Zp |
| BinaryOperator::!= | 18 | __extract_candidates |
| BinaryOperator::* | 16 | __ddf_Zp |
| CXXOperatorCallExpr::operator/ | 15 | __edf_Zp |
| UnaryOperator::- | 12 | __extract_candidates |
| CXXOperatorCallExpr::operator< | 12 | __factor_multivar |
| CXXOperatorCallExpr::operator-= | 12 | __hensel_lift |
| UnaryOperator::-- | 12 | __lll_reduce |
| BinaryOperator::/ | 11 | __build_cld_matrix |
| BinaryOperator::|| | 9 | __binomial |
| CXXOperatorCallExpr::operator% | 8 | __upoly_divmod_mod |
| CXXOperatorCallExpr::operator<= | 7 | __hensel_lift |
| BinaryOperator::% | 5 | __build_cld_matrix |
| CXXOperatorCallExpr::operator+= | 5 | __isqrt_ceil |
| CXXOperatorCallExpr::operator> | 5 | __mtshl_coeff_bound |
| CXXOperatorCallExpr::operator<< | 4 | __vanhoeij_recombine |
| CXXOperatorCallExpr::operator/= | 2 | __binomial |
| CXXOperatorCallExpr::operator>= | 2 | __isqrt_ceil |
| CompoundAssignOperator::*= | 2 | __select_eval_point |
| CompoundAssignOperator::+= | 2 | __vanhoeij_recombine |
| CXXOperatorCallExpr::operator<<= | 1 | __isqrt_ceil |
| CompoundAssignOperator::/= | 1 | __select_eval_point |
| CompoundAssignOperator::-= | 1 | __vanhoeij_recombine |