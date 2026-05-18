# AST Kind 直方图

> 生成时间：2026-05-19 00:58
> 覆盖 66 个函数

按出现次数降序。每 kind 列首次出现的函数名（样例）。

| Kind | Count | First Seen In |
|---|---|---|
| DeclRefExpr | 6423 | __assign_partial_zp |
| ImplicitCastExpr | 5472 | __assign_partial_zp |
| CXXOperatorCallExpr | 1347 | __assign_partial_zp |
| MemberExpr | 1187 | __assign_partial_zp |
| VarDecl | 1152 | __assign_partial_zp |
| DeclStmt | 1128 | __assign_partial_zp |
| CXXMemberCallExpr | 878 | __assign_partial_zp |
| CXXConstructExpr | 639 | __assign_partial_zp |
| IntegerLiteral | 603 | __assign_partial_zp |
| CXXBindTemporaryExpr | 558 | __assign_partial_zp |
| ExprWithCleanups | 547 | __assign_partial_zp |
| BinaryOperator | 541 | __assign_partial_zp |
| MaterializeTemporaryExpr | 536 | __assign_partial_zp |
| CallExpr | 515 | __assign_partial_zp |
| CompoundStmt | 388 | __assign_partial_zp |
| IfStmt | 282 | __binomial |
| ParmVarDecl | 268 | __assign_partial_zp |
| UnaryOperator | 262 | __assign_partial_zp |
| ReturnStmt | 168 | __assign_partial_zp |
| RecordType | 166 | __assign_partial_zp |
| CXXFunctionalCastExpr | 156 | __binomial |
| TemplateArgument | 146 | __assign_partial_zp |
| ForStmt | 146 | __assign_partial_zp |
| ElaboratedType | 99 | __assign_partial_zp |
| CXXForRangeStmt | 92 | __build_cld_matrix |
| CXXBoolLiteralExpr | 82 | __edf_Zp |
| TemplateSpecializationType | 71 | __assign_partial_zp |
| CStyleCastExpr | 70 | __build_cld_matrix |
| CXXDefaultArgExpr | 70 | __build_cld_matrix |
| SubstTemplateTypeParmType | 68 | __assign_partial_zp |
| FunctionDecl | 66 | __assign_partial_zp |
| ParenExpr | 60 | __build_cld_matrix |
| BindingDecl | 60 | __extract_monomial_content |
| ContinueStmt | 52 | __edf_Zp |
| StringLiteral | 46 | __cld_polys |
| SourceLocExpr | 46 | __cld_polys |
| CXXTemporaryObjectExpr | 46 | __lll_reduce |
| CXXMethodDecl | 41 | __build_cld_matrix |
| ConditionalOperator | 39 | __build_cld_matrix |
| DecompositionDecl | 30 | __extract_monomial_content |
| BreakStmt | 27 | __ddf_Zp |
| FieldDecl | 27 | __lll_reduce |
| CXXStaticCastExpr | 26 | __cld_polys |
| LambdaExpr | 24 | __build_cld_matrix |
| CXXRecordDecl | 24 | __build_cld_matrix |
| CXXDestructorDecl | 24 | __build_cld_matrix |
| PredefinedExpr | 23 | __cld_polys |
| TypeAliasDecl | 20 | __assign_partial_zp |
| WhileStmt | 20 | __edf_Zp |
| InitListExpr | 12 | __ddf_Zp |
| CXXConstructorDecl | 12 | __factor_Zp |
| ImplicitValueInitExpr | 12 | __hensel_tree_build |
| CXXConversionDecl | 11 | __build_cld_matrix |
| CXXStdInitializerListExpr | 9 | __ddf_Zp |
| CompoundAssignOperator | 6 | __select_eval_point |
| CXXCtorInitializer | 4 | __lll_reduce |
| FloatingLiteral | 3 | __heuristic_starting_precision |
| DoStmt | 2 | __wang_core |

## 罕见 Kind 样例（≤3 次出现）

### `DoStmt` （2 次，首现 `__wang_core`）

```json
{"id": "0x7be5fcc7bac0", "kind": "DoStmt", "range": {"line": 2375, "col": 29}}
```

### `FloatingLiteral` （3 次，首现 `__heuristic_starting_precision`）

```json
{"id": "0x7e2e3491bc40", "kind": "FloatingLiteral", "range": {"line": null, "col": 14}, "type": {"qualType": "double"}, "valueCategory": "prvalue", "value": "2.5"}
```
