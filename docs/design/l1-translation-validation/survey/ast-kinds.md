# AST Kind 直方图

> 生成时间：2026-04-21 11:07
> 覆盖 65 个函数

按出现次数降序。每 kind 列首次出现的函数名（样例）。

| Kind | Count | First Seen In |
|---|---|---|
| DeclRefExpr | 6354 | __assign_partial_zp |
| ImplicitCastExpr | 5419 | __assign_partial_zp |
| CXXOperatorCallExpr | 1332 | __assign_partial_zp |
| MemberExpr | 1177 | __assign_partial_zp |
| VarDecl | 1145 | __assign_partial_zp |
| DeclStmt | 1123 | __assign_partial_zp |
| CXXMemberCallExpr | 869 | __assign_partial_zp |
| CXXConstructExpr | 627 | __assign_partial_zp |
| IntegerLiteral | 596 | __assign_partial_zp |
| CXXBindTemporaryExpr | 542 | __assign_partial_zp |
| ExprWithCleanups | 535 | __assign_partial_zp |
| BinaryOperator | 535 | __assign_partial_zp |
| MaterializeTemporaryExpr | 522 | __assign_partial_zp |
| CallExpr | 507 | __assign_partial_zp |
| CompoundStmt | 381 | __assign_partial_zp |
| IfStmt | 275 | __binomial |
| ParmVarDecl | 267 | __assign_partial_zp |
| UnaryOperator | 260 | __assign_partial_zp |
| RecordType | 165 | __assign_partial_zp |
| ReturnStmt | 165 | __assign_partial_zp |
| CXXFunctionalCastExpr | 155 | __binomial |
| TemplateArgument | 145 | __assign_partial_zp |
| ForStmt | 145 | __assign_partial_zp |
| ElaboratedType | 99 | __assign_partial_zp |
| CXXForRangeStmt | 92 | __build_cld_matrix |
| CXXBoolLiteralExpr | 82 | __edf_Zp |
| TemplateSpecializationType | 71 | __assign_partial_zp |
| CStyleCastExpr | 70 | __build_cld_matrix |
| CXXDefaultArgExpr | 69 | __build_cld_matrix |
| SubstTemplateTypeParmType | 68 | __assign_partial_zp |
| FunctionDecl | 65 | __assign_partial_zp |
| ParenExpr | 60 | __build_cld_matrix |
| BindingDecl | 60 | __extract_monomial_content |
| ContinueStmt | 52 | __edf_Zp |
| StringLiteral | 46 | __cld_polys |
| SourceLocExpr | 46 | __cld_polys |
| CXXTemporaryObjectExpr | 46 | __lll_reduce |
| CXXMethodDecl | 41 | __build_cld_matrix |
| ConditionalOperator | 39 | __build_cld_matrix |
| DecompositionDecl | 30 | __extract_monomial_content |
| FieldDecl | 27 | __lll_reduce |
| CXXStaticCastExpr | 26 | __cld_polys |
| BreakStmt | 26 | __ddf_Zp |
| LambdaExpr | 24 | __build_cld_matrix |
| CXXRecordDecl | 24 | __build_cld_matrix |
| CXXDestructorDecl | 24 | __build_cld_matrix |
| PredefinedExpr | 23 | __cld_polys |
| TypeAliasDecl | 20 | __assign_partial_zp |
| WhileStmt | 20 | __edf_Zp |
| CXXConstructorDecl | 12 | __factor_Zp |
| ImplicitValueInitExpr | 12 | __hensel_tree_build |
| CXXConversionDecl | 11 | __build_cld_matrix |
| InitListExpr | 11 | __ddf_Zp |
| CXXStdInitializerListExpr | 8 | __ddf_Zp |
| CompoundAssignOperator | 6 | __select_eval_point |
| CXXCtorInitializer | 4 | __lll_reduce |
| FloatingLiteral | 3 | __heuristic_starting_precision |
| DoStmt | 2 | __wang_core |

## 罕见 Kind 样例（≤3 次出现）

### `DoStmt` （2 次，首现 `__wang_core`）

```json
{"id": "0x7abbf71ccde0", "kind": "DoStmt", "range": {"line": 2375, "col": 29}}
```

### `FloatingLiteral` （3 次，首现 `__heuristic_starting_precision`）

```json
{"id": "0x74e61b3f1a40", "kind": "FloatingLiteral", "range": {"line": null, "col": 14}, "type": {"qualType": "double"}, "valueCategory": "prvalue", "value": "2.5"}
```
