# AST Kind 直方图

> 生成时间：2026-05-12 01:45
> 覆盖 70 个函数

按出现次数降序。每 kind 列首次出现的函数名（样例）。

| Kind | Count | First Seen In |
|---|---|---|
| DeclRefExpr | 6534 | __assign_partial_zp |
| ImplicitCastExpr | 5548 | __assign_partial_zp |
| CXXOperatorCallExpr | 1366 | __assign_partial_zp |
| MemberExpr | 1230 | __assign_partial_zp |
| VarDecl | 1181 | __assign_partial_zp |
| DeclStmt | 1156 | __assign_partial_zp |
| CXXMemberCallExpr | 911 | __assign_partial_zp |
| CXXConstructExpr | 654 | __assign_partial_zp |
| IntegerLiteral | 608 | __assign_partial_zp |
| CXXBindTemporaryExpr | 560 | __assign_partial_zp |
| ExprWithCleanups | 552 | __assign_partial_zp |
| BinaryOperator | 543 | __assign_partial_zp |
| MaterializeTemporaryExpr | 542 | __assign_partial_zp |
| CallExpr | 523 | __assign_partial_zp |
| CompoundStmt | 397 | __assign_partial_zp |
| IfStmt | 288 | __binomial |
| ParmVarDecl | 278 | __assign_partial_zp |
| UnaryOperator | 268 | __assign_partial_zp |
| ReturnStmt | 176 | __assign_partial_zp |
| RecordType | 168 | __assign_partial_zp |
| CXXFunctionalCastExpr | 157 | __binomial |
| TemplateArgument | 148 | __assign_partial_zp |
| ForStmt | 147 | __assign_partial_zp |
| ElaboratedType | 99 | __assign_partial_zp |
| CXXForRangeStmt | 95 | __build_cld_matrix |
| CXXBoolLiteralExpr | 82 | __edf_Zp |
| TemplateSpecializationType | 71 | __assign_partial_zp |
| FunctionDecl | 70 | __assign_partial_zp |
| CStyleCastExpr | 70 | __build_cld_matrix |
| CXXDefaultArgExpr | 70 | __build_cld_matrix |
| SubstTemplateTypeParmType | 68 | __assign_partial_zp |
| ParenExpr | 62 | __build_cld_matrix |
| BindingDecl | 60 | __extract_monomial_content |
| ContinueStmt | 52 | __edf_Zp |
| StringLiteral | 48 | __cld_polys |
| SourceLocExpr | 48 | __cld_polys |
| CXXTemporaryObjectExpr | 46 | __lll_reduce |
| CXXMethodDecl | 41 | __build_cld_matrix |
| ConditionalOperator | 40 | __build_cld_matrix |
| DecompositionDecl | 30 | __extract_monomial_content |
| CXXStaticCastExpr | 27 | __cld_polys |
| BreakStmt | 27 | __ddf_Zp |
| FieldDecl | 27 | __lll_reduce |
| LambdaExpr | 24 | __build_cld_matrix |
| CXXRecordDecl | 24 | __build_cld_matrix |
| CXXDestructorDecl | 24 | __build_cld_matrix |
| PredefinedExpr | 24 | __cld_polys |
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
{"id": "0x7ada75153500", "kind": "DoStmt", "range": {"line": 2375, "col": 29}}
```

### `FloatingLiteral` （3 次，首现 `__heuristic_starting_precision`）

```json
{"id": "0x7eaf1dc07ef0", "kind": "FloatingLiteral", "range": {"line": null, "col": 14}, "type": {"qualType": "double"}, "valueCategory": "prvalue", "value": "2.5"}
```
