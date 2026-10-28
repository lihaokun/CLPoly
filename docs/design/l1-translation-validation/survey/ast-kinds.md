# AST Kind 直方图

> 生成时间：2026-05-24 20:33
> 覆盖 66 个函数

按出现次数降序。每 kind 列首次出现的函数名（样例）。

| Kind | Count | First Seen In |
|---|---|---|
| DeclRefExpr | 6442 | __assign_partial_zp |
| ImplicitCastExpr | 5581 | __assign_partial_zp |
| CXXOperatorCallExpr | 1347 | __assign_partial_zp |
| MemberExpr | 1183 | __assign_partial_zp |
| VarDecl | 1152 | __assign_partial_zp |
| DeclStmt | 1128 | __assign_partial_zp |
| CXXMemberCallExpr | 878 | __assign_partial_zp |
| IntegerLiteral | 650 | __assign_partial_zp |
| CXXConstructExpr | 639 | __assign_partial_zp |
| CXXBindTemporaryExpr | 558 | __assign_partial_zp |
| BinaryOperator | 541 | __assign_partial_zp |
| CallExpr | 538 | __assign_partial_zp |
| MaterializeTemporaryExpr | 536 | __assign_partial_zp |
| ExprWithCleanups | 530 | __assign_partial_zp |
| CompoundStmt | 376 | __assign_partial_zp |
| IfStmt | 282 | __binomial |
| ParmVarDecl | 262 | __assign_partial_zp |
| UnaryOperator | 262 | __assign_partial_zp |
| RecordType | 196 | __assign_partial_zp |
| TemplateArgument | 176 | __assign_partial_zp |
| ReturnStmt | 168 | __assign_partial_zp |
| ElaboratedType | 149 | __assign_partial_zp |
| ForStmt | 146 | __assign_partial_zp |
| CXXFunctionalCastExpr | 133 | __binomial |
| TemplateSpecializationType | 101 | __assign_partial_zp |
| SubstTemplateTypeParmType | 100 | __assign_partial_zp |
| CStyleCastExpr | 93 | __build_cld_matrix |
| CXXForRangeStmt | 92 | __build_cld_matrix |
| ParenExpr | 83 | __build_cld_matrix |
| CXXBoolLiteralExpr | 82 | __edf_Zp |
| StringLiteral | 69 | __cld_polys |
| FunctionDecl | 66 | __assign_partial_zp |
| BindingDecl | 60 | __extract_monomial_content |
| ContinueStmt | 52 | __edf_Zp |
| CXXTemporaryObjectExpr | 46 | __lll_reduce |
| ConditionalOperator | 39 | __build_cld_matrix |
| CXXMethodDecl | 35 | __build_cld_matrix |
| DecompositionDecl | 30 | __extract_monomial_content |
| BreakStmt | 27 | __ddf_Zp |
| FieldDecl | 27 | __lll_reduce |
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
| FloatingLiteral | 3 | __heuristic_starting_precision |
| DoStmt | 2 | __wang_core |
| CXXStaticCastExpr | 1 | __factor_multivar |
| CXXDefaultArgExpr | 1 | __lll_factorize |

## 罕见 Kind 样例（≤3 次出现）

### `CXXDefaultArgExpr` （1 次，首现 `__lll_factorize`）

```json
{"id": "0x879664030", "kind": "CXXDefaultArgExpr", "range": {"line": null, "col": null}, "type": {"qualType": "int"}, "valueCategory": "prvalue"}
```

### `CXXStaticCastExpr` （1 次，首现 `__factor_multivar`）

```json
{"id": "0x7ae9092e8", "kind": "CXXStaticCastExpr", "range": {"line": null, "col": 68}, "type": {"desugaredQualType": "unsigned long long", "qualType": "uint64_t", "typeAliasDeclId": "0x7b6f9d048"}, "valueCategory": "prvalue", "castKind": "NoOp"}
```

### `DoStmt` （2 次，首现 `__wang_core`）

```json
{"id": "0x892ec6240", "kind": "DoStmt", "range": {"line": 2375, "col": 29}}
```

### `FloatingLiteral` （3 次，首现 `__heuristic_starting_precision`）

```json
{"id": "0x7c43260b8", "kind": "FloatingLiteral", "range": {"line": null, "col": 14}, "type": {"qualType": "double"}, "valueCategory": "prvalue", "value": "2.5"}
```
