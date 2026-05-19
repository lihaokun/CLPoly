# AST Kind 直方图

> 生成时间：2026-05-19 01:18
> 覆盖 5 个函数

按出现次数降序。每 kind 列首次出现的函数名（样例）。

| Kind | Count | First Seen In |
|---|---|---|
| DeclRefExpr | 284 | __hensel_step |
| ImplicitCastExpr | 193 | __hensel_step |
| MemberExpr | 101 | __hensel_step |
| CXXOperatorCallExpr | 67 | __hensel_step |
| VarDecl | 56 | __hensel_step |
| CXXMemberCallExpr | 55 | __hensel_step |
| DeclStmt | 53 | __hensel_step |
| ExprWithCleanups | 28 | __hensel_step |
| CXXBindTemporaryExpr | 25 | __hensel_step |
| CallExpr | 24 | __hensel_step |
| CXXConstructExpr | 23 | __hensel_step |
| MaterializeTemporaryExpr | 14 | __hensel_step |
| CompoundStmt | 12 | __hensel_step |
| ParmVarDecl | 11 | __hensel_step |
| RecordType | 8 | factorize |
| CXXForRangeStmt | 7 | __hensel_step |
| TemplateArgument | 7 | factorize |
| IfStmt | 6 | __hensel_step |
| FunctionDecl | 5 | __hensel_step |
| IntegerLiteral | 5 | __hensel_step |
| ReturnStmt | 5 | __upoly_const_term |
| CXXFunctionalCastExpr | 4 | __hensel_step |
| ElaboratedType | 4 | factorize |
| SubstTemplateTypeParmType | 4 | factorize |
| TemplateSpecializationType | 3 | factorize |
| ForStmt | 2 | __hensel_step |
| BindingDecl | 2 | factorize |
| BinaryOperator | 1 | __upoly_const_term |
| TypeAliasDecl | 1 | factorize |
| DecompositionDecl | 1 | factorize |

## 罕见 Kind 样例（≤3 次出现）

### `BinaryOperator` （1 次，首现 `__upoly_const_term`）

```json
{"id": "0x7a18a461fa28", "kind": "BinaryOperator", "range": {"line": null, "col": 13}, "type": {"qualType": "bool"}, "valueCategory": "prvalue", "opcode": "=="}
```

### `BindingDecl` （2 次，首现 `factorize`）

```json
{"id": "0x7c38d741ed60", "kind": "BindingDecl", "range": {"line": null, "col": 21}, "isReferenced": true, "name": "fac"}
```

### `DecompositionDecl` （1 次，首现 `factorize`）

```json
{"id": "0x7c38d741edf0", "kind": "DecompositionDecl", "range": {"line": null, "col": 14}, "isUsed": true, "type": {"qualType": "std::pair<clpoly::basic_polynomial<clpoly::basic_monomial<clpoly::lex_<clpoly::less>>, clpoly::ZZ, clpoly::lex_<clpoly::less>>, unsigned long> &"}, "init": "c"}
```

### `ForStmt` （2 次，首现 `__hensel_step`）

```json
{"id": "0x7667f358a380", "kind": "ForStmt", "range": {"line": 425, "col": 9}}
```

### `TemplateSpecializationType` （3 次，首现 `factorize`）

```json
{"id": "0x7c38d74158c0", "kind": "TemplateSpecializationType", "type": {"qualType": "polynomial_<ZZ, clpoly::grlex_<clpoly::less>>"}, "isAlias": true, "templateName": "polynomial_"}
```

### `TypeAliasDecl` （1 次，首现 `factorize`）

```json
{"id": "0x7c38d7415988", "kind": "TypeAliasDecl", "range": {"line": null, "col": 9}, "isReferenced": true, "name": "Poly", "type": {"desugaredQualType": "clpoly::basic_polynomial<clpoly::basic_monomial<clpoly::grlex_<clpoly::less>>, clpoly::ZZ, clpoly::grlex_<clpoly::less>>", "qualType": "polynomial_<ZZ, grlex_<less>>"}}
```
