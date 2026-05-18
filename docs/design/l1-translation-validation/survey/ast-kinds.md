# AST Kind 直方图

> 生成时间：2026-05-19 01:11
> 覆盖 3 个函数

按出现次数降序。每 kind 列首次出现的函数名（样例）。

| Kind | Count | First Seen In |
|---|---|---|
| MemberExpr | 16 | __upoly_const_term |
| DeclRefExpr | 16 | __upoly_const_term |
| CXXMemberCallExpr | 14 | __upoly_const_term |
| ParmVarDecl | 7 | __upoly_const_term |
| CXXConstructExpr | 6 | __upoly_const_term |
| ReturnStmt | 4 | __upoly_const_term |
| ImplicitCastExpr | 4 | __upoly_const_term |
| FunctionDecl | 3 | __upoly_const_term |
| CompoundStmt | 3 | __upoly_const_term |
| IntegerLiteral | 3 | __upoly_const_term |
| IfStmt | 2 | __upoly_const_term |
| ExprWithCleanups | 2 | __upoly_const_term |
| CXXFunctionalCastExpr | 2 | __upoly_const_term |
| CXXBindTemporaryExpr | 2 | __upoly_const_term |
| CallExpr | 2 | __upoly_divmod |
| VarDecl | 2 | __upoly_mod |
| BinaryOperator | 1 | __upoly_const_term |
| DeclStmt | 1 | __upoly_mod |

## 罕见 Kind 样例（≤3 次出现）

### `BinaryOperator` （1 次，首现 `__upoly_const_term`）

```json
{"id": "0x7f74c89eba28", "kind": "BinaryOperator", "range": {"line": null, "col": 13}, "type": {"qualType": "bool"}, "valueCategory": "prvalue", "opcode": "=="}
```

### `CXXBindTemporaryExpr` （2 次，首现 `__upoly_const_term`）

```json
{"id": "0x7f74c89eb698", "kind": "CXXBindTemporaryExpr", "range": {"line": null, "col": 31}, "type": {"desugaredQualType": "clpoly::ZZ", "qualType": "ZZ"}, "valueCategory": "prvalue", "temp": "0x7f74c89eb690", "dtor": {"id": "0x5a88523d8b58", "kind": "CXXDestructorDecl", "name": "~ZZ", "type": {"qualType": "void () noexcept"}}}
```

### `CXXFunctionalCastExpr` （2 次，首现 `__upoly_const_term`）

```json
{"id": "0x7f74c89eb6b8", "kind": "CXXFunctionalCastExpr", "range": {"line": null, "col": 31}, "type": {"desugaredQualType": "clpoly::ZZ", "qualType": "ZZ"}, "valueCategory": "prvalue", "castKind": "ConstructorConversion", "conversionFunc": {"id": "0x5a88523d76e8", "kind": "CXXConstructorDecl", "name": "ZZ", "type": {"qualType": "void (int)"}}}
```

### `CallExpr` （2 次，首现 `__upoly_divmod`）

```json
{"id": "0x7ff3fc146d18", "kind": "CallExpr", "range": {"line": 55, "col": 9}, "type": {"qualType": "void"}, "valueCategory": "prvalue"}
```

### `CompoundStmt` （3 次，首现 `__upoly_const_term`）

```json
{"id": "0x7f74c89ebcd0", "kind": "CompoundStmt", "range": {"line": 744, "col": 5}}
```

### `DeclStmt` （1 次，首现 `__upoly_mod`）

```json
{"id": "0x7a43ebf000e0", "kind": "DeclStmt", "range": {"line": 43, "col": 9}}
```

### `ExprWithCleanups` （2 次，首现 `__upoly_const_term`）

```json
{"id": "0x7f74c89eb6e0", "kind": "ExprWithCleanups", "range": {"line": null, "col": 31}, "type": {"desugaredQualType": "clpoly::ZZ", "qualType": "ZZ"}, "valueCategory": "prvalue", "cleanupsHaveSideEffects": true}
```

### `FunctionDecl` （3 次，首现 `__upoly_const_term`）

```json
{"id": "0x7f74c89eb490", "kind": "FunctionDecl", "range": {"line": null, "col": 5}, "isUsed": true, "name": "__upoly_const_term", "mangledName": "_ZN6clpoly18__upoly_const_termERKNS_16basic_polynomialINS_9umonomialENS_2ZZENS_5ulessEEE", "type": {"qualType": "ZZ (const upolynomial_<ZZ> &)"}, "inline": true}
```

### `IfStmt` （2 次，首现 `__upoly_const_term`）

```json
{"id": "0x7f74c89eb708", "kind": "IfStmt", "range": {"line": 745, "col": 9}}
```

### `IntegerLiteral` （3 次，首现 `__upoly_const_term`）

```json
{"id": "0x7f74c89eb5f0", "kind": "IntegerLiteral", "range": {"line": null, "col": 34}, "type": {"qualType": "int"}, "valueCategory": "prvalue", "value": "0"}
```

### `VarDecl` （2 次，首现 `__upoly_mod`）

```json
{"id": "0x7a43ebefff50", "kind": "VarDecl", "range": {"line": null, "col": 9}, "isUsed": true, "name": "q", "type": {"desugaredQualType": "clpoly::basic_polynomial<clpoly::umonomial, clpoly::Zp, clpoly::uless>", "qualType": "upolynomial_<Zp>"}, "init": "call"}
```
