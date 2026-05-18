# AST 调研总览

> 生成时间：2026-05-19 01:11

## 覆盖情况

- 目标函数总数（TRANSLATION_SCOPE）: **66**
- 成功 dump AST 的函数数: **3**
- 失败/缺失的函数数: **0**

## 全局统计

- 不同的 AST kind 种数: **18**
- 不同的运算符种数: **1**
- 不同的类型（qualType）种数: **23**
- AST 节点总数（累计）: **90**

## 潜在难点提示

### 控制流构造
- ✅ `IfStmt`: 2
- — `WhileStmt`: 0
- — `ForStmt`: 0
- — `CXXForRangeStmt`: 0
- — `DoStmt`: 0
- — `BreakStmt`: 0
- — `ContinueStmt`: 0
- ✅ `ReturnStmt`: 4
- — `SwitchStmt`: 0
- — `CaseStmt`: 0
- — `DefaultStmt`: 0
- — `CXXTryStmt`: 0
- — `CXXThrowExpr`: 0
- — `CXXCatchStmt`: 0
- — `GotoStmt`: 0

### Lambda / 闭包
- — `LambdaExpr`: 0
- — `CXXRecordDecl`: 0
- ✅ `CallExpr`: 2
- — `CXXOperatorCallExpr`: 0

### 结构化绑定 / 迭代器
- — `DecompositionDecl`: 0
- — `BindingDecl`: 0
- ✅ `CXXBindTemporaryExpr`: 2
- — `MaterializeTemporaryExpr`: 0
- ✅ `CXXConstructExpr`: 6
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
- ✅ `ExprWithCleanups`: 2
- — `CXXThisExpr`: 0
- — `UnresolvedMemberExpr`: 0

## 每函数使用的 kind 数量

| Function | Distinct Kinds |
|---|---|
| __upoly_const_term | 15 |
| __upoly_mod | 12 |
| __upoly_divmod | 8 |