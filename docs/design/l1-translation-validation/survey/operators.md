# 运算符直方图

> 生成时间：2026-05-19 01:18

包括 `CXXOperatorCallExpr`（用户类型重载）、`UnaryOperator`（基本类型 `!` `-` `~` `++` 等）、`BinaryOperator`（`+ - * / == ` 等）、`CompoundAssignOperator`（`+= -= *=` 等）。

| Operator | Count | First Seen In |
|---|---|---|
| CXXOperatorCallExpr::operator* | 21 | __hensel_step |
| CXXOperatorCallExpr::operator!= | 11 | __hensel_step |
| CXXOperatorCallExpr::operator++ | 11 | __hensel_step |
| CXXOperatorCallExpr::operator= | 9 | __hensel_step |
| CXXOperatorCallExpr::operator+ | 6 | __hensel_step |
| CXXOperatorCallExpr::operator*= | 4 | __hensel_step |
| CXXOperatorCallExpr::operator- | 3 | __hensel_step |
| CXXOperatorCallExpr::operator-> | 2 | __hensel_step |
| BinaryOperator::== | 1 | __upoly_const_term |