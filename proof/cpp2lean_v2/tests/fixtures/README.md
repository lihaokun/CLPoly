# Pass 1 parse 测试 fixtures

5 个最简单的 CLPoly 函数，覆盖 Pass 1 需要处理的基础 AST 构造。

## Fixture 列表

| Fixture | 源位置 | AST 覆盖 |
|---|---|---|
| `make_zp.cc` | `polynomial_factorize_zp.hh:20` | FunctionDecl、ReturnStmt、CXXConstructExpr、ParmVarDecl |
| `upoly_mod.cc` | `polynomial_factorize_zp.hh:39-46` | 多语句、VarDecl + LetStmt、CallExpr、MemberExpr |
| `upoly_divmod.cc` | `polynomial_factorize_zp.hh:49-56` | void 函数、non-const ref 参数（输出参数）|
| `symmetric_mod.cc` | `polynomial_factorize_univar.hh:94-103` | if-else、BinaryOperator、CompoundAssignOperator（`-=`） |
| `upoly_const_term.cc` | `polynomial_factorize_univar.hh:743-748` | if-return 链、CXXMemberCallExpr (`.empty()`, `.back()`, `.deg()`) |

## 测试策略

对每个 fixture：

1. Pass 1 parse 生成 `HIRFunc`
2. assert HIR₀ 不变量（允许 `UnknownStmt`/`UnknownExpr`，但记录数量）
3. 预期 HIR 结构详见每 fixture 的 `expected_hir.md`

## 使用

```python
from pathlib import Path
from tests.fixtures import load_fixture
from passes.pass1_parse import parse_pass

ast = load_fixture("make_zp")
hir0 = parse_pass(ast)
assert hir0.base_name == "__make_zp"
assert hir0.params[0].name == "val"
assert hir0.params[0].ty == BaseType.INT64
```
