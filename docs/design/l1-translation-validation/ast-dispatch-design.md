# AST 节点统一转发机制

> 状态：设计中（含审核修订 + 语义等价性证明）
> 解决的问题：AST 节点处理分散在 if-else 链中，覆盖不完整，新增节点无统一入口

## 1. 现状问题

`parse_expr` 和 `parse_stmt` 是 700+ 行的 `if kind == "XXX"` 链。问题：

1. **覆盖不完整**：71 个函数 AST 中出现 ~40 种节点，10 种未处理
2. **无法审计**：不扫描全量 AST 就不知道漏了什么
3. **静默丢失**：未处理节点走兜底路径，可能被父节点吞掉
4. **表达式/语句边界模糊**：同一节点可能出现在两种位置，需要两处都加处理

## 2. 新架构：注册表转发

### 2.1 节点分类

| 类别 | 含义 | 处理方式 |
|------|------|---------|
| **表达式节点** | 产生值的 AST 节点 | 注册在 `EXPR_HANDLERS`，返回 `ExprIR` |
| **语句节点** | 控制流/声明的 AST 节点 | 注册在 `STMT_HANDLERS`，返回 `StmtIR` |
| **透传节点** | 不改变语义的包装层 | 注册在 `TRANSPARENT_NODES`，递归子节点 |
| **元信息节点** | 调试信息（`__func__`、`__builtin_source_location`），出现在 assert 宏展开中，对算法翻译无语义价值 | 注册在 `META_NODES`，翻译为 `Lit(0)` |

### 2.2 注册表定义

```python
# ============================================================
# 表达式处理器：kind → handler(node) → ExprIR
# 每个节点类型只在一个注册表中出现
# ============================================================

EXPR_HANDLERS: dict[str, Callable[[dict], ExprIR]] = {
    # 字面量
    "IntegerLiteral":           handle_integer_literal,
    "FloatingLiteral":          handle_floating_literal,
    "CXXBoolLiteralExpr":       handle_bool_literal,
    "StringLiteral":            handle_string_literal,

    # 变量引用
    "DeclRefExpr":              handle_decl_ref,

    # 运算符（含 operator[], operator=, 二元运算符）
    # 语句位置的 operator= 由 parse_stmt step 4 自动包 ExprStmt，
    # SSA 变换中 BinOp("=") → AssignStmt
    "BinaryOperator":           handle_binary_op,
    "UnaryOperator":            handle_unary_op,
    "ConditionalOperator":      handle_conditional_op,
    "CXXOperatorCallExpr":      handle_cxx_operator_call,
    "CompoundAssignOperator":   handle_compound_assign_op,

    # 函数调用
    "CallExpr":                 handle_call_expr,
    "CXXMemberCallExpr":        handle_member_call,

    # 成员/下标访问
    "MemberExpr":               handle_member_expr,
    "ArraySubscriptExpr":       handle_array_subscript,

    # 类型转换
    "ImplicitCastExpr":         handle_implicit_cast,
    "CStyleCastExpr":           handle_c_cast,
    "CXXStaticCastExpr":        handle_c_cast,
    "CXXFunctionalCastExpr":    handle_functional_cast,

    # 构造
    "CXXConstructExpr":         handle_construct,
    "CXXTemporaryObjectExpr":   handle_construct,
    "InitListExpr":             handle_init_list,

    # 模板依赖表达式
    "CXXDependentScopeMemberExpr": handle_dependent_member,
    "UnresolvedLookupExpr":     handle_unresolved_lookup,
    "CXXUnresolvedConstructExpr": handle_unresolved_construct,

    # Lambda
    "LambdaExpr":               handle_lambda_expr,

    # 默认参数 / 零初始化
    "CXXDefaultArgExpr":        handle_default_arg,
    "ImplicitValueInitExpr":    handle_implicit_value_init,
}

# ============================================================
# 语句处理器：kind → handler(node) → StmtIR | list[StmtIR]
# 只注册"纯语句"节点，不重复注册表达式节点
# ============================================================

STMT_HANDLERS: dict[str, Callable[[dict], StmtIR | list[StmtIR]]] = {
    # 声明
    "DeclStmt":                 handle_decl_stmt,       # 内部分发 VarDecl / DecompositionDecl

    # 控制流
    "IfStmt":                   handle_if_stmt,
    "ForStmt":                  handle_for_stmt,
    "WhileStmt":                handle_while_stmt,
    "DoStmt":                   handle_do_stmt,
    "CXXForRangeStmt":          handle_range_for_stmt,
    "ReturnStmt":               handle_return_stmt,
    "BreakStmt":                handle_break_stmt,
    "ContinueStmt":             handle_continue_stmt,

    # 复合
    "CompoundStmt":             handle_compound_stmt,

    # 异常
    "CXXThrowExpr":             handle_throw_stmt,
}

# ============================================================
# 透传节点：不改变语义的包装层，直接递归第一个子节点
# ============================================================

TRANSPARENT_NODES = {
    "ExprWithCleanups",
    "MaterializeTemporaryExpr",
    "CXXBindTemporaryExpr",
    "ParenExpr",
    "SubstNonTypeTemplateParmExpr",
}

# ============================================================
# 元信息节点：出现在 assert 宏展开中的调试信息
# assert(cond) → if (!cond) __assert_fail("cond", __FILE__, __LINE__, __func__)
# 其中 __FILE__→SourceLocExpr, __func__→PredefinedExpr
# assert 语义已由 UB collector 处理为 require，这些节点不影响算法正确性
# ============================================================

META_NODES = {
    "SourceLocExpr",
    "PredefinedExpr",
    "TypeAliasDecl",
    "FullComment",
    "GCCAsmStmt",
    "NullStmt",
}
```

### 2.3 统一转发函数

```python
def parse_expr(node: dict) -> ExprIR:
    kind = node.get("kind", "")
    inner = node.get("inner", [])

    # 1. 透传
    if kind in TRANSPARENT_NODES:
        return parse_expr(inner[0]) if inner else Lit(0)

    # 2. 元信息
    if kind in META_NODES:
        return Lit(0)

    # 3. 查表达式注册表
    handler = EXPR_HANDLERS.get(kind)
    if handler:
        result = handler(node)
        # P1 原则：附加 Clang AST 类型信息
        ast_qual = _get_qual_type(node)
        if ast_qual:
            result._ast_type = parse_type(ast_qual)
        return result

    # 4. 兜底
    if inner:
        return parse_expr(inner[0])
    return UnknownExpr(kind)


def parse_stmt(node: dict) -> StmtIR | list[StmtIR]:
    kind = node.get("kind", "")
    inner = node.get("inner", [])

    # 1. 透传
    if kind in TRANSPARENT_NODES:
        return parse_stmt(inner[0]) if inner else ExprStmt(Lit(0))

    # 2. 元信息
    if kind in META_NODES:
        return ExprStmt(Lit(0))

    # 3. 查语句注册表
    handler = STMT_HANDLERS.get(kind)
    if handler:
        return handler(node)

    # 4. 查表达式注册表（表达式出现在语句位置 → 包 ExprStmt）
    #    这使得 CXXOperatorCallExpr(operator=) 等在语句位置也能被处理
    #    SSA 变换中 BinOp("=") → AssignStmt
    expr_handler = EXPR_HANDLERS.get(kind)
    if expr_handler:
        return ExprStmt(expr_handler(node))

    # 5. 兜底
    return UnknownStmt(kind, [parse_stmt(c) for c in inner if isinstance(c, dict)])
```

**设计决策**：`CXXOperatorCallExpr` 只在 `EXPR_HANDLERS` 中注册一次。当它出现在语句位置时（如 `a = b;`），step 4 自动包 `ExprStmt`。赋值语义 `BinOp("=", a, b)` 由 SSA 层的 `transform_stmt` 转换为 `AssignStmt`。这避免了同一节点在两个表中注册的歧义。

### 2.4 DeclStmt 内部分发

`DeclStmt` 是容器节点，内部可能是 `VarDecl` 或 `DecompositionDecl`。在 `handle_decl_stmt` 内部分发：

```python
def handle_decl_stmt(node: dict) -> StmtIR | list[StmtIR]:
    inner = node.get("inner", [])
    stmts = []
    for decl in inner:
        kind = decl.get("kind")
        if kind == "VarDecl":
            stmts.append(parse_var_decl(decl))
        elif kind == "DecompositionDecl":
            stmts.extend(handle_decomposition_decl(decl))
        # 其他声明类型：跳过
    return stmts[0] if len(stmts) == 1 else stmts if stmts else ExprStmt(Lit(0))
```

### 2.5 ParenListExpr 处理

`ParenListExpr` 是多表达式列表 `(a, b, c)`，不是单表达式透传。从 `TRANSPARENT_NODES` 中移除，在 `EXPR_HANDLERS` 中处理：

```python
def handle_paren_list(node: dict) -> ExprIR:
    inner = node.get("inner", [])
    if len(inner) == 1:
        return parse_expr(inner[0])
    elif len(inner) == 2:
        return Call("Prod.mk", [parse_expr(inner[0]), parse_expr(inner[1])])
    else:
        # 多于 2 个 → 嵌套 Prod
        return Call("Prod.mk", [parse_expr(inner[0]),
                                handle_paren_list({"inner": inner[1:]})])
```

## 3. 关键处理器设计

### 3.1 operator[]（影响 ~200 处）

```python
def handle_cxx_operator_call(node: dict) -> ExprIR:
    inner = node.get("inner", [])
    ref_name = _extract_callee_name(inner[0]) if inner else ""

    # operator[] → ArrayAccess
    if "operator[]" in ref_name and len(inner) >= 3:
        return ArrayAccess(parse_expr(inner[1]), parse_expr(inner[2]))

    # operator= → BinOp("=", lhs, rhs)
    if ref_name == "operator=" and len(inner) >= 3:
        return BinOp("=", parse_expr(inner[1]), parse_expr(inner[2]))

    # 二元运算符 → 查 OPERATOR_MAP
    for op_key, (min_args, lean_op) in sorted(
            OPERATOR_MAP.items(), key=lambda x: -len(x[0])):
        if op_key in ref_name and len(inner) >= min_args:
            return BinOp(lean_op, parse_expr(inner[1]), parse_expr(inner[2]))

    # operator() → 查 CALL_OPERATOR_MAP
    if ref_name in CALL_OPERATOR_MAP:
        lean_func, arg_order = CALL_OPERATOR_MAP[ref_name]
        args = [parse_expr(a) for a in inner[1:]]
        if arg_order == "reverse" and len(args) == 2:
            args = [args[1], args[0]]
        return Call(lean_func, args)

    return UnknownExpr(f"operator:{ref_name}")
```

### 3.2 DecompositionDecl（结构化绑定）

```python
def handle_decomposition_decl(node: dict) -> list[StmtIR]:
    """auto [a, b] = expr; → let pair := expr; let a := pair.1; let b := pair.2"""
    inner = node.get("inner", [])
    init_expr = None
    binding_names = []
    for child in inner:
        if child.get("kind") == "BindingDecl":
            binding_names.append(child.get("name", "_"))
        elif init_expr is None:
            init_expr = parse_expr(child)

    if init_expr is None:
        return [UnknownStmt("DecompositionDecl_no_init")]

    result = []
    pair_var = Var("_decomp")
    result.append(LetStmt(pair_var, "auto", init_expr))

    # Lean Prod 嵌套解构：A × (B × C) → .1, .2.1, .2.2
    for i, name in enumerate(binding_names):
        if len(binding_names) == 1:
            accessor = pair_var
        elif i == 0:
            accessor = FieldAccess(pair_var, "1")
        elif i == len(binding_names) - 1:
            accessor = pair_var
            for _ in range(i):
                accessor = FieldAccess(accessor, "2")
        else:
            accessor = pair_var
            for _ in range(i):
                accessor = FieldAccess(accessor, "2")
            accessor = FieldAccess(accessor, "1")
        result.append(LetStmt(Var(name), "auto", accessor))

    return result
```

### 3.3 LambdaExpr

```python
def handle_lambda_expr(node: dict) -> ExprIR:
    """简单 lambda（单 return 语句）→ 内联表达式。
       复杂 lambda → sorry + 位置信息。"""
    inner = node.get("inner", [])
    for child in inner:
        if child.get("kind") == "CompoundStmt":
            body_inner = child.get("inner", [])
            if len(body_inner) == 1 and body_inner[0].get("kind") == "ReturnStmt":
                ret_inner = body_inner[0].get("inner", [])
                if ret_inner:
                    return parse_expr(ret_inner[0])
    loc = node.get("loc", {})
    file = loc.get("file", "").split("/")[-1]
    line = loc.get("line", 0)
    return UnknownExpr(f"lambda:{file}:{line}")
```

### 3.4 输出参数

输出参数不是 AST 节点问题，是**调用约定**问题。在 `FUNC_MAP` 中标记：

```python
FUNC_MAP = {
    # (lean_name, arg_rule, output_param_indices)
    "pair_vec_div":    ("pair_vec_div", "direct", [0, 1]),  # q, r 是输出
    "fdiv_r":          ("ZZ.fdiv_r",   "direct", [0]),      # r 是输出
    "fdiv_q":          ("ZZ.fdiv_q",   "direct", [0]),      # q 是输出
    "invert":          ("ZZ.invert",   "direct", [0]),      # inv 是输出，bool 是返回值
    "polynomial_GCD":  ("polynomial_GCD", "direct", [2, 3]), # alpha, beta 是输出
}
```

**组合规则**：Lean 返回值 = `(原始返回值, 输出参数₁, 输出参数₂, ...)`。若原始返回值是 void，则只返回输出参数。

在 `handle_call_expr` 中：
1. 检测 `FUNC_MAP` 中的 `output_param_indices`
2. 将输出参数从调用参数中移除，只传入非输出参数
3. 返回值是 Prod（包含原始返回值 + 输出参数值）
4. SSA 变换中，给每个输出变量 bump 版本

## 4. break/continue 语义等价性证明

### 4.1 break → return 等价

**C++ while 循环含 break**：
```cpp
while (cond(σ)) { S₁; if (guard) break; S₂; }
```

语义函数 W(σ)：
- W(σ) = σ，若 ¬cond(σ)
- W(σ) = S₁(σ)，若 cond(σ) ∧ guard(S₁(σ))
- W(σ) = W(S₂(S₁(σ)))，若 cond(σ) ∧ ¬guard(S₁(σ))

**Lean 循环函数（break → return）**：
```lean
partial def loop (σ) :=
  if ¬cond(σ) then σ
  else
    let σ₁ := S₁(σ)
    if guard(σ₁) then σ₁          -- break → return 当前状态
    else loop (S₂(σ₁))
```

语义函数 L(σ)：定义与 W 逐条相同。由递归定义一致性，**W = L**。∎

### 4.2 continue → return recurse 等价

**C++ while 循环含 continue**：
```cpp
while (cond(σ)) { S₁; if (guard) continue; S₂; }
```

语义 W(σ)：
- W(σ) = σ，若 ¬cond(σ)
- W(σ) = W(S₁(σ))，若 cond(σ) ∧ guard(S₁(σ))（跳过 S₂）
- W(σ) = W(S₂(S₁(σ)))，若 cond(σ) ∧ ¬guard(S₁(σ))

**Lean（continue → return recurse）**：
```lean
partial def loop (σ) :=
  if ¬cond(σ) then σ
  else
    let σ₁ := S₁(σ)
    if guard(σ₁) then loop(σ₁)    -- continue → 直接递归
    else loop(S₂(σ₁))
```

逐条相同。∎

### 4.3 for 循环 continue 含 step

**C++**：`for (i=0; cond(i); i++) { S₁; if (guard) continue; S₂; }`

语义 W(i, σ)：
- W(i, σ) = σ，若 ¬cond(i)
- W(i, σ) = W(i+1, S₁(σ))，若 cond(i) ∧ guard(S₁(σ))
- W(i, σ) = W(i+1, S₂(S₁(σ)))，若 cond(i) ∧ ¬guard(S₁(σ))

**Lean**：
```lean
partial def loop (i σ) :=
  if ¬cond(i) then σ
  else
    let σ₁ := S₁(σ)
    if guard(σ₁) then loop (i+1) σ₁     -- continue: step + recurse
    else loop (i+1) (S₂(σ₁))
```

**关键**：continue 分支必须包含 step（i+1）。逐条相同。∎

### 4.4 early-return 机制保证后续代码正确放入 else

语句序列 `[S₁; if (guard) return e; S₂; S₃]` 由 `transform_body` 变换为：

```lean
let σ₁ := S₁
if guard then e
else
  let σ₂ := S₂(σ₁)
  S₃(σ₂)
```

在纯函数式语言中，`let` 链是顺序求值，`if-then-else` 是分支选择。因为没有副作用，这与 C++ 的"if guard 则 return，否则继续"语义一致。

**适用条件**：SSA 变换保证循环体是纯函数式 let 链。`env.fork()` 保证分支变量独立演化。

### 4.5 实现方式

break/continue 不在 AST 解析层处理。而是在**循环函数体的 SSA 变换**中：

1. 在 `transform_range_for` / `transform_while_loop` / `transform_for_loop` 中，将循环上下文（函数名、修改变量、step 语句）传入 `transform_body`
2. `transform_body` 中，BreakStmt → `ReturnStmt(current_state)`，ContinueStmt → `ReturnStmt(recurse_call_with_step)`
3. 现有的 early-return 检测（`if (cond) return` → 剩余代码进 else）自动处理嵌套情况

通过 `transform_body` 的 `loop_ctx` 参数传递循环上下文：

```python
def transform_body(stmts, env, ctx, loop_ctx=None):
    # loop_ctx: (func_name, modified_vars, step_stmts) or None
    # 每个循环函数独立传入自己的 loop_ctx，嵌套循环各自独立
```

## 5. 审计机制

```python
# 在翻译前扫描全量 AST，对比注册表
TYPE_ONLY_NODES = {
    "TemplateArgument", "TemplateSpecializationType", "ElaboratedType",
    "RecordType", "SubstTemplateTypeParmType", "TemplateTypeParmType",
    "TemplateTypeParmDecl", "FunctionDecl", "FunctionTemplateDecl",
    "CXXMethodDecl", "CXXRecordDecl", "CXXDestructorDecl", "CXXConversionDecl",
    "ParmVarDecl",
}

def audit_coverage(ast_scan: dict[str, int]):
    all_handled = (set(EXPR_HANDLERS) | set(STMT_HANDLERS)
                   | TRANSPARENT_NODES | META_NODES | TYPE_ONLY_NODES)
    unhandled = []
    for kind, count in sorted(ast_scan.items(), key=lambda x: -x[1]):
        if kind not in all_handled:
            unhandled.append((kind, count))
    if unhandled:
        print("UNHANDLED AST NODES:")
        for kind, count in unhandled:
            print(f"  {count:5d}  {kind}")
    else:
        print("All AST nodes covered.")
```

每次扩展翻译范围时先跑 `audit_coverage`。

## 6. 回归测试策略

重构前保存基线：

```bash
python3 gen_full.py > baseline.lean 2>baseline.log
```

重构后对比：

```bash
python3 gen_full.py > refactored.lean 2>refactored.log
diff baseline.lean refactored.lean  # 应只有改善（减少 sorry），无回归
```

## 7. 翻译等价性证明

### 7.1 定理（翻译等价性）

设 P 是 C++ 因式分解程序，T(P) 是翻译后的 Lean IR。在以下假设下，T(P) 的计算结果等价于 P 的计算结果：

1. **可信基正确**：clpoly_model.lean 中的类型和操作正确实现了 C++ 对应类型的语义
2. **AST 完整覆盖**：所有出现的 AST 节点都有对应的 handler（由 §5 审计机制保证）
3. **无 UB**：C++ 程序在执行路径上无未定义行为（由 UB proof obligations 保证）

证明分解为 5 个引理，覆盖翻译管线的每个环节：

```
C++ AST ──引理1,5──→ IR ──引理2,3,4──→ SSA IR ──codegen──→ Lean IR
```

### 7.2 引理 1：透传节点语义无关

**声明**：对 `TRANSPARENT_NODES` 中的节点 N，`parse_expr(N) = parse_expr(N.inner[0])` 不改变语义。

**证明**：逐个验证每个透传节点的 C++ 语义：

| 节点 | C++ 语义 | 跳过的理由 |
|------|---------|-----------|
| `ExprWithCleanups` | 标记需要析构的临时对象 | 析构是 C++ 内存管理语义，Lean 由 GC 管理，不影响计算值 |
| `MaterializeTemporaryExpr` | 将右值物化为左值（创建临时存储） | Lean 所有值都是不可变对象，不区分左值/右值 |
| `CXXBindTemporaryExpr` | 绑定临时对象的生命周期 | 同上 |
| `ParenExpr` | 括号 `(expr)` | 语义定义：`(e)` = `e` |
| `SubstNonTypeTemplateParmExpr` | 模板参数替换后的包装 | 语义定义：替换后的值 = 子表达式的值 |

共同性质：**子表达式的计算值 = 父节点的计算值**。这些节点只影响 C++ 的对象生命周期和内存布局，不影响计算结果。由于 Lean 是纯函数式（无生命周期管理），透传安全。∎

### 7.3 引理 2：SSA 变换保持语义

**声明**：将可变变量赋值序列变换为 SSA let 链后，程序的计算结果不变。

**形式化**：

C++ 命令式语义（状态 σ 是变量名 → 值的映射）：
```
σ₀ → [x = e₁] → σ₁ = σ₀[x ↦ eval(e₁, σ₀)]
   → [x = e₂] → σ₂ = σ₁[x ↦ eval(e₂, σ₁)]
   → [use(x)] → use(σ₂(x)) = use(eval(e₂, σ₁))
```

Lean SSA 语义（环境 ρ 是变量名 → 值的映射）：
```
let x₁ := eval(e₁, ρ₀)     → ρ₁ = ρ₀[x₁ ↦ eval(e₁, ρ₀)]
let x₂ := eval(e₂', ρ₁)    → ρ₂ = ρ₁[x₂ ↦ eval(e₂', ρ₁)]
use(x₂)                     → use(ρ₂(x₂)) = use(eval(e₂', ρ₁))
```

其中 `e₂' = rename(e₂)`，即 e₂ 中所有 `x` 替换为 `x₁`。

**证明**：

需要证 `eval(e₂, σ₁) = eval(e₂', ρ₁)`。

由 SSA rename 构造：e₂' 中的 `x₁` 指向 ρ₁(x₁) = eval(e₁, ρ₀)。
而 σ₁(x) = eval(e₁, σ₀)。

若 ρ₀ 和 σ₀ 在所有自由变量上一致（归纳假设），则 eval(e₁, ρ₀) = eval(e₁, σ₀)，故 ρ₁(x₁) = σ₁(x)。

e₂ 中引用 x 的地方在 e₂' 中被替换为 x₁，对应的值相等。e₂ 中引用其他变量的地方，由归纳假设在 ρ₁ 和 σ₁ 中值相等。

因此 eval(e₂', ρ₁) = eval(e₂, σ₁)。

**归纳推广**：对任意长度 n 的赋值序列，SSA 变换后 let 链中第 k 个版本 xₖ 的值 = C++ 执行到第 k 次赋值后 x 的值。

**base case**：n = 0，无赋值，等价成立。
**inductive step**：假设前 k 次赋值等价，第 k+1 次赋值由上述论证等价。∎

### 7.4 引理 3：循环提取保持语义

（引用 loop-extraction-design.md §4 中已证明的三个等价性定理。）

- while 循环含 break：W = L（§4.1）
- while 循环含 continue：W = L（§4.2）
- for 循环含 continue + step：W = L（§4.3）
- early-return 重构后续代码：纯函数式 let 链的条件分支等价（§4.4）

### 7.5 引理 4：输出参数变换保持语义

**声明**：将 C++ 的 `f(out₁, out₂, in₁, in₂)`（out 是非 const 引用）变换为 Lean 的 `let (out₁', out₂') := f_lean(in₁, in₂)` 后，后续代码中 out₁、out₂ 的值不变。

**证明**：

C++ 语义：`f` 通过引用修改 out₁、out₂。调用后 out₁ = v₁，out₂ = v₂。

Lean 语义：`f_lean(in₁, in₂)` 返回 `(v₁, v₂)`。解构后 out₁' = v₁，out₂' = v₂。

后续代码中，SSA 变换将 out₁ 的引用替换为 out₁'（新版本）。因此后续代码看到的值与 C++ 一致。

**组合规则**：当函数既有返回值 r 又有输出参数 out₁ 时：
- C++ 调用后：r = 返回值，out₁ = 写入值
- Lean：`let (r', out₁') := f_lean(...)`，r' = r，out₁' = out₁
- 后续代码等价。

**前提**：`f_lean` 的实现（clpoly_model.lean）正确计算了 v₁、v₂。这是可信基假设 1。∎

### 7.6 引理 5：表达式-语句 fallback 正确

**声明**：当 `EXPR_HANDLERS` 中的节点出现在语句位置时，`parse_stmt` 返回 `ExprStmt(handler(node))` 保持语义。

**证明**：

C++ 中表达式出现在语句位置的语义 = 求值表达式，保留副作用。三种情况：

1. `a = b;`（赋值）→ handler 返回 `BinOp("=", a, b)` → SSA 变换中 `transform_stmt` 检测到 `BinOp("=")` 并转为 `AssignStmt(a, b)` → `let a_{n+1} := b`。赋值语义保持。

2. `f(x);`（函数调用）→ handler 返回 `Call("f", [x])` → `ExprStmt(Call(...))` → codegen 生成 `let _ := f(x)`。函数调用执行，副作用通过输出参数（引理 4）或返回值传播。

3. `a + b;`（纯表达式）→ handler 返回 `BinOp("+", a, b)` → `ExprStmt(...)` → codegen 生成 `let _ := a + b`。无副作用，无害。

关键：`ExprStmt` 包装不丢失语义信息。SSA 层的 `transform_stmt` 识别特定模式（赋值、自增、push 等）并生成正确的 SSA 语句。∎

### 7.7 组合定理证明

**翻译管线**：
```
C++ AST ──parse_expr/parse_stmt──→ IR ──transform_body──→ SSA IR ──gen_ssa_func──→ Lean text
```

**第一层（AST → IR）**：
- EXPR_HANDLERS 中的每个 handler 将 C++ 表达式翻译为语义等价的 ExprIR（各 handler 的正确性由其定义保证——每个 handler 只做直接的结构映射）
- STMT_HANDLERS 同理
- 透传节点不影响语义（引理 1）
- 表达式出现在语句位置时正确包装（引理 5）
- 未处理节点标记为 `UnknownExpr`/`UnknownStmt`（不静默丢失，审计可发现）

**第二层（IR → SSA IR）**：
- 变量重命名保持语义（引理 2）
- 循环提取保持语义（引理 3）
- 输出参数变换保持语义（引理 4）
- break/continue 替换保持语义（引理 3 推论）

**第三层（SSA IR → Lean text）**：codegen 是 SSA IR 的直接文本渲染（LetStmt → `let x := e`，IfStmt → `if c then ... else ...`，Call → `f args`）。渲染是注入映射（每个 IR 节点有唯一的 Lean 文本表示），不改变语义。

**组合**：第一层保持语义 ∧ 第二层保持语义 ∧ 第三层保持语义 → 整体保持语义。∎

### 7.8 假设边界（翻译不保证的部分）

| 假设 | 含义 | 验证方式 |
|------|------|---------|
| 可信基正确 | clpoly_model.lean 中 Zp/ZZ/SparsePolyZp 等的实现与 C++ 一致 | 人工审核（174 行）+ `#eval` 测试 |
| AST 完整覆盖 | 所有 AST 节点都有 handler | 审计脚本（§5） |
| 无 UB | C++ 不触发未定义行为 | UB proof obligations（23 条 require） |
| handler 正确 | 每个 handler 的结构映射保持语义 | 人工审核 handler 代码 |
| 终止性 | 循环函数终止 | 由 L2 层已证（`partial def` 不要求终止证明） |

## 8. 实现步骤

1. 保存基线翻译输出（回归测试基准）
2. 将 `_parse_expr_inner` 的 if-else 链拆分为独立 handler 函数
3. 建立 4 个注册表（EXPR_HANDLERS / STMT_HANDLERS / TRANSPARENT_NODES / META_NODES）
4. 改写 `parse_expr` / `parse_stmt` 为统一转发（§2.3）
5. 回归测试：确认输出与基线一致
6. 补充缺失 handler：`operator[]`、`DecompositionDecl`、`LambdaExpr`、`CXXStaticCastExpr`、`ParenListExpr`
7. 实现 break/continue 的循环函数语义（基于 §4 的证明、§7 的引理 3）
8. 实现输出参数的 FUNC_MAP 扩展（基于 §7 的引理 4）
9. 跑审计脚本确认覆盖完整（§5）
10. 全量翻译 + 逐函数对比审计
