# 循环提取设计：loop-as-function 架构

> 状态：设计中
> 解决的问题：循环后代码截断、嵌套循环作用域、`let rec` 内联导致的 codegen 复杂性

## 1. 问题

当前架构把循环内联为 `let rec`：

```lean
partial def foo_ir (...) :=
  let x := ...
  let rec while_123 (x) :=
    if done then x
    else while_123 (x + 1)
  while_123 x              -- 整个 let rec 块是一个表达式
  -- ↑ 这里是最终返回值，后续代码全部丢失
```

导致 6 个系统性 bug：
1. **循环后代码截断**（~30 函数）：`transform_body` 遇到 TailRec 就 break
2. **嵌套循环作用域**：内层 `let rec` 的初始调用截断外层 body
3. **输出参数丢失**：函数调用结果被 `let _ :=` 丢弃
4. **返回值错误**：循环返回单个累加器，丢失其他修改变量
5. **复杂 codegen**：`gen_tailrec` 需要处理 bind_var、latest 追踪等

## 2. 新架构：loop-as-function

**核心思想**：每个循环提取为独立的 `partial def`，主函数中只有调用。

### 2.1 翻译模型

```
C++:                              Lean:
┌─────────────────────┐           ┌─────────────────────────────┐
│ void foo() {        │           │ partial def while_123       │
│   int x = 0;        │           │     (x : Int) : Int :=     │
│   while (x < 10) {  │    ──→    │   if x >= 10 then x        │
│     x += 1;         │           │   else while_123 (x + 1)   │
│   }                 │           │                             │
│   return x + 1;     │           │ partial def foo_ir : Int := │
│ }                   │           │   let x := 0               │
└─────────────────────┘           │   let x_2 := while_123 x   │
                                  │   (x_2 + 1)                │
                                  └─────────────────────────────┘
```

主函数体是纯 `let` 链——没有 `let rec`，没有作用域问题。

### 2.2 循环函数的签名

循环函数的**参数** = 循环体内读写的外部变量（即循环携带的状态）。
循环函数的**返回值** = 循环退出时这些变量的最终值。

单变量：
```lean
partial def loop_123 (result : Array Zp) : Array Zp :=
  ...
```

多变量：
```lean
partial def loop_456 (x : Int) (y : Int) (found : Bool) : Int × Int × Bool :=
  if done then (x, y, found)
  else loop_456 x' y' found'
```

主函数中解构：
```lean
let (x_2, y_2, found_2) := loop_456 x_1 y_1 found_1
-- 后续代码用 x_2, y_2, found_2
```

### 2.3 三种循环的翻译

**while 循环**：
```cpp
while (cond) { body; }
```
```lean
partial def while_N (vars...) : RetType :=
  if ¬cond then vars        -- 退出
  else
    <body>
    while_N vars'            -- 递归
```

**for 循环**：
```cpp
for (init; cond; step) { body; }
```
```lean
partial def for_N (i : Int) (vars...) : RetType :=
  if ¬cond then vars
  else
    <body>
    for_N (i + step) vars'
```

**range-for 循环**：
```cpp
for (auto& elem : coll) { body; }
```
```lean
partial def range_N (idx : UInt64) (coll : Array T) (vars...) : RetType :=
  if idx >= coll.size then vars
  else
    let elem := coll[idx]!
    <body>
    range_N (idx + 1) coll vars'
```

### 2.4 嵌套循环

自然处理——每层循环各自一个 `partial def`，互相调用：

```cpp
for (auto& a : arr1) {
    for (auto& b : arr2) {
        result.push_back(a + b);
    }
}
```

```lean
partial def range_inner (idx : UInt64) (coll : Array T)
    (result : Array U) : Array U :=
  if idx >= coll.size then result
  else
    let b := coll[idx]!
    let result_2 := result.push (a + b)    -- a 来自闭包/参数
    range_inner (idx + 1) coll result_2

partial def range_outer (idx : UInt64) (coll : Array T)
    (result : Array U) (arr2 : Array T) : Array U :=
  if idx >= coll.size then result
  else
    let a := coll[idx]!
    let result_2 := range_inner 0 arr2 result  -- 内层循环 = 函数调用
    range_outer (idx + 1) coll result_2 arr2

partial def foo_ir (...) :=
  let result := #[]
  let result_2 := range_outer 0 arr1 result arr2
  result_2
```

**嵌套循环的作用域问题消失了**——因为没有 `let rec` 嵌套。

### 2.5 循环后代码

**截断问题消失了**——循环调用就是一个普通的 `let`：

```cpp
auto sqf = __squarefree_Zp(f);    // while loop inside
// 循环后代码
if (!c.empty()) {
    auto g = __extract_pth_root(c);
    // ...
}
return result;
```

```lean
-- while 循环提取为 while_N
let (result_2, c_2) := while_N result_1 c_1 w_1
-- 循环后代码自然保留
if ¬(Array.isEmpty c_2) then
  let g := __extract_pth_root_ir c_2
  ...
else result_2
```

## 3. 实现方案

### 3.1 数据结构变更

在 `SSAFunc` 中加 `aux_defs` 字段，存放提取出的循环函数：

```python
@dataclass
class SSAFunc:
    name: str
    params: list[ParamIR]
    ret_type: TypeIR
    body: list[StmtIR]
    aux_defs: list[SSAFunc] = field(default_factory=list)  # 新增
```

### 3.2 SSA 变换（ssa_transform.py）

**`transform_range_for` / `transform_for_loop` / `transform_while_loop`**：

不再返回 `[coll_init, idx_init, TailRec]`，改为：

1. 构造循环函数 `SSAFunc`（参数 = 循环变量，body = 循环体 + 递归调用）
2. 注册到当前函数的 `aux_defs`
3. 返回一个 `LetStmt` 或多个 `LetStmt`（调用循环函数 + 解构返回值）

```python
def transform_range_for(stmt, env, func_ctx):
    # ... 识别循环变量 ...
    
    # 1. 构造循环函数
    loop_func = SSAFunc(
        name=f"range_{loop_id}",
        params=[idx_param, coll_param] + modified_var_params,
        ret_type=return_type,  # 单变量或 Prod
        body=loop_body + [recursive_call],
    )
    func_ctx.aux_defs.append(loop_func)
    
    # 2. 返回调用语句
    call_expr = Call(f"range_{loop_id}", [Lit(0), coll_var] + current_vars)
    if len(modified_vars) == 1:
        new_var = env.bump(modified_vars[0])
        return [LetStmt(new_var, var_type, call_expr)]
    else:
        # 多变量：解构 Prod
        result_stmts = []
        pair_var = env.bump(f"_loop_result_{loop_id}")
        result_stmts.append(LetStmt(pair_var, "auto", call_expr))
        for i, name in enumerate(modified_vars):
            new_var = env.bump(name)
            result_stmts.append(LetStmt(new_var, env.get_type(name),
                                        FieldAccess(pair_var, str(i+1))))
        return result_stmts
```

**`transform_body`**：删除 TailRec break 逻辑——因为不再生成 TailRec，循环返回的是普通 LetStmt。

### 3.3 代码生成（lean_codegen.py）

**`gen_ssa_func`**：先输出 `aux_defs`（循环函数），再输出主函数：

```python
def gen_ssa_func(func: SSAFunc) -> str:
    parts = []
    # 先输出所有辅助循环函数
    for aux in func.aux_defs:
        parts.append(gen_ssa_func(aux))  # 递归：嵌套循环的循环函数
    # 再输出主函数
    parts.append(gen_main_func(func))
    return "\n\n".join(parts)
```

**删除 `gen_tailrec`**——不再需要。循环函数就是普通的 `partial def`。

### 3.4 多变量返回值

当循环修改多个变量时，返回 `Prod`：

```lean
partial def while_123 (x : Int) (y : Int) : Int × Int :=
  if done then (x, y)
  else while_123 (x + 1) (y - 1)

-- 调用方解构
let result := while_123 x_1 y_1
let x_2 := result.1
let y_2 := result.2
```

单变量时直接返回，不包 Prod。

### 3.5 闭包变量

内层循环可能引用外层变量（如嵌套 range-for 中 `a` 来自外层）。两种处理方式：

**方案 A**：将闭包变量作为循环函数的额外参数（纯函数，无闭包）：
```lean
partial def range_inner (idx : UInt64) (coll : Array T)
    (result : Array U) (a : Zp) : Array U :=  -- a 是闭包变量
```

**方案 B**：利用 Lean 的 `where` 或闭包（更简洁但可能影响 partial def）。

**选方案 A**——纯函数，无闭包依赖，更清晰。

## 4. 自然解决的问题

| 问题 | 旧架构 | 新架构 |
|------|--------|--------|
| 循环后代码截断 | `let rec` 后 break | 循环 = 普通 `let`，后续代码自然保留 |
| 嵌套循环作用域 | 内层 `let rec` 截断外层 body | 每层各自 `partial def`，互相调用 |
| 返回值不完整 | `_pick_accumulator` 只选一个变量 | 返回所有修改变量的 Prod |
| codegen 复杂度 | `gen_tailrec` 有 bind_var/latest/... | 循环函数 = 普通 `partial def`，无特殊处理 |
| 循环变量名冲突 | 嵌套时 `__coll`/`__idx` 同名 shadow | 每层在独立函数中，无冲突 |

## 5. 不变的部分

- Clang AST 解析（`clang_ast.py`）：不变
- 类型映射（`class_map.py`）：不变
- UB 收集（`ub_collector.py`）：不变
- `identify_loop_vars`：不变（仍用于确定循环参数）
- `_collect_local_decls`：不变

## 6. 实现步骤

1. 在 `SSAFunc` 加 `aux_defs` 字段
2. 重写 `transform_range_for`：生成 SSAFunc + LetStmt 调用
3. 重写 `transform_for_loop` / `transform_while_loop`：同上
4. 修改 `transform_body`：删除 TailRec break
5. 修改 `gen_ssa_func`：先输出 aux_defs
6. 删除 `gen_tailrec`、`_pick_accumulator`、`latest` 计算
7. 测试：71 函数全量翻译 + 逐函数审计

## 7. 闭包变量识别

循环体引用但未声明的变量 = 闭包变量，必须作为循环函数的额外参数。

### 7.1 算法

```python
def collect_free_vars(stmts: list[StmtIR]) -> set[Var]:
    """收集语句列表中的自由变量（被引用但未在本作用域声明）。"""
    referenced = set()   # (name, version) 被引用
    declared = set()     # (name, version) 在 body 内声明

    def scan_expr(expr):
        if isinstance(expr, Var):
            referenced.add((expr.name, expr.version))
        elif isinstance(expr, BinOp):
            scan_expr(expr.lhs); scan_expr(expr.rhs)
        elif isinstance(expr, Call):
            for a in expr.args: scan_expr(a)
        elif isinstance(expr, ArrayAccess):
            scan_expr(expr.arr); scan_expr(expr.idx)
        elif isinstance(expr, FieldAccess):
            scan_expr(expr.obj)
        elif isinstance(expr, ArrayPush):
            scan_expr(expr.arr); scan_expr(expr.elem)
        elif isinstance(expr, Cast):
            scan_expr(expr.expr)
        elif isinstance(expr, CondExpr):
            scan_expr(expr.cond); scan_expr(expr.then_e); scan_expr(expr.else_e)
        elif isinstance(expr, UnaryOp):
            scan_expr(expr.operand)

    def scan_stmt(stmt):
        if isinstance(stmt, LetStmt):
            declared.add((stmt.var.name, stmt.var.version))
            scan_expr(stmt.value)
        elif isinstance(stmt, IfStmt):
            scan_expr(stmt.cond)
            for s in stmt.then_body + stmt.else_body: scan_stmt(s)
        elif isinstance(stmt, ReturnStmt) and stmt.value:
            scan_expr(stmt.value)
        elif isinstance(stmt, ExprStmt):
            scan_expr(stmt.expr)

    for s in stmts: scan_stmt(s)
    return {Var(n, v) for n, v in referenced - declared}
```

### 7.2 使用

在 `transform_range_for` 等中，SSA 变换循环体后调用：

```python
body_stmts = transform_body(all_body_raw, loop_env)
free = collect_free_vars(body_stmts)
# free 中的变量 = 闭包变量，加入循环函数参数
```

闭包变量作为循环函数的**只读参数**（不出现在返回值中）。

## 8. break/continue 处理

**break**：在循环函数体中，`if (cond) break` → `if cond then (modified_vars...) else <continue body>`。即直接返回当前状态。

**continue**：`if (cond) continue` → `if cond then loop_func(idx+1, vars...) else <continue body>`。即跳过本轮 body，直接递归。

这与当前 TailRec 中的处理逻辑一致，只是发生在独立函数中而非 `let rec` 内。

## 9. Prod 解构

Lean 的 `A × B × C` = `A × (B × C)`（右结合）。解构：

| 变量数 | 返回类型 | 解构方式 |
|--------|---------|---------|
| 1 | `A` | `let a := loop(...)` |
| 2 | `A × B` | `let r := loop(...); let a := r.1; let b := r.2` |
| 3 | `A × B × C` | `let r := loop(...); let a := r.1; let b := r.2.1; let c := r.2.2` |
| n | 嵌套 Prod | `.1` 取第一个，`.2....2` 取最后一个 |

## 10. 不在本次范围

- **输出参数**（根因 B）：`pair_vec_div(q, r, f, g)` 的 q/r 是输出参数。需要单独设计，不在本次重构中。
- **operator[] / operator- 未翻译**（根因 C）：AST 解析层的问题，与循环架构无关。
- **结构化绑定**（根因 D）：AST 解析层的问题。
