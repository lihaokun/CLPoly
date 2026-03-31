# 表达式类型传播设计

> 解决翻译器的核心类型问题：Lean 编译时 55 个错误。

## 1. 问题

当前 `ExprIR` 节点不携带类型信息。翻译器在生成 Lean 代码时不知道每个表达式的类型，导致：

1. **Zp vs UInt64 混淆**：`f.front().second` 的类型是 `Zp`，但后续代码把它当 `UInt64` 用
2. **构造函数未识别**：`Zp(val, p)` 在 Clang 中是 `CXXConstructExpr`，翻译器没区分 pair 构造和自定义类型构造
3. **Array.size 类型**：Lean 返回 `Nat`，CLPoly 期望 `UInt64`

## 2. 根因

违反 P1 原则：**翻译器不做类型推断——所有类型信息从 Clang AST 传播。**

Clang AST 的每个表达式节点都有 `type.qualType`，但 `ExprIR` 没有存储。

## 3. 错误分类

55 个错误的实际构成：

| 错误类型 | 数量 | 类型传播能修？ |
|---------|------|-------------|
| Type mismatch | 14 | ✅ |
| Application type mismatch | 9 | ✅ |
| Function expected at | 3 | ❌ 结构问题 |
| unexpected token | 5 | ❌ 语法问题 |
| 级联错误 | 24 | ⚠️ 前面修好后大部分消失 |

**类型传播预期修 23 个（Type + Application mismatch），残留 ~5-8 个结构/语法问题需单独修。**

## 4. 方案

### 4.1 给 ExprIR 加类型字段（不用 monkey-patch）

在 `ir_types.py` 中，给**所有 ExprIR dataclass** 加 `_ast_type` 字段：

```python
@dataclass
class Var:
    name: str
    version: int = 0
    _ast_type: TypeIR = None    # Clang AST 类型（新增）

@dataclass
class BinOp:
    op: str
    lhs: ExprIR
    rhs: ExprIR
    _ast_type: TypeIR = None    # 新增

# ... 所有 ExprIR 类型都加
```

**理由**：monkey-patch（`result._ast_type = ...`）不类型安全、不可见于定义。dataclass 字段是正确做法。

### 4.2 parse_expr 附加类型

```python
def parse_expr(node: dict) -> ExprIR:
    result = _parse_expr_inner(node)
    ast_qual = node.get("type", {}).get("qualType", "")
    if ast_qual:
        result._ast_type = parse_type(ast_qual)
    return result
```

### 4.3 合成节点赋类型

`ssa_transform` 中创建的合成节点必须手动赋 `_ast_type`：

| 合成场景 | 类型来源 |
|---------|---------|
| phi 节点 `CondExpr(c, a, b)` | 从 then 分支表达式 `a._ast_type` 获取 |
| 步进 `BinOp("+", var, Lit(1))` | 从 `var._ast_type` 获取（同类型加法） |
| `Call("Array.set!", [arr, idx, elem])` | 从 `arr._ast_type` 获取 |
| `Call("_with", [root, field, val])` | 从 `root._ast_type` 获取 |
| `ArrayPush(arr, elem)` | 从 `arr._ast_type` 获取 |

规则：**合成节点的类型 = 操作数中最具体的类型**。

### 4.4 rename_expr 传播类型

```python
def rename_expr(expr, env):
    if isinstance(expr, Var):
        new_var = env.current(expr.name)
        new_var._ast_type = expr._ast_type    # 传播！
        return new_var
    if isinstance(expr, BinOp):
        result = BinOp(expr.op, rename_expr(expr.lhs, env), rename_expr(expr.rhs, env))
        result._ast_type = expr._ast_type     # 传播！
        return result
    # ... 所有 rename 分支都传播 _ast_type
```

### 4.5 gen_coercion 统一转换

```python
def gen_coercion(expr_str: str, source: TypeIR, target: TypeIR) -> str:
    """source 类型表达式 → target 类型。"""
    if source == target:
        return expr_str

    # BaseType → BaseType：查 CAST_TABLE
    if isinstance(source, BaseType) and isinstance(target, BaseType):
        key = (source, target)
        if key in CAST_TABLE:
            return CAST_TABLE[key].format(e=expr_str)

    # StructType → BaseType
    if isinstance(source, StructType) and isinstance(target, BaseType):
        if source.name == "Zp" and target == BaseType.UINT64:
            return f"{expr_str}.val"
        if source.name == "UMonomial" and target == BaseType.UINT64:
            return f"{expr_str}.deg"

    # BaseType → StructType
    if isinstance(source, BaseType) and isinstance(target, StructType):
        if target.name == "Zp":
            return f"(Zp.mk {expr_str} 0)"  # prime 需从上下文获取

    # fallback
    return f"({expr_str} : {gen_type(target)})"
```

### 4.6 gen_stmt(LetStmt) 自动 coercion

```python
def gen_stmt(stmt):
    if isinstance(stmt, LetStmt):
        declared = stmt.typ
        actual = getattr(stmt.value, '_ast_type', None)
        val = gen_expr(stmt.value)
        if actual and actual != declared:
            val = gen_coercion(val, actual, declared)
        ...
```

### 4.7 CXXConstructExpr 按目标类型区分

```python
if kind == "CXXConstructExpr":
    typ = parse_type(qualType)
    if isinstance(typ, StructType) and typ.name == "Zp" and len(inner) == 2:
        return Call("Zp.mk", [parse_expr(inner[0]), parse_expr(inner[1])])
    if isinstance(typ, PairType) and len(inner) == 2:
        return Call("Prod.mk", [parse_expr(inner[0]), parse_expr(inner[1])])
    if isinstance(typ, StructType) and typ.name == "SparsePolyZp":
        return Call("Array.empty", [])
    ...
```

## 5. 可行性论证

### 5.1 完整性

| 表达式来源 | 有 `_ast_type`？ | 保证机制 |
|-----------|-----------------|---------|
| `parse_expr`（原始节点） | ✅ | §4.2：从 Clang `type.qualType` |
| `rename_expr`（SSA 重命名） | ✅ | §4.4：从原节点传播 |
| `transform_if` phi | ✅ | §4.3：从操作数推断 |
| `transform_step` 步进 | ✅ | §4.3：从计数器类型 |
| 其他合成节点 | ✅ | §4.3：从操作数推断 |

**结论：所有 ExprIR 节点都有 `_ast_type`。**

### 5.2 正确性

**命题**：`gen_coercion(val, source, target)` 产出的 Lean 表达式类型为 `target`。

**证明**（逐 case）：
- `source == target`：返回 `val`（类型不变）✓
- `CAST_TABLE` 命中：CAST_TABLE 的每个条目已验证（§blueprint 引理 1）✓
- `Zp → UInt64`：`val.val` 是 `Zp` 结构体的 `UInt64` 字段 ✓
- `UMonomial → UInt64`：`val.deg` 是 `UMonomial` 的 `UInt64` 字段 ✓
- fallback `(val : target)`：Lean 类型标注，如果没有 Coe 实例则编译报错（不会静默错误）✓

**结论：gen_coercion 要么产出正确类型，要么 Lean 编译报错（不会产出错误类型的代码）。**

### 5.3 预期效果

| 错误类型 | 数量 | 修复 |
|---------|------|------|
| Type mismatch | 14 | → 0（gen_coercion） |
| Application type mismatch | 9 | → 0（参数 coercion） |
| Function expected at | 3 | 不在范围，需单独修 |
| unexpected token | 5 | 不在范围，需单独修 |
| 级联错误 | 24 | 大部分消失（前面错误修好后） |
| **总计** | **55** | **→ ~5-8 残留** |

### 5.4 不在范围的问题

以下问题需单独修，不属于类型传播：

1. **"Function expected at"**（3 个）：通常是 `let rec` 返回值结构问题
2. **"unexpected token"**（5+2 个）：if/else 嵌套或 let 作用域问题

## 6. 实施步骤

1. `ir_types.py`：所有 ExprIR dataclass 加 `_ast_type: TypeIR = None`
2. `clang_ast.py`：`parse_expr` 附加类型 + `CXXConstructExpr` 按类型区分
3. `ssa_transform.py`：`rename_expr` 传播 + 合成节点赋类型
4. `lean_codegen.py`：`gen_coercion` + `gen_stmt(LetStmt)` 自动 coercion
5. 编译验证：55 → ~5-8
6. 残留的结构/语法错误单独修

## 7. 前提条件

| 前提 | 验证 |
|------|------|
| Clang AST 每个表达式有 `type.qualType` | Clang 标准保证（编译通过的代码） |
| CAST_TABLE 覆盖所有 BaseType 对 | 已验证（blueprint §4.1） |
| Zp/UMonomial 结构体字段名正确 | 由 prelude 定义保证 |
| rename_expr 覆盖所有 ExprIR 类型 | 代码审查确认（12 种 ExprIR 全有分支） |
