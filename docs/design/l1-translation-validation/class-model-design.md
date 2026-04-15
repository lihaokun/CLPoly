# 可信基类建模 + 统一映射机制

> **状态**：已实现。13 个 C++ 函数翻译为 465 行 Lean，0 sorry，0 编译错误。

## 1. 设计原则

翻译器的所有 domain-specific 知识集中在 **`class_map.py`**（12 个配置表）和 **`clpoly_model.lean`**（可信基实现）中。翻译器逻辑代码（clang_ast / ssa_transform / lean_codegen / ub_collector）只查表，不包含任何 CLPoly 特有的硬编码。

**三层分离**：

```
class_map.py     ← 映射规则（C++ 操作 → Lean 函数名）
clpoly_model.lean ← 可信基实现（Lean 中的类型 + 正确语义）
翻译器逻辑        ← 只查表、只处理 AST 结构
```

新增类/方法/运算符只需改 `class_map.py` + `clpoly_model.lean`，不需改翻译器逻辑。

## 2. 配置表总览（class_map.py）

| 表名 | 用途 | 条目数 |
|------|------|--------|
| **CLASS_MAP** | 类方法/构造函数映射（Zp, ZZ, UMonomial, SparsePolyZp） | 4 类 × ~8 方法 |
| **FUNC_MAP** | 独立函数映射（derivative, GCD, pow 等） | 8 |
| **FIELD_MAP** | 成员字段名映射（first→fst, second→snd 等） | 5 |
| **OPERATOR_MAP** | C++ operator overload → Lean 运算符 | 13 |
| **CALL_OPERATOR_MAP** | 函数调用运算符 operator() → Lean 函数 | 1 |
| **STRUCT_COERCE_MAP** | StructType → BaseType 隐式转换 | 3 |
| **CAST_TABLE** | BaseType 间的类型转换（在 lean_codegen.py 中） | 13 |
| **UNSAFE_CAST_PAIRS** | 需要证明安全的类型转换 | 2 |
| **LEAN_BUILTINS** | Lean 标准库函数（不加 `_ir` 后缀） | ~20 |
| **TRANSLATION_SCOPE** | 翻译范围内的 C++ 函数（加 `_ir` 后缀） | 13 |
| **MUTATING_METHODS** | 修改 this 对象的方法（自动从 CLASS_MAP 派生） | 自动 |
| **EMPTY_CONTAINER_METHODS** | 空容器调用是 UB 的方法（自动从 CLASS_MAP 派生） | 自动 |
| **ASSERT_FAIL_NAMES** | assert 宏展开后的函数名 | 3 |

## 3. CLASS_MAP：类方法映射

```python
CLASS_MAP = {
    "Zp": {
        "lean_type": StructType("Zp", []),
        "constructors": {
            (BaseType.INT64, BaseType.UINT64): "Zp.ofInt",
            (BaseType.UINT64, BaseType.UINT64): "Zp.ofUInt64",
            (BaseType.INT64,): "Zp.ofInt",
            (): "default",
        },
        "methods": {
            "number": ("field", "val"),
            "prime": ("field", "prime"),
            "val": ("field", "val"),
            "inv": ("method", "Zp.inv"),
        },
        "operators": {
            "+": None, "-": None, "*": None,
            "/": "Zp.div",
            "==": None, "!=": None,
        },
    },
    "ZZ": { ... },        # 大整数 → Int
    "UMonomial": { ... },  # 单项式
    "SparsePolyZp": { ... }, # 稀疏多项式
}
```

### 3.1 方法类别

| 类别 | 含义 | Lean 生成 | SSA 处理 |
|------|------|----------|---------|
| `field` | 字段访问 | `obj.field_name` | 无 |
| `method` | 纯方法调用 | `lean_name obj` | 无 |
| `mutate` | 修改 this | `lean_name obj` | SSA 新版本 |
| `mutate_push` | push_back | `obj.push elem` | SSA 新版本 |
| `noop` | 无操作（reserve） | `obj` | 无 |
| `identity` | 透传（.data()） | `obj` | 无 |

### 3.2 构造函数匹配

CXXConstructExpr 的处理流程：
1. **Copy/move 构造检测**：1 个参数且类型 == 目标类型 → `parse_expr(arg)`（identity）
2. **CLASS_MAP 构造函数查找**：先按参数类型精确匹配，再按参数数量 fallback
3. **PairType** + 2 参数 → `Prod.mk`
4. **ArrayType** + 0 参数 → `#[]`
5. **BaseType** → 取最后一个参数
6. Fallback → `Prod.mk` 或 `parse_expr(arg)`

### 3.3 自动派生表

以下表从 CLASS_MAP 自动计算，不需要手动维护：

```python
NOOP_METHODS, MUTATING_METHODS = _derive_method_sets()
# 从 CLASS_MAP 的 methods 中按 category 分类

EMPTY_CONTAINER_METHODS = _derive_empty_container_methods()
# 从 CLASS_MAP 的 methods 中找 lean_name 含 "!" 的方法
```

## 4. FUNC_MAP：独立函数映射

```python
FUNC_MAP = {
    "derivative":     ("SparsePolyZp.derivative", "direct"),
    "polynomial_GCD": ("SparsePolyZp.gcd", "direct"),
    "pair_vec_div":   ("SparsePolyZp.divmod", "out2_drop1"),
    "get_deg":        ("SparsePolyZp.getDeg", "direct"),
    "move":           ("id", "identity"),
    "make_pair":      ("Prod.mk", "make_pair"),
    "pow":            ("HPow.hPow", "direct"),
}
```

参数规则：
- `"direct"` — 直接传参
- `"identity"` — 返回第一个实参（std::move）
- `"make_pair"` — Prod.mk
- `"out2_drop1"` — pair_vec_div 特殊处理（输出前 2 个，输入中间 2 个，丢弃最后 1 个）

## 5. OPERATOR_MAP：运算符映射

```python
OPERATOR_MAP = {
    "operator/": (3, "/"),  "operator*": (3, "*"),
    "operator+": (3, "+"),  "operator-": (3, "-"),
    "operator%": (3, "%"),
    "operator==": (3, "=="), "operator!=": (3, "!="),
    "operator<": (3, "<"),   "operator>": (3, ">"),
    "operator<=": (3, "<="), "operator>=": (3, ">="),
    "operator<<": (3, "<<"), "operator>>": (3, ">>"),
}
```

CXXOperatorCallExpr 处理流程：
1. `_extract_callee_name` 展开 ImplicitCastExpr 链提取 operator 名
2. OPERATOR_MAP 查表（按长度降序匹配，`operator==` 优先于 `operator=`）
3. CALL_OPERATOR_MAP 查表（`operator()` → `Rng.next`）
4. 复合赋值（`operator*=`）→ 展开为 `lhs = lhs * rhs`

## 6. STRUCT_COERCE_MAP：结构体类型转换

```python
STRUCT_COERCE_MAP = {
    ("Zp", BaseType.UINT64):     "{e}.val",
    ("Zp", BaseType.INT64):      "({e}.val.toNat : Int)",
    ("UMonomial", BaseType.UINT64): "{e}.deg",
}
```

当 Clang AST 中 ImplicitCastExpr 从 StructType 转换到 BaseType 时查此表。

## 7. 可信基实现（clpoly_model.lean）

### 7.1 文件结构

```
§1. Zp：Z/pZ 系数（算术 + extGcd 模逆）
§2. ZZ：大整数（abbrev Int）
§3. UMonomial：单变量单项式
§4. SparsePolyZp：稀疏多项式（derivative, gcd, divmod, comp）
§5. Array 辅助（front!, size_u64）
§5a. Rng：伪随机数生成器（xorshift64）
§6. 验证测试（#eval）
```

### 7.2 关键实现

- **Zp.ofInt**：正确处理负数（`v % p` 后检查 `< 0` 加 `p`）
- **Zp.modInv**：扩展欧几里得算法（`extGcdAux`）
- **SparsePolyZp.derivative**：正确的 Zp 系数求导
- **Rng.next**：xorshift64 伪随机（确定性，足够背靠背测试）
- **SparsePolyZp.comp**：占位（返回 0，pair_vec_div 的 Lean 实现不使用比较器）

### 7.3 隐式转换

```lean
instance : Coe Zp UInt64 where coe z := z.val
instance : Coe Zp Int where coe z := z.val.toNat
```

## 8. SSA 变换的统一机制

### 8.1 if/else phi 节点

当 if/else 两分支修改了同一个变量时，使用 **BlockExpr** 封装分支内的 let 链：

```lean
-- 正确：分支内变量在 Lean 的 if-then 作用域内可见
let result_4 := if cond then
    let result_2 := result_1 * b_1
    let result_3 := __upoly_mod_ir result_2 modpoly
    result_3
else result_1
```

BlockExpr（`ir_types.py`）= 语句列表 + 最终值表达式。codegen 输出为 Lean 的 let-in 链。

### 8.2 局部变量排除

循环变量识别（`identify_loop_vars`）和 phi 变量识别（`transform_if`）都排除**在当前作用域内首次声明的局部变量**：

```python
# 循环：排除循环体内 LetStmt 声明的变量
loop_vars = identify_loop_vars(body) - _collect_local_decls(body)

# if/else：只 phi 在 if 之前已存在的变量
diff = {k: v for k, v in diff.items() if k in pre_existing}
```

### 8.3 range-for 变换

`elem_let`（循环变量赋值 `let term := __coll[idx]!`）与循环体一起经过 `transform_body`，确保正确的 SSA 版本号。

### 8.4 成员赋值

`term.second *= lc_inv` 的处理：
1. `clang_ast.py`：检测 `operator*=` → 展开为 `term.second = term.second * lc_inv`
2. FieldAccess LHS 保留字段信息 → `ExprStmt(BinOp("=", FieldAccess(...), rhs))`
3. `ssa_transform.py`：`transform_member_assign` 生成 functional update + `Array.set!`

## 9. UB 证明目标（ub_collector.py）

### 9.1 检测规则

| UB 类型 | 检测位置 | 表驱动 |
|---------|---------|--------|
| DIV_BY_ZERO | BinOp `/` `%` | 运算符字面匹配 |
| SHIFT_OOB | BinOp `<<` `>>` | 运算符字面匹配 |
| SIGNED_OVERFLOW | BinOp `+` `-` `*` 且 signed 上下文 | `_is_signed_context` |
| ARRAY_OOB | ArrayAccess | 结构匹配 |
| EMPTY_CONTAINER | Call 到 EMPTY_CONTAINER_METHODS | 从 CLASS_MAP 派生 |
| UNSAFE_CAST | Cast 到 UNSAFE_CAST_PAIRS | 从 class_map.py 查表 |

### 9.2 注入过滤

`inject_obligations` 只将**引用参数名的 obligation** 提升到函数签名：

```python
def _all_vars_in_scope(lean_prop, param_names):
    # 提取标识符，去掉版本后缀，检查是否都是参数名
```

引用循环内部变量或函数体内 let 变量的 obligation 被过滤掉（它们应该作为循环不变式或局部 require 处理，不能提升到函数签名）。

## 10. 查找顺序

翻译器遇到 C++ 操作时，按以下统一顺序查找：

1. **Copy/move 构造** → identity（`_types_equal` 检测）
2. **CLASS_MAP** → 构造函数 / 方法 / 字段
3. **FUNC_MAP** → 独立函数
4. **OPERATOR_MAP** → 运算符重载（CXXOperatorCallExpr）
5. **CALL_OPERATOR_MAP** → 函数调用运算符（operator()）
6. **STRUCT_COERCE_MAP** → 结构体→基本类型隐式转换
7. **CAST_TABLE** → 基本类型间转换
8. **LEAN_BUILTINS** → Lean 标准库函数
9. **TRANSLATION_SCOPE** → 翻译范围内函数（加 `_ir` 后缀）
10. **都不匹配 → `sorry`** + 注释标注

## 11. 验证结果

| 指标 | 值 |
|------|---|
| 翻译函数数 | 13 |
| 生成 Lean 行数 | 465 |
| sorry 数 | 0 |
| UNKNOWN/CAST 标记数 | 0 |
| Lean 编译错误数 | 0 |
| 可信基 Lean 模型行数 | ~150 |
| class_map.py 配置表数 | 12 |
