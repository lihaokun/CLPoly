# HIR 数据结构与 Pass 1-5 设计

> Stage 1 Week 4 主产出
> 基于：`cpp-construct-catalog.md`（Week 1）+ `type-system.md`（Week 2）+ `cpp-subset-semantics.md`（Week 3）+ `v1-reuse-inventory.md`（Week 4 Agent）
> 日期：2026-04-21

---

## §0 导读

### §0.1 目标

定义 cpp2lean v2 的 HIR（High-level IR）数据结构 + 5 个 HIR Pass 的精确规格。每个 Pass 的输入/输出/不变量/转换规则都要明确，为 Stage 2 实现提供 1:1 可编程的蓝图。

### §0.2 非目标

- **不定义 MIR**：MIR 是 SSA 形式 + CFG，留给 Week 5
- **不写具体 Python 代码**：只定义 dataclass schema 和 Pass 算法；Stage 2 负责具体实现
- **不处理形式证明**：Week 3 的 `cpp-subset-semantics.md` 已给出非形式引理，Stage 4 (精化证明) 再考虑 Lean 化

### §0.3 整体架构回顾

```
AST (Clang JSON)
    ↓ Pass 1: parse
HIR₀: 原始 IR（允许一切）
    ↓ Pass 2: ref_elim
HIR₁: 无 T& 参数
    ↓ Pass 3: lambda_lift
HIR₂: 无 inline Lambda
    ↓ Pass 4: iter_recognize
HIR₃: 无裸迭代器
    ↓ Pass 5: operator_resolve
HIR₄: 所有 Call 的 callee 已解析
    ↓ (Week 5) Pass 6: ssa_build
MIR₀: SSA + phi nodes
    ↓ (Week 5) Pass 7: loop_lower
MIR₁: 循环已提取，break/continue 已下降
    ↓ (Week 5) Pass 8: codegen
Lean 源代码
```

每个 HIR₀-HIR₄ 共用 **同一个 dataclass 家族**，但**不变量**递进。每 Pass 的入口/出口有 runtime assert 验证不变量。

---

## §1 HIR 数据结构

### §1.1 类型系统（TypeIR）

```python
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Union

class BaseType(Enum):
    """基础数值类型，对应 type-system.md §1."""
    UINT64 = "UInt64"
    INT64 = "Int64"
    INT32 = "Int32"
    UINT32 = "UInt32"
    NAT = "Nat"           # size_t → Nat
    BOOL = "Bool"
    FLOAT = "Float"        # double
    UNIT = "Unit"          # void

@dataclass
class NamedType:
    """命名类型（含模板实例化后的完整名）。"""
    name: str              # "Zp", "ZZ", "SparsePolyZp", "MvPolyZZ", ...
    # 不展开模板参数：实例化后的类型名就是完整标识

@dataclass
class ArrayType:
    elem: 'TypeIR'

@dataclass
class PairType:
    fst: 'TypeIR'
    snd: 'TypeIR'

@dataclass
class TupleType:
    """多元组，>2 元素。"""
    elems: list['TypeIR']

@dataclass
class OptionType:
    """Option α = none | some α；主要由 `StdMap.find` 等返回。"""
    inner: 'TypeIR'

@dataclass
class StdMapType:
    key: 'TypeIR'
    value: 'TypeIR'

@dataclass
class UnknownType:
    """占位，`parse` Pass 偶尔产出；后续 Pass 不允许。"""
    raw: str

TypeIR = Union[BaseType, NamedType, ArrayType, PairType, TupleType,
               OptionType, StdMapType, UnknownType]
```

**设计要点**：
- 不区分 `const T` 和 `T`：Lean 值语义透明（`cpp-subset-semantics.md` §4.2）
- 不表示 `T&` / `T*`：`ref_elim` Pass 之前用 `RefType` 标记，Pass 2 后消除；本文件 §1.4 定义
- 不处理模板参数：实例化后的类型名已具体化（`upolynomial_<Zp>` → `NamedType("SparsePolyZp")`）

### §1.2 引用标记（仅 HIR₀）

```python
@dataclass
class RefType:
    """C++ T& 或 T&&，HIR₀ 允许，ref_elim Pass 后消除。"""
    inner: 'TypeIR'
    is_const: bool = False
    is_rvalue: bool = False
```

**不变量**：
- HIR₀：参数类型可出现 `RefType`
- HIR₁ 及之后：**无 `RefType`**（被 `ref_elim` 转为 tuple 返回值）

### §1.3 表达式（ExprIR）

```python
@dataclass
class Var:
    """变量引用。HIR 阶段 version 始终为 0；MIR 阶段才有 SSA 版本号。"""
    name: str
    version: int = 0
    ty: Optional['TypeIR'] = None  # 从 AST 传播

@dataclass
class Lit:
    """字面量。"""
    value: int | bool | float | str
    ty: 'TypeIR' = field(default=BaseType.INT32)

@dataclass
class BinOp:
    """二元运算（基础类型）。"""
    op: str                 # "+" "-" "*" "/" "%" "<" ">" "==" "!=" "&&" "||" ...
    lhs: 'ExprIR'
    rhs: 'ExprIR'
    ty: Optional['TypeIR'] = None

@dataclass
class UnaryOp:
    """一元运算。"""
    op: str                 # "!" "-" "++" "--" "~"
    operand: 'ExprIR'
    ty: Optional['TypeIR'] = None

@dataclass
class CondExpr:
    """三元 ? :"""
    cond: 'ExprIR'
    then_e: 'ExprIR'
    else_e: 'ExprIR'
    ty: Optional['TypeIR'] = None

@dataclass
class Call:
    """函数/方法调用。
    callee 形态：
      HIR₀-HIR₃: 可为未解析的 raw 字符串（如 "operator[]" / "std::max"）
      HIR₄: 必须是已解析的 Lean 函数名（如 "Array.get!" / "max"）
    """
    callee: str | 'Expr' | 'UnresolvedOp'
    args: list['ExprIR']
    ty: Optional['TypeIR'] = None

@dataclass
class UnresolvedOp:
    """HIR₀-HIR₃：未解析的运算符或方法调用。
    operator_resolve Pass 把它替换为具体 Call."""
    op_name: str           # "operator[]", "operator+", "size"
    receiver_ty: Optional['TypeIR'] = None

@dataclass
class ArrayAccess:
    """arr[i]"""
    arr: 'ExprIR'
    idx: 'ExprIR'
    ty: Optional['TypeIR'] = None

@dataclass
class FieldAccess:
    """obj.field_name"""
    obj: 'ExprIR'
    field_name: str
    ty: Optional['TypeIR'] = None

@dataclass
class Cast:
    """显式/隐式类型转换。"""
    expr: 'ExprIR'
    source_ty: 'TypeIR'
    target_ty: 'TypeIR'
    cast_kind: str          # "IntegralCast" / "NoOp" / "LValueToRValue" 等

@dataclass
class LambdaExpr:
    """HIR₀-HIR₁ 允许的内联 Lambda。lambda_lift Pass 之后消除。"""
    captures: list['Capture']      # 捕获变量
    params: list['HIRParam']
    body: list['StmtIR']
    ty: Optional['TypeIR'] = None

@dataclass
class Capture:
    """Lambda 捕获条目。"""
    name: str
    by_ref: bool            # [&x] True, [=x] False
    is_default: bool        # [&] True（不是具名），[x] False

@dataclass
class IteratorExpr:
    """HIR₀-HIR₂ 原始迭代器（`v.begin()`, `v.end()`, `*it` 等）。
    iter_recognize Pass 识别模式后消除。"""
    kind: str               # "begin" / "end" / "deref" / "increment"
    container: Optional['ExprIR'] = None
    operand: Optional['ExprIR'] = None

@dataclass
class BlockExpr:
    """语句块作为表达式值（phi 节点材料）。"""
    stmts: list['StmtIR']
    value: 'ExprIR'
    ty: Optional['TypeIR'] = None

@dataclass
class TupleExpr:
    """tuple 构造（ref_elim 后可能产生）。"""
    elems: list['ExprIR']
    ty: Optional['TypeIR'] = None

@dataclass
class ArrayLit:
    """Array 字面量 #[...]"""
    elems: list['ExprIR']
    elem_ty: 'TypeIR'

@dataclass
class UnknownExpr:
    """parse Pass 的未识别表达式，后续 Pass 必须在报错前消除或通报。"""
    kind: str
    children: list['ExprIR'] = field(default_factory=list)
    raw: str = ""

ExprIR = Union[Var, Lit, BinOp, UnaryOp, CondExpr, Call, UnresolvedOp,
               ArrayAccess, FieldAccess, Cast, LambdaExpr, Capture,
               IteratorExpr, BlockExpr, TupleExpr, ArrayLit, UnknownExpr]
```

**简化设计**（相对 v1 ir_types.py）：
- 删除 `ArrayPush`：v1 用于表示 `.push_back()`，v2 统一用 `Call(callee="Array.push", ...)` 表达
- 删除 `ExceptType`：CLPoly 不用异常（`cpp-subset-semantics.md` §7.3 确认）
- 新增 `UnresolvedOp` / `IteratorExpr` / `LambdaExpr` / `Capture`：各自 Pass 中转节点
- 新增 `TupleExpr` / `ArrayLit`：`ref_elim`、`initializer_list` 产物

### §1.4 语句（StmtIR）

```python
@dataclass
class LetStmt:
    """let x : T := e"""
    var: 'Var'
    ty: 'TypeIR'
    value: 'ExprIR'

@dataclass
class AssignStmt:
    """x = e（HIR 允许；MIR 阶段被 SSA 消除）"""
    target: 'ExprIR'        # 可以是 Var、ArrayAccess、FieldAccess
    value: 'ExprIR'

@dataclass
class CompoundAssignStmt:
    """x op= e （+=, -=, *=, /= 等；HIR₀-HIR₃ 允许）
    operator_resolve Pass 展开为 AssignStmt + BinOp。"""
    target: 'ExprIR'
    op: str                 # "+", "-", "*", "/", "%"
    value: 'ExprIR'

@dataclass
class IfStmt:
    cond: 'ExprIR'
    then_body: list['StmtIR']
    else_body: list['StmtIR'] = field(default_factory=list)

@dataclass
class WhileStmt:
    cond: 'ExprIR'
    body: list['StmtIR']

@dataclass
class ForStmt:
    """C++ for (init; cond; step) { body }"""
    init: list['StmtIR']
    cond: 'ExprIR'
    step: list['StmtIR']
    body: list['StmtIR']

@dataclass
class RangeForStmt:
    """C++ for (auto& x : container) { body }"""
    var: 'Var'
    var_ty: 'TypeIR'
    container: 'ExprIR'
    body: list['StmtIR']
    # 结构化绑定 [k, v]
    decomposition: Optional[list['Var']] = None

@dataclass
class DoWhileStmt:
    body: list['StmtIR']
    cond: 'ExprIR'

@dataclass
class BreakStmt:
    pass

@dataclass
class ContinueStmt:
    pass

@dataclass
class ReturnStmt:
    value: Optional['ExprIR'] = None

@dataclass
class RequireStmt:
    """assert / UB 前置条件。由 parse Pass (assert 识别) 或
    operator_resolve/ref_elim (UB) 生成。放到函数签名里。"""
    cond: 'ExprIR'
    name: str               # "hp" / "h_shift" / "hi"
    source: str             # "assert" / "div_by_zero" / "array_oob" / ...
    uid: int = 0            # 唯一标识，用于 def-use 跟踪

@dataclass
class ExprStmt:
    """表达式作为语句（丢弃结果），如 `f();` 或 `i++;`"""
    expr: 'ExprIR'

@dataclass
class BlockStmt:
    """显式的 { ... }，用于 scope 作用域。"""
    stmts: list['StmtIR']

@dataclass
class UnknownStmt:
    """parse Pass 的未识别语句。后续 Pass 遇到必须报错。"""
    kind: str
    children: list['StmtIR'] = field(default_factory=list)

StmtIR = Union[LetStmt, AssignStmt, CompoundAssignStmt,
               IfStmt, WhileStmt, ForStmt, RangeForStmt, DoWhileStmt,
               BreakStmt, ContinueStmt, ReturnStmt,
               RequireStmt, ExprStmt, BlockStmt, UnknownStmt]
```

**设计要点**：
- 保留 `C++ for/while/range-for/do-while` 四种循环：Pass 7 `loop_lower` 在 MIR 阶段统一下降为 `partial def` 尾递归，HIR 阶段仍然结构化
- `AssignStmt` 和 `CompoundAssignStmt` 在 HIR 合法；Pass 6 `ssa_build` 消除
- `RangeForStmt` 包含可选 `decomposition`（结构化绑定的 30 处）

### §1.5 顶层结构

```python
@dataclass
class HIRParam:
    name: str
    ty: 'TypeIR'
    is_ref: bool = False        # HIR₀: T&; HIR₁+: 必须 False
    is_const_ref: bool = False  # HIR₀: const T&; 后续保留
    is_output: bool = False     # HIR₀: 由 parse Pass 标注，非 const ref

@dataclass
class HIRFunc:
    """单个 Lean 定义对应的函数 IR。"""
    base_name: str              # 如 "factorize"
    instance_suffix: str        # "upoly" / "lex" / "grlex" / ""
    mangled_name: str           # Clang mangled name
    qual_type: str              # 完整 C++ 签名（debug 用）
    params: list['HIRParam']
    ret_ty: 'TypeIR'
    body: list['StmtIR']
    requires: list['RequireStmt'] = field(default_factory=list)
    # Pass 3 lambda_lift 后，该函数依赖的 lambda 独立 def
    aux_lambdas: list['HIRFunc'] = field(default_factory=list)

    @property
    def lean_name(self) -> str:
        """生成的 Lean 函数名。"""
        if self.instance_suffix:
            return f"{self.base_name}_{self.instance_suffix}_ir"
        else:
            return f"{self.base_name}_ir"

@dataclass
class HIRProgram:
    """完整程序，包含所有翻译的函数。"""
    funcs: list['HIRFunc']
```

**与 v1 `FuncIR` 对比**：
- v1 `FuncIR.params` 用 `ParamIR(is_output: bool)`，但缺 `instance_suffix`/`mangled_name`
- v1 `SSAFunc` 和 `FuncIR` 分离；v2 用单一 `HIRFunc`，SSA 形式通过 Var.version > 0 表达
- v1 `aux_defs` 是"提取的循环函数"，v2 移到 MIR Pass 7 产出；HIR 层只有 `aux_lambdas`

---

## §2 HIR 不变量阶梯

每个 HIR_i 共享同一套 dataclass，但不变量递进。入口/出口用 `assert_invariant_i()` 函数验证。

### §2.1 不变量表

| 阶段 | 不变量 | 由什么 Pass 保证 |
|---|---|---|
| **HIR₀** | 原始——允许一切。`UnknownStmt`/`UnknownExpr` 可能出现（记录未识别节点）| Pass 1 `parse` |
| **HIR₁** | **无 `T&`/`T&&` 参数**（`HIRParam.is_ref == False`）；所有 ref 参数已转 tuple | Pass 2 `ref_elim` |
| **HIR₂** | **无 inline `LambdaExpr`**；所有 lambda 已提升为 `HIRFunc` in `aux_lambdas` | Pass 3 `lambda_lift` |
| **HIR₃** | **无裸 `IteratorExpr`**；所有迭代器模式已转为高阶 Call（`Array.filter`/`Array.foreach`）或显式 index loop | Pass 4 `iter_recognize` |
| **HIR₄** | **所有 `Call.callee` 为具体函数名字符串**；无 `UnresolvedOp`；所有 `CompoundAssignStmt` 已展开；所有 `Cast` 的 source/target 可查 `CAST_TABLE` | Pass 5 `operator_resolve` |

### §2.2 runtime assert 函数

```python
def assert_hir1_invariant(func: HIRFunc):
    """HIR₁ 出口 assert："""
    for p in func.params:
        assert not p.is_ref, f"{func.base_name}: param {p.name} still has ref after ref_elim"

def assert_hir2_invariant(func: HIRFunc):
    """HIR₂ 出口 assert："""
    _no_lambda = True
    def walk(node):
        nonlocal _no_lambda
        if isinstance(node, LambdaExpr):
            _no_lambda = False
        if hasattr(node, '__dict__'):
            for v in node.__dict__.values():
                if isinstance(v, list):
                    for x in v: walk(x)
                elif hasattr(v, '__dict__'):
                    walk(v)
    for stmt in func.body:
        walk(stmt)
    assert _no_lambda, f"{func.base_name}: inline LambdaExpr after lambda_lift"

# HIR₃, HIR₄ 类似
```

### §2.3 违反不变量 → 报错

v1 的一大问题是**默认兜底**（不认识就 `let _ := expr`）。v2 改为**严格报错**：

```python
# 每个 Pass 遇到 UnknownStmt/UnknownExpr/未注册构造：
raise TranslationError(f"{pass_name}: unknown construct {kind} in {func.base_name}:{line}")
```

`UnknownStmt`/`UnknownExpr` 在 HIR₀ 是警告（parse Pass 记录），HIR₁ 开始是**错误**。

---

## §3 Pass 1: `parse` (AST → HIR₀)

### §3.1 职责

把 Clang AST JSON（已确保单态化，见 Week 1 Day 1）转为 HIR₀ 节点。

**前提**：AST 是经过 `TRANSLATION_INSTANCES` 配置（§5.4.5 `type-system.md`）挑选的**实例化 FunctionDecl**（mangledName ≠ None）。不处理模板定义。

### §3.2 输入/输出

- **输入**：`ast_json : dict`（单个 FunctionDecl 节点）+ `TRANSLATION_INSTANCES : dict`（全表）
- **输出**：`HIRFunc`（可能含 `UnknownStmt`/`UnknownExpr`）

### §3.3 转换规则（AST kind → HIR 节点）

来自 `cpp-construct-catalog.md` §13 的 Pass-AST 表，`parse` 需处理全部 62 种 kind：

**语句类**：
| AST kind | HIR 节点 |
|---|---|
| `CompoundStmt` | `BlockStmt` |
| `DeclStmt` + `VarDecl` | `LetStmt(var, ty, init)` |
| `BinaryOperator::=` | `AssignStmt(target, value)` |
| `CompoundAssignOperator::*=` | `CompoundAssignStmt(target, "*", value)` |
| `IfStmt` | `IfStmt(cond, then_body, else_body)` |
| `WhileStmt` | `WhileStmt(cond, body)` |
| `ForStmt` | `ForStmt(init, cond, step, body)` |
| `CXXForRangeStmt` | `RangeForStmt(var, ty, container, body, decomposition)` |
| `DoStmt` | `DoWhileStmt(body, cond)` |
| `BreakStmt` | `BreakStmt()` |
| `ContinueStmt` | `ContinueStmt()` |
| `ReturnStmt` | `ReturnStmt(value)` |
| `CallExpr(__assert_fail)` 在 if 内 | `RequireStmt(cond, ...)` |
| `ExprStmt` / 其他 | `ExprStmt(expr)` |
| `DecompositionDecl` | 作为 `RangeForStmt.decomposition` 字段 |
| 其他未识别 | `UnknownStmt(kind, ...)` |

**表达式类**：
| AST kind | HIR 节点 |
|---|---|
| `DeclRefExpr` | `Var(name)` |
| `IntegerLiteral` / `FloatingLiteral` / `CXXBoolLiteralExpr` | `Lit(value, ty)` |
| `BinaryOperator::<` 等 | `BinOp(op, lhs, rhs)` |
| `UnaryOperator::!` 等 | `UnaryOp(op, operand)` |
| `ConditionalOperator` | `CondExpr(cond, then, else)` |
| `CallExpr` | `Call(callee_str, args)` — callee 先存字符串（raw）|
| `CXXOperatorCallExpr` | `Call(UnresolvedOp(op_name), args)` |
| `CXXMemberCallExpr` | `Call(UnresolvedOp(f"{obj_ty}.{method}"), [obj] + args)` |
| `ArraySubscriptExpr` | `ArrayAccess(arr, idx)` |
| `MemberExpr` | `FieldAccess(obj, field_name)` |
| `ImplicitCastExpr` / `CStyleCastExpr` / `CXXStaticCastExpr` | `Cast(expr, source_ty, target_ty, cast_kind)` |
| `CXXConstructExpr` | `Call(UnresolvedOp(f"construct_{ty_name}"), args)` |
| `LambdaExpr` | `LambdaExpr(captures, params, body)` — captures 从源码正则提取 |
| `StringLiteral` / `SourceLocExpr` / `PredefinedExpr` | 忽略（assert 消息）|
| 其他 | `UnknownExpr(kind, children, raw)` |

### §3.4 类型传播（P1 原则）

**严格从 AST 读取，不推断**。每个 AST 节点有 `type.qualType`，直接转为 `TypeIR`（通过 `TYPE_PARSE_TABLE`）。

```python
def parse_type(qt: str) -> TypeIR:
    # 基础
    if qt in ("uint64_t", "unsigned long"): return BaseType.UINT64
    if qt in ("int64_t", "long"): return BaseType.INT64
    # ...
    # 数组
    if qt.startswith("std::vector<"):
        elem = parse_type(extract_template_arg(qt, 0))
        return ArrayType(elem)
    # pair / map
    if qt.startswith("std::pair<"): return PairType(...)
    if qt.startswith("std::map<"): return StdMapType(...)
    # CLPoly 类型
    if qt in ("upolynomial_<ZZ>", "basic_polynomial<umonomial, ZZ, uless>"):
        return NamedType("SparsePolyZZ")
    # ... 依 type-system.md 完整表
    return UnknownType(qt)
```

### §3.5 参数分类

对每个 `ParmVarDecl`，根据 qualType 判断：

```python
def parse_param(parm_json) -> HIRParam:
    qt = parm_json["type"]["qualType"]
    is_const_ref = qt.startswith("const ") and qt.endswith("&")
    is_ref = (qt.endswith("&") or qt.endswith("*")) and not is_const_ref
    inner = qt.replace("const ", "").rstrip("&*").strip()
    return HIRParam(
        name=parm_json["name"],
        ty=parse_type(inner),
        is_ref=is_ref,
        is_const_ref=is_const_ref,
        is_output=is_ref  # 初始猜测；ref_elim Pass 验证
    )
```

### §3.6 错误处理

- **允许 `UnknownStmt`/`UnknownExpr`**：记录 kind + 子节点，供人工检查
- 结束时打印 "N unknown nodes in func X"
- Pass 2-5 遇到 Unknown 节点**报错退出**

### §3.7 实现规模

预估 ~400 行 Python（参考 v1 `clang_ast.py` 1287 行中 ~400 行可搬）。

---

## §4 Pass 2: `ref_elim` (HIR₀ → HIR₁)

### §4.1 职责

消除 C++ `T&` / `T&&` 参数，转为 tuple 返回值。依据 `type-system.md` §7 + `cpp-subset-semantics.md` §4.3-4.4（引理 L4.1）。

### §4.2 输入/输出

- **输入**：`func : HIRFunc`（HIR₀）
- **输出**：同 func 改写签名 + body 的最后一条 `ReturnStmt`

### §4.3 算法

```
ref_elim(func):
    ref_params = [p for p in func.params if p.is_ref]
    if not ref_params:
        func.params[*].is_ref = False  # const ref 保留 is_const_ref
        return func
    
    # 1. 改写参数签名：移除 ref 参数，每个作为普通参数传入初值
    new_params = [p for p in func.params if not p.is_ref]
    for rp in ref_params:
        rp.is_ref = False
        new_params.append(rp)  # 放到参数末尾作为"初值"
    
    # 2. 改写返回类型：原返回 + ref 参数的 tuple
    ref_tys = [p.ty for p in ref_params]
    if func.ret_ty == BaseType.UNIT:
        new_ret_ty = single_or_tuple(ref_tys)
    else:
        new_ret_ty = TupleType([func.ret_ty] + ref_tys)
    
    # 3. 改写 body 末尾的 ReturnStmt
    new_body = rewrite_returns(func.body, ref_params)
    # 所有 return 点追加 ref_params 的当前值（HIR₀ 尚无 SSA，用原 name）
    
    # 4. 改写调用方（在 Pass 结束后对整个 HIRProgram 统一改写）
    #    但由于翻译单元 = 一个 HIRFunc，调用方改写放到 Pass 5 或单独的 call_rewrite 步骤
    
    func.params = new_params
    func.ret_ty = new_ret_ty
    func.body = new_body
    return func

rewrite_returns(body, ref_params):
    """在每个 return 点追加 ref 参数的当前值作为 tuple。"""
    result = []
    for stmt in body:
        if isinstance(stmt, ReturnStmt):
            # 构造 tuple
            orig_val = stmt.value or Lit(None, BaseType.UNIT)
            elems = [orig_val] + [Var(p.name) for p in ref_params]
            result.append(ReturnStmt(TupleExpr(elems)))
        elif isinstance(stmt, IfStmt):
            result.append(IfStmt(
                stmt.cond,
                rewrite_returns(stmt.then_body, ref_params),
                rewrite_returns(stmt.else_body, ref_params)
            ))
        elif isinstance(stmt, (WhileStmt, ForStmt, RangeForStmt, DoWhileStmt, BlockStmt)):
            # 递归到 body
            stmt.body = rewrite_returns(stmt.body, ref_params)
            result.append(stmt)
        else:
            result.append(stmt)
    
    # 若 body 末尾无 return（void 函数），追加一条 return (ref_params_tuple)
    if not result or not isinstance(result[-1], ReturnStmt):
        elems = [Var(p.name) for p in ref_params]
        result.append(ReturnStmt(
            elems[0] if len(elems) == 1 else TupleExpr(elems)
        ))
    return result
```

### §4.4 调用方改写

由于 HIR 阶段每个 `HIRFunc` 独立处理，**调用方改写推迟到 Pass 5 `operator_resolve`**（那里需要全局重新解析 Call）。

具体：在 Pass 5 中，看到 `Call(callee="foo_ir", args=[...])` 时，若 `foo_ir` 在 `TRANSLATION_INSTANCES` 表里有 `output_count > 0`，则：
```
old:  Call("foo_ir", args)
new:  LetStmt(tmp, ret_ty_tuple, Call("foo_ir", args_with_initial_ref_values))
      + 对 tmp 解构到各 output 变量
```

### §4.5 输出参数检测策略

**v2 改进**：不再依赖手工维护的 `TRANSLATION_SCOPE_OUTPUT_PARAMS`（v1 8 处配置错位）。每 `HIRFunc` 的参数表直接给出 ref 参数索引。

```python
# 运行时用的 table，从 HIRProgram 的所有 HIRFunc 动态构建
OUTPUT_PARAMS = {
    func.base_name + (f"_{func.instance_suffix}" if func.instance_suffix else ""): 
        [i for i, p in enumerate(func.params) if p.is_ref]
    for func in program.funcs
}
```

### §4.6 语义保持

见 `cpp-subset-semantics.md` §4.4 引理 L4.1。

### §4.7 实现规模

~200 行 Python（v1 的输出参数识别 + tuple 改写）。

---

## §5 Pass 3: `lambda_lift` (HIR₁ → HIR₂)

### §5.1 职责

把所有 `LambdaExpr` 提升为独立 `HIRFunc`（存入 `aux_lambdas`），调用点替换为对新函数的显式调用。依据 `cpp-subset-semantics.md` §5。

### §5.2 输入/输出

- **输入**：`func : HIRFunc`（HIR₁）
- **输出**：同 func + body 内所有 inline Lambda 被引用替换 + `aux_lambdas` 列表填充

### §5.3 算法

```
lambda_lift(func):
    lambda_counter = 0
    def lift(lambda_expr: LambdaExpr, parent_name: str) -> tuple[Call, HIRFunc]:
        nonlocal lambda_counter
        lambda_counter += 1
        lam_name = f"_lambda_{parent_name}_{lambda_counter}"
        
        # 1. 分析 Lambda body：哪些 captures 被修改？
        modified_caps = find_modified_captures(lambda_expr.body, lambda_expr.captures)
        
        # 2. 构造新 HIRFunc
        #    参数 = captures + lambda 原参数
        new_params = [
            HIRParam(c.name, lookup_type(c.name, func), is_const_ref=not c.by_ref)
            for c in lambda_expr.captures
        ] + lambda_expr.params
        
        # 3. 返回类型：若无 modified capture → 原 Lambda 返回
        #           若有 → (原返回, 修改后 capture1, ..., 修改后 capN) tuple
        if modified_caps:
            ret_ty = TupleType([lambda_expr.body_ret_ty] + [lookup_type(c) for c in modified_caps])
            # 改写 body 末尾：追加修改后的 capture 值
            new_body = append_modified_caps_to_return(lambda_expr.body, modified_caps)
        else:
            ret_ty = lambda_expr.body_ret_ty
            new_body = lambda_expr.body
        
        aux = HIRFunc(
            base_name=lam_name,
            instance_suffix="",
            mangled_name="",
            qual_type="<lifted lambda>",
            params=new_params,
            ret_ty=ret_ty,
            body=new_body,
        )
        
        # 4. 构造调用点替换
        call_args = [Var(c.name) for c in lambda_expr.captures] + \
                    [Var(p.name) for p in lambda_expr.params]
        call = Call(lam_name, call_args)
        
        return call, aux
    
    def walk_replace(node, parent_name):
        if isinstance(node, LambdaExpr):
            call, aux = lift(node, parent_name)
            func.aux_lambdas.append(aux)
            return call
        # 递归处理 children
        for field in node.__dict__:
            ...  # 递归替换
        return node
    
    func.body = [walk_replace(s, func.base_name) for s in func.body]
    return func
```

### §5.4 `modified_captures` 检测

扫描 lambda body，对每条语句算 `target_root(expr) -> str | None`，若目标根是某个 capture，就加入 `modified_caps`。

**`target_root` 递归规则**：
- `Var(x)` → `x`
- `ArrayAccess(arr, _)` → `target_root(arr)`
- `FieldAccess(obj, _)` → `target_root(obj)`
- `Cast(inner, _)` → `target_root(inner)`
- `Call(UnresolvedOp("operator[]"), [arr, _])` → `target_root(arr)`
  — Pass 1 对 `a[i]` 可能生成 `Call` 而非 `ArrayAccess`（例如嵌套 `M[i][k]` 或受 `operator-=` / `swap` 包裹时），必须递归到实参 0 上。
- 其他情况 → `None`

**判定为修改的 pattern**：
| 源码示例 | HIR 节点 | 被修改的根 |
|---|---|---|
| `x = e` | `AssignStmt(target, e)` | `target_root(target)` |
| `x += e` / `x *= e` / `x <<= e` ... | `CompoundAssignStmt(target, op, e)` | `target_root(target)` |
| `++x` / `--x` / `x++` / `x--` | `ExprStmt(Call(UnresolvedOp("operator++"), [x]))` 等 | `target_root(args[0])` |
| `a += b` / `a -= b` ...（作为 ExprStmt） | `ExprStmt(Call(UnresolvedOp("operator+="), [lhs, rhs]))` | `target_root(args[0])` |
| `swap(a, b)` / `std::swap(a, b)` | `ExprStmt(Call("swap", [a, b]))` — callee 是字符串 | `target_root(args[0])` **和** `target_root(args[1])` 都算 |

**遗漏案例（已知 TODO）**：
- mutator method call（`vec.push_back(x)`、`vec.clear()` 等）未识别为修改 receiver。当前 Pass 3 的 15 个 lifted lambda 中没有此模式；若后续引入需要扩展。
- 非 `swap` 的自由函数若修改 ref 实参（如自定义 out-param helper）同样未识别。保持关闭；Pass 4 的 SSA/alias 分析阶段会兜底。

**保守策略已弃用**：早期 `translation-fix-plan.md` M4 建议 `[&]` 全捕获全部视作修改，但实测会产生大量 spurious tuple 返回值，现改为上述**精确检测**；smoke `test_smoke_all_65` 覆盖所有 15 宿主 26 lambda，任何漏检会导致后续 Pass 失败。

**输出方式**：Pass 3 实现不新建字段，`modified_caps` 通过 aux `HIRFunc.qual_type` 字符串追加 `| modified_captures=[...]` 记录；同时把对应 capture param 标为 `is_ref=True, is_const_ref=False`（其余 capture 标 `is_const_ref=True`）。

### §5.5 调用方改写（modified capture 写回）

**设计目标**：调用点最终要变成：

```
old:  LAMBDA(...)  (inline LambdaExpr)
target:
  若无 modified → Call("_lambda_N", [cap1, cap2, x, y])
  若有 modified caps [cap1] →
    LetStmt(_lam_out, ret_ty, Call("_lambda_N", [cap1, cap2, x, y]))
    + AssignStmt(Var("cap1"), FieldAccess(Var("_lam_out"), "1"))
    + 若有原返回值则 AssignStmt(Var("result"), FieldAccess(Var("_lam_out"), "0"))
```

**Pass 3 当前实现**：仅做 `LambdaExpr → Var(lam_name)` 简单替换，保留宿主 body 里原 `Call(...)`/`sort(...,lam)` 结构。换言之：
- `sort(it1, it2, LAMBDA(...))` → `sort(it1, it2, Var("_lambda_N"))`
- `auto f = LAMBDA(...); f(x)` → `let f := Var("_lambda_N"); Var("f")(x)`

aux `HIRFunc.ret_ty` 不改写为 tuple；`modified_captures` 只保存在 `qual_type` 字符串里，供下游 Pass 读取。

**剩余工作（后续 Pass 认领）**：
- **tuple 返回类型改写**（ret_ty 追加 modified caps）—— 留给 Pass 5 `operator_resolve` 或 Pass 6 `ssa_build`，需先清楚调用点的表达式语境。
- **capture 写回展开**（`LetStmt + AssignStmt` 注入）—— 同上。
- **body 里 return 语句改写**（`return e` → `return (e, cap1, ...)`）—— 同上。
- 当前 smoke 可通过，因为 HIR₂→HIR₃ 不检查这三项，下游依赖 `qual_type` 的 `modified_captures=[...]` 发起改写。

### §5.6 归属（aux_lambdas 位置）

提升后的 lambda 作为 `HIRFunc.aux_lambdas` 字段。Pass 8 codegen 时**先输出 aux_lambdas 再输出主函数**（Lean 要求被调函数先定义）。

### §5.7 语义保持

见 `cpp-subset-semantics.md` §5.3 引理 L5.1 + 5.4（`[&]` capture 写回正确性）。

### §5.8 实现规模

~300 行 Python。

---

## §6 Pass 4: `iter_recognize` (HIR₂ → HIR₃)

### §6.1 职责

识别迭代器模式，把裸迭代器操作（`IteratorExpr`、classic iterator loop、compact-erase 双指针模式）转为高阶 Call 或显式索引循环。依据 `cpp-subset-semantics.md` §6 + `iterators.md`。

### §6.2 识别的 4 种模式

| 模式 | 频次 | 转换目标 |
|---|---|---|
| `RangeForStmt` | 92 | 保留（已是高阶形式）+ 结构化绑定 desugar |
| Structured-binding range-for | 24 | `RangeForStmt.decomposition` 字段处理 |
| **Compact-erase 双指针** | **4** | **`AssignStmt(Var(v), Call("Array.filter", [Var(v), pred_lambda]))`** |
| Classic iterator loop | 1 | 展开为显式 index-based `ForStmt` |

### §6.3 Compact-erase 模式识别

**输入模式**（HIR₂ 形式，抽象化后）：

```
LetStmt(it_0, _, IteratorExpr("begin", Var(v)))
LetStmt(out_0, _, Var(it_0))
WhileStmt(
    BinOp("!=", Var(it_0), IteratorExpr("end", Var(v))),
    body=[
        IfStmt(pred_expr_with_deref(Var(it_0)),
               then=[ ... AssignStmt(FieldAccess(Var(out_0), "deref"), IteratorExpr("deref", Var(it_0)))
                     + IteratorExpr("increment", Var(out_0)) ... ],
               else=[]),
        AssignStmt(Var(it_0), IteratorExpr("increment", Var(it_0)))
    ]
)
Call("erase", [Var(v), Var(out_0), IteratorExpr("end", Var(v))])
```

**匹配算法**：
```python
def match_compact_erase(stmts: list[StmtIR]) -> Optional[tuple[Var, LambdaExpr]]:
    """尝试识别 compact-erase 模式。
    成功返回 (容器 Var, 过滤谓词 Lambda)。"""
    # 1. 前 2 条 LetStmt 是 it = v.begin(); out = it;
    # 2. 主 WhileStmt with cond it != v.end()
    # 3. 末尾 erase(v, out, v.end())
    # 4. 构造 lambda_pred from IfStmt.cond
    ...
```

**转换**：

```
match [
    LetStmt(it, _, IteratorExpr("begin", Var(v))),
    LetStmt(out, _, Var(it)),
    WhileStmt(...cond_expr = *it...),
    ExprStmt(Call("erase", ...))
]
→
AssignStmt(Var(v), Call("Array.filter", [Var(v), LambdaExpr(pred_body)]))
```

### §6.4 Structured-binding range-for

在 HIR₀ → HIR₁ 时，`parse` Pass 已把 `DecompositionDecl` 放到 `RangeForStmt.decomposition`。Pass 4 将其 desugar 为 body 开头的 `LetStmt`：

```
RangeForStmt(
    var=_x,
    ty=PairType(K, V),
    container=map,
    body=...body_using_k_and_v...,
    decomposition=[Var("k"), Var("v")]
)
→
RangeForStmt(
    var=_x,
    ty=PairType(K, V),
    container=map,
    body=[
        LetStmt(Var("k"), K, FieldAccess(Var("_x"), "fst")),
        LetStmt(Var("v"), V, FieldAccess(Var("_x"), "snd")),
    ] + body,
    decomposition=None
)
```

### §6.5 Classic iterator loop

Pattern: `for (auto it = v.begin(); it != v.end(); ++it) { body }`

**转换**：展开为显式 `ForStmt`

```
ForStmt(
    init=[LetStmt(Var("__idx"), BaseType.NAT, Lit(0))],
    cond=BinOp("<", Var("__idx"), Call("Array.size", [Var(v)])),
    step=[CompoundAssignStmt(Var("__idx"), "+", Lit(1))],
    body=[
        LetStmt(Var("it"), elem_ty, ArrayAccess(Var(v), Var("__idx"))),
        ... body with *it replaced by Var("it"), ++it removed ...
    ]
)
```

然后 Pass 7 `loop_lower` 再处理 ForStmt。

### §6.6 并行双迭代器（1 处）

`__upoly_divmod_mod` 的 zip-walk 不识别，保留原始 HIR（由 `UnresolvedOp("operator++")` 等表示）。最终由 `operator_resolve` Pass 解析为普通 Call。

### §6.7 语义保持

见 `cpp-subset-semantics.md` §6.3 引理 L6.1。

### §6.8 实现规模

~350 行 Python。

---

## §7 Pass 5: `operator_resolve` (HIR₃ → HIR₄)

### §7.1 职责

把所有 `UnresolvedOp`、`CompoundAssignStmt`、`Cast` 查表解析为具体 Lean 函数调用。依据 `type-system.md` §8（CAST_TABLE）+ `class_map.py` CLASS_MAP/FUNC_MAP。

### §7.2 输入/输出

- **输入**：`func : HIRFunc`（HIR₃）
- **输出**：同 func，所有 `UnresolvedOp` 替换为具体 `Call`、`Cast` 替换为具体 Call 或 noop、`CompoundAssignStmt` 展开为 `AssignStmt + BinOp`

### §7.3 查表解析

从 v1 复用的 `class_map.py`：

```python
# CLASS_MAP: 类型 → 方法表
CLASS_MAP = {
    "SparsePolyZp": {
        "front!": ("direct", "SparsePolyZp.front!"),
        "back!": ("direct", "SparsePolyZp.back!"),
        "push_back": ("mutate", "Array.push"),  # mutate = 需要 AssignStmt
        "size": ("direct", "Array.size"),
        # ...
    },
    "Zp": {
        "prime": ("field", "prime"),     # FieldAccess
        "number": ("field", "val"),
        "inv": ("direct", "Zp.inv"),
        # ...
    },
    # ...
}

# FUNC_MAP: 全局函数 → Lean 函数
FUNC_MAP = {
    "std::sort": ("direct", "Array.qsort"),
    "std::max": ("direct", "max"),
    "std::iota": ("special", "iota_handler"),
    "poly_convert": ("direct", "poly_convert"),
    # ...
}

# CAST_TABLE: (source, target) → Lean 表达式模板
CAST_TABLE = {
    ("int", "size_type"): "{x}.toNat.toUInt64",
    ("int", "int64_t"): "{x}.toInt64",
    # ... from type-system.md §8.3
}
```

### §7.4 转换规则

**`UnresolvedOp`**：
```python
def resolve(op_call: Call) -> Call | ExprIR:
    op = op_call.callee  # UnresolvedOp
    receiver_ty = op_call.args[0].ty if op_call.args else None
    
    # 运算符类型
    if op.op_name.startswith("operator"):
        # 查 CLASS_MAP
        method_name = op.op_name[len("operator"):]  # "[]", "+", "!=", ...
        if isinstance(receiver_ty, NamedType):
            cls = CLASS_MAP.get(receiver_ty.name, {})
            if method_name in cls:
                disposition, target = cls[method_name]
                if disposition == "direct":
                    return Call(target, op_call.args)
                elif disposition == "mutate":
                    # 生成 AssignStmt 而非 Call
                    ...
    # 查 FUNC_MAP
    if op.op_name in FUNC_MAP:
        disposition, target = FUNC_MAP[op.op_name]
        return Call(target, op_call.args)
    
    raise TranslationError(f"operator_resolve: can't resolve {op.op_name} on {receiver_ty}")
```

**`Cast`**：

```python
def resolve_cast(cast: Cast) -> ExprIR:
    if cast.cast_kind in ("NoOp", "LValueToRValue", "FunctionToPointerDecay",
                          "BuiltinFnToFnPtr", "ArrayToPointerDecay", "ToVoid"):
        return cast.expr  # 消除 cast
    
    if cast.cast_kind == "IntegralCast":
        key = (str(cast.source_ty), str(cast.target_ty))
        template = CAST_TABLE.get(key)
        if template:
            return build_expr_from_template(template, cast.expr)
        raise TranslationError(f"operator_resolve: cast {key} not in CAST_TABLE")
    
    if cast.cast_kind == "ConstructorConversion":
        # 查 FUNC_MAP 或直接用构造器
        ...
    
    if cast.cast_kind == "IntegralToFloating":
        return Call("Int.toFloat", [cast.expr])
    if cast.cast_kind == "FloatingToIntegral":
        return Call("Float.toInt", [cast.expr])
    
    raise TranslationError(f"operator_resolve: unknown cast kind {cast.cast_kind}")
```

**`CompoundAssignStmt`**：

```python
def expand_compound(stmt: CompoundAssignStmt) -> AssignStmt:
    return AssignStmt(
        target=stmt.target,
        value=BinOp(stmt.op, stmt.target, stmt.value)
    )
```

### §7.5 调用方 output param 解构（Pass 2 遗留）

对每个 `Call(callee, args)`：

```python
def rewrite_output_param_call(call: Call, call_stmt: StmtIR) -> list[StmtIR]:
    fn_name = call.callee
    output_indices = OUTPUT_PARAMS.get(fn_name, [])
    if not output_indices:
        return [call_stmt]  # 无 output，不改
    # 生成 (tmp = call; 解构到各 output var)
    tmp = Var("_out_tmp")
    ...
```

### §7.6 `require` 生成（UB 点）

`operator_resolve` 遇到**会触发 UB 的运算符**时，**在同一语句前插入 `RequireStmt`**：

| 运算符 | 生成的 require |
|---|---|
| `BinOp("/", _, b)` for integer | `RequireStmt(BinOp("≠", b, Lit(0)), ..., "div_by_zero")` |
| `BinOp("%", _, b)` | 同上 |
| `Call("Array.get!", [arr, i])` | `RequireStmt(BinOp("<", i, Call("Array.size", [arr])), ..., "array_oob")` |
| `Call("Array.set!", [arr, i, _])` | 同上 |
| `Call("MvPoly.front!", [v])` | `RequireStmt(UnaryOp("!", Call("Array.isEmpty", ...)), ..., "empty_container")` |
| `Call("MvPoly.back!", [v])` | 同上 |

Signed overflow 单独扫描 `UB-6` 站点，批量生成 `Int.noOverflow` require。

### §7.7 实现规模

~250 行 Python（+ CLASS_MAP/FUNC_MAP/CAST_TABLE 的数据表，复用自 v1）。

---

## §8 手动走查（3 个复杂函数）

对 Week 1 发现的 3 个最复杂函数，手动走过 Pass 1-5，验证设计可行。

### §8.1 `__lll_reduce`（5 lambda、132 UB 站点，最复杂）

**HIR₀**（Pass 1 后，节选）：
```
HIRFunc("__lll_reduce", instance_suffix="",
  params=[
    HIRParam("M", NamedType("LLLMatrix"), is_ref=True),     # 输出！
    HIRParam("U", NamedType("LLLMatrix"), is_ref=True),     # 输出！
    HIRParam("n", BaseType.INT32, is_const_ref=False),      # 值
    ...
  ],
  ret_ty=BaseType.UNIT,  # void
  body=[
    LetStmt(Var("mu"), ArrayType(ArrayType(NamedType("QQ"))), Call(UnresolvedOp("construct"), ...)),
    LetStmt(Var("B_gs"), ArrayType(NamedType("QQ")), ...),
    LetStmt(Var("row_sub"), NamedType("Lambda"), LambdaExpr(
      captures=[Capture("M", by_ref=True), Capture("U", by_ref=True), ...],
      params=[HIRParam("i", INT32), HIRParam("j", INT32), HIRParam("c", NamedType("ZZ"))],
      body=[...ForStmt(... AssignStmt(ArrayAccess(ArrayAccess(Var("M"), Var("i")), Var("k")), ...))...],
    )),
    LetStmt(Var("row_swap"), NamedType("Lambda"), LambdaExpr(...)),
    # ... 3 more lambdas
    WhileStmt(cond=..., body=[...IfStmt(...)+ Call(Var("row_sub"), [k, k-1, mu_kk1])...]),
  ]
)
```

**Pass 2 ref_elim**：`M` 和 `U` 从参数移到返回 tuple。签名变为：
```
params=[n, ...non-ref params..., M_init, U_init]
ret_ty=TupleType([BaseType.UNIT, NamedType("LLLMatrix"), NamedType("LLLMatrix")])
```
所有 return 点（只有末尾一个隐式 return）包装为 `ReturnStmt(TupleExpr([Unit, M, U]))`。

**Pass 3 lambda_lift**：5 个 Lambda 提升：
- `_lambda___lll_reduce_1_ir`（row_sub，修改 M、U）
- `_lambda___lll_reduce_2_ir`（row_swap，修改 M、U）
- `_lambda___lll_reduce_3_ir`（dot）
- `_lambda___lll_reduce_4_ir`（round_qq）
- `_lambda___lll_reduce_5_ir`（row_sub 第二个 body）

每个带 `[&]` capture 的 lambda 返回 tuple 含修改后 M/U。

**Pass 4 iter_recognize**：内部有若干 range-for，全转为高阶 Array.foreach；无 compact-erase 模式。

**Pass 5 operator_resolve**：
- `M[i][j]` → `Array.get! (Array.get! M i) j`，前插 2 个 `RequireStmt`（外层和内层越界）
- `ZZ(0)` → `(0 : ZZ)`
- `QQ.abs` → `Rat.natAbs`
- Lambda 调用 `row_sub(k, k-1, mu_kk1)` → `_lambda___lll_reduce_1_ir M U k (k-1) mu_kk1`

**UB require 累积**：约 130 个，主要是 OOB（108）+ signed overflow（21）。

### §8.2 `__mtshl_step_j`（2 lambda，25 UB）

**特点**：2 个 lambda（`lc_correct`、`product_F`）都 `[&]` 捕获并修改 `F` 数组。

**Pass 2**：F 是 `non-const ref` 参数 → 转 tuple 返回。

**Pass 3**：2 lambda 提升，各带 modified captures 返回。

**Pass 5**：Call `_lambda___mtshl_step_j_1_ir F lc_tau_i` 的返回值解构到 F。

### §8.3 `__zassenhaus_recombine`（2 lambda，31 UB）

**特点**：包含 `next_combination` lambda 做子集枚举，核心算法 140+ 行。

**关键挑战**：do-while 循环（`do { ... } while (next_combination)`）。Pass 7 `loop_lower` 需要处理，但 HIR 层保留 `DoWhileStmt`。

**Pass 3**：next_combination lambda 提升为独立 def，body 含 ArrayAccess 和赋值操作。

**Pass 5**：`std::iota(T.begin(), T.end(), 0)` → `Array.range r`（r 是 T 的初始 size）。

---

## §9 实现计划（Stage 2 映射）

| Pass | 文件 | 预估代码量 | 复用来源 |
|---|---|---|---|
| `parse` | `pass_parse.py` | 400 行 | 30% 从 v1 `clang_ast.py` 复用 |
| `ref_elim` | `pass_ref_elim.py` | 200 行 | 新写（v1 手工 OUTPUT_PARAMS 废弃）|
| `lambda_lift` | `pass_lambda_lift.py` | 300 行 | 新写 |
| `iter_recognize` | `pass_iter_recognize.py` | 350 行 | 少量 v1 启发式参考 |
| `operator_resolve` | `pass_operator_resolve.py` | 250 行 | 查表逻辑新写；CLASS_MAP 复用 |
| **HIR 总计** | | **~1500 行** | |

加上数据结构定义：
- `ir_types.py`（HIR 节点 dataclass）：150 行
- `class_map.py`（查表）：500 行（v1 复用）
- `clang_hybrid.py`（frontend）：200 行（v1 复用）

**HIR 层总体：~2350 行**。加上 MIR 层（Week 5 设计）+ codegen，总体估 ~4250 行（与 v1 5956 行相比减少 29%）。

---

## §10 与 Week 5（MIR）衔接

HIR₄ → MIR 需要新的 dataclass 族：
- `MIRStmt`（含 `PhiStmt`、`TailCallStmt`、`ValueStmt`、**无** `WhileStmt`/`ForStmt`/`BreakStmt`/`ContinueStmt`/`ReturnStmt`）
- `MIRFunc`（含 `cfg : CFG`、`ssa_vars : dict[Var, Version]`）

**Pass 6 `ssa_build`**：
- 输入：HIR₄
- 输出：MIR₀（SSA + phi + CFG）
- 依据：`cpp-subset-semantics.md` §3 L3.1（Appel 1998）

**Pass 7 `loop_lower`**：
- 输入：MIR₀
- 输出：MIR₁（循环已提取为 partial def 尾递归 + break/continue/return 下降）
- 依据：`cpp-subset-semantics.md` §2 L2.2/L2.3（Winskel + Aeneas）

**Pass 8 `codegen`**：
- 输入：MIR₁
- 输出：Lean 源代码
- 依据：`type-system.md` §1-8 的具体映射

---

## §11 Week 4 验收

- [x] §1 HIR 数据结构完整定义（TypeIR、ExprIR、StmtIR、HIRParam、HIRFunc、HIRProgram）
- [x] §2 HIR 不变量阶梯（HIR₀-HIR₄ + runtime assert）
- [x] §3 Pass 1 `parse` 完整规格（AST kind → HIR 节点表）
- [x] §4 Pass 2 `ref_elim` 算法 + 调用方改写策略
- [x] §5 Pass 3 `lambda_lift` 算法 + modified capture 处理
- [x] §6 Pass 4 `iter_recognize` 4 种模式处理
- [x] §7 Pass 5 `operator_resolve` 查表解析 + UB require 生成
- [x] §8 手动走查 3 个复杂函数（`__lll_reduce` / `__mtshl_step_j` / `__zassenhaus_recombine`）
- [x] §9 实现计划（~1500 行 HIR Pass + 850 行数据结构/frontend）
- [x] §10 与 Week 5 (MIR) 的衔接接口

**本文件是 Stage 1 Week 4 的硬性产出**，为 Stage 2 的 HIR 层实现提供 1:1 可编程的蓝图。
