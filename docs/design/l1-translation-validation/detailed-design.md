# C++ → Lean IR 翻译器细化设计

> 基于 `blueprint.md` 架构阶段产出。不改变模块划分、接口规约和核心流程。
>
> **核心原则**：翻译器不做类型推断——所有类型信息从 Clang AST 传播（blueprint §2a P1）。

---

## 1. 数据结构定义

所有模块共享的 IR 类型定义，放在 `ir_types.py`：

```python
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional

# ============================================================
# 类型系统
# ============================================================

class BaseType(Enum):
    UINT64 = "UInt64"
    INT64 = "Int"        # Lean Int (arbitrary precision, but maps to int64)
    UINT32 = "UInt32"
    UINT128 = "UInt128"
    BOOL = "Bool"
    VOID = "Unit"

@dataclass
class ArrayType:
    elem: 'TypeIR'

@dataclass
class PairType:
    fst: 'TypeIR'
    snd: 'TypeIR'

@dataclass
class StructType:
    name: str                          # e.g. "SparsePolyZp"
    fields: list[tuple[str, 'TypeIR']] # [(field_name, type)]

@dataclass
class ExceptType:
    inner: 'TypeIR'                    # Except Error inner

TypeIR = BaseType | ArrayType | PairType | StructType | ExceptType | str
# str 用于未识别类型（保留原始 C++ 类型名）

# ============================================================
# 表达式
# ============================================================

@dataclass
class Var:
    name: str
    version: int = 0                   # SSA 版本号

@dataclass
class Lit:
    value: int
    typ: BaseType = BaseType.UINT64

@dataclass
class BinOp:
    op: str                            # "+", "-", "*", "/", "%", "<<", ">>",
                                       # "<", ">", "<=", ">=", "==", "!="
    lhs: 'ExprIR'
    rhs: 'ExprIR'

@dataclass
class UnaryOp:
    op: str                            # "!", "-" (negation)
    operand: 'ExprIR'

@dataclass
class CondExpr:
    cond: 'ExprIR'
    then_e: 'ExprIR'
    else_e: 'ExprIR'

@dataclass
class Call:
    func: str
    args: list['ExprIR']

@dataclass
class ArrayAccess:
    arr: 'ExprIR'
    idx: 'ExprIR'

@dataclass
class FieldAccess:
    obj: 'ExprIR'
    field_name: str

@dataclass
class Cast:
    expr: 'ExprIR'
    target_type: TypeIR

@dataclass
class ArrayPush:
    arr: 'ExprIR'
    elem: 'ExprIR'

@dataclass
class UnknownExpr:
    kind: str
    children: list['ExprIR'] = field(default_factory=list)
    raw: str = ""

ExprIR = Var | Lit | BinOp | UnaryOp | CondExpr | Call | ArrayAccess | \
         FieldAccess | Cast | ArrayPush | UnknownExpr

# ============================================================
# 语句
# ============================================================

@dataclass
class LetStmt:
    """let x_{n} : T := e"""
    var: Var
    typ: TypeIR
    value: ExprIR

@dataclass
class AssignStmt:
    """x = e（SSA 前，SSA 后变为 LetStmt）"""
    target: Var
    value: ExprIR

@dataclass
class IfStmt:
    cond: ExprIR
    then_body: list['StmtIR']
    else_body: list['StmtIR'] = field(default_factory=list)

@dataclass
class ReturnStmt:
    value: Optional[ExprIR] = None

@dataclass
class Require:
    """assert(cond) 或隐式 UB → require h : cond"""
    cond: ExprIR
    name: str                          # 证明参数名，如 "hp"
    source: str                        # "assert" | "div_by_zero" | "array_oob" | ...

@dataclass
class Throw:
    """throw error → Except.error"""
    error_tag: str
    message: str = ""

@dataclass
class TailRec:
    """循环 → 尾递归"""
    func_name: str                     # 生成的递归函数名
    params: list[tuple[Var, TypeIR]]   # 循环变量 + 累积状态
    exit_cond: ExprIR                  # 循环条件取反 = 退出条件
    break_cond: Optional[ExprIR]       # break 条件（无则 None）
    body: list['StmtIR']              # 循环体
    step: list['StmtIR']              # 步进（i++ 等）

@dataclass
class ExprStmt:
    """表达式语句（如纯函数调用）"""
    expr: ExprIR

@dataclass
class UnknownStmt:
    kind: str
    children: list['StmtIR'] = field(default_factory=list)

StmtIR = LetStmt | AssignStmt | IfStmt | ReturnStmt | Require | Throw | \
         TailRec | ExprStmt | UnknownStmt

# ============================================================
# 函数
# ============================================================

@dataclass
class ParamIR:
    name: str
    typ: TypeIR
    is_output: bool = False            # T& out 参数
    is_const_ref: bool = False         # const T& 参数

@dataclass
class FuncIR:
    """Clang AST 解析后的函数表示"""
    name: str
    params: list[ParamIR]
    ret_type: TypeIR
    body: list[StmtIR]
    source_file: str = ""
    source_line: int = 0

@dataclass
class SSAFunc:
    """SSA 变换后的函数表示"""
    name: str
    params: list[ParamIR]              # 输出参数已转为返回值
    ret_type: TypeIR                   # 可能是 ExceptType
    requires: list[Require]            # 前置条件（assert + 隐式 UB）
    body: list[StmtIR]                # 全部 LetStmt（单赋值）
    has_throw: bool = False            # 是否包含 throw

# ============================================================
# UB 证明目标
# ============================================================

class UBType(Enum):
    DIV_BY_ZERO = auto()
    ARRAY_OOB = auto()
    SIGNED_OVERFLOW = auto()
    SHIFT_OOB = auto()
    ASSERT_FAIL = auto()

@dataclass
class UBObligation:
    func_name: str
    source_line: int
    ub_type: UBType
    lean_prop: str                     # Lean Prop 表达式
    context: list[str]                 # 相关 let 绑定
```

## 2. 模块细化

### 2.1 clang_ast — Clang AST JSON 解析器

```
文件：clang_ast.py
依赖：ir_types.py, json (标准库)
```

| 函数 | 签名 | 功能 | 调用 |
|------|------|------|------|
| `parse_translation_unit` | `(json_str: str) → list[FuncIR]` | 入口：解析整个翻译单元 | `parse_func_decl` |
| `parse_func_decl` | `(node: dict) → FuncIR` | 解析单个函数声明 | `parse_param`, `parse_stmt` |
| `parse_param` | `(node: dict) → ParamIR` | 解析参数（含 const ref/out ref 判断） | `parse_type` |
| `parse_type` | `(qual_type: str) → TypeIR` | C++ 类型字符串 → TypeIR | — |
| `parse_stmt` | `(node: dict) → StmtIR` | 分发到具体语句解析 | `parse_expr`, `parse_if`, `parse_for`, ... |
| `parse_expr` | `(node: dict) → ExprIR` | 分发到具体表达式解析 | — |
| `parse_for` | `(node: dict) → list[StmtIR]` | for 循环 → 语句列表（含 AssignStmt） | `parse_stmt`, `parse_expr` |
| `parse_while` | `(node: dict) → list[StmtIR]` | while 循环 → 语句列表 | `parse_stmt`, `parse_expr` |

```python
def parse_translation_unit(json_str: str) -> list[FuncIR]:
    """解析 Clang AST JSON，提取所有有函数体的 FunctionDecl。

    跳过：无函数体的声明、隐式声明、系统头文件中的声明。
    """

def parse_type(qual_type: str) -> TypeIR:
    """C++ qualType 字符串 → TypeIR。

    映射规则：
    - "uint64_t" / "unsigned long" → BaseType.UINT64
    - "int" / "int32_t" → BaseType.INT64
    - "unsigned __int128" → BaseType.UINT128
    - "bool" / "_Bool" → BaseType.BOOL
    - "void" → BaseType.VOID
    - "std::vector<T>" / "vector<T>" → ArrayType(parse_type(T))
    - "std::pair<A,B>" → PairType(parse_type(A), parse_type(B))
    - 其他 → str（保留原始类型名，后续手动处理）
    """

def parse_param(node: dict) -> ParamIR:
    """解析参数声明。

    判断 is_output：qualType 含 "&" 且不含 "const" → is_output=True
    判断 is_const_ref：qualType 含 "const" 和 "&" → is_const_ref=True
    """

def parse_stmt(node: dict) -> StmtIR:
    """按 Clang AST kind 分发：

    DeclStmt → LetStmt（提取变量名+类型+初始值）
    BinaryOperator(=) → AssignStmt
    IfStmt → IfStmt
    ForStmt → parse_for()
    WhileStmt → parse_while()
    ReturnStmt → ReturnStmt
    CallExpr(assert) → Require(source="assert")
    CXXThrowExpr → Throw
    CompoundStmt → 递归 parse 每个子语句
    其他 → UnknownStmt(kind, children)
    """

def parse_expr(node: dict) -> ExprIR:
    """按 Clang AST kind 分发：

    IntegerLiteral → Lit
    DeclRefExpr → Var
    BinaryOperator → BinOp
    UnaryOperator → UnaryOp
    ConditionalOperator → CondExpr
    CallExpr → Call
    ArraySubscriptExpr → ArrayAccess
    MemberExpr → FieldAccess
    CStyleCastExpr / ImplicitCastExpr → Cast 或直接递归（无损转换时）
    CXXMemberCallExpr(push_back) → ArrayPush
    其他 → UnknownExpr(kind)
    """
```

**复用**：原型 `cpp2lean.py` 的 `translate_expr` 方法提供了表达式解析框架，可直接扩展。

**错误处理**：UnknownStmt/UnknownExpr 保留原始 AST kind 和子节点，不丢失信息。解析后统计未知节点，生成警告报告。

### 2.2 ssa_transform — SSA 变换引擎

```
文件：ssa_transform.py
依赖：ir_types.py
```

| 函数 | 签名 | 功能 | 调用 |
|------|------|------|------|
| `transform_func` | `(func: FuncIR) → SSAFunc` | 入口：整个函数的 SSA 变换 | `transform_params`, `transform_body` |
| `transform_params` | `(params: list[ParamIR]) → (new_params, out_vars)` | 分离输出参数 → 返回值 | — |
| `transform_body` | `(stmts: list[StmtIR], env: VarEnv) → list[StmtIR]` | 语句列表 SSA 变换 | `rename_vars`, `transform_loop` |
| `rename_vars` | `(expr: ExprIR, env: VarEnv) → ExprIR` | 表达式中变量引用 → 最新版本 | — |
| `new_version` | `(env: VarEnv, name: str) → Var` | 创建变量新版本 x_{n+1} | — |
| `transform_loop` | `(stmt: ForStmt/WhileStmt, env: VarEnv) → TailRec` | for/while → TailRec 节点 | `identify_loop_vars`, `transform_body` |
| `identify_loop_vars` | `(stmts: list[StmtIR]) → list[Var]` | 识别循环中被修改的变量（= 尾递归参数） | — |
| `transform_if` | `(stmt: IfStmt, env: VarEnv) → (list[StmtIR], VarEnv)` | if/else SSA：两分支可能产生不同版本 → phi | `transform_body` |
| `collect_requires` | `(stmts: list[StmtIR]) → list[Require]` | 提取所有 Require 节点到函数签名 | — |
| `detect_throws` | `(stmts: list[StmtIR]) → bool` | 检测是否包含 Throw 节点 | — |

```python
class VarEnv:
    """变量版本环境。跟踪每个变量的当前版本号。"""

    def __init__(self):
        self.versions: dict[str, int] = {}   # name → current version

    def current(self, name: str) -> Var:
        """获取变量当前版本。"""
        return Var(name, self.versions.get(name, 0))

    def bump(self, name: str) -> Var:
        """创建变量新版本。x_n → x_{n+1}。"""
        v = self.versions.get(name, 0) + 1
        self.versions[name] = v
        return Var(name, v)

    def fork(self) -> 'VarEnv':
        """创建分支副本（用于 if/else 两分支）。"""
        new = VarEnv()
        new.versions = dict(self.versions)
        return new

    def merge(self, other: 'VarEnv', changed: set[str]) -> dict[str, Var]:
        """合并两分支（if/else 后）。
        对版本不同的变量创建 phi 节点（在 Lean 中 = if 表达式选择）。
        """

def transform_loop(stmt, env: VarEnv) -> TailRec:
    """for/while → TailRec。

    步骤：
    1. identify_loop_vars：扫描循环体，找所有被赋值的变量
    2. 这些变量 + 循环计数器 = TailRec 的参数
    3. 循环条件取反 = exit_cond
    4. 扫描循环体中的 break → break_cond
    5. 循环体 SSA 变换 → TailRec.body
    6. 步进（i++ 等） → TailRec.step
    """

def transform_if(stmt: IfStmt, env: VarEnv) -> tuple[list[StmtIR], VarEnv]:
    """if/else 的 SSA 变换。

    问题：then 分支可能修改 x → x_2，else 分支可能修改 x → x_3。
    if 后需要 "phi"：x_4 = if cond then x_2 else x_3。

    实现：
    1. fork env → env_then, env_else
    2. transform_body(then, env_then), transform_body(else, env_else)
    3. merge：对每个版本不同的变量，生成 LetStmt(x_4 = if cond then x_2 else x_3)
    4. 返回合并后的 env
    """
```

**关键设计决策**：

1. **if/else phi 节点**：C++ 中 if/else 两分支可能对同一变量赋不同值。SSA 需要 phi 节点。在 Lean 中，phi = `let x_4 := if cond then x_2 else x_3`。`transform_if` 通过 fork/merge 实现。

2. **循环变量识别**：`identify_loop_vars` 扫描循环体中所有 AssignStmt 的目标变量。这些变量成为尾递归函数的参数，循环体中的赋值变为递归调用时的参数传递。

3. **输出参数**：`transform_params` 将 `T& out` 参数从参数列表移除，加入返回类型（多返回值时用 tuple）。函数体中对 out 的赋值在 SSA 中自然变为 let 链，最终版本作为返回值。

### 2.3 lean_codegen — Lean 4 代码生成器

```
文件：lean_codegen.py
依赖：ir_types.py
```

| 函数 | 签名 | 功能 | 调用 |
|------|------|------|------|
| `generate_file` | `(funcs: list[SSAFunc]) → str` | 入口：生成完整 .lean 文件 | `generate_func` |
| `generate_func` | `(func: SSAFunc) → str` | 生成单个 partial def | `gen_signature`, `gen_body` |
| `gen_signature` | `(func: SSAFunc) → str` | 生成函数签名（含 require 参数） | `gen_type` |
| `gen_type` | `(t: TypeIR) → str` | TypeIR → Lean 类型字符串 | — |
| `gen_body` | `(stmts: list[StmtIR], indent: int) → str` | 语句列表 → Lean let 链 | `gen_stmt`, `gen_expr` |
| `gen_stmt` | `(stmt: StmtIR, indent: int) → str` | 单条语句 → Lean 代码 | `gen_expr` |
| `gen_expr` | `(expr: ExprIR) → str` | 表达式 → Lean 表达式字符串 | — |
| `gen_tailrec` | `(tr: TailRec, indent: int) → str` | TailRec → partial def 尾递归函数 | `gen_body` |
| `gen_require` | `(req: Require) → str` | Require → 函数参数 `(h : prop)` | `gen_expr` |
| `gen_throw` | `(thr: Throw) → str` | Throw → `.error .tag` | — |
| `gen_imports` | `(funcs: list[SSAFunc]) → str` | 生成必要的 import 语句 | — |
| `gen_ub_comments` | `(obligations: list[UBObligation]) → str` | UB 目标以注释形式附加 | — |

```python
def generate_func(func: SSAFunc) -> str:
    """生成单个函数的 Lean 代码。

    格式：
    partial def {name}_ir {require_params} {value_params} : {ret_type} :=
      {body}

    规则：
    - require 参数放在值参数之前
    - 有 throw 时返回类型包 Except
    - TailRec 生成为嵌套的 partial def（where clause）
    """

def gen_expr(expr: ExprIR) -> str:
    """表达式翻译规则：

    Var(name, ver) → "{name}_{ver}"（或 "{name}" 如果 ver=0）
    Lit(value, UINT64) → f"({value} : UInt64)"
    BinOp("+", l, r) → f"({gen_expr(l)} + {gen_expr(r)})"
    BinOp("<<", l, r) → f"({gen_expr(l)} <<< {gen_expr(r)})"
    BinOp(">>", l, r) → f"({gen_expr(l)} >>> {gen_expr(r)})"
    BinOp("&&", l, r) → f"({gen_expr(l)} && {gen_expr(r)})"
    Call(func, args) → f"({func}_ir {' '.join(gen_expr(a) for a in args)} {ub_proofs})"
    ArrayAccess(arr, idx) → f"({gen_expr(arr)}[{gen_expr(idx)}]!)"
    FieldAccess(obj, field) → f"{gen_expr(obj)}.{field}"
    Cast(expr, UINT64) → f"({gen_expr(expr)}).toUInt64"
    ArrayPush(arr, elem) → f"({gen_expr(arr)}.push {gen_expr(elem)})"
    CondExpr(c, t, e) → f"(if {gen_expr(c)} then {gen_expr(t)} else {gen_expr(e)})"
    UnknownExpr(kind) → f"/- UNKNOWN: {kind} -/ sorry"
    """

def gen_tailrec(tr: TailRec, indent: int) -> str:
    """TailRec → partial def 尾递归。

    输出格式：
    partial def {tr.func_name} {params} : {ret_type} :=
      if {exit_cond} then {acc}
      else if {break_cond} then {acc}    -- 有 break 时
      else
        {body}
        {tr.func_name} {step} {new_acc}
    """
```

### 2.4 ub_collector — UB 证明目标收集器

```
文件：ub_collector.py
依赖：ir_types.py
```

| 函数 | 签名 | 功能 | 调用 |
|------|------|------|------|
| `collect_all` | `(func: SSAFunc) → list[UBObligation]` | 入口：收集函数中所有 UB 目标 | `scan_stmt`, `scan_expr` |
| `scan_stmt` | `(stmt: StmtIR, ctx: list[str]) → list[UBObligation]` | 扫描语句中的 UB 点 | `scan_expr` |
| `scan_expr` | `(expr: ExprIR, ctx: list[str]) → list[UBObligation]` | 扫描表达式中的 UB 点 | — |
| `make_div_check` | `(divisor: ExprIR, ctx) → UBObligation` | 除法 → b ≠ 0 | `gen_lean_prop` |
| `make_oob_check` | `(arr, idx: ExprIR, ctx) → UBObligation` | 数组 → i < size | `gen_lean_prop` |
| `make_shift_check` | `(amount: ExprIR, ctx) → UBObligation` | 移位 → n < 64 | `gen_lean_prop` |
| `gen_lean_prop` | `(ub_type, operands) → str` | 生成 Lean Prop 字符串 | — |

```python
def scan_expr(expr: ExprIR, ctx: list[str]) -> list[UBObligation]:
    """扫描表达式中的 UB 点。

    BinOp("/", _, rhs) → make_div_check(rhs)
    BinOp("%", _, rhs) → make_div_check(rhs)
    BinOp("<<", _, rhs) → make_shift_check(rhs)
    BinOp(">>", _, rhs) → make_shift_check(rhs)
    ArrayAccess(arr, idx) → make_oob_check(arr, idx)

    注意：
    - unsigned +/-/* 无 UB（mod 2^64），不生成
    - signed +/-/* 有溢出 UB，但 CLPoly 极少用 signed 算术
    - assert 已在 ssa_transform 中转为 Require，此处不重复
    """
```

### 2.5 refine_template — 精化定理骨架生成器

```
文件：refine_template.py
依赖：ir_types.py
```

| 函数 | 签名 | 功能 | 调用 |
|------|------|------|------|
| `generate_refinements` | `(l1_funcs, l2_mapping) → str` | 入口：生成精化定理 .lean 文件 | `gen_refinement` |
| `gen_refinement` | `(l1: SSAFunc, l2_name: str, bridge: str) → str` | 单个函数的精化定理骨架 | — |
| `load_l2_mapping` | `(config_path: str) → dict` | 加载 L1→L2 函数对应配置 | — |

```python
def gen_refinement(l1: SSAFunc, l2_name: str, bridge: str) -> str:
    """生成精化定理骨架。

    无 throw：
    theorem {l1.name}_refines {params} {requires} :
        {bridge}({l1.name}_ir {args}) = {l2_name} ({bridge_args}) := by
      sorry

    有 throw：
    theorem {l1.name}_refines {params} {requires} :
        {l1.name}_ir {args} = .ok result →
        {bridge}(result) = {l2_name} ({bridge_args}) := by
      sorry
    """
```

**L1→L2 函数对应配置**（`l1_l2_mapping.toml`）：

```toml
[nmod_mul]
l1 = "nmod_mul_ir"
l2 = "(· * ·) : ZMod p → ZMod p → ZMod p"
bridge = "UInt64.val % p.val"

[squarefree_Zp]
l1 = "squarefree_Zp_ir"
l2 = "sqfZp"
bridge = "SparsePolyZp.toPoly"
```

### 2.6 back2back — 背靠背测试框架

```
文件：back2back.py
依赖：ir_types.py, subprocess, json, tempfile
```

| 函数 | 签名 | 功能 | 调用 |
|------|------|------|------|
| `run_all_tests` | `(test_vectors, lean_file, cpp_exe) → Report` | 入口：运行全部测试 | `run_cpp`, `run_lean`, `compare` |
| `run_cpp` | `(cpp_exe: str, inputs: dict) → dict` | 执行 C++ 测试程序 | subprocess |
| `run_lean` | `(lean_file: str, func: str, inputs: dict) → dict` | 生成 #eval 文件并执行 | subprocess, tempfile |
| `compare` | `(cpp_out, lean_out) → bool` | 比较两个输出 | — |
| `gen_eval_file` | `(func_name: str, inputs: dict) → str` | 生成临时 .lean 文件（#eval 调用） | — |
| `gen_test_vectors` | `(func_name: str, n: int) → list[dict]` | 生成随机测试向量 | random |

```python
def gen_eval_file(func_name: str, inputs: dict) -> str:
    """生成可执行的 Lean 文件。

    import CLPoly.L1.{module}

    #eval do
      let result := {func_name}_ir {arg1} {arg2} ... {sorry_proofs}
      IO.println s!"{result}"

    注意：require 参数用 `by decide` / `by omega` 证明。
    对于具体测试值（如 p=13），这些命题是可判定的，Lean 可自动证。
    不需要 sorry。
    """

def compare(cpp_out: dict, lean_out: dict) -> bool:
    """比较输出。

    整数：精确比较
    多项式：按 (degree, coeff) 对列表比较
    数组：逐元素比较
    Except：.ok 时比较内容，.error 时比较 tag
    """
```

## 3. 调用关系总图

```
cpp2lean.py (主入口)
│
├── clang_ast.parse_translation_unit(json)
│   └─→ list[FuncIR]
│
├── ssa_transform.transform_func(func)  ×N
│   └─→ list[SSAFunc]
│
├── lean_codegen.generate_file(ssa_funcs)
│   └─→ output.lean（可编译+可执行）
│
├── ub_collector.collect_all(ssa_func)  ×N
│   └─→ list[UBObligation]（嵌入 output.lean 注释）
│
├── refine_template.generate_refinements(ssa_funcs, l2_mapping)
│   └─→ refinement_skeleton.lean（全 sorry）
│
└── back2back.run_all_tests(vectors, output.lean, cpp_exe)
    └─→ test_report.json
```

## 4. 与已有代码的复用

| 已有代码 | 复用方式 |
|---------|---------|
| `cpp2lean.py` 原型 | `translate_expr` 框架 → `clang_ast.parse_expr` 基础 |
| `cpp2lean.py` 原型 | `translate_stmt` 框架 → `clang_ast.parse_stmt` 基础 |
| `extract_test.cpp` | 背靠背测试的 C++ 测试桩 |
| CLPoly test suite | 测试向量来源（`test/test_factorize_*.cc`） |
| L2 Lean 文件 | `refine_template` 的 L2 函数签名来源 |

## 5. 错误处理策略

| 场景 | 处理方式 |
|------|---------|
| 未知 AST 节点 | `UnknownStmt`/`UnknownExpr` 保留，生成 `sorry` + 警告 |
| 未知类型 | 保留原始 C++ 类型名字符串，生成 `sorry` + 警告 |
| SSA 变换中遇到指针算术 | 生成 `UnknownStmt`，报错"超出 CLPoly C++ 子集" |
| Lean 编译失败 | 输出错误位置 + 对应 C++ 源码行号，辅助人工修复 |
| 背靠背测试不一致 | 输出差异详情（函数名、输入、C++ 输出、Lean 输出） |
| 背靠背测试中的 require 参数 | 不影响主体翻译。测试时用具体值 + `by decide` / `by omega` |

## 6. 开发顺序

按架构阶段确定的 Phase 0-4，细化为函数级实现顺序：

### Phase 0（基础设施）

1. `ir_types.py` — 全部数据结构
2. `clang_ast.py` — `parse_type`, `parse_expr`, `parse_param`, `parse_func_decl`, `parse_translation_unit`
3. `lean_codegen.py` — `gen_type`, `gen_expr`, `gen_signature`, `gen_body`（仅 LetStmt/ReturnStmt）
4. `cpp2lean.py` — 组合 1-3，端到端翻译纯算术函数
5. **验收**：`nmod_mul` 自动翻译 → `lake env lean` 编译通过

### Phase 1（控制流）

6. `clang_ast.py` — `parse_for`, `parse_while`, `parse_stmt`（IfStmt, break, continue）
7. `ssa_transform.py` — `VarEnv`, `transform_body`, `rename_vars`, `new_version`
8. `ssa_transform.py` — `transform_loop`, `identify_loop_vars`
9. `ssa_transform.py` — `transform_if`（fork/merge）
10. `ssa_transform.py` — `transform_params`, `collect_requires`, `detect_throws`
11. `lean_codegen.py` — `gen_tailrec`, `gen_require`, `gen_throw`
12. **验收**：`upoly_make_monic` 自动翻译

### Phase 2（数据结构）

13. `clang_ast.py` — 模板单态化处理、struct/class 成员访问
14. `ir_types.py` — `SparsePolyZp` 等 CLPoly 特有类型的 Lean 定义
15. `lean_codegen.py` — struct 生成、ArrayPush 处理
16. **验收**：`__squarefree_Zp` 翻译骨架

### Phase 3（完整管线+测试）

17. `ub_collector.py` — 全部 UB 检查函数
18. `refine_template.py` — 精化骨架生成
19. `back2back.py` — 测试框架
20. `cpp2lean.py` — 跨函数调用图排序
21. **验收**：`polynomial_factorize_zp.hh` 全翻译 + 背靠背测试
