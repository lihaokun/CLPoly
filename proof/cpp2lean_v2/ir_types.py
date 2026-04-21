"""
cpp2lean v2 — HIR + MIR 数据结构

基于 hir-design.md §1 和 mir-design.md §1。
每个 HIRFunc 经过 Pass 1-5 变为 HIR₄，再经过 Pass 6-7 变为 MIR₁。
Pass 8 codegen 把 MIR₁ 转为 Lean 源码字符串。
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Union


# ============================================================
# §1 类型系统
# ============================================================

class BaseType(Enum):
    """基础数值类型（type-system.md §1）。"""
    UINT64 = "UInt64"
    INT64 = "Int64"
    INT32 = "Int32"
    UINT32 = "UInt32"
    NAT = "Nat"              # size_t → Nat
    BOOL = "Bool"
    FLOAT = "Float"           # double
    UNIT = "Unit"             # void


@dataclass(frozen=True)
class NamedType:
    """命名类型（模板实例化后的完整名）。"""
    name: str                 # "Zp", "ZZ", "SparsePolyZp", "MvPolyZZ", ...


@dataclass(frozen=True)
class ArrayType:
    elem: 'TypeIR'


@dataclass(frozen=True)
class PairType:
    fst: 'TypeIR'
    snd: 'TypeIR'


@dataclass(frozen=True)
class TupleType:
    """多元组，3+ 元素。2 元组用 PairType。"""
    elems: tuple['TypeIR', ...]


@dataclass(frozen=True)
class OptionType:
    inner: 'TypeIR'


@dataclass(frozen=True)
class StdMapType:
    key: 'TypeIR'
    value: 'TypeIR'


@dataclass(frozen=True)
class RefType:
    """C++ T& / T&&（仅 HIR₀ 允许）。"""
    inner: 'TypeIR'
    is_const: bool = False
    is_rvalue: bool = False


@dataclass(frozen=True)
class UnknownType:
    """parse Pass 偶尔产出；后续 Pass 不允许。"""
    raw: str


TypeIR = Union[
    BaseType, NamedType, ArrayType, PairType, TupleType,
    OptionType, StdMapType, RefType, UnknownType,
]


# ============================================================
# §2 表达式 (ExprIR)
# ============================================================

@dataclass
class Var:
    """变量引用。HIR: version=0; MIR: SSA 版本号 ≥ 0。"""
    name: str
    version: int = 0
    ty: Optional[TypeIR] = None

    def lean_name(self) -> str:
        return self.name if self.version == 0 else f"{self.name}_{self.version}"


@dataclass
class Lit:
    """字面量。"""
    value: object              # int / bool / float / str
    ty: TypeIR = BaseType.INT32


@dataclass
class BinOp:
    op: str                    # "+" "-" "*" "/" "%" "<" ">" "==" ...
    lhs: 'ExprIR'
    rhs: 'ExprIR'
    ty: Optional[TypeIR] = None


@dataclass
class UnaryOp:
    op: str                    # "!" "-" "++" "--" "~"
    operand: 'ExprIR'
    ty: Optional[TypeIR] = None


@dataclass
class CondExpr:
    """C++ 三元 ? :"""
    cond: 'ExprIR'
    then_e: 'ExprIR'
    else_e: 'ExprIR'
    ty: Optional[TypeIR] = None


@dataclass
class UnresolvedOp:
    """HIR₀-HIR₃ 中的未解析运算符/方法调用。
    operator_resolve Pass 替换为具体 Call。"""
    op_name: str                # "operator[]", "operator+", "size"
    receiver_ty: Optional[TypeIR] = None


@dataclass
class Call:
    """函数/方法调用。
    HIR₀-HIR₃: callee 可为 str 或 UnresolvedOp
    HIR₄+: callee 必为具体 Lean 函数名（str）"""
    callee: 'str | UnresolvedOp'
    args: list['ExprIR']
    ty: Optional[TypeIR] = None


@dataclass
class ArrayAccess:
    arr: 'ExprIR'
    idx: 'ExprIR'
    ty: Optional[TypeIR] = None


@dataclass
class FieldAccess:
    obj: 'ExprIR'
    field_name: str
    ty: Optional[TypeIR] = None


@dataclass
class Cast:
    """类型转换（显式或隐式）。"""
    expr: 'ExprIR'
    source_ty: TypeIR
    target_ty: TypeIR
    cast_kind: str              # "IntegralCast" / "NoOp" / ...


@dataclass
class Capture:
    """Lambda 捕获条目。"""
    name: str
    by_ref: bool                # [&x] True, [=x] False
    is_default: bool = False    # [&]/[=] 默认捕获


@dataclass
class LambdaExpr:
    """HIR₀-HIR₁ 内联 Lambda；lambda_lift Pass 之后消除。"""
    captures: list[Capture]
    params: list['HIRParam']
    body: list['StmtIR']
    ty: Optional[TypeIR] = None


@dataclass
class IteratorExpr:
    """HIR₀-HIR₂ 原始迭代器；iter_recognize Pass 之后消除。"""
    kind: str                   # "begin" / "end" / "deref" / "increment"
    container: Optional['ExprIR'] = None
    operand: Optional['ExprIR'] = None


@dataclass
class BlockExpr:
    """语句块作为表达式值（phi 节点材料）。"""
    stmts: list['StmtIR']
    value: 'ExprIR'
    ty: Optional[TypeIR] = None


@dataclass
class TupleExpr:
    """tuple 构造。"""
    elems: list['ExprIR']
    ty: Optional[TypeIR] = None


@dataclass
class ArrayLit:
    """Array 字面量 #[...]"""
    elems: list['ExprIR']
    elem_ty: TypeIR


@dataclass
class UnknownExpr:
    """parse Pass 的未识别表达式。"""
    kind: str
    children: list['ExprIR'] = field(default_factory=list)
    raw: str = ""


ExprIR = Union[
    Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
    ArrayAccess, FieldAccess, Cast,
    LambdaExpr, IteratorExpr,
    BlockExpr, TupleExpr, ArrayLit,
    UnknownExpr,
]


# ============================================================
# §3 语句 (StmtIR) — HIR 层
# ============================================================

@dataclass
class LetStmt:
    """let x : T := e"""
    var: Var
    ty: TypeIR
    value: ExprIR


@dataclass
class AssignStmt:
    """HIR: x = e (mutation)。Pass 6 ssa_build 消除。"""
    target: ExprIR              # Var / ArrayAccess / FieldAccess
    value: ExprIR


@dataclass
class CompoundAssignStmt:
    """HIR₀-HIR₃: x op= e。Pass 5 operator_resolve 展开为 AssignStmt + BinOp。"""
    target: ExprIR
    op: str                     # "+" "-" "*" "/" "%"
    value: ExprIR


@dataclass
class IfStmt:
    cond: ExprIR
    then_body: list['StmtIR']
    else_body: list['StmtIR'] = field(default_factory=list)


@dataclass
class WhileStmt:
    cond: ExprIR
    body: list['StmtIR']


@dataclass
class ForStmt:
    """C++ for (init; cond; step) { body }"""
    init: list['StmtIR']
    cond: ExprIR
    step: list['StmtIR']
    body: list['StmtIR']


@dataclass
class RangeForStmt:
    """C++ for (auto& x : container) { body }"""
    var: Var
    var_ty: TypeIR
    container: ExprIR
    body: list['StmtIR']
    decomposition: Optional[list[Var]] = None    # 结构化绑定 [k, v]


@dataclass
class DoWhileStmt:
    body: list['StmtIR']
    cond: ExprIR


@dataclass
class BreakStmt:
    pass


@dataclass
class ContinueStmt:
    pass


@dataclass
class ReturnStmt:
    value: Optional[ExprIR] = None


@dataclass
class RequireStmt:
    """assert / UB 前置条件。放到函数签名中。"""
    cond: ExprIR
    name: str                   # "hp" / "h_shift" / "hi"
    source: str                 # "assert" / "div_by_zero" / "array_oob" / ...
    uid: int = 0


@dataclass
class ExprStmt:
    """表达式作为语句（结果丢弃）。"""
    expr: ExprIR


@dataclass
class BlockStmt:
    """显式 { ... } 作用域。"""
    stmts: list['StmtIR']


@dataclass
class UnknownStmt:
    """parse Pass 的未识别语句；后续 Pass 遇到必须报错。"""
    kind: str
    children: list['StmtIR'] = field(default_factory=list)


StmtIR = Union[
    LetStmt, AssignStmt, CompoundAssignStmt,
    IfStmt, WhileStmt, ForStmt, RangeForStmt, DoWhileStmt,
    BreakStmt, ContinueStmt, ReturnStmt,
    RequireStmt, ExprStmt, BlockStmt,
    UnknownStmt,
]


# ============================================================
# §4 顶层结构
# ============================================================

@dataclass
class HIRParam:
    name: str
    ty: TypeIR
    is_ref: bool = False            # HIR₀: T&; HIR₁+: 必须 False
    is_const_ref: bool = False      # HIR₀: const T&; 保留
    is_output: bool = False         # 初始 = is_ref


@dataclass
class HIRFunc:
    """单个 Lean 定义对应的函数 IR。"""
    base_name: str                  # "factorize"
    instance_suffix: str             # "upoly" / "lex" / "grlex" / ""
    mangled_name: str                # Clang mangled name
    qual_type: str                   # 完整 C++ 签名（debug）
    params: list[HIRParam]
    ret_ty: TypeIR
    body: list[StmtIR]
    requires: list[RequireStmt] = field(default_factory=list)
    aux_lambdas: list['HIRFunc'] = field(default_factory=list)

    @property
    def lean_name(self) -> str:
        if self.instance_suffix:
            return f"{self.base_name}_{self.instance_suffix}_ir"
        return f"{self.base_name}_ir"


@dataclass
class HIRProgram:
    funcs: list[HIRFunc]


# ============================================================
# §5 MIR 层（Pass 6 之后）— 暂占位，Stage 2 Week 4 实现
# ============================================================

@dataclass
class PhiStmt:
    """SSA phi 节点：target := phi({pred_bb_id: src_var})"""
    target: Var
    ty: TypeIR
    sources: dict[int, Var]


@dataclass
class BasicBlock:
    """基本块：非跳转语句 + 末尾 terminator。"""
    bb_id: int
    stmts: list['MIRStmt']
    terminator: Optional['Terminator'] = None


@dataclass
class Terminator:
    pass


@dataclass
class JumpTerm(Terminator):
    target: int


@dataclass
class CondJumpTerm(Terminator):
    cond: ExprIR
    then_bb: int
    else_bb: int


@dataclass
class ReturnTerm(Terminator):
    value: Optional[ExprIR] = None


@dataclass
class TailCallTerm(Terminator):
    target_func: str
    args: list[ExprIR]


@dataclass
class CFG:
    entry: int
    blocks: dict[int, BasicBlock]
    preds: dict[int, list[int]] = field(default_factory=dict)
    succs: dict[int, list[int]] = field(default_factory=dict)


MIRStmt = Union[PhiStmt, LetStmt, RequireStmt]


@dataclass
class MIRFunc:
    """Pass 6/7 产物。"""
    base_name: str
    instance_suffix: str
    mangled_name: str
    params: list[HIRParam]
    ret_ty: TypeIR
    requires: list[RequireStmt]
    cfg: CFG
    aux_defs: list['MIRFunc'] = field(default_factory=list)

    @property
    def lean_name(self) -> str:
        if self.instance_suffix:
            return f"{self.base_name}_{self.instance_suffix}_ir"
        return f"{self.base_name}_ir"


# ============================================================
# §6 错误类型
# ============================================================

class TranslationError(Exception):
    """翻译失败。附上上下文方便 debug。"""

    def __init__(self, pass_name: str, func_name: str, reason: str, line: int = 0):
        self.pass_name = pass_name
        self.func_name = func_name
        self.reason = reason
        self.line = line
        super().__init__(f"[{pass_name}] {func_name}:{line}: {reason}")
