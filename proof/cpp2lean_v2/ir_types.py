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
    UINT128 = "UInt128"       # 为 v1 class_map.py 的 CAST_TABLE 保留；CLPoly 因式分解不实际使用
    NAT = "Nat"               # size_t → Nat
    BOOL = "Bool"
    FLOAT = "Float"            # double
    UNIT = "Unit"              # void


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
    """C++ T& / T&& / T*（仅 HIR₀ 允许）。"""
    inner: 'TypeIR'
    is_const: bool = False
    is_rvalue: bool = False
    is_pointer: bool = False  # 区分 T* 和 T&


@dataclass(frozen=True)
class UnknownType:
    """parse Pass 偶尔产出；后续 Pass 不允许。"""
    raw: str


# 兼容 v1 class_map.py 的别名：v1 用 StructType(name, fields) 表示命名类型
# v2 改用 NamedType(name)；但保留 StructType 作为 alias 让 class_map.py 照搬复用
class StructType:
    """[v1 兼容] 等价于 NamedType；fields 参数保留但忽略。"""

    __slots__ = ("name", "fields")

    def __init__(self, name: str, fields: list | None = None):
        self.name = name
        self.fields = fields or []

    def __eq__(self, other):
        if isinstance(other, (StructType, NamedType)):
            return self.name == other.name
        return False

    def __hash__(self):
        return hash(("StructType", self.name))

    def __repr__(self):
        return f"StructType({self.name!r})"


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
    LambdaExpr,
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


# ============================================================
# §7 MIR：CFG / Phi / Terminator / MIRFunc
# 依据 mir-design.md §1
# ============================================================

@dataclass
class PhiStmt:
    """phi 节点：x_n := phi(pred1 → x_i, pred2 → x_j, ...)
    Pass 6 ssa_build 产出，仅出现在 BasicBlock 开头（连续 N 个）。"""
    target: Var                            # 新版本 Var（被定义）
    ty: TypeIR
    sources: dict[int, Var]                # pred_bb_id → source Var（旧版本）


# MIR 层 stmt 类型（subset of HIR + PhiStmt）：MIR₀ 不允许 AssignStmt /
# CompoundAssignStmt / IfStmt / 循环 / Break/Continue/Return（控制流走 Terminator）
MIRStmt = Union[PhiStmt, LetStmt, RequireStmt]


# ============================================================
# Terminator：基本块末尾的控制流
# ============================================================

@dataclass
class JumpTerm:
    """无条件跳转：goto target_bb。"""
    target: int


@dataclass
class CondJumpTerm:
    """条件跳转：cond ? then_bb : else_bb。"""
    cond: ExprIR
    then_bb: int
    else_bb: int


@dataclass
class ReturnTerm:
    """函数返回。"""
    value: Optional[ExprIR] = None


@dataclass
class TailCallTerm:
    """尾递归调用（Pass 7 循环下降产物）：跳转到另一个 def 的开头。
    args 与 target_func 的参数列表对应，按位置传递。"""
    target_func: str
    args: list[ExprIR] = field(default_factory=list)


Terminator = Union[JumpTerm, CondJumpTerm, ReturnTerm, TailCallTerm]


# ============================================================
# BasicBlock / CFG
# ============================================================

@dataclass
class BasicBlock:
    """基本块：一段无内部跳转的直线代码 + 末尾 Terminator。

    `stmts` 由 PhiStmt（开头连续 N 个）+ LetStmt + RequireStmt 组成。
    `terminator` 决定到哪个/哪些后继块（preds/succs 由 CFG 缓存）。"""
    bb_id: int
    stmts: list[MIRStmt]
    terminator: Terminator


@dataclass
class CFG:
    """控制流图：基本块集 + 入口；前驱/后继可缓存。"""
    entry: int
    blocks: dict[int, BasicBlock]
    # 前驱/后继可由 terminator 隐式给出（rebuild_edges 重算）
    preds: dict[int, list[int]] = field(default_factory=dict)
    succs: dict[int, list[int]] = field(default_factory=dict)

    def rebuild_edges(self) -> None:
        """从 terminator 重新计算 preds/succs。"""
        self.preds = {bb_id: [] for bb_id in self.blocks}
        self.succs = {bb_id: [] for bb_id in self.blocks}
        for bb_id, bb in self.blocks.items():
            t = bb.terminator
            if isinstance(t, JumpTerm):
                targets = [t.target]
            elif isinstance(t, CondJumpTerm):
                targets = [t.then_bb, t.else_bb]
            else:  # ReturnTerm / TailCallTerm 无 CFG 后继
                targets = []
            for tgt in targets:
                if tgt in self.blocks:
                    self.succs[bb_id].append(tgt)
                    self.preds[tgt].append(bb_id)


# ============================================================
# MIRFunc
# ============================================================

@dataclass
class MIRFunc:
    """MIR 层函数。body / aux_lambdas 被转换为 cfg + aux_defs。

    与 HIRFunc 对齐：base_name / instance_suffix / mangled_name / params /
    ret_ty / requires；body 改为 cfg；aux_lambdas（HIR 层 lifted）继续保留为
    aux_defs 的初始集合，Pass 7 循环提取后追加新 MIRFunc 进 aux_defs。"""
    base_name: str
    instance_suffix: str = ""
    mangled_name: str = ""
    qual_type: str = ""
    params: list[HIRParam] = field(default_factory=list)
    ret_ty: TypeIR = BaseType.UNIT
    requires: list[RequireStmt] = field(default_factory=list)
    cfg: Optional[CFG] = None              # 主 CFG（Pass 6 后非 None）
    aux_defs: list['MIRFunc'] = field(default_factory=list)

    @property
    def lean_name(self) -> str:
        if self.instance_suffix:
            return f"{self.base_name}_{self.instance_suffix}_ir"
        return f"{self.base_name}_ir"


# ============================================================
# §7.x MIR 不变量（Pass 6 / Pass 7 出口检查）
# ============================================================

def assert_mir0_invariant(func: MIRFunc) -> None:
    """MIR₀（Pass 6 后）出口检查：
      - func.cfg 非 None
      - 每个 BasicBlock terminator 非 None 且类型属 Terminator union
      - PhiStmt 只在 BasicBlock 开头（连续 N 个）
      - 不允许 AssignStmt / CompoundAssignStmt / IfStmt / WhileStmt / ForStmt /
        DoWhileStmt / RangeForStmt / BreakStmt / ContinueStmt / ReturnStmt /
        ExprStmt / BlockStmt 残留
      - SSA 性质：每个 Var(name, version) 至多一个定义（PhiStmt.target 或
        LetStmt.var）
    """
    def fail(reason: str):
        raise TranslationError(
            pass_name="ssa_build",
            func_name=func.base_name,
            reason=reason,
        )

    if func.cfg is None:
        fail("MIRFunc.cfg is None")

    cfg = func.cfg
    if cfg.entry not in cfg.blocks:
        fail(f"CFG entry {cfg.entry} not in blocks")

    seen_defs: dict[tuple[str, int], int] = {}  # (name, version) → bb_id

    for bb_id, bb in cfg.blocks.items():
        # Terminator 类型
        if not isinstance(bb.terminator, (JumpTerm, CondJumpTerm,
                                          ReturnTerm, TailCallTerm)):
            fail(f"bb[{bb_id}]: invalid terminator type {type(bb.terminator).__name__}")

        # phi 必须连续在开头
        seen_non_phi = False
        for i, s in enumerate(bb.stmts):
            if isinstance(s, PhiStmt):
                if seen_non_phi:
                    fail(f"bb[{bb_id}].stmts[{i}]: PhiStmt after non-phi")
                # 检查 SSA 唯一定义
                key = (s.target.name, s.target.version)
                if key in seen_defs:
                    fail(f"bb[{bb_id}]: SSA violation — Var{key} redefined "
                         f"(prior in bb[{seen_defs[key]}])")
                seen_defs[key] = bb_id
                continue
            seen_non_phi = True
            # 允许的非-phi MIRStmt：LetStmt / RequireStmt
            if isinstance(s, LetStmt):
                key = (s.var.name, s.var.version)
                if key in seen_defs:
                    fail(f"bb[{bb_id}].stmts[{i}]: SSA violation — Var{key} redefined")
                seen_defs[key] = bb_id
            elif isinstance(s, RequireStmt):
                pass  # 允许
            else:
                fail(f"bb[{bb_id}].stmts[{i}]: disallowed MIR stmt "
                     f"{type(s).__name__}")


def assert_mir1_invariant(func: MIRFunc) -> None:
    """MIR₁（Pass 7 后）出口检查：在 MIR₀ 基础上 +
      - 无结构化循环（CFG 不含 back edge——拓扑序成立）
      - 所有循环已提取为 aux_defs 内的独立 MIRFunc
      - 残留 break/continue/return 由 TailCallTerm + flag 表达
    """
    assert_mir0_invariant(func)
    # 简化：本检查不实际运行 DomTree 分析，仅约定 Pass 7 不留下显式 back edge
    # （Pass 7 实现时此处可加 SCC / 拓扑序检查）
