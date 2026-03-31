"""
CLPoly C++ → Lean IR 翻译器：共享数据结构

所有模块共享的 IR 类型定义。
"""

from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional, Union

# ============================================================
# 类型系统
# ============================================================

class BaseType(Enum):
    UINT64 = "UInt64"
    INT64 = "Int"
    UINT32 = "UInt32"
    UINT128 = "UInt128"
    BOOL = "Bool"
    VOID = "Unit"

@dataclass
class ArrayType:
    elem: TypeIR

@dataclass
class PairType:
    fst: TypeIR
    snd: TypeIR

@dataclass
class StructType:
    name: str
    fields: list[tuple[str, TypeIR]]

@dataclass
class ExceptType:
    inner: TypeIR

TypeIR = Union[BaseType, ArrayType, PairType, StructType, ExceptType, str]

# ============================================================
# 表达式
# ============================================================

@dataclass
class Var:
    name: str
    version: int = 0

    def lean_name(self) -> str:
        if self.version == 0:
            return self.name
        return f"{self.name}_{self.version}"

@dataclass
class Lit:
    value: int
    typ: BaseType = BaseType.UINT64

@dataclass
class BinOp:
    op: str
    lhs: ExprIR
    rhs: ExprIR

@dataclass
class UnaryOp:
    op: str
    operand: ExprIR

@dataclass
class CondExpr:
    cond: ExprIR
    then_e: ExprIR
    else_e: ExprIR

@dataclass
class Call:
    func: str
    args: list[ExprIR]

@dataclass
class ArrayAccess:
    arr: ExprIR
    idx: ExprIR

@dataclass
class FieldAccess:
    obj: ExprIR
    field_name: str

@dataclass
class Cast:
    expr: ExprIR
    target_type: TypeIR
    source_type: TypeIR = None  # 原始类型（如果已知）

@dataclass
class ArrayPush:
    arr: ExprIR
    elem: ExprIR

@dataclass
class UnknownExpr:
    kind: str
    children: list[ExprIR] = field(default_factory=list)
    raw: str = ""

ExprIR = Union[Var, Lit, BinOp, UnaryOp, CondExpr, Call, ArrayAccess,
               FieldAccess, Cast, ArrayPush, UnknownExpr]

# ============================================================
# 语句
# ============================================================

@dataclass
class LetStmt:
    var: Var
    typ: TypeIR
    value: ExprIR

@dataclass
class AssignStmt:
    target: Var
    value: ExprIR

@dataclass
class IfStmt:
    cond: ExprIR
    then_body: list[StmtIR]
    else_body: list[StmtIR] = field(default_factory=list)

@dataclass
class ReturnStmt:
    value: Optional[ExprIR] = None

@dataclass
class Require:
    cond: ExprIR
    name: str
    source: str  # "assert" | "div_by_zero" | "array_oob" | "shift_oob"

@dataclass
class Throw:
    error_tag: str
    message: str = ""

@dataclass
class TailRec:
    func_name: str
    params: list[tuple[Var, TypeIR]]
    exit_cond: ExprIR
    break_cond: Optional[ExprIR]
    body: list[StmtIR]
    step: list[StmtIR]

@dataclass
class ExprStmt:
    expr: ExprIR

@dataclass
class UnknownStmt:
    kind: str
    children: list[StmtIR] = field(default_factory=list)

StmtIR = Union[LetStmt, AssignStmt, IfStmt, ReturnStmt, Require, Throw,
               TailRec, ExprStmt, UnknownStmt]

# ============================================================
# 函数
# ============================================================

@dataclass
class ParamIR:
    name: str
    typ: TypeIR
    is_output: bool = False
    is_const_ref: bool = False

@dataclass
class FuncIR:
    name: str
    params: list[ParamIR]
    ret_type: TypeIR
    body: list[StmtIR]
    source_file: str = ""
    source_line: int = 0

@dataclass
class SSAFunc:
    name: str
    params: list[ParamIR]
    ret_type: TypeIR
    requires: list[Require]
    body: list[StmtIR]
    has_throw: bool = False

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
    lean_prop: str
    context: list[str]
