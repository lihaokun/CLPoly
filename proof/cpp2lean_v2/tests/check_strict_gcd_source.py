#!/usr/bin/env python3
"""Admission gate for the raw Euclidean specialization of C++ ``gcd``.

The strict Lean implementation may replace object moves by raw-buffer pointer
rotation only while this exact source control-flow shape remains unchanged.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))
sys.path.insert(0, str(V2_ROOT / "tests"))

from build_strict_gcd import dense_method_ast
from check_strict_divrem import stable_ast
from ir_types import Call, ExprStmt, IfStmt, WhileStmt
from pass1_parse import parse_pass


AST_SHA256 = "1b6f7b4200ddd82c6c7a632c2b89c87cf7fd8781e1eee8816597db630b305656"


def walk_statements(stmts):
    for stmt in stmts:
        yield stmt
        if isinstance(stmt, IfStmt):
            yield from walk_statements(stmt.then_body)
            yield from walk_statements(stmt.else_body)
        elif isinstance(stmt, WhileStmt):
            yield from walk_statements(stmt.body)


def main() -> None:
    ast = dense_method_ast("gcd")
    digest = hashlib.sha256(json.dumps(
        stable_ast(ast), sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if digest != AST_SHA256:
        raise SystemExit(
            f"dense_upoly_zp::gcd AST drift: expected {AST_SHA256}, got {digest}")

    hir = parse_pass(ast)
    loops = [stmt for stmt in walk_statements(hir.body)
             if isinstance(stmt, WhileStmt)]
    if len(loops) != 1:
        raise SystemExit(f"gcd must contain exactly one Euclid while loop: {len(loops)}")
    loop = loops[0]
    calls = [stmt.expr for stmt in loop.body
             if isinstance(stmt, ExprStmt) and isinstance(stmt.expr, Call)]
    direct = [call.callee for call in calls if isinstance(call.callee, str)]
    if direct.count("divrem") != 1:
        raise SystemExit(f"gcd Euclid loop no longer has one divrem call: {direct}")

    rendered = repr(loop.body)
    required_moves = (
        "Var(name='a'", "Call(callee='move', args=[Var(name='b'",
        "Var(name='b'", "Call(callee='move', args=[Var(name='r'",
    )
    missing = [fragment for fragment in required_moves if fragment not in rendered]
    if missing:
        raise SystemExit(f"gcd Euclid buffer rotation drift: missing {missing}")
    if "<method>._gcd_hgcd" not in repr(hir.body):
        raise SystemExit("gcd no longer exposes the HGCD branch separately")
    print("PASS: gcd AST has the admitted Euclid divrem/move rotation and HGCD split")


if __name__ == "__main__":
    main()
