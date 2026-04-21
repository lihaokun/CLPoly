"""
Pass 1 parse 单元测试。

对 5 个 fixture 函数验证 Pass 1 的基础正确性。
"""

from __future__ import annotations
import sys
import json
import subprocess
from pathlib import Path

# Add v2 to path
V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, RefType, ArrayType, PairType, UnknownType,
    Var, Lit, BinOp, UnaryOp, Call, UnresolvedOp, FieldAccess,
    LetStmt, IfStmt, ReturnStmt, RequireStmt, ExprStmt, AssignStmt,
    HIRFunc, HIRParam,
    LambdaExpr, IteratorExpr,  # 不应在 HIR₀ 里出现于这些简单函数
    UnknownStmt, UnknownExpr,
)
from pass1_parse import parse_pass, assert_hir0_invariant, parse_type


PROJECT_ROOT = V2_ROOT.parent.parent
FIXTURES_DIR = V2_ROOT / "tests" / "fixtures"


def _system_includes() -> list[str]:
    """与 survey_ast.py 一致：抓系统 include。"""
    try:
        r = subprocess.run(
            ["g++", "-E", "-Wp,-v", "-x", "c++", "-"],
            input="", capture_output=True, text=True,
        )
        paths = []
        cap = False
        for line in r.stderr.splitlines():
            if "#include <...>" in line:
                cap = True
                continue
            if line.startswith("End of search"):
                break
            if cap and line.startswith(" "):
                paths.append(f"-I{line.strip()}")
        import glob
        for d in glob.glob("/usr/lib/llvm-*/lib/clang/*/include"):
            paths.append(f"-I{d}")
        return paths
    except Exception:
        return []


def dump_ast(fixture_cc: Path, func_name: str) -> dict | None:
    """Dump 单个函数的 AST JSON。"""
    cmd = [
        "clang++", "-std=c++17", f"-I{PROJECT_ROOT}",
    ] + _system_includes() + [
        "-Xclang", "-ast-dump=json",
        "-Xclang", f"-ast-dump-filter={func_name}",
        "-fsyntax-only",
        str(fixture_cc),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    if not r.stdout.strip():
        return None

    decoder = json.JSONDecoder()
    raw = r.stdout.strip()
    pos = 0
    candidates = []
    while pos < len(raw):
        try:
            obj, end = decoder.raw_decode(raw[pos:])
            candidates.append(obj)
            pos += end
            while pos < len(raw) and raw[pos] in " \t\n\r":
                pos += 1
        except json.JSONDecodeError:
            break

    def has_body(n: dict) -> bool:
        if n.get("kind") == "FunctionDecl":
            return any(
                c.get("kind") == "CompoundStmt"
                for c in n.get("inner", []) if isinstance(c, dict)
            )
        return False

    def expand(n: dict) -> list[dict]:
        if n.get("kind") == "FunctionTemplateDecl":
            return [
                c for c in n.get("inner", [])
                if isinstance(c, dict) and c.get("kind") == "FunctionDecl"
            ]
        return [n]

    exact = [c for c in candidates if c.get("name") == func_name]
    expanded = []
    for c in exact:
        expanded.extend(expand(c))

    mangled_with_body = [n for n in expanded if n.get("mangledName") and has_body(n)]
    if mangled_with_body:
        return mangled_with_body[-1]
    with_body = [n for n in expanded if has_body(n)]
    if with_body:
        return with_body[-1]
    return expanded[0] if expanded else None


# ============================================================
# 测试：parse_type 基础
# ============================================================

def test_parse_type_basic():
    assert parse_type("uint64_t") == BaseType.UINT64
    assert parse_type("int64_t") == BaseType.INT64
    assert parse_type("int") == BaseType.INT32
    assert parse_type("bool") == BaseType.BOOL
    assert parse_type("double") == BaseType.FLOAT
    assert parse_type("void") == BaseType.UNIT


def test_parse_type_ref():
    r = parse_type("const ZZ &")
    assert isinstance(r, RefType)
    assert r.is_const
    assert r.inner == NamedType("ZZ")

    r2 = parse_type("upolynomial_<ZZ> &")
    assert isinstance(r2, RefType)
    assert not r2.is_const
    assert r2.inner == NamedType("SparsePolyZZ")


def test_parse_type_clpoly():
    assert parse_type("ZZ") == NamedType("ZZ")
    assert parse_type("Zp") == NamedType("Zp")
    assert parse_type("upolynomial_<Zp>") == NamedType("SparsePolyZp")
    assert parse_type("upolynomial_<ZZ>") == NamedType("SparsePolyZZ")


def test_parse_type_stl():
    assert parse_type("std::vector<int>") == ArrayType(BaseType.INT32)
    assert parse_type("std::pair<int, long>") == PairType(BaseType.INT32, BaseType.INT64)


# ============================================================
# 测试：fixture 文件逐个通过 Pass 1
# ============================================================

def _test_fixture(func_name: str, fixture_file: str,
                  expected_params: list[tuple[str, type]],
                  min_body_stmts: int = 1):
    """通用 fixture 测试辅助：
    - 检查函数名、参数名/类型、body 非空。
    - 断言 HIR₀ 不变量（允许 Unknown）。
    """
    fixture_cc = FIXTURES_DIR / fixture_file
    assert fixture_cc.exists(), f"fixture missing: {fixture_cc}"

    ast = dump_ast(fixture_cc, func_name)
    assert ast is not None, f"Clang dump returned nothing for {func_name}"
    assert ast.get("kind") == "FunctionDecl", f"top-level not FunctionDecl: {ast.get('kind')}"

    hir = parse_pass(ast)
    assert hir.base_name == func_name
    assert len(hir.params) == len(expected_params), \
        f"{func_name}: expected {len(expected_params)} params, got {len(hir.params)}"
    for i, (pname, ptype_cls) in enumerate(expected_params):
        actual = hir.params[i]
        assert actual.name == pname, f"param {i}: expected {pname}, got {actual.name}"
        # 类型检查（按 class）
        ty = actual.ty
        assert isinstance(ty, ptype_cls) or (isinstance(ptype_cls, type) and isinstance(ty, ptype_cls)) \
            or ty == ptype_cls, f"param {i} type: expected {ptype_cls}, got {ty}"
    assert len(hir.body) >= min_body_stmts, \
        f"{func_name}: expected at least {min_body_stmts} body stmts, got {len(hir.body)}"

    assert_hir0_invariant(hir)
    return hir


def test_fixture_make_zp():
    hir = _test_fixture(
        func_name="__make_zp",
        fixture_file="make_zp.cc",
        expected_params=[("val", BaseType), ("p", BaseType)],
        min_body_stmts=1,
    )
    # body 应是 ReturnStmt（可能先被 BlockStmt 包裹；parse_pass 已扁平化）
    stmts = hir.body
    ret_stmts = [s for s in stmts if isinstance(s, ReturnStmt)]
    assert len(ret_stmts) == 1, f"expected 1 return, got {len(ret_stmts)}"
    assert ret_stmts[0].value is not None
    # 返回类型应为 Zp
    assert hir.ret_ty == NamedType("Zp"), f"ret_ty: {hir.ret_ty}"


def test_fixture_upoly_mod():
    hir = _test_fixture(
        func_name="__upoly_mod",
        fixture_file="upoly_mod.cc",
        expected_params=[("f", NamedType), ("g", NamedType)],  # const T& → NamedType
        min_body_stmts=2,
    )
    # 应有 return
    ret_stmts = [s for s in hir.body if isinstance(s, ReturnStmt)]
    assert len(ret_stmts) == 1


def test_fixture_upoly_divmod():
    hir = _test_fixture(
        func_name="__upoly_divmod",
        fixture_file="upoly_divmod.cc",
        expected_params=[
            ("q", NamedType), ("r", NamedType),
            ("f", NamedType), ("g", NamedType),
        ],
        min_body_stmts=1,
    )
    # q 和 r 应标为 ref
    assert hir.params[0].is_ref, "q should be ref"
    assert hir.params[1].is_ref, "r should be ref"
    # f 和 g 应标为 const_ref
    assert hir.params[2].is_const_ref, "f should be const_ref"
    assert hir.params[3].is_const_ref, "g should be const_ref"


def test_fixture_symmetric_mod():
    hir = _test_fixture(
        func_name="__symmetric_mod",
        fixture_file="symmetric_mod.cc",
        expected_params=[("a", NamedType), ("m", NamedType)],  # const ZZ&
        min_body_stmts=3,
    )
    # 应含 IfStmt
    if_stmts = [s for s in hir.body if isinstance(s, IfStmt)]
    assert len(if_stmts) >= 1, f"expected IfStmt, got {[type(s).__name__ for s in hir.body]}"


def test_fixture_upoly_const_term():
    hir = _test_fixture(
        func_name="__upoly_const_term",
        fixture_file="upoly_const_term.cc",
        expected_params=[("f", NamedType)],
        min_body_stmts=2,
    )
    # 应含多个 IfStmt
    if_stmts = [s for s in hir.body if isinstance(s, IfStmt)]
    assert len(if_stmts) >= 2, f"expected ≥2 IfStmt, got {len(if_stmts)}"


if __name__ == "__main__":
    # pytest-less 简易 runner
    tests = [
        test_parse_type_basic,
        test_parse_type_ref,
        test_parse_type_clpoly,
        test_parse_type_stl,
        test_fixture_make_zp,
        test_fixture_upoly_mod,
        test_fixture_upoly_divmod,
        test_fixture_symmetric_mod,
        test_fixture_upoly_const_term,
    ]
    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL {t.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
            failed += 1
    print(f"\n=== {passed}/{passed+failed} passed ===")
    sys.exit(0 if failed == 0 else 1)
