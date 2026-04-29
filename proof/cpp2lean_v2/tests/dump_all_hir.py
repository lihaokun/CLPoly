"""
Pass 1 HIR₀ dumper — 把 65 函数的 HIR 转为可读文本供等价性审核。

输出到 /tmp/hir0_dump/<func_name>.txt
"""

from __future__ import annotations
import sys
import json
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, RefType, ArrayType, PairType, TupleType,
    StdMapType, OptionType, UnknownType, TypeIR,
    Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
    ArrayAccess, FieldAccess, Cast, Capture, LambdaExpr,
    BlockExpr, TupleExpr, ArrayLit, UnknownExpr,
    LetStmt, AssignStmt, CompoundAssignStmt, IfStmt, WhileStmt, ForStmt,
    RangeForStmt, DoWhileStmt, BreakStmt, ContinueStmt, ReturnStmt,
    RequireStmt, ExprStmt, BlockStmt, UnknownStmt,
    HIRFunc, HIRParam,
)
from pass1_parse import parse_pass
from class_map import TRANSLATION_SCOPE

AST_CACHE_DIR = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
OUT_DIR = Path("/tmp/hir0_dump")


def fmt_type(ty: TypeIR, depth=0) -> str:
    if isinstance(ty, BaseType):
        return ty.value
    if isinstance(ty, NamedType):
        return ty.name
    if isinstance(ty, ArrayType):
        return f"Array<{fmt_type(ty.elem)}>"
    if isinstance(ty, PairType):
        return f"Pair<{fmt_type(ty.fst)}, {fmt_type(ty.snd)}>"
    if isinstance(ty, TupleType):
        return f"Tuple<{', '.join(fmt_type(e) for e in ty.elems)}>"
    if isinstance(ty, StdMapType):
        return f"StdMap<{fmt_type(ty.key)}, {fmt_type(ty.value)}>"
    if isinstance(ty, OptionType):
        return f"Option<{fmt_type(ty.inner)}>"
    if isinstance(ty, RefType):
        cst = "const " if ty.is_const else ""
        suffix = "&&" if ty.is_rvalue else "&"
        return f"{cst}{fmt_type(ty.inner)}{suffix}"
    if isinstance(ty, UnknownType):
        return f"?[{ty.raw[:30]}]"
    return str(ty)


def fmt_expr(e, depth=0) -> str:
    if isinstance(e, Var):
        return e.name if e.version == 0 else f"{e.name}_{e.version}"
    if isinstance(e, Lit):
        return repr(e.value)
    if isinstance(e, BinOp):
        return f"({fmt_expr(e.lhs)} {e.op} {fmt_expr(e.rhs)})"
    if isinstance(e, UnaryOp):
        return f"({e.op}{fmt_expr(e.operand)})"
    if isinstance(e, CondExpr):
        return f"({fmt_expr(e.cond)} ? {fmt_expr(e.then_e)} : {fmt_expr(e.else_e)})"
    if isinstance(e, UnresolvedOp):
        return f"<{e.op_name}>"
    if isinstance(e, Call):
        callee = e.callee if isinstance(e.callee, str) else fmt_expr(e.callee)
        args = ", ".join(fmt_expr(a) for a in e.args)
        return f"{callee}({args})"
    if isinstance(e, ArrayAccess):
        return f"{fmt_expr(e.arr)}[{fmt_expr(e.idx)}]"
    if isinstance(e, FieldAccess):
        return f"{fmt_expr(e.obj)}.{e.field_name}"
    if isinstance(e, Cast):
        return f"cast<{fmt_type(e.target_ty)}>({fmt_expr(e.expr)})[{e.cast_kind}]"
    if isinstance(e, LambdaExpr):
        return f"LAMBDA(captures={e.captures}, params={len(e.params)})"
    if isinstance(e, BlockExpr):
        return f"BLOCK(...;{fmt_expr(e.value)})"
    if isinstance(e, TupleExpr):
        return f"({', '.join(fmt_expr(x) for x in e.elems)})"
    if isinstance(e, ArrayLit):
        return f"#[{', '.join(fmt_expr(x) for x in e.elems)}]"
    if isinstance(e, UnknownExpr):
        return f"UNKNOWN_EXPR<{e.kind}>"
    return str(e)


def fmt_stmt(s, indent=0) -> str:
    pad = "  " * indent
    if isinstance(s, LetStmt):
        return f"{pad}let {s.var.name} : {fmt_type(s.ty)} := {fmt_expr(s.value)}"
    if isinstance(s, AssignStmt):
        return f"{pad}{fmt_expr(s.target)} = {fmt_expr(s.value)}"
    if isinstance(s, CompoundAssignStmt):
        return f"{pad}{fmt_expr(s.target)} {s.op}= {fmt_expr(s.value)}"
    if isinstance(s, IfStmt):
        out = [f"{pad}if {fmt_expr(s.cond)} {{"]
        for t in s.then_body:
            out.append(fmt_stmt(t, indent + 1))
        if s.else_body:
            out.append(f"{pad}}} else {{")
            for t in s.else_body:
                out.append(fmt_stmt(t, indent + 1))
        out.append(f"{pad}}}")
        return "\n".join(out)
    if isinstance(s, WhileStmt):
        out = [f"{pad}while {fmt_expr(s.cond)} {{"]
        for t in s.body:
            out.append(fmt_stmt(t, indent + 1))
        out.append(f"{pad}}}")
        return "\n".join(out)
    if isinstance(s, ForStmt):
        out = [f"{pad}for (init;cond={fmt_expr(s.cond)};step) {{"]
        for t in s.init:
            out.append(f"{pad}  [init] {fmt_stmt(t, 0)}")
        for t in s.body:
            out.append(fmt_stmt(t, indent + 1))
        for t in s.step:
            out.append(f"{pad}  [step] {fmt_stmt(t, 0)}")
        out.append(f"{pad}}}")
        return "\n".join(out)
    if isinstance(s, RangeForStmt):
        decomp = f" [decomp={[v.name for v in s.decomposition]}]" if s.decomposition else ""
        out = [f"{pad}for ({s.var.name} : {fmt_type(s.var_ty)} in {fmt_expr(s.container)}){decomp} {{"]
        for t in s.body:
            out.append(fmt_stmt(t, indent + 1))
        out.append(f"{pad}}}")
        return "\n".join(out)
    if isinstance(s, DoWhileStmt):
        out = [f"{pad}do {{"]
        for t in s.body:
            out.append(fmt_stmt(t, indent + 1))
        out.append(f"{pad}}} while {fmt_expr(s.cond)}")
        return "\n".join(out)
    if isinstance(s, BreakStmt):
        return f"{pad}break"
    if isinstance(s, ContinueStmt):
        return f"{pad}continue"
    if isinstance(s, ReturnStmt):
        v = fmt_expr(s.value) if s.value else ""
        return f"{pad}return {v}"
    if isinstance(s, RequireStmt):
        return f"{pad}require {s.name} : {fmt_expr(s.cond)}  [{s.source}]"
    if isinstance(s, ExprStmt):
        return f"{pad}(expr) {fmt_expr(s.expr)}"
    if isinstance(s, BlockStmt):
        out = [f"{pad}{{"]
        for t in s.stmts:
            out.append(fmt_stmt(t, indent + 1))
        out.append(f"{pad}}}")
        return "\n".join(out)
    if isinstance(s, UnknownStmt):
        return f"{pad}UNKNOWN_STMT<{s.kind}> ({len(s.children)} children)"
    return f"{pad}???({type(s).__name__})"


def fmt_func(func: HIRFunc) -> str:
    out = []
    out.append(f"# {func.base_name}" + (f" [{func.instance_suffix}]" if func.instance_suffix else ""))
    out.append(f"mangled: {func.mangled_name}")
    out.append(f"qualType: {func.qual_type}")
    out.append("")
    out.append("## Params")
    for i, p in enumerate(func.params):
        flags = []
        if p.is_ref: flags.append("REF")
        if p.is_const_ref: flags.append("CONST-REF")
        if p.is_output: flags.append("OUTPUT")
        out.append(f"  [{i}] {p.name} : {fmt_type(p.ty)} {' '.join(f'[{f}]' for f in flags)}")
    out.append("")
    out.append(f"## Return type: {fmt_type(func.ret_ty)}")
    out.append("")
    out.append("## Body")
    for s in func.body:
        out.append(fmt_stmt(s, 0))
    return "\n".join(out)


def _dump_factorize_instances(out_dir: Path):
    """对 factorize 特殊处理：dump 全部 3 个实例（upoly/lex/grlex）。"""
    import subprocess
    PROJECT_ROOT = V2_ROOT.parent.parent

    # 内联 system_includes：
    def sys_incs():
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

    cmd = [
        "clang++", "-std=c++17", f"-I{PROJECT_ROOT}",
    ] + sys_incs() + [
        "-Xclang", "-ast-dump=json",
        "-Xclang", "-ast-dump-filter=factorize",
        "-fsyntax-only",
        str(V2_ROOT.parent / "cpp2lean" / "instantiate.cc"),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

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

    # 展开 FunctionTemplateDecl 的实例化
    def expand(n: dict) -> list[dict]:
        if n.get("kind") == "FunctionTemplateDecl":
            return [c for c in n.get("inner", [])
                    if isinstance(c, dict) and c.get("kind") == "FunctionDecl"]
        if n.get("kind") == "FunctionDecl":
            return [n]
        return []

    def has_body(n: dict) -> bool:
        return any(c.get("kind") == "CompoundStmt"
                   for c in n.get("inner", []) if isinstance(c, dict))

    # 只取 name == "factorize" 且有 mangled + body
    instances = []
    for cand in candidates:
        if cand.get("name") == "factorize":
            for inst in expand(cand):
                if inst.get("mangledName") and has_body(inst):
                    instances.append(inst)

    # 按 qualType 分类生成后缀
    def suffix_of(qt: str) -> str:
        if "upolynomial_<ZZ>" in qt: return "upoly"
        if "grlex_<less>" in qt: return "grlex"
        if "lex_<less>" in qt: return "lex"
        if "QQ," in qt: return "qq"
        return "unknown"

    dumped = 0
    for inst in instances:
        qt = inst.get("type", {}).get("qualType", "")
        sfx = suffix_of(qt)
        try:
            hir = parse_pass(inst)
            out_file = out_dir / f"factorize_{sfx}.txt"
            out_file.write_text(fmt_func(hir))
            dumped += 1
        except Exception as e:
            print(f"FAIL factorize_{sfx}: {type(e).__name__}: {e}")
    return dumped


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    count = 0
    for func_name in sorted(TRANSLATION_SCOPE):
        if func_name == "factorize":
            n = _dump_factorize_instances(OUT_DIR)
            print(f"Dumped {n} factorize instances")
            count += n
            # 也写一个 factorize.txt 保持兼容（取 grlex 实例）
            cache_file = AST_CACHE_DIR / f"{func_name}.json"
            if cache_file.exists():
                try:
                    with open(cache_file) as f:
                        ast = json.load(f)
                    hir = parse_pass(ast)
                    (OUT_DIR / "factorize.txt").write_text(fmt_func(hir))
                except Exception:
                    pass
            continue

        cache_file = AST_CACHE_DIR / f"{func_name}.json"
        if not cache_file.exists():
            print(f"SKIP {func_name}: no AST cache")
            continue
        try:
            with open(cache_file) as f:
                ast = json.load(f)
            hir = parse_pass(ast)
            out_file = OUT_DIR / f"{func_name}.txt"
            out_file.write_text(fmt_func(hir))
            count += 1
        except Exception as e:
            print(f"FAIL {func_name}: {type(e).__name__}: {e}")
    print(f"\nDumped {count} items to {OUT_DIR}")


if __name__ == "__main__":
    main()
