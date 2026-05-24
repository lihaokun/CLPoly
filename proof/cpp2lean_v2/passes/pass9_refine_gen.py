"""
Pass 9: refine_gen — 从 MIRFunc + REFINEMENT_MAP 生成精化定理骨架

输入：MIRFunc 列表（Pass 8 已生成 L1 corpus）
输出：proof/lean/CLPoly/Refinement/*.lean（全 sorry 骨架）

四类定理形态（由 result_kind 决定）：

  direct_eq    : l1_ir x = l2 x                        （纯数值函数）
  map_eq       : toPolyList (l1_ir f) p = l2 (f.toPoly p)  （列表输出，Array→List）
  map_eq_pair  : toPoly (l1_ir f).snd = l2 (f.toPoly p)    （pair 输出，取第二个分量）
  pair_eq      : toPoly (l1_ir ...).1 = (l2 ...).1 ∧ ...  （多返回值）
  conditional  : l1_ir f = .ok result → toPoly result = ...（有 throw/Except）

参数信息的三个来源（优先级从高到低）：
  1. MIRFunc（来自 Pass 1-8 流水线）
  2. REFINEMENT_MAP 的 l1_params / l2_call 字段
  3. _FIRST_PARAM / _HAS_SPARSE_POLY 固定表
"""

from __future__ import annotations
import sys
from pathlib import Path
from typing import Optional

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from class_map import REFINEMENT_MAP
from ir_types import (
    BaseType,
    MIRFunc,
    NamedType,
    ArrayType,
    PairType,
)


# ============================================================
# 输出目录
# ============================================================

REFINEMENT_DIR = Path(__file__).resolve().parent.parent.parent / "lean" / "CLPoly" / "Refinement"


# ============================================================
# 片段生成
# ============================================================

LEAN_KEYWORDS = frozenset({
    "let", "fun", "match", "with", "do", "if", "then", "else",
    "where", "by", "and", "or", "not", "in", "of", "from",
    "def", "theorem", "sorry",
})


def _safe_name(s: str) -> str:
    if s in LEAN_KEYWORDS:
        return f"«{s}»"
    return s


def _to_lean_type(ty) -> str:
    if isinstance(ty, BaseType):
        return ty.value
    if isinstance(ty, NamedType):
        return ty.name
    if isinstance(ty, ArrayType):
        return f"Array ({_to_lean_type(ty.elem)})"
    if isinstance(ty, PairType):
        return f"({_to_lean_type(ty.fst)} × {_to_lean_type(ty.snd)})"
    return " sorry /* type */ "


# 固定参数表（当 MIRFunc 和 REFINEMENT_MAP.l1_params 都不可用时使用）
_PARAM_TABLE = {
    "__make_zp":        [("val", "Int64"), ("modulus", "UInt64")],
    "__upoly_make_monic": [("f", "SparsePolyZp")],
    "__upoly_divmod":   [("q", "SparsePolyZp"), ("r", "SparsePolyZp"),
                         ("f", "SparsePolyZp"), ("g", "SparsePolyZp")],
    "__squarefree_Zp": [("f", "SparsePolyZp")],
    "__ddf_Zp": [("f", "SparsePolyZp")],
    "__symmetric_mod": [("a", "ZZ"), ("m", "ZZ")],
    "__binomial": [("n", "Int64"), ("k", "Int64")],
    "__isqrt_ceil": [("n", "ZZ")],
    "__upoly_mod":      [("f", "SparsePolyZp"), ("g", "SparsePolyZp")],
    "__edf_Zp":         [("result", "Array SparsePolyZp"), ("f", "SparsePolyZp"),
                         ("d", "UInt64"), ("rng", "Rng")],
    "__factor_Zp":      [("f", "SparsePolyZp")],
    "factorize":        [("F", "SparsePolyZZ")],
}

# 哪些 base_name 的参数包含 SparsePolyZp
_SPARSE_POLY_FUNCS = {
    "factorize", "__squarefree_Zp", "__ddf_Zp", "__edf_Zp",
    "__factor_Zp", "__upoly_make_monic", "__upoly_divmod",
}


def _get_params(base_name: str) -> list[tuple[str, str]]:
    return _PARAM_TABLE.get(base_name, [("x", "ZZ")])


def _has_sparse_poly(base_name: str) -> bool:
    return base_name in _SPARSE_POLY_FUNCS


def _first_sparse_poly_param(base_name: str) -> str:
    for name, ty in _PARAM_TABLE.get(base_name, []):
        if "SparsePolyZp" in ty:
            return name
    return "x"


def _emit_theorem_body(l1_name: str, info: dict, f: Optional[MIRFunc] = None,
                       base_name: str = "") -> str:
    """生成定理结论（= 号右侧的表达式）。"""
    kind = info["result_kind"]
    l2_name = info["l2_name"]
    l2_call = info.get("l2_call", l2_name)
    bn = base_name or l1_name.replace("_ir", "")

    # L1 调用参数（Lean 空格分隔函数参数）
    if f is not None:
        l1_call_params = " ".join(_safe_name(p.name) for p in f.params)
    else:
        l1_call_params = " ".join(_safe_name(p[0]) for p in _get_params(bn))
    l1_call = f"Generated.{l1_name} {l1_call_params}"

    if kind == "direct_eq":
        if _has_sparse_poly(bn):
            return f"  SparsePolyZp.toPoly p ({l1_call}) = {l2_call}"
        else:
            return f"  {l1_call} = {l2_call}"
    elif kind == "map_eq":
        return f"  toPolyList ({l1_call}) p = {l2_call}"
    elif kind == "map_eq_pair":
        return f"  SparsePolyZp.toPoly p ({l1_call}).snd = SparsePolyZp.toPoly p ({l2_call})"
    elif kind == "pair_eq":
        return (
            f"  SparsePolyZp.toPoly p ({l1_call}).1 = SparsePolyZp.toPoly p ({l2_call}).1 ∧\n"
            f"  SparsePolyZp.toPoly p ({l1_call}).2 = SparsePolyZp.toPoly p ({l2_call}).2"
        )
    elif kind == "conditional":
        return (
            f"  {l1_call} = .ok result →\n"
            f"  SparsePolyZp.toPoly p result = {l2_call}"
        )
    else:
        return f"  sorry  /* unknown result_kind: {kind} */"


def _emit_refinement_theorem(f: Optional[MIRFunc], info: dict, base_name: str = "") -> str:
    if f is not None:
        l1_lean_name = f.lean_name
        bn = f.base_name
    else:
        l1_lean_name = info.get("l1_name", f"{base_name}_ir")
        bn = base_name

    cpp_source = info.get("cpp_source", "unknown")
    doc = info.get("doc", "")
    l2_name = info["l2_name"]

    lines: list[str] = []
    # docstring
    lines.append("/--")
    lines.append(f"  L1 `{l1_lean_name}` (C++: `{cpp_source}`) → L2 `{l2_name}`")
    if doc:
        lines.append(f"")
        lines.append(f"  {doc}")
    lines.append(f"")
    lines.append(f"  C++ source: `{cpp_source}` — Clang AST → cpp2lean v2 Pass 1-8 → Corpus.lean")
    lines.append(f"  L2 model : Algorithm/{info.get('refinement_file', '?')}.lean — hand-written, proven correct")
    lines.append(f"  Bridge   : SparsePolyZp.toPoly (see Math/Univariate.lean)")
    lines.append(f"")
    lines.append(f"  Proof status: Skeleton (sorry) — fill to complete L1→L2 verification chain")
    lines.append("-/")

    # 参数列表
    if f is not None:
        params = [(p.name, _to_lean_type(p.ty)) for p in f.params]
    else:
        params = info.get("l1_params", _get_params(bn))

    # 定理签名
    lines.append(f"theorem {l1_lean_name}_refines (p : ℕ) [hp : Fact (Nat.Prime p)]")
    for pname, ptype in params:
        lines.append(f"    ({_safe_name(pname)} : {ptype})")

    # 前提条件（SparsePolyZp 参数）
    has_spzp = any("SparsePolyZp" in t for _, t in params)
    for pname, ptype in params:
        if "SparsePolyZp" in ptype:
            lines.append(f"    (hwf_{pname} : SparsePolyZp.WellFormed p {_safe_name(pname)})")
            lines.append(f"    (hred_{pname} : SparsePolyZp.AllReduced p {_safe_name(pname)}.toList)")
    if has_spzp:
        lines.append(f"    (hp_size : 2 * p ≤ UInt64.size)")

    # 结论
    body = _emit_theorem_body(l1_lean_name, info, f, base_name=bn)
    lines.append(f"    : {body} :=")
    lines.append(f"  by")
    lines.append(f"  sorry")
    lines.append("")

    return "\n".join(lines)


def _emit_refinement_file(funcs: list[Optional[MIRFunc]], base_names: list[str], module: str) -> str:
    lines: list[str] = []
    lines.append("-- Auto-generated by cpp2lean v2 Pass 9 (refine_gen)")
    lines.append("-- L1 → L2 refinement theorem skeletons")
    lines.append("")

    seen_imports: set[str] = set()
    for bn in base_names:
        imp = REFINEMENT_MAP.get(bn, {}).get("l2_import", "")
        if imp and imp not in seen_imports:
            lines.append(f"import {imp}")
            seen_imports.add(imp)
    lines.append("import CLPoly.Model")
    lines.append("import CLPoly.Generated.Corpus")
    lines.append("import CLPoly.Refinement.Basic")
    lines.append("import CLPoly.Math.Univariate")
    lines.append("")
    lines.append("set_option autoImplicit false")
    lines.append("")
    lines.append("open Polynomial")
    lines.append("open CLPoly.Math")
    lines.append("")
    lines.append("namespace Refinement")
    lines.append("")
    lines.append("variable {p : ℕ} [hp : Fact (Nat.Prime p)]")
    lines.append("")

    for i, bn in enumerate(base_names):
        f = funcs[i] if i < len(funcs) else None
        info = REFINEMENT_MAP.get(bn, {})
        theorem = _emit_refinement_theorem(f, info, base_name=bn)
        lines.append(theorem)

    lines.append("end Refinement")
    lines.append("")
    return "\n".join(lines)


# ============================================================
# Pass 9 入口
# ============================================================

def refine_gen_pass(top_funcs: Optional[list[MIRFunc]] = None) -> None:
    mir_index: dict[str, MIRFunc] = {}
    if top_funcs is not None:
        for f in top_funcs:
            mir_index[f.base_name] = f

    file_groups: dict[str, dict] = {}
    for base_name, info in REFINEMENT_MAP.items():
        rfile = info["refinement_file"]
        if rfile not in file_groups:
            file_groups[rfile] = {"info": info, "funcs": [], "base_names": []}
        file_groups[rfile]["base_names"].append(base_name)
        if base_name in mir_index:
            file_groups[rfile]["funcs"].append(mir_index[base_name])
        else:
            file_groups[rfile]["funcs"].append(None)

    REFINEMENT_DIR.mkdir(parents=True, exist_ok=True)
    for rfile, group in sorted(file_groups.items()):
        src = _emit_refinement_file(
            group["funcs"], group["base_names"], rfile)
        out_path = REFINEMENT_DIR / f"{rfile}.lean"
        out_path.write_text(src)
        n_thms = len(group["base_names"])
        print(f"  Refinement: {out_path} ({n_thms} theorems)", file=sys.stderr)


def main():
    import json
    from pass1_parse import parse_pass
    from pass2_ref_elim import ref_elim_pass
    from pass2b_callsite_ref_elim import callsite_ref_elim_pass
    from pass3_lambda_lift import lambda_lift_pass
    from pass3b_lambda_ref_elim import lambda_ref_elim_pass
    from pass4_iter_recognize import iter_recognize_pass
    from pass5_operator_resolve import operator_resolve_pass
    from pass6_ssa_build import ssa_build_pass
    from pass7_loop_lower import loop_lower_pass

    ast_cache = Path(__file__).resolve().parent.parent / "tests" / "ast_cache"
    targets = [k for k in REFINEMENT_MAP if not k.startswith("_loop_")]
    mirs = []
    for fn in targets:
        cache = ast_cache / f"{fn}.json"
        if not cache.exists():
            print(f"  skip: {fn} (no AST cache)", file=sys.stderr)
            continue
        try:
            with open(cache) as f:
                ast = json.load(f)
            h = lambda_ref_elim_pass(lambda_lift_pass(callsite_ref_elim_pass(
                ref_elim_pass(parse_pass(ast)))))
            h3 = iter_recognize_pass(h)
            h4, _ = operator_resolve_pass(h3)
            h4 = callsite_ref_elim_pass(h4, callee_filter={"Rng.next_advance"})
            mir0 = ssa_build_pass(h4)
            mir1 = loop_lower_pass(mir0)
            mirs.append(mir1)
        except Exception as e:
            print(f"  error: {fn} ({type(e).__name__}: {e})", file=sys.stderr)
    refine_gen_pass(mirs)


if __name__ == "__main__":
    main()
