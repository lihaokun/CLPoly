"""
Pass 6 ssa_build 全量烟测 — 65 函数（factorize 展开 3 实例）。

流程：Pass 1-6；验证 MIR₀ 不变量 + 统计 phi 节点数 + CFG 复杂度。

输出：
- stdout: OK / FAIL + phi 总数 + 平均 BB 数
- `docs/.../survey/pass6-smoke.md`: 详细报告
- `/tmp/mir0_dump/`: 每函数 MIR 文本 dump（基本块 + terminator）
"""

from __future__ import annotations
import json
import sys
import traceback
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, ArrayType, RefType, PairType, StdMapType, TypeIR,
    Var, Lit, BinOp, UnaryOp, CondExpr, UnresolvedOp, Call,
    ArrayAccess, FieldAccess, Cast, LambdaExpr, BlockExpr, TupleExpr, ArrayLit,
    LetStmt, RequireStmt,
    PhiStmt, JumpTerm, CondJumpTerm, ReturnTerm, TailCallTerm,
    BasicBlock, CFG, MIRFunc, HIRFunc, TranslationError,
    assert_mir0_invariant,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass
from pass3_lambda_lift import lambda_lift_pass
from pass4_iter_recognize import iter_recognize_pass
from pass5_operator_resolve import operator_resolve_pass
from pass6_ssa_build import ssa_build_pass
from class_map import TRANSLATION_SCOPE

sys.path.insert(0, str(V2_ROOT / "tests"))
from smoke_pass2_full import _get_factorize_instances

AST_CACHE_DIR = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
OUT_MD = V2_ROOT.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "pass6-smoke.md"
DUMP_DIR = Path("/tmp/mir0_dump")


def _fmt_var(v: Var) -> str:
    return f"{v.name}_{v.version}" if v.version > 0 else v.name


def _fmt_expr(e) -> str:
    if isinstance(e, Var): return _fmt_var(e)
    if isinstance(e, Lit): return str(e.value)
    if isinstance(e, BinOp):
        return f"({_fmt_expr(e.lhs)} {e.op} {_fmt_expr(e.rhs)})"
    if isinstance(e, UnaryOp): return f"({e.op}{_fmt_expr(e.operand)})"
    if isinstance(e, Call):
        cs = e.callee if isinstance(e.callee, str) else f"<op>{e.callee.op_name}"
        return f"{cs}({', '.join(_fmt_expr(a) for a in e.args)})"
    if isinstance(e, Cast):
        return f"cast<{e.target_ty}>({_fmt_expr(e.expr)})"
    if isinstance(e, FieldAccess): return f"{_fmt_expr(e.obj)}.{e.field_name}"
    if isinstance(e, ArrayAccess): return f"{_fmt_expr(e.arr)}[{_fmt_expr(e.idx)}]"
    if isinstance(e, CondExpr):
        return f"({_fmt_expr(e.cond)} ? {_fmt_expr(e.then_e)} : {_fmt_expr(e.else_e)})"
    if isinstance(e, LambdaExpr): return "<lambda>"
    if isinstance(e, TupleExpr):
        return f"({', '.join(_fmt_expr(x) for x in e.elems)})"
    if isinstance(e, ArrayLit):
        return f"#[{', '.join(_fmt_expr(x) for x in e.elems)}]"
    return str(e)[:80]


def _fmt_term(t) -> str:
    if isinstance(t, JumpTerm): return f"goto bb{t.target}"
    if isinstance(t, CondJumpTerm):
        return f"if {_fmt_expr(t.cond)} then bb{t.then_bb} else bb{t.else_bb}"
    if isinstance(t, ReturnTerm):
        return f"return {_fmt_expr(t.value)}" if t.value else "return"
    if isinstance(t, TailCallTerm):
        return f"tailcall {t.target_func}({', '.join(_fmt_expr(a) for a in t.args)})"
    return f"???({type(t).__name__})"


def fmt_mir(func: MIRFunc) -> str:
    lines = [f"# {func.lean_name}", ""]
    lines.append(f"params: {[(p.name, str(p.ty)) for p in func.params]}")
    lines.append(f"ret: {func.ret_ty}")
    if func.cfg is None:
        lines.append("(no cfg)")
        return "\n".join(lines)
    lines.append(f"entry: bb{func.cfg.entry}")
    lines.append(f"blocks: {len(func.cfg.blocks)}")
    for bb_id in sorted(func.cfg.blocks):
        bb = func.cfg.blocks[bb_id]
        preds = func.cfg.preds.get(bb_id, [])
        lines.append(f"\nbb{bb_id}  preds={preds}")
        for s in bb.stmts:
            if isinstance(s, PhiStmt):
                srcs = ", ".join(f"bb{p}: {_fmt_var(v)}" for p, v in s.sources.items())
                lines.append(f"  {_fmt_var(s.target)} := phi({srcs})")
            elif isinstance(s, LetStmt):
                lines.append(f"  let {_fmt_var(s.var)} := {_fmt_expr(s.value)}")
            elif isinstance(s, RequireStmt):
                lines.append(f"  require {s.name}: {_fmt_expr(s.cond)}  [{s.source}]")
            else:
                lines.append(f"  ???({type(s).__name__})")
        lines.append(f"  → {_fmt_term(bb.terminator)}")
    return "\n".join(lines)


def main():
    DUMP_DIR.mkdir(parents=True, exist_ok=True)
    targets = sorted(TRANSLATION_SCOPE)
    print(f"Pass 1+2+3+4+5+6 smoke on {len(targets)} functions...", file=sys.stderr)

    ok = 0; fail = 0; fail_details = []
    total_phi = 0
    total_blocks = 0
    # 残余 silent-bug 计数器（每项越接近 0 越好，0 = clean）
    metrics = {
        "phi_undef_ver0": 0,    # B2 残：phi.sources Var.version=0 引用非-param
        "sideeff_write": 0,     # B7 残：__sideeff_* := __write__(...) 而 root 是 param
        "sideeff_other": 0,     # 其他 __sideeff_*（多为合法 discarded return）
        "aux_dropped": 0,       # B6 残：HIR.aux_lambdas 数 != MIR.aux_defs 数
        "rangefor_no_writeback": 0,  # P0-1：ranged-for 有 cont latch __write__ 但无 exit 写回
        "lambda_ref_callsite_unwritten": 0,  # P0(lambda)：含 ref 的 lifted lambda 调用点未 destructure
    }

    def _walk_mir(mir, hir_aux_count, hir_aux_with_ref=None):
        from ir_types import LetStmt, Call, Var
        import re
        n_phi = 0
        param_names = {p.name for p in mir.params}
        if len(mir.aux_defs) != hir_aux_count:
            metrics["aux_dropped"] += 1
        # 收集本函数所有 lifted lambda 名字（用于检测调用点是否 destructure）
        lambda_names = {aux.base_name for aux in mir.aux_defs}
        # 含 modified ref 的 lambda（调用点必须 destructure）
        # 由于 Pass 3b 后 aux_defs 已经 ref-elim 过，看 ret_ty 不再有 is_ref；
        # 改用 hir_aux_with_ref（外部传入的修复前快照）作为参考。若无快照，跳过此检测。
        ref_lambda_names = hir_aux_with_ref or set()
        # P0-1：每个 cont 实例（按 N 分组）应有 latch `let cont_N_X := __write__(cont_N[idx],...)`
        # 同时 exit 应有写回（一个 LetStmt 把 cont_N_M 赋给某个变量/字段）。
        # 检测：某个 cont_N 在 latch 出现 __write__ 写回但全 dump 没有任何 LetStmt
        # 把 cont_N_* 作为 RHS 给到非 cont 变量。
        cont_writes_in_latch: dict[str, bool] = {}      # "cont_4" -> True if latch __write__
        cont_used_as_rhs: dict[str, bool] = {}          # "cont_4" -> True if used as rhs to non-cont LetStmt
        for bb in mir.cfg.blocks.values():
            for s in bb.stmts:
                if isinstance(s, PhiStmt):
                    n_phi += 1
                    for src in s.sources.values():
                        if isinstance(src, Var) and src.version == 0 \
                                and src.name not in param_names:
                            metrics["phi_undef_ver0"] += 1
                elif isinstance(s, LetStmt) and s.var.name.startswith("__sideeff_"):
                    if isinstance(s.value, Call) and s.value.callee == "__write__":
                        metrics["sideeff_write"] += 1
                    else:
                        metrics["sideeff_other"] += 1
                        # P0 lambda：含 ref 的 lifted lambda 被 sideeff（未 destructure）
                        if isinstance(s.value, Call) and isinstance(s.value.callee, str):
                            if s.value.callee in ref_lambda_names:
                                metrics["lambda_ref_callsite_unwritten"] += 1
                # P0-1 cont latch __write__ 检测
                if isinstance(s, LetStmt):
                    m = re.match(r'__rangefor_cont_(\d+)_\d+', s.var.name)
                    if m and isinstance(s.value, Call) and s.value.callee == "__write__":
                        cont_writes_in_latch[m.group(1)] = True
                    # cont_*_* 作为 RHS 被赋给非 cont 变量（exit 写回）
                    if not s.var.name.startswith("__rangefor_cont_"):
                        if isinstance(s.value, Var):
                            mr = re.match(r'__rangefor_cont_(\d+)_\d+', s.value.name)
                            if mr:
                                cont_used_as_rhs[mr.group(1)] = True
                        elif isinstance(s.value, Call) and s.value.callee == "__write__" \
                                and len(s.value.args) >= 2 and isinstance(s.value.args[1], Var):
                            mr = re.match(r'__rangefor_cont_(\d+)_\d+', s.value.args[1].name)
                            if mr:
                                cont_used_as_rhs[mr.group(1)] = True
        # 凡 latch 写但 exit 没回写的 cont = 漏写回
        for cn, has_latch in cont_writes_in_latch.items():
            if has_latch and not cont_used_as_rhs.get(cn, False):
                metrics["rangefor_no_writeback"] += 1
        # 递归 aux_defs
        for aux in mir.aux_defs:
            _walk_mir(aux, len(aux.aux_defs))  # aux 没有 HIR 比对源
        return n_phi

    def process(tag, ast, dump_name):
        nonlocal ok, fail, total_phi, total_blocks
        try:
            from pass2b_callsite_ref_elim import callsite_ref_elim_pass
            from pass3b_lambda_ref_elim import lambda_ref_elim_pass
            # 修复前快照：在 Pass 3b 之前记录哪些 lifted lambda 含 ref param
            hir_pre_3b = lambda_lift_pass(callsite_ref_elim_pass(ref_elim_pass(parse_pass(ast))))
            ref_lambda_names = {aux.base_name for aux in hir_pre_3b.aux_lambdas
                                if any(p.is_ref for p in aux.params)}
            hir3 = iter_recognize_pass(lambda_ref_elim_pass(hir_pre_3b))
            hir4, _ = operator_resolve_pass(hir3)
            mir = ssa_build_pass(hir4)
            assert_mir0_invariant(mir)
            (DUMP_DIR / f"{dump_name}.txt").write_text(fmt_mir(mir))
            n_phi = _walk_mir(mir, len(hir4.aux_lambdas), ref_lambda_names)
            total_phi += n_phi
            total_blocks += len(mir.cfg.blocks)
            ok += 1
        except TranslationError as e:
            fail += 1
            fail_details.append((tag, f"TranslationError: {e.reason}"))
        except Exception as e:
            fail += 1
            tb = traceback.format_exc()[:300]
            fail_details.append((tag, f"{type(e).__name__}: {str(e)[:150]}"))

    for fn in targets:
        if fn == "factorize":
            instances = _get_factorize_instances()
            for suffix, ast in instances:
                process(f"factorize_{suffix}", ast, f"factorize_{suffix}")
            continue
        cache = AST_CACHE_DIR / f"{fn}.json"
        if not cache.exists():
            fail += 1; fail_details.append((fn, "no AST cache")); continue
        with open(cache) as f: ast = json.load(f)
        process(fn, ast, fn)

    n_total = ok + fail
    avg_phi = total_phi / max(ok, 1)
    avg_blocks = total_blocks / max(ok, 1)

    lines = []
    lines.append("# Pass 1-6 全量烟测")
    lines.append("")
    lines.append(f"- 目标 HIRs：**{n_total}**（factorize 展开 3 实例）")
    lines.append(f"- OK: **{ok}** / FAIL: **{fail}**")
    lines.append(f"- 平均 phi 节点数: **{avg_phi:.1f}** / 函数")
    lines.append(f"- 平均 BasicBlock 数: **{avg_blocks:.1f}** / 函数")
    lines.append(f"- 总 phi: {total_phi} / 总 BB: {total_blocks}")
    lines.append("")
    lines.append("## Silent-bug 残余指标")
    lines.append("")
    lines.append(f"- B2 phi.sources ver=0 非-param: **{metrics['phi_undef_ver0']}**")
    lines.append(f"- B6 aux_dropped (HIR.aux_lambdas != MIR.aux_defs): **{metrics['aux_dropped']}**")
    lines.append(f"- B7 __sideeff_*=__write__(...) (root 非 def): **{metrics['sideeff_write']}**")
    lines.append(f"- P0-1 ranged-for 漏写回原 container: **{metrics['rangefor_no_writeback']}**")
    lines.append(f"- P0(lambda) 含 ref 的 lifted lambda 调用点未 destructure: **{metrics['lambda_ref_callsite_unwritten']}**")
    lines.append(f"- 其他 __sideeff_* (多为 discarded return): {metrics['sideeff_other']}")
    lines.append("")
    if fail_details:
        lines.append("## FAIL")
        lines.append("")
        for fn, r in fail_details:
            lines.append(f"- `{fn}`: {r}")
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines))

    print(f"\n=== {ok} OK / {fail} FAIL ===", file=sys.stderr)
    print(f"  total phi: {total_phi} (avg {avg_phi:.1f}/func)", file=sys.stderr)
    print(f"  total BB:  {total_blocks} (avg {avg_blocks:.1f}/func)", file=sys.stderr)
    print(f"  silent-bug residue:", file=sys.stderr)
    for k, v in metrics.items():
        print(f"    {k}: {v}", file=sys.stderr)
    print(f"\nReport: {OUT_MD}", file=sys.stderr)
    print(f"MIR₀ dumps: {DUMP_DIR}", file=sys.stderr)
    # 关键：B2/B6/P0-1 残余非 0 → smoke 失败（B7 sideeff_write 接受少量残余）
    if metrics["aux_dropped"] > 0:
        print(f"FAIL: aux_dropped > 0 — B6 regression", file=sys.stderr)
        return 1
    if metrics["rangefor_no_writeback"] > 0:
        print(f"FAIL: rangefor_no_writeback > 0 — P0-1 regression", file=sys.stderr)
        return 1
    if metrics["lambda_ref_callsite_unwritten"] > 0:
        print(f"FAIL: lambda_ref_callsite_unwritten > 0 — P0(lambda) regression", file=sys.stderr)
        return 1
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
