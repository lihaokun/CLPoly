"""
Pass 2 ref_elim 全量烟测 — 65 函数。

流程：Pass 1 parse → Pass 2 ref_elim；验证 HIR₁ 不变量 + 统计转换效果。

输出：
- stdout: OK / WARN / FAIL 统计
- `docs/.../survey/pass2-smoke.md`: 详细报告
- `/tmp/hir1_dump/`: 每函数的 HIR₁ dump（供人工审核）

可以随时重跑作为回归测试。
"""

from __future__ import annotations
import json
import sys
import traceback
from pathlib import Path
from collections import Counter, defaultdict

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))
sys.path.insert(0, str(V2_ROOT / "passes"))

from ir_types import (
    BaseType, NamedType, PairType, TupleType, RefType, UnknownType,
    HIRFunc, HIRParam, TranslationError,
)
from pass1_parse import parse_pass
from pass2_ref_elim import ref_elim_pass, assert_hir1_invariant
from class_map import TRANSLATION_SCOPE

sys.path.insert(0, str(V2_ROOT / "tests"))
from dump_all_hir import fmt_func

AST_CACHE_DIR = V2_ROOT.parent / "cpp2lean" / "_ast_cache"
OUT_MD = V2_ROOT.parent.parent / "docs" / "design" / "l1-translation-validation" / "survey" / "pass2-smoke.md"
DUMP_DIR = Path("/tmp/hir1_dump")


def fmt_type(ty) -> str:
    """复用 dump_all_hir 的类型格式化。"""
    from dump_all_hir import fmt_type as _fmt
    return _fmt(ty)


def _get_factorize_instances() -> list[tuple[str, dict]]:
    """抓 factorize 全部 3 实例的 AST，返回 [(suffix, ast_dict), ...]。

    复用 dump_all_hir._dump_factorize_instances 的 clang 调用模式，
    但直接返回 ast 供 Pass 2 处理（而非立即 dump）。
    """
    import subprocess
    import glob
    PROJECT_ROOT = V2_ROOT.parent.parent

    def sys_incs() -> list[str]:
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

    def suffix_of(qt: str) -> str:
        if "upolynomial_<ZZ>" in qt: return "upoly"
        if "grlex_<less>" in qt: return "grlex"
        if "lex_<less>" in qt: return "lex"
        if "QQ," in qt: return "qq"
        return "unknown"

    instances: list[tuple[str, dict]] = []
    for cand in candidates:
        if cand.get("name") == "factorize":
            for inst in expand(cand):
                if inst.get("mangledName") and has_body(inst):
                    qt = inst.get("type", {}).get("qualType", "")
                    instances.append((suffix_of(qt), inst))
    return instances


def main():
    DUMP_DIR.mkdir(parents=True, exist_ok=True)

    targets = sorted(TRANSLATION_SCOPE)
    print(f"Pass 1+2 smoke on {len(targets)} functions...", file=sys.stderr)

    ok_count = 0
    transformed_count = 0  # 有 ref 参数被消除的函数数
    passthrough_count = 0   # 无 ref 参数的函数数
    fail_count = 0
    fail_details = []

    # 统计每个函数的转换前/后
    per_func = {}  # func_name → dict

    def process_ast(tag: str, ast: dict, dump_name: str):
        """对单个 AST 跑 Pass 1+2，更新统计 + dump HIR₁。
        tag: 报告用的显示名（如 "__upoly_divmod" 或 "factorize_upoly"）
        dump_name: 输出文件名（不含 .txt）"""
        nonlocal ok_count, transformed_count, passthrough_count, fail_count
        try:
            hir0 = parse_pass(ast)
            pass1_ref_params = [p for p in hir0.params if p.is_ref]
            pass1_ref_count = len(pass1_ref_params)
            pass1_ret_ty = hir0.ret_ty

            hir1 = ref_elim_pass(hir0)
            assert_hir1_invariant(hir1)

            pass2_has_ref = any(p.is_ref for p in hir1.params)
            assert not pass2_has_ref, "should have no is_ref after Pass 2"

            if pass1_ref_count > 0:
                transformed_count += 1
            else:
                passthrough_count += 1
            ok_count += 1

            per_func[tag] = {
                "status": "OK",
                "ref_count": pass1_ref_count,
                "ref_names": [p.name for p in pass1_ref_params],
                "pass1_ret_ty": str(pass1_ret_ty),
                "pass2_ret_ty": str(hir1.ret_ty),
                "n_params_before": len(hir0.params),
                "n_params_after": len(hir1.params),
            }
            (DUMP_DIR / f"{dump_name}.txt").write_text(fmt_func(hir1))

        except TranslationError as e:
            fail_count += 1
            fail_details.append((tag, f"TranslationError: {e.reason}"))
            per_func[tag] = {"status": "FAIL", "reason": e.reason}
        except AssertionError as e:
            fail_count += 1
            fail_details.append((tag, f"AssertionError: {e}"))
            per_func[tag] = {"status": "FAIL", "reason": str(e)}
        except Exception as e:
            fail_count += 1
            tb = traceback.format_exc()[:500]
            fail_details.append((tag, f"{type(e).__name__}: {e}"))
            per_func[tag] = {"status": "FAIL", "reason": f"{type(e).__name__}: {e}", "traceback": tb}

    for func_name in targets:
        # factorize 特殊处理：抓全 3 实例
        if func_name == "factorize":
            instances = _get_factorize_instances()
            if not instances:
                fail_count += 1
                fail_details.append((func_name, "clang returned 0 factorize instances"))
                per_func[func_name] = {"status": "FAIL", "reason": "0 instances"}
                continue
            for suffix, ast in instances:
                tag = f"factorize_{suffix}"
                process_ast(tag, ast, tag)
            # 同时保留 factorize.txt 作兼容（用 grlex 或 lex 实例）
            for suffix, ast in instances:
                if suffix in ("grlex", "lex"):
                    # 再 dump 一次作为 factorize.txt 兼容
                    try:
                        hir0 = parse_pass(ast)
                        hir1 = ref_elim_pass(hir0)
                        (DUMP_DIR / "factorize.txt").write_text(fmt_func(hir1))
                    except Exception:
                        pass
                    break
            continue

        cache_file = AST_CACHE_DIR / f"{func_name}.json"
        if not cache_file.exists():
            fail_count += 1
            fail_details.append((func_name, "no AST cache"))
            continue
        with open(cache_file) as f:
            ast = json.load(f)
        process_ast(func_name, ast, func_name)

    # ======== 输出报告 ========
    lines = []
    lines.append("# Pass 1 + Pass 2 全量烟测（65 函数）")
    lines.append("")
    lines.append(f"- 总函数数：**{len(targets)}**")
    lines.append(f"- OK: **{ok_count}**")
    lines.append(f"  - 有 ref 参数被消除: **{transformed_count}**")
    lines.append(f"  - 无 ref 参数（passthrough）: **{passthrough_count}**")
    lines.append(f"- FAIL: **{fail_count}**")
    lines.append("")

    if fail_details:
        lines.append("## FAIL 列表")
        lines.append("")
        for fn, reason in fail_details:
            lines.append(f"- `{fn}`: {reason}")
        lines.append("")

    # 按 ref_count 分布
    lines.append("## 转换后函数详情（按 ref 参数数分组）")
    lines.append("")
    by_refcount = defaultdict(list)
    for fn, d in per_func.items():
        if d.get("status") == "OK":
            by_refcount[d["ref_count"]].append(fn)
    for rc in sorted(by_refcount.keys()):
        funcs = by_refcount[rc]
        lines.append(f"### {rc} ref 参数（{len(funcs)} 函数）")
        lines.append("")
        if rc == 0:
            # 简短列出
            for fn in funcs:
                lines.append(f"- `{fn}`")
        else:
            lines.append("| 函数 | ref 参数名 | Pass 1 返回类型 → Pass 2 返回类型 |")
            lines.append("|---|---|---|")
            for fn in sorted(funcs):
                d = per_func[fn]
                refs = ", ".join(d["ref_names"])
                lines.append(f"| `{fn}` | {refs} | `{d['pass1_ret_ty']}` → `{d['pass2_ret_ty']}` |")
        lines.append("")

    # 核验：与 Week 2 Day 1 scan_ref_params.md 的预期值对照
    lines.append("## 与 `ref_params.md` 预期值核对")
    lines.append("")
    total_ref_params = sum(d["ref_count"] for d in per_func.values()
                           if d.get("status") == "OK")
    lines.append(f"- 本轮测得 ref 参数总数: **{total_ref_params}**")
    lines.append(f"- `ref_params.md` 预期（Week 2 Day 1）: 26")
    lines.append("")
    if total_ref_params == 26:
        lines.append("✅ 与预期值完全一致")
    else:
        lines.append(f"⚠️ 不一致，差异 {abs(total_ref_params - 26)}")
    lines.append("")

    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(lines))

    # stdout 总结
    print(f"\n=== {ok_count} OK / {fail_count} FAIL ===", file=sys.stderr)
    print(f"  transformed (had ref): {transformed_count}", file=sys.stderr)
    print(f"  passthrough (no ref):  {passthrough_count}", file=sys.stderr)
    print(f"  total ref params: {total_ref_params} (expected: 26)", file=sys.stderr)
    print(f"\nReport: {OUT_MD}", file=sys.stderr)
    print(f"HIR₁ dumps: {DUMP_DIR}", file=sys.stderr)

    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
