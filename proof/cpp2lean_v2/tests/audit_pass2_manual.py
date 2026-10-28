"""
手动 Pass 2 审核（代替 agent 审核，使用 diff + 结构检查）。

验证：
1. passthrough 函数：HIR₀ 和 HIR₁ 除 gensym 名差异外字节级一致
2. transformed 函数：
   - params 无 [REF] 标记（HIR₁ dump 里不出现 "REF"）
   - 返回类型与预期一致
   - body 所有 return 都被包装或有 trailing return
"""

from __future__ import annotations
import re
import sys
import subprocess
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
HIR0 = Path("/tmp/hir0_dump")
HIR1 = Path("/tmp/hir1_dump")


def normalize_gensym(text: str) -> str:
    """把 `__decomp_<digits>` 的数字标准化为 `<id>`，消除 Python id() 差异。"""
    return re.sub(r"__decomp_\d+", "__decomp_<id>", text)


def diff_passthrough(name: str) -> tuple[bool, str]:
    """比对 HIR₀ 和 HIR₁ 文件（归一 gensym 后）。
    返回 (is_equal, reason)。"""
    f0 = HIR0 / f"{name}.txt"
    f1 = HIR1 / f"{name}.txt"
    if not f0.exists():
        return False, f"missing HIR0: {f0}"
    if not f1.exists():
        return False, f"missing HIR1: {f1}"
    t0 = normalize_gensym(f0.read_text())
    t1 = normalize_gensym(f1.read_text())
    if t0 == t1:
        return True, "byte-identical (after gensym normalization)"
    # 差异情况，返回 diff 摘要
    diff = subprocess.run(
        ["diff", "--brief", "-"],
        input=t0, capture_output=True, text=True,
    )
    # 写临时文件做正式 diff
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as a, \
         tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as b:
        a.write(t0); a.flush()
        b.write(t1); b.flush()
        r = subprocess.run(["diff", a.name, b.name],
                           capture_output=True, text=True)
    Path(a.name).unlink(missing_ok=True)
    Path(b.name).unlink(missing_ok=True)
    return False, r.stdout[:500]


def check_transformed(name: str, expected_refs: list[str],
                      expected_ret: str) -> tuple[bool, list[str]]:
    """验证 transformed 函数。
    - HIR₁ 的 params 段应无 [REF] / [OUTPUT] 标记
    - HIR₁ Return type 应为 expected_ret
    - HIR₁ 末尾应有 return + expected_refs 的 tuple（若原 void）
    """
    f1 = HIR1 / f"{name}.txt"
    if not f1.exists():
        return False, [f"missing HIR1 dump: {f1}"]
    text = f1.read_text()
    errors = []

    # 检查 params 段无 [REF] / [OUTPUT]
    params_section = text.split("## Return type:")[0].split("## Params")[1]
    if "[REF]" in params_section:
        errors.append(f"param 仍有 [REF] 标记")
    if "[OUTPUT]" in params_section:
        errors.append(f"param 仍有 [OUTPUT] 标记")

    # 检查返回类型
    ret_match = re.search(r"## Return type: (.+)", text)
    if not ret_match:
        errors.append("no Return type line")
    else:
        actual_ret = ret_match.group(1).strip()
        if actual_ret != expected_ret:
            errors.append(f"ret 类型不符：expected {expected_ret!r}, got {actual_ret!r}")

    # 统计 return 点数量：每个 return 都该包装成 tuple（含 ref 变量名）
    body_section = text.split("## Body")[1] if "## Body" in text else ""
    return_lines = [l for l in body_section.splitlines() if l.strip().startswith("return")]
    if not return_lines:
        errors.append("body 无 return")
    for rl in return_lines:
        # 每个 return 应含所有 ref 参数名
        for ref_name in expected_refs:
            # 允许 ref_name 在 tuple 内以 "name," 或 "name)" 形式出现
            if ref_name not in rl:
                # 可能在复杂 tuple 里，放宽：检查 `, <ref_name>)` 或 `(<ref_name>`
                if f", {ref_name})" not in rl and f", {ref_name}," not in rl \
                   and f"({ref_name}," not in rl and f"({ref_name})" not in rl:
                    errors.append(f"return 未包装 {ref_name!r}: `{rl.strip()[:120]}`")
                    break

    return len(errors) == 0, errors


# ============================================================
# Passthrough 列表（HIR₁ 应与 HIR₀ 一致）
# ============================================================

PASSTHROUGH = [
    # Zp (9)
    "__make_zp", "__upoly_mod", "__upoly_powmod", "__extract_pth_root",
    "__squarefree_Zp", "__upoly_subtract_x", "__upoly_subtract_one",
    "__ddf_Zp", "__factor_Zp",
    # Univar (21 + factorize 3 实例)
    "__symmetric_mod", "__upoly_symmetric_mod", "__upoly_norm_l2_sq",
    "__binomial", "__isqrt_ceil", "__mignotte_bound",
    "__upoly_mul_mod", "__hensel_tree_build", "__hensel_lift",
    "__heuristic_starting_precision", "__upoly_norm_l1",
    "__upoly_primitive", "__subset_product_mod", "__upoly_const_term",
    "__zassenhaus_recombine", "__cld_polys", "__extract_candidates",
    "__vanhoeij_recombine", "__factor_recombine", "__select_prime",
    "__lll_factorize", "__factor_squarefree_primitive_ZZ",
    "__upoly_to_poly",
    "factorize_upoly", "factorize_lex", "factorize_grlex",
    # Wang (10)
    "__wang_leading_coeff", "__taylor_coeff_zp", "__polynomial_to_zp",
    "__symmetric_mod_poly", "__assign_partial_zp", "__factor_multivar",
    "__select_eval_point", "__mtshl_lift", "__mtshl_coeff_bound",
    "__wang_core",
]


# ============================================================
# Transformed 列表（HIR₁ 应剥离 ref 并包装 return）
# ============================================================

TRANSFORMED = [
    # (name, [ref_names], expected_ret_type_literal)
    # Zp (4)
    ("__upoly_make_monic", ["f"], "Pair<Zp, SparsePolyZp>"),
    ("__upoly_random", ["rng"], "Pair<SparsePolyZp, Rng>"),
    ("__upoly_divmod", ["q", "r"], "Pair<SparsePolyZp, SparsePolyZp>"),
    ("__edf_Zp", ["result", "rng"], "Pair<Array<SparsePolyZp>, Rng>"),
    # Univar (10)
    ("__upoly_mod_coeff", ["f"], "SparsePolyZZ"),
    ("__upoly_divmod_mod", ["q", "r"], "Pair<SparsePolyZZ, SparsePolyZZ>"),
    ("__hensel_tree_build_recursive", ["nodes"], "Array<HenselNode>"),
    ("__hensel_step", ["node"], "HenselNode"),
    ("__hensel_extract_factors", ["factors"], "Array<SparsePolyZZ>"),
    ("__hensel_lift_recursive", ["nodes"], "Array<HenselNode>"),
    ("__hensel_step_linear", ["node"], "HenselNode"),
    ("__hensel_lift_linear_recursive", ["nodes"], "Array<HenselNode>"),
    ("__build_cld_matrix", ["M"], "Pair<Int32, LLLMatrix>"),
    ("__lll_reduce", ["M", "U"], "Tuple<Array<Int32>, LLLMatrix, LLLMatrix>"),
    # Wang (8)
    ("__si_vandermonde_solve", ["coeffs"], "Pair<Bool, Array<Zp>>"),
    ("__si_theta_array_eval", ["images"], "Array<SparsePolyZp>"),
    ("__mtshl_zp_univar_mdp", ["sigma"], "Pair<Bool, Array<SparsePolyZp>>"),
    ("__mtshl_multi_bdp", ["result"], "Pair<Bool, Array<MvPolyZp>>"),
    ("__mtshl_sparse_int", ["result"], "Pair<Bool, Array<MvPolyZp>>"),
    ("__mtshl_wmds", ["result"], "Pair<Bool, Array<MvPolyZp>>"),
    ("__mtshl_step_j", ["F"], "Pair<Bool, Array<MvPolyZp>>"),
    ("__extract_monomial_content", ["var_factors"],
     "Pair<MvPolyZZ, Array<Pair<Variable, Int64>>>"),
]


def main():
    print(f"=== Pass 2 手动审核 ===\n")

    # 1. Passthrough 审核
    print(f"## Passthrough ({len(PASSTHROUGH)} 函数)\n")
    pt_ok, pt_fail = 0, 0
    pt_fails = []
    for name in PASSTHROUGH:
        ok, reason = diff_passthrough(name)
        if ok:
            pt_ok += 1
        else:
            pt_fail += 1
            pt_fails.append((name, reason))
    print(f"OK: {pt_ok} / FAIL: {pt_fail}")
    for name, reason in pt_fails:
        print(f"  FAIL {name}:")
        for line in reason.splitlines()[:5]:
            print(f"    {line}")

    # 2. Transformed 审核
    print(f"\n## Transformed ({len(TRANSFORMED)} 函数)\n")
    tr_ok, tr_fail = 0, 0
    tr_fails = []
    for name, refs, expected_ret in TRANSFORMED:
        ok, errors = check_transformed(name, refs, expected_ret)
        if ok:
            tr_ok += 1
        else:
            tr_fail += 1
            tr_fails.append((name, errors))
    print(f"OK: {tr_ok} / FAIL: {tr_fail}")
    for name, errors in tr_fails:
        print(f"  FAIL {name}:")
        for e in errors:
            print(f"    - {e}")

    total_ok = pt_ok + tr_ok
    total = len(PASSTHROUGH) + len(TRANSFORMED)
    print(f"\n=== 合计 {total_ok}/{total} PASS ===")
    return 0 if (pt_fail + tr_fail) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
