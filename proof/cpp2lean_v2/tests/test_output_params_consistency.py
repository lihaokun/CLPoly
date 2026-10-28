"""
test_output_params_consistency.py — 机械校验 TRANSLATION_SCOPE_OUTPUT_PARAMS
注册的 ref-out 索引与 C++ 源签名一致。

Why this exists:
    Phase 6 第七轮审视发现 4 个 __mtshl_* 函数注册索引全错（如 [3, 8] 撞 6-arg
    函数），且因 Pass 2b silent skip 越界 + 之前不处理 LetStmt RHS 而潜伏 6+
    轮。手写注册表是"无类型外部世界"，必须机械化校验。

How it works:
    1) 扫 clpoly/*.hh / *.cc 找每个注册函数的定义
    2) 解析 params（处理嵌套模板 <...>）→ 分类 const-ref / non-const-ref / 值
    3) 把 non-const-ref 索引列表与注册值对比
    4) 报告 mismatch / 缺源 / OK

Limitation:
    - 仅查 clpoly/ 目录的 source；STL/GMP 函数（std::swap、fdiv_q 等）跳过
      （标 ⚠ external，不视为错）
    - 不解析 default args（罕见且不影响 ref 类别）
    - 模板特化 / 多 overload：用 #N arity suffix 区分（class_map.py 已用此约定）
"""

from __future__ import annotations
import re
import sys
from pathlib import Path

V2_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(V2_ROOT))

from class_map import TRANSLATION_SCOPE_OUTPUT_PARAMS

CLPOLY_SRC = V2_ROOT.parent.parent / "clpoly"


def _split_top_level_commas(s: str) -> list[str]:
    """按顶层 `,` 切分（跳过 `<...>` `(...)` 内部）。"""
    parts: list[str] = []
    depth = 0
    cur: list[str] = []
    for ch in s:
        if ch in "<(":
            depth += 1
            cur.append(ch)
        elif ch in ">)":
            depth -= 1
            cur.append(ch)
        elif ch == "," and depth == 0:
            parts.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    if cur:
        parts.append("".join(cur).strip())
    return [p for p in parts if p]


def _classify_param(p: str) -> str:
    """分类 param 字符串：'const-ref' / 'ref-out' / 'value'。"""
    # 简化：去 default arg
    if "=" in p:
        p = p.split("=")[0].strip()
    has_ref = "&" in p
    # 检查 const：仅当 const 在 `&` 之前才算 const-ref（const T& vs T const&）
    # 简单规则：含 "const " 且含 "&" → const-ref（容忍 `T& const` 极其罕见）
    is_const = bool(re.search(r"\bconst\b", p))
    if has_ref and not is_const:
        return "ref-out"
    if has_ref and is_const:
        return "const-ref"
    return "value"


_TYPE_KEYWORDS = re.compile(
    r"\b(void|bool|int|int8_t|int16_t|int32_t|int64_t|uint8_t|uint16_t|uint32_t|"
    r"uint64_t|size_t|char|short|long|float|double|auto|signed|unsigned|"
    r"static|inline|template|polynomial_|upolynomial_|std::|ZZ|QQ|Zp|"
    r"variable|umonomial|basic_monomial|HenselNode|HenselTree|"
    r"PrimeSelectionResult|WangLcResult|MultivariateFactor|"
    r"Factorization|polynomial)\b"
)
_INVALID_HEAD_CHARS = set(";{}")  # 不含 `,` —— template `<class A, class B>` 与注释中
                                  # 出现的 `,` 是合法的；只用 `;` `{` `}` 区分语句边界


def _find_function_defs(name: str) -> list[tuple[Path, int, list[str]]]:
    """在 CLPOLY_SRC 下找 name 的**定义**（含 body `{...}`，排除调用点）。

    判定定义：
      1) name 前的 line head 段含**返回类型/修饰关键字**（void/bool/int/template/...）
      2) line head 段不含 `;` `,` `{` `}` 等终结字符
      3) 匹配的 `)` 之后下一个非空字符必须是 `{`（必须有 body）
    """
    defs: list[tuple[Path, int, list[str]]] = []
    name_re = re.compile(r"\b" + re.escape(name) + r"\s*\(")
    for f in list(CLPOLY_SRC.glob("*.hh")) + list(CLPOLY_SRC.glob("*.cc")):
        raw_text = f.read_text()
        # 预剥行注释 //... 到行尾（保留换行以保持 line_no 一致）；
        # 否则 `// factors[start..end)` 注释中的 `)` 会干扰 paren counter
        text = re.sub(r"//[^\n]*", "", raw_text)
        for m in name_re.finditer(text):
            paren_open = m.end() - 1  # `(` 位置
            # 1) 找 line head（向前找 `\n`，跨多行模板支持）
            #    简化：取 name 之前直至最近的 `;` 或 `}` 或文件头
            line_head_start = max(
                text.rfind(";", 0, m.start()),
                text.rfind("}", 0, m.start()),
            ) + 1
            head = text[line_head_start:m.start()]
            # 剥行注释 // ... 到行尾（避免注释中 `;` `{` `}` 干扰）
            head_clean = re.sub(r"//[^\n]*", "", head)
            if not _TYPE_KEYWORDS.search(head_clean):
                continue
            if any(c in head_clean for c in _INVALID_HEAD_CHARS):
                continue

            # 2) 匹配 `)`
            depth = 1
            i = paren_open + 1
            while i < len(text) and depth > 0:
                ch = text[i]
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                i += 1
            if depth != 0:
                continue
            params_str = text[paren_open + 1:i - 1].strip()

            # 3) `)` 后下一个非空字符必须是 `{`（有 body）
            tail = text[i:i + 200].lstrip()
            if not tail.startswith("{"):
                continue

            params = _split_top_level_commas(params_str)
            line_no = text[:m.start()].count("\n") + 1
            defs.append((f, line_no, params))
    return defs


def check_one(key: str, registered: list[int]) -> tuple[str, str]:
    """返回 (status, detail) — status ∈ {OK, MISMATCH, EXTERNAL, NOT_FOUND}。"""
    if "#" in key:
        name, arity_s = key.split("#", 1)
        arity = int(arity_s)
    else:
        name = key
        arity = None

    defs = _find_function_defs(name)
    if not defs:
        return "EXTERNAL", f"<no definition in clpoly/>: {key}"

    # 按 arity 过滤
    if arity is not None:
        defs_arity = [d for d in defs if len(d[2]) == arity]
        if not defs_arity:
            # clpoly 有同名但不同 arity → 多半是 STL 版本（std::swap 等）
            return "EXTERNAL", (f"{key}: clpoly 仅有 arity {[len(d[2]) for d in defs]} "
                                 f"重载（无 arity={arity}）；目标实为 STL/外部")
        defs = defs_arity

    # 用第一个匹配；多重载时若 arity 都同则只校第一个（同 arity 通常是
    # 模板/特化重复，参数类别相同）
    f, line, params = defs[0]
    actual_ref_out: list[int] = []
    classes: list[str] = []
    for i, p in enumerate(params):
        c = _classify_param(p)
        classes.append(c)
        if c == "ref-out":
            actual_ref_out.append(i)

    if actual_ref_out == sorted(registered):
        return "OK", f"{key}: ref-out at {actual_ref_out} ✓"
    return "MISMATCH", (
        f"{key}: registered={sorted(registered)}, actual={actual_ref_out}; "
        f"def at {f.name}:{line}, params: "
        + ", ".join(f"{i}={c}" for i, c in enumerate(classes))
    )


def test_all_registrations():
    """校验所有 TRANSLATION_SCOPE_OUTPUT_PARAMS 条目与 C++ 签名一致。"""
    mismatches: list[str] = []
    externals: list[str] = []
    okays: list[str] = []
    for key, registered in TRANSLATION_SCOPE_OUTPUT_PARAMS.items():
        status, detail = check_one(key, registered)
        if status == "OK":
            okays.append(detail)
        elif status == "EXTERNAL":
            externals.append(detail)
        else:  # MISMATCH or NOT_FOUND
            mismatches.append(detail)

    print(f"\n=== TRANSLATION_SCOPE_OUTPUT_PARAMS 校验 ===")
    print(f"OK: {len(okays)}")
    for s in okays:
        print(f"  ✓ {s}")
    if externals:
        print(f"\nEXTERNAL（clpoly/ 外，跳过）: {len(externals)}")
        for s in externals:
            print(f"  ⚠ {s}")
    if mismatches:
        print(f"\nMISMATCH: {len(mismatches)}")
        for s in mismatches:
            print(f"  ✗ {s}")
        assert False, f"{len(mismatches)} 个注册与 C++ 签名不一致"
    print("\n全部一致（external 除外）。")


if __name__ == "__main__":
    try:
        test_all_registrations()
    except AssertionError as e:
        print(f"\nFAIL: {e}", file=sys.stderr)
        sys.exit(1)
    print("PASS")
