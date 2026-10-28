"""B2B 单条结果对比

接受三个 JSON：cpp_resp / lean_resp / expected (可选)。
返回 (status, detail)：
  status ∈ {"PASS", "FAIL_CPP", "FAIL_LEAN", "FAIL_BOTH",
            "FAIL_DIFF" (cpp 与 lean 不一致), "FAIL_EXPECTED"}
  detail: 字符串描述

per-type NORMALIZERS 钩子：先 normalize 后比较。
"""
from __future__ import annotations

from typing import Any, Callable

# type_name → normalize 函数；normalize 输入是 {"type":..., "val":...} 节点的 val
NORMALIZERS: dict[str, Callable[[Any], Any]] = {
    # 默认无；按需添加，如：
    # "SparsePolyZp": lambda v: sorted(v, key=lambda t: -t[0]),
}


def _normalize_node(node: Any) -> Any:
    """递归遍历 JSON 树，对每个 {"type":..., "val":...} 节点应用 NORMALIZERS。"""
    if isinstance(node, dict):
        if "type" in node and "val" in node and node["type"] in NORMALIZERS:
            return {
                "type": node["type"],
                "val": NORMALIZERS[node["type"]](node["val"]),
            }
        return {k: _normalize_node(v) for k, v in node.items()}
    if isinstance(node, list):
        return [_normalize_node(x) for x in node]
    return node


def _eq(a: Any, b: Any) -> bool:
    return _normalize_node(a) == _normalize_node(b)


def diff_one(cpp_resp: dict, lean_resp: dict, expected: dict | None = None
             ) -> tuple[str, str]:
    """对比单条响应。"""
    cpp_ok = cpp_resp.get("ok", False)
    lean_ok = lean_resp.get("ok", False)

    if not cpp_ok and not lean_ok:
        if expected is not None and "error_contains" in expected:
            needle = expected["error_contains"]
            cpp_err = str(cpp_resp.get("err", ""))
            lean_err = str(lean_resp.get("err", ""))
            if needle in cpp_err and needle in lean_err:
                return ("PASS", "")
            return ("FAIL_EXPECTED",
                    f"expected both errors to contain {needle!r}; "
                    f"cpp_err={cpp_err}; lean_err={lean_err}")
        return ("FAIL_BOTH",
                f"cpp_err={cpp_resp.get('err')}; lean_err={lean_resp.get('err')}")
    if not cpp_ok:
        return ("FAIL_CPP", f"cpp_err={cpp_resp.get('err')}")
    if not lean_ok:
        return ("FAIL_LEAN", f"lean_err={lean_resp.get('err')}")

    cpp_ret = cpp_resp.get("ret")
    lean_ret = lean_resp.get("ret")

    if not _eq(cpp_ret, lean_ret):
        return ("FAIL_DIFF",
                f"cpp_ret={cpp_ret!r} != lean_ret={lean_ret!r}")

    if expected is not None and not _eq(cpp_ret, expected):
        return ("FAIL_EXPECTED",
                f"cpp/lean agree but != expected={expected!r}; got={cpp_ret!r}")

    return ("PASS", "")
