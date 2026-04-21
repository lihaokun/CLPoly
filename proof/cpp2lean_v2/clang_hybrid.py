"""
混合解析器：libclang 定位函数 + Clang JSON dump 翻译单个函数。

libclang 用于：
- 编译文件（检查错误）
- 列出目标文件中的函数
- 获取函数的行号范围

Clang JSON 用于：
- 对单个函数做详细 AST 翻译（使用已有的 clang_ast.py）
  通过 -ast-dump-filter 只 dump 目标函数，避免 773MB 全量 JSON
"""

from __future__ import annotations
import subprocess
import json
from clang.cindex import Index, CursorKind
from ir_types import FuncIR
from clang_ast import parse_func_decl


def _get_system_includes() -> list[str]:
    try:
        result = subprocess.run(
            ["g++", "-E", "-Wp,-v", "-x", "c++", "-"],
            input="", capture_output=True, text=True
        )
        paths = []
        capture = False
        for line in result.stderr.splitlines():
            if "#include <...>" in line:
                capture = True
                continue
            if line.startswith("End of search"):
                break
            if capture and line.startswith(" "):
                paths.append(f"-I{line.strip()}")
        # libclang 的内建头文件（stddef.h 等）
        import glob
        for d in glob.glob("/usr/lib/llvm-*/lib/clang/*/include"):
            paths.append(f"-I{d}")
        return paths
    except Exception:
        return []


SYSTEM_INCLUDES = _get_system_includes()


def list_functions(filename: str, args: list[str] = None,
                   target_file: str = None) -> list[str]:
    """用 libclang 列出目标文件中的函数名。

    target_file 可以是逗号分隔的多个过滤词，如 "factorize_zp,factorize_univar,factorize_wang"。
    """
    if args is None:
        args = ["-std=c++17", "-I./"]

    index = Index.create()
    tu = index.parse(filename, args=args + SYSTEM_INCLUDES)

    errors = [d for d in tu.diagnostics if d.severity >= 3]
    if errors:
        print(f"Warning: {len(errors)} compilation errors")
        for e in errors[:3]:
            print(f"  {e}")

    # 支持多个 target_file 过滤词
    filters = [f.strip() for f in target_file.split(",")] if target_file else []

    funcs = []
    seen = set()
    template_names = set()
    for cursor in tu.cursor.walk_preorder():
        if cursor.kind in (CursorKind.FUNCTION_DECL, CursorKind.FUNCTION_TEMPLATE) \
                and cursor.is_definition():
            loc = cursor.location
            if loc.file:
                if filters and not any(f in loc.file.name for f in filters):
                    continue
                name = cursor.spelling
                if name not in seen:
                    seen.add(name)
                    funcs.append((name, loc.file.name))
                if cursor.kind == CursorKind.FUNCTION_TEMPLATE:
                    template_names.add(name)
    return funcs, template_names


def parse_function_json(filename: str, func_name: str,
                        args: list[str] = None) -> FuncIR | None:
    """用 Clang JSON dump 翻译单个函数（通过 -ast-dump-filter）。"""
    if args is None:
        args = ["-std=c++17", "-I./"]

    cmd = ["clang++"] + args + SYSTEM_INCLUDES + [
        "-Xclang", "-ast-dump=json",
        "-Xclang", f"-ast-dump-filter={func_name}",
        "-fsyntax-only", filename
    ]

    timeout_sec = 30 if "instantiate" in filename else 15
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec)
    except subprocess.TimeoutExpired:
        print(f"Timeout dumping {func_name}")
        return None

    if result.returncode != 0 and not result.stdout.strip():
        print(f"Error dumping {func_name}: {result.stderr[:200]}")
        return None

    try:
        ast = json.loads(result.stdout)
    except json.JSONDecodeError:
        # -ast-dump-filter 子串匹配可能输出多个顶层 JSON 对象
        # 逐个解析，选名字精确匹配的
        decoder = json.JSONDecoder()
        raw = result.stdout.strip()
        pos = 0
        ast = None
        while pos < len(raw):
            try:
                obj, end = decoder.raw_decode(raw[pos:])
                if obj.get("name") == func_name:
                    ast = obj
                    break
                pos += end
                # 跳过空白
                while pos < len(raw) and raw[pos] in ' \t\n\r':
                    pos += 1
            except json.JSONDecodeError:
                break
        if ast is None:
            print(f"JSON parse error for {func_name}")
            return None

    # ast 是单个 FunctionDecl 节点
    if ast.get("kind") == "FunctionDecl":
        return parse_func_decl(ast)

    # 模板函数：FunctionTemplateDecl 包含 FunctionDecl 子节点
    # 优先选有 CompoundStmt body 的实例化版本
    if ast.get("kind") == "FunctionTemplateDecl":
        inner = ast.get("inner", [])
        # 先找实例化（后面的 FunctionDecl 通常是实例化）
        for child in reversed(inner):
            if child.get("kind") == "FunctionDecl":
                has_body = any(c.get("kind") == "CompoundStmt"
                              for c in child.get("inner", []) if isinstance(c, dict))
                if has_body:
                    return parse_func_decl(child)
        # fallback：取第一个 FunctionDecl
        for child in inner:
            if child.get("kind") == "FunctionDecl":
                return parse_func_decl(child)

    return None


def parse_file(filename: str, args: list[str] = None,
               target_file: str = None) -> list[FuncIR]:
    """混合解析：libclang 列函数 + JSON dump 逐个翻译。

    对 instantiate.cc 模式：libclang 发现函数（含模板实例化），
    但 JSON dump 用原始头文件（更快，因为不需要编译全量 STL）。
    """
    import sys as _sys
    func_entries, template_names = list_functions(filename, args, target_file)
    print(f"Found {len(func_entries)} functions ({len(template_names)} templates) in {filename}",
          file=_sys.stderr)

    funcs = []
    for name, source_file in func_entries:
        print(f"  Translating {name}...", end=" ", flush=True, file=_sys.stderr)
        if name in template_names and filename != source_file:
            # 模板函数：用 instantiate.cc dump（含实例化 AST，auto 已展开）
            func = parse_function_json(filename, name, args)
        else:
            # 普通函数：用原始头文件（更快）
            func = parse_function_json(source_file, name, args)
        if func is None and filename != source_file:
            # fallback：尝试另一个源
            other = filename if name not in template_names else source_file
            func = parse_function_json(other, name, args)
        if func:
            funcs.append(func)
            print(f"OK ({len(func.body)} stmts)", file=_sys.stderr)
        else:
            print("FAILED", file=_sys.stderr)

    return funcs
