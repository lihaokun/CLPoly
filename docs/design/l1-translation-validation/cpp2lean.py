#!/usr/bin/env python3
"""
CLPoly C++ → Lean IR 翻译器原型

用法: clang++ -Xclang -ast-dump=json -fsyntax-only file.cpp | python3 cpp2lean.py

输出: Lean 4 def + UB 证明目标注释
"""

import json
import sys

class LeanTranslator:
    def __init__(self):
        self.indent = 0
        self.ub_obligations = []  # UB 证明目标收集
        self.var_types = {}  # 变量 → 类型映射

    def emit(self, s):
        return "  " * self.indent + s

    def translate_type(self, qual_type):
        """C++ 类型 → Lean 类型"""
        t = qual_type.replace("const ", "").strip()
        mapping = {
            "uint64_t": "Nat",
            "uint32_t": "Nat",
            "int": "Int",
            "int64_t": "Int",
            "unsigned __int128": "Nat",  # 建模为 Nat（无溢出）
            "void": "Unit",
        }
        for cpp_t, lean_t in mapping.items():
            if cpp_t in t:
                return lean_t
        return f"/{t}/"  # 未知类型标记

    def translate_expr(self, node):
        """AST 表达式 → Lean 表达式"""
        kind = node.get("kind", "")

        if kind == "IntegerLiteral":
            return str(node.get("value", "0"))

        if kind == "DeclRefExpr":
            ref = node.get("referencedDecl", {})
            return ref.get("name", "?")

        if kind == "BinaryOperator":
            op = node.get("opcode", "?")
            inner = node.get("inner", [])
            if len(inner) >= 2:
                lhs = self.translate_expr(inner[0])
                rhs = self.translate_expr(inner[1])
                lean_op = {
                    "+": "+", "-": "-", "*": "*",
                    "/": "/", "%": "%",
                    "<<": "* 2^", ">>": "/ 2^",
                    "<": "<", ">": ">", ">=": ">=", "<=": "<=",
                    "==": "==", "!=": "!=",
                }.get(op, f" /{op}/ ")

                # UB 检查
                if op == "/" and rhs != "0":
                    self.ub_obligations.append(f"  -- UB-check: {rhs} ≠ 0")

                return f"({lhs} {lean_op} {rhs})"
            return f"?binop_{op}"

        if kind == "UnaryOperator":
            op = node.get("opcode", "")
            inner = node.get("inner", [])
            if inner:
                operand = self.translate_expr(inner[0])
                if op == "!":
                    return f"(¬{operand})"
                return f"({op}{operand})"

        if kind == "ConditionalOperator":
            inner = node.get("inner", [])
            if len(inner) >= 3:
                cond = self.translate_expr(inner[0])
                then_e = self.translate_expr(inner[1])
                else_e = self.translate_expr(inner[2])
                return f"(if {cond} then {then_e} else {else_e})"

        if kind == "CStyleCastExpr" or kind == "ImplicitCastExpr":
            inner = node.get("inner", [])
            if inner:
                return self.translate_expr(inner[0])

        if kind == "ParenExpr":
            inner = node.get("inner", [])
            if inner:
                return self.translate_expr(inner[0])

        if kind == "CallExpr":
            inner = node.get("inner", [])
            if inner:
                fn_name = self.translate_expr(inner[0])
                args = [self.translate_expr(a) for a in inner[1:]]
                return f"({fn_name} {' '.join(args)})"

        if kind == "ArraySubscriptExpr":
            inner = node.get("inner", [])
            if len(inner) >= 2:
                arr = self.translate_expr(inner[0])
                idx = self.translate_expr(inner[1])
                self.ub_obligations.append(f"  -- UB-check: {idx} < {arr}.size")
                return f"({arr}[{idx}])"

        if kind == "MemberExpr":
            inner = node.get("inner", [])
            name = node.get("name", "?")
            if inner:
                base = self.translate_expr(inner[0])
                return f"{base}.{name}"

        if kind == "CompoundAssignOperator":
            op = node.get("opcode", "?")
            inner = node.get("inner", [])
            if len(inner) >= 2:
                lhs = self.translate_expr(inner[0])
                rhs = self.translate_expr(inner[1])
                base_op = op.replace("=", "")
                return f"({lhs} {base_op} {rhs})"

        # 默认：返回 kind 作为占位
        inner = node.get("inner", [])
        if inner:
            return self.translate_expr(inner[0])
        return f"?{kind}"

    def translate_stmt(self, node):
        """AST 语句 → Lean let/if/match"""
        kind = node.get("kind", "")
        lines = []

        if kind == "DeclStmt":
            inner = node.get("inner", [])
            for decl in inner:
                if decl.get("kind") == "VarDecl":
                    name = decl.get("name", "?")
                    typ = self.translate_type(decl.get("type", {}).get("qualType", "?"))
                    self.var_types[name] = typ
                    init = decl.get("inner", [])
                    if init:
                        val = self.translate_expr(init[0])
                        lines.append(self.emit(f"let {name} : {typ} := {val}"))
                    else:
                        lines.append(self.emit(f"let {name} : {typ} := 0"))

        elif kind == "ReturnStmt":
            inner = node.get("inner", [])
            if inner:
                val = self.translate_expr(inner[0])
                lines.append(self.emit(f"return {val}"))

        elif kind == "IfStmt":
            inner = node.get("inner", [])
            if len(inner) >= 2:
                cond = self.translate_expr(inner[0])
                lines.append(self.emit(f"if {cond} then"))
                self.indent += 1
                then_lines = self.translate_stmt(inner[1])
                lines.extend(then_lines)
                self.indent -= 1
                if len(inner) >= 3:
                    lines.append(self.emit(f"else"))
                    self.indent += 1
                    else_lines = self.translate_stmt(inner[2])
                    lines.extend(else_lines)
                    self.indent -= 1

        elif kind == "CompoundStmt":
            inner = node.get("inner", [])
            for stmt in inner:
                lines.extend(self.translate_stmt(stmt))

        elif kind == "BinaryOperator" and node.get("opcode") == "=":
            inner = node.get("inner", [])
            if len(inner) >= 2:
                lhs = self.translate_expr(inner[0])
                rhs = self.translate_expr(inner[1])
                lines.append(self.emit(f"let {lhs} := {rhs}"))

        elif kind == "ForStmt":
            lines.append(self.emit("-- FOR loop (requires manual translation to Nat.fold)"))
            inner = node.get("inner", [])
            if inner:
                body = inner[-1] if inner else None
                if body:
                    lines.append(self.emit("sorry -- loop body:"))

        elif kind == "CompoundAssignOperator":
            inner = node.get("inner", [])
            op = node.get("opcode", "?")
            if len(inner) >= 2:
                lhs = self.translate_expr(inner[0])
                rhs = self.translate_expr(inner[1])
                base_op = op.replace("=", "")
                lines.append(self.emit(f"let {lhs} := ({lhs} {base_op} {rhs})"))

        else:
            # 尝试作为表达式处理
            expr = self.translate_expr(node)
            if expr and not expr.startswith("?"):
                lines.append(self.emit(f"let _ := {expr}"))

        return lines

    def translate_function(self, node):
        """FunctionDecl → Lean def"""
        name = node.get("name", "?")
        ret_type = self.translate_type(node.get("type", {}).get("qualType", "").split("(")[0])

        # 提取参数
        params = []
        body = None
        for child in node.get("inner", []):
            if child.get("kind") == "ParmVarDecl":
                pname = child.get("name", "?")
                ptype = self.translate_type(child.get("type", {}).get("qualType", "?"))
                params.append(f"({pname} : {ptype})")
                self.var_types[pname] = ptype
            elif child.get("kind") == "CompoundStmt":
                body = child

        lines = []
        lines.append(f"/-- C++ {name} → Lean IR (auto-translated) -/")
        lines.append(f"def {name}_ir {' '.join(params)} : {ret_type} :=")

        self.indent = 1
        self.ub_obligations = []

        if body:
            body_lines = self.translate_stmt(body)
            lines.extend(body_lines)

        # 附加 UB 证明目标
        if self.ub_obligations:
            lines.append("")
            lines.append(f"/- UB-freedom proof obligations for {name}:")
            for ob in self.ub_obligations:
                lines.append(ob)
            lines.append("-/")

        return "\n".join(lines)


def main():
    ast = json.load(sys.stdin)
    translator = LeanTranslator()

    print("-- Auto-generated Lean IR from C++ via Clang AST")
    print("-- Generated by cpp2lean.py")
    print("")

    for node in ast.get("inner", []):
        if node.get("kind") == "FunctionDecl" and node.get("inner"):
            # 只翻译有函数体的声明
            has_body = any(c.get("kind") == "CompoundStmt" for c in node.get("inner", []))
            if has_body:
                result = translator.translate_function(node)
                print(result)
                print("")

if __name__ == "__main__":
    main()
