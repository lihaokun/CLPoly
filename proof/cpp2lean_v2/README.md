# cpp2lean v2

> C++ → Lean 4 翻译器，Stage 1 设计稿已完成（见 `docs/design/l1-translation-validation/survey/`）。

## 架构

```
Clang AST (JSON)
  ↓ Pass 1: parse
HIR₀ (原始)
  ↓ Pass 2: ref_elim
HIR₁ (无 T&)
  ↓ Pass 3: lambda_lift
HIR₂ (无 inline lambda)
  ↓ Pass 4: iter_recognize
HIR₃ (无裸 iterator)
  ↓ Pass 5: operator_resolve
HIR₄ (Call callee 已解析)
  ↓ Pass 6: ssa_build (Cytron)
MIR₀ (SSA + phi + CFG)
  ↓ Pass 7: loop_lower
MIR₁ (循环已提取为 partial def)
  ↓ Pass 8: codegen
Lean 源代码
  ↓ back-to-back 测试
验证 vs C++ 输出
```

## 文件组织

```
cpp2lean_v2/
├── README.md               # 本文件
├── ir_types.py            # HIR + MIR 数据类
├── clang_hybrid.py        # [从 v1 搬] Clang AST 解析 frontend
├── class_map.py           # [从 v1 搬] CLASS_MAP / FUNC_MAP / CAST_TABLE
├── main.py                # Pass 管道入口 (编排)
├── codegen.py             # Pass 8: MIR → Lean 源代码
├── passes/
│   ├── pass1_parse.py       # AST → HIR₀
│   ├── pass2_ref_elim.py    # HIR₀ → HIR₁
│   ├── pass3_lambda_lift.py # HIR₁ → HIR₂
│   ├── pass4_iter_recognize.py  # HIR₂ → HIR₃
│   ├── pass5_operator_resolve.py # HIR₃ → HIR₄
│   ├── pass6_ssa_build.py    # HIR₄ → MIR₀
│   └── pass7_loop_lower.py   # MIR₀ → MIR₁
├── b2b/
│   ├── harness.cc            # C++ 端 test runner
│   ├── lean_runner.py        # Lean 端 #eval 执行
│   └── diff.py               # JSON diff
└── tests/
    ├── fixtures/              # 测试 C++ 代码片段
    └── test_passN.py          # 每 Pass 的单元测试
```

## 快速开始

```bash
# 翻译单个函数
python main.py --func __make_zp

# 完整翻译
python main.py --all --output /tmp/clpoly_ir_v2.lean

# 单元测试
python -m pytest tests/

# back-to-back
python b2b/run_b2b.py --test-vectors tests/vectors/factor_Zp.json
```

## 设计参考

| 文档 | 内容 |
|---|---|
| `translator-v2-plan.md` | 方案概览 |
| `hir-design.md` | HIR 数据结构 + Pass 1-5 规格 |
| `mir-design.md` | MIR 数据结构 + Pass 6-8 规格 + b2b 框架 |
| `type-system.md` | 类型映射表 + CAST_TABLE |
| `cpp-subset-semantics.md` | 语义保持论证 |
| `cpp-construct-catalog.md` | C++ 构造全量清单 |

## 状态

- ✅ **Stage 1 (设计)**：2026-04-21 完成
- 🚧 **Stage 2 (实现)**：进行中
  - [x] 目录骨架
  - [x] v1 复用（clang_hybrid, class_map）
  - [ ] ir_types.py
  - [ ] Pass 1 parse
  - [ ] Pass 2-8
  - [ ] b2b 框架
