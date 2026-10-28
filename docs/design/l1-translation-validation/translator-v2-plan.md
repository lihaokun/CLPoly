# cpp2lean 翻译器 v2 重构方案

> 制定日期：2026-04-20
> 背景：v6 审核暴露翻译器架构级缺陷（PASS 率从 v5 的 91% 跌至 12%），6 轮 R/S/M/N/P 根因迭代仍无法收敛，需要从零重新设计

## 1. 目标与范围

### 1.1 目标

将 CLPoly 的 C++ 因式分解代码（`clpoly/polynomial_factorize*.hh`）翻译为 Lean 4 L1 IR，严格 1:1 对应 C++ 语义，用于 L1→L2 精化证明。

### 1.2 范围

**翻译目标**：67 个活跃函数（模板实例化后），分布：
- Zp 模块：13 函数
- Univar 模块：34 函数（含顶层 `factorize`）
- Wang 模块：20 函数（不含 5 个死代码）

**已知死代码（不在范围）**：`__multivar_hensel_lift`、`__hensel_lift_one_var`、`__multivar_diophantine`、`__hensel_lc_correct`、`__pseudo_remainder_x1`（经典 Wang Hensel 路线，已被 MTSHL 替代，无调用者）

## 2. 当前翻译器为什么必须重构

v1-v6 经历 6 轮 "发现根因 → 修 → 暴露新根因" 的循环：

| 版本 | "根因"类别 | PASS 率 |
|---|---|---|
| v1 | R1-R6 | 38% |
| v2 | S1-S5 | 42% |
| v3 | M1-M5 | 64% |
| v4 | M1 回归 | 53% |
| v5 | N1-N5 | 91% |
| v6 | P1 pop-and-return | **12%** |

不存在"可枚举的根因"，是架构层面的病。深度根因：

### R0.1 没有形式语义，只有模式匹配
`ssa_transform.py` 2344 行全是 `isinstance` + 字符串匹配的启发式。没有定义过 C++ 子集 semantics、Lean IR invariant、翻译的保持性质。

### R0.2 SSA 变换是手搓的、ad hoc 的
没用标准 Cytron 算法，`VarEnv.versions` 手动 bump，phi 靠 fork env + merge，循环提取靠手动维护 `modified_names`。

### R0.3 没有 IR 分层
AST → SSA 一步到位，一个 Pass 同时做变量版本化、引用消除、循环提取、类型推断、输出参数重写、Lambda 提取、迭代器模式识别、错误路径包装。一锅炖。

### R0.4 没有机械化验证
每次改完 `python gen_full.py > /tmp/ir.lean; grep sorry; 读一下`。`0 sorry` 只证明"生成的代码能 parse"，不能证明语义。

### R0.5 三个目标互相冲突，隐性妥协失控
1:1 忠实 / Lean 可执行 / 0 sorry 冲突时，代码用"静默兜底"解决（`let _ := expr`、裸 `__coll`、`/- unary () -/`、`(id (0 : UInt64))` 当 ZZ 常量）。

## 3. 架构：3 IR + 8 Pass

### 3.1 IR 层次

| IR | 所处阶段 | 关键特性 |
|---|---|---|
| **AST** | Clang 原始输出 | JSON，含完整 C++ 语义细节（不由我们设计）|
| **HIR** | 引用/Lambda/迭代器消除后 | 纯函数式输入输出（无 `T&`）、Lambda 已提升为独立 def、迭代器模式已转为高阶函数、运算符已解析。**可变赋值仍保留**。控制流仍是 if/while/for 结构化形式 |
| **MIR** | SSA + 循环下降后 | 每变量单赋值、phi 节点显式、循环已提取为 `partial def` 尾递归、break/continue/return 下降为 flag+返回 |

### 3.2 Pass 列表

```
AST ──parse──▶ HIR₀ ──ref_elim──▶ HIR₁ ──lambda_lift──▶ HIR₂
                                                           │
                                            iter_recognize ▼
MIR₁ ◀──loop_lower── MIR₀ ◀──ssa_build── HIR₄ ◀──operator_resolve── HIR₃
 │
 ▼ codegen
Lean
```

| # | Pass | IR 变化 | 工作内容 |
|---|---|---|---|
| 1 | `parse` | AST → HIR₀ | Clang AST → 我们定义的 HIR 节点 |
| 2 | `ref_elim` | HIR₀ → HIR₁ | `T& out` 参数转为返回值 tuple；函数签名改写 |
| 3 | `lambda_lift` | HIR₁ → HIR₂ | Lambda 提为 `_lambda_N`；capture 作为参数；modified capture 加入返回值 tuple |
| 4 | `iter_recognize` | HIR₂ → HIR₃ | `begin/end/++` 模式识别为 `filter`/`map`/`compact` 等高阶操作 |
| 5 | `operator_resolve` | HIR₃ → HIR₄ | 运算符重载根据 Clang 类型信息解析到具体函数名（`Zp::operator+` → `Zp.add`）|
| 6 | `ssa_build` | HIR₄ → MIR₀ | CFG 构造 + dominance frontier + phi 放置 + 变量重命名（Cytron 算法）|
| 7 | `loop_lower` | MIR₀ → MIR₁ | `while`/`for` 提取为 `partial def` 尾递归；break/continue/return 下降为 flag 机制 |
| 8 | `codegen` | MIR₁ → Lean | 类型映射 + 字符串生成 |

### 3.3 HIR 内部 Pass 的数据结构

HIR₀-HIR₄ 共用数据结构，但不变量递进：

- **HIR₀**：原始，允许一切（`T&` 参数、inline Lambda、裸迭代器、未解析的运算符）
- **HIR₁**：保证无 `T&` 参数
- **HIR₂**：保证无 inline Lambda
- **HIR₃**：保证无裸 iterator
- **HIR₄**：保证所有 `Call` 的 callee 已解析为具体函数名

每个 Pass 入口 assert 前置条件，出口 assert 后置条件（~20 行/Pass 的 runtime 检查）。

### 3.4 HIR vs MIR 的数据结构差异

MIR 和 HIR 是不同的 dataclass 家族，Python type hint 静态可查：

- MIR 有 `PhiStmt` 节点，HIR 没有
- HIR 有 `WhileStmt`/`ForStmt`/`RangeForStmt`，MIR 没有（已被 `TailCallStmt` 替代）
- HIR 的 `LetStmt` 允许同名变量重复出现（mutation），MIR 不允许
- HIR 有 `BreakStmt`/`ContinueStmt`/`ReturnStmt`，MIR 只有 `TailCallStmt` 和 `ValueStmt`

## 4. 四个 Stage

### Stage 1：完整设计（5 周）

**不碰代码，只产文档。** 设计必须一次性覆盖全部 67 个目标函数的构造。

| 周 | 工作 | 产出 |
|---|---|---|
| 1 | 67 函数的 C++ 构造全量调研 | `cpp-construct-catalog.md` |
| 2 | 类型系统清查 | `type-system.md` |
| 3 | C++ 子集形式语义 | `cpp-subset-semantics.md` |
| 4 | HIR 数据结构 + Pass 1-5 规格 | `hir-design.md` |
| 5 | MIR 数据结构 + Pass 6-7 规格 + codegen + back-to-back 框架 | `mir-design.md`、`codegen-design.md`、`b2b-design.md` |

**Stage 1 验收标准**：
- 任一 C++ 构造出现在 67 函数中，都能在语义文档找到定义
- 任一构造都能追踪到"由 HIR 的 X 节点表示 / 由 MIR 的 Y 变换处理"
- 手动走一遍 `__lll_reduce`、`__mtshl_step_j`、`__zassenhaus_recombine` 的设计管线，能得出正确的 Lean IR 框架

### Stage 2：实现（4-5 周）

按 Stage 1 设计实现翻译器。新代码放 `proof/cpp2lean_v2/`，旧的保留作参考。

| 工作 | 估时 |
|---|---|
| HIR 数据结构 + parse + ref_elim + lambda_lift + iter_recognize + operator_resolve | 1.5 周 |
| MIR 数据结构 + ssa_build（CFG + dominance + Cytron） | 1.5 周 |
| loop_lower（循环提取 + break/continue/return 下降）| 0.5 周 |
| codegen + back-to-back 测试框架 | 1 周 |
| 单元测试 + 回归框架（跨 Stage 2 全程）| — |

**Stage 2 验收标准**：至少 1 个 trivial 函数（`__make_zp`）走完全管线，back-to-back 测试通过。

### Stage 3：67 函数 bring-up（3-4 周）

按依赖拓扑序批量推进。每函数 back-to-back 通过即结束。

**推进顺序建议**：
1. 叶节点工具函数（`__make_zp`、`__binomial`、`__isqrt_ceil` 等）
2. 基础多项式运算（`__upoly_mod`、`__upoly_divmod`、`__upoly_mul_mod` 等）
3. Zp 模块（`__squarefree_Zp` → `__ddf_Zp` → `__edf_Zp` → `__factor_Zp`）
4. Univar Hensel 模块
5. Univar 重组（Zassenhaus → Van Hoeij + LLL）
6. 顶层 `factorize`（单变量路径）
7. Wang MTSHL 内部
8. Wang 管线（`__wang_core` → `__factor_multivar`）
9. 顶层 `factorize` 多变量路径

**Stage 3 验收标准**：
- 67 函数全部 back-to-back 通过
- `factorize` 对三种输入（`upolynomial_<Zp>`、`upolynomial_<ZZ>`、`polynomial_<ZZ, grlex>`）端到端通过 10+ 测试用例

### Stage 4：精化证明（3-4 周，可与 Stage 3 并行）

每个 pipeline 函数证 L1 ≈ L2。L2 已有 6276 行 0 sorry 模型，精化就是桥接两层。

## 5. 总投入估计

| Stage | 周数 |
|---|---|
| 1. 完整设计 | 5 |
| 2. 实现 | 4-5 |
| 3. Bring-up | 3-4 |
| 4. 精化（可并行到 3）| 3-4 |
| **合计** | **13-16 周** |

## 6. 风险与缓解

### 6.1 单 IR 退化风险
风险：Pass 之间的"阶段约束"靠文档维护，纪律不严时可能退化（如 SSA 后又引入重复变量名）。
缓解：每 Pass 入口/出口 runtime assert；CI 跑 invariant 检查脚本。

### 6.2 Stage 1 设计遗漏
风险：调研时漏掉某个 C++ 构造，Stage 3 bring-up 时被迫回来改设计。
缓解：Stage 1 Week 1 的调研必须**机械扫描**全 67 函数的 Clang AST JSON，输出所有 `kind` 字段的直方图。人工审查不依赖记忆。

### 6.3 STL 原语形式化
风险：`std::sort`、`std::random_device`、`std::map` 各自需要 Lean 端的 opaque model 或实现。
缓解：Stage 1 Week 2 的类型系统清查显式列出所有 STL 依赖，作为 `clpoly_model.lean` 的 TODO。

### 6.4 模板实例化爆炸
风险：`upolynomial_<ZZ>`、`upolynomial_<Zp>`、`polynomial_<ZZ, lex>`、`polynomial_<ZZ, grlex>` 翻译后是独立定义还是共用。
缓解：Stage 1 Week 2 的类型系统清查明确模板单态化策略（推荐：每个实例化版本独立 Lean 类型别名）。

## 7. 与现有代码的关系

### 7.1 保留
- `proof/cpp2lean/clang_hybrid.py` + `clang_ast.py` 的 AST 解析部分（可复用到新 parse Pass）
- `proof/cpp2lean/class_map.py` 的 `FUNC_MAP`、`CLASS_MAP`、`TRANSLATION_SCOPE`、`clpoly_model.lean` 声明
- `proof/cpp2lean/lean_codegen.py` 的 prelude 生成（可复用到新 codegen）
- `test/`、`bench/`、相关 C++ 辅助

### 7.2 废弃（不删，保留作对比）
- `proof/cpp2lean/ssa_transform.py`（2344 行）
- `proof/cpp2lean/ir_types.py`
- `proof/cpp2lean/ub_collector.py`
- `proof/cpp2lean/gen_full.py`

### 7.3 新建
- `proof/cpp2lean_v2/hir.py` — HIR 数据结构
- `proof/cpp2lean_v2/mir.py` — MIR 数据结构
- `proof/cpp2lean_v2/passes/` — 8 个 Pass 各一个文件
- `proof/cpp2lean_v2/codegen.py` — Lean 代码生成
- `proof/cpp2lean_v2/b2b_runner.py` — back-to-back 测试框架
- `proof/cpp2lean_v2/main.py` — 主入口（替代 `gen_full.py`）
- `proof/cpp2lean_v2/tests/` — 单元测试 + 回归测试

## 8. 工作流与开发日志

按 `docs/workflow.md` 规范：
- 每 Stage 完成后在 `docs/devlog/` 写一条开发日志
- 每 Pass 完成后在 `docs/devlog/` 写一条开发日志
- 每个重要决策（如 STL 原语映射方式）写入 `docs/design/l1-translation-validation/` 对应文档
