# v1 翻译器代码复用盘点（Agent 调研）

> 由 Agent 调研 v1 cpp2lean 10 个模块（5956 行）的复用性，为 v2 HIR/MIR 设计分类。
> 日期：2026-04-21

## §1 文件分类表

| v1 模块 | 行数 | 分类 | v2 对应 Pass | 复用占比 |
|---|---|---|---|---|
| `clang_ast.py` | 1287 | **PARTIAL** | Pass 1 (parse) | 30% |
| `clang_hybrid.py` | 192 | **COPY** | Frontend | 100% |
| `ir_types.py` | 238 | **PARTIAL** | 数据结构 | 50% |
| `ssa_transform.py` | 2344 | **REWRITE** | Pass 2-7 拆分 | 0%（思路参考 10%）|
| `lean_codegen.py` | 619 | **PARTIAL** | Pass 8 (codegen) | 48% |
| `class_map.py` | 531 | **COPY** | 全局查表 | 100% |
| `ub_collector.py` | 311 | **REWRITE** | Pass 2 & 5 & 6 | 30% |
| `gen_full.py` | 217 | **REFERENCE** | Orchestration | 新写 |
| `gen_b2b_lean.py` | 145 | **REFERENCE** | Test only | 保留 |
| `cpp2lean.py` | 72 | **DISCARD** | 新 CLI | 0% |

**合计**：v1 5956 行 → v2 预计 ~4250 行（减少 29%）

## §2 v2 各 Pass 代码量预估

| Pass | 新代码 | 复用 | 备注 |
|---|---|---|---|
| Pass 1 `parse` | 400 | 30% v1 clang_ast | AST→HIR₀ 基础映射 |
| Pass 2 `ref_elim` | 200 | 0%（v1 手工 OUTPUT_PARAMS 废弃）| 从 AST 自动推导 |
| Pass 3 `lambda_lift` | 300 | 少 | Lambda 提升 + modified capture 写回 |
| Pass 4 `iter_recognize` | 350 | 参考 v1 启发式 | 4 种模式机械化 |
| Pass 5 `operator_resolve` | 250 | 大量 class_map | 查表解析 + require 生成 |
| Pass 6 `ssa_build` | 600 | 新写（Cytron 算法）| CFG + dominance + phi |
| Pass 7 `loop_lower` | 450 | 参考 v1 循环模式 | break/continue/return flag 下降 |
| Pass 8 `codegen` | 350 | 48% v1 lean_codegen | MIR→Lean 字符串 |
| **8 Pass 小计** | **2900** | | |
| 数据结构 `ir_types.py` | 150 | 50% v1 | HIR + MIR dataclass |
| Frontend `clang_hybrid.py` | 200 | 100% v1 | 直接复用 |
| 查表 `class_map.py` | 500 | 100% v1 | 直接复用 |
| Orchestration | 200 | 新写 | Pass 串联 + CLI |
| 测试 `pass_tests/` | 300 | 新写 | 每 Pass 单元测试 |
| **总计** | **4250** | | |

## §3 关键风险区域

### §3.1 Pass 6 (ssa_build) — Cytron 算法（600 行）
- 风险：CFG 构造、支配边界计算、phi 放置错误 → 整个 SSA 化失败
- 缓解：提前编写单元测试，参考《编译器：原理、技术与工具》第 2 版 §9 章

### §3.2 Pass 4 (iter_recognize) — 迭代器模式（350 行）
- 风险：v1 启发式遗漏 compact-erase 等关键模式
- 缓解：`iterators.md` 已列举全部 4 种模式（92 range-for + 4 compact-erase + 1 classic + 1 parallel），Pass 4 逐个实现

### §3.3 Pass 2 ↔ Pass 3 交界 — lambda 内的 ref capture（200+300 行）
- 风险：Lambda 的 `[&]` capture 如果是 C++ `T&` 参数，需要嵌套 tuple 返回
- 缓解：Pass 2 先跑（消除顶层 ref），Pass 3 再跑（Lambda 内 capture 只可能是 local var 或已消除的 ref）

## §4 立即可搬的代码

1. **`clang_hybrid.py`**（192 行）：libclang + JSON dump 的 frontend — 整体搬
2. **`class_map.py`**（531 行）：知识库 CLASS_MAP / FUNC_MAP / CAST_TABLE — 整体搬
3. **`ir_types.py` 的基础节点**（100 行）：`BaseType`、`Var`、`Lit`、`BinOp`、`UnaryOp`、`CondExpr`、`Call`、`LetStmt`、`AssignStmt`、`IfStmt`、`ReturnStmt` — 搬
4. **`clang_ast.py` 的 `parse_type()`**（150 行）：类型字符串解析 — 搬
5. **`lean_codegen.py` 的基础函数**（300 行）：`gen_type()`, `gen_coercion()`, `gen_expr()` 的大部分 — 搬
6. **`ub_collector.py` 的 `_scan_expr()`**（100 行）：UB 点识别 — 改造后搬

## §5 彻底重写的代码

**`ssa_transform.py`（2344 行）整体重写**，因为它混合了 6 个独立任务：
- 参数重写 → Pass 2
- 变量版本化 → Pass 6
- 引用消除 → Pass 2
- Lambda 提升 → Pass 3
- 迭代器识别 → Pass 4
- 循环提取 → Pass 7

混在一起是 v6 崩盘的核心原因。v2 强制分离。

**`cpp2lean.py`（72 行）CLI 入口废弃**，v2 用 Pass 管道式编排，新写 CLI。

## §6 最终建议

**优先级**：
1. 先实施 Pass 1-5（HIR 层，相对独立）
2. 再做 Pass 6-8（MIR 层，依赖 Cytron）
3. class_map.py 和 clang_hybrid.py 可立即搬用无风险
4. 每 Pass 实施后立即写单元测试（`pass_tests/test_passN.py`）

**Stage 2 里程碑**：
- Week 1：基础设施（clang_hybrid 搬 + ir_types 定义 + class_map 搬）
- Week 2-3：Pass 1 + Pass 2 + 单元测试
- Week 4：Pass 3 + Pass 4
- Week 5：Pass 5
- Week 6：Pass 6（Cytron）
- Week 7：Pass 7 + Pass 8
- Week 8：集成测试 + back-to-back
