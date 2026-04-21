# cpp2lean v2 重构 TODO

> 追踪 `translator-v2-plan.md` 的执行进度。每项完成打勾 + 注日期。

## Stage 1：完整设计（5 周）

### Week 1：67 函数的 C++ 构造全量调研

目标：输出 `cpp-construct-catalog.md`，列出 67 函数用到的所有 C++ 构造 + 出现位置。

- [x] **Day 1 (2026-04-21)**：写调研脚本 `scripts/survey_ast.py`
  - 输入：`proof/cpp2lean/instantiate.cc` 的 Clang AST JSON
  - 输出：`survey/ast-kinds.md` / `operators.md` / `types.md` / `summary.md`
  - 发现并修：脚本挑错版本（诊断 B）+ `__taylor_coeff` 死代码（诊断 A）+ 4 处 generic lambda（诊断 C）
  - 结果：AST 完全单态化（101 dependent 节点 → 0），24 测试全通过
- [x] **Day 2 (2026-04-21)**：5 个专题扫描脚本
  - `scan_lambdas.py` → `lambdas.md`：26 in-scope lambda，0 generic
  - `scan_iterators.py` → `iterators.md`：92 range-for、4 compact-erase 关键模式
  - `scan_stl.py` → `stl.md`：std::vector 2066、pair 685、map 94、RNG 三件套
  - `scan_decomposition.py` → `decompositions.md`：30 处 DecompositionDecl，全 std::pair
  - `scan_ref_params.py` → `ref_params.md`：26 个输出参数，**8 函数配置错位**
  - 附带发现：补漏第 5 个 generic lambda（`polynomial_factorize.hh:108`）
- [x] **Day 3 (2026-04-21)**：整合 `cpp-construct-catalog.md`
  - 16 节 397 行：范围、类型、控制流、Lambda、迭代器、STL、运算符、Cast、内存、模板、函数调用、字面量、Pass 覆盖表、设计决策、缺口、衔接
  - 每 AST kind 有对应翻译策略；每 Pass 列出必须处理的 kind
- [x] **Week 1 回顾 (2026-04-21)**：调研完成 — 65 函数 + 3 实例 = 67 Lean 定义，62 AST kind，45 种运算符，1201 类型，全部登记

### Week 2：类型系统清查

目标：输出 `type-system.md`，包含完整 C++ → Lean 类型映射表。

- [x] **Day 1 (2026-04-21)**：§1 基础数值类型映射 + §2 UB 分析（685 UB 站点枚举）
  - 新增脚本：`scan_numeric_types.py`、`scan_casts.py`、`scan_ub_sites.py`
  - Agent 调研：`clpoly-model-inventory.md`（380 行）、`lean-stdlib-catalog.md`（622 行）
- [x] **Day 2 (2026-04-21)**：§3 CLPoly 核心类型（ZZ/QQ/Zp/Variable/Monomial/SparsePoly/MvPoly/Factorization + 辅助 struct）
  - Agent 调研：`mathlib-poly-types.md`（511 行）
  - Gap 分析：clpoly_model.lean 需新增 MonomialOrder enum、3 辅助 struct、补齐 Zp.inv
- [x] **Day 3 (2026-04-21)**：§4 STL 容器 shim
  - vector→Array、pair→Prod、map→StdMap（60 行完整实现）
  - sort→Array.qsort（stdlib 直接对应）+ Agent `lean-sort-api.md`
  - RNG 公理化、move/swap/iota/max/min 映射
- [x] **Day 4 (2026-04-21)**：§5 模板实例化 + §6（并入 §3.3.9）
  - `factorize` 3 实例本质是 3 个独立源函数（非同模板不同类型）
  - 命名方案：factorize_{upoly, lex, grlex}_ir
  - HIR 节点扩展：base_name + instance_suffix + mangled_name
- [x] **Day 5 (2026-04-21)**：§7 引用指针消除 + §8 CAST_TABLE
  - 26 non-const ref → tuple 返回；1 pointer 个案；ref_elim Pass 从 AST 自动推导
  - 5670 cast 全覆盖：709 显式查表、4803 noop、158 FUNC_MAP
  - CAST_TABLE Python 结构 + IntegralCast 源/目对完整枚举
- [x] **Week 2 回顾 (2026-04-21)**：type-system.md 1286 行完成所有 8 节

### Week 3：C++ 子集形式语义

目标：输出 `cpp-subset-semantics.md`，覆盖调研清查出的所有构造的 denotational semantics。

- [x] **基本运算语义 (2026-04-21)**：§1 整数算术 UB + bool 短路 + 引理 L1.1/L1.2
- [x] **控制流语义 (2026-04-21)**：§2 if/while/for/range-for/break/continue/return/do-while + L2.1-L2.3
- [x] **变量/赋值语义 (2026-04-21)**：§3 SSA 变换 + L3.1（引用 Appel 1998）
- [x] **函数调用语义 (2026-04-21)**：§4 值/const ref/ref_elim + L4.1/L4.2（引用 Leroy 2009）
- [x] **Lambda 语义 (2026-04-21)**：§5 闭包/capture/lifting + L5.1
- [x] **迭代器语义 (2026-04-21)**：§6 Range-for/compact-erase/classic + L6.1
- [x] **assert / throw 语义 (2026-04-21)**：§7 require 映射（throw 已确认不使用）
- [x] **UB 催化点 (2026-04-21)**：§8 685 站点 → require 精确对应 + L8.1/L8.2
- [x] **Week 3 回顾 (2026-04-21)**：cpp-subset-semantics.md 630 行 + semantic-references.md 167 行（Agent）
  - 主定理 L9.1 + 8 个引理分别对应 8 个 HIR/MIR Pass
  - 引用基础：Appel 1998 / Aeneas ICFP 2022 / Winskel 1993 / Leroy 2009

### Week 4：HIR 设计

目标：输出 `hir-design.md`，包含 HIR 数据结构定义 + Pass 1-5 规格。

- [x] **HIR 节点定义 (2026-04-21)**：§1 Stmt/Expr/Type/HIRFunc 完整 dataclass
- [x] **不变量阶梯 (2026-04-21)**：§2 HIR₀-HIR₄ 阶梯 + runtime assert 函数
- [x] **Pass 1 规格 (2026-04-21)**：§3 AST kind → HIR 节点映射表（覆盖 62 种 kind）
- [x] **Pass 2 规格 (2026-04-21)**：§4 ref_elim 算法 + 从 AST 自动推导（废弃 v1 OUTPUT_PARAMS 手工表）
- [x] **Pass 3 规格 (2026-04-21)**：§5 lambda_lift + modified capture 写回（保守策略 [&] 全部返回）
- [x] **Pass 4 规格 (2026-04-21)**：§6 iter_recognize 4 种模式（range-for 92、compact-erase 4、classic 1、parallel 1）
- [x] **Pass 5 规格 (2026-04-21)**：§7 operator_resolve + CAST_TABLE 查表 + UB require 生成
- [x] **手动走查 (2026-04-21)**：§8 __lll_reduce / __mtshl_step_j / __zassenhaus_recombine 三函数走过 Pass 1-5
- [x] **Week 4 回顾 (2026-04-21)**：hir-design.md 1212 行 + v1-reuse-inventory.md
  - HIR 层预估 ~2350 行代码；v2 总量 ~4250 行（v1 5956 → 减 29%）

### Week 5：MIR + codegen + 测试框架设计（2026-04-21 完成）

目标：输出 `mir-design.md`、`codegen-design.md`、`b2b-design.md`。

- [x] **MIR 节点定义 (2026-04-21)**：§1 CFG/BasicBlock/Terminator/PhiStmt/MIRFunc
- [x] **Pass 6 规格 (2026-04-21)**：§2 Cytron 4 阶段（CFG + DomTree + DF + Phi+Rename）
  - 选 Cooper-Harvey-Kennedy 迭代 DomTree（简单 20 行）
  - 选 Minimal SSA（无需活跃性分析）
  - 预估 ~600 行 Python
- [x] **Pass 7 规格 (2026-04-21)**：§3 循环提取 + flag 下降
  - 自然循环识别 + live_vars 分析
  - Flag 方案（_ret_flag + _ret_val），不引入 Except
- [x] **Pass 8 规格 (2026-04-21)**：§4 CFG 线性化 + Lean 源代码生成
  - 按支配关系 DFS 嵌套 let 链
  - aux_defs 拓扑排序输出
- [x] **back-to-back 测试框架 (2026-04-21)**：§5
  - C++ 侧：~500 行 harness 手写 emitter（零第三方依赖）
  - Lean 侧：`Lean.Data.Json` + ToJson 实例（~30 行）
  - Diff 工具：~50 行 Python（`math.isclose` 浮点近似）
  - JSON Lines 格式（升级 v1 的空格分隔）
  - **警告**：Stage 2 必须替换 v1 的伪原语（`polynomial_GCD_ir = fun f g => f`）
- [x] **单元测试 + 回归框架 (2026-04-21)**：§7 ~300 行跨 8 个 pass
- [x] **Week 5 回顾 (2026-04-21)**：mir-design.md 1008 行 + mir-references.md 240 行
  - v2 总量从 4250 修正为 ~5300 行（包含 b2b 框架 + 单元测试）
  - Lean 无现成 CFG/DomTree 基础设施，全部 Python 内部实现
  - phi 节点在 codegen 下降为 let / 尾递归参数，Lean 不暴露 phi

## Stage 2：实现（4-5 周）

> Stage 1 完成后再规划细节

- [ ] 创建 `proof/cpp2lean_v2/` 目录结构
- [ ] HIR 数据结构实现（`hir.py`）
- [ ] Pass 1-5 实现 + 单元测试
- [ ] MIR 数据结构实现（`mir.py`）
- [ ] Pass 6（ssa_build）实现 + 单元测试
- [ ] Pass 7（loop_lower）实现 + 单元测试
- [ ] Pass 8（codegen）实现
- [ ] back-to-back 测试框架实现
- [ ] 回归测试基建
- [ ] Stage 2 验收：`__make_zp` 走完全管线 + back-to-back 通过

## Stage 3：67 函数 bring-up（3-4 周）

> Stage 2 完成后再规划细节

按依赖拓扑序推进，每函数 back-to-back 通过即结束。

- [ ] 批次 1：叶节点工具函数（~10 个）
- [ ] 批次 2：基础多项式运算（~10 个）
- [ ] 批次 3：Zp 模块 5 函数
- [ ] 批次 4：Univar Hensel（~10 个）
- [ ] 批次 5：Univar 重组（~6 个）
- [ ] 批次 6：Univar 顶层 + factorize 单变量路径
- [ ] 批次 7：Wang MTSHL 内部（~10 个）
- [ ] 批次 8：Wang 管线（~6 个）
- [ ] 批次 9：factorize 多变量路径
- [ ] Stage 3 验收：67 函数全通过 + factorize 三种输入端到端

## Stage 4：精化证明（3-4 周，可并行到 Stage 3）

> Stage 2 完成后再规划细节

- [ ] `__squarefree_Zp` L1 ≈ L2
- [ ] `__ddf_Zp` L1 ≈ L2
- [ ] `__edf_Zp` L1 ≈ L2
- [ ] `__factor_Zp` L1 ≈ L2
- [ ] `__select_prime` L1 ≈ L2
- [ ] Hensel 链精化
- [ ] `__factor_recombine` L1 ≈ L2
- [ ] `__wang_core` L1 ≈ L2
- [ ] `__mtshl_lift` L1 ≈ L2
- [ ] `__factor_multivar` L1 ≈ L2
- [ ] 顶层 `factorize` 端到端精化定理

## 当前状态

**今天（2026-04-20）**：方案已敲定，准备开始 Stage 1 Week 1 Day 1 — 写调研脚本 `scripts/survey_ast.py`。
