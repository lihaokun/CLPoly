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

- [ ] **基本类型映射**：`uint64_t` / `int64_t` / `int` / `bool` / `double` / `size_t`
- [ ] **复合类型映射**：`std::vector<T>` / `std::pair<A,B>` / `std::map<K,V>` / `std::tuple`
- [ ] **CLPoly 类型**：`ZZ` / `QQ` / `Zp` / `upolynomial_<T>` / `polynomial_<T, Order>` / `basic_monomial` / `Variable` 等
- [ ] **模板实例化策略**：每个实例化生成独立 Lean 类型别名还是共用
- [ ] **引用/指针处理**：`T&`、`const T&`、`T*`、`T**` 的映射约定
- [ ] **STL 原语 Lean 端 shim**：`std::sort` → `Array.qsortWith`、`std::map` → `StdMap`（如何实现）等
- [ ] **Week 2 回顾**：类型表能覆盖调研出的全部类型

### Week 3：C++ 子集形式语义

目标：输出 `cpp-subset-semantics.md`，覆盖调研清查出的所有构造的 denotational semantics。

- [ ] **基本运算语义**：整数算术（含 UB）、布尔、比较
- [ ] **控制流语义**：if/else、while、for、range-for、break、continue、return
- [ ] **变量/赋值语义**：声明、mutation、作用域
- [ ] **函数调用语义**：值传递、引用传递、输出参数
- [ ] **Lambda 语义**：闭包、capture-by-value、capture-by-reference
- [ ] **迭代器语义**：`begin/end/++/*` 的抽象模型
- [ ] **assert / throw 语义**：require / Except.error 的映射理由
- [ ] **UB 催化点**：除零、数组越界、移位越界、有符号溢出等
- [ ] **Week 3 回顾**：每构造语义都能追踪到 Lean 层

### Week 4：HIR 设计

目标：输出 `hir-design.md`，包含 HIR 数据结构定义 + Pass 1-5 规格。

- [ ] **HIR 节点定义**：Stmt 家族、Expr 家族的 dataclass
- [ ] **不变量阶梯**：HIR₀ → HIR₁ → HIR₂ → HIR₃ → HIR₄ 每阶段的不变量
- [ ] **Pass 1 规格**：`parse`（AST → HIR₀）的输入输出契约 + 特殊 case 处理
- [ ] **Pass 2 规格**：`ref_elim`（HIR₀ → HIR₁）的重写规则 + 输出参数 tuple 约定
- [ ] **Pass 3 规格**：`lambda_lift`（HIR₁ → HIR₂）的 capture 处理 + modified capture 写回
- [ ] **Pass 4 规格**：`iter_recognize`（HIR₂ → HIR₃）支持的迭代器模式枚举
- [ ] **Pass 5 规格**：`operator_resolve`（HIR₃ → HIR₄）的运算符解析查表
- [ ] **手动走查**：选 3 个复杂函数（`__lll_reduce`、`__mtshl_step_j`、`__zassenhaus_recombine`）手动走过 Pass 1-5，确认 HIR₄ 形式正确
- [ ] **Week 4 回顾**：HIR 能完整表达 67 函数所有构造

### Week 5：MIR + codegen + 测试框架设计

目标：输出 `mir-design.md`、`codegen-design.md`、`b2b-design.md`。

- [ ] **MIR 节点定义**：Stmt 家族（含 `PhiStmt`、`TailCallStmt`、`ValueStmt`）、Expr 家族
- [ ] **Pass 6 规格**：`ssa_build`（HIR₄ → MIR₀）
  - CFG 构造（basic block 划分 + edges）
  - Dominator tree 计算
  - Dominance frontier 计算
  - Phi 节点放置（Cytron 算法）
  - 变量重命名
- [ ] **Pass 7 规格**：`loop_lower`（MIR₀ → MIR₁）
  - 循环识别（natural loops）
  - 循环提取为 `partial def` 尾递归
  - break → `_break_flag` + 尾返回
  - continue → 循环顶部尾调用
  - return in loop → `_ret_flag` + `_ret_val` + 尾返回
- [ ] **Pass 8 规格**：`codegen`（MIR₁ → Lean）
  - 类型字符串生成
  - `partial def` 生成
  - Phi 节点 → `if-then-else` 或 `match`
  - Require 参数生成
- [ ] **back-to-back 测试框架**：
  - C++ 执行流程（编译 harness + 输入 JSON → 输出 JSON）
  - Lean 执行流程（生成 `#eval` 文件 + `lake env lean` → 输出 JSON）
  - 双向 JSON 序列化约定
  - diff 工具
- [ ] **单元测试 + 回归框架**：每 Pass 的测试脚手架
- [ ] **手动走查**：选 2 个复杂函数（`__lll_reduce`、`__zassenhaus_recombine`）手动走过 Pass 6-8，确认 Lean 代码能编译
- [ ] **Week 5 回顾**：设计文档齐全；Stage 1 验收通过

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
