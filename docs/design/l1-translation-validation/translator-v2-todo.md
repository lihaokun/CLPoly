# cpp2lean v2 重构 TODO

> 追踪 `translator-v2-plan.md` 的执行进度。每项完成打勾 + 注日期。

## Stage 1：完整设计（5 周）

### Week 1：67 函数的 C++ 构造全量调研

目标：输出 `cpp-construct-catalog.md`，列出 67 函数用到的所有 C++ 构造 + 出现位置。

- [ ] **Day 1**：写调研脚本 `scripts/survey_ast.py`
  - 输入：`proof/cpp2lean/instantiate.cc` 的 Clang AST JSON
  - 输出：
    - 所有 `kind` 字段的直方图（哪些 AST 节点出现过、各出现多少次）
    - 每个 `kind` 列出首次出现的函数名 + 行号（便于查看样例）
    - 所有 `operator*` 类型的直方图（一元/二元/赋值/下标/调用）
    - 所有模板实例化类型的清单
- [ ] **Day 2-3**：人工补全调研（脚本抓不到的语义细节）
  - Lambda capture 类型（`[&]`、`[=]`、混合捕获）
  - 迭代器模式（`for (auto it = c.begin(); it != c.end(); ++it)`、range-for、双指针 compact）
  - 结构化绑定（`for (auto& [k, v] : m)`）
  - STL 依赖列表（`std::sort`、`std::swap`、`std::max/min`、`std::random_device`、`std::mt19937`、`std::uniform_int_distribution`、`std::map`、`std::pair`、`std::vector`）
  - 模板实例化策略（`upolynomial_<T>` 各实例）
  - 引用参数模式（const ref、non-const ref、output-only）
- [ ] **Day 4-5**：起草 `cpp-construct-catalog.md`
  - 按类别组织：控制流、类型、运算符、STL、Lambda、迭代器、模板
  - 每构造列出：语义、出现函数、典型 AST 样例
- [ ] **Week 1 回顾**：与用户对齐调研覆盖完整

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
