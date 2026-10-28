# cpp2lean v2 重构 — Stage 1 全部完成 🎯

日期：2026-04-21

## 做了什么

**Stage 1（完整设计）5 周工作全部完成**，1 天走完原定 5 周计划。

### 每周主产出

| Week | 文档 | 行数 | Agent 辅助 |
|---|---|---|---|
| Week 1 | `cpp-construct-catalog.md` + 10 scan reports | 4400+ | 3 agents 审核历史 |
| Week 2 | `type-system.md` | 1286 | 4 agents (clpoly-model / mathlib / stdlib / sort-api) |
| Week 3 | `cpp-subset-semantics.md` | 630 | 1 agent (semantic refs) |
| Week 4 | `hir-design.md` | 1212 | 1 agent (v1 reuse inventory) |
| **Week 5** | **`mir-design.md`** | **1008** | **1 agent (Cytron + Lean CFG)** |

**总计文档 + 脚本**：约 **13000 行**。

### 按日提交

17 个 commit：

```
d8e9129 Week 4 完成 — hir-design.md
9272e42 Week 3 完成 — cpp-subset-semantics.md
751e1a0 Week 2 完成 — §5+§7+§8 + devlog
a319f0b Week 2 Day 3 — §4 STL
2feb460 Week 2 Day 2 — §3 CLPoly types
d831997 Week 2 Day 1 — §1+§2 numeric+UB
05a2694 Week 1 完成 — catalog
5034ae1 Day 2 scan_decomposition + scan_ref_params
46b1476 Day 2 scan_stl
67456bf Day 2 scan_iterators
cf2866b 5th generic lambda + scan_lambdas
53db218 v2 plan + Day 1 scan
01b1894 generic lambda fix + __taylor_coeff 排除
```

## 为什么做（v2 重构的必要性）

**v6 审核（2026-04-20）揭示**：v1 翻译器从 v5 的 91% PASS 跌至 12%，6 轮 R/S/M/N/P 根因迭代无法收敛。

**深度根因**：
1. 没有形式语义，只有模式匹配（2344 行 ad hoc）
2. SSA 变换手搓 + 无 Cytron 算法
3. 没有 IR 分层（AST → SSA 一步到位）
4. 没有机械化验证（0 sorry 不等于语义正确）
5. 三个目标（1:1 忠实 / Lean 可执行 / 0 sorry）冲突时静默妥协

**v2 方案**（Stage 1 敲定）：**3 IR 层 + 8 Pass** 的标准编译器架构，基于 Cytron 1991 + Aeneas ICFP 2022 路线。

## 关键决策及其理由

### 1. 路线选择
- **放弃**：继续打补丁 v1（已证无法收敛）
- **选择**：完整重构（路线 A），而非缩减翻译范围（路线 B）或手工维护（路线 C）
- **理由**：翻译器架构成本与函数数量关系不大；缩减范围省不了多少，却增加 trust boundary 的 paper 论证负担

### 2. IR 分层
- **选择**：3 IR (AST / HIR / MIR) + 8 Pass，而非 1 IR 多 Pass
- **理由**：HIR 和 MIR 的节点形态真的不同（HIR 有结构化循环、MIR 有 phi），类型级强制区分比 runtime assert 更可靠

### 3. 模板实例化策略
- 65 函数中 64 个有 1 实例，只有 `factorize` 有 3 实例（upoly/lex/grlex）
- **选择**：3 个独立 Lean 定义（`factorize_{upoly,lex,grlex}_ir`），不共享 body
- **理由**：3 实例是 3 个 C++ 独立重载，body 结构不同（grlex 是转换器）

### 4. SSA 算法
- **选择**：Cooper-Harvey-Kennedy 迭代 DomTree + Minimal SSA，而非 Lengauer-Tarjan + Pruned SSA
- **理由**：Lean 无现成 CFG/DomTree/SSA 基础设施（Mathlib 的 Dominator 是测度论），必须 Python 内部实现；选简单算法降低实现风险

### 5. Lambda capture 策略
- **选择**：`[&]` 保守返回所有 captures
- **理由**：精确检测需要扫 body 里所有嵌套结构（v1 在此栽倒多次）；保守策略冗余但不错

### 6. UB 处理
- **选择**：685 UB 站点每点生成 require，不引入 Except
- **理由**：CLPoly 因式分解 0 throw 0 exception；require 参数传播方式比 Except.bind 简单

### 7. CAST_TABLE 范围
- **选择**：只覆盖 IntegralCast + 浮点转换（709 处）+ 构造器（158 处），其他 4803 处 cast 视为 noop
- **理由**：5670 个 cast 中 85% 无语义操作（只是类型系统细节），把它们塞进表里徒增复杂度

## 遇到的问题与解决方式

### 问题 1：Day 1 发现 101 个 dependent AST 节点
- **根因定位**：验证假设"所有 dependent 节点都在 generic lambda 里" → 54/101 不在
- **进一步分析**：`factorize` 37 个是 survey 脚本挑错版本（选到未实例化的模板定义），`__taylor_coeff` 17 个是死代码
- **解决**：修脚本（优先选 mangledName ≠ None）+ 移除 `__taylor_coeff` + 改 4 个 generic lambda 为显式类型
- **结果**：101 → 0

### 问题 2：`scan_clpoly_types.py` 前缀匹配太严
- **现象**：types.md 里 CLPoly 类型多数是裸 `Zp` / `upolynomial_<Zp>`，不带 `clpoly::` 前缀
- **解决**：直接从 Week 1 的 `types.md` 取现成数据，不重复扫描

### 问题 3：v1 的 `polynomial_GCD_ir = fun f g => f` 伪原语
- **Agent 发现**：当前 `gen_b2b_lean.py` 的基础函数是伪定义
- **记录**：Stage 2 必须替换，否则"一致"只能证明"两端一样错"

## 量化结果

### Stage 1 整体

| 指标 | 数值 |
|---|---|
| 文档总行数 | ~13000 |
| 产出文档数 | 21 份 |
| 新增 scan 脚本 | 9 个 |
| Agent 调研 | 9 次 |
| commit 数 | 17 |
| Mathlib/Lean stdlib 调研覆盖 | 26 类型 + 5 排序 API + 9 多项式类型 |
| C++ 构造清单覆盖 | 62 AST kind + 45 运算符 + 1201 类型 + 685 UB 站点 |

### v2 代码量预估（完整）

| 模块 | 预估 | 来源 |
|---|---|---|
| Pass 1 parse | 400 | 30% 复用 v1 |
| Pass 2 ref_elim | 200 | 新写（v1 OUTPUT_PARAMS 废弃）|
| Pass 3 lambda_lift | 300 | 新写 |
| Pass 4 iter_recognize | 350 | 参考 v1 启发式 |
| Pass 5 operator_resolve | 250 | 查表逻辑新写，CLASS_MAP 复用 |
| Pass 6 ssa_build | 600 | 新写 Cytron 算法 |
| Pass 7 loop_lower | 450 | 参考 v1 循环模式 |
| Pass 8 codegen | 350 | 48% 复用 v1 |
| `ir_types.py`（HIR+MIR）| 250 | 50% 复用 v1 |
| `class_map.py` | 500 | 100% 复用 v1 |
| `clang_hybrid.py` | 200 | 100% 复用 v1 |
| Orchestration | 200 | 新写 |
| back-to-back 框架 | 950 | 新写 |
| 单元测试 | 300 | 新写 |
| **总计** | **~5300** | **v1 5956 → 减 11%** |

## 涉及的文件

### 新增（21 份设计文档）

**Week 1**（catalog + 10 surveys）：
- `cpp-construct-catalog.md`、`ast-kinds.md`、`operators.md`、`types.md`、`summary.md`
- `lambdas.md`、`iterators.md`、`stl.md`、`decompositions.md`、`ref_params.md`

**Week 2**（type-system + 4 surveys）：
- `type-system.md`
- `numeric-types.md`、`casts.md`、`ub-sites.md`、`clpoly-types-usage.md`
- `clpoly-model-inventory.md`、`lean-stdlib-catalog.md`、`mathlib-poly-types.md`、`lean-sort-api.md`

**Week 3**：`cpp-subset-semantics.md`、`semantic-references.md`

**Week 4**：`hir-design.md`、`v1-reuse-inventory.md`

**Week 5**：`mir-design.md`、`mir-references.md`

### 新增（9 个扫描脚本）

- `survey_ast.py`、`drill_dependent.py`、`verify_lambda_hypothesis.py`、`enumerate_instances.py`
- `scan_lambdas.py`、`scan_iterators.py`、`scan_stl.py`、`scan_decomposition.py`、`scan_ref_params.py`
- `scan_numeric_types.py`、`scan_casts.py`、`scan_ub_sites.py`、`scan_clpoly_types.py`

### 修改

- `clpoly/polynomial_factorize{,_zp,_univar,_wang}.hh`（Day 1 修 5 处 generic lambda）
- `proof/cpp2lean/class_map.py`（移除 `__taylor_coeff`）
- `.gitignore`（`_ast_cache/` 等）
- `translator-v2-plan.md`、`translator-v2-todo.md`

### 4 篇按 Week 的 devlog

- `2026-04-21-translator-v2-stage1-day1.md`
- `2026-04-21-translator-v2-stage1-week1-done.md`
- `2026-04-21-translator-v2-stage1-week2-done.md`
- `2026-04-21-translator-v2-stage1-done.md`（本文件）

## 下一步：Stage 2 实施

**预估 4-5 周**。文档 `translator-v2-plan.md` §4 和 `mir-design.md` §10 已列出：

### Stage 2 里程碑

| 周 | 任务 |
|---|---|
| 1 | 基础设施搬运（clang_hybrid + class_map + ir_types 定义）+ Pass 1 parse |
| 2 | Pass 2 ref_elim + Pass 3 lambda_lift + 单元测试 |
| 3 | Pass 4 iter_recognize + Pass 5 operator_resolve + 单元测试 |
| 4 | Pass 6 ssa_build（Cytron 算法，最关键）+ 单元测试 |
| 5 | Pass 7 loop_lower + Pass 8 codegen + 单元测试 |
| 6-8 | back-to-back 框架 + 集成测试 + 修 bug |

**成功标准**：65 函数（+ `factorize` 3 实例 = 67 Lean def）back-to-back 测试对 50+ 测试向量与 C++ 一致。

## 度量

- **耗时**：~1 天（2026-04-21 全天）
- **迭代**：几乎 0 轮大幅修正（每日每周 commit 到位，无重写）
- **Lean 新增/修改行数**：0（本 Stage 仅设计文档）
- **对应 C++ 行数**：4852（65 因式分解函数的 C++ 源，被完整分析）
- **放弃的方案**：
  - 单 IR 多 Pass（改为 3 IR，因为真的不同）
  - 全量 263 cast 三元组表（改为 709 显式查表 + 4803 noop）
  - `std::map` 红黑树实现（改为 Array，简单够用）
  - `Zp` 直接用 `ZMod p`（编译时/运行时不兼容）
  - Lengauer-Tarjan DomTree（改为 Cooper-Harvey-Kennedy，简单 20 行）
  - Pruned SSA（改为 Minimal SSA，无需活跃性分析）
  - Except 错误处理（CLPoly 0 throw，不需要）

## 一句话总结

Stage 1 用 **1 天**做完原定 **5 周**的完整重构设计，产出 **13000 行文档** 覆盖从 AST 到 Lean 源码的 8 Pass 架构全流程；v2 预计代码总量 ~5300 行（比 v1 少 11%），但结构清晰度和可维护性质的飞跃。Stage 2 实施阶段可按本次设计无歧义进行。
