# cpp2lean v2 重构 — Stage 1 Week 2 完成

日期：2026-04-21

## 做了什么

Stage 1 Week 2（类型系统清查）全部完成，同 Week 1 一样 1 天走完原 5 天计划。

### Day 1 — §1 基础数值 + §2 UB 分析
- 3 路并行：写 scan_numeric_types / scan_casts / scan_ub_sites + 2 个 Agent（clpoly_model 盘点 + Lean stdlib 目录）
- 产出：type-system.md §1+§2、6 份 survey
- **关键数据**：int32_t 3017、uint64_t 471、685 UB 站点枚举（UB-2 OOB 468、UB-6 signed 113）

### Day 2 — §3 CLPoly 核心类型
- Agent 调研 Mathlib 多项式类型（511 行 mathlib-poly-types.md）+ 我写 scan_clpoly_types.py
- 产出：§3（~340 行）覆盖 ZZ/QQ/Zp/Variable/Monomial/SparsePoly/MvPoly/Factorization
- **关键决策**：Zp 保留运行时 prime（ZMod p 编译时不兼容）；MvPoly 加 order 字段区分 lex/grlex

### Day 3 — §4 STL 容器 shim
- Agent 调研 Lean sort API（437 行 lean-sort-api.md）
- 产出：§4（~200 行）
- **关键**：StdMap 完整 Lean 实现（~60 行）、Array.qsort 直接对应 std::sort

### Day 4 — §5 模板实例化
- 查源码理清 factorize 3 实例的本质（3 个独立重载/特化）
- 产出：§5 + §6（已在 §3.3.9 覆盖）
- **关键发现**：factorize_grlex_ir 是转换器（把 grlex 转 lex 后调用 lex 版），因此 __factor_multivar 等下游函数只需 1 个 `<less>` 实例

### Day 5 — §7 引用指针 + §8 CAST_TABLE
- 产出：§7（~80 行）+ §8（~150 行）
- **§7 关键**：ref_elim Pass 从 AST 自动推导输出参数，不再依赖手工维护的 OUTPUT_PARAMS 表（v1 有 8 处配置错位）
- **§8 关键**：5670 cast 分 11 种 castKind，只 709 处需 CAST_TABLE 查表（IntegralCast + IntegralToFloating + FloatingToIntegral），4803 处 noop，158 处 FUNC_MAP

## 为什么做

v2 重构 Stage 1 要把所有 C++ 类型/转换/引用/模板实例的 Lean 映射**一次性敲定**，避免实现阶段返工。

Week 1 产物 `cpp-construct-catalog.md` 列出"有哪些构造"，Week 2 产物 `type-system.md` 列出"每构造映射到什么 Lean 类型/函数"。两份文档合起来给 Week 4-5 的 HIR/MIR 设计提供完整可查的输入。

## 关键决策及其理由

1. **`ZZ` → `Int`（直接 abbrev）而非自定义 struct**：GMP 的任意精度与 Lean `Int` 语义完全一致，L1 层不需要暴露 GMP 的 small-int 优化路径。若未来需精细模型可加 `structure ZZ { ... }` 但暂不必要。

2. **`Zp` 保留运行时 prime 字段**：Mathlib `ZMod p` 要求 p 编译期已知（`[Fact (p.Prime)]` 实例）。CLPoly 的 Zp 是运行时对象（prime 存在字段里）。L1 用自定义 struct，L2 精化时通过 `Zp.toZMod` 映射。

3. **`MvPoly` 加 `order : MonomialOrder` 字段**：factorize 的 3 实例中有 lex 和 grlex 两种多变量版本，不能在 L1 静态区分（都是 `polynomial_<ZZ, ...>`）。用运行时字段保留 order 信息。

4. **`StdMap` 用按 key 排序的 Array 实现**：CLPoly 的 map 都是小尺寸（度数 0..n 作 key），Array + 线性/二分查找性能足够，实现简单且 computable。不引入红黑树。

5. **`factorize` 3 实例独立 Lean 定义**：3 实例源是 3 个 C++ 独立重载（非同模板不同具体化），body 结构差异大，不能共享翻译代码。生成 `factorize_{upoly,lex,grlex}_ir` 3 个独立 def。

6. **`ref_elim` 从 AST 自动推导**：v1 手工维护 `TRANSLATION_SCOPE_OUTPUT_PARAMS` 有 8 处配置错位（scan_ref_params 发现）。v2 直接从 Clang `ParmVarDecl.type.qualType` 读 `&` 修饰，零维护成本。

7. **CAST_TABLE 只覆盖必需**：5670 cast 中 85% 是 noop（NoOp/LValueToRValue/FunctionToPointerDecay 等），HIR 阶段整体消除。只有 IntegralCast（704）+ 浮点转换（5）需显式查表（709 处），加 FUNC_MAP 覆盖 ConstructorConversion（158）。避免 263 unique 三元组的"全量表"陷阱。

## 遇到的问题与解决方式

- **问题 1**：scan_clpoly_types.py 的 qualType 前缀匹配太严（期望 `clpoly::Zp` 但 AST 里多数裸 `Zp`），数据不准
- **解决**：改用 `types.md` 的现成 qualType 直方图（Day 1 已扫），不重复扫描

- **问题 2**：第一眼以为 factorize 3 实例是"同模板不同具体化"，但 __factor_multivar 只有 1 实例令人困惑
- **定位**：查源码发现 polynomial_factorize.hh §8.5（lex 主版）+ §8.6（grlex 转换器）+ polynomial_factorize_univar.hh §8.8（univariate 非模板）是**3 个独立定义**，不是同一模板
- **结论**：`factorize_grlex_ir` 只是转换器，下游实际走 lex 版，所以 `__factor_multivar` 只需 1 个 `<less>` 实例

## 量化结果

| 指标 | 数值 |
|---|---|
| type-system.md 总行数 | **1286** |
| 8 节覆盖 | §1 基础数值 + §2 UB + §3 CLPoly + §4 STL + §5 模板 + §6 struct + §7 ref/ptr + §8 CAST_TABLE |
| 新增 scan 脚本 | 4 个（numeric_types、casts、ub_sites、clpoly_types） |
| Agent 调研文档 | 4 份（clpoly-model-inventory、lean-stdlib-catalog、mathlib-poly-types、lean-sort-api），合计 1950 行 |
| UB 站点 | 685（7 类，UB-2 OOB 468 主导）|
| CAST_TABLE | 709 处显式查表 + 4803 noop + 158 FUNC_MAP = 5670 全覆盖 |
| clpoly_model.lean 需新增代码 | ~60-80 行（StdMap 8 方法 + 3 辅助 struct + MonomialOrder + Zp.inv）|
| 测试 | Day 1 的 4 lambda edit 后 24 test files 全通过 |

## 涉及的文件

### 修改
- `docs/design/l1-translation-validation/translator-v2-todo.md`（Week 2 全打勾）

### 新增
- `docs/design/l1-translation-validation/survey/type-system.md`（1286 行，Week 2 主产出）
- `docs/design/l1-translation-validation/survey/numeric-types.md`
- `docs/design/l1-translation-validation/survey/casts.md`
- `docs/design/l1-translation-validation/survey/ub-sites.md`
- `docs/design/l1-translation-validation/survey/clpoly-types-usage.md`
- `docs/design/l1-translation-validation/survey/clpoly-model-inventory.md`（Agent）
- `docs/design/l1-translation-validation/survey/lean-stdlib-catalog.md`（Agent）
- `docs/design/l1-translation-validation/survey/mathlib-poly-types.md`（Agent）
- `docs/design/l1-translation-validation/survey/lean-sort-api.md`（Agent）
- `proof/cpp2lean/scripts/scan_numeric_types.py`
- `proof/cpp2lean/scripts/scan_casts.py`
- `proof/cpp2lean/scripts/scan_ub_sites.py`
- `proof/cpp2lean/scripts/scan_clpoly_types.py`

## 下一步

**Stage 1 Week 3：C++ 子集形式语义** — 基于 `cpp-construct-catalog.md` + `type-system.md` 输出 `cpp-subset-semantics.md`：
- 每构造的 denotational semantics
- 每 UB 点的触发条件与 Lean `require` 的精确对应
- 语义保持论证的核心引理（基础类型算术、SSA、循环、lambda、iterator）
- 为 Week 4-5 的 HIR/MIR 设计提供形式化依据

## 度量

- 耗时：~半天（3 路并行 + Agent 调研 + 集中写作）
- 迭代：0 轮修正（每节写完即 commit，无重写）
- Lean 新增/修改行数：0（仅文档；clpoly_model.lean 代码改动留到 Week 4 实现阶段）
- 对应 C++ 行数：4852（65 因式分解函数的 C++ 源）
- 放弃的方案：
  - `std::map` 用红黑树实现（改为 Array，避免不必要的复杂度）
  - CAST_TABLE 全量 263 三元组（改为只覆盖 709 关键 cast）
  - `Zp` 直接用 `ZMod p`（因 prime 编译时/运行时差异放弃）
