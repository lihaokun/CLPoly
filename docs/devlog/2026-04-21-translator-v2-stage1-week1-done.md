# cpp2lean v2 重构 — Stage 1 Week 1 完成

日期：2026-04-21

## 做了什么

Stage 1 Week 1（67 函数 C++ 构造全量调研）全部完成。1 天走完原 5 天计划：

### Day 1（上午）
- 写 `survey_ast.py` 主调研脚本 + `drill_dependent.py` / `verify_lambda_hypothesis.py` / `enumerate_instances.py` 辅助脚本
- 扫出 101 个 dependent AST 节点 → 逐步归零
  - 诊断 B（脚本挑错版本）：`survey_ast.py` 优先选 `mangledName ≠ None` 的 FunctionDecl
  - 诊断 A（死代码）：从 `TRANSLATION_SCOPE` 移除 `__taylor_coeff`
  - 诊断 C（generic lambda）：4 处 `const auto&` → 具体类型（`__factor_Zp`、`__wang_leading_coeff`、`__factor_multivar`、单变量 `factorize`）
- commit `01b1894`、`53db218`

### Day 2（下午）
- 5 个专题扫描脚本 + 对应 survey 报告：
  - `scan_lambdas.py` + `lambdas.md`（26 in-scope lambda 细节）
  - `scan_iterators.py` + `iterators.md`（92 range-for + 4 compact-erase 关键模式）
  - `scan_stl.py` + `stl.md`（全 STL 依赖清单）
  - `scan_decomposition.py` + `decompositions.md`（30 处 DecompositionDecl）
  - `scan_ref_params.py` + `ref_params.md`（26 输出参数 + 8 配置错位）
- 补漏第 5 个 generic lambda（`polynomial_factorize.hh:108`，Day 1 grep glob 漏掉无下划线的文件名）
- commit `cf2866b`、`67456bf`、`46b1476`、`5034ae1`

### Day 3（晚上）
- 整合 `cpp-construct-catalog.md`（397 行 16 节）：
  - §1 翻译范围（65+3=67 Lean 定义）
  - §2 类型系统（基础数值 + CLPoly + STL 容器）
  - §3 控制流（IfStmt 275、ForStmt 145、CXXForRangeStmt 92 等；确认无 switch/goto/try/throw）
  - §4 Lambda 全清单（26 个，0 generic）
  - §5 迭代器模式（92 range-for、4 compact-erase）
  - §6 STL 算法（6 种 90 调用）
  - §7 运算符（45 种）
  - §8 Cast / §9 内存（全零）/ §10 模板（无 dependent）/ §11 函数调用 / §12 字面量
  - §13 Pass-AST 覆盖表（各 Pass 必须处理的 kind）
  - §14 已知设计决策 / §15 设计缺口（交 Week 3）/ §16 与后续 Week 衔接
- 打勾 `translator-v2-todo.md` 的 Week 1 条目

## 为什么做

v6 审核暴露翻译器架构级缺陷，PASS 率从 91% 跌至 12%，6 轮 R/S/M/N/P 根因迭代无法收敛。Stage 1 的核心价值：**一次性机械化调研全 65 函数**，让 IR 和 Pass 设计在纸面上完整覆盖所有 C++ 构造，避免"先搞定简单函数再扩展"导致的返工。

Week 1 的产物（`cpp-construct-catalog.md`）是**硬性验收**：任一 C++ 构造未被登记视为设计缺口，必须在 Week 3（语义）和 Week 4-5（IR 设计）之前全部覆盖。

## 关键决策及其理由

1. **Day 2 五个脚本并行做而非串行**：每脚本 ~200-350 行独立实现，共享已缓存的 AST JSON，互不干扰，压缩迭代周期到 1 天
2. **发现的 8 个 `TRANSLATION_SCOPE_OUTPUT_PARAMS` 配置错位不修**：
   - v2 的 `ref_elim` Pass 从 AST 自动推导，不再依赖手工维护表
   - 当前翻译器仍在用这份表，但我们的目标是重写，不是修补
3. **`cpp-construct-catalog.md` 只列策略不写实现**：Week 4-5 的 HIR/MIR 设计会把策略具体化为代码

## 遇到的问题与解决方式

- **问题**：Day 1 首扫显示 101 个 dependent 节点，违背"Clang 单态化 AST 是输入假设"
- **定位**：`verify_lambda_hypothesis.py` 检查每节点的 LambdaExpr 祖先 → 54 个在 lambda 外
- **根因**：不是 Clang 问题，是脚本 bug（挑到模板定义而非实例化）+ 1 个死代码 + 5 处 generic lambda（而非 4 处，第 5 个漏在无下划线文件名）
- **结果**：修完后 0 dependent，完全单态化假设成立，大幅降低 Week 3 语义文档复杂度

- **问题**：C++ 源码没有 generic lambda 的专题扫描，只能靠 grep
- **解决**：Day 2 的 `scan_lambdas.py` 正则扫 `[...]\(...auto` 模式 + 与 AST 交叉验证

## 量化结果

| 指标 | 数值 |
|---|---|
| TRANSLATION_SCOPE | 66 → **65**（移除 `__taylor_coeff`）|
| 实际 Lean 定义数 | **67**（`factorize` 有 3 实例）|
| AST kind 种数 | **62** |
| 不同运算符 | **45** |
| 不同 qualType | **1201** |
| Lambda（in-scope）| **26**（0 generic、13 无捕获、10 `[&]`、3 具名）|
| Range-for | **92**（24 含结构化绑定）|
| Compact-erase 双指针 | **4**（关键模式）|
| 输出参数 | **26**（8 函数配置错位）|
| Dependent AST 节点 | 101 → **0** |
| 测试结果 | 24 test files 全通过 |
| 代码改动 | 5 处 CLPoly lambda + 2 处 class_map.py + 新增 5 scan + 4 diag 脚本（共 9 个）+ 10 份 survey 文件 |

## 涉及的文件

### 修改
- `clpoly/polynomial_factorize.hh`（Edit 5：Day 2 补漏的 generic lambda）
- `clpoly/polynomial_factorize_zp.hh` / `polynomial_factorize_univar.hh` / `polynomial_factorize_wang.hh`（Day 1 的 4 个 lambda edit）
- `proof/cpp2lean/class_map.py`（移除 `__taylor_coeff`）
- `.gitignore`（`_ast_cache/` 等）
- `docs/design/l1-translation-validation/translator-v2-todo.md`（打勾 Week 1）

### 新增
- `docs/design/l1-translation-validation/translator-v2-plan.md`
- `docs/design/l1-translation-validation/translator-v2-todo.md`
- `docs/design/l1-translation-validation/survey/{ast-kinds,operators,types,summary,lambdas,iterators,stl,decompositions,ref_params,cpp-construct-catalog}.md`（10 份）
- `proof/cpp2lean/scripts/{survey_ast,drill_dependent,verify_lambda_hypothesis,enumerate_instances,scan_lambdas,scan_iterators,scan_stl,scan_decomposition,scan_ref_params}.py`（9 个脚本）

## 下一步

**Stage 1 Week 2：类型系统清查** — 基于 `cpp-construct-catalog.md` §2 展开为完整的 `type-system.md`：
- 每个 C++ 类型 → Lean 类型映射的详细表达
- CLPoly 自定义类型（`polynomial_<T, order>` 等）的 Lean 结构化定义
- 模板实例化处理策略（尤其 `factorize` 3 实例）
- STL 原语的 Lean 端 shim 设计（含 `StdMap`、`Rng`、`Array.qsortWith` 等）
- 引用/指针消除约定
- 类型转换表（`CAST_TABLE`）枚举
