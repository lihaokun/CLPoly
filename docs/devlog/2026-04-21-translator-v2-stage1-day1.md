# cpp2lean v2 重构 — Stage 1 Week 1 Day 1

日期：2026-04-21

## 做了什么

1. **写了重构方案文档**：`docs/design/l1-translation-validation/translator-v2-plan.md`（3 IR + 8 Pass 架构 + 4 Stage 时间表）和 `translator-v2-todo.md`（分周分日的执行清单）。
2. **Stage 1 Week 1 Day 1：写 AST 调研脚本并执行**：
   - `proof/cpp2lean/scripts/survey_ast.py`：扫 TRANSLATION_SCOPE 里每个函数的 Clang AST，输出 kind 直方图、运算符直方图、类型直方图
   - `proof/cpp2lean/scripts/drill_dependent.py`：深挖 dependent AST 节点分布
   - `proof/cpp2lean/scripts/verify_lambda_hypothesis.py`：验证"所有 dependent 节点都在 generic lambda 里"的假设
   - `proof/cpp2lean/scripts/enumerate_instances.py`：枚举每个函数名对应的全部 template 实例化
3. **发现并修复 3 个问题**：
   - 诊断 B（翻译器挑版本 bug）：`survey_ast.py` 原来用 `has_body + [-1]` 挑候选，会选到 FunctionTemplateDecl 的未实例化**模板定义**（`mangledName=None`、含 `var_order`/`comp` 等模板参数）。改为**优先选 `mangledName ≠ None` 的 FunctionDecl**。
   - 诊断 A（死代码）：`__taylor_coeff` (ZZ 版) 只被已排除的 `__multivar_diophantine` 和 `__hensel_lift_one_var` 调用（5 个经典 Wang Hensel 死代码函数的一部分），自身也是死代码。从 `class_map.py` 的 `TRANSLATION_SCOPE` 移除。
   - 诊断 C（generic lambda）：4 处 `[](const auto& a, const auto& b) { ... }` 导致 lambda `operator()` 本身是 function template，body 里 `a.first` / `get_deg(a)` 未单态化。替换为具体类型。

## 为什么做

**v6 审核（2026-04-20）揭示翻译器架构级缺陷**——PASS 率从 v5 的 91% 跌至 12%，6 轮 R/S/M/N/P 根因迭代无法收敛。决定重构翻译器（路线 A）。Stage 1 首要任务是对全 65 个目标函数做机械化的 C++ 构造调研，为 IR 和 Pass 设计提供**一次性完整**的输入清单，避免"先搞定简单函数再扩展"导致的返工。

调研里发现 Clang 的 AST dump 混着模板定义和实例化，需要区分；也发现 CLPoly 源里有 4 处 generic lambda 阻止完整单态化。直接修掉这些是为了让 Stage 1 Week 3 的语义文档**不需要包含模板/dependent 规则**，大幅降低 IR 设计复杂度。

## 关键决策及其理由

1. **选 Option 1（改源码）而不是 Option 2（翻译器处理 generic lambda）**：
   - Option 1 零运行时影响（显式类型就是 `auto` 解析出的类型），5 处小改
   - Option 2 需要翻译器理解模板实例化机制（~200 行代码 + 语义文档增章节）
   - 权衡：5 行源码改动 vs 几百行翻译器复杂度 —— 明显前者划算

2. **从 TRANSLATION_SCOPE 剔除 `__taylor_coeff`**：
   - 它的调用者（`__multivar_diophantine`、`__hensel_lift_one_var`）本身已作为死代码排除
   - 编译期 Clang 也确认无实例化（`mangledName=None`，仅模板定义）
   - 无实际调用者 = 无精化证明义务

3. **保留 `factorize` 3 个实例化**（`upolynomial_<ZZ>` / `polynomial_<ZZ, lex>` / `polynomial_<ZZ, grlex>`）：
   - 这是唯一一个同名多实例化的函数，其他 64 个函数各只有 1 个实例化
   - 翻译器架构需要能处理这种情况，后续 Stage 1 Week 4 HIR 设计要覆盖
   - 3 个实例化映射到 3 个 Lean `_ir` 定义（后缀区分），总计 67 个 Lean 定义

## 遇到的问题与解决方式

- **问题 1**：首次 survey 显示 101 个 dependent AST 节点（CXXDependentScopeMemberExpr 55、UnresolvedLookupExpr 19、TemplateTypeParmDecl 18、FunctionTemplateDecl 9），推翻了"Clang 完全单态化"的初步判断
- **解决**：写 `verify_lambda_hypothesis.py` 检查每个节点的 LambdaExpr 祖先。第一轮假设失败（54/101 不在 lambda 内，主要集中在 `factorize` 和 `__taylor_coeff`）
- **根因**：
  - 54 个节点中 37 个来自 `factorize`（因为我的脚本挑错了版本 —— 挑了模板定义）
  - 17 个来自 `__taylor_coeff`（本身是死代码）
- **最终**：修完 3 个问题（脚本 bug + 死代码 + generic lambda），**101 dependent 节点全部归零**，假设成立

## 量化结果

| 指标 | 数值 |
|---|---|
| TRANSLATION_SCOPE 调整 | 66 → 65（移除 `__taylor_coeff`）|
| 实例化分布 | 64 函数 × 1 + 1 函数（`factorize`） × 3 = 67 个 Lean 定义 |
| AST dependent 节点 | 101 → 0 |
| AST kind 种数 | 62（覆盖全部 65 函数的 AST）|
| 不同运算符种数 | 45 |
| 不同类型种数 | 1201 |
| 测试结果 | 24 test files 全通过，0 failure |
| 代码改动 | 4 处 CLPoly lambda + 1 处 class_map.py + 新增 4 个调研脚本 |

## 涉及的文件

### 修改
- `clpoly/polynomial_factorize_zp.hh`（`__factor_Zp` lambda）
- `clpoly/polynomial_factorize_wang.hh`（`__wang_leading_coeff`、`__factor_multivar` lambda）
- `clpoly/polynomial_factorize_univar.hh`（单变量 `factorize` lambda）
- `proof/cpp2lean/class_map.py`（TRANSLATION_SCOPE 移除 `__taylor_coeff`）
- `.gitignore`（新增 `_ast_cache/`、`stderr.log`）

### 新增
- `docs/design/l1-translation-validation/translator-v2-plan.md`
- `docs/design/l1-translation-validation/translator-v2-todo.md`
- `docs/design/l1-translation-validation/survey/{ast-kinds,operators,types,summary}.md`
- `proof/cpp2lean/scripts/survey_ast.py`
- `proof/cpp2lean/scripts/drill_dependent.py`
- `proof/cpp2lean/scripts/verify_lambda_hypothesis.py`
- `proof/cpp2lean/scripts/enumerate_instances.py`

## 下一步

**Stage 1 Week 1 Day 2-3**：人工补全调研脚本抓不到的语义细节：
- Lambda capture 类型（`[&]`/`[=]`/混合）
- 迭代器模式（range-for、begin/end/++、双指针 compact）
- 结构化绑定（`for (auto& [k, v] : m)`）
- STL 依赖完整清单（`std::sort`、`std::swap`、`std::random_device` 等）
- 引用参数模式（const ref、non-const ref、output-only）
- 模板实例化策略（已知：`factorize` 有 3 个）
