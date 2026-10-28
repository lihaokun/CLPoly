# CLPoly 因式分解 C++ 子集构造完整清单

> 生成日期：2026-04-21
> 范围：`clpoly/polynomial_factorize{,_zp,_univar,_wang}.hh` 的 65 个 TRANSLATION_SCOPE 函数
> 源数据：`survey/` 下 9 份调研报告（ast-kinds、operators、types、summary、lambdas、iterators、stl、decompositions、ref_params）

**本文件是 Stage 1 Week 1 的硬性验收产出**：列出 cpp2lean v2 翻译器 HIR / MIR / 各 Pass 必须覆盖的所有 C++ 构造，及其初步翻译策略。任一构造在此文件未被登记即视为**设计缺口**。

---

## 1. 翻译范围与模板实例化

### 1.1 数量总览

| 类别 | 数量 | 说明 |
|---|---|---|
| `TRANSLATION_SCOPE` 函数名 | **65** | 不含 6 个已排除的死代码 |
| 实际要生成的 Lean 定义（含多实例化）| **67** | 64 个 1:1 + 1 个 `factorize` 3 实例 |
| AST node 种类 | **62** | 从最常见 DeclRefExpr(6354) 到最少的 DoStmt(2) |
| 不同运算符 | **45** | 下节详列 |
| 不同 qualType | **1201** | 模板实例化爆炸，但去重后可控 |

### 1.2 已排除的死代码（6 个）

| 函数 | 理由 |
|---|---|
| `__multivar_hensel_lift` | 经典 Wang Hensel 路线，被 MTSHL 替代，生产代码无调用 |
| `__hensel_lift_one_var` | 同上 |
| `__hensel_lc_correct` | 仅被 `__hensel_lift_one_var` 调用 |
| `__multivar_diophantine` | 仅被 `__hensel_lift_one_var` 调用 |
| `__pseudo_remainder_x1` | 仅被 `__multivar_diophantine` 调用 |
| `__taylor_coeff`（ZZ 版）| 仅被上述已排除函数调用（`__taylor_coeff_zp` 保留，Zp 版有生产调用）|

### 1.3 `factorize` 的 3 个实例化

| 实例化 | mangled qualType | 对应 Lean 定义名（提议）|
|---|---|---|
| `factorize<upolynomial_<ZZ>>` | `factorization<upolynomial_<ZZ>> (const upolynomial_<ZZ> &)` | `factorize_upoly_ir` |
| `factorize<polynomial_<ZZ, lex_<less>>>` | `factorization<polynomial_<ZZ, lex_<less>>> (const polynomial_<ZZ, lex_<less>> &)` | `factorize_lex_ir` |
| `factorize<polynomial_<ZZ, grlex_<less>>>` | `factorization<polynomial_<ZZ, grlex_<less>>> (const polynomial_<ZZ, grlex_<less>> &)` | `factorize_grlex_ir` |

**设计意义**：HIR 节点设计必须允许"同一源函数名 → 多 Lean 定义"。其他 64 个 TRANSLATION_SCOPE 函数都是 1:1 映射，此为唯一例外。

---

## 2. 类型系统

### 2.1 基础数值类型

| C++ | Lean | 使用频次（约）| 备注 |
|---|---|---|---|
| `uint64_t` | `UInt64` | 高 | 系数、度数、哈希 |
| `int64_t` | `Int64` | 中 | 普通整型（配对用 int 更多）|
| `int` | `Int` | 高 | 索引、循环计数器 |
| `uint32_t` | `UInt32` | 低 | 特定字段 |
| `int32_t` | `Int32` | 低 | 特定字段 |
| `size_t` | `Nat` or `USize` | 12 | 容器大小 |
| `bool` | `Bool` | 高 | 逻辑值 |
| `double` | `Float` | 3 | 仅 `__heuristic_starting_precision` 用 `std::log`、`std::ceil` |

### 2.2 CLPoly 核心类型

| C++ | Lean | 使用频次 | 备注 |
|---|---|---|---|
| `clpoly::ZZ` | `ZZ`（已在 `clpoly_model.lean`）| 高 | 任意精度整数（GMP wrapper）|
| `clpoly::QQ` | `QQ` | 中 | 有理数（两个 ZZ）|
| `clpoly::Zp` | `Zp`（结构体：`val : UInt64`、`prime : UInt64`）| 高 | 素数模 p 元素 |
| `clpoly::variable` | `Variable` | 中 | 变量名占位 |
| `clpoly::basic_monomial<order>` | `Monomial`（按序号序列化） | 中 | 多变量单项式 |
| `clpoly::polynomial_<T, order>` | `MvPoly T`（统一抽象）| 高 | 多变量多项式 |
| `clpoly::upolynomial_<T>` | `SparsePoly T` | 高 | 单变量多项式 |
| `clpoly::factorization<Poly>` | `Factorization Poly` | 中 | `{content, factors : Array (Poly × Nat)}` |

### 2.3 STL 容器（已覆盖）

| C++ | Lean | 频次 | 设计要点 |
|---|---|---|---|
| `std::vector<T>` | `Array T` | **2066** | `.push_back → .push`，`.erase(it, end) → .take n` |
| `std::pair<A,B>` | `A × B` | 685 | 原生 product |
| `std::tuple_element` | （不直接翻译）| 235 | 仅类型特征，SFINAE 辅助 |
| `std::map<K,V>` | `StdMap K V`（自定义）| 94 | 按 key 排序 array + find/insert |
| `std::list<T>` | `Array T` | 34 | 不需要 linked list 语义 |
| `std::set<K>` | `StdSet K` / `Array K` | 6 | 仅 `__extract_monomial_content` 用 |
| `std::initializer_list` | `Array` | 3 | 列表初始化 |

### 2.4 STL 随机数三件套

| C++ | Lean 对应 | 频次 |
|---|---|---|
| `std::mt19937` | `Rng`（`clpoly_model.lean` 已定义）| 27 |
| `std::uniform_int_distribution` | 公理化：`Rng.next : Rng → Nat → (Nat × Rng)` | 6 |
| `std::random_device` | `axiom Rng.seed : Nat`（仅 `__mtshl_sparse_int` 用）| 3 |

### 2.5 引用与指针

全 192 参数分类：

| 传递方式 | 数量 | 比例 | 翻译策略 |
|---|---|---|---|
| `T` (value) | 40 | 21% | 不变 |
| `const T&` | 125 | 65% | 不变（Lean 值语义，引用透明）|
| `T&` (non-const ref) | 26 | **14%** | **HIR `ref_elim` 转为返回值 tuple** |
| `T*` (pointer) | 1 | 0.5% | 仅 1 处，个案处理 |
| `T&&` (rvalue ref) | 0 | — | 不用 |

**输出参数汇总**：26 个 non-const ref 参数分布在 **22 个函数**中。8 个函数的旧 `TRANSLATION_SCOPE_OUTPUT_PARAMS` 配置**错位**，v2 的 `ref_elim` Pass 从 AST 自动推导，不再依赖手工表。

---

## 3. 控制流

全覆盖且范围可控：

| 构造 | AST kind | Count | HIR/MIR 策略 |
|---|---|---|---|
| if-else | `IfStmt` | 275 | MIR 原样保留，phi 在 SSA build 阶段产生 |
| while 循环 | `WhileStmt` | 20 | MIR `loop_lower` 提取为 `partial def` 尾递归 |
| for 循环 | `ForStmt` | 145 | 同上 |
| range-for | `CXXForRangeStmt` | **92** | `iter_recognize` Pass 识别为高阶操作（详 §5）|
| do-while | `DoStmt` | **2** | 仅 `__wang_core`、`__zassenhaus_recombine` — `loop_lower` 特殊处理 |
| break | `BreakStmt` | 26 | `loop_lower` 下降为 `_break_flag + tail-return` |
| continue | `ContinueStmt` | 52 | `loop_lower` 下降为"跳过 body 尾部 + 递归" |
| return（含循环内）| `ReturnStmt` | 165 | 循环内 return 用 `_ret_flag + _ret_val` 机制 |
| 三元运算 | `ConditionalOperator` | 39 | 直接映射为 `if-then-else` |

**未用到的控制流**（已确认）：`switch` / `case` / `default` / `goto` / `try` / `throw` / `catch` — 简化 HIR 语义。

### 3.1 结构化绑定范围

30 处 `DecompositionDecl`，其中 25 处在 range-for 内（std::map 或 std::vector<std::pair> 遍历），5 处独立。全部绑定 `std::pair`（无 tuple 或 custom struct）。

典型模式：`[mono, coeff]`、`[var, deg]`、`[fi, ei]`、`[gk, mk]`、`[fac, mult]`（5 次为同一 pattern）

---

## 4. Lambda 全清单

### 4.1 概览

26 个 in-scope lambda 分布在 14 个函数里：

| 宿主 | lambda 数 |
|---|---|
| `__lll_reduce` | 5 |
| `__wang_core` | 3 |
| `__factor_multivar` / `__mtshl_sparse_int` / `__mtshl_step_j` / `__vanhoeij_recombine` / `__wang_leading_coeff` / `__zassenhaus_recombine` | 各 2 |
| 其他 7 函数 | 各 1 |

### 4.2 按捕获模式分布

| 模式 | 数量 | 说明 |
|---|---|---|
| `[]`（无捕获）| 13 | std::sort 比较器 / 纯函数谓词 |
| `[&]`（默认引用）| 10 | 修改外层变量，`lambda_lift` 要 capture 加入返回值 tuple |
| `[x]`（具名 by-value）| 3 | `[j]`、`[p]`、`[use_large_prime]` 各 1 |
| `[=]`（默认值）/ `[&, x]`/`[=, &x]` 混合 | 0 | **不使用**，简化 `lambda_lift` |

### 4.3 按 body 大小

| 大小 | 数量 | 翻译难度 |
|---|---|---|
| 1 行 | 2 | trivial |
| 2-3 行 | 6 | 容易 |
| 4-10 行 | 12 | 中等 |
| 11-30 行 | 6 | 复杂（`__lll_reduce` 的 row_sub、row_swap；`__mtshl_step_j` 的 lc_correct、product_F）|

### 4.4 HIR `lambda_lift` Pass 设计要点

1. **Lambda 统一提升**：每 lambda 生成独立 `_lambda_N_ir` 定义
2. **`[&]` 捕获**：调用点改为显式传参；修改后的变量作为返回 tuple 的后半部分
3. **`[x]` 具名**：同 `[&]` 处理，但语义上是 by-value（lambda 内部不修改）
4. **Generic lambda**：**源已确保无**（5 处改为显式类型）
5. **lambda 作为 std::sort 比较器**：特殊识别为 `Array.qsortWith cmp_fn`

---

## 5. 迭代器模式

### 5.1 优先级排序（翻译难度）

| 模式 | 数量 | 策略 |
|---|---|---|
| Range-for `for (auto& x : v)` | 92 | 直接 HIR for 节点保留，`iter_recognize` 后转 `v.foldl`/`v.foreach` |
| Structured-binding range-for `for (auto& [k, v] : m)` | 24 | 对 map 遍历专用分支：`m.foldl (fun acc (k, v) => ...)` |
| **Compact-erase 双指针** | **4** | **关键**：模式匹配为 `v.filter predicate` |
| 并行双迭代器（zip-walk）| 1 | `__upoly_divmod_mod`，无抽象模式，保留为显式 loop |
| Classic iterator loop | 1 | 罕见，同上 |

### 5.2 Compact-erase 4 处定位

| 宿主 | 文件:行 | 对象 |
|---|---|---|
| `__upoly_mod_coeff` | univar.hh:193-205 | `f.data()` — 模系数后去零项 |
| `__hensel_step` | univar.hh:423-429 + 469-475 | `e.data()` / `ep.data()` — Hensel 步后 compact（2 处）|
| `__hensel_step_linear` | univar.hh:649-656 | `e.data()` — 线性 Hensel 后 compact |

---

## 6. STL 算法使用（6 种，90 次调用）

| 算法 | 次数 | Lean 对应 | 备注 |
|---|---|---|---|
| `std::move` | 63 | `id` | 值语义下是 noop |
| `std::sort` | 8 | `Array.qsortWith cmp` | 5 处用 lambda 比较器（已改具体类型）|
| `std::max` | 7 | `max` | Mathlib/stdlib |
| `std::iota` | 5 | `Array.range` | 生成 `[start, start+1, ..., end-1]` |
| `std::swap` | 4 | 值重绑 | `(b, a) := (a, b)`，Lean 无 in-place swap |
| `std::min` | 3 | `min` | 同 max |

未使用：`std::find*`、`std::remove*`、`std::unique`、`std::copy*`、`std::reverse`、`std::rotate`、`std::for_each` 等。

---

## 7. 运算符（45 种）

### 7.1 `CXXOperatorCallExpr`（用户类型重载，共 ~1332 次）

| Operator | 次数 | Lean 映射（示例：Zp/ZZ）|
|---|---|---|
| `operator[]` | 468 | `a[i]!`（Array）/ `m.find k` (StdMap) |
| `operator*` | 207 | `a * b` or `Zp.mul a b` |
| `operator=` | 167 | SSA 重绑 |
| `operator!=` | 129 | `!(a == b)` |
| `operator++` | 110 | 迭代器推进 / UInt64 递增 |
| `operator()` | 40 | lambda 调用 / 构造函数（需区分）|
| `operator-` (二元) | 38 | 减法 |
| `operator+` | 30 | 加法 |
| `operator->` | 29 | 迭代器解引用 |
| `operator*=` | 23 | `a := a * b`（SSA） |
| `operator==` | 23 | 等比较 |
| `operator-=` | 12 | `a := a - b` |
| `operator/` | 11 | 除法 |
| `operator<` | 11 | 比较 |
| `operator%` | 8 | 取模 |
| `operator<=` | 7 | 比较 |
| `operator+=` | 5 | `a := a + b` |
| `operator>` | 5 | 比较 |
| `operator<<` | 4 | 左移（或流输出）|
| `operator/=` / `operator>=` / `operator<<=` | 各 1-2 | |

### 7.2 `BinaryOperator`（基础类型，~535 次）

| Operator | 次数 | Lean |
|---|---|---|
| `<` / `=` / `-` / `==` / `+` / `>` | 146 / 74 / 62 / 40 / 40 / 37 | 对应 Lean 原生 |
| `&&` / `\|\|` | 34 / 9 | `Bool.and` / `Bool.or` |
| `>=` / `<=` / `!=` / `*` / `/` / `%` | 24 / 19 / 18 / 16 / 11 / 5 | 原生 |

### 7.3 `UnaryOperator`（~260 次）

| Operator | 次数 | Lean |
|---|---|---|
| `++` | 144 | `i := i + 1`（SSA）|
| `!` | 69 | `!` |
| `__extension__` | 23 | GCC 扩展（含 signed-integer 操作），视作 noop |
| `-` | 12 | 取负 |
| `--` | 12 | `i := i - 1` |

### 7.4 `CompoundAssignOperator`（基础类型 `+=` 等，共 6 次）

| Operator | 次数 | Lean |
|---|---|---|
| `*=` | 2 | `a := a * b` |
| `+=` | 2 | `a := a + b` |
| `/=` / `-=` | 各 1 | |

**注意**：绝大多数 compound assign 实际是用户类型（`operator*=`/`operator+=` 等，见 §7.1），基础类型的 CompoundAssignOperator 仅 6 次。

---

## 8. 类型转换（Cast）

| Kind | Count | 策略 |
|---|---|---|
| `ImplicitCastExpr` | 5419 | 最多，Clang 已标注 source/target 类型 → 查 `CAST_TABLE` |
| `CXXConstructExpr` | 627 | 构造函数调用，`operator_resolve` 映射到 Lean 构造器 |
| `CXXFunctionalCastExpr` | 155 | `Type(value)` 显式转换 → 同 ImplicitCast |
| `CStyleCastExpr` | 70 | `(Type)value` → 同上 |
| `CXXStaticCastExpr` | 26 | `static_cast<T>(...)` → 同上 |

**P1 原则**：翻译器**不推断类型**，全部从 AST 读取（参见 `blueprint.md` §2b）。

---

## 9. 内存与生命周期（全部确认不使用）

| 构造 | Count | 结论 |
|---|---|---|
| `CXXNewExpr` | 0 | **无动态分配** |
| `CXXDeleteExpr` | 0 | 无 |
| `CXXThisExpr` | 0 | 无显式 `this`（方法调用通过 MemberExpr 自动携带）|
| `CXXTryStmt` / `CXXThrowExpr` / `CXXCatchStmt` | 0 / 0 / 0 | **无异常机制** |
| `GotoStmt` | 0 | 无 |
| `CXXBindTemporaryExpr` | 542 | 临时对象绑定，SSA 化时变为 let 绑定 |
| `MaterializeTemporaryExpr` | 522 | 临时物化，同上 |
| `ExprWithCleanups` | 535 | 需要运行析构的表达式，Lean 值语义下可忽略 |

---

## 10. 模板相关 AST（受控范围内）

| Kind | Count | 说明 |
|---|---|---|
| `FunctionTemplateDecl` | 9 | **仅在 Lambda 闭包内出现**（operator() 原型 + 转换运算符），不用作翻译单元 |
| `TemplateTypeParmDecl` | 18 | 同上（Lambda 闭包内的 auto 参数）|
| `TemplateArgument` | 145 | 模板参数具现化信息，用于 `SubstTemplateTypeParmType` 替换 |
| `TemplateSpecializationType` | 71 | 模板类型具现化（如 `polynomial_<ZZ, lex>`）|
| `SubstTemplateTypeParmType` | 68 | 已替换的模板类型参数（Clang 保留记录）|
| `TypeAliasDecl` | 20 | `using Poly = ...` 等内部 typedef |

**关键**：在经过 Day 1 的 3 项修复（survey 脚本 bug + `__taylor_coeff` 死代码 + 5 处 generic lambda → 显式类型）后，**本设计所需处理的 C++ 函数完全单态化**（`CXXDependentScopeMemberExpr` / `UnresolvedLookupExpr` / `DependentScopeDeclRefExpr` 总数为 0）。模板内的 FunctionTemplateDecl/TemplateTypeParmDecl 仅残留在 Lambda closure 的内部，不影响主体翻译。

---

## 11. 函数调用与成员访问

| Kind | Count | 说明 |
|---|---|---|
| `DeclRefExpr` | 6354 | 变量/函数引用，最常见 |
| `MemberExpr` | 1177 | 成员访问 `.`/`->` |
| `CXXMemberCallExpr` | 869 | `obj.method(...)` |
| `CallExpr` | 507 | 普通函数调用 |
| `VarDecl` | 1145 | 变量声明（含 let-binding 源）|
| `ParmVarDecl` | 267 | 参数声明 |
| `DeclStmt` | 1123 | 声明语句容器 |
| `CompoundStmt` | 381 | `{...}` 块 |

**SSA 所涉**：每个 `VarDecl`（声明 + 可能初始化）+ 每个 `BinaryOperator::=` 赋值 + 每个 compound assign / `++` / `--` + 每个被 ref 参数调用（`operator_resolve` 识别 ref callee 后）都是一个 mutation 点，需要 `ssa_build` 递增版本号。

---

## 12. 字面量

| Kind | Count | Lean |
|---|---|---|
| `IntegerLiteral` | 596 | `(N : UInt64/Int/...)`（类型从 AST）|
| `CXXBoolLiteralExpr` | 82 | `true` / `false` |
| `StringLiteral` | 46 | `String`（日志用，大概率可忽略）|
| `SourceLocExpr` | 46 | `__builtin_source_location()` — assert 消息用，可忽略 |
| `PredefinedExpr` | 23 | `__PRETTY_FUNCTION__` 等 — assert 消息用 |
| `FloatingLiteral` | 3 | `Float`（`__heuristic_starting_precision` 的 `2.5`、`1.0`、`2.0`）|

---

## 13. 各 Pass 的 AST 构造覆盖表

| Pass | 必须处理的 AST kind |
|---|---|
| `parse` | 全部 62 种（AST → HIR 直接映射）|
| `ref_elim` | `ParmVarDecl`（ref 修饰）+ `ReturnStmt` + `CallExpr`（调用方）|
| `lambda_lift` | `LambdaExpr` + `CXXRecordDecl`（闭包类）+ `CXXMethodDecl`（operator()）|
| `iter_recognize` | `CXXForRangeStmt` + iterator init 模式（源码辅助）|
| `operator_resolve` | `CXXOperatorCallExpr` + `ImplicitCastExpr` + `CXXConstructExpr` + `CXXMemberCallExpr` |
| `ssa_build` | `VarDecl` + `BinaryOperator::=` + Compound assign + `UnaryOperator::++/--` + `IfStmt`（phi 放置）|
| `loop_lower` | `ForStmt` / `WhileStmt` / `DoStmt` / `CXXForRangeStmt` + `BreakStmt` + `ContinueStmt` + `ReturnStmt`（循环内）|
| `codegen` | 全部 HIR/MIR 节点 → Lean 字符串 |

---

## 14. 已知设计决策

1. **模板实例化**：把 `factorize` 的 3 实例化生成为 3 个独立 Lean 定义，其他 64 函数各 1 定义
2. **输出参数**：从 AST 的 ParmVarDecl 类型自动识别 non-const ref，`ref_elim` Pass 统一改为 tuple 返回值
3. **Lambda**：全部提升为 `_lambda_N_ir` 独立定义；`[&]` capture 作参数，修改值作返回 tuple
4. **迭代器 compact 模式**：4 处双指针模式识别为 `Array.filter`
5. **STL 容器**：vector→Array、pair→×、map→StdMap（自定义）、list→Array
6. **RNG**：axiom `Rng.next` 形式，与已有 L2 一致
7. **除法/乘法/GCD 等多项式原语**：信任 `clpoly_model.lean` 声明 + Mathlib 等价，不展开逐行翻译（翻译时查 FUNC_MAP）
8. **类型转换**：全部 cast kind 映射到 `CAST_TABLE` 查表，不做推断
9. **异常/Goto/动态分配**：零出现，HIR 不包含对应节点
10. **控制流下降**：所有循环统一 `partial def` 尾递归，break/continue/return 用 flag 机制（Aeneas 风格）

---

## 15. 设计缺口（需 Stage 1 Week 3 详细讨论）

- [ ] **`operator()` 区分**：40 次 `CXXOperatorCallExpr::operator()` 里哪些是 lambda 调用、哪些是 RNG `dist(rng)`、哪些是用户自定 functor？需要在 `operator_resolve` Pass 里判别
- [ ] **`operator->` 29 次**：只在迭代器路径出现？需确认
- [ ] **`__extension__` 23 次**：GCC 扩展，当 noop 处理还是精确建模？
- [ ] **`CXXBindTemporaryExpr` / `MaterializeTemporaryExpr`**：在哪些场景必须保留（析构副作用）？CLPoly 的 C++ 子集里是否全部可 SSA 化为 let？
- [ ] **`StringLiteral`/`SourceLocExpr`/`PredefinedExpr`**：assert 消息用途，翻译器是否完全丢弃（让 `require` 不带消息）？
- [ ] **`do-while`（2 处）**：`loop_lower` 的 `while-loop` 变体处理方式
- [ ] **`initializer_list`（3 处）**：展开为数组字面量？
- [ ] **`CXXConstructExpr`（627 次）**：哪些是真正的对象构造（有副作用）、哪些是类型转换（应按 cast 处理）？

这些缺口由 Week 3 `cpp-subset-semantics.md` 的详细语义章节覆盖。

---

## 16. 与 Week 2、3、4、5 的衔接

- **Week 2**（类型系统清查）：细化本文件 §2 的映射表，特别是模板实例化策略与 `StdMap`/`Factorization` 的 Lean 结构定义
- **Week 3**（C++ 子集形式语义）：为 §3-7 的每个构造写 denotational semantics
- **Week 4**（HIR 设计）：基于 §13 的 Pass-AST 表，定义 HIR 数据结构和 Pass 1-5 规格
- **Week 5**（MIR + codegen 设计）：同上，Pass 6-8

**Stage 1 Week 1 验收**：本文件完成即达标。每个 §13 列出的 AST kind 都能在 §3-12 找到对应条目与翻译策略。
