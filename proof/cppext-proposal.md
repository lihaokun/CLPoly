# C++ 到 Lean4 翻译拓展提案

> **目标读者**：C→Lean 翻译工具的开发者，需要将其工具拓展为支持 C++，覆盖 CLPoly 项目。
> **分析对象**：CLPoly C++ 多项式计算库（~17.6K 行 C++，~25K Clang AST 节点，66 个翻译目标函数）。
> **产物对接目标**：生成的 Lean 代码需无缝接入现有 L2 算法模型精化层（`CLPoly/Refinement/`）。
> **参考实现**：本项目 `proof/cpp2lean_v2/` 已有 C++→Lean 翻译器（Clang AST → Lean4），可作为特性映射的参考。
> **项目仓库**：<https://github.com/lihaokun/CLPoly>
> **文档位置**：`proof/cppext-proposal.md`

---

## 目录

1. [项目背景与验证架构](#1-项目背景与验证架构)
2. [代码全景：规模与构成](#2-代码全景规模与构成)
3. [C++ 特性全景清单](#3-c-特性全景清单)
4. [C→Lean 与 C++→Lean 的差距分析](#4-cleanc-lean-的差距分析)
5. [C++ 特性翻译策略（核心）](#5-c-特性翻译策略核心)
6. [生成的 Lean 代码与 L2 的对接设计](#6-生成的-lean-代码与-l2-的对接设计)
7. [分阶段实施路线图](#7-分阶段实施路线图)
8. [关键挑战与风险](#8-关键挑战与风险)
9. [附录](#9-附录)

---

## 1. 项目背景与验证架构

### 1.1 CLPoly 是什么

CLPoly 是一个多变量多项式计算 C++ 库，支持 Z[x₁,…,xₙ] 和 Q[x₁,…,xₙ] 上的算术、GCD、因式分解、结式、特征列、Gröbner 基、实根隔离。泛型模板设计，支持不同单项式序和系数域。

### 1.2 三层验证架构

```
C++ 实现 (clpoly/ — 17.6K 行)
    ↑  翻译器 + 人工审核，确保语义对应
L1: 实现模型 (Model.lean + Generated/Corpus.lean — 1:1 with C++)
    ↑  精化证明 (CLPoly/Refinement/*.lean)
L2: 算法模型 (Algorithm/*.lean — Mathlib 多项式类型, 抽象数据结构)
    ↑  正确性证明
L3: 数学基础 (Math/*.lean — 有限域理论,  Hensel 引理, 不可约判定)
    ↑
Lean 4 定理证明器 + Mathlib4
```

- **L3 + L2**: 已完成（0 sorry, `lake build` 3071/3071 全绿）
- **L1 Model.lean**: 已完成（含 Zp、ZZ、QQ、SparsePolyZp、SparsePolyZZ、MvPolyZZ 等核心类型定义）
- **L1 Corpus.lean**: 翻译器产物（C++ 函数的 1:1 Lean 翻译）
- **L1 Refinement**: 精化证明层，正在开发中（证明 Corpus.lean 中函数的行为等价于 Algorithm/*.lean 中的算法）

### 1.3 翻译产物在验证链中的位置

翻译器将 C++ 函数翻译为 Lean 的 `partial def`（带 `_ir` 后缀），这些定义在 `Corpus.lean` 中。精化层 `Refinement/*.lean` 提供形式化定理，将每个 `_ir` 函数与其在 `Algorithm/*.lean` 中的数学模型关联起来：

```
Corpus.lean                    Algorithm/*.lean
  ddf_Zp_ir (f : SparsePolyZp)   DDF.ddf (f : Polynomial (ZMod p))
     ↓  refinement theorem          ↓
  Refinement/DDF.lean:
    theorem ddf_Zp_ir_refines (hf : ...) :
      toSparsePolyZp (DDF.ddf (fromSparsePolyZp f)) = ddf_Zp_ir f
```

因此，翻译器不只需要生成语法正确的 Lean 代码，还需要：
1. **类型兼容**：函数签名使用 `Model.lean` 中定好的类型（`Zp`, `SparsePolyZp`, `MvPolyZZ` 等），而非 C++ 原始类型
2. **命名一致**：函数名符合 `_ir` 后缀约定，便于精化层通过 `REFINEMENT_MAP` 定位
3. **结构透明**：控制流结构保持直接可读，便于手写精化证明

### 1.4 推荐翻译流水线

```
Clang AST (JSON)
  ↓ Phase 1: parse → HIR（原始 AST 映射，56 种节点）
  ↓ Phase 2: ref_elim → HIR₁（消除 T& 参数，转为 tuple 返回值）
  ↓ Phase 3: lambda_lift → HIR₂（lambda 提升 + 捕获分析）
  ↓ Phase 4: iter_recognize → HIR₃（迭代器 → 高阶函数/索引）
  ↓ Phase 5: operator_resolve → HIR₄（运算符/方法/构造器 → 已知函数）
  ↓ Phase 6: ssa_build → MIR（SSA + CFG, Cytron 算法）
  ↓ Phase 7: loop_lower → MIR'（循环 → partial def 尾递归）
  ↓ Phase 8: codegen → Lean 4 代码（使用 Model 类型, _ir 后缀）
  ↓ Phase 9: refine_gen → 精化定理骨架（对接 Algorithm/*.lean）
```

---

## 2. 代码全景：规模与构成

### 2.1 源文件统计

| 类别 | 文件数 | 行数（约） |
|------|--------|-----------|
| 核心库头文件 `.hh` | 27 | 13,200 |
| 核心库实现 `.cc` | 7 | 1,300 |
| 辅助库（number） | 2 个 `.hh` | 1,100 |
| 解析库 | 1 个 `.hh` | 3,000+ |
| **总计** | 34 个头文件 + 7 源文件 | **~17,646** |

### 2.2 翻译目标规模（TRANSLATION_SCOPE）

- **66 个函数**（67 个 Lean 定义，`factorize` 有 3 种实例化）
- **~25,059 AST 节点**（累计）
- **56 种不同 AST 节点类型**
- **44 种不同运算符**
- **1,199 种不同 qualType**

### 2.3 控制流分布

| 构造 | 出现次数 | 翻译策略 |
|------|---------|---------|
| if-else | 282 | MIR 保留 + phi 节点 |
| for 循环 | 146 | partial def 尾递归 |
| range-for | 92 | desugar → 索引循环 / filter |
| while | 20 | partial def 尾递归 |
| do-while | 2 | 展开为 while |
| break | 27 | flag + tail-return |
| continue | 52 | skip body + recurse |
| return | 168 | ReturnTerm |
| 三元表达式 | 39 | if-then-else |

**已确认零出现**（可无需支持）：`switch`/`case`/`default`/`goto`/`try`/`catch`/`throw`/`CXXNewExpr`/`CXXDeleteExpr`/`CXXThisExpr`。

---

## 3. C++ 特性全景清单

所有 C++ 特性按"是否存在于 C"和"使用频率"分类。

### 3.1 与 C 共享的特性

| 特性 | C 等价物 | CLPoly 使用频次 | 备注 |
|------|----------|----------------|------|
| `if`/`else` | 相同 | 282 | |
| `for` | 相同 | 146 | |
| `while` | 相同 | 20 | |
| `do-while` | 相同 | 2 | |
| `break`/`continue` | 相同 | 79 | |
| `return` | 相同 | 168 | |
| `? :` 三元 | 相同 | 39 | |
| `int64_t`/`uint64_t` | `<stdint.h>` | 全库 | |
| `bool` | `_Bool` | 全库 | |
| `double` | 相同 | 仅 `__heuristic_starting_precision` | |
| `sizeof` | 相同 | `ZZ.hh` 平台检测 | |
| `assert` | `<assert.h>` | ~50+ | 含 `NDEBUG` 保护 |
| `#define`/`#ifdef` | 相同 | 条件编译 | |
| `nullptr` | `NULL` | 43 | |
| 整数算术 | 相同 | 全库 | |
| 指针 | 相同 | VHC 堆结构、GMP | |

### 3.2 C++ 独有特性（按翻译难度分组）

#### 组 A：核心不可绕过（高频率，必须优先支持）

| 特性 | CLPoly 使用频率 | 产生翻译难点 |
|------|----------------|-------------|
| **模板**（class/function template） | **全库骨架**，全部 .hh 均为模板 | 实例化选择、类型映射 |
| **类（class）** | 10+ 核心类 | 构造函数、成员函数、字段访问 |
| **运算符重载** | 45 种运算符，~1,332 次 `CXXOperatorCallExpr` | 解析到类型类/trait |
| **引用 `T&` / `const T&`** | 192 个参数，65% const ref，14% 非 const ref | ref-elim pass、别名分析 |
| **`std::vector`** | 2,066 次使用，核心数据结构 | → `Array T` |
| **`std::pair`** | 685 次使用，多项式项表示 | → `A × B` |
| **范围 for 循环** | 92 次 | desugar 为索引循环 |
| **auto** | 全库普遍 | 类型从 Clang AST 传播 |
| **constexpr** | 147 次 | 展开为运行时值 |

#### 组 B：重要但可组合处理

| 特性 | CLPoly 使用频率 | 翻译策略 |
|------|----------------|---------|
| **移动语义（`std::move`）** | 63+ 次 | 语义等同拷贝，去除 move |
| **`std::map`** | 94 次 | → `StdMap K V` |
| **Lambda 表达式** | 26 次（14 个函数中） | 提升为独立函数 + 捕获参数 |
| **结构化绑定** | 30 次 | 展开为 `let` 链 |
| **`if constexpr`** | 6 次 | 编译期已定，单分支保留 |
| **`using` 别名** | 42 次 | 类型展开 |
| **`std::sort`** | 8 次 | → `Array.qsort` |
| **`noexcept`** | 8 次 | 可忽略 |

#### 组 C：低频率但存在

| 特性 | 出现位置 | 次数 |
|------|---------|------|
| `std::initializer_list` | 构造函数/赋值 | 5 |
| `std::list` | 变量跟踪 | 34 |
| `std::set` | 图算法 | 6 |
| `__builtin_clzll` 等 GCC 内建 | 算术优化 | 9 |
| `__int128` / `unsigned __int128` | Barrett 乘、溢出检测 | 6 |
| 内联汇编（`__asm__`） | `dense_upoly_zp.hh` 进位加法 | 2 |
| `friend` 函数 | 运算符重载、输出 | 69 |
| `new[]`/`delete[]` | VHC 堆结构 | 多处 |
| `throw` | 错误报告 | 15+ |
| `explicit` | 构造函数 | 8 |
| `static` 局部/成员 | 全局常量、注册表 | 68 |
| `mutable` | Ideal 类缓存 | 3 |
| `std::atomic` | profile 统计 | 微量 |

#### 组 D：本项目不使用（无需支持）

继承/虚函数/多态、`override`/`final`、`enum class`、智能指针、`std::variant`/`optional`/`any`、`std::function`、C++20 concepts、`static_assert`、`decltype`、RTTI、`volatile`、`std::thread`/`mutex`、`extern template`、协程、`goto`、`switch`/`case`。

---

## 4. C→Lean 与 C++→Lean 的差距分析

### 4.1 AST 节点差异

C 翻译工具应已处理 C 的 AST 节点。C++ 新增主要节点如下：

| AST 节点 | C 中 | C++ 中频次 | 差异性质 |
|----------|------|-----------|---------|
| `CXXOperatorCallExpr` | ✗ | 1,347 | **全新**：运算符重载调用 |
| `CXXMemberCallExpr` | ✗ | 878 | **全新**：方法调用 |
| `CXXConstructExpr` | ✗ | 639 | **全新**：构造函数调用 |
| `MemberExpr` | ✗（仅 struct 字段 `.`） | 1,183 | **语义扩展**：含方法、嵌套类型 |
| `TemplateSpecializationType` | ✗ | 101 | **全新**：模板实例化类型 |
| `SubstTemplateTypeParmType` | ✗ | 100 | **全新**：模板参数替换 |
| `CXXForRangeStmt` | ✗ | 92 | **全新**：范围 for |
| `DecompositionDecl` | ✗ | 30 | **全新**：结构化绑定 |
| `LambdaExpr` | ✗ | 24 | **全新**：Lambda |
| `CXXTemporaryObjectExpr` | ✗ | 46 | **全新**：临时对象 |
| `CXXFunctionalCastExpr` | ✗ | 133 | **语义扩展**：函数式转型 |
| `CXXBindTemporaryExpr` | ✗ | 558 | **全新**：临时对象绑定，可 strip |
| `MaterializeTemporaryExpr` | ✗ | 536 | **全新**：临时物化，可 strip |
| `ExprWithCleanups` | ✗ | 530 | **全新**：析构清理，可 strip |

### 4.2 类型系统差距

| 维度 | C | C++ | 对翻译的影响 |
|------|---|-----|------------|
| 自定义类型 | `struct`（仅数据） | `class`（方法/构造/访问控制） | 方法需独立函数派发 |
| 泛型 | 无 | 模板 | 实例化后为具体类型，需类型名解析 |
| 引用 | 仅指针 | `T&`, `const T&`, `T&&` | ref-elim pass 消除引用参数 |
| 命名空间 | 无 | `namespace clpoly` | strip 前缀 |
| 容器 | 数组 | `vector`, `pair`, `map` | 映射到 Array/Pair/StdMap |

### 4.3 C 工具已有但需调整的功能

| C→Lean 已有功能 | C++ 需要调整 |
|-----------------|-------------|
| `DeclRefExpr` → `Var` | 兼容，无需调整 |
| `IntegerLiteral` → 字面量 | 新增 `CXXBoolLiteralExpr` |
| `BinaryOperator` → `BinOp` | 兼容 |
| `UnaryOperator` → `UnaryOp` | C++ 前缀 `++`/`--` 需展开为 `x := x + 1` |
| `ImplicitCastExpr` → `Cast` | C++ 更多 cast 类型需扩展 cast table |
| `CallExpr` → `Call` | 新增 `CXXOperatorCallExpr`/`CXXMemberCallExpr`/`CXXConstructExpr` |
| `IfStmt` → `CondJump` | 兼容 |
| `ForStmt`/`WhileStmt`/`DoStmt` → 循环提取 | 新增 `CXXForRangeStmt` |
| `DeclStmt` + `VarDecl` → 变量声明 | C++ 变量可含构造函数调用 |

---

## 5. C++ 特性翻译策略（核心）

每个特性标注：**出现频率** / **语义差异** / **Lean 映射** / **复杂度** ★→★★★★★ / **前置依赖**。

### 5.1 模板（Templates）

**频率**：全库骨架，27+ `.hh` 全部为模板。

**语义差异**：C 无等价物。CLPoly 使用模板实现泛型多态（替代 OOP 继承）。

**关键简化**：翻译器不需要实现模板实例化——**Clang AST dump 已包含完全实例化后的代码**。模板参数在 AST 中以 `TemplateSpecializationType` 和 `SubstTemplateTypeParmType` 呈现，翻译器直接处理具体类型。

**多实例化**：`factorize` 有 3 种实例化（`lex`/`grlex`/`grevlex`），通过后缀区分。

**策略**：
1. 维护一个类型映射表（`CLASS_MAP`），将实例化的 C++ 类型名映射到 `Model.lean` 中的 Lean 类型
2. 维护一个函数映射表（`FUNC_MAP`），将模板函数实例化名映射到 `_ir` 后缀的 Lean 函数名
3. 同一函数的不同实例化使用不同后缀（如 `factorize_lex_ir` / `factorize_grlex_ir`）

```
C++ 模板类型            Lean 类型（Model.lean）
─────────────────────────────────────────────
Zp                   → Zp
ZZ                   → ZZ
QQ                   → QQ
basic_monomial<grlex> → Monomial Grlex
basic_monomial<lex>   → Monomial Lex
polynomial_<ZZ,grlex> → MvPolyZZ
upolynomial_<ZZ>      → SparsePolyZZ
upolynomial_<Zp>      → SparsePolyZp
factorization<Poly>   → Factorization Poly
```

**复杂度**：★★☆☆☆ — 模板已实例化，翻译器只需类型名映射。

### 5.2 类与成员函数（Classes & Member Functions）

**频率**：10+ 核心类，878 次 `CXXMemberCallExpr`，1,183 次 `MemberExpr`。

**语义差异**：C `struct` 只有字段；C++ `class` 有成员函数、访问控制。CLPoly 不使用虚函数。

**策略**：
1. `obj.method(args)` → `method_ir(obj, args)`（成员函数变为自由函数，`this` 成为首参）
2. `obj.field` → `obj.field`（字段映射）
3. `T obj(args)` 构造 → `construct_T(args)`（通过 `CONSTRUCTOR_MAP`）
4. 所有方法通过 `CLASS_MAP` 注册

**CLASS_MAP 数据结构**（需由翻译器维护）：

```python
CLASS_MAP = {
    "Zp": {
        "lean_type": "Zp",
        "constructors": {
            ("Zp", "Zp(uint64_t,uint64_t)"): "Zp.mk",
            ("Zp", "Zp()"):                  "default",
        },
        "methods": {
            "prime": ("field", "prime"),       # → obj.prime
            "inv":   ("method", "Zp.inv"),     # → Zp.inv obj
            "val":   ("field", "val"),
        },
        "operators": {
            "+": None,      # None → Lean typeclass Add
            "-": None,      # Lean typeclass Sub
            "*": None,      # Lean typeclass Mul
            "/": ("method", "Zp.div"),
        },
    },
}
```

**复杂度**：★★★☆☆ — 映射规则固定，但条目众多（~20 个方法/字段/运算符/类）。

### 5.3 运算符重载（Operator Overloading）

**频率**：45 种运算符，~1,332 次 `CXXOperatorCallExpr` + ~535 次 `BinaryOperator` + ~260 次 `UnaryOperator`。

**语义差异**：C 运算符只对内建类型有效；C++ 可重载。

**运算符 → Lean 映射**：

| 运算符 | C++ 行为（CLPoly） | Lean 映射 |
|--------|-------------------|-----------|
| `+` `-` `*` `/` | 多项式/数值类重载 | typeclass `Add`/`Sub`/`Mul`/`Div` |
| `==` `!=` `<` `>` `<=` `>=` | 比较重载 | typeclass `BEq`/`LT`/`LE` |
| `+=` `-=` `*=` `/=` | 复合赋值 | 展开为 `x := x + rhs` |
| `++` `--`（前缀） | 迭代器递增 | `x := x + 1` |
| `[]` | 多项式项访问 | `Array.get` / `Array.set!` |
| `()` | lambda/RNG/仿函数 | 按接收者类型区分 |
| `<<` `>>` | ostream/istream | I/O，翻译中忽略 |
| `->` | 迭代器解引用 | `__deref__` 模式 |

**策略**：
1. 运算符 `+` / `-` / `*` → 直接映射到 Lean 同名运算符（利用 typeclass 重载）
2. 复合赋值 `+=` → 展开为 `x := x + rhs`
3. `operator[]` → `Array.get a i`（需生成边界检查 `require`）
4. `operator()` → 检查接收者类型：lambda→提升函数；RNG→`Rng.next_advance`；仿函数→`__call__`
5. `operator<<` / `>>` → 丢弃（I/O 不在验证范围内）

**复杂度**：★★★☆☆ — 规则明确但数量大。`operator()` 需要类型敏感的上下文区分。

### 5.4 引用（References: T&, const T&, T&&）

**频率**：`const T&`: 125 参数（65%），`T&`: 26 参数（14%），`T&&`: 构造函数。

**语义差异**：C 只有指针。C++ 引用是别名、非空、不可重绑定。

**策略**：
- **`const T&`**: 等价于传值，直接映射为值参数
- **`T&`（可变输出参数）**: 需要 **ref-elim**：
  - 函数返回 `void` + N 个 ref-out → 新返回值为 `Tuple(ref_types...)`
  - 函数返回 `R` + N 个 ref-out → 新返回值为 `Tuple(R, ref_types...)`
  - 调用点 destructure

```
C++: void divmod(A a, B b, C& q, C& r);
Lean: divmod_ir (a : A) (b : B) : C × C
```

- **`T&&`（右值引用）**: 在 move 构造/赋值中出现，翻译为传值。Lean 无别名问题。

**复杂度**：★★★☆☆ — ref-elim 算法经典，需处理嵌套调用和别名链。

### 5.5 移动语义（Move Semantics）

**频率**：63+ 次 `std::move`。

**语义差异**：`std::move(x)` 只是将 `x` 转为右值引用以触发移动构造。移动后源不再使用。

**策略**：`std::move(x)` → `x`（恒等映射）。Lean 是纯函数式语言，传值即拷贝，不存在移动销毁。

```
C++: return std::move(result);         → Lean: result
C++: v.push_back(std::move(item));     → Lean: v := Array.push v item
```

**复杂度**：★☆☆☆☆

### 5.6 STL 容器与算法

#### 容器映射

| C++ 容器 | 频次 | Lean 映射 |
|----------|------|----------|
| `std::vector<T>` | 2,066 | `Array T` |
| `std::pair<A,B>` | 685 | `A × B` |
| `std::map<K,V>` | 94 | `StdMap K V` |
| `std::list<T>` | 34 | `Array T` |
| `std::set<K>` | 6 | `StdSet K` |
| `std::initializer_list` | 5 | `Array T` 字面量 |
| `std::tuple<...>` | 3 | nested `Prod` |

#### 算法映射

| C++ 算法 | 频次 | Lean 映射 |
|----------|------|----------|
| `std::move` | 63 | 恒等（见 §5.5） |
| `std::sort` | 8 | `Array.qsort` |
| `std::max`/`std::min` | 10 | `max`/`min` |
| `std::iota` | 5 | `Array.range` |
| `std::swap` | 4 | 直接展开 |
| `std::lower_bound` | 1 | `Array.findIdx?` |

#### 迭代器模式（5 种）

| 模式 | 频次 | 策略 | 示例 |
|------|------|------|------|
| 0: 范围 for | 92 | 索引化 | `for (auto& x : v)` → 索引 for |
| 1: 结构化绑定 range-for | 24 | 索引 + pair 解构 | `for (auto& [f, e] : v)` |
| 2: 紧凑擦除双指针 | 4 | `Array.filter` | 双指针→ `v.filter(pred)` |
| 3: 并行双迭代器 | 1 | 手动索引 | `__upoly_divmod_mod` |
| 4: 经典迭代器循环 | 1 | 索引化 | `while (it != v.end())` |

**复杂度**：★★★☆☆ — 模式 2 需要语句级模式识别。

### 5.7 范围 for 循环（Range-based for）

**频率**：92 次。

**策略**：desugar 为索引 for。

```
C++: for (auto& term : gz.data())
Lean: for __i : [0:gz.data.size) do
        let term := gz.data[__i]
        ...
```

**复杂度**：★★☆☆☆

### 5.8 Lambda 表达式

**频率**：26 个 lambda（14 个函数），均为排序比较器。捕获模式：`[]` 13 个 / `[&]` 10 个 / `[x]` 3 个。

**策略**：
1. **无捕获** `[]` → 提升为独立函数
2. **引用捕获** `[&]` → 提升 + 捕获变量作为函数参数，被修改者通过返回值传回
3. **按值捕获** `[x]` → 提升 + 捕获值副本作为参数

```
C++: std::sort(factors.begin(), factors.end(),
                [](auto& a, auto& b) { return degree(a.first) < degree(b.first); });
Lean: def _lambda_N_ir (a : Poly × UInt64) (b : Poly × UInt64) : Bool :=
        degree a.fst < degree b.fst
      Array.qsort factors _lambda_N_ir
```

**复杂度**：★★★☆☆

### 5.9 结构化绑定（Structured Bindings）

**频率**：30 次 `DecompositionDecl`（25 次在 range-for 中，5 次独立）。

**策略**：展开为 `let` 链。

```
C++: auto [c_g, pp_g] = __upoly_primitive(g);
Lean: let __tmp := upoly_primitive_ir g
      let c_g := __tmp.fst
      let pp_g := __tmp.snd
```

**复杂度**：★★☆☆☆

### 5.10 constexpr / if constexpr

**频率**：`constexpr` 147 次，`if constexpr` 6 次。

**策略**：
- `constexpr` 函数 → 普通 Lean 函数（标记无意义）
- `if constexpr (false)` → 删除分支
- `if constexpr (true)` → 保留分支，删除 else

**复杂度**：★☆☆☆☆

### 5.11 异常（Exceptions）

**频率**：15+ 处 `throw`，无 `try`/`catch`。

**策略**：`throw` → `require precondition`。不引入 `Except` monad。

```
C++: throw std::domain_error("division by zero");
Lean: 函数签名的 requires 中加入 h : divisor ≠ 0
```

**复杂度**：★☆☆☆☆

### 5.12 auto 类型推导

**频率**：全库普遍。

**策略**：翻译器不做类型推导——Clang AST 的每个 `VarDecl` 已标注完整 `qualType`，翻译器直接读取。

```
C++: auto ptr = F.begin();
Clang AST: VarDecl 'auto' → qualType = '__gnu_cxx::__normal_iterator<...>'
Lean: let ptr : Iterator Poly := Iterator.begin F
```

**复杂度**：★☆☆☆☆

### 5.13 namespace 与 using

**策略**：
- 在 AST parse 时 strip `clpoly::` 前缀
- `using` 类型别名展开为具体类型，不保留别名

**复杂度**：★☆☆☆☆

### 5.14 friend 声明

**频率**：69 次。

**策略**：`friend` 函数翻译为普通自由函数。Lean 无私有成员。

**复杂度**：★☆☆☆☆

### 5.15 static 成员与变量

**频率**：68 处。

**策略**：
- `static constexpr` 常量 → Lean `def`
- `static` 局部变量（GMP 一次性初始化） → 顶部 `let`
- 非常量静态成员（如 `variable` 注册表） → `axiom` 或参数化

**复杂度**：★★☆☆☆

### 5.16 explicit 构造函数与转换

**频率**：8 处。

**策略**：翻译无特殊处理（Lean 无隐式转换）。

**复杂度**：★☆☆☆☆

### 5.17 noexcept 与移动操作

**频率**：8 处（全部 move 构造/赋值 + swap）。

**策略**：忽略 `noexcept`。

**复杂度**：★☆☆☆☆

### 5.18 mutable 成员

**频率**：3 处（`Ideal` 类的缓存）。

**策略**：外化为函数参数 + 返回值：

```
C++: mutable std::vector<Poly> gb_cache_;
     mutable bool gb_computed_ = false;
Lean: def computeGroebnerBasis (gens : Array Poly) (cache : GroebnerCache)
        : GroebnerCache × Array Poly := ...
```

**复杂度**：★★★☆☆

### 5.19 std::initializer_list

**频率**：5 处声明 + 大量 `{{m, c}}` 字面量。

**策略**：翻译为 `Array` 字面量 `#[...]`。

**复杂度**：★★☆☆☆

### 5.20 GCC/Clang 扩展

**频率**：
- `__builtin_clzll`: 9 次 → `UInt64.clz`
- `__builtin_add_overflow` 等: 6 次 → 建模为 spec+pure fn
- `unsigned __int128`: 6 次 → `UInt128`（现有 Lean 类型）
- 内联汇编 `_add_carry3`: 2 处 → `axiom` + spec

```
// 内联汇编策略：受信任原语
axiom add_carry3 (a b c : UInt64) : UInt64 × UInt64 × UInt64
// axiom_spec: add_carry3 a b c = (lo, carry1, carry2)
//   where lo + carry1 * 2^64 + carry2 * 2^128 = a + b + c
```

**复杂度**：★★★★☆

### 5.21 new/delete 与手动内存管理

**频率**：VHC 算法堆分配 + GMP 包装。

**策略**：`new T[n]` → `Array.mkEmpty n`。GMP 已由 `ZZ` 类封装，翻译器不处理。

**复杂度**：★★☆☆☆

### 5.22 宏与条件编译

**频率**：`CLPOLY_PROFILE`、`NDEBUG`、`DEBUG`、`__SIZEOF_INT128__` 等。

**策略**：翻译器只处理当前平台预处理后的代码。性能计时宏丢弃。需确保 Clang AST dump 时使用与构建一致的 flag。

**复杂度**：★☆☆☆☆

### 5.23 继承与虚函数（本项目不使用）

CLPoly 零次使用 `virtual`/`override`/`dynamic_cast`/`typeid`，全部通过模板静态多态实现。不涉及。

**如果其他项目需要**（仅供扩展参考）：需虚函数表 → 代数 sum type 或显式 vtable。★★★★★ 极高复杂度。

---

## 6. 生成的 Lean 代码与 L2 的对接设计

这是翻译器设计中最关键的部分——翻译产物的最终目标是装配到现有精化层中，而不仅仅是生成语法正确的 Lean 代码。

### 6.1 当前已有基础设施

**Model.lean**（`CLPoly/Model.lean`）定义了所有 L1 类型：

| Lean 类型 | 含义 | 对应 C++ 类 |
|-----------|------|-------------|
| `Zp` | `val: UInt64`, `prime: UInt64` | `clpoly::Zp` |
| `ZZ` | 小整数优化 + GMP fallback | `clpoly::ZZ` |
| `QQ` | 分子分母约分形式 | `clpoly::QQ` |
| `Variable` | `serial: UInt64`, `name: String` | `clpoly::variable` |
| `UMonomial` | 单变量单项式 | `clpoly::umonomial` |
| `Monomial o` | 多变量单项式 | `clpoly::basic_monomial<o>` |
| `SparsePolyZp` | 稀疏 Zp 多项式 | `clpoly::upolynomial_<Zp>` |
| `SparsePolyZZ` | 稀疏 ZZ 多项式 | `clpoly::upolynomial_<ZZ>` |
| `MvPolyZZ` | 多变量多项式 | `clpoly::polynomial_<ZZ, o>` |
| `MvPolyZp` | 多变量 Zp 多项式 | `clpoly::polynomial_<Zp, o>` |
| `Factorization α` | 因式分解结果 | `clpoly::factorization<Poly>` |

**Algorithm/*.lean** 定义了 L2 算法模型，使用 Mathlib 类型（`Polynomial (ZMod p)` 等）。

**Refinement/*.lean** 当前包含精化定理骨架（全部 `sorry`），形式如：

```lean
theorem ddf_Zp_ir_refines (args...) :
    toSparsePolyZp (DDF.ddf (fromSparsePolyZp f)) = ddf_Zp_ir f
```

### 6.2 翻译产物必须满足的接口契约

#### 契约 1：类型一致性

翻译器输出的函数签名**必须使用 `Model.lean` 中已有的类型**，不能自行生成新的类型定义。

```
✓ 正确：def ddf_Zp_ir (f : SparsePolyZp) : Array (SparsePolyZp × UInt64) := ...
✗ 错误：def ddf_Zp_ir (f : MyCustomPolyType) : ... := ...
```

#### 契约 2：命名约定

- 所有翻译函数须有 `_ir` 后缀（如 `poly_gcd_ir` → 非 `poly_gcd`）
- 局部辅助函数和 lambdas 用 `_lambda_N_ir` 命名
- 构造函数用 `construct_T_ir` 命名

#### 契约 3：`require` 前置条件

所有 UB 来源（除零、越界、空容器、有符号溢出等）必须在 `require` 中声明，以便 L2 精化证明者能看到约束。

```lean
def poly_div_ir (a b : SparsePolyZZ) (h : b ≠ default) : SparsePolyZZ := ...
                                                  ^^^^^^^^^^^^^^^^
                                                  除零保护，L2 精化也需此约束
```

#### 契约 4：`partial def` vs `def`

- 含循环的函数 → `partial def`（终止性由 L2 保证）
- 无循环的纯计算函数 → `def`（全函数）

```lean
partial def ddf_Zp_ir (f : SparsePolyZp) : ... :=   -- 循环→partial
def make_zp_ir (v : UInt64) (p : UInt64) : Zp :=    -- 无循环→def
```

### 6.3 精化定理自动生成

翻译器的 Phase 9（`refine_gen`）应为每个翻译函数生成精化定理骨架：

```lean
-- 自动生成模板（入 proof/skeleton）：

theorem ddf_Zp_ir_refines (f : SparsePolyZp)
    (hf : f.WellFormed) :   -- L1 不变式
    toArray (DDF.ddf (fromSparsePolyZp f)) = ddf_Zp_ir f :=
by
  -- TODO: 精化证明
  sorry
```

生成的模板应包含：
1. **定理签名**：参数、前置条件（与 `require` 对应的 `h : ...`）
2. **等式结论**：L2 函数（经桥接转换） = L1 `_ir` 函数
3. **桥接函数调用**：`toSparsePolyZp` / `fromSparsePolyZp` 等类型转换
4. **L1 不变式**（如 `SparsePolyZp` 的 `WellFormed` 性质）：作为额外的假设

**Phase 9 的 `REFINEMENT_MAP` 结构**：

```python
REFINEMENT_MAP = {
    "ddf_Zp_ir": {
        "l2_model": "DDF.ddf",
        "bridge": {
            "to": "toSparsePolyZp",
            "from": "fromSparsePolyZp",
        },
        "theorem_shape": "map_eq",  # 5 种模板之一
        "invariants": ["WellFormed"],
    },
    # ...
}
```

### 6.4 验证方法

| 验证层次 | 方法 | 覆盖 |
|---------|------|------|
| 类型检查 | `lake build` 通过 | 100% |
| 语义等价 | B2B 测试：C++ / Lean 翻译 对同一输入运行，输出一致 | 功能级别 |
| L2 精化 | 精化定理（`refines`）填补 `sorry` 后验证 | 数学级别 |

**B2B 测试流程**：
```
C++ 程序计算 f(input) → 输出 v_cpp
Lean 翻译 f_ir(input) → 输出 v_lean
断言：v_cpp 与 v_lean 在类型映射下等价
```

目前已有 138 个 B2B 测例覆盖 Zp/ZZ/SparsePolyZp/SparsePolyZZ/MvPoly 等模块。

---

## 7. 分阶段实施路线图

### 阶段 1：基础扩展（覆盖 ~60%）

**目标**：能处理与 C 接近的 C++ 代码，翻译核心函数。

| 优先级 | 特性 | 工作量 | 前提 |
|--------|------|--------|------|
| P0 | namespace strip（`clpoly::`） | 2 天 | C→Lean AST parse |
| P0 | 类 → struct 降级 | 3 天 | C 的 struct 支持 |
| P0 | `std::vector<T>` → `Array T` | 2 天 | 类型映射扩展 |
| P0 | `std::pair<A,B>` → `A × B` | 1 天 | |
| P1 | `CXXConstructExpr` 处理 | 3 天 | 类的降级 |
| P1 | `CXXMemberCallExpr` 处理 | 5 天 | CLASS_MAP |
| P1 | 范围 for desugar | 3 天 | |
| P1 | `const T&` 传值化 | 1 天 | |
| P1 | 结构化绑定展开 | 2 天 | |
| P1 | auto 类型处理（Clang AST 已标注） | 1 天 | |

**阶段 1 交付**：可翻译 ~20 个核心函数（`__make_zp`, `__upoly_subtract_x` 等），`lake build` 通过，类型正确。

**验证**：与已有 `Corpus.lean` 中对应函数的类型签名和输出结构逐一对比。

### 阶段 2：C++ 特有语义（覆盖 ~90%）

**目标**：处理重载、引用、移动语义、Lambda，翻译全部 66 个目标函数。

| 优先级 | 特性 | 工作量 | 前提 |
|--------|------|--------|------|
| P0 | 运算符重载解析（算术/比较/下标） | 5 天 | CLASS_MAP |
| P1 | `T&` ref-elim pass | 5 天 | |
| P1 | `std::move` 消除 | 1 天 | |
| P1 | FUNC_MAP 模板实例化映射 | 3 天 | |
| P2 | lambda 提升（无捕获+引用捕获） | 5 天 | 引用处理 |
| P2 | `std::map` → `StdMap` | 2 天 | |
| P2 | `std::sort` → `Array.qsort` | 1 天 | |
| P2 | 紧凑擦除 → `Array.filter` | 5 天 | |
| P2 | `if constexpr` 分支消除 | 1 天 | |
| P2 | static 成员外部化 | 3 天 | |

**阶段 2 交付**：全部 66 个函数可通过 `lake build`。

**验证**：与 `Corpus.lean` 逐函数输出对比 + B2B 回归（138 个测例）。

### 阶段 3：高级特性与 L2 对接（覆盖 ~99%）

**目标**：处理平台扩展、手动内存、生成精化定理骨架。

| 优先级 | 特性 | 工作量 | 前提 |
|--------|------|--------|------|
| P2 | `__builtin_clzll` → `UInt64.clz` | 1 天 | |
| P2 | `__builtin_*_overflow` 映射 | 2 天 | |
| P2 | `unsigned __int128` → `UInt128` | 3 天 | |
| P3 | `explicit` / `noexcept` 忽略 | 0.5 天 | |
| P3 | `friend` 声明忽略 | 0.5 天 | |
| P3 | `std::initializer_list` → `#[]` | 2 天 | |
| P3 | `new[]` → Array 分配 | 2 天 | |
| P3 | `mutable` 成员参数化 | 3 天 | |
| P3 | Phase 9 refine_gen 实现 | 5 天 | 全部 Phase 8 |
| P4 | 内联汇编 → `axiom` | 3 天 | |
| P4 | `throw` → `require` | 2 天 | |
| P4 | `__extension__` 忽略 | 0.5 天 | |
| P4 | `std::iota` → `Array.range` | 1 天 | |

**阶段 3 交付**：
- 翻译器可处理整个 CLPoly 库
- Phase 9 为关键函数生成精化定理骨架（含 L2 函数引用、桥接调用、L1 不变式假设）
- `lake build` 全绿
- 全部 138 个 B2B 测试通过

### 总工作量

| 阶段 | 工作量 | 覆盖 |
|------|--------|------|
| 阶段 1 | ~23 人天 | ~60% |
| 阶段 2 | ~30 人天 | ~90% |
| 阶段 3 | ~25 人天 | ~99% |
| **总计** | **~78 人天** | **~99%** |

> 注：假设 C→Lean 工具已有完整 C AST 解析、类型映射、控制流翻译基础。若从零搭建则需额外增加基础工作。

---

## 8. 关键挑战与风险

### 挑战 1：`operator()` 的上下文区分

40 次 `operator()`，根据接收者类型语义不同：

| 接收者 | 语义 | Lean 映射 |
|--------|------|----------|
| `std::mt19937` | RNG | `Rng.next_advance` + 状态传递 |
| `std::uniform_int_distribution` | 分布 | 组合 RNG 调用 |
| lambda 对象 | lambda | 提升后的函数调用 |
| 自定义仿函数 | 自定义操作 | 自由函数调用 |

**方案**：按接收者类型派发，维护 `CALL_OPERATOR_MAP`。

### 挑战 2：迭代器紧凑擦除模式

4 处模式 2 循环，翻译为 `Array.filter` 需要模式匹配：

```cpp
out = v.begin();
for (i = v.begin(); i != v.end(); ++i)
    if (pred(*i)) { *out = std::move(*i); ++out; }
v.erase(out, v.end());
```

**方案**：Phase 4 中实现模式匹配。若失败可回退为显式索引循环。

### 挑战 3：内联汇编建模

`_add_carry3` 的 x86_64/ARM64 内联汇编需要受信任原语。

**方案**：`axiom + spec`。精化证明中验证 spec 与 L2 行为一致。

### 挑战 4：mutable 缓存状态

`Ideal` 类的惰性缓存（`gb_cache_`, `gb_computed_`）在纯函数模型中状态需外化。

**方案**：函数签名增加 cache 参数，通过返回值传递更新后的 cache。

### 挑战 5：精化定理的自动化

Phase 9 生成的定理骨架目前全部 `sorry`。

**方案**：
- 生成完整的定理签名（参数 + 前置条件 + 结论）
- 为每类函数形状（`map_eq`/`pair_eq`/`direct_eq` 等 5 种模板）预置证明框架
- 等完成 B2B 覆盖后，可以考虑使用 Lean 的 `rfl`/`simp` 对部分简单函数自动证明

### 挑战 6：翻译精度验证

如何确保 C++ 语义在 Lean 中完全保留？

**方案**：
1. **B2B 测试**（已建）：C++ vs Lean 翻译同一输入的输出对比
2. **类型系统验证**：`lake build` 确保类型使用正确
3. **精化定理**：L2 已知正确的算法模型作为参考标准

---

## 9. 附录

### A.1 所有 C++ 特性一览

| # | 特性 | CLPoly 使用 | C 等价 | 复杂度 | 阶段 |
|---|------|------------|--------|--------|------|
| 1 | 类模板 | 全库 | ✗ | ★★ | 1 |
| 2 | 函数模板 | 全库 | ✗ | ★★ | 1 |
| 3 | 成员函数 | 878 次 | ✗ | ★★★ | 1 |
| 4 | 运算符重载 | 1,332 次 | ✗ | ★★★ | 2 |
| 5 | `const T&` 参数 | 125 次 | ✗ | ★ | 1 |
| 6 | `T&` 输出参数 | 26 次 | ✗ | ★★★ | 2 |
| 7 | `T&&` 右值引用 | 构造函数 | ✗ | ★★ | 2 |
| 8 | `std::move` | 63 次 | ✗ | ★ | 2 |
| 9 | `std::vector<T>` | 2,066 次 | 数组 | ★★ | 1 |
| 10 | `std::pair<A,B>` | 685 次 | struct | ★ | 1 |
| 11 | `std::map<K,V>` | 94 次 | ✗ | ★★ | 2 |
| 12 | `std::sort` | 8 次 | qsort | ★★ | 2 |
| 13 | `auto` | 全库 | ✗ | ★ | 1 |
| 14 | 范围 for | 92 次 | ✗ | ★★ | 1 |
| 15 | 结构化绑定 | 30 次 | ✗ | ★★ | 1 |
| 16 | Lambda | 26 次 | ✗ | ★★★ | 2 |
| 17 | `constexpr` | 147 次 | ✗ | ★ | 1 |
| 18 | `if constexpr` | 6 次 | ✗ | ★ | 2 |
| 19 | `throw`（无 try） | 15+ 次 | ✗ | ★ | 3 |
| 20 | `explicit` | 8 次 | ✗ | ★ | 3 |
| 21 | `noexcept` | 8 次 | ✗ | ★ | 3 |
| 22 | `friend` | 69 次 | ✗ | ★ | 3 |
| 23 | `static` 成员 | 68 次 | ✗ | ★★ | 2 |
| 24 | `mutable` | 3 次 | ✗ | ★★★ | 3 |
| 25 | `using` 别名 | 42 次 | ✗ | ★ | 1 |
| 26 | `namespace` | 全库 | ✗ | ★ | 1 |
| 27 | `new[]`/`delete[]` | VHC 堆 | malloc | ★★ | 3 |
| 28 | `__builtin_*` | 9 次 | ✗ | ★★ | 3 |
| 29 | `__int128` | 6 次 | ✗ | ★★★ | 3 |
| 30 | 内联汇编 | 2 处 | ✗ | ★★★★ | 3 |
| 31 | 继承/虚函数 | **不使用** | — | N/A | — |
| 32 | `enum class` | **不使用** | — | N/A | — |
| 33 | 智能指针 | **不使用** | — | N/A | — |

### A.2 AST 节点快速参考

| 节点 | C++ 特有 | 频次 | 翻译行为 |
|------|---------|------|---------|
| `DeclRefExpr` | ✗ | 6,442 | → Var |
| `ImplicitCastExpr` | ✗ | 5,581 | → Cast（传播类型） |
| `CXXOperatorCallExpr` | ✓ | 1,347 | → Call（经 CLASS_MAP 解析） |
| `MemberExpr` | ✓ | 1,183 | → field access / method call |
| `VarDecl` | ✗ | 1,152 | → let / phi |
| `DeclStmt` | ✗ | 1,128 | → let |
| `CXXMemberCallExpr` | ✓ | 878 | → Call（经 CLASS_MAP） |
| `CXXConstructExpr` | ✓ | 639 | → construct_T（经 CONSTRUCTOR_MAP） |
| `CXXBindTemporaryExpr` | ✓ | 558 | → strip（透明包装） |
| `BinaryOperator` | ✗ | 541 | → BinOp |
| `ExprWithCleanups` | ✓ | 530 | → strip |
| `MaterializeTemporaryExpr` | ✓ | 536 | → strip |
| `IfStmt` | ✗ | 282 | → CondJump |
| `ForStmt` | ✗ | 146 | → loop_lower |
| `CXXFunctionalCastExpr` | ✓ | 133 | → Cast |
| `CXXForRangeStmt` | ✓ | 92 | → desugar |
| `CXXBoolLiteralExpr` | ✓ | 82 | → Lit |
| `LambdaExpr` | ✓ | 24 | → lambda_lift |
| `DecompositionDecl` | ✓ | 30 | → let 链 |
| `WhileStmt` | ✗ | 20 | → loop_lower |
| `DoStmt` | ✗ | 2 | → while 展开 |

### A.3 关键类型映射

| C++ 类型 | Lean 类型（Model.lean） |
|----------|-----------------------|
| `Zp` | `Zp` |
| `ZZ` | `ZZ` |
| `QQ` | `QQ` |
| `variable` | `Variable` |
| `umonomial` | `UMonomial` |
| `basic_monomial<comp>` | `Monomial comp` |
| `basic_polynomial<umonomial, Tc, uless>` | `SparsePoly Tc` |
| `basic_polynomial<basic_monomial<comp>, Tc, comp>` | `MvPoly Tc comp` |
| `factorization<Poly>` | `Factorization` |
| `dense_upoly_zp` | `DensePoly64` |
| `Ideal<Tc, comp>` | `Ideal` |

---

> **关于此文档**：本文基于 CLPoly 代码库（2026-05-26）和已有 Clang AST 调研数据（`docs/design/l1-translation-validation/survey/`）分析编写。
> 各特性使用量来自全库 grep + AST survey 统计，翻译策略参考了现有生成产物的接口约定（`_ir` 后缀、`require` 模式、REFINEMENT_MAP 结构）。
