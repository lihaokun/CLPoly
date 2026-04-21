# CLPoly C++ → Lean 类型系统映射

> 生成日期：2026-04-21
> 覆盖：65 TRANSLATION_SCOPE 函数（`factorize` 3 实例 → 67 Lean 定义）
> 基于：`cpp-construct-catalog.md`、`numeric-types.md`、`casts.md`、`stl.md`、`clpoly-model-inventory.md`、`lean-stdlib-catalog.md`

**Stage 1 Week 2 的硬性验收产出**：给每种 C++ 类型、每种转换、每种 STL 原语定义权威的 Lean 对应，供 Week 4-5 HIR/MIR 设计直接引用。

---

## §1 基础数值类型

### §1.1 使用频次（实测）

| C++ 类型 | 总出现 | 变量 | 参数 | 返回 | 典型用途 |
|---|---|---|---|---|---|
| `int` / `int32_t` | **3017** | 172 | 30 | 1 | 循环计数器、索引、degree 计数 |
| `bool` | 1052 | 20 | 1 | 6 | 控制流、flag |
| `uint64_t` | 471 | 23 | 21 | 0 | Zp 系数、UMonomial.deg、hash、prime |
| `int64_t` | 352 | 25 | 5 | 0 | 符号位敏感的度数计算 |
| `size_t` | 322 | 33 | 0 | 0 | 容器大小、索引 |
| `char` | 69 | 0 | 0 | 0 | 字符串字面量（assert 消息用，可忽略）|
| `double` | 30 | 2 | 0 | 0 | `__heuristic_starting_precision` 专用 |
| `uint32_t` | 25 | 0 | 0 | 0 | UMonomial 字段 |

### §1.2 映射表（精确模型）

| C++ | Lean 精确模型 | 语义模型 | 备注 |
|---|---|---|---|
| `uint64_t` | `UInt64` | `Fin (2^64)` / mod 2^64 | 核心类型，与 C++ unsigned 完全对应 |
| `int64_t` | `Int64` / `Int` | `Fin [-2^63, 2^63)` 或无界 `Int` | L1 用 `Int64`（精确），L2 可抽象为 `Int` |
| `uint32_t` | `UInt32` | `Fin (2^32)` | 仅 UMonomial 字段用，整型运算罕见 |
| `int32_t` / `int` | `Int32` / `Int` | 同上 | L1 `Int32`，L2 `Int` |
| `size_t` | `USize` 或 `Nat` | `Nat` | 容器索引；L1 可用 `USize`（平台相关），L2 直接 `Nat` |
| `ptrdiff_t` | `Int` | `Int` | 未实际使用（0 次）|
| `bool` | `Bool` | `Prop` 或 `Bool` | 原生 |
| `char` | `Char` 或 `UInt8` | `Fin 256` | 仅字符串字面量，翻译时可整体忽略 |
| `float` | `Float` | `Float` | **0 使用** |
| `double` | `Float` | `Float` | 仅 `__heuristic_starting_precision`（3 处 FloatingLiteral）|
| `unsigned __int128` | `UInt128`（自定义）| `Nat` | **CLPoly 因式分解不使用**（Barrett reduction 用；但因式分解 TRANSLATION_SCOPE 无）|

### §1.3 Lean 4 原生运算与 C++ 语义对比

参考 `lean-stdlib-catalog.md`。关键一致/不一致：

| 运算 | C++ uint64_t 语义 | Lean UInt64 行为 | 一致？ |
|---|---|---|---|
| `a + b` | mod 2^64 | mod 2^64 | ✅ |
| `a - b` | mod 2^64 | mod 2^64（不 panic） | ✅ |
| `a * b` | mod 2^64 | mod 2^64 | ✅ |
| `a / b` (b≠0) | `⌊a/b⌋` | `⌊a/b⌋` | ✅ |
| `a / 0` | **UB** | **返回 0（不 panic）** | ❌ 需 `require b ≠ 0` |
| `a % 0` | **UB** | **返回 a（不 panic）** | ❌ 需 `require b ≠ 0` |
| `a << n` (n<64) | `a × 2^n` mod 2^64 | 相同 | ✅ |
| `a << n` (n≥64) | **UB** | **wrap 到 n mod 64** | ❌ 需 `require n < 64` |
| `a >> n` (n<64) | `⌊a / 2^n⌋` | 相同 | ✅ |
| `a >> n` (n≥64) | **UB** | 相同 wrap | ❌ 需 `require n < 64` |
| `a == b` / `a != b` / `a < b` / ... | bool 比较 | 相同 | ✅ |
| `~a` | bitwise not | `UInt64.complement a` | ✅ |
| `a & b` / `a \| b` / `a ^ b` | bitwise | `.land` / `.lor` / `.xor` | ✅ |

**signed 整数**（`int`/`int64_t`）的 C++ 语义严格：**溢出是 UB**，但 Lean `Int32`/`Int64` 也以 mod wrap 实现。需要 `require` 明确排除溢出输入。

### §1.4 Bool 短路求值

| C++ | Lean | 一致？ |
|---|---|---|
| `a && b` 短路（b 不求值当 a=false） | `a && b` 短路（pattern match）| ✅ |
| `a \|\| b` 短路（b 不求值当 a=true） | 短路 | ✅ |
| `!a` | `!a` | ✅ |

**Lean 4 的 `Bool.and`/`Bool.or` 在编译器层面保证短路**（通过 pattern match）。这保证了含副作用运算（本项目无副作用，但原则上成立）的语义一致。

### §1.5 运算符的结果类型直方图（实测）

从 `numeric-types.md` 摘抄，用于验证 `operator_resolve` Pass 的类型映射完整：

| Operator | 结果类型 | 次数 | Lean 映射 |
|---|---|---|---|
| `BinaryOperator::<` | bool | 146 | `decide (a < b) : Bool`（L2 可 `a < b : Prop`）|
| `UnaryOperator::++` | int32 | 106 | `i := i + 1` |
| `UnaryOperator::!` | bool | 69 | `!a` |
| `BinaryOperator::-` | int32 | 50 | `a - b` |
| `BinaryOperator::==` | bool | 40 | `a == b` |
| `BinaryOperator::>` | bool | 37 | `a > b` |
| `BinaryOperator::&&` | bool | 34 | `a && b` |
| `BinaryOperator::+` | int32 | 31 | `a + b` |
| `BinaryOperator::=` | bool | 29 | SSA rebind |
| ... 详见 `numeric-types.md` | | | |

**全 40 种 (op, 结果类型) 对**都有 Lean 对应，无遗漏。

### §1.6 double / float 的处理

`__heuristic_starting_precision`（univar.hh:607）是**唯一**用 `double` 的函数。3 处 `FloatingLiteral`：`2.5`、`1.0`、`2.0`。计算 `std::log(p)`、`std::ceil(x)`。

**翻译策略**：
- L1 IR 用 `Float`（Lean 原生）+ `Float.log`、`Float.ceil` — 精确建模
- L2 抽象掉浮点（可用近似的 `Nat.log` 或作为黑盒公理）

此函数用于启发式选择 Hensel lifting 初始精度，**不影响正确性**（大精度总是安全的，只影响性能）。L1 翻译策略可以简化为"直接调用 `Float.log`/`Float.ceil` opaque"。

---

## §2 UB 分析与 `require` 生成

L1 IR 的一个核心工作：**把 C++ 的隐式 UB 点转为 Lean 的显式 `require` 前置条件**，否则 Lean 会静默执行（如 `UInt64 / 0 = 0`），和 C++ UB 不一致。

UB 分类的权威来源：`ub-catalog.md`（UB-1 到 UB-8）。实测站点：`ub-sites.md`（扫 65 函数 AST）。

### §2.1 UB 分类与实测站点（扫 65 函数）

| UB 类型 | C++ 触发 | Lean 行为 | 生成的 require | **实测次数** |
|---|---|---|---|---|
| **UB-1** 除以零 | `a / b`, `a % b`（整数）| Lean: `x/0→0`、`x%0→x` 不 panic | `require hb : b ≠ 0` | **31** |
| **UB-2** 数组越界 | `arr[i]`（`operator[]`）| Lean `.get!` 返回 default 不 panic | `require hi : i < arr.size` | **468** |
| **UB-3** 空容器 back/front | `v.back()`/`v.front()` on empty | CLPoly 用 `.front!` 返回 default | `require v.size > 0` | **42** |
| **UB-4** 移位越界 | `a << n`, `a >> n`（n≥位宽）| Lean wrap 到 `n mod bits` | `require hn : n < 64` | **0**（因式分解不用位移）|
| **UB-6** signed 溢出（+/-/*）| C++ `int + int` 溢出 | Lean `Int32/Int64` mod wrap | `require ¬overflow(a, b)` | **113** |
| **UB-7** unsigned→signed 截断 | `(int64_t)big_uint` | Lean `toInt64` mod wrap | `require x.val ≤ INT64_MAX` | **8** |
| **UB-8** assert 失败 | `assert(cond)` | Lean 不执行 | `require h : cond` | **23** |
| **合计** | | | | **685** |

### §2.2 按函数分布（前 10）

| 函数 | UB 总数 | 主导类型 |
|---|---|---|
| `__lll_reduce` | **132** | OOB 108（密集矩阵访问）、signed overflow 21 |
| `__wang_core` | 53 | OOB 35、signed 15 |
| `__wang_leading_coeff` | 53 | OOB 31、Div0 7、Empty 8 |
| `__si_vandermonde_solve` | 43 | OOB 38、signed 4 |
| `__mtshl_sparse_int` | 38 | OOB 38 |
| `__mtshl_wmds` | 33 | OOB 31、Empty 2 |
| `__zassenhaus_recombine` | 31 | OOB 14、signed 14 |
| `__mtshl_lift` | 30 | OOB 23、signed 7 |
| `__mtshl_step_j` | 25 | OOB 25 |
| `__vanhoeij_recombine` | 24 | OOB 9、signed 11 |

完整表见 `ub-sites.md`。**一个函数平均 ~10 个 UB 站点**，**平均一个 UB 点需一个 `require` 参数 + 一个精化证明义务**。

### §2.3 关键 Lean vs C++ 语义差异（必须生成 require）

从 `lean-stdlib-catalog.md` 摘抄。这些是 Lean **不主动** panic 但 C++ 是 UB 的点：

| 操作 | C++ | Lean (stdlib) | 翻译必须 |
|---|---|---|---|
| `UInt64 / 0` | UB | 返回 0 | `require hb : b ≠ 0` |
| `UInt64 % 0` | UB | 返回 a | `require hb : b ≠ 0` |
| `a <<< n` (n≥64) | UB | wrap to `n mod 64` | `require hn : n < 64` |
| `Array.get! i`（i ≥ size）| UB（vector 无检查）| 返回 `default` | `require hi : i < arr.size` |
| `Array.set! i v`（i ≥ size）| UB | silent no-op（return unchanged）| `require hi : i < arr.size` |
| `Array.pop empty` | 未定义 | 返回空 array | `require ¬arr.isEmpty` |
| signed 加减乘溢出 | UB | mod wrap | `require Int.noOverflow(...)` |

**所以"Lean 不 panic = 翻译必须显式 require"**。这是为什么 685 个 UB 站点都要生成义务。

### §2.4 有符号溢出（UB-6）的处理策略

113 个 signed 溢出站点是最大的挑战。CLPoly 的 signed 用途：

| 用途 | 占比 | 策略 |
|---|---|---|
| **循环计数器** `for (int i = 0; i < n; ++i)` | ~60% | 通过**循环不变量** `0 ≤ i < n < 2^31` 自动消化 |
| **度数计算** `deg(f) - deg(g)` | ~20% | 通过 `deg` 的自然范围（≤ 2^32）保证 |
| **子集枚举索引** | ~10% | 小数字（`< 2^20`），trivial |
| **符号位计算** | ~10% | `if (x < 0) ...`，结果是 bool 不是算术 |

**实际需要人工证明的溢出 require 预估：< 10 个**（大部分靠循环不变量自动化）。

### §2.5 assert(...) 的模式识别

从 AST 识别 `assert(cond)`：

C++ 预处理展开后（Debug 编译）：
```c
if (!cond) __assert_fail("cond", "file.cc", line, __PRETTY_FUNCTION__);
```

AST 模式：
```
IfStmt
├── UnaryOperator(!) (条件)
└── CallExpr(__assert_fail, ...)
```

翻译器识别此模式 → 消除 IfStmt，提取 `cond` 作为 `require h : cond` 放入函数签名。23 处 assert 站点。

### §2.6 throw / Except 确认不使用

`CXXThrowExpr` 和 `CXXTryStmt` 在 AST 中均为 **0**。因式分解代码的所有错误路径通过 return 值表达（如 `empty result`、`false flag`）。

**L1 IR 不引入 `Except`**，所有函数返回普通类型 / tuple。这简化了 HIR 节点定义和 `ref_elim` Pass。

### §2.7 require 总量预测

| UB 类型 | 站点数 | 平均复用率 | 预估 require 参数数 |
|---|---|---|---|
| UB-1 除零 | 31 | ~80%（同一变量多次除法）| ~10 独立 require |
| UB-2 OOB | 468 | ~70%（循环内索引复用）| ~150 独立 require |
| UB-3 Empty | 42 | ~60% | ~15 独立 |
| UB-6 Signed | 113 | ~90%（循环不变量大覆盖）| ~10 独立 |
| UB-7 U→S | 8 | 0% | 8 |
| UB-8 Assert | 23 | 0% | 23 |
| **合计** | 685 | | **~200 独立 require 参数** |

平均每函数 **~3 个显式 require 参数**。精化证明时通过 Lean 的自动决策过程（`omega`、`decide`）消化大部分。

### §2.2 signed 溢出的特殊处理

C++ `int` / `int64_t` 的溢出是 **UB**。Lean `Int32` / `Int64` 以 mod wrap 实现，不 panic。

CLPoly 中 signed 整数的主要用途：
1. **循环计数器**（3017 次 `int32_t`）：`for (int i = 0; i < n; ++i)` — 只要 `n ≤ INT_MAX` 就不溢出。通常通过函数前置条件（`n < 2^31`）保证
2. **度数计算**：可能计算 `deg1 - deg2`，需要 `require deg1 ≥ deg2` 或用 unsigned
3. **符号信息**：`sign`、`comparison result`，不做算术

**策略**：翻译器**保守地**对每个 signed 算术生成 `require -2^31 ≤ result < 2^31`（或 `-2^63` 等）。这使义务数激增但语义精确。后续可以在**精化证明**阶段用上下文信息（循环不变量、参数界）来消化义务。

### §2.3 assert(...) 的处理

`numeric-types.md` 没统计 assert，但从 `cpp-construct-catalog.md` 可知 CLPoly 用 `assert()` 做**开发时不变量检查**（Release 编译 no-op）。

AST 层面，`assert(cond)` 展开为：
```
if (!cond) __assert_fail(...)
```
`__assert_fail` 是 `noreturn` 标注，Clang 知道。

**翻译策略**：
- 识别模式 `if (!cond) __assert_fail(...)` → 生成 `require h : cond`
- 不翻译 `__assert_fail` 的具体参数（字符串消息忽略）

### §2.4 throw / Except 的处理

CLPoly 因式分解代码**完全不使用 throw**（`CXXThrowExpr` 和 `CXXTryStmt` 都是 0）—— 这是原 `blueprint.md` 曾预留的错误处理，但实际代码里错误路径用 `return` 返回空值或特殊值表示。

**翻译策略**：**不引入 `Except`**，所有 L1 函数返回普通类型。这比 `blueprint.md` §3.2 的预设更简洁。

### §2.5 require 参数的分布预测

基于实测 AST 统计：

| require 类型 | 预估生成次数 | 说明 |
|---|---|---|
| 数组越界 | 200-400 | 每次 `operator[]` 生成 1 个 |
| 除以零 | 7-10 | 每次 `/` / `%` 一次 |
| 移位越界 | 4-5 | 较少 |
| signed 溢出 | **若保守** 数百 | 每次 signed 算术；可通过循环不变量化简 |
| assert | ~50 | `cpp-construct-catalog.md` 未精确统计，估计这个量级 |
| `.front()`/`.back()` empty | ~30 | `MvPolyZZ.front!` 等方法 |

**总计**：一个函数平均 **5-15 个 require**，全项目 **~500-1500 require 参数**。

这是 L1→L2 精化证明的主要义务来源。大部分由循环不变量、范围推理自动消化；少数需要人工提示。

### §2.6 require 的 Lean 语法

```lean
partial def nmod_mul_ir (a b p : UInt64) (norm : UInt32)
    (hp : p ≠ 0)                    -- assert(p != 0) → require
    (h_shift : norm.val < 64)       -- 隐式 UB: shift → require
    : UInt64 :=
  let pn := p <<< norm.val          -- norm.val 已有 h_shift 前提
  let q := a / p                    -- p 已有 hp 前提
  ...
```

- `require` 参数紧跟普通参数，在函数签名里
- 调用方必须在调用点提供证明
- 函数体内可直接使用被保证的值

**`ref_elim` Pass 之后**，require 和输出参数的 tuple 一起组织：
```
signature: (value_params) (require_params) → (ret_tuple)
```

---

## §3 CLPoly 核心类型（Day 2 起草）

*（本节留到 Week 2 Day 2 完成。计划覆盖 `ZZ`、`QQ`、`Zp`、`Variable`、`basic_monomial`、`polynomial_`、`upolynomial_`、`factorization`）*

---

## §4 STL 容器 shim（Day 3）

*（计划覆盖 `vector` → `Array`、`pair` → `×`、`map` → `StdMap`、`sort` → `qsortWith`、RNG 三件套公理化）*

---

## §5 模板实例化策略（Day 4）

*（计划：`factorize` 3 实例命名与 HIR 处理）*

---

## §6 Factorization / WangLcResult / PrimeSelectionResult struct（Day 4）

*（计划：CLPoly 结果类型的 Lean 结构化定义）*

---

## §7 引用与指针（Day 5）

*（计划：`T&` → tuple 返回；`T*` 的 1 处特殊处理）*

---

## §8 CAST_TABLE 完整枚举（Day 5）

*（计划：所有 263 unique cast 三元组 → Lean 表达式映射）*

---

## 当前进度

- [x] **Day 1 (2026-04-21)**：§1 基础数值 + §2 UB 分析 完成
- [ ] **Day 2**：§3 CLPoly 核心类型
- [ ] **Day 3**：§4 STL 容器 shim
- [ ] **Day 4**：§5 + §6 模板实例化 + struct 定义
- [ ] **Day 5**：§7 + §8 引用指针 + CAST_TABLE
