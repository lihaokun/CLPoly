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

## §3 CLPoly 核心类型

### §3.1 使用频次（实测，from `types.md`）

| C++ 类型（简化）| 带 const 出现 | 无 const 出现 | 合计 | 宿主数 |
|---|---|---|---|---|
| `upolynomial_<Zp>` | 112 | 340 | **452** | Zp 模块主导 |
| `upolynomial_<ZZ>` | 103 | 325 | **428** | Univar 模块主导 |
| `Zp` | 81 | 243 | **324** | 22 个函数 |
| `basic_polynomial<...Zp...>` | 79 | 207 | **286** | polynomial_<Zp,...> 底层 |
| `PolyZp`（typedef 别名）| — | 255 | **255** | Wang 模块内部 |
| `polynomial_<ZZ, lex_<less>>` | 79 | 115 | **194** | 多变量 ZZ |
| `basic_monomial<lex_<less>>` | — | 79 | **79** | 16 函数 |
| `polynomial_<Zp, lex_<less>>` | 66 | 80 | **146** | 多变量 Zp 中间 |
| `basic_polynomial<umonomial,Zp,uless>` | — | 79 | **79** | upolynomial_<Zp> 底层 |
| `clpoly::ZZ` | — | 107 | **107** | 25 函数 |
| `std::vector<upolynomial_<ZZ>>` | — | 99 | **99** | 因子列表 |
| `std::vector<PolyZp>` | — | 81 | **81** | — |
| `std::vector<Zp>` | — | 68 | **68** | 系数序列 |

### §3.2 类型映射决策表（总览）

来自 `mathlib-poly-types.md` 调研结论：

| CLPoly C++ 类型 | L1 映射 | L2 精化 | 决策理由 |
|---|---|---|---|
| `ZZ`（GMP 封装）| `ZZ : structure { val : Int }` 或直接 `Int` | `Int` | L1 可完全用 `Int`（GMP 精确到 Int 一致）|
| `QQ`（有理数）| 直接 `Rat` | `Rat` | Mathlib 已完整 |
| `Zp`（运行时素数）| **自定义** `Zp { val : UInt64, prime : UInt64 }` | `ZMod p.toNat` | `ZMod p` 要求 p 编译时，CLPoly 是运行时 |
| `variable` | `Variable = Nat` 或 `Fin n`（简化）| `Fin n`（严格）| L1 直接用 UInt64 作索引 |
| `basic_monomial<order>` | **自定义** `Monomial { exps : Array UInt64, order : MonomialOrder }` | `Finsupp (Fin n) ℕ` | order 信息必须保留于 L1 |
| `umonomial`（单变量单项式）| **自定义** `UMonomial { deg : UInt64 }` | `ℕ` | 只需一个度数 |
| `upolynomial_<T>` | **自定义** `SparsePoly T = Array (UMonomial × T)` | `Polynomial T`（经 `toPoly`）| 稀疏 + 降序，和 C++ 一致 |
| `polynomial_<T, order>` | **自定义** `MvPoly T = Array (Monomial × T)` | `MvPolynomial σ T` | 稀疏 + 单调序 |
| `factorization<Poly>` | **自定义** `Factorization Poly = { content : T, factors : Array (Poly × UInt64) }` | `Multiset (Poly × ℕ)`（L2 已用）| 结构直接对应 |

### §3.3 每个类型的详细 Lean 定义

#### §3.3.1 `ZZ`

**C++**：GMP `mpz_t` 封装，`clpoly::ZZ` struct 带 `int64_t _val` + `mpz_ptr _mpz`（小整数路径 + mpz 任意精度路径）。

**L1 Lean**：

```lean
/-- L1 精确模型：保留双路径可能但简化为 Int -/
abbrev ZZ := Int

-- 关键操作（映射到 Int stdlib）
-- ZZ.abs : ZZ → ZZ := Int.natAbs |> Int.ofNat
-- ZZ.add, ZZ.sub, ZZ.mul, ZZ.div, ZZ.mod : 已在 Int
-- ZZ.gcd : ZZ → ZZ → ZZ := Int.gcd (Nat)
```

**AST 实测调用**（from `clpoly-types-usage.md`）：
- `.operator bool()`（12 次）—— 判零：`ZZ != 0` → `x ≠ 0`
- `.sizeinbase(base)`（3 次）—— 位数估计：`Int.log2 x + 1` 或 Mathlib `Nat.log`
- `.fdiv_ui(d)`（3 次）—— 向下除以 unsigned：`Int.fdiv x d`

**L2 精化**：`Int`。映射函数：`toInt : ZZ → Int := id`。

#### §3.3.2 `QQ`

**C++**：两个 `ZZ`（分子 + 分母），自动化简。

**L1 Lean**：

```lean
abbrev QQ := Rat

-- QQ.num, QQ.den : 已在 Rat
-- 自动化简：Rat.mk 默认 reduce
```

**实测使用**：极少（0 次 as variable / parameter type）。仅 `__lll_reduce` 和 `__heuristic_starting_precision` 间接用（通过 `std::vector<QQ>`）。

**L2 精化**：`Rat` 直接。

#### §3.3.3 `Zp`

**C++**：
```cpp
struct Zp {
  int64_t _val;
  uint64_t _p;
  // ... Barrett reduction 辅助字段 ...
};
```

**L1 Lean**（已在 `clpoly_model.lean` 部分定义）：

```lean
structure Zp where
  val : UInt64
  prime : UInt64
deriving Repr, Inhabited

def Zp.mk (v : Int) (p : UInt64) : Zp :=
  let n := (v % p.toNat).toNat
  ⟨n.toUInt64, p⟩

/-- 加/减/乘：Barrett 简化为 `val_x op val_y mod prime` -/
def Zp.add (a b : Zp) : Zp := ⟨(a.val + b.val) % a.prime, a.prime⟩
def Zp.sub (a b : Zp) : Zp := ⟨(a.val + a.prime - b.val) % a.prime, a.prime⟩
def Zp.mul (a b : Zp) : Zp := ⟨(a.val * b.val) % a.prime, a.prime⟩
def Zp.neg (a : Zp) : Zp := ⟨(a.prime - a.val) % a.prime, a.prime⟩
def Zp.inv (a : Zp) : Zp := sorry  -- 扩展欧几里得

/-- 相等仅比较值；prime 必须相同（若不同是 bug）-/
instance : BEq Zp where
  beq a b := a.val == b.val && a.prime == b.prime
```

**AST 实测方法**（from `clpoly-types-usage.md`）：
- `.prime()`（9 次）—— 返回 `prime` 字段
- `.inv()`（2 次）—— 模乘逆元
- `.number()`（1 次）—— 返回值字段

**L2 精化**：
```lean
def Zp.toZMod : (z : Zp) → ZMod z.prime.toNat :=
  fun z => ZMod.mk z.prime.toNat z.val.toNat
```
**需前置条件** `Fact (z.prime.toNat.Prime)`，由 `__select_prime` 建立。

#### §3.3.4 `Variable`

**C++**：`clpoly::variable` — 简单封装（可能带 name string）。

**L1 Lean**：

```lean
abbrev Variable := UInt64  -- 变量索引
```

**实测**：13 次出现，仅作索引。不访问字符串名字段。

#### §3.3.5 `basic_monomial<order>` 与 `umonomial`

**C++**：
```cpp
template<class order_> struct basic_monomial {
  std::vector<std::pair<variable, uint32_t>> _data;
  // order_ 是 tag 类型，决定 operator< 的实现
};

struct umonomial {
  uint32_t deg;
};
```

**L1 Lean**：

```lean
/-- 多变量单项式：变量→指数 的数组，保持 order 不变量 -/
structure Monomial where
  /-- 按变量 id 排序的 (var, exp) 对；exp > 0 -/
  data : Array (Variable × UInt32)
deriving Repr, Inhabited

/-- Monomial order 作为参数传递（lex / grlex / etc.）-/
inductive MonomialOrder where
  | lex
  | grlex
  | degrevlex

/-- 单变量单项式 -/
structure UMonomial where
  deg : UInt64
deriving Repr, Inhabited
```

**AST 方法**：
- `basic_monomial.begin()`/`.end()`（各 4 次）—— range-for 触发的迭代器
- `umonomial.deg`（1 次显式访问）—— 度数字段

**L2 精化**：
- `Monomial.toFinsupp : Monomial → (Fin n →₀ ℕ)`
- `UMonomial.toNat : UMonomial → ℕ := .deg.toNat`

#### §3.3.6 `upolynomial_<T>`

**C++**：
```cpp
template<class T> using upolynomial_ = basic_polynomial<umonomial, T, uless>;
// basic_polynomial 是 std::vector<std::pair<Mon, T>>
```

**L1 Lean**：

```lean
/-- 单变量稀疏多项式：降序 (UMonomial, coeff) 对列表 -/
abbrev SparsePoly (T : Type) := Array (UMonomial × T)

abbrev SparsePolyZZ := SparsePoly Int
abbrev SparsePolyZp := SparsePoly Zp

-- 不变量：degrees 严格降序，coeffs 非零
-- 这些不变量在 L1 不强制（partial def 不验证），
-- 在 L2 精化中作为 `toPoly` 的前置条件
```

**L2 精化**：

```lean
def SparsePoly.toPolyZp (p : SparsePolyZp) : Polynomial (ZMod p_val) :=
  p.foldl (fun acc (u, c) => acc + Polynomial.C c.toZMod * Polynomial.X ^ u.deg.toNat) 0
```

#### §3.3.7 `polynomial_<T, order>`

**C++**：
```cpp
template<class T, class order> using polynomial_ =
  basic_polynomial<basic_monomial<order>, T, order>;
```

**L1 Lean**：

```lean
/-- 多变量稀疏多项式：按 order 排序的 (Monomial, coeff) 对列表 -/
structure MvPoly (T : Type) where
  terms : Array (Monomial × T)
  order : MonomialOrder
deriving Repr, Inhabited

abbrev MvPolyZZ := MvPoly Int
abbrev MvPolyZp := MvPoly Zp
```

**L2 精化**：`MvPolynomial (Fin n) T`，通过 `toMvPolynomial` 折叠。

#### §3.3.8 `factorization<Poly>`

**C++**：
```cpp
template<class Poly> struct factorization {
  typename Poly::coeff_type content;
  std::vector<std::pair<Poly, uint64_t>> factors;
};
```

**L1 Lean**：

```lean
structure Factorization (Poly : Type) (Coeff : Type) where
  content : Coeff
  factors : Array (Poly × UInt64)
deriving Repr, Inhabited
```

**实例化**：
- `Factorization SparsePolyZZ Int`（`factorize(upolynomial_<ZZ>)` 返回）
- `Factorization MvPolyZZ Int`（多变量）
- `Factorization SparsePolyZp Zp`（Zp 因式分解）

**L2 精化**：`Multiset (Poly × ℕ)`（数学乘积）。

#### §3.3.9 辅助 struct（`__prime_selection_result`、`__wang_lc_result`、`__hensel_node`）

从 C++ 源：

```cpp
struct __prime_selection_result {
  uint64_t prime;
  std::vector<upolynomial_<Zp>> factors;
  bool irreducible;
};

struct __wang_lc_result {
  bool success;
  std::vector<polynomial_<ZZ, lex>> scaled_factors;
  std::vector<polynomial_<ZZ, lex>> lc_targets;
  ZZ delta;
};

struct __hensel_node {
  /* 树节点：g, h, s, t, parent_idx, children ... */
};
```

**L1 Lean**：

```lean
structure PrimeSelectionResult where
  prime : UInt64
  factors : Array SparsePolyZp
  irreducible : Bool

structure WangLcResult where
  success : Bool
  scaled_factors : Array MvPolyZZ
  lc_targets : Array MvPolyZZ
  delta : Int

structure HenselNode where
  g : SparsePolyZZ
  h : SparsePolyZZ
  s : SparsePolyZZ
  t : SparsePolyZZ
  parent_idx : Int
  left_idx : Int
  right_idx : Int
```

**L2 精化**：这些是算法内部结构，L2 不直接建模；精化证明时"解包"成每个字段。

### §3.4 Gap 分析：`clpoly_model.lean` 现状 vs 需要

基于 `clpoly-model-inventory.md`（380 行，14 types、1 axiom、~40 functions、6 TODOs）：

| 需要 | 已有 | 状态 |
|---|---|---|
| `ZZ`（别名 `Int`）| ✓ 已定义 `abbrev ZZ := Int` | OK |
| `QQ`（别名 `Rat`）| ✓ 已定义 | OK |
| `Zp` 结构体 | ✓ 已定义 `Zp { val, prime }` + `Zp.add/sub/mul/neg` | OK |
| `Zp.inv` | ❌ **TODO**（扩展欧几里得未写）| 补齐 |
| `Variable` | ✓ `abbrev Variable := UInt64` | OK |
| `Monomial` | ❌ 简化为"指数向量数组"，缺 `MonomialOrder` 标记 | **重做** |
| `UMonomial` | ✓ 已定义 | OK |
| `SparsePoly` / `SparsePolyZZ` / `SparsePolyZp` | ✓ 已定义 | OK |
| `MvPoly` / `MvPolyZZ` / `MvPolyZp` | ✓ 基础结构，缺 order 字段 | **重做** |
| `Factorization` | ✓ 已定义 | OK |
| `PrimeSelectionResult` / `WangLcResult` / `HenselNode` | ❌ 未定义 | **新增** |
| `StdMap` | ✓ 基础，但迭代器 `.end/.begin` 是占位 | 翻译时改用 filter/find，不依赖裸迭代器 |
| `Rng`（公理化）| ✓ 已定义 | OK |
| `__upoly_divmod_mod` Lean 实现 | ❌ TODO（占位 `(f, #[])`）| 补齐 |
| `__upoly_gcd` Lean 实现 | ❌ TODO | 补齐或用 `EuclideanDomain.gcd` |

### §3.5 需要 `clpoly_model.lean` 新增的内容

```lean
-- §3.5.1 MonomialOrder 枚举
inductive MonomialOrder where
  | lex
  | grlex
deriving Repr, Inhabited, BEq

-- §3.5.2 MvPoly 改版（加 order 字段）
structure MvPoly (T : Type) where
  terms : Array (Monomial × T)
  order : MonomialOrder
deriving Repr, Inhabited

-- §3.5.3 辅助 struct
structure PrimeSelectionResult where
  prime : UInt64
  factors : Array SparsePolyZp
  irreducible : Bool
deriving Repr, Inhabited

structure WangLcResult where
  success : Bool
  scaled_factors : Array MvPolyZZ
  lc_targets : Array MvPolyZZ
  delta : Int
deriving Repr, Inhabited

structure HenselNode where
  g : SparsePolyZZ
  h : SparsePolyZZ
  s : SparsePolyZZ
  t : SparsePolyZZ
  parent_idx : Int
  left_idx : Int
  right_idx : Int
deriving Repr, Inhabited

-- §3.5.4 Zp.inv（扩展欧几里得）
partial def Zp.inv (a : Zp) : Zp :=
  -- 扩展欧几里得：a * inv_a ≡ 1 mod prime
  sorry  -- Week 2 Day 3 补齐，或用 Mathlib 的 ZMod.inv 的公理对应
```

### §3.6 命名约定

从 `clpoly_model.lean` 现状 + `cpp-construct-catalog.md` 的实例化发现：

| C++ 类型族 | L1 命名约定 | 示例 |
|---|---|---|
| 有具体类型参数 | `<Family><TypeSuffix>` | `SparsePolyZp`, `MvPolyZZ` |
| 泛型多项式族 | `<Family>` | `SparsePoly`, `MvPoly` |
| 辅助 struct | `<C++Name 驼峰>` | `PrimeSelectionResult`, `WangLcResult` |
| 枚举 | `<Family><Kind>` | `MonomialOrder` |

### §3.7 与 `factorize` 3 实例的对接

`factorize` 3 实例对应 3 种输入类型：

| 实例 | 输入类型 | Lean 命名 | 返回类型 |
|---|---|---|---|
| `factorize<upolynomial_<ZZ>>` | 单变量 ZZ | `factorize_upoly_ir` | `Factorization SparsePolyZZ Int` |
| `factorize<polynomial_<ZZ, lex>>` | 多变量 lex | `factorize_lex_ir` | `Factorization MvPolyZZ Int` |
| `factorize<polynomial_<ZZ, grlex>>` | 多变量 grlex | `factorize_grlex_ir` | `Factorization MvPolyZZ Int` |

多变量的 `lex` vs `grlex` 仅在 `MvPoly.order` 字段区分，其他逻辑共享。

---

## §4 STL 容器 shim

### §4.1 `std::vector<T>` → `Array T`

**使用频次**：2066 次（`stl.md` 数据），最常见容器。

**操作映射表**：

| C++ 操作 | Lean 对应 | 注意事项 |
|---|---|---|
| `v.push_back(x)` | `v.push x` | 返回新 Array（SSA 化） |
| `v.pop_back()` | `v.pop` | 返回新 Array，**空 Array 时 no-op** |
| `v.back()` | `v.back!` | 空 Array 时 `default` — 需 `require v.size > 0` |
| `v.front()` | `v.front!` | 同上 |
| `v[i]` 读 | `v[i]!` / `v.get! i` | 越界 `default` — 需 `require i < v.size` |
| `v[i] = x` | `v.set! i x` | 越界 silent no-op — 需 `require i < v.size` |
| `v.size()` | `v.size` | 返回 `Nat` |
| `v.empty()` | `v.isEmpty` | `v.size == 0` |
| `v.reserve(n)` | `id`（noop） | Lean Array 无容量预留 |
| `v.clear()` | `Array.empty` | 新建空 Array |
| `v.begin()`, `v.end()` | 不直接翻译 | 迭代器模式由 `iter_recognize` 处理 |
| `v.erase(it, v.end())` | `v.take k`（k 由 it 推导） | compact 模式专用 |
| `v.erase(v.begin() + i)` | `v.eraseIdx! i`（自定义） | 删除指定位置 |
| `v.insert(v.begin() + i, x)` | `v.insertAt! i x`（自定义） | 插入 |
| `std::vector<T>(n)` | `Array.mkArray n default` | 大小 n 默认值 |
| `std::vector<T>(n, v)` | `Array.mkArray n v` | 大小 n 值 v |

**注**：Lean `Array.pop` 对空 Array 返回空 Array（safe），但 C++ `pop_back` 是 UB — 翻译器应在 `require ¬v.isEmpty` 保证。

### §4.2 `std::pair<A, B>` → `A × B`

**使用频次**：685 次。

**操作映射**：

| C++ 操作 | Lean 对应 |
|---|---|
| `std::pair<A, B>(a, b)` | `(a, b)` 或 `Prod.mk a b` |
| `p.first` | `p.1` 或 `p.fst` |
| `p.second` | `p.2` 或 `p.snd` |
| `std::make_pair(a, b)` | `(a, b)` |
| 结构化绑定 `auto [x, y] = p` | `let (x, y) := p` |

**Lean 4 语法**：Prod 直接使用元组语法 `(a, b) : α × β`。

### §4.3 `std::map<K, V>` → 自定义 `StdMap K V`

**使用频次**：94 次（主要在 `__extract_monomial_content`、`__mtshl_sparse_int`、`__select_eval_point` 等 7 处）。

**C++ 语义**：有序关联容器，按 key 排序，支持 `O(log n)` find/insert/erase。

**Lean 端设计：按 key 排序的 Array**

```lean
/-- StdMap：按 key 严格升序的 Array；仿 std::map 语义 -/
structure StdMap (K V : Type) where
  entries : Array (K × V)
  -- 不变量（不强制）：entries[i].1 < entries[i+1].1
deriving Repr, Inhabited

namespace StdMap

/-- 空 map -/
def empty {K V : Type} : StdMap K V := ⟨#[]⟩

/-- 查找 key，返回 Option V（CLPoly 用 `.find` 返回迭代器，shim 用 Option） -/
def find [BEq K] [Ord K] (m : StdMap K V) (k : K) : Option V :=
  -- 二分查找或线性（CLPoly 的 map size 很小，线性足够）
  m.entries.findSome? (fun (k', v) => if k' == k then some v else none)

/-- 插入或更新 -/
def insert [BEq K] [Ord K] (m : StdMap K V) (k : K) (v : V) : StdMap K V :=
  -- 线性扫描插入点（保持有序）
  go m.entries 0
where
  go (arr : Array (K × V)) (i : Nat) : StdMap K V :=
    if h : i < arr.size then
      let (k', _) := arr[i]
      if k' == k then ⟨arr.set! i (k, v)⟩
      else if compare k k' = .lt then ⟨arr.insertAt! i (k, v)⟩  -- 需自定义 insertAt!
      else go arr (i + 1)
    else ⟨arr.push (k, v)⟩

/-- 删除 key -/
def erase [BEq K] (m : StdMap K V) (k : K) : StdMap K V :=
  ⟨m.entries.filter (fun (k', _) => !(k' == k))⟩

/-- 是否包含 key -/
def contains [BEq K] (m : StdMap K V) (k : K) : Bool :=
  m.entries.any (fun (k', _) => k' == k)

/-- 大小 -/
def size (m : StdMap K V) : Nat := m.entries.size

/-- 转 List 方便 foldl/遍历 -/
def toList (m : StdMap K V) : List (K × V) := m.entries.toList

/-- 从 Array 构造（假设已排序） -/
def fromSortedArray (arr : Array (K × V)) : StdMap K V := ⟨arr⟩

end StdMap
```

**为什么用 Array 而非红黑树**：
- CLPoly 的 map 都是小尺寸（度数 0..n 作 key，n < 1000）
- Array + 二分查找 O(log n)，性能足够
- 实现简单，`#eval` 可执行
- 精化证明时只关心"key→value 映射"语义，实现细节无所谓

**已知的 C++ `std::map` 操作及其 Lean 对应**：

| C++ | Lean | 备注 |
|---|---|---|
| `m[k]` 读 | `(m.find k).getD default` | 键不存在返回 default |
| `m[k] = v` | `m.insert k v` | 插入或更新 |
| `m.find(k)` → 迭代器 | `m.find k : Option V` | 简化为 Option |
| `m.erase(k)` | `m.erase k` | |
| `m.erase(it)` | `m.erase key_of_it` | shim：先从 it 拿 key |
| `m.empty()` | `m.size == 0` | |
| `m.size()` | `m.size` | |
| `for (auto& [k, v] : m)` | `m.toList.forM (fun (k, v) => ...)` | 结构化绑定对应 §4.2 |
| `m.begin()` / `m.end()` | 由 range-for 吸收 | 裸迭代器不翻译 |

### §4.4 `std::sort` → `Array.qsort`

**使用频次**：8 处。详见 `lean-sort-api.md`（Agent 调研）。

**直接对应**：Lean 4 stdlib `Array.qsort` 签名

```
Array.qsort : {α : Type} → Array α → (α → α → Bool) → Array α
```

- 比较器风格：**Bool less-than**（返回 `true` 表示 `a < b`），与 C++ `std::sort` 完全一致
- Computable ✅（可 `#eval`）
- 非稳定排序（quicksort）；CLPoly 对稳定性无要求
- 最坏 O(n²)，平均 O(n log n)

CLPoly 的 5 处 std::sort 的比较器都已改为显式类型（Day 1 修复）：

```cpp
std::sort(v.begin(), v.end(),
    [](const Pair& a, const Pair& b) { return a.first.deg() < b.first.deg(); });
```

**Lean 翻译**：

```lean
v.qsort (fun a b => a.fst.deg < b.fst.deg)
```

一行对应，无需自实现。Lean 4 stdlib 无 `qsortWith`/`sortBy` 变体，但用 Bool 比较器可表达所有情况。

### §4.5 RNG 三件套

**使用频次**：`std::mt19937` 27、`std::uniform_int_distribution` 6、`std::random_device` 3。
**宿主函数**：`__upoly_random`、`__edf_Zp`、`__factor_Zp`、`__mtshl_sparse_int`。

**Lean 端公理化**（已在 `clpoly_model.lean` 部分定义）：

```lean
/-- Rng: 伪随机数生成器状态，公理化（类似 State monad）-/
axiom Rng : Type
axiom Rng.mk : UInt64 → Rng
axiom Rng.next : Rng → UInt64 → (UInt64 × Rng)
-- 语义：Rng.next rng max 返回 ((rand in [0, max]), new_rng)

/-- std::random_device shim -/
axiom Rng.seed : UInt64 := default  -- 或由 __mtshl_sparse_int 直接用
```

**C++ → Lean 映射**：

| C++ | Lean | 备注 |
|---|---|---|
| `std::mt19937 rng(seed)` | `let rng := Rng.mk seed` | |
| `std::mt19937 rng(rd())` | `let rng := Rng.mk Rng.seed` | `rd()` 公理化 seed |
| `std::uniform_int_distribution<uint64_t> dist(0, n-1)` | 融入 `Rng.next` | dist 不是独立对象 |
| `dist(rng)` | `let (r, rng') := Rng.next rng n` | 返回 `(sample, new_rng)` |

**约束**：`Rng` 是 opaque；精化证明不展开，只证"分裂函数返回合法因子"的概率性质。这与 L2 的 `exists_nonQR_poly`（EDF 中 AdjoinRoot 的有限域计数）对接。

### §4.6 其他 STL 工具

| C++ | 次数 | Lean 对应 | 备注 |
|---|---|---|---|
| `std::move(x)` | **63** | `x`（id） | Lean 值语义，move 无实际效果 |
| `std::swap(a, b)` | 4 | `let (a', b') := (b, a); ...` | Lean 无 in-place swap，用 let-rebind；SSA 化 |
| `std::iota(begin, end, v)` | 5 | `Array.range n \| >.map (fun i => i + v)` | 生成 `[v, v+1, ..., v+n-1]` |
| `std::max(a, b)` | 7 | `max a b` | 原生 |
| `std::min(a, b)` | 3 | `min a b` | 原生 |

**用例示例 — `std::iota`**：

```cpp
std::vector<int> T(r);
std::iota(T.begin(), T.end(), 0);  // T = [0, 1, ..., r-1]
```

```lean
let T := Array.range r  -- 返回 #[0, 1, ..., r-1]
-- 若 std::iota 起始值非 0：
-- let T := (Array.range r).map (· + start)
```

### §4.7 `initializer_list`

**使用频次**：3 处（`__ddf_Zp`、`__mtshl_zp_univar_mdp`、`__upoly_powmod`）。

**C++ 用法**：`std::vector<int>{1, 2, 3}` 或 `func({a, b, c})`。

**Lean 对应**：`#[a, b, c]`（Array 字面量）或 `[a, b, c]`（List 字面量）。

### §4.8 不使用的 STL

根据 `stl.md`，CLPoly 因式分解**不使用**以下 STL（全 0 次）：

- `std::unique_ptr` / `std::shared_ptr` / `std::weak_ptr`（智能指针）
- `std::thread` / `std::mutex` / `std::atomic`（并发）
- `std::exception` / `std::runtime_error`（异常类）
- `std::cout` / `std::cerr`（IO）
- `std::string` / `std::string_view`（仅字符串字面量用于 assert 消息，翻译时忽略）

简化 HIR：不用定义对应的 Lean shim。

### §4.9 总结：需要在 `clpoly_model.lean` 新增的容器 shim

```lean
-- §4.9.1 StdMap（见 §4.3 完整实现）
structure StdMap (K V : Type) where ... deriving ...
namespace StdMap
  def empty, find, insert, erase, contains, size, toList, fromSortedArray
end StdMap

-- §4.9.2 Array 辅助（CLPoly 专用）
def Array.eraseIdx! {α : Type} [Inhabited α] (a : Array α) (i : Nat) : Array α :=
  -- 删除位置 i 的元素（无 copy 原生性能差，但语义清晰）
  (a.toList.eraseIdx i).toArray  -- 占位实现

def Array.insertAt! {α : Type} (a : Array α) (i : Nat) (x : α) : Array α :=
  (a.toList.insertIdx i x).toArray  -- 占位实现

-- §4.9.3 Rng（已公理化，保持）
axiom Rng : Type
axiom Rng.mk : UInt64 → Rng
axiom Rng.next : Rng → UInt64 → (UInt64 × Rng)
axiom Rng.seed : UInt64
```

需要新增的代码量约 **60-80 行 Lean**，主要是 `StdMap` 的 8 个方法 + Array 辅助 + RNG 的 seed 公理。

---

## §5 模板实例化策略

### §5.1 规模

65 TRANSLATION_SCOPE 函数中：
- **64 个** 恰好 1 个 mangled 实例化 → 1 Lean 定义
- **1 个**（`factorize`）有 3 个 mangled 实例化 → 3 Lean 定义

**总计 67 个 Lean `_ir` 定义。**

### §5.2 `factorize` 3 实例的本质：**3 个独立源函数**

关键认识：**"factorize 有 3 个实例"不是"一个模板被 3 种类型具体化"**，而是**同一函数名 `factorize` 有 3 个互相独立的 C++ 源定义**（C++ 函数重载）：

| 实例 | 源位置 | 源形式 | Lean 命名 |
|---|---|---|---|
| `factorize<upolynomial_<ZZ>>` | `polynomial_factorize_univar.hh:1564` | **非模板**（重载） | `factorize_upoly_ir` |
| `factorize<polynomial_<ZZ, lex_<less>>>` | `polynomial_factorize.hh:22` | `template<class var_order>`（具体化到 `var_order=less`） | `factorize_lex_ir` |
| `factorize<polynomial_<ZZ, grlex_<less>>>` | `polynomial_factorize.hh:133`（§8.6） | `template<class comp>`（具体化到 `comp=grlex_<less>`） | `factorize_grlex_ir` |

**3 个源 body 结构不同**：
- `factorize_upoly_ir`：直接处理单变量 ZZ
- `factorize_lex_ir`：多变量 lex 序的完整 pipeline（sqf → primitive → 若单变量则 delegate 给 upoly，否则调 `__factor_multivar`）
- `factorize_grlex_ir`：**转换器**——把输入从 `grlex` 转到 `lex`，调用 `factorize_lex_ir`，再把结果转回 `grlex`

### §5.3 其他 64 函数的单实例化

**都是因为其他函数只被 1 个具体的类型路径调用**。示例（均 1 实例）：

| 函数 | 仅存在的实例 | 实例化原因 |
|---|---|---|
| `__factor_multivar` | `<less>`（即 lex_<less>）| 仅被 `factorize_lex_ir` 调用 |
| `__wang_core` | `<less>` | 仅被 `__factor_multivar<less>` 调用 |
| `__mtshl_lift` | `<less>` | 同上链 |
| `__factor_squarefree_primitive_ZZ` | 非模板 | 单变量 ZZ 专用 |
| `__hensel_lift` | `<less>` | |
| `__ddf_Zp`、`__squarefree_Zp`、`__factor_Zp` | 非模板 | Zp 专用 |

**`factorize_grlex_ir` 不创建新的 `__factor_multivar<grlex>` 实例**，因为它**先把 grlex 转 lex**，然后才调 `__factor_multivar<lex>`。所以没有 `grlex` 版本的下游链。

### §5.4 HIR 处理方案

#### §5.4.1 实例识别

`parse` Pass 对每个函数：
1. 从 Clang AST JSON 找所有 mangledName 非空的 `FunctionDecl`（`enumerate_instances.py` 已示范）
2. 对 TRANSLATION_SCOPE 的每个函数名：
   - 若恰 1 个实例 → 生成 1 个 `_ir` Lean def
   - 若多个实例（仅 `factorize` 一例）→ 按实例的 qualType 签名生成多个 `_ir`，命名用类型后缀

#### §5.4.2 命名后缀规则

```python
def instance_suffix(qualType: str) -> str:
    """根据实例的 C++ 类型签名生成 Lean 名称后缀。"""
    # factorize<upolynomial_<ZZ>>
    if "upolynomial_" in qualType: return "upoly"
    # factorize<polynomial_<ZZ, lex_<less>>>
    if "lex_<less>" in qualType: return "lex"
    # factorize<polynomial_<ZZ, grlex_<less>>>
    if "grlex_<less>" in qualType: return "grlex"
    # 其他情况（当前未出现）
    return hash_suffix(qualType)
```

生成 Lean 定义名：`{base_name}_ir` 若单实例，`{base_name}_{suffix}_ir` 若多实例。

#### §5.4.3 子调用解析

在实例化版本的 AST 里，所有子函数调用**已经指向具体的实例化版本**（Clang 完成名字查找 + 实例化）。HIR 的 `parse` Pass 直接用：

```python
# 在 factorize_lex_ir 的 body 里，Clang AST 的 CallExpr 指向 __factor_multivar<less>
# mangledName: _ZN6clpoly17__factor_multivarINS_4lessEE...
# parse Pass 根据 mangledName 映射到 Lean 名 __factor_multivar_ir
# （因为只有 1 个实例，不用后缀）

# 在 factorize_grlex_ir 的 body 里，AST 的 CallExpr 指向 factorize<less>
# （即 lex 版本的递归调用）
# parse Pass 映射到 factorize_lex_ir
```

**关键机制**：翻译器用 **mangledName → Lean 名** 的 dict，构建时一次性填充，查找一次 O(1)。

#### §5.4.4 HIR 节点扩展

`HIRFunc` 增加一个字段：

```python
@dataclass
class HIRFunc:
    base_name: str          # 如 "factorize"
    instance_suffix: str     # 如 "lex" / "grlex" / "upoly" / ""
    mangled_name: str
    qual_type: str          # 完整签名字符串（用于 debug / 保留）
    lean_name: str          # 最终生成 "factorize_lex_ir"
    params: list[HIRParam]
    body: HIRBlock
```

#### §5.4.5 实例列表生成

翻译器**运行前**用 `enumerate_instances.py` 脚本产出完整实例表：

```
TRANSLATION_INSTANCES = {
    "factorize": [
        {"suffix": "upoly",  "mangled": "_Z...", "qualType": "..."},
        {"suffix": "lex",    "mangled": "_Z...", "qualType": "..."},
        {"suffix": "grlex",  "mangled": "_Z...", "qualType": "..."},
    ],
    "__factor_multivar": [
        {"suffix": "", "mangled": "_Z...", "qualType": "..."},
    ],
    # ... 其他 63 个 ...
}
```

这张表由调研脚本生成，翻译器查表确定每个函数的实例数与 Lean 名。

### §5.5 Mathlib 精化映射（§8.5-8.8 对齐）

| C++ 实例 | L1 Lean | L2 精化目标 |
|---|---|---|
| `factorize<upolynomial_<ZZ>>` | `factorize_upoly_ir` | 单变量 `Polynomial ℤ` 的 Mathlib 因式分解 |
| `factorize<polynomial_<ZZ, lex>>` | `factorize_lex_ir` | `MvPolynomial (Fin n) ℤ` with `lex` 序 |
| `factorize<polynomial_<ZZ, grlex>>` | `factorize_grlex_ir` | `MvPolynomial (Fin n) ℤ` with `grlex` 序 |

但 L2 模型（已完成 6276 行 0 sorry）**只关心数学正确性**（不可约因子的乘积等于输入），与 monomial order 无关。因此：
- L2 的 `factorize_correct` 定理对 `factorize_lex_ir` / `factorize_grlex_ir` **通用**
- 精化证明需要证明 "grlex 版本 = lex 版本 ∘ 转换" 即可直接继承正确性

---

## §6 辅助 struct（已在 §3.3.9 覆盖）

见 §3.3.9：`PrimeSelectionResult`、`WangLcResult`、`HenselNode`、`Factorization` 的 Lean struct 定义。

这 4 个结构是算法内部状态容器，由 `parse` Pass 直接生成对应 Lean struct。不需要精化到 Mathlib 类型（L2 直接用相同 struct 或拆解字段）。

---

## §7 引用与指针消除

### §7.1 参数模式统计（`ref_params.md`）

| 传递方式 | 数量 | 比例 | v2 策略 |
|---|---|---|---|
| `T` (value) | 40 | 21% | 保留（Lean 值传递） |
| `const T&` | 125 | 65% | 保留为值参数（Lean 值语义，引用透明） |
| **`T&`** (non-const ref) | **26** | **14%** | **`ref_elim` Pass：提升为返回值 tuple** |
| `T*` (pointer) | 1 | 0.5% | 个案处理 |
| `T&&` (rvalue ref) | 0 | — | 不使用 |

### §7.2 `ref_elim` Pass 规则

对每个 C++ 函数：
1. 扫描参数表，识别所有 `T&`（含 `const T&` 排除）
2. 对每个 ref 参数，将其从参数列表**移到返回值 tuple**
3. 在函数体末尾自动追加该变量的最终 SSA 版本

**转换前**：
```cpp
void foo(ZZ& result, SparsePoly& work, int n) {
    // 修改 result 和 work
    result = 42;
    work.push_back(ZZ(1));
}
```

**转换后**（HIR₁）：
```
def foo_ir (result_0 : ZZ) (work_0 : SparsePoly) (n : Int) : ZZ × SparsePoly :=
    let result_1 := 42
    let work_1 := work_0.push (ZZ.ofInt 1)
    (result_1, work_1)
```

### §7.3 自动推导 vs 配置表

**v1 问题**：`TRANSLATION_SCOPE_OUTPUT_PARAMS` 是手工维护的 dict，8 个函数配置错位（`ref_params.md` §4 发现）。

**v2 方案**：`ref_elim` Pass 从 AST 的 `ParmVarDecl.type.qualType` **直接读取**是否带 `&`，不再依赖配置。零手工维护错误源。

### §7.4 调用方重写

对调用 `foo(result, work, n)` 的每个调用点：

**转换前**（C++）：
```cpp
ZZ result;
SparsePoly work;
foo(result, work, n);
use(result, work);
```

**转换后**（HIR₁）：
```
let result_0 : ZZ := default
let work_0 : SparsePoly := #[]
let (result_1, work_1) := foo_ir result_0 work_0 n
use(result_1, work_1)
```

注意：调用方必须初始化 ref 参数的初值（或用 `default`）。

### §7.5 1 个指针参数的特殊处理

唯一的 `T*` 参数位于（需 AST 精确定位）：

从 `ref_params.md` 数据推断，最可能是：
- `__factor_squarefree_primitive_ZZ` 的内部某个 `mpz_t*` 传参（GMP 原生）
- 或 `poly_convert` 的 `comp_ptr()`（类型参数）

**策略**：把 `T*` 当作 `T`（Lean 值传递，ptr 语义不翻译）。如果实际行为依赖于指针别名，手工在 `clpoly_model.lean` 里用 opaque 函数包装。

### §7.6 `std::move` 的交互

CLPoly 代码常写 `foo(std::move(x), ...)`，但对**输出参数**（ref）不会 move（move 是对 rvalue 的优化）。对**输入 const ref 参数**的 `std::move` 在 Lean 中等同于 `id`。

---

## §8 CAST_TABLE：完整类型转换表

### §8.1 Cast 种类分布（`casts.md`）

| castKind | 次数 | 说明 | Lean 处理 |
|---|---|---|---|
| `FunctionToPointerDecay` | 1786 | `func → &func` | **noop**（Lean 不区分函数和函数指针）|
| `NoOp` | 1630 | 类型同质 | **noop** |
| `LValueToRValue` | 1263 | 读取左值 | **noop**（Lean 值语义）|
| **`IntegralCast`** | **704** | **整数 ↔ 整数** | **查表** |
| `ConstructorConversion` | 142 | 构造函数转换 | 映射到 Lean 构造器 |
| `BuiltinFnToFnPtr` | 53 | builtin 到 fn ptr | noop |
| `ArrayToPointerDecay` | 46 | 数组到指针 | noop（CLPoly 不用指针算术）|
| `ToVoid` | 25 | 丢弃值 | `let _ := ...` |
| `UserDefinedConversion` | 16 | 用户定义 | 映射到 Lean 函数 |
| `IntegralToFloating` | 4 | 整数 → 浮点 | `Float.ofNat` / `Int.toFloat` |
| `FloatingToIntegral` | 1 | 浮点 → 整数 | `Float.toInt` / `Float.floor.toInt` |
| **合计** | **5670** | | |

**实际需要 CAST_TABLE 查表的只有**：`IntegralCast` (704) + `IntegralToFloating` (4) + `FloatingToIntegral` (1) + `UserDefinedConversion` (16) + `ConstructorConversion` (142，但多数是 `T(x)` 格式平凡) = **~870 处**。

其余 ~4800 处 cast 在 HIR 阶段**整体消除**（翻译为恒等映射）。

### §8.2 IntegralCast 完整映射表

基于 `casts.md` 的前 80 三元组，IntegralCast 的主要源/目对：

| Source | Target | 次数 | Lean 表达式 | UB 义务 |
|---|---|---|---|---|
| `int` | `size_type` (unsigned long) | 465 | `i.toNat.toUInt64` | `require i ≥ 0` |
| `int` | `int64_t` | 61 | `i.toInt64` | 无（int64 容纳 int）|
| `size_type` | `int` | 44 | `s.toInt32` | `require s ≤ INT_MAX` |
| `int` | `uint64_t` | 37 | `i.toNat.toUInt64` | `require i ≥ 0` |
| `int` | `size_t` | 26 | `i.toNat.toUInt` | `require i ≥ 0` |
| `uint64_t` | `int` | ~20 | `u.toInt32` | `require u ≤ INT_MAX`（UB-7）|
| `uint64_t` | `int64_t` | ~15 | `u.toInt64` | `require u ≤ INT64_MAX`（UB-7）|
| `int64_t` | `int` | ~10 | `i.toInt32` | `require -2^31 ≤ i < 2^31` |
| `int64_t` | `uint64_t` | ~5 | `i.toNat.toUInt64` | `require i ≥ 0` |
| `bool` | `int` | ~5 | `if b then 1 else 0` | 无 |
| `int` | `bool` | ~5 | `i ≠ 0` | 无 |

### §8.3 CAST_TABLE Python 结构

```python
CAST_TABLE: dict[tuple[str, str], str] = {
    # (source_type, target_type) → Lean 表达式模板
    #   "{x}" 是源值占位符
    #   (需前置条件 require)
    ("int", "size_type"):     "({x}).toNat.toUInt64",      # require x ≥ 0
    ("int", "int64_t"):       "({x}).toInt64",             # 无 UB
    ("size_type", "int"):     "({x}).toInt32",             # require x ≤ INT_MAX
    ("int", "uint64_t"):      "({x}).toNat.toUInt64",      # require x ≥ 0
    ("int", "size_t"):        "({x}).toNat.toUInt",        # require x ≥ 0
    ("uint64_t", "int"):      "({x}).toInt32",             # UB-7 require
    ("uint64_t", "int64_t"):  "({x}).toInt64",             # UB-7 require
    ("int64_t", "int"):       "({x}).toInt32",             # UB-6 require
    ("int64_t", "uint64_t"):  "({x}).toNat.toUInt64",      # require x ≥ 0
    ("bool", "int"):          "(if {x} then 1 else 0)",
    ("int", "bool"):          "({x} ≠ 0)",

    # IntegralToFloating (4 处，全在 __heuristic_starting_precision)
    ("int", "double"):        "Int.toFloat ({x})",
    ("uint64_t", "double"):   "Nat.toFloat ({x}).toNat",

    # FloatingToIntegral (1 处)
    ("double", "int"):        "({x}).toInt",  # 截断（C++ 是 truncate 不是 floor）

    # NoOp 类（映射恒等）
    ("int", "int"):           "{x}",
    ("uint64_t", "uint64_t"): "{x}",
    # ... 所有同型对 ...
}
```

### §8.4 UserDefinedConversion（16 处）

从 `casts.md` 的前 80 看，UserDefinedConversion 主要是：
- `ZZ(int literal)` — `ZZ` 构造函数：`Int.ofNat` / 直接字面量
- `Zp(val, prime)` — `Zp` 构造函数：`⟨val, prime⟩`
- `QQ(int)` — `Rat` 构造函数

**处理**：映射到 `FUNC_MAP`（CLASS_MAP）里的 Lean 构造器名。不需专门 CAST_TABLE 条目。

### §8.5 ConstructorConversion（142 处）

绝大多数是 `T(x)` 的显式构造，本质是 function call。映射同 UserDefinedConversion：查 `CLASS_MAP` → Lean 构造器。

### §8.6 CAST_TABLE 完整性验证

**covered**：IntegralCast（704）+ IntegralToFloating（4）+ FloatingToIntegral（1）= 709 处显式整数/浮点 cast，全部覆盖。

**noop**：FunctionToPointerDecay + NoOp + LValueToRValue + BuiltinFnToFnPtr + ArrayToPointerDecay + ToVoid = 4803 处 自动归零。

**查表**（`FUNC_MAP`）：ConstructorConversion + UserDefinedConversion = 158 处。

**总计 5670 cast 全部覆盖**。

### §8.7 v1 翻译器的不足

v1 的 `ssa_transform.py` 对 cast 处理分散在多处：
- `ImplicitCastExpr` 在 `clang_ast.py` 里 unwrap（大部分情况 OK）
- `CStyleCastExpr` / `CXXStaticCastExpr` 零散处理
- 没有 CAST_TABLE，整数类型转换靠运行时字符串拼接

**v2 改为统一 CAST_TABLE**，`operator_resolve` Pass 查表一次完成。消除"不同 cast kind 在不同代码路径被处理"的混乱。

---

## 当前进度（Week 2 完成）

---

## §6 Factorization / WangLcResult / PrimeSelectionResult struct（Day 4）

*（计划：CLPoly 结果类型的 Lean 结构化定义）*

---

## 当前进度（Week 2 完成）

- [x] **Day 1 (2026-04-21)**：§1 基础数值 + §2 UB 分析
- [x] **Day 2 (2026-04-21)**：§3 CLPoly 核心类型
- [x] **Day 3 (2026-04-21)**：§4 STL 容器 shim
- [x] **Day 4 (2026-04-21)**：§5 模板实例化 + §6（已并入 §3.3.9）
- [x] **Day 5 (2026-04-21)**：§7 引用指针消除 + §8 CAST_TABLE

## Week 2 验收

**本文件（`type-system.md`）作为 Week 2 的硬性产出**，覆盖：
- §1 所有基础数值类型 + §2 685 UB 站点 + §3 所有 CLPoly 核心类型（14 个）
- §4 全部 STL 容器 shim + §5 模板实例化方案 + §6 辅助 struct
- §7 引用/指针消除 + §8 CAST_TABLE 覆盖 5670 cast 全部

**给 Week 3（语义）/ Week 4-5（HIR/MIR 设计）提供完整输入**。
