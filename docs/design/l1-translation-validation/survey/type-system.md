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
- [x] **Day 2 (2026-04-21)**：§3 CLPoly 核心类型完成
- [ ] **Day 3**：§4 STL 容器 shim
- [ ] **Day 4**：§5 + §6 模板实例化 + struct 定义（§6 已在 §3.3.9 部分覆盖）
- [ ] **Day 5**：§7 + §8 引用指针 + CAST_TABLE
