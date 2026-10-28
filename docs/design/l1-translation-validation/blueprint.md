# C++ → Lean IR 翻译器蓝图

## 0. 目标

构建一个自动翻译器，将 CLPoly 的 C++ 因式分解代码（4,617 行）翻译为 Lean 4 定义，同时生成所有 UB-freedom 证明目标。翻译后的 Lean IR 与已有的 L2 算法模型之间的精化关系可以用 Lean 本身证明，从而完成 L1→L2 的验证链。

```
C++ 源码 ──Clang──→ AST JSON ──翻译器──→ Lean IR 定义
                                          ↓
                              UB proof obligations
                                          ↓
                              refinement proof: L1 IR ≈ L2 Algorithm
```

## 1. 架构

```
┌─────────────────────────────────────────────────────┐
│                    cpp2lean                          │
│                                                     │
│  Phase 1: Clang AST 提取                            │
│  clang++ -Xclang -ast-dump=json -fsyntax-only       │
│                                                     │
│  Phase 2: AST → SSA IR                              │
│  ├─ 变量重命名（消除 mutation → let 链）            │
│  ├─ 循环 → 尾递归                                   │
│  ├─ assert → Require 节点                           │
│  └─ 类型推断 + 溢出标注                             │
│                                                     │
│  Phase 3: SSA IR → Lean 4 代码生成                  │
│  ├─ 类型映射（uint64_t → UInt64）                   │
│  ├─ partial def 生成（循环翻译，不需要终止推断）    │
│  ├─ require 生成（UB + assert 前置条件）            │
│  └─ UB proof obligation 收集                        │
│                                                     │
│  Phase 4: 精化模板生成                              │
│  ├─ L1 IR 函数签名 → L2 函数签名对齐               │
│  └─ refinement theorem skeleton                     │
└─────────────────────────────────────────────────────┘
```

## 2. C++ → Lean 类型映射

### 2.1 基本类型

| C++ 类型 | Lean 精确模型 | Lean 抽象模型 | 说明 |
|---------|-------------|-------------|------|
| `uint64_t` | `UInt64` | `Nat` | 精确：mod 2^64 语义；抽象：无溢出假设 |
| `int64_t` | `Int64` | `Int` | 精确：补码语义 |
| `uint32_t` | `UInt32` | `Nat` | |
| `unsigned __int128` | `UInt128` (自定义) | `Nat` | Barrett 乘法核心 |
| `bool` | `Bool` | `Prop` | |

### 2.2 复合类型

| C++ 类型 | Lean 模型 | UB 检查 |
|---------|----------|--------|
| `T*` / `T&` | 不翻译（内联化） | — |
| `vector<T>` | `Array T` | 越界：`i < arr.size` |
| `pair<A,B>` | `A × B` | — |
| CLPoly `upolynomial_<Zp>` | `SparsePolyZp` (自定义) | 有序不变量 |
| CLPoly `basic_monomial` | `Monomial` (自定义) | — |

### 2.3 CLPoly 特有类型

```lean
/-- CLPoly 稀疏多项式：降序 (degree, coeff) 对列表 -/
structure SparsePolyZp where
  terms : Array (Nat × UInt64)
  -- 不变量：度数严格降序，系数非零
  inv_sorted : ∀ i j, i < j → j < terms.size →
    (terms[i]!).1 > (terms[j]!).1
  inv_nonzero : ∀ i, i < terms.size → (terms[i]!).2 ≠ 0

/-- 到 Mathlib Polynomial 的精化映射 -/
def SparsePolyZp.toPoly (s : SparsePolyZp) (p : Nat) :
    Polynomial (ZMod p) := ...
```

## 2a. 翻译器输入

翻译器有四个输入：

| 输入 | 性质 | 内容 |
|------|------|------|
| C++ 源码文件 | 翻译对象 | `polynomial_factorize_zp.hh` 等 |
| `clpoly_model.lean` | 可信基定义 | 类型+方法的 Lean 模型（手写，Lean 编译验证） |
| `class_map.py` | 操作映射（固定） | C++ 类操作 → Lean 函数名 |
| 翻译范围配置 | 每次翻译的输入 | 本次翻译哪些函数 |

**翻译范围**：显式列出被翻译的 C++ 函数集合。翻译器据此决定函数调用的链接方式：

```
调用 f()：
  f 在翻译范围内 → f_ir（调用翻译后的版本）
  f 在 FUNC_MAP 中 → 映射到 clpoly_model.lean 的函数
  f 在 LEAN_BUILTINS 中 → 直接使用 Lean 标准库
  都不在 → sorry
```

翻译范围是配置（每次可变），不是 class_map 的一部分（class_map 是固定知识）。

## 2b. 核心设计原则

**P0: 翻译器只查表，不硬编码。** 类操作查 CLASS_MAP，独立函数查 FUNC_MAP，函数链接查翻译范围。不在表中 → sorry。

**P1: 翻译器不做类型推断——所有类型信息从 Clang AST 传播。**

Clang AST 的每个表达式节点都有精确的 `type.qualType`。每个隐式类型转换都表示为显式的 `ImplicitCastExpr` 节点（含源类型和目标类型）。翻译器**直接使用这些信息**，不自己猜测或推断类型。

这意味着：
- 变量类型：从 `VarDecl` 的 `type.qualType` 获取，不从 RHS 推断
- 类型转换：从 `ImplicitCastExpr` / `CStyleCastExpr` 获取，不特判运算符
- 函数返回类型：从 `FunctionDecl` 的 `type.qualType` 获取
- UInt128 ↔ UInt64 转换：由 Clang 的 Cast 节点决定，翻译器只做 Cast → Lean 表达式的映射

**违反此原则的代码是 bug。** 如果翻译器需要"猜"一个类型，说明 Clang AST 解析不完整——应该修解析，不是加启发式。

所有隐式类型转换的处理流程：(1) 从 Clang AST 的 `ImplicitCastExpr` 读取 `(source_type, target_type)`；(2) 查 `CAST_TABLE` 映射到对应的 Lean 转换表达式。`CAST_TABLE` 是 C++ 类型转换语义到 Lean 表达式的静态映射（如 `(UINT128, UINT64) → uint128_lo`），不做推断。未在表中的转换 → `sorry`。

**P2: 翻译规则与 C++ 构造一一对应，不合并不拆分。**

每种 C++ 语句/表达式类型对应一条翻译规则。翻译器不优化、不合并、不重排。生成的 Lean 代码可能冗余（如 `(0 : UInt64).toNat.toUInt64`），但**忠实反映 C++ 语义**。优化是后续可选步骤。

**P3: 未知构造保留为 sorry，不丢弃不猜测。**

翻译器遇到未识别的 AST 节点时，生成 `UnknownStmt`/`UnknownExpr` → Lean 中输出 `sorry` + 注释。不尝试"近似翻译"。

## 3. 语句翻译规则

### 3.1 变量声明 + 赋值 → let 链（SSA 变换）

```cpp
// C++
uint64_t x = a + b;
x = x * c;           // mutation!
uint64_t y = x + 1;
```

```lean
-- Lean IR (SSA: 每次赋值创建新变量)
let x_0 : UInt64 := a + b
let x_1 : UInt64 := x_0 * c    -- x_0 renamed to x_1
let y : UInt64 := x_1 + 1
```

**规则**：每个赋值 `x = e` 生成 `let x_{n+1} := e[x_n/x]`。最终版本号在后续引用中使用。

### 3.2 三种 C++ 构造的翻译

C++ 代码中有三种需要区分的构造：

| C++ 构造 | 性质 | Lean 翻译 | 义务 |
|---------|------|----------|------|
| 隐式 UB（`a/b`, `arr[i]`, `a<<n`） | 无检查，静默出错 | `require` 前置条件 | **必须证明** |
| `assert(cond)` | 开发者认为必成立 | `require` 前置条件 | **必须证明** |
| `throw` / 显式错误 | 对非法输入的处理 | `Except.error` | 不证；条件正确性 |

```cpp
assert(p != 0);           // → require hp : p ≠ 0（必须证明）
if (f.empty()) throw ...;  // → if f.empty then .error（显式错误路径）
uint64_t r = a / b;       // → require hb : b ≠ 0（必须证明）
```

```lean
-- 函数签名：隐式 UB + assert 全部变成 require
-- 返回类型：有 throw 时用 Except，无 throw 时直接返回
partial def foo_ir (p : UInt64)
    (hp : p ≠ 0)                         -- assert(p != 0)
    (h_shift : norm.val < 64)            -- 隐式 UB: shift
    : Except Error UInt64 :=             -- 有 throw
  if f.empty then .error .emptyInput     -- throw
  else .ok (...)                         -- 正常路径，UB-free
```

### 3.3 if/else → Lean if

```cpp
if (r >= pn) r -= pn;
```

```lean
let r_1 := if r_0 ≥ pn then r_0 - pn else r_0
-- 注：uint64 减法无 UB（mod 2^64）；此处 if guard 已保证 r_0 ≥ pn
```

### 3.4 所有循环 → `partial def` 尾递归（统一规则）

**所有** `for`、`while`、含 `break` 的循环统一翻译为 `partial def` 尾递归。不使用 `Nat.fold`、不使用 `fuel`、不使用 `termination_by`。

**理由**：
- 终止性在 L2 已证明（6,276 行，0 sorry），L1 只关心与 L2 一致性
- `partial def` 可 `#eval`（背靠背测试）
- 翻译器不需要推断终止度量（参考 Aeneas, ICFP 2022）
- 一种规则覆盖所有循环模式，语义保持论证只做一次

**翻译规则**：

```cpp
// for 循环
for (int i = 0; i < n; i++) { body; }
// while 循环
while (cond) { body; }
// for + break
for (int i = 0; i < n; i++) { if (done) break; body; }
```

全部翻译为同一形式：

```lean
partial def loop_ir (i : Nat) (acc : State) : State :=
  if i ≥ n then acc                    -- 循环条件不满足 → 退出
  else if done then acc                -- break → 返回当前状态
  else loop_ir (i + 1) (body acc i)    -- 继续
```

嵌套循环 = 嵌套翻译，每层独立的 `partial def`：

```lean
partial def inner_ir (j : Nat) (acc : State) : State :=
  if j ≥ m then acc
  else if found j then acc             -- inner break
  else inner_ir (j + 1) (step acc j)

partial def outer_ir (s : Nat) (acc : State) : State :=
  if s > limit then acc
  else outer_ir (s + 1) (inner_ir 0 acc)
```

### 3.5 return → 最终表达式

```cpp
return r >> norm;
```

```lean
.ok (r >>> norm)    -- 有 throw 的函数返回 Except
-- 或
r >>> norm          -- 无 throw 的函数直接返回
```

### 3.6 函数调用 → 直接调用

```cpp
poly[i].coeff = nmod_mul(poly[i].coeff, lc_inv, p, ninv, norm);
```

```lean
let coeff_new := nmod_mul_ir (poly[i]!.coeff) lc_inv p ninv norm h_ub...
-- require: i < poly.size
```

## 4. UB 证明目标系统

### 4.1 分类

| 类型 | C++ 触发 | Lean 处理 | 义务 |
|-----|---------|----------|------|
| 除以零 | `a / b` | `require hb : b ≠ 0` | 必须证明 |
| 数组越界 | `arr[i]` | `require hi : i < arr.size` | 必须证明 |
| 无符号下溢 | `a - b` (unsigned) | 无 UB（mod 2^64 语义） | 不需要（但精化证明可能需要 `b ≤ a` 以匹配 L2） |
| 有符号溢出 | `a + b` (signed) | `require h : INT_MIN ≤ a+b ≤ INT_MAX` | 必须证明 |
| 移位越界 | `a << n` | `require h : n < 64` | 必须证明 |
| assert | `assert(cond)` | `require h : cond` | 必须证明 |
| throw | `throw error(...)` | `Except.error` | 不证明 |

### 4.2 输出格式

```lean
/-- C++ nmod_mul → Lean IR -/
partial def nmod_mul_ir (a b p ninv : UInt64) (norm : UInt32)
    (hp : p ≠ 0)                -- assert(p != 0) → 必须证明
    (h_shift : norm.val < 64)   -- 隐式 UB → 必须证明
    : UInt64 :=                 -- 无 throw → 直接返回
  let pn := p <<< norm.val
  let a_shifted := a <<< norm.val
  ...
```

## 5. 精化证明

### 5.1 正确性陈述

精化定理的一般形式：**如果没有 throw（.ok），则 L1 结果与 L2 一致。**

```lean
-- 无 throw 的函数：直接等式
theorem nmod_mul_refines (a b p ninv : UInt64) (norm : UInt32)
    (hp : p ≠ 0) (h_shift : norm.val < 64) :
    (nmod_mul_ir a b p ninv norm hp h_shift).val % p.val =
    (a.val * b.val) % p.val := by
  sorry -- Barrett reduction 正确性

-- 有 throw 的函数：条件等式
theorem factorize_refines (f : SparsePolyZp) (p : UInt64) (h_ub...) :
    factorize_ir f p h_ub... = .ok result →
    result.map toPoly = factorize_l2 (f.toPoly p.val) := by
  sorry
```

### 5.2 终止性

L1 不证终止。终止性由 L2 已有证明保证：
- `mtshl_invariant_terminates`（MTSHL 循环）
- `multi_bdp_terminates`（Taylor 循环）
- `ddf_correct`（DDF 循环，6 不变量）
- 等（L2 共 6,276 行，0 sorry）

L1 的 `partial def` 在精化证明中通过 L2 的终止性间接获得终止保证。

### 5.3 端到端

```lean
/-- 完整验证链：L1 → L2 → L3 → Spec -/
theorem factorize_end_to_end (f : SparsePolyZp) (p : UInt64) (h_ub...) :
    factorize_ir f p h_ub... = .ok result →
    FactorZpCorrect (f.toPoly p.val) ... (result.map toPoly) := by
  -- 步骤 1: L1 → L2 精化（factorize_refines）
  -- 步骤 2: L2 正确性（factor_Zp_instantiate，已证 0 sorry）
  sorry
```

## 5a. 语义保持论证

### 5a.1 CLPoly C++ 子集

CLPoly 因式分解代码使用的 C++ 子集：

| 用到的 | 没用到的 |
|-------|---------|
| `uint64_t`/`int` 算术 | 异常（除显式 throw） |
| `vector` 顺序访问 | `goto` / `setjmp` |
| `for`/`while`/`if`/`break`/`return` | 虚函数 / 函数指针 |
| 函数调用（静态分发） | 多线程 / 原子操作 |
| `assert` | 指针别名（值传递+const引用） |
| 显式 `throw` | 动态内存（内循环无 new） |

### 5a.2 翻译正确性证明

#### 定理（语义保持）

设 P 是 CLPoly C++ 子集（§5a.1）中的函数，T(P) 是翻译后的 Lean IR。设 σ 是满足所有 `require` 前置条件的输入状态（即 UB-free 且 assert 成立）。则：

1. 若 C++ 执行 P(σ) 正常返回值 v，则 T(P)(σ) 返回 `.ok v'`，且 v' 与 v 在类型映射下一致
2. 若 C++ 执行 P(σ) 触发 throw，则 T(P)(σ) 返回 `.error e`
3. 若 C++ 执行 P(σ) 触发 UB，则该执行不满足 `require` 前置条件（矛盾，不在定理范围内）

#### 证明

对翻译规则逐条论证。定义 C++ 语义 ⟦·⟧_C 和 Lean 语义 ⟦·⟧_L，证明 ⟦T(S)⟧_L = ⟦S⟧_C 对每种语句 S。

**引理 1（类型映射保持运算）：** 对 CLPoly 使用的每种运算，C++ 语义和 Lean 语义在 UB-free 条件下一致。

| C++ 运算 | C++ 语义 | Lean 翻译 | Lean 语义 | UB-free 条件 | 一致性 |
|---------|---------|----------|----------|------------|--------|
| `a + b` (uint64) | (a+b) mod 2^64 | `a + b` (UInt64) | (a+b) mod 2^64 | 无（uint64 加法无 UB） | 定义相同 |
| `a - b` (uint64) | (a-b) mod 2^64 | `a - b` (UInt64) | (a-b) mod 2^64 | 无（uint64 减法无 UB） | 定义相同 |
| `a * b` (uint64) | (a*b) mod 2^64 | `a * b` (UInt64) | (a*b) mod 2^64 | 无 | 定义相同 |
| `a / b` (uint64) | a ÷ b | `a / b` (UInt64) | a ÷ b | b ≠ 0 | 定义相同 |
| `a << n` (uint64) | a × 2^n（n < 64 时） | `a <<< n` (UInt64) | (a × 2^n) mod 2^64 | n < 64 | C++ 标准一致 |
| `a >> n` (uint64) | a ÷ 2^n | `a >>> n` (UInt64) | a ÷ 2^n | n < 64 | 定义相同 |
| `(unsigned __int128)a * b` | a×b（精确，128 位） | `uint128_mul a b` | a×b（精确，UInt128） | 无 | 需证 uint128_mul_correct |
| `arr[i]` | arr 的第 i 个元素 | `arr[i]!` (Array) | arr 的第 i 个元素 | i < arr.size | 定义相同 |
| `a < b` (uint64) | 布尔比较 | `a < b` (UInt64) | 布尔比较 | 无 | 定义相同 |
| `a == b` | 布尔等于 | `a == b` | 布尔等于 | 无 | 定义相同 |
| `a && b` | 短路与 | `a && b` (Bool) | 短路与 | 无 | Lean Bool.and 短路语义一致 |
| `!a` | 逻辑非 | `!a` (Bool) | 逻辑非 | 无 | 定义相同 |
| `a ? b : c` | 条件表达式 | `if a then b else c` | 条件表达式 | 无 | 语义一致 |
| `(uint64_t)x` (from int) | 截断/扩展 | `x.toUInt64` | 截断到 mod 2^64 | 无 UB（定义行为） | C++ 标准一致 |
| `s.field` | 成员访问 | `s.field` (Lean 结构体) | 投影 | 无 | 定义相同 |

**证明**：Lean 4 的 `UInt64` 定义为 `Fin (2^64)`，其加/减/乘/除/移位/比较运算与 C++ uint64_t 的定义完全一致（均为 mod 2^64 算术）。Bool 的 `&&`/`||` 在 Lean 中也是短路求值（`Bool.and` 先求左，左为 false 则不求右），与 C++ 一致。条件表达式 `a ? b : c` = `if a then b else c`，语义相同。整数类型转换（unsigned→unsigned）在 C++ 中是定义行为（mod 2^N 截断），Lean 的 `toUInt64` 等函数有相同语义。结构体成员访问是纯投影操作。UInt128 需要独立证明 `uint128_mul_correct`。数组访问在 `i < arr.size` 条件下语义一致。 ∎

**引理 2（SSA 变换保持语义）：** 在无别名条件下，变量 mutation 重命名为 let 链不改变程序语义。

设 C++ 代码为：
```
T x = e1;
x = e2(x);    // mutation
... use(x) ...
```

翻译后：
```
let x_0 = e1
let x_1 = e2(x_0)
... use(x_1) ...
```

**证明**：

**无别名前提**：CLPoly 因式分解代码通过值传递和 const 引用传参。在无别名条件下，对变量 x 的赋值 `x = e2(x)` 只影响 x 本身，不通过指针/引用影响其他变量。

**语义等价**：在 C++ 中，`x = e2(x)` 后 x 的值为 e2(旧 x)。在 Lean 中，`let x_1 = e2(x_0)` 后 x_1 的值为 e2(x_0)。由于 x_0 = 旧 x，两者相等。后续所有对 x 的引用在 C++ 中读取最新值 = e2(旧 x)，在 Lean 中引用 x_1 = e2(x_0) = e2(旧 x)。语义一致。

**别名情况**：若存在别名（`T& y = x; y = 5;` 也改变 x），则 SSA 变换不保持语义。CLPoly 因式分解代码不存在这种模式（已通过代码审计确认：所有函数参数为值传递或 const 引用，无可变引用别名）。

**形式化**：此性质等价于"在单线程、无别名的命令式程序中，SSA 变换是语义保持的"，这是 Appel 1998 定理的直接实例。 ∎

**引理 3（循环→尾递归保持语义）：** 结构化循环翻译为 `partial def` 尾递归函数保持语义。

设 C++ for 循环为：
```
State acc = init;
for (int i = 0; i < n; i++) {
    if (break_cond(i, acc)) break;
    acc = body(i, acc);
}
// 此处 acc = 最终值
```

翻译后：
```lean
partial def loop_ir (i : Nat) (acc : State) : State :=
  if i ≥ n then acc
  else if break_cond i acc then acc
  else loop_ir (i + 1) (body i acc)
```

**证明**：对循环执行步数 k 归纳。

**Base（k=0）**：循环条件 `i ≥ n` 立即成立（i = n），C++ 返回 acc。Lean `if i ≥ n then acc` 返回 acc。一致。

**Step（k→k+1）**：假设循环在第 k 步后的状态为 (i_k, acc_k)，且 loop_ir(i_k, acc_k) = C++ 从 (i_k, acc_k) 开始执行的结果（归纳假设）。

在第 k+1 步：
- **C++**：检查 `i_k < n`（是），检查 `break_cond(i_k, acc_k)`：
  - 若 true：退出循环，返回 acc_k
  - 若 false：执行 `acc_{k+1} = body(i_k, acc_k)`，继续 (i_k+1, acc_{k+1})
- **Lean**：计算 `loop_ir(i_k, acc_k)`：
  - `i_k ≥ n`？否 → `break_cond i_k acc_k`？
  - 若 true：返回 acc_k ✓
  - 若 false：返回 `loop_ir(i_k + 1, body i_k acc_k)` = 由归纳假设 = C++ 结果 ✓

**break**：break 在 C++ 中立即退出循环，返回当前 acc。在 Lean 中，`if break_cond then acc` 立即返回当前 acc。语义一致。

**嵌套循环**：内层循环翻译为独立的 `partial def`，外层调用内层。由引理 3 对内层成立 + 引理 3 对外层成立（将内层调用视为普通函数调用），嵌套保持语义。

**while 循环**：`while(cond(acc)){body}` 翻译为：
```lean
partial def while_ir (acc : State) : State :=
  if ¬cond acc then acc
  else while_ir (body acc)
```
与 for 循环的归纳论证完全一致：cond 扮演循环条件，body 扮演步进。归纳基础：cond 为 false 时两者都返回 acc。归纳步：cond 为 true 时两者都执行 body 然后继续。 ∎

**引理 3a（continue→循环顶部跳转保持语义）：**

C++ `for(i=0;i<n;i++){if(cond) continue; body}` 等价于 `for(i=0;i<n;i++){if(!cond) body}`。

翻译为：
```lean
partial def loop_ir (i : Nat) (acc : State) : State :=
  if i ≥ n then acc
  else if cond i acc then loop_ir (i + 1) acc    -- continue = 跳过 body，直接 i+1
  else loop_ir (i + 1) (body i acc)
```

**证明**：continue 跳过当前迭代的剩余部分，进入下一迭代。翻译为 `if cond then loop_ir (i+1) acc`（不执行 body，直接递归到 i+1）。与 C++ 的 continue 语义一致：跳过 body，执行 step（i++），检查条件。 ∎

**引理 3b（if/else 保持语义）：**

C++ `if(cond) S1 else S2` 翻译为 Lean `if cond then T(S1) else T(S2)`。

**证明**：C++ 和 Lean 的 if/else 都是：求值 cond → 若 true 执行 S1/T(S1)，若 false 执行 S2/T(S2)。由归纳假设 T(S1) 和 T(S2) 分别保持语义，组合保持语义。 ∎

**引理 3c（return 保持语义）：**

C++ `return e` 翻译为 Lean 表达式 `T(e)`（函数体的最终表达式）或 `.ok T(e)`（有 Except 时）。

**证明**：return 终止函数执行并返回 e 的值。Lean 函数体的最终表达式即为返回值。由引理 1，T(e) 与 e 的值在类型映射下一致。 ∎

**引理 4（assert→require 保持语义）：**

C++ 代码 `assert(P); body` 在 P 成立时的行为 = body。

Lean 翻译 `(h : P) → body` 在提供 P 的证明后的行为 = body。

**证明**：assert 在 P 成立时是 no-op（release 编译直接忽略）。`require h : P` 将 P 作为前置条件，在 P 成立时函数体直接执行。两者在 P 成立时行为一致。在 P 不成立时：C++ assert 触发 abort（不在 UB-free 定理范围内，因为我们要求证明 assert 成立），Lean require 无法被调用。 ∎

**引理 5（throw→Except 保持语义）：**

C++ 代码 `if (error_cond) throw E; body` 的行为：
- error_cond 为 true：抛出异常 E，程序终止
- error_cond 为 false：执行 body

Lean 翻译 `if error_cond then .error E else .ok (body)` 的行为：
- error_cond 为 true：返回 `.error E`
- error_cond 为 false：返回 `.ok (body 的结果)`

**证明**：两者的分支条件相同（error_cond）。正常路径（.ok）的值 = C++ body 的返回值（由引理 1-4 对 body 成立）。错误路径的标记一致。 ∎

**引理 6（输出参数→返回值保持语义）：**

C++ 函数通过引用参数返回结果的两种模式：

**模式 A（单次赋值）**：`void f(T& out, ...) { ... out = result; }`

翻译为 `def f_ir (...) : T := result`。在无别名条件下，out 的最终值 = result，等于 Lean 返回值。

**模式 B（累积修改）**：`void f(vector<T>& out, ...) { ... out.push_back(x); ... out.push_back(y); }`

翻译为：out 初始值作为输入参数，每次 `push_back` 在 SSA 中生成新版本：
```lean
partial def f_ir (out_0 : Array T) ... : Array T :=
  let out_1 := out_0.push x
  ...
  let out_2 := out_1.push y
  out_2    -- 返回最终版本
```

**证明**：模式 A 是模式 B 的特例（只有一次赋值）。模式 B 由引理 2（SSA 变换）直接覆盖：每次 `push_back` 是对 out 的 mutation，SSA 变换将其重命名为 out_0, out_1, out_2...，最终版本作为返回值。在无别名条件下（out 不与其他参数共享），语义一致。 ∎

**引理 7（跨函数调用保持语义）：**

设函数 g 调用函数 f：`result = f(args)`。若引理 1-6 对 f 的翻译成立（即 `f_ir` 与 C++ 的 f 语义一致），则 g 中调用 f 的翻译也保持语义。

**证明**：g 的翻译中，`f(args)` 被翻译为 `f_ir args h_ub...`。由 f 的语义保持（归纳假设），`f_ir args` 的返回值与 C++ 中 `f(args)` 的返回值一致。g 的后续代码使用该返回值，由引理 2（SSA）语义保持。

**require 传播**：f 的 `require` 前置条件成为 g 调用 f 时的证明义务。翻译器在 g 中生成对应的 `h_ub` 参数或 `have` 声明。若 g 的 `require` 前置条件蕴含 f 的 `require`，则 g 调用 f 时前置条件自动满足。

**Except 传播**：若 f 返回 `Except`，g 中的 `let result ← f_ir args` 在 f 返回 `.error` 时传播错误（do notation 的 bind 语义），在 `.ok` 时继续。这与 C++ 中 throw 沿调用栈传播的语义一致。 ∎

**主定理证明**：

对程序的**调用图**（call graph）拓扑序归纳。叶函数（不调用其他翻译函数）的正确性由引理 1-6 保证。非叶函数的正确性由引理 7（函数调用）+ 引理 1-6（函数体内语句）组合保证。

对每个函数体的**语句结构**归纳。语句序列 `S1; S2` 翻译为 `let x = T(S1); T(S2)[x/result_of_S1]`：
- 若 T(S1) 保持语义（归纳假设），则 x 与 C++ S1 的结果一致
- T(S2) 中引用 x（= S1 的结果），由引理 2（SSA）语义保持
- 组合：整个序列保持语义

由调用图归纳 + 语句归纳的组合，整个程序的翻译保持语义。 ∎

#### 前提条件总结

翻译正确性依赖以下前提，均需通过代码审计或工具验证：

| 前提 | 验证方式 |
|------|---------|
| P1: CLPoly 因式分解代码无指针别名 | 代码审计：所有参数为值传递或 const 引用；无可变全局变量。可用 Clang `-Wstrict-aliasing` 辅助检查 |
| P2: CLPoly 因式分解代码无 goto/setjmp | 代码审计 + `grep -r "goto\|setjmp\|longjmp" *.hh`（已确认无结果） |
| P3: CLPoly 因式分解代码无隐式异常 | 代码审计：无 `new`（不触发 `bad_alloc`）；无 dynamic_cast；vector 使用 `[]` 不是 `.at()`（不抛异常）。显式 throw 保留为 Except.error |
| P4: CLPoly 因式分解代码无多线程 | 代码审计：因式分解函数不使用 `std::thread`/`std::mutex`/`std::atomic` |
| P5: CLPoly 因式分解代码无 `continue` 或仅在简单 for 中使用 | 代码审计（若有 continue 则由引理 3a 覆盖） |
| P6: Lean UInt64 = C++ uint64_t（mod 2^64 语义） | Lean 4 源码定义 `UInt64 := Fin (2^64)` + C++17 标准 §6.8.1（unsigned 运算 mod 2^N） |
| P7: uint128_mul 实现正确 | 独立 Lean 证明（uint128_mul_correct）+ 背靠背测试 |
| P8: Clang AST JSON 忠实反映 C++ 语义 | 信任 Clang 编译器（工业级工具，非我们验证范围） |

### 5a.3 参考文献

- Appel, "SSA is Functional Programming", 1998 — SSA→函数式的理论基础
- Ho & Protzenko, "Aeneas: Rust Verification by Functional Translation", ICFP 2022 — Rust MIR→Lean 4 翻译器（循环用 partial fixpoint，我们采用相同策略）

## 6. 分阶段实施计划

### Phase 0: 基础设施（1-2 周）

- [ ] Clang AST JSON 解析器（Python）
- [ ] SSA 变换引擎（mutation → let 链）
- [ ] 基本类型映射（uint64_t, int, bool）
- [ ] 简单函数翻译（纯算术，无循环）
- [ ] **验收**：`nmod_mul` 自动翻译，编译通过

### Phase 1: 控制流（2-3 周）

- [ ] if/else → Lean if
- [ ] 所有循环 → partial def 尾递归（统一规则）
- [ ] break → 提前返回
- [ ] assert → require
- [ ] throw → Except.error
- [ ] 嵌套函数调用
- [ ] **验收**：`upoly_make_monic` 自动翻译

### Phase 2: CLPoly 数据结构（2-3 周）

- [ ] `SparsePolyZp` 类型定义 + `toPoly` 精化映射
- [ ] vector 操作 → Array 操作 + 越界检查
- [ ] pair/tuple 翻译
- [ ] 模板单态化（`upolynomial_<Zp>` → `SparsePolyZp`）
- [ ] **验收**：`__squarefree_Zp` 自动翻译骨架

### Phase 3: 完整管线（3-4 周）

- [ ] 跨函数调用图处理
- [ ] UB 证明目标批量生成
- [ ] 精化定理骨架自动生成
- [ ] 背靠背测试框架
- [ ] **验收**：`polynomial_factorize_zp.hh` 全文件翻译 + 背靠背测试通过

### Phase 4: 精化证明（持续）

- [ ] Barrett reduction 精化（nmod_mul → Zp 乘法）
- [ ] 稀疏多项式精化（SparsePolyZp → Polynomial）
- [ ] 各算法子过程精化（SQF, DDF, EDF, Hensel, ...）
- [ ] 端到端精化（L1 factorize = L2 factorize）

## 7. 关键难点分析

### 7.1 128 位乘法建模

CLPoly 的 Barrett reduction 使用 `unsigned __int128`。Lean 没有原生 128 位整数。

**方案**：定义 `UInt128` 为 `(hi : UInt64) × (lo : UInt64)`，实现乘法/移位操作，证明与 `Nat` 乘法一致。

```lean
structure UInt128 where
  hi : UInt64
  lo : UInt64

def UInt128.toNat (x : UInt128) : Nat :=
  x.hi.toNat * 2^64 + x.lo.toNat

def uint128_mul (a b : UInt64) : UInt128 := ...

theorem uint128_mul_correct (a b : UInt64) :
    (uint128_mul a b).toNat = a.toNat * b.toNat := by
  sorry -- 需要位运算推理
```

### 7.2 指针/引用消除

CLPoly 大量使用引用参数和指针。翻译器需要：
1. **引用参数**：转为返回值（输出参数 → 多返回值）
2. **数组指针**：转为 Array + index
3. **this 指针**：转为显式 self 参数

```cpp
void foo(vector<int>& out, const vector<int>& in) { ... }
```
```lean
def foo_ir (input : Array Int) : Array Int := ...
```

### 7.3 模板单态化

CLPoly 使用模板泛型（`upolynomial_<T>`），但因式分解只用 `Zp` 和 `ZZ` 实例。

**方案**：Clang AST 已经单态化模板。直接从单态化后的 AST 翻译。

### 7.4 循环终止性

**不是翻译器的职责。** 所有循环统一用 `partial def`，不需要终止推断。终止性由 L2 层的已有证明保证（6,276 行，0 sorry）。L1 精化证明中，通过 L2 的终止性间接获得 L1 的终止保证。

### 7.5 副作用和全局状态

CLPoly 的因式分解模块基本无全局状态（函数式风格）。少数例外：
- `random()` 调用 → 参数化（与 EDF splits 处理相同）
- 输出参数 → 返回值

## 8. 工具链

```
Python 3.10+
├── clang_ast.py       -- Clang AST JSON 解析
├── ssa_transform.py   -- C++ AST → SSA IR
├── lean_codegen.py    -- SSA IR → Lean 4 代码
├── ub_collector.py    -- UB 证明目标收集
├── refine_template.py -- 精化定理骨架生成
└── cpp2lean.py        -- 主入口（组合上述模块）

依赖：
- clang++ 18+ (JSON AST dump)
- Python 3.10+ (match/case)
- Lean 4 + Mathlib (输出验证)
```

## 9. 成功标准

| 阶段 | 标准 |
|------|------|
| Phase 0 | `nmod_mul` 自动翻译 → Lean 编译通过 |
| Phase 1 | 10+ 个 CLPoly 函数自动翻译 → Lean 编译通过 |
| Phase 2 | `__squarefree_Zp` 完整翻译骨架 → Lean 编译通过（带 sorry） |
| Phase 3 | 整个 `polynomial_factorize_zp.hh` 翻译 → 编译通过 |
| Phase 4 | 至少 1 个子过程端到端精化证明（0 sorry） |

## 10. 与现有 L2 验证的集成

翻译器产出的 L1 IR 与现有 L2 模型之间的精化关系：

```
C++ 源码 ──翻译器──→ L1 Lean IR (SparsePolyZp, UInt64)
                          ↕ 精化证明
                     L2 Algorithm (Polynomial (ZMod p), Nat)
                          ↕ 已完成 (6,276 行, 0 sorry)
                     L3 Math (Mathlib)
                          ↕ Lean 内核
                     Spec (FactorZpCorrect, FactorZZCorrect)
```

精化证明的核心定理：
```lean
theorem l1_refines_l2 :
    ∀ f : SparsePolyZp, ∀ p : UInt64,
      (factorize_zp_ir f p).map (·.toPoly p.val) =
      factorize_zp_l2 (f.toPoly p.val) := by ...
```

这将闭合 C++ → Lean → 数学规约 的完整验证链。

---

## 附录 A：模块功能规约（§3.1 格式）

### A.1 clang_ast — Clang AST JSON 解析器

```
模块名称：clang_ast

功能描述：解析 Clang 输出的 AST JSON，提取函数声明、
         参数列表、函数体语句树，构建内部 IR 表示。

前置条件（Requires）：
  - 输入是合法 JSON（由 clang++ -Xclang -ast-dump=json 产生）
  - C++ 源码已通过 clang 编译（-fsyntax-only 无错误）

后置条件（Ensures）：
  - 输出 FuncIR 列表，每个 FuncIR 包含：
    - name : str（函数名）
    - params : List[ParamIR]（参数名+类型）
    - ret_type : TypeIR（返回类型）
    - body : StmtIR（函数体语句树）
  - 所有 C++ 语句/表达式 kind 被映射到 StmtIR/ExprIR 枚举
  - 未识别的 kind → UnknownStmt(kind, children)（不丢失信息）

不变式：无（无状态，纯函数）
副作用：无
```

### A.2 ssa_transform — SSA 变换引擎

```
模块名称：ssa_transform

功能描述：将 clang_ast 输出的 IR 变换为 SSA 形式：
         消除变量 mutation（每次赋值创建新版本），
         将循环转为尾递归，将 assert 转为 Require 节点。

前置条件（Requires）：
  - 输入是 clang_ast 产出的 FuncIR
  - 函数体内无 goto、setjmp/longjmp、无隐式异常（显式 throw 保留为 Throw 节点）

后置条件（Ensures）：
  - 输出 SSAFunc，每个变量只赋值一次（let 绑定）
  - 每个变量有唯一版本号（x_0, x_1, ...）
  - 所有后续引用使用最新版本号
  - assert(cond) → Require(cond)（保留，不过滤）
  - for(init;cond;step){body} → TailRec(init, cond, step, body)
  - while(cond){body} → TailRec(unit, cond, unit, body)
  - 输出参数（T& out）→ 返回值

不变式：
  - 语义保持：SSA 形式与原 C++ 在无 UB 时行为一致

副作用：无
```

### A.3 lean_codegen — Lean 4 代码生成器

```
模块名称：lean_codegen

功能描述：将 SSA IR 翻译为 Lean 4 源代码。
         生成可执行（computable）的 def 定义。

前置条件（Requires）：
  - 输入是 ssa_transform 产出的 SSAFunc
  - 所有类型已解析（无 UnknownType）

后置条件（Ensures）：
  - 输出 .lean 文件内容（字符串）
  - 每个函数生成 `partial def <name>_ir (params) : RetType := ...`
  - 所有 def 可 #eval（partial def 可执行），不使用 noncomputable
  - 类型映射：uint64_t → UInt64, int → Int, bool → Bool
  - Require(cond) → 函数参数中的 (h : cond) 前置条件（UB + assert）
  - throw → Except.error（返回类型包 Except）
  - TailRec → partial def 尾递归（不需要 termination_by）
  - UB 检查点以 require 参数形式生成

不变式：
  - 生成的 Lean 代码在 Lean 4 + Mathlib 下可编译（无语法错误）
  - 可执行：所有函数可用 #eval 运行

副作用：无
```

### A.4 ub_collector — UB 证明目标收集器

```
模块名称：ub_collector

功能描述：遍历 SSA IR，识别所有可能的未定义行为（UB）点，
         生成对应的 Lean 证明目标。

前置条件（Requires）：
  - 输入是 ssa_transform 产出的 SSAFunc

后置条件（Ensures）：
  - 输出 List[UBObligation]，每个包含：
    - location : (func_name, line_no)
    - ub_type : enum (DivByZero | ArrayOOB | UnsignedUnderflow |
                      SignedOverflow | ShiftOOB | AssertFail)
    - lean_prop : str（Lean Prop 表达式）
    - context : List[str]（相关 let 绑定上下文）
  - 完整性：每个 C++ UB 点对应恰好一个 obligation
  - 输出可直接嵌入 Lean 文件作为 theorem 或 have 声明

不变式：无
副作用：无
```

### A.5 refine_template — 精化定理骨架生成器

```
模块名称：refine_template

功能描述：对比 L1 IR 函数签名和 L2 算法函数签名，
         生成精化定理（refinement theorem）的 Lean 骨架。

前置条件（Requires）：
  - L1 IR 函数列表（来自 lean_codegen 输出）
  - L2 算法函数列表（来自 CLPoly/Algorithm/*.lean）
  - 类型映射表（SparsePolyZp.toPoly → Polynomial (ZMod p) 等）

后置条件（Ensures）：
  - 输出 List[RefinementTheorem]，每个包含：
    - l1_func : str（L1 函数名）
    - l2_func : str（L2 函数名）
    - theorem_stmt : str（Lean theorem 签名 + sorry 体）
    - type_bridge : str（toPoly 等精化映射）
  - 骨架可直接放入 .lean 文件编译（所有 sorry）

不变式：无
副作用：无
```

### A.6 back2back — 背靠背测试框架

```
模块名称：back2back

功能描述：对翻译后的 Lean IR 和原始 C++ 执行相同输入，
         比较输出一致性。轻量级翻译正确性验证，
         在精化证明之前捕获翻译 bug。

前置条件（Requires）：
  - L1 Lean IR 函数可执行（computable def，可 #eval）
  - C++ 函数可编译执行
  - 测试向量集（手动 + 随机生成）

后置条件（Ensures）：
  - 对每个测试向量 (input, expected_output)：
    - C++ 执行结果 == expected_output
    - Lean #eval 结果 == expected_output
    - 两者一致 ↔ 翻译正确（在该输入上）
  - 输出测试报告：通过/失败/不一致的用例列表

不变式：无
副作用：文件 I/O（生成临时 .lean 文件供 #eval）
```

---

## 附录 B：模块间接口规约（§3.2 格式）

### B.1 Clang → clang_ast

```
接口：Clang CLI → clang_ast

输入数据：JSON 字符串（Clang AST dump，根节点 kind=TranslationUnitDecl）
输出数据：List[FuncIR]（内部 IR 表示）

协议约定：
  - 调用方责任：确保输入是合法 Clang AST JSON（clang 无编译错误）
  - 被调用方责任：解析所有 FunctionDecl 节点；
    未识别的 AST kind 包装为 UnknownStmt/UnknownExpr（不丢弃）
```

### B.2 clang_ast → ssa_transform

```
接口：clang_ast → ssa_transform

输入数据：FuncIR（函数名 + 参数 + 类型 + 语句树）
输出数据：SSAFunc（SSA 变换后的函数，变量单赋值）

协议约定：
  - 调用方责任：FuncIR 中的语句树结构完整（无 None 子节点）
  - 被调用方责任：
    - 变量 mutation 消除（x = e → let x_{n+1} = e）
    - 循环转尾递归（for → TailRec，while → TailRec）
    - assert 保留为 Require 节点
    - 输出参数转返回值
```

### B.3 ssa_transform → lean_codegen

```
接口：ssa_transform → lean_codegen

输入数据：SSAFunc（SSA 形式函数）
输出数据：str（Lean 4 源代码）

协议约定：
  - 调用方责任：SSAFunc 中所有变量单赋值，类型已标注
  - 被调用方责任：
    - 生成可编译的 Lean 4 def
    - 所有 def 可执行（computable，可 #eval）
    - Require → 函数参数前置条件（require）
    - throw → Except.error
    - TailRec → partial def 尾递归（不需要 termination_by）
```

### B.4 ssa_transform → ub_collector

```
接口：ssa_transform → ub_collector

输入数据：SSAFunc（SSA 形式函数）
输出数据：List[UBObligation]（UB 证明目标）

协议约定：
  - 调用方责任：SSAFunc 结构完整
  - 被调用方责任：
    - 每个 UB 点（除法、数组访问、减法、移位、assert）恰好一个 obligation
    - obligation 的 lean_prop 是合法 Lean Prop 表达式
    - 不遗漏（完整性）、不重复
```

### B.5 lean_codegen + L2 → refine_template

```
接口：lean_codegen 输出 + L2 文件 → refine_template

输入数据：
  - L1 函数签名列表（从生成的 .lean 文件提取）
  - L2 函数签名列表（从 CLPoly/Algorithm/*.lean 提取）
  - 类型映射配置（TOML/JSON）

输出数据：Lean theorem 骨架（.lean 文件，全 sorry）

协议约定：
  - 调用方责任：L1 和 L2 函数有明确对应关系（配置文件指定）
  - 被调用方责任：
    - 生成的 theorem 签名类型正确
    - sorry 体可编译
    - 包含 toPoly 等精化映射的 import
```

### B.6 lean_codegen + C++ → back2back

```
接口：lean_codegen 输出 + C++ 源码 → back2back

输入数据：
  - Lean IR .lean 文件（可执行 def）
  - C++ 可执行文件（编译后）
  - 测试向量（JSON 格式：输入参数 + 期望输出）

输出数据：测试报告（JSON：每个用例的 pass/fail + diff）

协议约定：
  - 调用方责任：
    - Lean IR 可 #eval（computable）
    - C++ 可编译执行
    - 测试向量在 UB-free 范围内（不触发未定义行为）
  - 被调用方责任：
    - 对每个测试向量，分别执行 C++ 和 Lean
    - 比较输出（按类型：整数精确比较，多项式按系数比较）
    - 输出不一致时报告差异详情
```

---

## 附录 C：背靠背测试设计

### C.1 可执行性要求

翻译出的 Lean IR **必须**是 computable：

```lean
-- ✅ 正确：使用 UInt64，可 #eval
def nmod_mul_ir (a b p ninv : UInt64) (norm : UInt32) : UInt64 := ...

#eval nmod_mul_ir 7 11 13 ...  -- 可执行

-- ❌ 错误：使用 Nat/noncomputable，不可执行
noncomputable def nmod_mul_abstract (a b p : Nat) : Nat := ...
```

**设计决策**：L1 IR 使用具体类型（`UInt64`, `Array`, `UInt128`），不使用抽象数学类型。这保证可执行性，且与 C++ 语义精确对应。精化证明连接到使用 `Nat`/`Polynomial` 的 L2 模型。

### C.2 测试流程

```
1. 生成测试向量
   ├── 手动：从 CLPoly test suite 提取关键用例
   ├── 随机：生成随机多项式 + 随机素数
   └── 边界：p=3, p=2^63-25, deg=0, deg=1000

2. C++ 执行
   ├── 编译 CLPoly 测试程序
   ├── 输入测试向量 → 输出结果
   └── 序列化为 JSON

3. Lean 执行
   ├── 生成 #eval 调用的 .lean 文件
   ├── lake env lean eval_test.lean → 输出结果
   └── 解析输出

4. 比较
   ├── 逐字段比较 C++ 和 Lean 输出
   ├── 不一致 → 翻译 bug（定位到具体函数）
   └── 全部一致 → 翻译在该输入集上正确
```

### C.3 测试层次

| 层次 | 粒度 | 目的 |
|-----|------|------|
| **单元** | 单个函数（nmod_mul, upoly_add） | 验证基本运算翻译 |
| **集成** | 子管线（SQF, DDF, EDF） | 验证函数间交互 |
| **端到端** | 完整因式分解 | 验证整体行为 |

### C.4 与精化证明的关系

```
背靠背测试（轻量、自动、不完备）
  → 捕获大多数翻译 bug
  → 建立翻译正确性的信心

精化证明（重量、手动/半自动、完备）
  → 数学保证翻译正确
  → 对所有输入成立

推荐工作流：
  1. 翻译函数
  2. 背靠背测试通过（快速反馈）
  3. 然后做精化证明（高保证）
```

### C.5 示例：nmod_mul 背靠背测试

```python
# test_vectors.json
[
  {"a": 7, "b": 11, "p": 13, "ninv": ..., "norm": 0, "expected": 12},
  {"a": 0, "b": 100, "p": 17, "ninv": ..., "norm": 0, "expected": 0},
  {"a": 18446744073709551557, "b": 2, "p": 18446744073709551557,
   "ninv": ..., "norm": 0, "expected": 2},
  ...
]
```

```lean
-- eval_nmod_mul.lean (auto-generated)
import CLPoly.L1.NmodMul

#eval nmod_mul_ir 7 11 13 <ninv> 0    -- expect 12
#eval nmod_mul_ir 0 100 17 <ninv> 0   -- expect 0
```

```
$ lake env lean eval_nmod_mul.lean
12
0
2
$ diff <(lean_output) <(cpp_output)
(empty = all match)
```
