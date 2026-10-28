# Mathlib4 多项式类型调研与 L1 映射策略

**日期**：2026-04-19  
**作者**：代码调研  
**目标**：为 CLPoly L1 IR 提供 Mathlib 类型映射决策表

---

## 执行摘要

| CLPoly C++ 类型 | Mathlib 现成支持 | L1 映射方案 | 说明 |
|----------------|-----------------|-----------|------|
| `polynomial_<T, order>` 多变量稀疏多项式 | ❌ MvPolynomial 用 Finsupp，结构不同 | **自定义 struct** | 保留 Array 稀疏表示 + 单调序；L2 精化到 MvPolynomial |
| `upolynomial_<T>` 单变量稀疏多项式 | ✓ 部分支持 | **L1 保留，L2 精化** | L1 用 List<(deg, coeff)>；toPoly → Polynomial (ZMod p) |
| `ZZ`（GMP） | ❌ Int 功能足但不是 GMP | **shim → Int** | L1 用 struct{ptr, type}；L2 精化到 Int |
| `QQ`（有理数） | ✓ Rat 直接可用 | **直接用 Rat** | L1 可直接映射；含 normalize |
| `Zp`（运行时素数） | ❌ ZMod p 素数在编译时 | **自定义 Zp struct** | L1 {coeff: Nat, prime: Nat}；L2 精化到 ZMod p.toNat |
| `basic_monomial<order>` | ❌ 无对应（Finsupp 是 σ →₀ ℕ） | **自定义 Monomial** | L1 保留 order 信息；L2 映射到 Finsupp (Fin n) ℕ |
| `variable` | ✓ Fin n 或 σ | **直接用 Fin n 或 Fin (n.succ)** | 变量索引 |
| `factorization<Poly>` | ✓ 部分 Multiset | **List (Poly × ℕ) 保留** | C++ Array 直接对应 List；L2 精化到 Multiset |

---

## 第 1 部分：单变量多项式

### 1.1 Mathlib `Polynomial R` 概览

**定义**（`Mathlib/Algebra/Polynomial/Basic.lean`）:
```lean
structure Polynomial (R : Type*) [Semiring R] where
  toFinsupp : AddMonoidAlgebra R ℕ
```

本质上是 `Finsupp ℕ R`，即 `ℕ →₀ R`（有限支持函数）。

**关键 API**：
- `coeff n p`：系数访问（支持多态）
- `natDegree p`：最高非零项的度数（ZMod p 时精确度数）
- `leadingCoeff p`：首项系数
- `C a`：常数多项式
- `X`：变量 X
- `degree p`：有序偶（∈ WithBot ℕ）

### 1.2 算术与除法

**加法、乘法**：环运算，非交换多项式直接支持。

**除法与模**（`Mathlib/Algebra/Polynomial/Div.lean`）:
- `divByMonic f g`（g monic）：精确商（Ring 足够）
- `modByMonic f g`：精确余项
- `modByMonic_add_div : f %ₘ g + g * (f /ₘ g) = f`

**GCD**（`Mathlib/RingTheory/...`）：
- 不存在 `Polynomial.gcd`
- 用 `EuclideanDomain.gcd`（对 `ZMod p[X]` 工作）
- 返回结果**不一定 monic**，需 `normalize` 后得到 monic
- `normalize_dvd_iff : normalize a ∣ b ↔ a ∣ b`（关键桥接）

### 1.3 与 L2 现用情况

**L2 用法**（见 `CLPoly/Algorithm/SquarefreeZp.lean`）:
- `Polynomial (ZMod p)`：直接用于 Yun 算法、DDF 等
- `EuclideanDomain.gcd` + `normalize`：SQF 分解循环
- `modByMonic_add_div`：精确除法
- `derivative`：导数（自动简化）
- `Monic` 谓词：检查首项系数 = 1

**L2 没有的**：多变量系数（那时用 `Polynomial (Polynomial ℤ)` 嵌套）

### 1.4 L1 映射决策

| C++ 概念 | L1 类型 | 精化到 L2 |
|---------|--------|----------|
| `upolynomial_<ZMod p>` | `struct { terms: List (Nat × Nat), prime }` | `Polynomial (ZMod p)` via `toPoly` |
| `degree` | `Nat`（简化） | `natDegree` |
| `leadingCoeff` | `Nat` | `leadingCoeff` |
| 多项式除法 | 分别实现除/余 | `divByMonic` / `modByMonic` |
| 最大公因子 | 单独 GCD 子程序 | `EuclideanDomain.gcd` + `normalize` |
| 求导 | 单独导数子程序 | `derivative` |

**建议**：
- L1 保留稀疏表示 + 单调序性质，用 `List (Nat × Nat)` 而非 `Array`（Lean 更自然）
- `toPoly` 函数桥接到 Mathlib，证明不变量守持

---

## 第 2 部分：多变量多项式

### 2.1 Mathlib `MvPolynomial σ R` 概览

**定义**（`Mathlib/Algebra/MvPolynomial/Basic.lean`）:
```lean
def MvPolynomial (σ : Type*) (R : Type*) [CommSemiring R] :=
  AddMonoidAlgebra R (σ →₀ ℕ)
```

其中 `σ →₀ ℕ` 是 monomial 索引（变量指数向量，有限支持）。

**与 CLPoly 的数据结构差异**：

| 维度 | CLPoly `polynomial_<T, order>` | Mathlib MvPolynomial σ R |
|-----|--------|---------|
| 单项式索引 | C++ `monomial_t` struct + order 信息 | `σ →₀ ℕ`（无 order） |
| 存储 | `Array<{monomial, coeff}>` 按 order 排序 | `Finsupp (σ →₀ ℕ) R`（支持集 + 函数） |
| order 消息 | 硬编码（grlex/grevlex）及算法决策 | **无对应**（Mathlib 无 term order） |
| 计算性 | 稀疏 Array 高效 | Finsupp 严格有限支持，数学性强 |

### 2.2 Mathlib 多变量 API

**关键函数**：
- `MvPolynomial.eval α`：代入求值（`σ → R` → `R`）
- `MvPolynomial.rename σ ρ`：变量重索引
- `MvPolynomial.degreeOf i`：关于第 i 个变量的次数
- `MvPolynomial.induction_on`：归纳法则
- `MvPolynomial.finSuccEquiv`：n+1 变量 ↔ Polynomial(n 变量多项式)

**L2 用法**（`CLPoly/Algorithm/Wang.lean`）:
```lean
open MvPolynomial
-- 用 finSuccEquiv 将多变量问题化为单变量多项式
def mvContentX0 (f : MvPolynomial (Fin (n+1)) ℤ) :=
  (finSuccEquiv ℤ n f).content
```

### 2.3 Finsupp 与 稀疏表示

**Finsupp** (`Mathlib/Data/Finsupp/Defs.lean`):
```lean
structure Finsupp (α : Type*) (M : Type*) [Zero M] where
  support : Finset α
  toFun : α → M
  mem_support_toFun : ∀ a, a ∈ support ↔ toFun a ≠ 0
```

**与 CLPoly Array 的对比**：
- Finsupp 保证有限支持，数学上清晰
- 但 Finsupp 不强制 order（Mathlib 无 term order 概念）
- CLPoly Array + order 信息是**计算优化**

### 2.4 L1 映射决策

**结论**：多变量多项式**需要自定义 shim**。

| 决策 | 原因 |
|------|------|
| ❌ 不直接用 MvPolynomial σ R | Mathlib 无 term order；稀疏 Array + order 是算法必需 |
| ✓ L1 自定义 struct | 1:1 对应 C++ 的 `polynomial_<T, order>`：Array sparse terms + order 标志 |
| ✓ L2 精化到 MvPolynomial | 证明 L1 稀疏和与 L2 多项式等值；order 分开处理 |

**L1 struct 建议**：
```lean
structure SparsePolyMv (σ : Type*) (R : Type*) (order : MonoOrder σ) where
  terms : List ({m : σ →₀ ℕ} × R)  -- 单项式 + 系数对，按 order 排序
  h_sorted : terms.Chain' (fun a b => order a.1 b.1)
  h_nonzero : ∀ t ∈ terms, t.2 ≠ 0

-- L2 精化
def SparsePolyMv.toMvPoly : MvPolynomial σ R :=
  terms.foldr (fun ⟨m, c⟩ acc => acc + C c * X^m) 0
```

---

## 第 3 部分：有限域与数域

### 3.1 整数 `Int` 与 `Rat`

**Int**：Lean core，直接可用。
- GCD：`Int.gcd`（Euclidean 实现）
- 因式分解：`UniqueFactorizationMonoid.normalizedFactors : Int → Multiset ℤ`
- L2 已用于多变量系数

**Rat**：Lean core，Mathlib 扩展丰富。
- 分子/分母：`.num`、`.den`
- 规范化：自动在构造时完成
- L2 已用（见 `Spec.lean` 的 `FactorQQ`）

**L1 映射**：
- `ZZ` struct（GMP 指针） → **L1 保留 struct 包装，L2 精化到 `Int`**
- `QQ` 直接可用 → **L2 直接用 `Rat`，L1 对应 C++ Rational 结构**

### 3.2 ZMod p：编译时素数 vs 运行时素数

**Mathlib ZMod p**：
```lean
def ZMod (n : ℕ) : Type  -- n 在编译时固定
-- p 是素数时的实例
instance [Fact (Nat.Prime p)] : Field (ZMod p)
```

**特点**：
- ✓ 素数时是域（自动获得 Field 实例）
- ✓ Ring homomorphism：`ZMod.cast : ℤ → ZMod p`
- ✓ 可计算的代表：`ZMod.val a : ℕ` 返回 `[0, p)` 中的代表
- ❌ 素数在编译时，不可动态改变
- ❌ 不能在运行时分支 `if p.Prime then ... else ...`

**CLPoly Zp 的特点**：
- 值 + 素数对（`{value: uint64, prime: uint64}`）
- 运行时多项式系数来自动态的 p
- 需要动态检查 p 是否素数

### 3.3 L1 映射决策

**不能直接用 ZMod p**，需要自定义 shim。

**L1 选项**：

| 选项 | 优缺点 |
|------|--------|
| 自定义 `Zp` struct | ✓ 1:1 对应 C++；✓ 运行时 p；❌ Lean 中 p.Prime 是 Prop，不自动可计算 |
| 用 Nat 表示元素（隐含 p） | ✓ 简单；❌ 失去 Field 实例，所有运算手写 |
| 条件依赖 ZMod | ❌ L1 代码复杂度高 |

**推荐方案**：

L1 层：
```lean
structure Zp where
  val : Nat
  prime : Nat
  h_prime : Nat.Prime prime
  h_range : val < prime

-- 运算内联（对应 C++ Zp 类）
def Zp.add (a b : Zp) : Zp where
  val := if a.val + b.val < a.prime then a.val + b.val else a.val + b.val - a.prime
  prime := a.prime
  ...
```

L2 精化：
```lean
-- 给定 h : Nat.Prime p，将 L1 Zp 转换到 ZMod p
def Zp.toZMod (a : Zp) (hp : Nat.Prime p) (ha_eq_p : a.prime = p) :
    ZMod p := ⟨a.val, by omega⟩

-- 证明 L1 运算与 ZMod 运算一致
theorem Zp.add_spec (a b : Zp) (hp : Nat.Prime a.prime) :
    (a.add b).toZMod = a.toZMod + b.toZMod := by ...
```

---

## 第 4 部分：稀疏表示策略

### 4.1 List vs Array vs Finsupp

| 工具 | 定义 | 计算性 | 有序性 | L1 适配 |
|------|------|--------|--------|---------|
| `List` | 递归数据 | ✓ computable | ❌ 需手动维护 | ✓ 自然，pattern match |
| `Array` | 可变向量（实现细节） | ❌ noncomputable | ✓ 内在有序 | ✓ 对应 C++ std::vector |
| `Finsupp` | support + 函数 | ❌ noncomputable | ❌ 无序（纯数学） | ❌ 不适合 L1 |

### 4.2 CLPoly 稀疏多项式的特征

- **Array 存储**：高效迭代、缓存友好
- **单调序**：算法决策依赖（grevlex → 最高次幂优先）
- **运行时 order**：动态选择 term order（grlex/grevlex 等）

### 4.3 L1/L2 表示策略

**L1 建议**：用 `List (Nat × Nat)` 而非 `Array`
- ✓ Lean 中更易形式化（递归定义、pattern match）
- ✓ 定义精确的 sorted invariant：`List.Chain' (fun a b => a.1 > b.1) terms`
- ✓ 避免数组越界证明的复杂性
- ❌ 性能略低，但 L1 不关心效率

**L2 精化**：
- 证明 L1 List 的不变量对应 Mathlib Finsupp
- 若需要 order，单独定义 `TermOrder σ` trait，不在多项式结构中

---

## 第 5 部分：因式分解结果

### 5.1 Mathlib 因式分解 API

**normalizedFactors**（`Mathlib/RingTheory/.../NormalizedFactors.lean`）:
```lean
noncomputable def normalizedFactors (a : α) : Multiset α
  -- UFD α 中 a 的唯一不可约因子分解
  -- 返回 Multiset（支持重数）
```

**特点**：
- ✓ Multiset 自动计重数（无需 `(Poly, ℕ)` 对）
- ✓ 只需证明 UFD，自动获得存在性
- ❌ 非构造性（无算法）

### 5.2 CLPoly factorization 结构

```cpp
struct factorization<Poly> {
  content : Poly,  // 常数倍数（单位）
  factors : Array<(Poly, uint64)>  // (因子, 重数) 对
};
```

### 5.3 L1/L2 表示

**L1**：保留 C++ 结构
```lean
structure Factorization (Poly : Type*) where
  content : Poly
  factors : List (Poly × ℕ)
  h_no_unit : ∀ f ∈ factors.map Prod.fst, ¬IsUnit f
```

**L2**：精化到 Multiset
```lean
def Factorization.toMultiset (f : Factorization (Polynomial ℤ)) :
    Multiset (Polynomial ℤ) :=
  (f.factors.map (fun ⟨p, k⟩ => replicate k p)).join
```

---

## 第 6 部分：L2 现用调查

### 6.1 L2 多项式类型统计

基于 `CLPoly/Algorithm/` 和 `CLPoly/Pipeline/` 的 grep 分析：

| 类型 | 出现次数 | 主要文件 | 用途 |
|-----|---------|--------|------|
| `Polynomial (ZMod p)` | 120+ | SquarefreeZp, DDF, EDF, Hensel | 单变量因式分解 |
| `MvPolynomial σ ℤ` | 80+ | Wang, MultivarFactor | 多变量系数/求值 |
| `Polynomial (Polynomial ℤ)` | 30+ | Hensel, Recombine | 多变量系数嵌套 |
| `ZMod p` | 50+ | 各算法模块 | 有限域计算 |
| `List (Polynomial ...)` | 70+ | Spec, 各算法 | 因子列表 |

### 6.2 关键 API 依赖

**GCD 相关**：
- `EuclideanDomain.gcd`：30+ 使用点
- `normalize`：15+ 使用点（与 gcd 协作）
- `Monic`：25+ 使用点（gcd 结果正规化）

**多项式运算**：
- `divByMonic`, `modByMonic`：20+ 使用点
- `derivative`：10+ 使用点
- `natDegree`, `leadingCoeff`：50+ 使用点

**多变量**：
- `MvPolynomial.eval`：30+ 使用点
- `MvPolynomial.finSuccEquiv`：8 使用点
- `MvPolynomial.induction_on`：5 使用点

### 6.3 L2 不需要的 Mathlib

| API | 原因 |
|-----|------|
| `Polynomial.gcd` | 不存在（后续作为独立模块，不用 Mathlib 定义） |
| `Polynomial.factors` | L2 用自定义算法，不用 normalizedFactors |
| Hensel 提升 | L2 自定义算法，Mathlib 无标准库支持 |
| MDP 求解器 | L2 自定义（稀疏插值/BDP/WMDS），Mathlib 无 |

---

## 第 7 部分：最终映射决策表

### 7.1 CLPoly C++ 类型 → L1/L2

| C++ 类型 | L1 类型 | L2 类型 | 备注 |
|---------|--------|--------|------|
| `polynomial_<ZZ, grlex>` | `struct SparsePolyZZ { terms: List, prime_func }` | `MvPolynomial σ ℤ` | 自定义 shim；L1 保留稀疏 + order |
| `polynomial_<Zp, grevlex>` | `struct SparsePolyZp { terms: List, prime: Nat }` | `MvPolynomial σ (ZMod p)` | 自定义；order 分离 |
| `upolynomial_<ZZ>` | `struct SparsePolyUZZ { terms: List (Nat × ZZ) }` | `Polynomial ℤ` | L1 保留稀疏；L2 精化 |
| `upolynomial_<Zp>` | `struct SparsePolyZp { terms: List (Nat × Nat), prime: Nat }` | `Polynomial (ZMod p)` | **现有 L1 Types.lean** |
| `ZZ` | `struct ZZ { val: Nat, ... }` | `Int` | L1 保留 GMP 模型；L2 精化 |
| `QQ` | `struct Rat { num, den: ZZ, ... }` | `Rat` | 直接可用 |
| `Zp` | `struct Zp { coeff: Nat, prime: Nat, h_prime }` | `ZMod p` 条件投影 | 自定义；运行时 prime |
| `variable` | `Fin n` 或索引值 | `Fin n` | 直接 |
| `basic_monomial<grlex>` | `σ →₀ ℕ` with order info | `σ →₀ ℕ` | Finsupp 无 order；order 分离 |
| `factorization<Poly>` | `struct { content, factors: List }` | `Multiset` 或 `List` | L1 保留结构；L2 精化 |

### 7.2 需要自定义 shim 的列表

| CLPoly 类型 | Mathlib 缺口 | shim 位置 |
|------------|-----------|---------|
| `polynomial_<T, order>` 多变量 | 无 term order 概念 | L1 Impl/Types.lean |
| `Zp` 运行时素数 | ZMod p 素数在编译时 | L1 Impl/Types.lean |
| `basic_monomial<order>` | Finsupp 无 order | L1 Impl/Types.lean |
| `ZZ` GMP 包装 | Int 无 GMP 语义 | L1 Impl/ZZImpl.lean |

### 7.3 可直接复用 Mathlib

| Mathlib 类型/函数 | 用途 | L2 文件 |
|-----------------|------|--------|
| `Polynomial (ZMod p)` | 单变量有限域多项式 | SquarefreeZp, DDF, EDF, ... |
| `MvPolynomial σ R` | 多变量多项式（无 order） | Wang, MultivarFactor |
| `EuclideanDomain.gcd` | GCD 算法 | SquarefreeZp yunLoop |
| `divByMonic`, `modByMonic` | 精确多项式除法 | SQF, 各算法 |
| `normalize` | 多项式正规化（monic） | 与 gcd 协作 |
| `ZMod p` Field 实例 | 有限域性质 | Math/FiniteFieldFact |
| `Int`, `Rat` | 整数/有理数 | 核心数据类型 |
| `Multiset` | 重数集合 | Spec 定义因式分解 |
| `MvPolynomial.eval` | 多变量求值 | Wang 算法 |
| `MvPolynomial.finSuccEquiv` | 多变量层化 | Wang content 计算 |

---

## 第 8 部分：特殊问题解答

### 8.1 ZMod p 的 Computability

**问题**：ZMod p 在 Mathlib 中是 noncomputable（Fintype 基础），如何在 L1 中处理可计算性？

**答**：
- L1 不需要 Lean 可计算性（L1 是形式化模型，对应 C++ 可执行代码）
- L1 中 `Zp` struct 的运算（`add`, `mul` 等）定义为 noncomputable，因为它们操作的是数学对象
- 可计算性在 L0（C++ 执行层）保证；L1 只需形式验证正确性
- L2 可能需要 noncomputable 实例（Field 实例需要 ZMod p）

### 8.2 GCD 的 Monic 问题

**问题**：EuclideanDomain.gcd 不保证 monic，每次 GCD 都要 normalize，很麻烦。

**答**：这是 Mathlib 的设计（一般环不有 monic 概念），不可回避。
- L2 中 Polynomial (ZMod p) 有 Monic 定义
- 每次 GCD 后加一行 `let y := normalize (gcd w c)` + `have hy_monic := Polynomial.monic_normalize ...`
- SQF 已验证这个模式可行（0 sorry）

### 8.3 Order 的处理

**问题**：CLPoly 的 term order（grlex/grevlex）在 Mathlib 中无对应。

**答**：
- **L1 层**：order 作为 `TermOrder σ` class 或 `enum` 参数，编码在稀疏多项式结构中
- **L2 层**：抽象掉 order，只用 `MvPolynomial σ R`（数学对象）；算法正确性不依赖具体 order（只要一致）
- **Spec**：多变量因式分解 spec 不提及 order，只保证最终分解的正确性

### 8.4 Multiset vs List for 因式分解

**问题**：Mathlib 用 Multiset，CLPoly 用 Array/List。两者区别？

**答**：
- **Multiset**：无序，允许重复；数学上代表"集合带重数"；`normalizedFactors` 返回 Multiset
- **List**：有序，允许重复；对应 C++ Array
- **L1**：保留 `List (Poly × ℕ)` 对（对应 C++ factorization.factors）
- **L2**：精化关系 `factors_list_to_multiset : List.prod ≈ Multiset.prod`（顺序不影响乘积）

---

## 第 9 部分：实施建议

### 9.1 L1 Types 文件结构

```lean
-- CLPoly/Impl/Types.lean

-- 1. Zp 自定义（对应 C++ Zp class）
structure Zp where ...
def Zp.add, Zp.mul, ... : Zp → Zp → Zp

-- 2. 稀疏单变量多项式（对应 upolynomial_<Zp>）
structure SparsePolyZp where ...
  h_sorted : List.Chain' (fun a b => a.1 > b.1) terms

-- 3. 稀疏多变量多项式（对应 polynomial_<Zp, order>）
structure SparsePolyMv (σ : Type*) where ...
  h_sorted : List.Chain' (fun a b => order a.1 b.1) terms

-- 4. L2 精化桥接
noncomputable def SparsePolyZp.toPoly : Polynomial (ZMod p)
noncomputable def SparsePolyMv.toMvPoly : MvPolynomial σ (ZMod p)

-- 5. 精化定理（关键不变量）
theorem SparsePolyZp.toPoly_add : (a.add b).toPoly = a.toPoly + b.toPoly
...
```

### 9.2 L1 验证策略

**阶段 1**（types + 基本运算）：
- 定义 Zp, SparsePolyZp, SparsePolyMv 结构
- 证明运算的基本性质（封闭性、结合律等）
- 0 sorry

**阶段 2**（L2 精化）：
- 定义 `toPoly` / `toMvPoly` 函数
- 证明 L1 运算与 L2 运算对应（`add_spec`, `mul_spec` 等）
- 依赖 L2 的 Polynomial/MvPolynomial 定理
- 0 sorry

**阶段 3**（算法正确性）：
- 对各个算法（SQF, DDF, EDF, Hensel, Wang）建立 L1 模型
- 证明 L1 模型与 L2 算法一致
- 继承 L2 的 0-sorry 证明

### 9.3 注意事项

1. **避免过早精化**：L1 应直接对应 C++ 控制流，不应引入 Mathlib 复杂性
2. **分离关切**：order、GCD、多项式除法 分别定义，不混杂在主体中
3. **不变量文档**：每个 struct 的 invariant（sorted, nonzero, range）都要注释清楚
4. **测试精化**：小规模 L1 to L2 的精化定理先单独验证，再组装

---

## 参考资料

- `Mathlib/Algebra/Polynomial/Basic.lean`：Polynomial 定义 + API
- `Mathlib/Algebra/MvPolynomial/Basic.lean`：MvPolynomial 定义
- `Mathlib/Data/Finsupp/Defs.lean`：Finsupp（有限支持函数）
- `Mathlib/Data/ZMod/Defs.lean`：ZMod p 定义
- `Mathlib/RingTheory/UniqueFactorizationDomain/NormalizedFactors.lean`：因式分解
- `Mathlib/Algebra/Polynomial/Div.lean`：divByMonic / modByMonic
- `CLPoly/CLAUDE.md`：L2 开发经验（especially §关键 API 陷阱）
- `CLPoly/Impl/Types.lean`：已有 L1 Zp + SparsePolyZp 模型

