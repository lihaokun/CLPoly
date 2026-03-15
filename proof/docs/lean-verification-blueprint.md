# 因式分解模块 Lean 4 形式化验证蓝图

## 1. 目标

对 CLPoly 的 Zp[x] → Z[x] 因式分解管线进行**机器检查的形式化验证**，覆盖：

1. **数学正确性**：算法输出确实是输入的不可约分解
2. **算法不变量**：循环不变量、终止性、中间状态的数学性质
3. **表示层安全**：整数溢出、截断、模运算正确性
4. **内存操作安全**：数组越界、别名冲突、move-after-use

验证架构遵循 `docs/pre-verification-architecture.md` 的分支二（完整数学对象）路径。

## 2. 整体架构

```
┌─────────────────────────────────────────────────┐
│               Lean 4 定理证明器                    │
│                                                   │
│  L3: 数学基石                                      │
│      有限域理论、多项式环、不可约性、Hensel 引理       │
│      ↕ 正确性证明                                  │
│  L2: 算法模型                                      │
│      DDF/EDF/SQF/Hensel/LLL 的算法逻辑（Mathlib 类型）│
│      ↕ 精化证明                                    │
│  L1: 实现模型                                      │
│      1:1 对应 C++（uint64 语义、数组越界、move）      │
│                                                   │
└───────────────────┬─────────────────────────────┘
                    │ 人工审查：Lean 模型 ↔ C++ 代码对应
                    ▼
              C++ 实现（现有代码）
```

**信任假设 (TCB)**：
- Lean 4 核心 + Mathlib 正确
- L1 实现模型忠实反映 C++ 语义（人工审查保证，可通过 crosscheck 增强信心）
- 编译器 + 硬件正确

## 3. L3：数学基础层

### 3.1 Mathlib 已有（可直接使用）

| 概念 | Mathlib 位置 | 用于 |
|------|------------|------|
| `Polynomial R` | `Mathlib.RingTheory.Polynomial.Basic` | 多项式环 |
| `ZMod p` | `Mathlib.Data.ZMod.Basic` | 有限域 Zp |
| `Irreducible` | `Mathlib.Algebra.Prime.Basic` | 不可约性定义 |
| `IsCoprime` | `Mathlib.Algebra.IsCoprime` | 互素 |
| `Polynomial.derivative` | `Mathlib.RingTheory.Polynomial.Basic` | 导数 |
| `Squarefree` | `Mathlib.Algebra.Squarefree.Basic` | 无平方 |
| `GCDMonoid` | `Mathlib.Algebra.GCDMonoid.Basic` | GCD 存在性 |
| `Fintype (ZMod p)` | `Mathlib.Data.ZMod.Basic` | |ZMod p| = p |

### 3.2 需要新建/扩展

| 概念 | 数学内容 | 难度 | 用于 |
|------|---------|------|------|
| **定理 2.1** | $x^{p^d} - x = \prod\{g : \text{irred}, \deg g \mid d\}$ | 中 | DDF 正确性核心 |
| **推论 2.2** | $\gcd(x^{p^d}-x, f) = \prod\{f_i : \deg f_i \mid d\}$ | 低（2.1 直推） | DDF |
| **Hensel 引理** | p-adic 提升：$f \equiv g \cdot h \pmod{p} \Rightarrow$ 提升到 $\pmod{p^k}$ | 高 | Hensel 提升 |
| **Mignotte 界** | $\|g\|_\infty \leq \binom{n}{k} \|f\|_2$ | 中 | 精度估算 |
| **Cantor-Zassenhaus 概率** | 分裂概率 $\geq 1 - 2^{1-k}$ | 中 | EDF 终止性 |
| **Fermat 小定理（多项式版）** | $a^{p^d-1} = 1$ in $\mathbb{F}_{p^d}^*$ | 低（Mathlib 有标量版） | EDF |

### 3.3 定理 2.1 的证明策略

这是 DDF 正确性的数学基石。证明路线：

```
Mathlib: ZMod p 是域
  → Polynomial (ZMod p) 是 Euclidean domain
  → roots of x^{p^d} - x = 全部 F_{p^d} 元素（Mathlib 有 Frobenius）
  → 不可约多项式 g, deg g | d ↔ g 的根全在 F_{p^d} 中
  → g | (x^{p^d} - x)
  → 无平方 + UFD → 精确分解
```

## 4. L2：算法模型

### 4.1 需要建模的函数

建模 `polynomial_factorize_zp.hh` 的算法逻辑，操作 Mathlib 数学类型（`Polynomial (ZMod p)` 等），抽象掉 C++ 实现细节（整数表示、数组布局、内存管理）：

#### 4.1.1 `__squarefree_Zp`（无平方分解）

```lean
/-- Zp[x] 无平方分解：返回 (因子, 重数) 列表 -/
def squarefree_Zp (f : ZpPoly) : List (ZpPoly × Nat) := ...
```

**证明义务**：
- `∀ (g, e) ∈ result, Squarefree g`
- `∀ (g, e) ∈ result, e ≥ 1`
- `f = unit * ∏ (g, e) ∈ result, g ^ e`
- `∀ i ≠ j, IsCoprime result[i].1 result[j].1`

#### 4.1.2 `__ddf_Zp`（按度分组）

```lean
/-- DDF：输入首一无平方多项式，返回 (度d的不可约因子之积, d) 列表 -/
def ddf_Zp (f : ZpPoly) : List (ZpPoly × Nat) := ...
```

**证明义务**：
- 循环不变量：`h ≡ x^{p^d} (mod f_star)`（引理 3.1）
- 分裂正确性：`gd = ∏{f_i : irred, f_i | f_star, deg f_i = d}`（定理 3.2）
- 提前终止：`deg f_star < 2d → f_star 至多一个不可约因子`（命题 3.3）
- 终止性：`deg f_star` 严格递减或 `d` 递增至 `deg f_star < 2d`

#### 4.1.3 `__edf_Zp`（等度分裂，Cantor-Zassenhaus）

```lean
/-- EDF：输入所有不可约因子度为 d 的首一无平方多项式，返回全部不可约因子 -/
partial def edf_Zp (f : ZpPoly) (d : Nat) (rng : RngState) : List ZpPoly := ...
```

**证明义务**：
- CRT 分解：`Zp[x]/(f) ≅ ∏ Zp[x]/(f_i)`（引理 4.1）
- 分裂概率：`P[非平凡分裂] ≥ 1 - 2^{1-k}`（定理 4.4）
- `partial` 标注：概率终止，附加定理证明期望迭代次数有界

#### 4.1.4 `__upoly_powmod`（模幂）

```lean
/-- 二进制模幂：base^exp mod modpoly -/
def upoly_powmod (base : ZpPoly) (exp : Nat) (modpoly : ZpPoly) : ZpPoly := ...
```

**证明义务**：
- `upoly_powmod base exp modpoly ≡ base ^ exp (mod modpoly)`
- 终止性：`exp` 每次右移，well-founded on `Nat`

#### 4.1.5 `__upoly_subtract_x` / `__upoly_subtract_one`（辅助函数）

```lean
def upoly_subtract_x (h : ZpPoly) : ZpPoly := ...
def upoly_subtract_one (h : ZpPoly) : ZpPoly := ...
```

**证明义务**：
- `upoly_subtract_x h = h - X`
- `upoly_subtract_one h = h - 1`
- **B3 对应**：建模中显式要求 `-1 ≡ p-1 (mod p)` 的表示正确性

#### 4.1.6 `__factor_Zp`（编排）

```lean
/-- Zp[x] 完整不可约分解 -/
def factor_Zp (f : ZpPoly) : Zp × List (ZpPoly × Nat) := ...
```

**证明义务（顶层定理）**：
```lean
theorem factor_Zp_correct (f : ZpPoly) (hf : f ≠ 0) :
    let (lc, factors) := factor_Zp f
    -- 1. 乘积还原
    lc • ∏ (g, e) ∈ factors, g ^ e = f
    -- 2. 每个因子不可约
    ∧ ∀ (g, e) ∈ factors, Irreducible g ∧ Monic g ∧ e ≥ 1
    -- 3. 因子两两互素
    ∧ ∀ i j, i ≠ j → IsCoprime factors[i].1 factors[j].1 := by
  sorry
```

### 4.2 第二阶段：Z[x] 因式分解

按 `polynomial_factorize_univar.hh` 建模（在 Zp 层完成后）：

| 函数 | 证明义务核心 |
|------|------------|
| `__hensel_lift` | Hensel 引理：提升后因子模 $p^k$ 仍是分解 |
| `__heuristic_starting_precision` | Mignotte 界的正确性 |
| `__zassenhaus_recombine` | 子集乘积 = 真因子（mod → Z 的提升） |
| `__lll_factorize` | Phase 1 + Phase 2 的完整性：不丢因子 |

## 5. L1：实现模型

### 5.1 uint64_t 模型

```lean
/-- 64-bit 无符号整数，模 2^64 算术 -/
abbrev U64 := Fin (2^64)

/-- 64-bit 有符号整数，补码表示 -/
def I64 := { n : Int // -2^63 ≤ n ∧ n < 2^63 }

/-- C++ 的 (int64_t)(uint64_t v) 强制转换 -/
def cast_u64_to_i64 (v : U64) : I64 := ...
-- 关键性质：v ≥ 2^63 时发生回绕

/-- Zp(int64_t val, uint64_t p) 构造函数语义 -/
def zp_from_i64 (val : I64) (p : U64) : U64 := ...

/-- Zp(uint64_t val, uint64_t p) 构造函数语义 -/
def zp_from_u64 (val : U64) (p : U64) : U64 := val % p
```

**B3 重现与证明**：
```lean
-- 旧代码（有 bug）的模型：
def old_neg_one (p : U64) : U64 :=
  zp_from_i64 (cast_u64_to_i64 (p - 1)) p

-- 新代码（修复后）的模型：
def new_neg_one (p : U64) : U64 :=
  zp_from_u64 (p - 1) p

-- 证明旧代码对 p > 2^63 错误：
theorem old_neg_one_wrong (p : U64) (hp : p.val > 2^63) :
    old_neg_one p ≠ p - 1 := by ...

-- 证明新代码对任意 p 正确：
theorem new_neg_one_correct (p : U64) (hp : p.val ≥ 2) :
    new_neg_one p = p - 1 := by ...
```

### 5.2 Vector 模型

```lean
/-- C++ std::vector<T> 的 Lean 模型 -/
structure Vec (α : Type) where
  data : Array α
  -- 不变量：data.size 是实际元素数

/-- 带越界检查的索引访问 -/
def Vec.get (v : Vec α) (i : Nat) (h : i < v.data.size) : α :=
  v.data.get ⟨i, h⟩

/-- push_back：返回新 Vec，旧引用不再有效 -/
def Vec.push_back (v : Vec α) (x : α) : Vec α :=
  ⟨v.data.push x⟩
```

### 5.3 Move 语义追踪

```lean
/-- 资源状态：可用 或 已移走 -/
inductive Ownership (α : Type)
  | alive : α → Ownership α
  | moved : Ownership α

/-- 使用资源：需证明未被 move -/
def use {α : Type} : Ownership α → (h : ¬ is_moved o) → α
  | .alive v, _ => v

/-- 移走资源：消耗并标记 -/
def move_from {α : Type} : Ownership α → (h : ¬ is_moved o) → α × Ownership α
  | .alive v, _ => (v, .moved)
```

在算法模型中，每个 `std::move` 对应一次 `move_from`，后续访问需提供 `¬ is_moved` 证明。若代码正确（move 后不再使用），证明自动成立。

### 5.4 迭代器失效模型

```lean
/-- Vector 引用：绑定到特定版本 -/
structure VecRef (α : Type) where
  vec_version : Nat    -- 创建时的 vector 版本号
  index : Nat

/-- push_back 递增版本号 -/
def Vec.push_back_versioned (v : Vec α) (x : α) : Vec α :=
  { data := v.data.push x, version := v.version + 1 }

/-- 通过引用访问：需版本号匹配 -/
def Vec.deref (v : Vec α) (ref : VecRef α)
    (h_ver : ref.vec_version = v.version)  -- 版本匹配 = 未失效
    (h_idx : ref.index < v.data.size) : α := ...
```

## 6. 实施计划

实施计划的权威版本在 `proof/docs/implementation-roadmap.md`（v2 自顶向下）。以下仅列出 Phase 划分概览，详细任务清单和验收标准见路线图。

| Phase | 内容 | 对应层 |
|-------|------|--------|
| Phase 0 | 环境搭建 + Mathlib API 验证 | — |
| Phase 1 | 顶层正确性骨架 + 接口规约锁定 | Spec + Pipeline |
| Phase 2 | L3 数学基石（Thm 2.1、Cor 2.2、CRT、概率分析） | L3 Math |
| Phase 3 | L2 算法模型 — Zp[x] 管线（SQF/DDF/EDF） | L2 Algorithm |
| Phase 4 | L2 算法模型 — Z[x] 管线（Hensel/重组/LLL） | L2 Algorithm + L3 Math |
| Phase 5 | L1 实现模型（U64/Vec/Move，1:1 对应 C++，精化 L2） | L1 Impl |

**Phase 0 已完成**：E1-E4 实验通过，Mathlib Gap Report 完成。

## 7. Path C 执行流程

### 7.1 单函数验证的完整工作流

以 `__ddf_Zp` 为例，每个函数走完以下 6 步：

```
Step 1: 提取           C++ 源码 → 标注版 C++（标记控制流 + 变量生命期）
Step 2: L2 建模        提取算法逻辑 → Lean 4 函数定义（Mathlib 类型，抽象实现细节）
Step 3: 陈述           写出正确性定理（前条件 + 后条件 + 不变量）
Step 4: 证明（L2）     证明算法正确性（利用 L3 数学定理）
Step 5: L1 建模+证明   1:1 对应 C++ 控制流（U64/Vec/Ownership），精化证明 L1 行为 = L2
Step 6: 审查           人工核对 L1 模型 ↔ C++ 代码的逐行对应
```

### 7.2 Step 1：提取（C++ → 标注版 C++）

在原始 C++ 代码上标注，不改动代码本身。标注内容为 L1（1:1 对应）所需的信息：

```cpp
// [LEAN: loop_var h, f_star, d]
// [LEAN: invariant h ≡ x^{p^d} mod f_star]
// [LEAN: decreasing deg(f_star) or terminates when deg(f_star) < 2d]
for (uint64_t d = 1; ; ++d)
{
    if (get_deg(f_star) < (int64_t)(2 * d))  // [LEAN: early_exit]
        break;

    h = __upoly_powmod(h, ZZ(p), f_star);    // [LEAN: call powmod]
    auto h_minus_x = __upoly_subtract_x(h, p); // [LEAN: call subtract_x]
    auto gd = polynomial_GCD(h_minus_x, f_star); // [LEAN: call gcd]

    if (!gd.empty() && get_deg(gd) > 0)       // [LEAN: branch nontrivial_gd]
    {
        // [LEAN: f_star_new = f_star / gd, deg decreases]
        // [LEAN: move f_new → f_star, f_new is consumed]
        ...
    }
}
```

标注的目的：
- 明确哪些是循环变量、不变量、递减量
- 标记每个 `std::move` 的消耗点
- 标记每个数组访问的越界条件
- 为 Step 2 翻译提供精确映射

### 7.3 Step 2：L2 建模（C++ 算法逻辑 → Lean 4 算法模型）

L2 建模提取 C++ 的算法逻辑，操作 Mathlib 数学类型，抽象掉实现细节：

| C++ 构造 | L2 算法模型 | L1 实现模型（Phase 5） |
|---------|------------|---------------------|
| `for` 循环 | 数学归纳 / 递归 | 顶层递归 + `termination_by`（1:1 控制流） |
| `while(true)` | 存在性论证 | `partial def` 或 fuel 参数 |
| `polynomial_GCD(a, b)` | `GCDMonoid.gcd a b` | 调用 Euclid GCD 的 L1 模型 |
| `(int64_t)(p - 1)` | `p - 1`（自然数） | `cast_u64_to_i64 (p - 1)`（显式溢出语义） |
| `vec[i]` | `list.get i` | `Vec.get i (proof : i < size)` |
| `std::move(v)` | 无对应（纯函数式） | `Ownership.move_from v (proof : ¬moved)` |

**L2 建模原则**：Lean 函数捕获 C++ 的算法逻辑（DDF 循环、EDF 随机分裂、Hensel 提升），但操作 Mathlib 的 `Polynomial (ZMod p)` 等数学类型，抽象掉整数表示、数组布局、内存管理。

DDF 的 L2 算法模型：

```lean
/-- __ddf_Zp 的算法模型（L2）：操作 Mathlib Polynomial (ZMod p) -/
def ddfLoop (h f_star : Polynomial (ZMod p)) (d : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    : List (Polynomial (ZMod p) × ℕ) :=
  if f_star.natDeg < 2 * d then
    if 0 < f_star.natDeg then acc ++ [(f_star, f_star.natDeg)] else acc
  else
    let h' := (h ^ p) %ₘ f_star           -- powmod
    let gd := gcd (h' - X) f_star          -- gcd 步
    if 0 < gd.natDeg then
      ddfLoop (h' %ₘ (f_star /ₘ gd)) (f_star /ₘ gd) (d + 1) (acc ++ [(gd, d)])
    else
      ddfLoop h' f_star (d + 1) acc
termination_by f_star.natDeg + 1 - 2 * d
```

### 7.4 Step 3：陈述（写正确性定理）

对每个函数，陈述两类定理（L2 和 L1 分别陈述）：

**A. 功能正确性**（L2 算法模型）
```lean
theorem ddf_Zp_correct (f : ZpPoly) (hm : Monic f) (hsq : Squarefree f) :
    let result := ddf_Zp f p
    -- 每个 (gd, d) 中 gd 恰为 f 的全部 d 次不可约因子之积
    ∀ (gd, d) ∈ result,
      gd ∣ f
      ∧ (∀ q, Irreducible q → q ∣ gd → q.natDeg = d)
      ∧ (∀ q, Irreducible q → q ∣ f → q.natDeg = d → q ∣ gd)
    -- 不遗漏：f 的所有不可约因子都出现在某个 gd 中
    ∧ (∀ q, Irreducible q → q ∣ f →
         ∃ (gd, d) ∈ result, q ∣ gd ∧ q.natDeg = d) := by
  sorry
```

**B. 循环不变量**（L2 算法模型，辅助引理）
```lean
theorem ddf_loop_invariant (h f_star : ZpPoly) (d : Nat) :
    -- 在第 d 次迭代开始时（powmod 之后）
    h ≡ X ^ (p ^ d) [MOD f_star]
    -- f_star 的所有不可约因子度 ≥ d
    ∧ (∀ q, Irreducible q → q ∣ f_star → q.natDeg ≥ d) := by
  sorry
```

**C. 实现模型安全**（L1 实现模型，Phase 5）
```lean
/-- L1 精化：实现模型行为 = L2 算法模型 -/
theorem ddf_Zp_impl_refines (f : ZpPoly) (p : U64) (hp : p.val ≥ 2) :
    -- L1 实现模型的输出与 L2 算法模型一致
    ddf_Zp_impl f p = ddf_Zp_algo f p
    -- 且所有 U64 运算无溢出、数组访问在界内
    ∧ no_overflow (ddf_Zp_impl f p)
    ∧ all_array_accesses_in_bounds (ddf_Zp_impl f p) := by
  sorry
```

### 7.5 Step 4-5：证明

证明顺序：**自底向上**。

```
L3 数学定理（定理 2.1）
    ↓ 作为引理使用
L2 辅助函数正确性（powmod, subtract_x）
    ↓ 组合
L2 循环不变量（ddf_loop_invariant）
    ↓ 归纳
L2 功能正确性（ddf_Zp_correct）
    ↓ 并行
L1 表示层安全（ddf_Zp_no_overflow）
```

L2 证明策略：
- 循环不变量用**归纳法**：基础 `d=1` + 归纳步 `d → d+1`
- 归纳步的核心是定理 3.2（`gcd(h-x, f*) = 全部 d 次因子`）
- 终止性用 `f_star.natDeg` 的严格递减（提取非平凡 gd 时）或 `2d` 超过 `deg f_star`

L1 精化证明策略（Phase 5）：
- 精化关系：证明 L1（U64/Vec 操作）的输出与 L2（Mathlib 类型操作）一致
- `p - 1` 不溢出：`p ≥ 2 → p - 1 ≥ 1`，U64 减法安全
- 数组越界：追踪 `_coeffs.size()` = `deg + 1`，每次访问的 index ≤ deg

### 7.6 Step 6：人工审查

审查分两层进行：

**L2 审查**（Phase 3-4）：确认算法模型正确捕获了 C++ 的算法逻辑

| # | 检查项 | 方法 |
|---|--------|------|
| 1 | Lean 算法模型的**逻辑结构**与 C++ 算法一致 | 对照算法分支/循环/递归结构 |
| 2 | 数学操作（gcd、divByMonic 等）正确对应 C++ 调用 | 调用映射表 |
| 3 | 边界条件（提前终止、空输入等）处理一致 | 逐条对比 |

**L1 审查**（Phase 5）：确认实现模型 1:1 对应 C++ 代码

| # | 检查项 | 方法 |
|---|--------|------|
| 1 | L1 Lean 函数与 C++ 函数的**控制流同构** | 逐行对照 |
| 2 | 每个 C++ 变量在 L1 中有**同名对应** | 变量映射表 |
| 3 | C++ 的隐式转换在 L1 中**显式建模** | 检查所有 int/uint 赋值 |
| 4 | `std::move` 在 L1 中有 `Ownership.move_from` 对应 | move 点列表比对 |
| 5 | 数组访问的越界证明的**前提条件**与 C++ 的**运行时保证**一致 | 逐个检查 |

审查产物：`CLPoly/Review/DDF_review.md`，记录每个检查项的结论。

### 7.7 持续同步机制

当 C++ 代码修改时：

```
1. CI 检查：diff 涉及已建模函数 → 标记 "lean-model-sync-needed"
2. 开发者更新 Lean 模型（同一 PR 中）
3. Lean build 通过（证明仍成立）→ CI green
4. 重新执行 Step 6 审查 → 更新 review 文档
```

可通过 `git diff --name-only` + 函数名匹配自动触发。

## 8. 文件结构

```
lean/
├── CLPoly.lean                    -- 主入口
├── CLPoly/
│   ├── Math/                      -- L3: 数学基础
│   │   ├── FiniteField.lean       -- 定理 2.1, 推论 2.2
│   │   ├── HenselLemma.lean       -- Hensel 引理
│   │   ├── MignotteBound.lean     -- Mignotte 界
│   │   └── CantorZassenhaus.lean  -- EDF 概率分析
│   ├── Algorithm/                 -- L2: 算法模型
│   │   ├── ZpPoly.lean            -- ZpPoly 类型 + 基本操作
│   │   ├── Powmod.lean            -- upoly_powmod
│   │   ├── SquarefreeZp.lean      -- squarefree_Zp
│   │   ├── DDF.lean               -- ddf_Zp + 正确性证明
│   │   ├── EDF.lean               -- edf_Zp + 分裂概率
│   │   ├── FactorZp.lean          -- factor_Zp 编排
│   │   ├── HenselLift.lean        -- Hensel 提升
│   │   └── LLLFactorize.lean      -- Z[x] 完整流程
│   └── Impl/                      -- L1: 实现模型（1:1 对应 C++）
│       ├── UInt64.lean            -- U64/I64 模型 + 转换语义
│       ├── Barrett.lean           -- Barrett 模乘正确性
│       ├── Vector.lean            -- Vec 模型 + 越界安全
│       └── Ownership.lean         -- Move 语义追踪
└── lakefile.lean                  -- 构建配置
```

## 8. 验证强度总结

| bug 类别 | 覆盖层 | 机制 | 已知实例 |
|---------|-------|------|---------|
| 算法逻辑错误 | L2 | 循环不变量 + 正确性定理 | B2 Phase 2 触发 |
| 数学定理误用 | L3 + L2 | 数学证明 + 精化 | — |
| 整数溢出/截断 | L1 | U64/I64 显式建模 | B3 `(int64_t)(p-1)` |
| 数组越界 | L1 | Vec.get 需 `i < size` 证明 | — |
| 除零 | L1 + L2 | 前提条件 `divisor ≠ 0` | nmod_inv |
| 指针别名 | L1 | 函数参数 aliasing 显式建模 | divrem aliasing guard |
| use-after-move | L1 | Ownership 状态追踪 | — |
| 迭代器失效 | L1 | Vec 版本号匹配 | — |
| 概率终止 | L2 | 期望迭代次数有界证明 | EDF while(true) |

**唯一不覆盖**：编译器行为（UB 利用、优化）、硬件错误、OOM。这些由 UBSan + 测试兜底。

## 9. 风险与缓解

| 风险 | 影响 | 缓解 |
|------|------|------|
| Mathlib 缺定理 2.1 | Phase 1 阻塞 | Phase 0 先评估；最坏情况自行证明（~200 行 Lean） |
| Hensel 引理形式化过难 | Phase 3 延期 | Phase 1-2 独立有价值；Hensel 可后置 |
| Lean 模型与 C++ 漂移 | 验证失去意义 | CI 中加 `lean-model-review` 检查点；C++ 改动同步更新 Lean |
| 证明工程量超预期 | 进度延迟 | Phase 1 是试探性的，完成后重新评估 |

## 10. 成功标准

**Phase 1 完成时**：Lean 4 中有一个经过 `#check` 通过的 `ddf_correct` 定理，从定理 2.1 出发，经过循环不变量，证明 DDF 输出的每个 `(gd, d)` 确实是 `f` 中所有 `d` 次不可约因子之积。

**全部完成时**：`factor_Zp_correct` 和 `factor_ZZ_correct` 两个顶层定理通过 Lean 4 机器检查，覆盖从数学定义到算法实现的完整正确性链条。
