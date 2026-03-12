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
│  L3: 数学基础层                                    │
│      有限域理论、多项式环、不可约性、Hensel 引理       │
│      ↕ 正确性证明                                  │
│  L2: 算法层                                       │
│      DDF/EDF/SQF/Hensel/LLL 的 1:1 函数模型        │
│      ↕ 精化证明                                    │
│  L1: 表示层                                       │
│      uint64_t 语义、vector 模型、move 追踪           │
│                                                   │
└───────────────────┬─────────────────────────────┘
                    │ 人工审查：Lean 模型 ↔ C++ 代码对应
                    ▼
              C++ 实现（现有代码）
```

**信任假设 (TCB)**：
- Lean 4 核心 + Mathlib 正确
- L1 模型忠实反映 C++ 语义（人工审查保证，可通过 crosscheck 增强信心）
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

## 4. L2：算法层

### 4.1 需要建模的函数

按 `polynomial_factorize_zp.hh` 的结构 1:1 建模：

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

## 5. L1：表示层

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

### Phase 0：环境搭建（1 周）

- [ ] 创建 `lean/` 子目录，初始化 Lean 4 + Mathlib 项目
- [ ] 确认 Mathlib 的 `Polynomial (ZMod p)` API 可用性
- [ ] 定义 `ZpPoly` 类型别名和基本操作封装

### Phase 1：DDF 概念验证（2-3 周）

目标：端到端验证 `__ddf_Zp` 的正确性。

- [ ] L3: 形式化定理 2.1（或在 Mathlib 中找到等价结论）
- [ ] L3: 形式化推论 2.2
- [ ] L2: 建模 `upoly_powmod`，证明 `result ≡ base^exp (mod m)`
- [ ] L2: 建模 `upoly_subtract_x`，证明 `result = h - X`
- [ ] L2: 建模 `ddf_Zp`，证明循环不变量（引理 3.1）
- [ ] L2: 证明 DDF 分裂正确性（定理 3.2）
- [ ] L2: 证明提前终止正确性（命题 3.3）
- [ ] L1: 建模 `U64`/`I64`，证明 B3 bug 的 old-wrong/new-correct
- [ ] 评估报告：Mathlib 缺口、证明难度、后续工作量估算

### Phase 2：Zp 完整管线（3-4 周）

- [ ] L2: 建模 `squarefree_Zp`，证明输出无平方 + 乘积还原
- [ ] L2: 建模 `edf_Zp`，证明分裂概率下界
- [ ] L2: 建模 `factor_Zp`，证明顶层正确性定理
- [ ] L1: 完成 Vec/Move/迭代器模型
- [ ] L1: 对所有函数验证数组越界安全

### Phase 3：Z[x] Hensel + 重组（4-6 周）

- [ ] L3: 形式化 Hensel 引理（多项式 p-adic 提升）
- [ ] L3: 形式化 Mignotte 界
- [ ] L2: 建模 `__hensel_lift`
- [ ] L2: 建模 `__zassenhaus_recombine`
- [ ] L2: 建模 `__lll_factorize`，证明 Phase 1 + Phase 2 完整性
- [ ] L2: 顶层定理：`factor_ZZ_correct`

### Phase 4：LLL 验证（可选，高难度）

- [ ] L3: 形式化 LLL 算法正确性（格基约化保证）
- [ ] L2: 建模 `__lll_reduce`
- [ ] L2: 证明 van Hoeij 重组的因子绑定正确性

## 7. Path C 执行流程

### 7.1 单函数验证的完整工作流

以 `__ddf_Zp` 为例，每个函数走完以下 6 步：

```
Step 1: 提取           C++ 源码 → 标注版 C++（标记控制流 + 变量生命期）
Step 2: 翻译           标注版 C++ → Lean 4 函数定义（1:1 对应）
Step 3: 陈述           写出正确性定理（前条件 + 后条件 + 不变量）
Step 4: 证明（L2）     证明算法正确性（假设底层操作正确）
Step 5: 证明（L1）     证明表示层安全（整数溢出、数组越界、move）
Step 6: 审查           人工核对 Lean 模型 ↔ C++ 代码的逐行对应
```

### 7.2 Step 1：提取（C++ → 标注版 C++）

在原始 C++ 代码上标注，不改动代码本身。标注内容：

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

### 7.3 Step 2：翻译（标注版 C++ → Lean 4）

翻译规则表：

| C++ 构造 | Lean 4 对应 | 注意事项 |
|---------|------------|---------|
| `for (init; cond; step) { body }` | `let rec loop (vars) := if cond then body; loop (step vars) else vars` | 需 well-founded 递减证明 |
| `while (true) { ... break; }` | `partial def` 或 `loop` + 燃料参数 | EDF 用 `partial` |
| `uint64_t x = expr;` | `let x : U64 := expr` | 模 2^64 语义 |
| `(int64_t)(p - 1)` | `cast_u64_to_i64 (p - 1)` | 显式转换，可证明溢出 |
| `vec.push_back(x)` | `let vec' := vec.push_back x` | 旧 `vec` 不再使用 |
| `vec[i]` | `vec.get i (by ...)` | 需提供 `i < vec.size` 证明 |
| `auto x = std::move(y)` | `let (x, y') := Ownership.move_from y (by ...)` | `y'` 标记为 moved |
| `pair_vec_div(q, f, g, ...)` | `let (q, r) := poly_divmod f g` | 精确整除时 `r = 0` |
| `polynomial_GCD(a, b)` | `Polynomial.gcd a b` | 对接 Mathlib 或自定义 |

**1:1 对应原则**：Lean 函数的控制流结构必须与 C++ 一一对应。不允许"重写为更优雅的 Lean 风格"——否则 Step 6 的审查无法进行。

DDF 的翻译结果：

```lean
/-- __ddf_Zp 的 1:1 Lean 模型 -/
def ddf_Zp (f : ZpPoly) (p : U64) : List (ZpPoly × Nat) :=
  let h₀ : ZpPoly := X                       -- h = x
  let f_star₀ : ZpPoly := f
  let rec loop (h : ZpPoly) (f_star : ZpPoly) (d : Nat)
               (acc : List (ZpPoly × Nat))
               (h_deg_bound : f_star.deg ≥ 0 → h.deg < f_star.deg)
               : List (ZpPoly × Nat) :=
    if h_exit : f_star.deg < 2 * d then
      -- 提前终止：f_star 本身是不可约（或为 1）
      if f_star.deg > 0 then acc ++ [(f_star.monic, f_star.natDeg)]
      else acc
    else
      let h' := upoly_powmod h p f_star       -- h = h^p mod f*
      let h_minus_x := upoly_subtract_x h' p  -- h - x
      let gd := poly_gcd h_minus_x f_star     -- gcd(h-x, f*)
      if gd.deg > 0 then
        let f_star' := poly_exact_div f_star gd  -- f* / gd
        let h'' := upoly_mod h' f_star'           -- h mod new_f*
        loop h'' f_star' (d + 1) (acc ++ [(gd.monic, d)])
             (by ...)   -- 证明 h''.deg < f_star'.deg
      else
        loop h' f_star (d + 1) acc
             (by ...)   -- 不变量保持
  termination_by f_star.natDeg - 2 * d   -- well-founded: gap 递减或 d 递增
  loop h₀ f_star₀ 1 [] (by ...)
```

### 7.4 Step 3：陈述（写正确性定理）

对每个函数，陈述三类定理：

**A. 功能正确性**（L2 层）
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

**B. 循环不变量**（L2 层，辅助引理）
```lean
theorem ddf_loop_invariant (h f_star : ZpPoly) (d : Nat) :
    -- 在第 d 次迭代开始时（powmod 之后）
    h ≡ X ^ (p ^ d) [MOD f_star]
    -- f_star 的所有不可约因子度 ≥ d
    ∧ (∀ q, Irreducible q → q ∣ f_star → q.natDeg ≥ d) := by
  sorry
```

**C. 表示层安全**（L1 层）
```lean
theorem ddf_Zp_no_overflow (f : ZpPoly) (p : U64) (hp : p.val ≥ 2) :
    -- subtract_x 中的 p-1 不溢出
    ∀ h encountered in ddf_Zp execution,
      upoly_subtract_x h p 计算中无 U64 溢出
    -- 所有数组访问都在界内
    ∧ all_array_accesses_in_bounds (ddf_Zp f p) := by
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

L1 证明策略：
- `p - 1` 不溢出：`p ≥ 2 → p - 1 ≥ 1`，U64 减法安全
- 数组越界：追踪 `_coeffs.size()` = `deg + 1`，每次访问的 index ≤ deg

### 7.6 Step 6：人工审查

审查清单（每个函数必须通过）：

| # | 检查项 | 方法 |
|---|--------|------|
| 1 | Lean 函数与 C++ 函数的**控制流同构** | 逐行对照，确认分支/循环结构一致 |
| 2 | 每个 C++ 变量在 Lean 中有**同名对应** | 变量映射表 |
| 3 | 每个 C++ 函数调用在 Lean 中有**对应调用** | 调用图比对 |
| 4 | C++ 的隐式转换在 Lean 中**显式建模** | 检查所有 int/uint 赋值 |
| 5 | `std::move` 在 Lean 中有 `Ownership.move_from` 对应 | move 点列表比对 |
| 6 | 数组访问的越界证明的**前提条件**与 C++ 的**运行时保证**一致 | 逐个检查 |

审查产物：`lean/CLPoly/Review/DDF_review.md`，记录每个检查项的结论。

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
│   └── Repr/                      -- L1: 表示层
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
