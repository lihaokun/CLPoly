# L1 实现模型：设计与精化方案

> 状态：nl-proof v1（架构设计 + SQF 示例）
> 目标：1:1 对应 C++ 控制流，精化证明 L1 → L2

---

## 0. 总体架构

### 0.1 分层关系

```
L2 算法模型（已完成，3257 行，0 sorry）
  操作对象：Polynomial (ZMod p)，EuclideanDomain.gcd，Mathlib 类型
  控制流：数学归纳/递归，抽象掉 C++ 细节

    ↓ 精化证明（refine_*：L1 行为 = L2 行为）

L1 实现模型（本文档规划）
  操作对象：U64, SparsePolyZp, 模仿 C++ 数据结构
  控制流：1:1 对应 C++ 的 for/while/if
```

### 0.2 文件结构

```
CLPoly/Impl/
  Types.lean         — U64, Zp, SparsePolyZp 类型定义
  ZpArith.lean       — Zp 算术（add, mul, inv, Barrett 模乘）
  PolyArith.lean     — 稀疏多项式算术（mul, div, gcd, derivative）
  SquarefreeZp.lean  — __squarefree_Zp 的 L1 模型
  DDF.lean           — __ddf_Zp 的 L1 模型
  EDF.lean           — __edf_Zp 的 L1 模型
  Hensel.lean        — __hensel_lift 的 L1 模型
  Recombine.lean     — __zassenhaus_recombine 的 L1 模型
  Refinement/        — 精化证明（L1 → L2）
```

---

## 1. C++ 类型 → Lean L1 类型

### 1.1 设计决策

**简化**：不单独建模 `Zp` 类。C++ 的 `Zp` 包含值 + 素数 + Barrett 预计算，但每个多项式的所有系数共享同一个素数。在 L1 中，素数由 `SparsePolyZp` 统一管理，系数用裸 `UInt64` 值表示。

**度数表示**：C++ 用 `uint64_t`（通过 `umonomial`），L1 用 `Nat`。简化理由：CLPoly 不处理度数 > 2^64 的多项式，用 `Nat` 避免溢出证明的复杂性。此简化不影响算法正确性验证。

**Barrett 预计算**：省略。L1 只建模模运算的数学语义（`a * b % p`），不验证 Barrett 算法内部。精化时只需证 `impl_mul(a, b, p) = (a * b) % p`。

### 1.2 SparsePolyZp（核心类型）

C++:
```cpp
// upolynomial_<Zp> = vector<pair<umonomial, Zp>>
// 降幂排列，(degree, coefficient) 对
// Zp 包含 _i (值), _p (素数), _ninv, _norm
// 不变量：度数严格递减，系数非零，所有系数共享同一素数
```

Lean L1:
```lean
/-- 稀疏多项式 over Z/pZ，1:1 对应 C++ 的 upolynomial_<Zp>。
    系数用裸 UInt64 值（∈ [0, p)），素数在结构级共享。 -/
structure SparsePolyZp where
  terms : List (Nat × UInt64)  -- (度数, 系数值) 降序排列
  prime : Nat                   -- 素数 p
  h_prime : Nat.Prime prime     -- p 是素数
  h_sorted : terms.Chain' (fun a b => a.1 > b.1)  -- 度数严格降序
  h_nonzero : ∀ t ∈ terms, t.2.toNat ≠ 0           -- 系数值非零
  h_range : ∀ t ∈ terms, t.2.toNat < prime          -- 系数值 ∈ [0, p)
```

**便利方法**：
```lean
def SparsePolyZp.deg (sp : SparsePolyZp) : Nat :=
  sp.terms.head?.map Prod.fst |>.getD 0

def SparsePolyZp.isEmpty (sp : SparsePolyZp) : Bool :=
  sp.terms.isEmpty

def SparsePolyZp.leadCoeff (sp : SparsePolyZp) : UInt64 :=
  sp.terms.head?.map Prod.snd |>.getD 0
```

### 1.3 精化关系

```lean
/-- 将 L1 稀疏多项式转换为 L2 数学多项式 -/
noncomputable def SparsePolyZp.toPoly (sp : SparsePolyZp) :
    Polynomial (ZMod sp.prime) :=
  sp.terms.foldl
    (fun acc (d, c) => acc + Polynomial.C (c.toNat : ZMod sp.prime) * Polynomial.X ^ d)
    0
```

**精化定理模板**：对每个 L1 函数 `f_impl`，证明：
```lean
theorem refine_f (input : SparsePolyZp) (h_inv : input.isValid) :
    (f_impl input).toPoly = f_l2 input.toPoly
```

### 1.4 Zp 算术（内联操作）

不单独定义 `ZpImpl` 类型。Zp 算术作为 `SparsePolyZp` 上操作的内部步骤：

```lean
/-- 模加（对应 C++ Zp::operator+） -/
def zp_add (a b p : Nat) (ha : a < p) (hb : b < p) : Nat :=
  if a + b < p then a + b else a + b - p

/-- 模乘（对应 C++ Zp::__nmod_mul，Barrett 省略） -/
def zp_mul (a b p : Nat) (ha : a < p) (hb : b < p) : Nat :=
  (a * b) % p

/-- 模逆（对应 C++ inv_prime，扩展欧几里得） -/
def zp_inv (a p : Nat) (ha : a ≠ 0) (hp : Nat.Prime p) : Nat :=
  -- 扩展欧几里得算法
  sorry -- 定义在 ZpArith.lean 中
```

**精化定理**：
```lean
theorem refine_zp_add : (zp_add a b p ha hb : Nat) = (a + b : ZMod p).val
theorem refine_zp_mul : (zp_mul a b p ha hb : Nat) = (a * b : ZMod p).val
theorem refine_zp_inv : (zp_inv a p ha hp : Nat) = (a⁻¹ : ZMod p).val
```

---

## 2. SQF L1 模型（`__squarefree_Zp`）

### 2.1 C++ 结构（line 122-180）

```
__squarefree_Zp(f):
  if f' = 0:
    g = extract_pth_root(f)
    make_monic(g)
    sub = __squarefree_Zp(g)
    return [(s, e*p) for (s,e) in sub]

  c = GCD(f, f')
  w = f / c
  normalize(w)

  i = 1
  while deg(w) > 0:
    y = GCD(w, c)
    z = w / y; normalize(z)
    if deg(z) > 0:
      make_monic(z)
      result.push(z, i)
    c = c / y; normalize(c)
    w = y
    i++

  if deg(c) > 0:
    g = extract_pth_root(c)
    make_monic(g)
    sub = __squarefree_Zp(g)
    result.append([(s, e*p) for (s,e) in sub])

  return result
```

### 2.2 Lean L1 模型

```lean
-- Yun 内循环（对应 C++ while loop）
def yunLoop_impl (w c : SparsePolyZp) (i : Nat)
    (acc : List (SparsePolyZp × Nat)) :
    List (SparsePolyZp × Nat) × SparsePolyZp :=
  if w.deg = 0 then (acc, c)
  else
    let y := gcd_impl w c
    let z := div_impl w y |>.normalize
    let acc' := if z.deg > 0 then acc ++ [(z.makeMonic, i)] else acc
    let c' := div_impl c y |>.normalize
    yunLoop_impl y c' (i + 1) acc'
termination_by w.deg + c.deg  -- 同 L2

-- 顶层 SQF（对应 C++ __squarefree_Zp）
def squarefree_impl (f : SparsePolyZp) :
    List (SparsePolyZp × Nat) :=
  let f_deriv := derivative_impl f
  if f_deriv.isEmpty then
    let g := extractPthRoot_impl f |>.makeMonic
    (squarefree_impl g).map (fun (s, e) => (s, e * f.prime.toNat))
  else
    let c := gcd_impl f f_deriv
    let w := div_impl f c |>.normalize
    let (yun_result, c_rem) := yunLoop_impl w c 1 []
    if c_rem.deg > 0 then
      let g := extractPthRoot_impl c_rem |>.makeMonic
      yun_result ++ (squarefree_impl g).map (fun (s, e) => (s, e * f.prime.toNat))
    else
      yun_result
```

### 2.3 精化定理

```lean
theorem refine_squarefree (f : SparsePolyZp) (hf : f.isValid) :
    (squarefree_impl f).map (fun (s, e) => (s.toPoly, e)) =
    sqfZp f.toPoly  -- L2 算法
```

**前提**：需要先证明子操作的精化：
- `refine_gcd`：`gcd_impl(a, b).toPoly = EuclideanDomain.gcd a.toPoly b.toPoly`
- `refine_div`：`div_impl(a, b).toPoly = a.toPoly /ₘ b.toPoly`（when b monic）
- `refine_derivative`：`derivative_impl(f).toPoly = derivative f.toPoly`
- `refine_extractPthRoot`：`extractPthRoot_impl(f).toPoly = Polynomial.contract p f.toPoly`

---

## 3. 基本操作精化（依赖链）

### 3.1 依赖图

```
refine_squarefree
  ← refine_yunLoop
  ← refine_gcd         ← refine_divmod
  ← refine_div          ← refine_divmod
  ← refine_derivative   ← 直接系数操作
  ← refine_extractPthRoot ← 直接系数操作
  ← refine_makeMonic    ← refine_zp_inv
```

### 3.2 GCD 实现与精化（关键难点）

**C++ GCD**（polynomial_gcd.hh）：标准 Euclidean 算法 + normalize（首一化）。

**L1 GCD**：
```lean
/-- Euclidean GCD for sparse polynomials over Zp.
    1:1 对应 C++ polynomial_GCD。 -/
def poly_gcd_impl (a b : SparsePolyZp) : SparsePolyZp :=
  if b.isEmpty then a.makeMonic
  else
    let (_, r) := poly_divmod_impl a b
    poly_gcd_impl b r
termination_by b.deg
```

**精化难点**：Mathlib 的 `EuclideanDomain.gcd` 是**抽象定义**（通过 well-founded recursion），
不直接等于 Euclidean 算法的具体执行。

**解决路径**：不直接证 `gcd_impl.toPoly = EuclideanDomain.gcd a.toPoly b.toPoly`。
而是证 GCD 的**特征性质**：
```lean
theorem refine_poly_gcd (a b : SparsePolyZp) :
    let g := poly_gcd_impl a b
    -- g 整除 a 和 b
    g.toPoly ∣ a.toPoly ∧ g.toPoly ∣ b.toPoly
    -- g 是最大公因子
    ∧ ∀ d : Polynomial (ZMod p), d ∣ a.toPoly → d ∣ b.toPoly → d ∣ g.toPoly
    -- g 是 monic（normalize 后）
    ∧ Monic g.toPoly
```

这避免了与 Mathlib 抽象 GCD 的精确匹配问题。L2 证明中使用 GCD 的地方只依赖
整除性和互素性，不依赖 GCD 的具体值（`normalize` 保证唯一到 unit）。

### 3.3 多项式除法

**L1 实现**（对应 C++ `pair_vec_div`）：
```lean
/-- 多项式带余除法。b 必须 monic（或首项系数可逆）。
    返回 (商, 余式)。 -/
def poly_divmod_impl (a b : SparsePolyZp) :
    SparsePolyZp × SparsePolyZp :=
  -- 标准多项式长除法：
  -- 每步：q_coeff = lc(remainder) / lc(b)
  --       remainder -= q_coeff * x^(deg(r)-deg(b)) * b
  sorry -- 具体实现在 PolyArith.lean
```

**精化定理**：
```lean
theorem refine_poly_divmod (a b : SparsePolyZp) (hb : b.leadCoeff ≠ 0) :
    let (q, r) := poly_divmod_impl a b
    a.toPoly = b.toPoly * q.toPoly + r.toPoly
    ∧ r.deg < b.deg
```

### 3.4 导数和 p 次根

**导数**（直接系数操作）：
```lean
/-- 形式导数：对每项 (d, c) 映射为 (d-1, d*c mod p)，去掉 d=0 项。 -/
def derivative_impl (f : SparsePolyZp) : SparsePolyZp :=
  { terms := f.terms.filterMap (fun (d, c) =>
      if d = 0 then none else some (d - 1, zp_mul d.toUInt64 c f.prime ...))
    ... }
```

**p 次根提取**（直接系数操作）：
```lean
/-- 提取 p 次根：f(x) = g(x^p) → g。对每项 (d, c) 映射为 (d/p, c)。 -/
def extractPthRoot_impl (f : SparsePolyZp) : SparsePolyZp :=
  { terms := f.terms.map (fun (d, c) => (d / f.prime, c))
    ... }
```

这两个的精化比较直接：逐系数对应 Mathlib 的 `Polynomial.derivative` 和 `Polynomial.contract p`。

---

## 4. 实施策略

### 4.1 自底向上

1. **Phase L1.0**：类型定义（Types.lean，~100 行）
2. **Phase L1.1**：Zp 算术（ZpArith.lean，~150 行）+ 精化
3. **Phase L1.2**：多项式算术（PolyArith.lean，~300 行）+ 精化
4. **Phase L1.3**：GCD（~150 行）+ 精化
5. **Phase L1.4**：SQF L1（~200 行）+ 精化
6. **Phase L1.5**：DDF/EDF/Hensel/Recombine L1（~500 行总计）

### 4.2 估计

| Phase | 行数 | 依赖 |
|-------|------|------|
| L1.0 类型 | ~100 | — |
| L1.1 Zp 算术 | ~150 | L1.0 |
| L1.2 多项式算术 | ~300 | L1.1 |
| L1.3 GCD | ~150 | L1.2 |
| L1.4 SQF | ~200 | L1.3 |
| L1.5 其余算法 | ~500 | L1.3, L1.4 |
| **总计** | **~1400** | |

### 4.3 关键决策

1. **Barrett 模乘**：L1 中建模但精化时只证 `val(a*b mod p) = a*b mod p`，不验证 Barrett 算法内部
2. **稀疏表示不变量**：排序 + 非零系数。每个操作必须维护。
3. **溢出安全**：`uint64_t` 加法/乘法的溢出行为。CLPoly 的 Zp 算术已处理（`nmod_add`/`nmod_mul` 不溢出）。但 L1 需要显式证明。
4. **终止性**：L1 函数的终止性必须从 C++ 的循环不变量推导，对应 L2 的 `termination_by`。

---

## 5. Phase L1.0 内容（Types.lean）

与 §1.2 一致的最终定义：

```lean
/-- 稀疏多项式 over Z/pZ。1:1 对应 C++ upolynomial_<Zp>。
    系数用 UInt64 值（∈ [0, p)），素数在结构级共享。
    不变量：度数严格降序 + 系数非零 + 系数在范围内。 -/
structure SparsePolyZp where
  terms : List (Nat × UInt64)
  prime : Nat
  h_prime : Nat.Prime prime
  h_sorted : terms.Chain' (fun a b => a.1 > b.1)
  h_nonzero : ∀ t ∈ terms, t.2.toNat ≠ 0
  h_range : ∀ t ∈ terms, t.2.toNat < prime

-- 便利方法
def SparsePolyZp.deg (sp : SparsePolyZp) : Nat :=
  sp.terms.head?.map Prod.fst |>.getD 0

def SparsePolyZp.isEmpty (sp : SparsePolyZp) : Bool :=
  sp.terms.isEmpty

def SparsePolyZp.leadCoeff (sp : SparsePolyZp) : UInt64 :=
  sp.terms.head?.map Prod.snd |>.getD 0

-- 精化到 Polynomial (ZMod p)
noncomputable def SparsePolyZp.toPoly (sp : SparsePolyZp) :
    Polynomial (ZMod sp.prime) :=
  sp.terms.foldl
    (fun acc (d, c) => acc + Polynomial.C (c.toNat : ZMod sp.prime) * Polynomial.X ^ d)
    0

-- Zp 算术（内联辅助）
def zp_add (a b p : Nat) : Nat := if a + b < p then a + b else a + b - p
def zp_mul (a b p : Nat) : Nat := (a * b) % p
```
