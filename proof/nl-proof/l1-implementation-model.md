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

### 1.1 Zp 类（number.hh:87-200）

C++:
```cpp
class Zp {
    uint64_t _i;      // 值 ∈ [0, p)
    uint64_t _p;      // 素数
    uint64_t _ninv;   // Barrett 预计算
    uint32_t _norm;   // clz(p)
};
```

Lean L1:
```lean
structure ZpImpl where
  val : UInt64       -- 值 ∈ [0, p)
  prime : UInt64     -- 素数 p
  -- Barrett 预计算省略（不影响数学正确性）
  h_val_lt : val.toNat < prime.toNat  -- 不变量
```

精化关系：`ZpImpl` 对应 `ZMod p` 中的元素。
```lean
def ZpImpl.toZMod (z : ZpImpl) : ZMod z.prime.toNat :=
  (z.val.toNat : ZMod z.prime.toNat)
```

### 1.2 upolynomial_<Zp>（稀疏多项式）

C++:
```cpp
// upolynomial_<Zp> = vector<pair<umonomial, Zp>>
// 降幂排列，(degree, coefficient) 对
// 不变量：度数严格递减，系数非零
```

Lean L1:
```lean
structure SparsePolyZp where
  terms : List (Nat × ZpImpl)  -- (度数, 系数) 对
  prime : UInt64
  -- 不变量
  h_sorted : terms.Chain' (fun a b => a.1 > b.1)  -- 严格降序
  h_nonzero : ∀ t ∈ terms, t.2.val ≠ 0             -- 系数非零
```

精化关系：`SparsePolyZp` 对应 `Polynomial (ZMod p)`。
```lean
def SparsePolyZp.toPoly (sp : SparsePolyZp) : Polynomial (ZMod sp.prime.toNat) :=
  sp.terms.foldl (fun acc (d, c) => acc + Polynomial.C c.toZMod * Polynomial.X ^ d) 0
```

### 1.3 精化定理模板

对每个 L1 函数 `f_impl`，证明：
```lean
theorem refine_f (input : L1_type) (h_inv : L1_invariant input) :
    f_impl(input).toPoly = f_l2(input.toPoly)
```

即：L1 函数的结果转换为数学对象 = L2 函数的结果。

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
  ← refine_gcd         ← refine_mul, refine_divmod
  ← refine_div          ← refine_divmod
  ← refine_derivative   ← refine_sub, refine_scalar_mul
  ← refine_extractPthRoot ← 直接构造
  ← refine_makeMonic    ← refine_scalar_mul, refine_inv
```

### 3.2 最底层：Zp 算术精化

```lean
-- Zp 加法
theorem refine_zp_add (a b : ZpImpl) :
    (zp_add_impl a b).toZMod = a.toZMod + b.toZMod

-- Zp 乘法（Barrett 模乘）
theorem refine_zp_mul (a b : ZpImpl) :
    (zp_mul_impl a b).toZMod = a.toZMod * b.toZMod

-- Zp 逆元
theorem refine_zp_inv (a : ZpImpl) (ha : a.val ≠ 0) :
    (zp_inv_impl a).toZMod = (a.toZMod)⁻¹
```

### 3.3 多项式算术精化

```lean
-- 稀疏多项式乘法
theorem refine_poly_mul (a b : SparsePolyZp) :
    (poly_mul_impl a b).toPoly = a.toPoly * b.toPoly

-- 稀疏多项式除法（带余）
theorem refine_poly_divmod (a b : SparsePolyZp) (hb : b.isMonic) :
    (poly_div_impl a b).toPoly = a.toPoly /ₘ b.toPoly ∧
    (poly_mod_impl a b).toPoly = a.toPoly %ₘ b.toPoly

-- Euclidean GCD
theorem refine_poly_gcd (a b : SparsePolyZp) :
    (poly_gcd_impl a b).toPoly = normalize (EuclideanDomain.gcd a.toPoly b.toPoly)
```

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

## 5. 从 Phase L1.0 开始

### 5.1 Types.lean 内容

```lean
-- U64 语义（Lean 的 UInt64 已提供）

-- Zp 实现类型
structure ZpImpl where
  val : UInt64
  prime : UInt64
  h_prime : Nat.Prime prime.toNat
  h_val_lt : val.toNat < prime.toNat

-- 精化到 ZMod
def ZpImpl.toZMod (z : ZpImpl) : ZMod z.prime.toNat :=
  (z.val.toNat : ZMod z.prime.toNat)

-- 稀疏多项式
structure SparsePolyZp where
  terms : List (Nat × UInt64)  -- (度数, 系数值) 降序
  prime : UInt64
  h_prime : Nat.Prime prime.toNat
  h_sorted : List.Chain' (fun a b => a.1 > b.1) terms
  h_nonzero : ∀ t ∈ terms, t.2 ≠ 0

-- 精化到 Polynomial (ZMod p)
noncomputable def SparsePolyZp.toPoly (sp : SparsePolyZp) :
    Polynomial (ZMod sp.prime.toNat) :=
  sp.terms.foldl
    (fun acc (d, c) => acc + Polynomial.C (c.toNat : ZMod sp.prime.toNat) * Polynomial.X ^ d)
    0
```
