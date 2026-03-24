# Zassenhaus Recombination 算法 L2 模型

> 对应 C++: `__zassenhaus_recombine`（polynomial_factorize_univar.hh:750-882）

---

## 0. C++ 算法

```
input: f ∈ Z[x], lifted[0..r-1] (Hensel 因子, monic), m (精度)
output: irreducible factors of f

T = {0,...,r-1}    // 可用因子索引
f* = f              // 剩余多项式
result = []

for s = 1 to ⌊|T|/2⌋:
  for each s-subset S of T:
    g = (lc(f*) · ∏_{i∈S} lifted[i]) mod_sym m
    if g | f*:                    // trial division
      result.push(pp(g))
      f* = pp(f* / g)
      T = T \ S
      reset s = 1; break

if deg(f*) > 0: result.push(f*)
return result
```

## 1. L2 模型：循环不变量

### 1.1 不变量定义

```lean
/-- Zassenhaus 循环不变量：
    f ∼ remaining × ∏extracted，每个 extracted 不可约。-/
structure ZassenhausInvariant (f remaining : Polynomial ℤ)
    (extracted : List (Polynomial ℤ)) : Prop where
  prod_eq : Associated f (remaining * extracted.prod)
  all_irred : ∀ g ∈ extracted, Irreducible g
```

### 1.2 初始化

```lean
theorem zassenhaus_init (f : Polynomial ℤ) :
    ZassenhausInvariant f f [] :=
  ⟨by simp, fun _ h => absurd h (List.not_mem_nil _)⟩
```

**证明**：`f ∼ f * [].prod = f * 1 = f`。空列表无元素。

### 1.3 因子提取步

对应 C++ 循环体中 "g | f* → 提取 pp(g)" 的逻辑。

```lean
theorem zassenhaus_extract
    (f remaining : Polynomial ℤ) (extracted : List (Polynomial ℤ))
    (h_inv : ZassenhausInvariant f remaining extracted)
    (g : Polynomial ℤ) (hg_irred : Irreducible g)
    (remaining' : Polynomial ℤ) (h_div : remaining = g * remaining') :
    ZassenhausInvariant f remaining' (g :: extracted)
```

**证明**：
- `f ∼ remaining * extracted.prod = (g * remaining') * extracted.prod`
- `= remaining' * (g * extracted.prod) = remaining' * (g :: extracted).prod`
- 新 extracted 的每个元素：g 不可约 + 旧 extracted 不可约。

### 1.4 终止（remaining 不可约）

对应 C++ 循环结束后 `if deg(f*) > 0: result.push(f*)`。

```lean
theorem zassenhaus_terminate_irred
    (f remaining : Polynomial ℤ) (extracted : List (Polynomial ℤ))
    (h_inv : ZassenhausInvariant f remaining extracted)
    (h_irred : Irreducible remaining) :
    RecombineCorrect f (remaining :: extracted)
```

**证明**：
- `f ∼ remaining * extracted.prod = (remaining :: extracted).prod`
- 所有元素：remaining 不可约 + extracted 不可约。

### 1.5 终止（remaining 是 unit）

```lean
theorem zassenhaus_terminate_unit
    (f remaining : Polynomial ℤ) (extracted : List (Polynomial ℤ))
    (h_inv : ZassenhausInvariant f remaining extracted)
    (h_unit : IsUnit remaining) :
    RecombineCorrect f extracted
```

**证明**：
- `f ∼ remaining * extracted.prod ∼ extracted.prod`（乘 unit）
- 所有元素不可约。

## 2. recombine_correct 的证明

```lean
theorem recombine_correct (f : Polynomial ℤ) (hf : f ≠ 0)
    : ∃ result, RecombineCorrect f result := by
  -- 初始化
  have h0 := zassenhaus_init f
  -- 迭代提取：对 remaining 使用 WfDvdMonoid 归纳
  -- （Z[x] 是 Noetherian → 不可约因子链有限）
  -- 每次提取减少 remaining 的因子数
  -- 终止时 remaining 不可约或 unit
```

**具体证明**：用 `WfDvdMonoid.exists_irreducible_factor` 证明每个非 unit、非不可约的 remaining 都有不可约因子（即 C++ 的子集枚举一定找到一个）。提取后 remaining 严格变小（在 dvd 良基序下）。终止时应用 `zassenhaus_terminate_irred` 或 `zassenhaus_terminate_unit`。

## 3. 行数估计

| 定理 | 行数 |
|------|------|
| `ZassenhausInvariant` 定义 | ~8 |
| `zassenhaus_init` | ~3 |
| `zassenhaus_extract` | ~10 |
| `zassenhaus_terminate_irred` | ~8 |
| `zassenhaus_terminate_unit` | ~8 |
| 更新 `recombine_correct` | ~20 |
| **总计** | ~57 |
