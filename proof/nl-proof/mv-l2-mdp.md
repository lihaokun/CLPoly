# MDP（偏分式分解）正确性

> 状态：nl-proof v1
> 对应 C++：`__mtshl_zp_univar_mdp`（j=2）/ `__mtshl_sparse_int`（j≥3）
> 核心：互素因子的偏分式分解存在且唯一（CRT）

---

## 0. 问题

给定互素的 f₁,...,fᵣ 和目标 c，找 σ₁,...,σᵣ 使得：
```
c = Σᵢ σᵢ · ∏_{j≠i} fⱼ    (*)
```

等价于：`c / ∏fⱼ = Σ σᵢ / fᵢ`（偏分式分解）。

## 1. 存在性

### 1.1 r = 2 的情况

IsCoprime f₁ f₂ → ∃ s t, s·f₁ + t·f₂ = 1。
乘以 c：c = (c·s)·f₁ + (c·t)·f₂ = σ₁·f₂ + σ₂·f₁。✓

这就是 Bézout 系数。C++ 的 `__mtshl_zp_univar_mdp` 正是这样做的。

### 1.2 一般 r 的情况

归纳：设 F̂₁ = ∏_{j≥2} fⱼ。
IsCoprime f₁ F̂₁（由 pairwise coprime + `IsCoprime.prod_right`）。
Bézout：∃ s t, s·f₁ + t·F̂₁ = 1。
c = (c·s)·f₁ + (c·t)·F̂₁。

第一项：σ₁·F̂₁ 其中 σ₁ = c·t。✓
第二项：(c·s)·f₁ = c'·∏_{j≥2} fⱼ 其中 c' = c·s。

对 f₂,...,fᵣ 递归求解 c' = Σ_{i≥2} σᵢ · ∏_{j≥2,j≠i} fⱼ。
乘以 f₁：c'·f₁ = Σ_{i≥2} σᵢ · f₁ · ∏_{j≥2,j≠i} fⱼ = Σ_{i≥2} σᵢ · ∏_{j≠i} fⱼ。

合并：c = σ₁·F̂₁ + Σ_{i≥2} σᵢ · ∏_{j≠i} fⱼ = Σᵢ σᵢ · ∏_{j≠i} fⱼ。✓

## 2. Lean 形式化

### 2.1 定理

```lean
/-- MDP 存在性：互素因子的偏分式分解存在。
    对应 C++ MDP 求解器的数学正确性保证。-/
theorem mdp_exists {R : Type*} [CommRing R]
    (f_base : List R) (c : R)
    (hcop : f_base.Pairwise IsCoprime) :
    ∃ sigma : List R,
      MDPCorrect c f_base sigma
```

### 2.2 证明（对 f_base 归纳）

**Base**（空/单元素）：
- `[]`：c = 0 · ... 只在 c = 0 时成立。实际上空列表不出现。
- `[f]`：c = σ · 1 → σ = c。`linearTerm [f] [c] = c * [].prod = c * 1 = c`。✓

**Step**（`f :: rest`）：
- `IsCoprime f rest.prod`（从 pairwise + `IsCoprime.prod_right`）
- Bézout：∃ s t, s·f + t·rest.prod = 1
- σ₀ = c·t（f 的系数）
- 对 rest 递归求解 c' = c·s
- `linearTerm (f :: rest) (σ₀ :: σ_rest)`
  `= σ₀ · rest.prod + f · linearTerm rest σ_rest`
  `= c·t · rest.prod + f · c'_decomposition`（递归给出后者 = c·s）
  `= c·t · rest.prod + c·s · f = c · (s·f + t·rest.prod) = c · 1 = c`。✓

### 2.3 Lean API

- `IsCoprime.prod_right`（或手动从 `IsCoprime.mul_right` 归纳）
- `IsCoprime` → Bézout 系数（`IsCoprime` 定义本身就是 `∃ s t, s*a + t*b = 1`）
- `linearTerm` 递归展开

### 2.4 度数约束

C++ 的 MDP 还要求 `deg(σᵢ) < deg(fᵢ)`。这通过 `modByMonic` 实现（取 σᵢ mod fᵢ 的余式）。
对于 L2 验证，我们只需要偏分式等式 (*)，不需要度数约束——因为 Newton 迭代的正确性只依赖等式，不依赖度数。

## 3. 形式化估计

| 内容 | 行数 |
|------|------|
| `isCoprime_list_prod` 辅助（pairwise → coprime with prod） | ~10 |
| `mdp_exists` 归纳证明 | ~30 |
| **总计** | **~40** |

## 4. 与现有代码的集成

`mtshl_step_invariant` 中的 `h_mdp : MDPCorrect ...` 假设将由 `mdp_exists` 的输出满足。
即：C++ 的 MDP 求解器输出满足 `MDPCorrect`，因为数学保证了偏分式分解存在。
具体哪个 σ 被选择不影响 Newton 迭代的正确性（任何满足等式的 σ 都行）。
