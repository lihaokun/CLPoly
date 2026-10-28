# 多变量 L2 模块 2：LC 分配（Valuation 提取）

> 状态：nl-proof v1
> 对应 C++：`__wang_leading_coeff`（lines 1325-1489）
> 5/11 算法 bug 出在此模块

---

## 0. C++ 算法逻辑

```
输入：f, univar_factors [u₁,...,uᵣ], eval_point α, uni_content cs
输出：σ₁,...,σᵣ ∈ Z[x₂,...,xₙ]，满足 ∏σᵢ ~ L = lc(f, x₁)

算法：
1. L = lc(f, x₁), δ = L(α)
2. factorize(L) → γ · ∏ lⱼ^eⱼ（递归分解 L）
3. Eⱼ = |lⱼ(α)|（求值每个 LC 因子）
   检查：Eⱼ ≥ 2（否则换点）
4. Non-divisor 检查：每个 Eⱼ 剥离与 cs·γ 及之前 E 的共享素因子后 > 1
5. wᵢ = |lc(uᵢ) · cs|（单变量因子的 LC × content）
6. Valuation 提取：
   for j = m-1 downto 0:
     for i = 0 to r-1:
       kᵢⱼ = 0
       while Eⱼ | wᵢ:
         wᵢ = wᵢ / Eⱼ
         kᵢⱼ++
       σᵢ *= lⱼ^kᵢⱼ
7. 守恒验证：∀j, Σᵢ kᵢⱼ = eⱼ（否则换点）
8. γ 吸收到 σ₀
```

---

## 1. L2 模型

### 1.1 Valuation 提取函数

```lean
/-- 整数 valuation：v_E(w) = max{k : E^k | w}。
    对应 C++ 的 while(R % Ej == 0) 循环。-/
noncomputable def intValuation (E w : ℤ) (hE : 1 < E.natAbs) : ℕ :=
  -- 存在最大 k 使得 E^k | w（因 E > 1，|E^k| 严格递增）
  Nat.find (exists_max_pow_dvd E w hE)
```

实际上用 `multiplicity E w` 更好（Mathlib 已有）。

```lean
/-- LC 分配核心循环：对每个 LC 因子 (lⱼ, eⱼ)，对每个单变量因子 uᵢ，
    提取 kᵢⱼ = v_{Eⱼ}(wᵢ) 并累积 σᵢ *= lⱼ^kᵢⱼ。
    对应 C++ 的双重 for 循环（lines 1431-1460）。-/
def lcDistribLoop
    -- LC 因子和它们的求值值
    (lc_factors : List (MvPolynomial (Fin n) ℤ × ℕ))  -- (lⱼ, eⱼ)
    (lc_evals : List ℤ)  -- Eⱼ = |lⱼ(α)|
    -- 单变量因子的 LC 值
    (w : List ℤ)  -- wᵢ = |lc(uᵢ) · cs|
    : List (MvPolynomial (Fin n) ℤ)  -- σ₁,...,σᵣ
```

### 1.2 输出规约

```lean
/-- LC 分配结果满足的条件。对应 C++ __wang_lc_result 的隐含不变量。-/
structure LCDistribCorrect {n : ℕ}
    (L : MvPolynomial (Fin n) ℤ)       -- lc(f, x₁)
    (γ : ℤ)                             -- content(L 的因式分解)
    (σ : List (MvPolynomial (Fin n) ℤ)) -- 输出
    (α : Fin n → ℤ)                     -- 求值点
    (univar_lcs : List ℤ)               -- lc(uᵢ)
    (cs : ℤ)                            -- uni_content
    : Prop where
  -- 1. 乘积条件：γ · ∏σᵢ ∼ L
  prod_eq : Associated L (C γ * σ.prod)
  -- 2. 求值一致：σᵢ(α) | wᵢ（wᵢ = lc(uᵢ)·cs）
  eval_dvd : List.Forall₂ (fun s w => MvPolynomial.eval α s ∣ w)
      σ (univar_lcs.map (· * cs))
  -- 3. 长度一致
  length_eq : σ.length = univar_lcs.length
```

---

## 2. 正确性证明

### 2.1 前置条件

```lean
structure LCDistribPrecond {n : ℕ}
    (lc_factors : List (MvPolynomial (Fin n) ℤ × ℕ))
    (lc_evals : List ℤ)
    (γ : ℤ) (cs : ℤ)
    : Prop where
  -- (C1) Eⱼ 两两 coprime
  evals_coprime : lc_evals.Pairwise (fun a b => Int.gcd a b = 1)
  -- (C2) Eⱼ 与 cs·γ 无共享素因子（non-divisor 检查后）
  evals_no_shared : ∀ E ∈ lc_evals, Int.gcd E (cs * γ) = 1
  -- (C3) Eⱼ ≥ 2
  evals_ge_two : ∀ E ∈ lc_evals, 1 < E.natAbs
  -- Eⱼ = |lⱼ(α)|
  evals_correct : List.Forall₂
      (fun lj E => E = (MvPolynomial.eval α lj.1).natAbs) lc_factors lc_evals
```

### 2.2 核心正确性定理

```lean
/-- Valuation 提取算法正确：
    在前置条件 (C1)-(C3) 下，若守恒验证 (C4) 通过，
    则输出 σ 满足 LCDistribCorrect。

    对应 C++ __wang_leading_coeff 的数学正确性。-/
theorem lc_distrib_correct
    (lc_factors : List (MvPolynomial (Fin n) ℤ × ℕ))
    (lc_evals : List ℤ)
    (γ : ℤ) (cs : ℤ) (α : Fin n → ℤ)
    (w : List ℤ)  -- wᵢ 初始值
    (hpre : LCDistribPrecond lc_factors lc_evals γ cs)
    -- L = γ · ∏ lⱼ^eⱼ
    (hL : Associated L (C γ * (lc_factors.map (fun pr => pr.1 ^ pr.2)).prod))
    -- 守恒 (C4)：Σᵢ kᵢⱼ = eⱼ
    (hconserv : ∀ j, j < lc_factors.length →
        (List.range w.length).sum (fun i => valuation (lc_evals.get j) (w.get i)) =
        (lc_factors.get j).2)
    : LCDistribCorrect L γ (lcDistribLoop lc_factors lc_evals w) α
        (w.map (· / cs)) cs
```

### 2.3 证明思路

**乘积条件**：
- σᵢ = ∏ⱼ lⱼ^kᵢⱼ
- ∏ᵢ σᵢ = ∏ᵢ ∏ⱼ lⱼ^kᵢⱼ = ∏ⱼ lⱼ^(Σᵢ kᵢⱼ) = ∏ⱼ lⱼ^eⱼ
- 最后一步用守恒条件 (C4)：Σᵢ kᵢⱼ = eⱼ
- 所以 γ · ∏σᵢ = γ · ∏ lⱼ^eⱼ ∼ L ✓

**求值一致**：
- σᵢ(α) = ∏ⱼ lⱼ(α)^kᵢⱼ = ∏ⱼ (±Eⱼ)^kᵢⱼ
- 由 valuation 提取定义：Eⱼ^kᵢⱼ | wᵢ 且 Eⱼ^(kᵢⱼ+1) ∤ wᵢ
- 因 (C1) Eⱼ 两两 coprime：∏ⱼ Eⱼ^kᵢⱼ | wᵢ
- 所以 |σᵢ(α)| | wᵢ ✓

---

## 3. Non-divisor 检查模型

### 3.1 C++ 逻辑

```cpp
// lines 1410-1428
vector<ZZ> d = { abs(cs * gamma) };
for j = 0 to m-1:
    q = Eⱼ
    for k = j downto 0:
        r = d[k]
        while r != 1:
            r = gcd(r, q)
            q = q / r
            if q == 1: return fail  // Eⱼ 无独有素因子
    d.push_back(q)
```

### 3.2 数学意义

每次循环从 Eⱼ 中剥离与 d[0],...,d[j] 共享的所有素因子。
d[0] = |cs·γ|，d[k] (k ≥ 1) = Eₖ 剥离后的残余。

循环后 q > 1 说明 Eⱼ 有一个"独有素因子"——不被任何 d[k] 整除。
这保证 valuation 提取时 Eⱼ^k | wᵢ 能唯一确定 lⱼ 的贡献。

### 3.3 L2 模型

```lean
/-- Non-divisor 检查：每个 Eⱼ 在剥离与前序值的共享因子后有残余 > 1。
    对应 C++ lines 1410-1428。-/
def nonDivisorCheck (evals : List ℤ) (base : ℤ) : Bool :=
  -- 递归剥离共享素因子
  ...
```

不需要证明这个检查的**完备性**（即好的 α 总能通过检查）。
只需证明：**若检查通过，则条件 (C1)+(C2) 成立**（soundness）。

---

## 4. 形式化估计

| 内容 | 行数 |
|------|------|
| `intValuation` / `multiplicity` 使用 | ~10 |
| `lcDistribLoop` 模型 | ~30 |
| `LCDistribCorrect` 规约 | ~15 |
| `LCDistribPrecond` 规约 | ~15 |
| `nonDivisorCheck` 模型 + soundness | ~40 |
| `lc_distrib_correct` 证明 | ~60 |
| **总计** | **~170** |

---

## 5. 依赖的 Mathlib API

| 需要 | Mathlib 路径 |
|------|------------|
| `multiplicity` | `Mathlib.RingTheory.Multiplicity` |
| `Int.gcd` | `Mathlib.Data.Int.GCD` |
| Coprime integers valuation independence | 需自证（~20 行） |
| List product manipulation | `List.prod_map`, `List.Forall₂` |
