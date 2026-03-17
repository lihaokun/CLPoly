# T2.2: 不可约多项式整除 X^{p^d} - X

> 状态：已形式化，0 sorry（无需 Monic 假设）

---

## 定理陈述

设 p 为素数，g ∈ F_p[x] 不可约，deg(g) | d。
则 g | (X^{p^d} - X)。

```lean
theorem irreducible_dvd_X_pow_sub_X
    (g : Polynomial (ZMod p)) (hg : Irreducible g)
    (d : ℕ) (hdvd : g.natDegree ∣ d) :
    g ∣ (X ^ (p ^ d) - X : Polynomial (ZMod p))
```

注：最终形式化去掉了 Monic 假设（`Irreducible.ne_zero` 足够构造 PowerBasis）。

## 证明

### Step 1: 构造扩域 K = F_p[x]/(g)

设 K = AdjoinRoot g，α = AdjoinRoot.root g。
由 g 不可约，K 是域（`AdjoinRoot.instField [Fact (Irreducible g)]`）。

### Step 2: card(K) = p^k

K 作为 F_p 上的向量空间，维数 = deg(g) = k（`AdjoinRoot.powerBasis`）。
因此 card(K) = card(F_p)^k = p^k。

**Mathlib 路径**：
- `AdjoinRoot.powerBasis (hg_ne : g ≠ 0)` 给出 PowerBasis，dim = natDegree g = k
- `Module.card_fintype` 或类似引理：`card K = (card (ZMod p)) ^ finrank`
- `ZMod.card p`：card(ZMod p) = p
- 组合得 card(K) = p^k

### Step 3: α^{p^d} = α

由 `FiniteField.pow_card_pow`：对有限域 K 中任意 a 和任意 n ∈ ℕ，
  a ^ (card K) ^ n = a

取 a = α, n = d/k（k | d 保证 d/k 是自然数），得：
  α ^ (p^k) ^ (d/k) = α

又 (p^k) ^ (d/k) = p ^ (k · (d/k)) = p^d（因为 k · (d/k) = d），所以：
  α ^ (p^d) = α

**Mathlib 路径**：
- `FiniteField.pow_card_pow (d/k) α`：`α ^ (card K) ^ (d/k) = α`
- `Nat.div_mul_cancel hdvd`：`k * (d / k) = d`（或 `d / k * k = d`）
- `pow_mul`：`p ^ (k * m) = (p^k)^m`
- 替换 card K = p^k，整理指数

### Step 4: g | X^{p^d} - X

α^{p^d} - α = 0（Step 3），即 α 是 X^{p^d} - X 的根。

在 AdjoinRoot g 中，mk g 是满射环同态 F_p[x] → K，核 = (g)。
所以 mk g (X^{p^d} - X) = aeval α (X^{p^d} - X) = α^{p^d} - α = 0。
由 mk g h = 0 ↔ g | h，得 g | X^{p^d} - X。

**Mathlib 路径**：
- `AdjoinRoot.aeval_eq p`：`aeval (root g) p = mk g p`
- 计算 `aeval α (X^{p^d} - X) = α^{p^d} - α = 0`（`map_sub`, `map_pow`, `aeval_X`）
- `AdjoinRoot.mk_eq_zero`：`mk g h = 0 ↔ g ∣ h`

## 预期难度

~50-100 行。主要工作：
1. 建立 AdjoinRoot 的有限域实例 + 计算 card（可能需要几个中间引理）
2. 指数算术：p^d = (p^k)^(d/k)（需要 Nat 算术引理）
3. aeval 计算（可能需要 simp 引理链）

## 潜在风险

- card(AdjoinRoot g) = p^k 的 Mathlib 路径可能需要多步中转
- `Fintype` 实例的自动合成可能有问题（AdjoinRoot 的 Fintype 依赖 powerBasis + ZMod 的 Fintype）
- `Fact (Irreducible g)` 需要手动构造：`haveI : Fact (Irreducible g) := ⟨hg⟩`
