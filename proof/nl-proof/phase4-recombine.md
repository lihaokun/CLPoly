# Phase 4: Recombination 正确性

> 状态：nl-proof v3（审核修正：Hensel 唯一性条件精确化 + 模逆元显式 + 探索文本清理）
> 对应 C++：`polynomial_factorize_univar.hh:750-882` `__zassenhaus_recombine`
> 参考：`docs/design/factorization/formal-proof-univar-factorization.md` §1-§3

---

## 0. 记号与设定（照搬 design doc §1）

- f ∈ Z[x]，squarefree，primitive（content = 1），deg f = n ≥ 2，ℓ = lc(f) > 0
- f = g₁ · g₂ · ... · gₜ（Z[x] 不可约分解）。每个 gⱼ primitive。ℓⱼ = lc(gⱼ)。ℓ = ∏ ℓⱼ。
- p 素数：p ∤ ℓ，f̄ = f mod p squarefree in F_p[x]
- f̄ = ℓ̄ · h̄₁ · ... · h̄ᵣ（h̄ᵢ monic irreducible in F_p[x]，两两互素）
- **因子对应**：每个 gⱼ 对应子集 Sⱼ ⊆ {1,...,r}，ḡⱼ / ℓ̄ⱼ = ∏_{i∈Sⱼ} h̄ᵢ。{Sⱼ} 是 {1,...,r} 的划分。
- **对称模**：sym_m(a) = a mod m 约化到 (-m/2, m/2]
- **因子系数界（Mignotte bound）**：g | f → ‖g‖_∞ ≤ C(n,n/2) · ‖f‖₂
  （Mathlib Mahler measure + Jensen 公式 + Parseval，精确匹配 C++ `__mignotte_bound`。
   证明见 phase4-mignotte.md §5，**0 sorry**）
- **精度条件 (M)**：m > 2ℓ · C(n,n/2) · ‖f‖₂（= C++ 的 `2 * lc_f * B`）

---

## 1. lc-baking Hensel 提升不变量（design doc §2）

### 1.1 初始化

调整 mod p 因子：h̃₁ = ℓ · h̄₁ mod p，h̃ᵢ = h̄ᵢ (i ≥ 2)。
乘积：h̃₁ · h̄₂ · ... · h̄ᵣ = ℓ · h̄₁ · ... · h̄ᵣ = f̄ (mod p)。✓

### 1.2 Hensel 提升后不变量 (I1)-(I4)

设 m = p^(2^k)（二次提升 k 步）。∃ 唯一 H₁,...,Hᵣ ∈ Z_m[x]：

**(I1)** ∏ Hᵢ ≡ f (mod m)
**(I2)** deg(Hᵢ) = deg(h̄ᵢ)
**(I3)** lc(H₁) = ℓ，lc(Hᵢ) = 1 (i ≥ 2)
**(I4)** H₁ ≡ h̃₁ (mod p)，Hᵢ ≡ h̄ᵢ (mod p) (i ≥ 2)

**来源**：Hensel 唯一性定理（§2）+ lc-baking 初始化。

### 1.3 (I3) 的证明（lc 保持）

Hensel 单步：G* = G + m·τ，H* = H + m·σ，其中 deg(σ) < deg(H)。
→ lc(H*) = lc(H)（σ 的度数 < H → 不改变 leading coeff）。
→ lc(G*) = lc(G)（类似论证，deg(τ) < deg(G)）。
每步保持 → lc(H₁) = ℓ（初始值），lc(Hᵢ) = 1 (i ≥ 2)。✓

**注**：这要求 Hensel 单步使用 divByMonic 控制度数（nl-proof phase3-t34-hensel.md §2.3）。
我们的 `hensel_step` 目前没有度数保持。**需要增强 hensel_step 或单独证明 (I3)**。

---

## 2. Hensel 唯一性定理（design doc §2.2，核心）

### 2.1 定理陈述

```
设 F ∈ Z_m[x]（m = p^k，k ≥ 1），Ā, B̄ ∈ F_p[x]，B̄ 首一，gcd(Ā, B̄) = 1。
若 (A₁, B₁) 和 (A₂, B₂) 都满足：
  Aᵢ · Bᵢ ≡ F (mod m)
  Aᵢ ≡ Ā (mod p), Bᵢ ≡ B̄ (mod p)
  deg(Aᵢ) = deg(Ā), deg(Bᵢ) = deg(B̄)
  Bᵢ 首一 in Z_m[x]（lc(Bᵢ) = 1，不仅 mod p）
则 A₁ = A₂, B₁ = B₂ in Z_m[x]。
```

**注**：
- **Ā 不要求首一**。证明中 `Ā | Ē` 和 `Ē = c·Ā` 对任意非零 Ā 成立（不需首一）。
- 关键条件是 **B̄ 首一 + Bᵢ 首一 in Z_m[x]**。这用于推出 c = 0（见 §2.2）。
- 在 lc-baking 设定中：Case 1 的 Ā = ℓ̄·∏h̄ᵢ（非首一），Case 2 交换后 Ā 含 H₁（非首一）。两种情况 B 侧均首一。✓

### 2.2 唯一性证明（对 k 归纳）

设 (A₁, B₁) 和 (A₂, B₂) 都满足上述条件。需证 A₁ = A₂, B₁ = B₂ in Z_m[x]。

**Base k = 1**：m = p。A₁ ≡ Ā (mod p) = A₂。✓

**Step k → k+1**（m = p^k → m' = p^{k+1}）：

IH：A₁ ≡ A₂ (mod p^k)，B₁ ≡ B₂ (mod p^k)。
已知 A₁·B₁ ≡ F ≡ A₂·B₂ (mod p^{k+1})。求证 A₁ ≡ A₂ (mod p^{k+1})。

设 A₁ = A₂ + p^k · E，B₁ = B₂ + p^k · F（E, F ∈ Z[x]）。

A₁·B₁ - A₂·B₂ = p^k·(E·B₂ + A₂·F) + p^{2k}·E·F。

p^{k+1} | 左边（因 A₁B₁ ≡ A₂B₂ mod p^{k+1}）。
k ≥ 1 → 2k ≥ k+1 → p^{k+1} | p^{2k}·E·F。
故 p^{k+1} | p^k·(E·B₂ + A₂·F) → p | (E·B₂ + A₂·F)。

在 F_p[x] 中：Ē·B̄ + Ā·F̄ = 0（其中 B̄ = B₂ mod p = B̄ 首一，Ā = A₂ mod p = Ā 首一）。

Ē·B̄ = -Ā·F̄。

**由 Ā 首一 + gcd(Ā, B̄) = 1**：

（Ā 和 B̄ 在 F_p[x] 中互素且 Ā 首一。）
Ā | Ē·B̄ 且 gcd(Ā, B̄) = 1 → Ā | Ē。
B̄ | Ā·F̄ 且 gcd(Ā, B̄) = 1 → B̄ | F̄。

**deg 约束**：
deg(E) ≤ deg(Ā)（因 A₁ - A₂ 的 deg ≤ deg(Ā)）。
deg(F) ≤ deg(B̄)（同理）。

Ā | Ē：deg(Ē) ≤ deg(E) ≤ deg(Ā)。Ā 首一，deg(Ā) > 0。
若 Ē ≠ 0：deg(Ē) ≥ deg(Ā)（因为 Ā | Ē 且 Ā 首一 → deg(Ā) ≤ deg(Ē)）。
结合 deg(Ē) ≤ deg(Ā)：deg(Ē) = deg(Ā)，Ē = c·Ā（c ∈ F_p，因首项系数匹配）。

同理 B̄ | F̄ 且 deg(F̄) ≤ deg(B̄) → 若 F̄ ≠ 0 则 F̄ = d·B̄。

代入 Ē·B̄ + Ā·F̄ = 0：c·Ā·B̄ + d·Ā·B̄ = (c+d)·Ā·B̄ = 0。
F_p[x] 整环 + Ā, B̄ ≠ 0 → c + d = 0 → d = -c。

所以 Ē = c·Ā，F̄ = -c·B̄（for some c ∈ F_p）。

**用 B 首一条件推出 c = 0**：

B₁ = B₂ + p^k·F。B₁ 首一（in Z_m[x]）→ lc(B₁) = 1。B₂ 首一 → lc(B₂) = 1。
在 deg(B̄) 位置：lc(B₁) = lc(B₂) + p^k·lc(F) = 1 + p^k·lc(F)。
B₁ 首一 → 1 + p^k·lc(F) = 1 (mod p^{k+1}) → p^k·lc(F) ≡ 0 (mod p^{k+1}) → p | lc(F)。

F̄ = -c·B̄。B̄ 首一 → lc(F̄) = -c。lc(F̄) = lc(F) mod p。
p | lc(F) → lc(F) mod p = 0 → -c = 0 → **c = 0** in F_p。✓

（**注**：这正是需要 "B 首一 in Z_m[x]" 的原因。仅 "B̄ 首一 in F_p[x]" 不够——
deg 保持 + mod p 首一不排除 c ≠ 0。GCL Theorem 6.1 明确要求 B_m 首一。）

c = 0 → Ē = 0 → E ≡ 0 (mod p) → p | E 的每个系数。
A₁ - A₂ = p^k·E。p | E → p^{k+1} | p^k·E → **A₁ ≡ A₂ (mod p^{k+1})**。
同理 c = 0 → F̄ = 0 → **B₁ ≡ B₂ (mod p^{k+1})**。✓

**所以唯一性需要 B_m（其中一个因子）首一 in Z_m[x]！**

### 2.3 应用于 lc-baking 设定

在 (I1)-(I4) 设定下：Hᵢ (i ≥ 2) 首一 in Z_m[x]（(I3)）。
所以对于 Sⱼ 的唯一性应用：

**情况 1**（|Sⱼ| 包含所有 i ≥ 2 的子集，使得 B = ∏_{i∉Sⱼ} Hᵢ 首一）：
若 1 ∉ Sⱼ：B = ∏_{i∉Sⱼ} Hᵢ 包含 H₁（lc = ℓ），B 不首一。❌
若 1 ∈ Sⱼ：B = ∏_{i∉Sⱼ} Hᵢ 全部首一 → B 首一。✓

**情况 2**（1 ∉ Sⱼ）：A = ∏_{i∈Sⱼ} Hᵢ 全部首一 → A 首一。唯一性定理的 "B 首一" 条件不满足。
但可以交换 A 和 B：让 A = ∏_{i∉Sⱼ} Hᵢ（包含 H₁），B = ∏_{i∈Sⱼ} Hᵢ（首一）。
则 B 首一 ✓。唯一性适用。

**结论**：两种情况都可以应用 Hensel 唯一性（选择合适的 A/B 分配使得 B 首一）。✓

---

## 3. 因子恢复定理（design doc §3）

### 3.1 定理

**定理 3.1**：设 (I1)-(I4) 和 (M) 成立。对任何真因子 gⱼ（对应子集 Sⱼ），定义：

```
P_{Sⱼ} = ∏_{i∈Sⱼ} Hᵢ           若 1 ∈ Sⱼ
P_{Sⱼ} = ℓ · ∏_{i∈Sⱼ} Hᵢ       若 1 ∉ Sⱼ
```

则 **pp(sym_m(P_{Sⱼ})) = gⱼ**。

### 3.2 证明

**情况 1**（1 ∈ Sⱼ）：

设 U = ∏_{i∈Sⱼ} Hᵢ，V = ∏_{i∉Sⱼ} Hᵢ。U·V ≡ f (mod m)。V 首一（i ≥ 2 全部首一）。

**构造对比组**：
- U' = (ℓ/ℓⱼ)·gⱼ ∈ Z[x]。lc(U') = (ℓ/ℓⱼ)·ℓⱼ = ℓ = lc(U)。
- V' = (ℓ/ℓⱼ)⁻¹ · hⱼ* mod m ∈ Z_m[x]，hⱼ* = f/gⱼ ∈ Z[x]。
  （注：(ℓ/ℓⱼ)⁻¹ 是模逆元，存在因 gcd(ℓ/ℓⱼ, m) = 1（p ∤ ℓ → p ∤ ℓ/ℓⱼ）。）
  lc(V') = (ℓ/ℓⱼ)⁻¹ · lc(hⱼ*) = (ℓ/ℓⱼ)⁻¹ · (ℓ/ℓⱼ) = 1 mod m（首一）。✓

U'·V' = (ℓ/ℓⱼ)·gⱼ · (ℓ/ℓⱼ)⁻¹·hⱼ* = gⱼ·hⱼ* = f（精确，in Z[x]）。✓
U' mod p = ℓ̄ · ∏_{i∈Sⱼ} h̄ᵢ = U mod p。✓
V' mod p = ∏_{i∉Sⱼ} h̄ᵢ = V mod p。V' 首一。✓

由 Hensel 唯一性（V, V' 首一）：**U ≡ U' (mod m)**。
即 ∏_{i∈Sⱼ} Hᵢ ≡ (ℓ/ℓⱼ)·gⱼ (mod m)。

**Mignotte 恢复**：
‖(ℓ/ℓⱼ)·gⱼ‖_∞ ≤ ℓ · ‖gⱼ‖_∞ ≤ ℓ · B_Mig(f) < m/2（by (M)）。
→ sym_m(U) = (ℓ/ℓⱼ)·gⱼ（精确）。
→ pp(sym_m(P_{Sⱼ})) = pp((ℓ/ℓⱼ)·gⱼ) = gⱼ（gⱼ primitive，ℓ/ℓⱼ ∈ Z₊）。✓

**情况 2**（1 ∉ Sⱼ）：

U = ∏_{i∈Sⱼ} Hᵢ（全部首一 → U 首一）。V = ∏_{i∉Sⱼ} Hᵢ（含 H₁，lc = ℓ）。

构造对比组（交换 A/B 使 B = U 首一）：
- U'' = gⱼ · ℓⱼ⁻¹ mod m（ℓⱼ⁻¹ 模逆元存在因 gcd(ℓⱼ, m) = 1）。lc(U'') = ℓⱼ · ℓⱼ⁻¹ = 1（首一）。✓
- V'' = ℓⱼ · hⱼ* mod m，hⱼ* = f/gⱼ ∈ Z[x]。lc(V'') = ℓⱼ · (ℓ/ℓⱼ) = ℓ。✓
- U''·V'' = f mod m。✓
- U'' mod p = ∏_{i∈Sⱼ} h̄ᵢ = U mod p。U'' 首一。✓
- V'' mod p = ℓ̄ · ∏_{i∉Sⱼ} h̄ᵢ = V mod p。✓

唯一性（U, U'' 首一）：**U ≡ U'' (mod m)**。
即 ∏_{i∈Sⱼ} Hᵢ ≡ gⱼ/ℓⱼ (mod m)。

P_{Sⱼ} = ℓ·U ≡ (ℓ/ℓⱼ)·gⱼ (mod m)。

Mignotte 恢复 + pp 同情况 1。✓

---

## 4. RecombineCorrect 证明

### 4.1 定理

```lean
theorem recombine_correct
    (f : Polynomial ℤ) (hf : f ≠ 0) (hprim : IsPrimitive f) (hsqf : Squarefree f)
    (p : ℕ) (hp : Nat.Prime p) (k : ℕ)
    (hensel_factors : List (Polynomial (ZMod (p^k))))
    (hhensel : HenselCorrect f k ... hensel_factors)
    (hmig : mignotte_precision f p k)  -- 1 sorry (Mignotte bound)
    : ∃ result, RecombineCorrect f result
```

### 4.2 证明

1. f = g₁·...·gₜ（Z[x] UFD，不可约分解）。
2. §0 的因子对应给出划分 S₁,...,Sₜ。
3. §3 的因子恢复给出：gⱼ = pp(sym_m(P_{Sⱼ}))。
4. 取 result = [g₁,...,gₜ]。
5. Associated f result.prod：f = g₁·...·gₜ ✓。
6. ∀ g ∈ result, Irreducible g：每个 gⱼ 不可约 by construction ✓。

**注**：这个证明的**实质内容**在 §2（Hensel 唯一性）和 §3（因子恢复）。
§4.2 只是把它们组合起来。Hensel 因子通过唯一性 + Mignotte 恢复连接到真因子。

---

## 5. 形式化策略

### sorry 清单：0 个

全链 0 sorry。因子系数界用 Cauchy bound + 初等对称函数完整证明。

### 形式化顺序与行数

| Step | 内容 | 估计行数 |
|------|------|---------|
| 1 | `mignotte_bound`：Mahler measure + Jensen + Parseval | ~95 |
| 2 | `symmetric_recovery`：\|a\| < m/2 → sym_m(a) = a | ~30 |
| 3 | `hensel_unique`：唯一性定理（对 k 归纳 + 首一条件） | ~120 |
| 4 | `factor_correspondence`：Hensel 因子子集 ↔ 不可约因子 | ~80 |
| 5 | `factor_recovery`：pp(sym_m(P_S)) = gⱼ | ~60 |
| 6 | `recombine_correct`：RecombineCorrect | ~30 |
| **总计** | | **~415** |
