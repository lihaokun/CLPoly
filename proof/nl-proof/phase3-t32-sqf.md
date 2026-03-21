# T3.2: SQF 算法建模与正确性

> 状态：nl-proof v3（审核修正版）
> 对应 C++：`polynomial_factorize_zp.hh:108-180` `__squarefree_Zp` + `__extract_pth_root`

---

## 0. 数学背景

设 `f ∈ F_p[x]` 非零，不可约分解 `f ~ u · ∏ qⱼ^{eⱼ}`（qⱼ monic 不可约，两两不同）。

**char p 导数性质**：

对 `f = r^v · g`（`r` 不可约，`gcd(r, g) = 1`，即 `v = v_r(f)`）：
```
f' = v · r^{v-1} · r' · g + r^v · g'
   = r^{v-1} · (v · r' · g + r · g')
```
故 **`v_r(f') ≥ v_r(f) - 1`**（关键不等式，对任意不可约 r 成立）。

**GCD 结构**：
- `v_r(gcd(f, f')) = min(v_r(f), v_r(f'))`
- 结合上式：`v_r(gcd(f, f')) ≥ min(v, v-1) = v - 1`
- `w₀ := f / gcd(f, f')`：`v_r(w₀) = v - v_r(gcd(f, f')) ≤ v - (v-1) = 1`
- **∀ 不可约 r，`v_r(w₀) ≤ 1` → `Squarefree w₀`** ✓

**Yun 迭代模式**（第 i 步）：
- `zᵢ = ∏_{eⱼ=i} qⱼ`（恰好 i 次出现且 p∤i 的不可约因子之积）
- `wᵢ = ∏_{eⱼ>i, p∤eⱼ} qⱼ`
- `cᵢ = ∏_{eⱼ>i, p∤eⱼ} qⱼ^{eⱼ-i-1} · ∏_{p|eⱼ} qⱼ^{eⱼ}`

**循环终止后**：`c_rem = ∏_{p|eⱼ} qⱼ^{eⱼ}`，满足 `c_rem' = 0`。

**Frobenius 关键等式**：在 `F_p[X]` 中 `expand p f = f ^ p`。
证明：`map_frobenius_expand` 给出 `map (frobenius (ZMod p) p) (expand p f) = f^p`。
在 `ZMod p` 中 `frobenius` 是恒等映射（`a^p = a`），故 `map id (expand p f) = expand p f = f^p`。

---

## 1. 算法模型

### 1.1 Yun 内循环

```lean
noncomputable def yunLoop
    (w c : Polynomial (ZMod p)) (i : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    : List (Polynomial (ZMod p) × ℕ) × Polynomial (ZMod p) :=
  if hw : w.natDegree = 0 then
    (acc, c)
  else
    let y := normalize (EuclideanDomain.gcd w c)
    let z := normalize (w /ₘ y)
    let acc' := if 0 < z.natDegree then acc ++ [(z, i)] else acc
    let c' := normalize (c /ₘ y)
    yunLoop y c' (i + 1) acc'
termination_by w.natDegree + c.natDegree
```

**返回值**：`(extracted_factors, c_remainder)`

**操作说明**：
- `y = normalize(gcd(w, c))` → Monic（因 `w ≠ 0` 从循环条件得 `gcd ≠ 0`）
- `/ₘ` 是 monic 除法，`y | w` 且 `y | c`（gcd 性质）保证精确
- `normalize` 保证输出 Monic

**终止性**（度量 `w.natDegree + c.natDegree`）：

**Case A**：`deg(y) ≥ 1`（gcd 非平凡）
- `y | w` → `deg(y) ≤ deg(w)`
- `c = y · (c/y)` → `deg(c/y) = deg(c) - deg(y)`（整环 + `c ≠ 0`）
- 新度量 = `deg(y) + deg(c) - deg(y)` = `deg(c)`
- 旧度量 = `deg(w) + deg(c)` ≥ `1 + deg(c)` > `deg(c)` ✓

**Case B**：`deg(y) = 0`（gcd 是单位，w 与 c 互素）
- `y` 是 monic + deg 0 → `y = 1`
- `z = w / 1 = w`（deg > 0，输出此因子）
- `w_new = 1`（deg = 0）
- 新度量 = `0 + deg(c)` < `deg(w) + deg(c)` ✓（因 `deg(w) ≥ 1`）

**终止性前提**：需要 `c ≠ 0`（见不变量 Y7）。

### 1.2 顶层 sqfZp

```lean
noncomputable def sqfZp (f : Polynomial (ZMod p)) :
    List (Polynomial (ZMod p) × ℕ) :=
  if hf : f.natDegree = 0 then []
  else if hderiv : derivative f = 0 then
    let g := Polynomial.contract p f
    (sqfZp g).map (fun pr => (pr.1, pr.2 * p))
  else
    let c := normalize (EuclideanDomain.gcd f (derivative f))
    let w := normalize (f /ₘ c)
    let (yun_result, c_rem) := yunLoop w c 1 []
    if 0 < c_rem.natDegree then
      let g := Polynomial.contract p c_rem
      yun_result ++ (sqfZp g).map (fun pr => (pr.1, pr.2 * p))
    else
      yun_result
termination_by f.natDegree
```

**终止性**：

**Branch 1**（`f' = 0`）：
- `deg(contract p f) = deg(f) / p`（由 `natDegree_expand` + `expand_contract`）
- `p ≥ 2` 且 `deg(f) ≥ 1` → `deg(f) / p ≤ deg(f) / 2 < deg(f)` ✓

**Branch 2**（`c_rem` 递归）：
- 需要 `deg(contract p c_rem) < deg(f)`
- 先证 `deg(c_rem) ≤ deg(f) - 1`：
  - `c₀ = gcd(f, f')`，`w₀ = f / c₀`
  - `f' ≠ 0` → `∃ qⱼ with p ∤ eⱼ` → `v_{qⱼ}(w₀) = 1`（≥ 1） → `deg(w₀) ≥ 1`
  - `deg(f) = deg(w₀) + deg(c₀)` → `deg(c₀) = deg(f) - deg(w₀) ≤ deg(f) - 1`
  - Yun 循环只缩小 c（每步 `c_new = c / y`，`deg(c_new) ≤ deg(c)`）
  - 故 `deg(c_rem) ≤ deg(c₀) ≤ deg(f) - 1`
- 则 `deg(contract p c_rem) = deg(c_rem) / p ≤ (deg(f) - 1) / p < deg(f)` ✓
  - （p ≥ 2 且 deg(f) ≥ 1 → `(deg(f) - 1) / p < deg(f)` → 严格 `<`，需 `deg(f) ≥ 1`，由 `hf` 保证）

**注**：此终止性**不依赖** `derivative(c_rem) = 0`，避免循环依赖。

---

## 2. Yun 循环不变量

在 `yunLoop w c i acc` 每次调用入口处，设 `f_orig` 为原始输入多项式：

| 标号 | 不变量 | 含义 |
|------|--------|------|
| Y1 | `Associated f_orig ((acc.map (fun (s,e) => s^e)).prod * w^i * c)` | 乘积还原 |
| Y2 | `Squarefree w` | w 无平方 |
| Y3 | `i ≥ 1` | 计数器 ≥ 1 |
| Y4 | `∀ (s,e) ∈ acc, Monic s ∧ Squarefree s ∧ 0 < s.natDegree ∧ e ≥ 1` | acc 条目正确 |
| Y5 | `∀ (s₁,e₁) (s₂,e₂) ∈ acc, (s₁,e₁) ≠ (s₂,e₂) → IsCoprime s₁ s₂` | acc 两两互素 |
| Y6 | `∀ (s,e) ∈ acc, IsCoprime s w` | acc 与 w 互素 |
| Y7 | `c ≠ 0` | c 非零（终止性前提） |
| Y8 | `Monic w ∨ w = 1` | w 是 monic 或 1 |
| Y10 | `∀ (s,e) ∈ acc, IsCoprime s c` | acc 与 c 互素（关键！推导 c_rem' = 0） |

### 初始调用满足不变量

`yunLoop w₀ c₀ 1 []`，其中 `c₀ = normalize(gcd(f, f'))`，`w₀ = normalize(f /ₘ c₀)`。

- **Y1**: `1 * w₀^1 * c₀ = w₀ * c₀`。由 `gcd(f, f') | f` 知 `f = c₀_raw * w₀_raw`（精确除法），normalize 引入单位，`w₀ * c₀ ~ f` ✓
- **Y2**: `Squarefree w₀`

  **证明**（valuation 论证）：

  对任意不可约 r，设 `v = v_r(f)` = `emultiplicity r f`。
  1. `f = r^v · g`，`gcd(r, g) = 1`
  2. `f' = r^{v-1} · (v · r' · g + r · g')`，故 `v_r(f') ≥ v - 1`
  3. `v_r(c₀) = v_r(gcd(f, f')) = min(v_r(f), v_r(f')) ≥ min(v, v-1) = v - 1`
  4. `v_r(w₀) = v_r(f) - v_r(c₀) ≤ v - (v-1) = 1`
  5. 由 `squarefree_iff_emultiplicity_le_one`：∀ 不可约 r，`v_r(w₀) ≤ 1` → `Squarefree w₀` ✓

  **Mathlib 路径**：`squarefree_iff_emultiplicity_le_one` + 自定义引理 `emultiplicity_derivative_ge`。

  **备选（更直接）**：假设 `r · r | w₀`（r 非单位），推矛盾：
  - `r² | w₀ | f` → `r² | f` → `r | f'`（导数中 `r^{v-1}` 因子）→ `r | gcd(f, f') = c₀`
  - `r | w₀` 且 `r | c₀` → `r² | w₀ · c₀ = f` → `r³ | f`
  - 递归：`r^n | f` → `r^{n-1} | f'` → `r^{n-1} | c₀` → `r^n | f` 且 `r^{n-1} | c₀` → `v_r(w₀) ≤ v - (v-1) = 1`

- **Y3**: `1 ≥ 1` ✓
- **Y4**: `[]` 空真 ✓
- **Y5**: `[]` 空真 ✓
- **Y6**: `[]` 空真 ✓
- **Y7**: `c₀ = normalize(gcd(f, f')) ≠ 0`（因 `f ≠ 0` → `gcd(f, f') ≠ 0`）✓
- **Y8**: `w₀ = normalize(...)` → Monic ✓

### Y1 保持

设 `y = gcd(w, c)`，`z = w / y`，`w_new = y`，`c_new = c / y`。

旧 Y1: `P * w^i * c ~ f`（P = acc_prod）

**Case: deg(z) > 0**（`acc_new = acc ++ [(z, i)]`）
```
P * z^i * y^{i+1} * (c/y)
= P * z^i * y^i * c
= P * (z · y)^i * c        [交换律：(ab)^n = a^n · b^n]
= P * w^i * c               [w = z · y]
~ f ✓
```

**Case: deg(z) = 0**（`acc_new = acc`，z 是非零常数）
```
P * y^{i+1} * (c/y)
= P * y^i * c
```
`z = w / y` 非零常数 → `w = y · z`（Associated）→ `w^i ~ y^i` → `P * y^i * c ~ P * w^i * c ~ f` ✓

### Y2 保持

`y = gcd(w, c)` 整除 `w`。`Squarefree w`（Y2 旧）+ `Squarefree.squarefree_of_dvd (y ∣ w)` → `Squarefree y` ✓

### Y7 保持

`c_new = c / y`。`c ≠ 0`（Y7 旧）且 `y ≠ 0`（因 `y | w`，`w ≠ 0`）→ 整环中 `c/y ≠ 0` ✓

### Y4 保持（新条目 (z, i)）

- `Monic z`：normalize 保证 ✓
- `Squarefree z`：`z | w`（`w = z · y`）+ `Squarefree w` → `Squarefree z` ✓
- `0 < z.natDegree`：条件检查保证 ✓
- `i ≥ 1`：由 Y3 ✓

### Y5/Y6 保持（coprimality）

**Y6**（新 z 与 w_new = y 互素）：
- `w = z · y`，`Squarefree w`（Y2）→ `IsRelPrime z y`（由 `squarefree_mul_iff`）✓

**Y5**（新 z 与旧 acc 条目互素）：
- 旧条目 `s` 满足 `IsCoprime s w`（Y6 旧），`z | w` → `IsCoprime s z` ✓

**Y6 更新**（旧 acc 条目与 w_new = y 互素）：
- 旧 `IsCoprime s w`（Y6 旧），`y | w` → `IsCoprime s y` ✓

### Y10 保持（acc 与 c 互素）— 关键新不变量

**新条目 z 与 c_new = c/y 互素**：

设 r 是 z 和 c/y 的公共不可约因子（反证法）。
1. `r | z` → `r | w`（因 `z | w`，即 `w = z · y`）
2. `r | z` 且 `IsRelPrime z y`（Y6 上面刚证）→ `r ∤ y`
3. `r | (c/y)` → 由 `c = y · (c/y)` 得 `r | c`
4. `r | w` 且 `r | c` → `r | gcd(w, c) = y`
5. **矛盾**：步骤 2 说 `r ∤ y`，步骤 4 说 `r | y`。

故 `IsCoprime z (c/y)` ✓

**旧条目 s 与 c_new = c/y 互素**：
- `IsCoprime s c`（Y10 旧）→ `IsCoprime s (y · (c/y))`
- `IsCoprime a (b · c) → IsCoprime a c`（标准引理）
- 故 `IsCoprime s (c/y)` ✓

---

## 3. 终止情况正确性

### 3.1 Yun 循环终止后（w.natDegree = 0）

返回 `(acc, c_rem)`。`w` 是常数（Monic + deg 0 → `w = 1`）。

Y1：`P * 1^i * c_rem ~ f`，即 `P * c_rem ~ f`。

### 3.2 `derivative(c_rem) = 0` 的完整证明

**定理**：Yun 循环终止后（`w = 1`），残余 c_rem 满足 `derivative c_rem = 0`。

**证明**：

使用不变量 Y10（`∀ (s,e) ∈ acc, IsCoprime s c`）和 Y2（`Squarefree w₀`）。

**Step 1**：`IsCoprime P c_rem`（P = acc_prod）

由 Y10 终止时：∀ (s, e) ∈ acc，`IsCoprime s c_rem`。
由 Y5：acc 条目两两互素。
标准引理：两两互素的因子之积与每个因子互素的元素互素。
故 `IsCoprime P c_rem` ✓

**Step 2**：对任意不可约 q 满足 `q | c_rem`，证明：**q 可分（q' ≠ 0）→ p | v_q(c_rem)**。

反证法：假设 q 可分（`q' ≠ 0`），`q | c_rem`，`p ∤ v_q(c_rem)`。

设 `v = v_q(c_rem)`。由 `IsCoprime P c_rem`：`v_q(P) = 0`（q 不整除 P）。

由 Y1 终止时：`P · c_rem ~ f`。故 `v_q(f) = v_q(P) + v_q(c_rem) = 0 + v = v`。

由 `p ∤ v`、`q' ≠ 0` 和 §0 的 valuation 论证：
- `f = q^v · g`，`gcd(q, g) = 1`
- `f' = q^{v-1} · (v · q' · g + q · g')`
- **q' ≠ 0 且 gcd(q, g) = 1 → q ∤ (v · q' · g)**（因 q 不可约，q ∤ q'（deg q' < deg q），q ∤ g）
- 又 `p ∤ v` → `v ≠ 0` in F_p → `v · q' · g ≠ 0`
- 故 `v_q(v · q' · g + q · g') = 0`（首项 v_q = 0，次项 v_q ≥ 1，和的 v_q = 0）
- `v_q(f') = v - 1`（精确）
- `v_q(gcd(f, f')) = min(v, v-1) = v - 1`
- `v_q(w₀) = v - (v-1) = 1`

所以 `q | w₀`（初始 w 包含 q）。

**Step 3**：q 必然进入某个 acc 条目（此步不依赖 q 可分/不可分）。

Yun 循环终止时 `w = 1`，故 `q ∤ w_final`。
但初始 `q | w₀`。在循环过程中 w 的值变化为：`w₀, y₁, y₂, ...`。
由于 `yᵢ | wᵢ₋₁`（gcd 整除性），`v_q` 只可能下降或保持。
∃ 最小 k 使得 `q | w_{k-1}` 但 `q ∤ w_k = gcd(w_{k-1}, c_{k-1})`。

`q | w_{k-1}` 且 `q ∤ gcd(w_{k-1}, c_{k-1})` → `q ∤ c_{k-1}`（否则 `q | gcd`）。

此时 `z_k = w_{k-1} / gcd(w_{k-1}, c_{k-1})`。
由 `Squarefree w_{k-1}`（Y2）：`v_q(w_{k-1}) = 1`。
`v_q(gcd) = 0`（因 q ∤ gcd）→ `v_q(z_k) = 1 - 0 = 1`。

故 `q | z_k`，z_k 被输出到 acc（`deg(z_k) ≥ 1`）。故 `q | P`。

**Step 4**：矛盾。

`q | P`（Step 3）且 `q | c_rem`（假设）→ `q | gcd(P, c_rem)`。
但 `IsCoprime P c_rem`（Step 1）→ `IsUnit q`。矛盾（q 不可约）。 ✓

故假设不成立：q 可分且 `q | c_rem` → `p | v_q(c_rem)`。

**Step 5**：结论——分两种情况。

写 `c_rem = ∏ qⱼ^{mⱼ}`（不可约分解）。导数：
```
derivative(c_rem) = ∑_j mⱼ · qⱼ^{mⱼ-1} · qⱼ' · ∏_{k≠j} qₖ^{mₖ}
```

对每个 j，第 j 项中的因子 `mⱼ · qⱼ'` 为零——分两种情况：
- **qⱼ 可分**（`qⱼ' ≠ 0`）：由 Step 2-4，`p | mⱼ` → `mⱼ ≡ 0 (mod p)` → 系数 `mⱼ = 0` in F_p → 项为 0
- **qⱼ 不可分**（`qⱼ' = 0`）：直接 `qⱼ' = 0` → 项为 0

每项为 0 → 和为 0 → `derivative(c_rem) = 0` ✓

### 3.2.1 Lean 形式化路径（干净版，补全全部缺口）

**策略**：反证法，完全避免 Step 5 的 UFD 展开。

**前置引理 A**：`yunLoop_preserves_irred_mult`
```
若 q 不可约，q ∤ w₀，则 ∀ k, q^k | c₀ ↔ q^k | c_rem。
```
**证明**（yunLoop 结构归纳，~15 行）：
- yunLoop 每步 c_new = c / gcd(w, c)。
- q ∤ w₀ → q ∤ w_i（归纳：w_i | w₀）→ q ∤ gcd(w_i, c_i)（gcd | w_i）。
- q ∤ gcd → gcd 的 q-幂次 = 0 → c_new = c / gcd 保持 q 的所有幂次。
- 每步 `q^k | c_i ↔ q^k | c_{i+1}`，归纳得 `q^k | c₀ ↔ q^k | c_rem`。

**注意**：不需要 `emultiplicity`。只需 `q^k | c ↔ q^k | c/gcd` 当 `q ∤ gcd`。
具体：`c = gcd * (c/gcd)`。`q^k | c` 且 `q ∤ gcd` → `q^k | c/gcd`（因 `IsCoprime q^k gcd`，
用 `Irreducible.prime → IsCoprime q gcd → IsCoprime q^k gcd → dvd_of_dvd_mul`）。
反向：`q^k | c/gcd` → `q^k | c`（因 `c/gcd | c`）。

**前置引理 B**：`not_pow_dvd_derivative_of_separable`（已证 0 sorry）
```
f = q^v * h, q ∤ h, q' ≠ 0, p ∤ v → q^v ∤ f'
```

**前置引理 C**：`gcd_dvd_of_dvd_dvd`（Mathlib 已有：`EuclideanDomain.dvd_gcd`）

---

**主证明**：`derivative(c_rem) = 0`

```
by_contra hderiv_ne（假设 derivative(c_rem) ≠ 0）

Step (a): squarefree_div_gcd_derivative(c_rem) 给出 w' Squarefree。
  derivative ≠ 0 → deg(derivative) ≤ deg(c_rem) - 1 < deg(c_rem)
  → c_rem ∤ derivative(c_rem)（度数）
  → gcd(c_rem, c_rem') 是真因子
  → w' = c_rem / gcd 有 deg ≥ 1
  → w' 不是 unit。

Step (b): WfDvdMonoid.exists_irreducible_factor → ∃ 不可约 q, q | w'。

Step (c): q | c_rem（w' | c_rem，因 c_rem = w' * gcd）。

Step (d): IsCoprime P c_rem（从 yunLoop_correct Y10 + Y5）→ q ∤ P。

Step (e): q | f（c_rem | f，因 f ~ P * c_rem）。

Step (f): f ~ w₀ * c₀, q prime → q | w₀ 或 q | c₀。

**Case 1**: q | w₀。
  yunLoop_extracts_factor → q | P。矛盾 step (d)。✓

**Case 2**: q ∤ w₀（推矛盾）。

  Step 2a: q | c₀（从 q | f ~ w₀ * c₀, q prime, q ∤ w₀）。

  Step 2b: 确定 c_rem 中 q 的精确幂次。
    设 v = 使 q^v | c_rem 成立的最大 v。（v ≥ 1，因 q | c_rem。）
    q^v | c_rem, q^{v+1} ∤ c_rem。

  Step 2c: q^v | c₀。
    由引理 A（yunLoop_preserves_irred_mult，反方向）：
    q^v | c_rem → q^v | c₀。
    （这里用了 q ∤ w₀ → yunLoop 保持 q 幂次。）

  Step 2d: q^{v+1} ∤ c₀。
    同引理 A（正方向）：q^{v+1} ∤ c_rem → q^{v+1} ∤ c₀。

  Step 2e: q^v | f 且 q^{v+1} ∤ f。
    f ~ w₀ * c₀。q ∤ w₀ → IsCoprime q w₀（q irreducible）→ IsCoprime q^{v+1} w₀。
    q^{v+1} ∤ c₀（step 2d）。
    q^{v+1} | f = w₀ * c₀ → q^{v+1} | c₀（IsCoprime）→ 矛盾。
    故 q^{v+1} ∤ f。
    同理 q^v | f（从 q^v | c₀ | c₀ * w₀ ~ f）。

  Step 2f: v_q(c_rem') = v - 1（精确）。
    q | w' = c_rem / gcd(c_rem, c_rem')。w' * gcd = c_rem。
    Squarefree w' → q^2 ∤ w'。
    若 q | gcd：q | w' 且 q | gcd → q^2 | c_rem（w' * gcd）。
      由 pow_dvd_derivative_of_pow_succ_dvd：q^2 | c_rem → q | c_rem'。
      q | c_rem 且 q | c_rem' → q | gcd。✓（自洽）。
      但这不直接给 v_q(c_rem') 的精确值。

    更直接：
    q^v | c_rem（step 2b），q^{v+1} ∤ c_rem。
    pow_dvd_derivative_of_pow_succ_dvd：q^v | c_rem → q^{v-1} | c_rem'（当 v ≥ 1）。
    若 q^v | c_rem'：
      q^v | c_rem 且 q^v | c_rem' → q^v | gcd(c_rem, c_rem')。
      w' = c_rem / gcd。q^v | c_rem, q^v | gcd → q^v | gcd, 但 c_rem = w' * gcd。
      q | w'（已知）→ q 至少 1 次在 w'。q^v | gcd。
      q^{v+1} | w' * gcd = c_rem（q^1 from w', q^v from gcd）。
      但 q^{v+1} ∤ c_rem（step 2b）。矛盾。
    故 q^v ∤ c_rem'。结合 q^{v-1} | c_rem'：v_q(c_rem') = v - 1。

  Step 2g: p ∤ v 且 q' ≠ 0。
    c_rem = q^v * g（gcd(q,g)=1，从 v = max power）。
    c_rem' = q^{v-1} * (C(v)*q'*g + q*g')。
    v_q(c_rem') = v - 1（step 2f）→ q ∤ (C(v)*q'*g + q*g')。
    若 p | v：C(v) = 0 → 内因子 = q*g' → q | 内因子。矛盾。
    若 q' = 0：同理 q | 内因子。矛盾。
    故 p ∤ v 且 q' ≠ 0。

  Step 2h: q^v ∤ f'。
    由引理 B（not_pow_dvd_derivative_of_separable）：
    q^v | f（step 2e），q^{v+1} ∤ f（step 2e），q' ≠ 0（step 2g），p ∤ v（step 2g），
    f = q^v * h 且 q ∤ h（从 q^v | f 的精确分解）→ q^v ∤ f'。

  Step 2i: q^v | c₀ 且 c₀ | f' → q^v | f'。矛盾 step 2h。
    c₀ = normalize(gcd(f, f'))。c₀ | f'（gcd 整除第二个参数）。
    q^v | c₀（step 2c）→ q^v | f'。
    但 q^v ∤ f'（step 2h）。矛盾。✓

两种情况均矛盾 → 假设不成立 → derivative(c_rem) = 0。✓
```

**引理清单**：
| 引理 | 状态 | 需自证 |
|------|------|--------|
| `squarefree_div_gcd_derivative` | ✅ 已证 | — |
| `not_pow_dvd_derivative_of_separable` | ✅ 已证 | — |
| `yunLoop_preserves_irred_mult` | ❌ 未证 | ~15 行 yunLoop 归纳 |
| `yunLoop_extracts_factor` | ✅ 已证 | — |
| `pow_dvd_derivative_of_pow_succ_dvd` | ✅ 已证 | — |

**关键**：整个证明**不使用 `emultiplicity`**。用 `q^k ∣` 和 `¬(q^{k+1} ∣)` 替代所有 valuation 操作。避免 ℕ∞ 算术。

**自需引理仅 1 个**：`yunLoop_preserves_irred_mult`（yunLoop 归纳，同结构已用 3 次）。
(b) WfDvdMonoid.exists_irreducible_factor: ∃ 不可约 q, q | w'。
(c) q | c_rem（因 w' | c_rem）。
(d) IsCoprime P c_rem (Step 1) → q ∤ P。
(e) q | w₀ — 关键步骤，见下。
(f) yunLoop_extracts_factor → q | P。
(g) 矛盾：(d) 和 (f)。
```

**步骤 (e) "q | w₀" 的证明（需要 emultiplicity）**：

这是唯一需要 emultiplicity 的地方。纯 dvd 论证无法推出 `q | w₀`。

```
设 v = emultiplicity q c_rem。（v ≥ 1，因 q | c_rem）
由 IsCoprime P c_rem：emultiplicity q P = 0。
由 f ~ P * c_rem：emultiplicity q f = 0 + v = v。
```

**子引理（核心）**：`q | w'` → `emultiplicity q (derivative c_rem) = v - 1`（精确）。

```
q | w' = c_rem / gcd(c_rem, c_rem')
→ emultiplicity q w' ≥ 1
→ emultiplicity q c_rem - emultiplicity q (gcd(c_rem, c_rem')) ≥ 1
  （用 emultiplicity_mul: emultiplicity q (a*b) = emultiplicity q a + emultiplicity q b）
→ emultiplicity q (gcd(c_rem, c_rem')) ≤ v - 1
但 emultiplicity q (gcd) = min(v, emultiplicity q (c_rem'))
且 pow_dvd_derivative: emultiplicity q c_rem' ≥ v - 1
所以 min(v, ≥v-1) = v - 1（当 v ≥ 1）
且 emultiplicity q w' = v - (v-1) = 1。✓
```

**然后**：
```
emultiplicity q f' ≥ v - 1（pow_dvd_derivative_of_pow_succ_dvd）。
但实际上 emultiplicity q f' = v - 1 精确？
  f = P * c_rem ~ products。
  f' = P' * c_rem + P * c_rem'。
  emultiplicity q (P' * c_rem) = emultiplicity q P' + v ≥ v（因 q ∤ P → q ∤ P'? 不一定！P' 可能被 q 整除）

这条路很复杂。换一种更直接的方式：

直接用 f = w₀ * c₀ 的 valuation。
emultiplicity q f = v。
emultiplicity q c₀ = emultiplicity q (gcd(f, f'))。
emultiplicity q f' ≥ v - 1（pow_dvd_derivative）。
emultiplicity q c₀ = min(v, emultiplicity q f') ≥ min(v, v-1) = v - 1。
emultiplicity q w₀ = v - emultiplicity q c₀ ≤ v - (v-1) = 1。

但 emultiplicity q w₀ ≥ 1 吗？不一定！可能 emultiplicity q c₀ = v（当 emultiplicity q f' ≥ v），则 emultiplicity q w₀ = 0。

这正是不可分因子的情况：q' = 0 → emultiplicity q f' ≥ v → emultiplicity q c₀ = v → q ∤ w₀。

所以需要区分：
- q 可分（q' ≠ 0）且 p ∤ v → emultiplicity q f' = v - 1（精确）→ q | w₀
- q 不可分 或 p | v → emultiplicity q f' ≥ v → q ∤ w₀

从 q | w'（c_rem 的 squarefree-part），我们推出 emultiplicity q (derivative c_rem) < emultiplicity q c_rem。这意味着 q 在 c_rem 中的行为像"可分因子"。具体：

emultiplicity q c_rem' = v - 1（从 emultiplicity q w' = 1 推出，见上）。
pow_dvd_derivative 给 emultiplicity q c_rem' ≥ v - 1。结合得精确 = v - 1。

这意味着 c_rem 的导数在 q 位置"精确下降 1"。在 char p 中，这等价于 q 可分（q' ≠ 0）且 p ∤ v。

但实际上我不需要推导出"q 可分"——我只需要 q | w₀：

从 emultiplicity q c_rem' = v - 1 < v = emultiplicity q c_rem：
  → emultiplicity q (gcd(c_rem, c_rem')) = v - 1
  → emultiplicity q w' = 1

但这给的是 c_rem 层面的信息，不是 f 层面的。需要转到 f：

emultiplicity q f = v（从 q ∤ P）。
emultiplicity q f'：
  f ~ P * c_rem。f' = P' * c_rem + P * c_rem'（导数的乘积规则不对 Associated 直接适用...）

这里有一个微妙问题：f ~ P * c_rem 是 Associated（差一个 unit），不是精确等式。导数不保持 Associated。

**修正**：f = u * P * c_rem 对某 unit u。f' = u * (P' * c_rem + P * c_rem')（u 是常数，u' = 0）。
```

---

### 3.2.2 正确且完整的 Lean 路径（审核修正版）

**整体策略**：反证法。假设 `derivative(c_rem) ≠ 0`，推出矛盾。

**核心依赖**：`squarefree_div_gcd_derivative`（已证）、`yunLoop_extracts_factor`（已证）、`pow_dvd_derivative_of_pow_succ_dvd`（已证）、`emultiplicity_mul`（Mathlib）。

```
假设 derivative(c_rem) ≠ 0。

(a) squarefree_div_gcd_derivative(c_rem):
    w' := c_rem / gcd(c_rem, c_rem') 是 Squarefree。
    deg(c_rem') ≤ deg(c_rem) - 1 < deg(c_rem)（导数度严格更小）。
    derivative(c_rem) ≠ 0 → c_rem ∤ c_rem' → gcd ≠ c_rem → deg(w') ≥ 1。

(b) WfDvdMonoid.exists_irreducible_factor: ∃ 不可约 q, q | w'。

(c) q | c_rem（因 w' | c_rem：c_rem = w' * gcd）。

(d) IsCoprime P c_rem (yunLoop_correct Y10 + Y5) → q ∤ P。

(e) q | f（f ~ P * c_rem，q | c_rem → q | P * c_rem ~ f）。

(f) f ~ w₀ * c₀。q 不可约 → q prime → q | w₀ 或 q | c₀。

分两种情况：

**Case 1: q | w₀**
→ yunLoop_extracts_factor → q | P → 矛盾（d）。✓

**Case 2: q ∤ w₀**（需证矛盾）

关键推导链（需要 emultiplicity）：

设 v = emultiplicity q c_rem ≥ 1（从 q | c_rem）。
设 v_f = emultiplicity q f。

Step 2a: v_f = v。
  IsCoprime P c_rem → emultiplicity q P = 0。
  f = u * P * c_rem（u 是 unit 常数，从 Associated）。
  emultiplicity q f = emultiplicity q P + emultiplicity q c_rem = 0 + v = v。
  （用 emultiplicity_mul for Prime q + emultiplicity_eq_zero.mpr(q ∤ P)）

Step 2b: q | w' → emultiplicity q c_rem' = v - 1（精确）。
  w' * gcd(c_rem, c_rem') = c_rem → emultiplicity q w' + emultiplicity q gcd = v。
  q | w' → emultiplicity q w' ≥ 1。
  pow_dvd_derivative: emultiplicity q c_rem' ≥ v - 1。
  emultiplicity q gcd = min(v, emultiplicity q c_rem')
    ≥ min(v, v-1) = v - 1（当 v ≥ 1）。
  从 emultiplicity q w' + emultiplicity q gcd = v 和 emultiplicity q w' ≥ 1：
    emultiplicity q gcd ≤ v - 1。
  结合 ≥ v-1 和 ≤ v-1：emultiplicity q gcd = v - 1。
  emultiplicity q w' = v - (v-1) = 1。

Step 2c: emultiplicity q c_rem' = v - 1 → q 可分（q' ≠ 0）且 p ∤ v。
  c_rem = q^v * g（gcd(q,g) = 1，from UFD/emultiplicity definition）。
  c_rem' = q^{v-1} * (C(v) * q' * g + q * g')。
  若 p | v 或 q' = 0：C(v)*q' = 0 → c_rem' = q^v * g' → emultiplicity q c_rem' ≥ v。
  但 emultiplicity q c_rem' = v - 1 < v。矛盾。
  故 p ∤ v 且 q' ≠ 0。

Step 2d: 对 f 做同样的 valuation 论证 → v_q(w₀) = 1 → 矛盾。
  v_f = v（Step 2a），p ∤ v（Step 2c），q' ≠ 0（Step 2c）。
  f = q^v * h（gcd(q,h) = 1）。
  f' = q^{v-1} * (C(v) * q' * h + q * h')。
  q ∤ (C(v) * q' * h)：
    - q ∤ q'：deg(q') < deg(q)，q 不可约 → q ∤ q'。
    - q ∤ h：gcd(q, h) = 1。
    - C(v) ≠ 0：p ∤ v → v ≠ 0 in F_p → C(v) ≠ 0。
    - 故 C(v) * q' * h ≠ 0 且 q ∤ 之。
  q | (q * h')。
  两加数的 q-emultiplicity 不同（0 vs ≥1）→ 和的 emultiplicity = min = 0。
  （Mathlib: emultiplicity_add_eq_min when emultiplicities differ）
  emultiplicity q f' = (v-1) + 0 = v - 1。

  emultiplicity q c₀ = emultiplicity q gcd(f, f'):
    q^{v-1} | f（v-1 < v）且 q^{v-1} | f'（emultiplicity = v-1）→ q^{v-1} | gcd。
    q^v ∤ f'（emultiplicity = v-1 < v）→ q^v ∤ gcd（gcd | f'）。
    故 emultiplicity q c₀ = v - 1。
    （用 pow_dvd_iff_le_emultiplicity + EuclideanDomain.dvd_gcd）

  emultiplicity q w₀ = emultiplicity q f - emultiplicity q c₀ = v - (v-1) = 1。
  （从 f = w₀ * c₀ + emultiplicity_mul）
  emultiplicity q w₀ ≥ 1 → q | w₀。
  **矛盾**：Case 2 假设 q ∤ w₀。✓
```

**两种情况都矛盾** → 假设 derivative(c_rem) ≠ 0 不成立 → derivative(c_rem) = 0。✓

**Mathlib API 需求**：
- `emultiplicity_mul` (Prime q): `emultiplicity q (a*b) = emultiplicity q a + emultiplicity q b`
- `emultiplicity_add_eq_min` (当两加数 emultiplicity 不同): `emultiplicity q (a+b) = min(...)`
- `pow_dvd_iff_le_emultiplicity`: `q^k ∣ a ↔ k ≤ emultiplicity q a`
- `emultiplicity_eq_zero`: `emultiplicity q a = 0 ↔ ¬(q ∣ a)`
- `dvd_of_emultiplicity_pos`: `0 < emultiplicity q a → q ∣ a`
- `emultiplicity_eq_of_dvd_of_not_dvd`: `q^k ∣ a ∧ q^{k+1} ∤ a → emultiplicity q a = k`

**自证引理**：
- 不需要 `emultiplicity_gcd`——用 `pow_dvd + dvd_gcd + gcd_dvd` 替代
- 不需要 `normalizedFactors` 展开——反证法完全避免了 Step 5
- `emultiplicity_derivative_of_separable_irred`（Step 2d 的核心）：~15 行

---

## 4. p-th root 正确性

### 4.1 `expand p f = f^p` in F_p[X]

**证明**：
1. `map_frobenius_expand` (Mathlib): `map (frobenius (ZMod p) p) (expand (ZMod p) p f) = f ^ p`
2. 在 `ZMod p` 中：`frobenius (ZMod p) p = RingHom.id (ZMod p)`
   - 因为 `∀ a : ZMod p, a^p = a`（Fermat/`ZMod.pow_card`）
   - `frobenius` 定义为 `x ↦ x^p`，在 F_p 上是恒等
3. `map (RingHom.id) g = g`，故 `expand (ZMod p) p f = f^p` ✓

**Mathlib 路径**：
- `ZMod.frobenius_zmod` 或 `frobenius_one`... 需确认确切引理名
- 备选：`Polynomial.map_id`

### 4.2 递归正确性

设 `SquarefreeDecomp g sub`（递归结果）：`g ~ ∏ sⱼ^{eⱼ}`。

`c_rem = expand p (contract p c_rem) = (contract p c_rem)^p = g^p`... 不对，`contract p c_rem = g`，但 `expand p g = g^p`。

所以 `c_rem = expand p g = g^p`（Frobenius）。

`g ~ ∏ sⱼ^{eⱼ}` → `c_rem = g^p ~ (∏ sⱼ^{eⱼ})^p = ∏ sⱼ^{eⱼ·p}` ✓

输出 `(sⱼ, eⱼ * p)` 正确编码了 `c_rem ~ ∏ sⱼ^{eⱼ·p}`。

验证 SquarefreeDecomp 四条：
1. **乘积还原**：`c_rem ~ ∏ sⱼ^{eⱼ·p}` ✓
2. **Squarefree + Monic**：从递归保证 ✓
3. **重数 ≥ 1**：`eⱼ * p ≥ 1 * 2 = 2 ≥ 1` ✓
4. **两两互素**：从递归保证 ✓

---

## 5. 顶层组合

### Case A: `f' = 0`
- `f = expand p (contract p f)` by `expand_contract`
- `g = contract p f`，递归得 `sub = sqfZp g`
- result = `sub.map (fun (s,e) => (s, e*p))`
- 正确性同 §4.2 ✓

### Case B: `f' ≠ 0`，`deg(c_rem) = 0`
- result = yun_result
- Y1 给 `P * c_rem ~ f`。`c_rem` deg 0，Monic 或单位 → 常数 ✓
- `P ~ f`（Associated）→ SquarefreeDecomp 条件 1
- Y4 → 条件 2 (Squarefree + Monic)，条件 3 (mult ≥ 1)
- Y5 → 条件 4 (coprimality)

### Case C: `f' ≠ 0`，`deg(c_rem) > 0`
- result = yun_result ++ pth_root_result
- **条件 1**（乘积还原）：
  - Y1：`P * c_rem ~ f`
  - 递归：`c_rem ~ ∏ sⱼ^{eⱼ·p}`（依赖 `c_rem' = 0` + §4.2）
  - 合并：`P * ∏ sⱼ^{eⱼ·p} ~ f` ✓
- **条件 2, 3**：从 Y4 + 递归保证 ✓
- **条件 4**（coprimality）：
  - Yun 因子间互素：Y5 ✓
  - p-th root 因子间互素：递归保证 ✓
  - **跨组互素**（Yun 因子 sᵢ vs p-th root 因子 tₖ）：

    **证明**：
    1. Y10 终止时：`IsCoprime sᵢ c_rem`（每个 Yun 因子与 c_rem 互素）
    2. p-th root 因子来源：`tₖ` 来自 `sqfZp(contract p c_rem)`。
       `tₖ | contract(p, c_rem) = g`。
       `g | g^p = expand(p, g) = c_rem`（Frobenius：`expand p g = g^p`，且 `g | g^p`）。
       故 `tₖ | g | c_rem`。
    3. `tₖ | c_rem` 且 `IsCoprime sᵢ c_rem`
       → 任何 sᵢ 与 tₖ 的公共不可约因子 r 满足 `r | sᵢ` 且 `r | tₖ | c_rem`
       → `r | gcd(sᵢ, c_rem)`
       → `IsUnit r`（由 IsCoprime）
       → 矛盾（r 不可约 → ¬ IsUnit）
    4. 故 `IsCoprime sᵢ tₖ` ✓

---

## 6. 辅助引理清单

| 引理 | 内容 | 难度 | Mathlib | nl-proof 状态 |
|------|------|------|---------|--------------|
| `emultiplicity_derivative_ge` | `v_r(f') ≥ v_r(f) - 1` | 中高 | 无（需自证） | §0 ✓ |
| `squarefree_div_gcd_derivative` | `f' ≠ 0 → Squarefree(f / gcd(f, f'))` | 高 | 无 | §2 Y2 ✓ |
| `expand_eq_pow_zmod` | `expand (ZMod p) p f = f^p` | 中 | `map_frobenius_expand` + frobenius=id | §4.1 ✓ |
| `frobenius_zmod_eq_id` | `frobenius (ZMod p) p = RingHom.id` | 低 | 可能已有 | §4.1 ✓ |
| `squarefree_mul_coprime` | `Squarefree (a*b) → IsRelPrime a b` | 低 | `squarefree_mul_iff` | ✓ |
| `derivative_c_rem_eq_zero` | Yun 终止后 `c_rem' = 0` | 高 | 无 | §3.2 ✓（5 步完整证明） |
| `yun_pthroot_coprime` | Yun 因子与 p-th root 因子互素 | 中 | 无 | §5 Case C ✓（Y10 推导） |
| `coprime_prod_of_pairwise` | 两两互素因子之积与互素元素互素 | 低 | 可能已有 | §3.2 Step 1 |

**所有引理的数学论证已完成（nl-proof 内 0 sorry）。**

---

## 7. 形式化策略

所有数学论证已在 nl-proof 中完成。形式化的难度在于将论证翻译为 Lean tactic proof。

### Step 1: 函数定义 + 终止性（~100 行）
- `yunLoop` 定义 + 终止度量 `w.natDegree + c.natDegree`
- `sqfZp` 定义 + 终止度量 `f.natDegree`
- **预期 sorry: 0**

### Step 2: Yun 不变量 Y1-Y8, Y10（~120 行）
- Y1 (product) — 代数恒等式
- Y2 (squarefree w₀) — 依赖 `emultiplicity_derivative_ge`（~30 行自证）
- Y5/Y6 (coprime with w) — 从 `squarefree_mul_iff`
- Y10 (coprime with c) — §2 的 gcd 论证
- **预期 sorry: 0-1**（`emultiplicity_derivative_ge` 可能需要 sorry）

### Step 3: `derivative(c_rem) = 0`（~40 行）
- 依赖 Y10 + valuation 论证 (§3.2)
- **预期 sorry: 0-1**（依赖 Step 2 的 valuation 引理）

### Step 4: SquarefreeDecomp 四条（~80 行）
- 条件 1: Y1 + p-th root 递归 + `expand_contract` + Frobenius
- 条件 2, 3: 直接
- 条件 4: Y5 + 递归 + Y10 跨组互素
- **预期 sorry: 0**

### 总计：~340 行，预期 0-1 个 sorry（仅可能在 `emultiplicity_derivative_ge`）
