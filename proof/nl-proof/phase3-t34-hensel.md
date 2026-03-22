# T3.4: Hensel 提升算法建模与正确性

> 状态：nl-proof v4（divByMonic 度数控制 + ker(π) 保持证明）
> 对应 C++：`polynomial_factorize_univar.hh:307-597` `__hensel_step` + `__hensel_lift`

---

## 0. 数学背景

设 `f ∈ Z[x]`，`p` 素数，`f ≡ g₁ · g₂ · ... · gᵣ (mod p)`（gᵢ monic 两两互素 in F_p[x]）。
Hensel 提升：∃ h₁,...,hᵣ with `f ≡ h₁·...·hᵣ (mod p^k)`，hᵢ ≡ gᵢ (mod p)，deg(hᵢ) = deg(gᵢ)。

---

## 1. 证明策略

**方案**：存在性证明。HenselCorrect spec 只要求 ∃ lifted factors。

**记号**：
- `φ_n := Int.castRingHom (ZMod n)`，`map_n := Polynomial.map φ_n`
- `π := ZMod.castHom (m ∣ m²) (ZMod m)` — 自然投射 ZMod(m²) → ZMod m
- `map_π := Polynomial.map π` — 多项式环上的投射

**关键性质**：`map_m = map_π ∘ map_{m²}`（`Polynomial.map_map`）。

**三步结构**：
1. `hensel_step`：2-factor 单步 mod m → mod m²（含度数保持）
2. `hensel_two_factor`：2-factor 迭代 mod p → mod p^{2^n}，投射到 p^k
3. `hensel_multifactor`：r-factor 归纳 → HenselCorrect

---

## 2. hensel_step（核心）

### 2.1 定理陈述

```lean
theorem hensel_step
    (m : ℕ) (hm : 1 < m)
    (f g h : Polynomial ℤ)
    (hprod : map_m f = map_m g * map_m h)
    (hcop : IsCoprime (map_m g) (map_m h))
    (hh_monic : Monic (map_m h))        -- h 在 ZMod m 下 monic
    : ∃ g' h' : Polynomial ℤ,
        -- (H1) f ≡ g'·h' (mod m²)
        map_{m²} f = map_{m²} g' * map_{m²} h'
        -- (H2) g' ≡ g (mod m)
        ∧ map_m g' = map_m g
        -- (H3) h' ≡ h (mod m)
        ∧ map_m h' = map_m h
        -- (H4) 度数保持
        ∧ natDegree (map_{m²} h') = natDegree (map_m h)
```

**注 1**：只要求 h 一侧的度数保持（`natDegree(map_{m²} h')` = `natDegree(map_m h)`）。
g 侧的度数保持通过乘积度数约束在多因子组合时推导。

**注 2**：需要 `hh_monic` 来做 divByMonic 并控制度数。

### 2.2 证明

在 ZMod(m²)[x] 中工作。设 ḡ = map_{m²} g, h̄ = map_{m²} h, f̄ = map_{m²} f。

**Step 1：ker(π) 的结构。**

π : ZMod(m²) →+* ZMod m（`ZMod.castHom`）。
ker(π) = 理想 (m) in ZMod(m²)。

**关键性质**：ker(π)² = 0。
证明：∀ a, b ∈ ker(π)，a = m·a'，b = m·b'（for some a', b' ∈ ZMod(m²)）。
a·b = m²·a'·b' = 0（因 m² = 0 in ZMod(m²)，by `ZMod.natCast_pow_eq_zero_of_le`）。✓

**推论**：若 p, q ∈ Polynomial(ZMod(m²)) 且 p, q 的每个系数 ∈ ker(π)，则 p·q = 0。
（p·q 的系数 = ∑ pᵢ·qⱼ，每项 pᵢ·qⱼ ∈ ker(π)² = 0。）✓

**Step 2：误差项 e。**

设 e := f̄ - ḡ·h̄ ∈ Polynomial(ZMod(m²))。

map_π(e) = map_π(f̄) - map_π(ḡ)·map_π(h̄) = map_m f - map_m g · map_m h = 0（by hprod）。

故 e 的每个系数 ∈ ker(π)。特别地 **e² = 0**（by 推论）。✓

**Step 3：Bézout 提升。**

`hcop` 给出 ∃ s̄, t̄ ∈ ZMod m[x]，`s̄ · map_m g + t̄ · map_m h = 1`。
（Lean: `hcop.exists` 或解构。）

用 `Polynomial.map_surjective`（`ZMod.intCast_surjective` → `map φ_m` 满射）：
∃ s, t ∈ Z[x]，`map_m s = s̄`，`map_m t = t̄`。

设 s̃ = map_{m²} s, t̃ = map_{m²} t。在 ZMod(m²)[x] 中：
`map_π(s̃·ḡ + t̃·h̄) = s̄·map_m g + t̄·map_m h = 1`。
故 `s̃·ḡ + t̃·h̄ = 1 + δ`，其中 δ ∈ ker(π)[x]。

**关键**：`e·δ = 0`（e 系数 ∈ ker(π)，δ 系数 ∈ ker(π)，ker(π)² = 0）。
故 `e·(1 + δ) = e + e·δ = e`。✓

即 `s̃·e·ḡ + t̃·e·h̄ = e·(s̃·ḡ + t̃·h̄) = e·(1+δ) = e`。✓

**Step 4：divByMonic 控制度数。**

h̄ = map_{m²} h。由 hh_monic 和 `map_m h = map_π h̄`：map_π h̄ monic in ZMod m[x]。

h̄ monic in ZMod(m²)[x]?
`leadingCoeff(h̄) = φ_{m²}(leadingCoeff(h))`。`π(leadingCoeff(h̄)) = φ_m(leadingCoeff(h)) = leadingCoeff(map_m h) = 1`（hh_monic）。
故 `leadingCoeff(h̄) = 1 + δ₀` with δ₀ ∈ ker(π)。δ₀ nilpotent → 1 + δ₀ unit（`IsNilpotent.isUnit_one_add`）。

但 Monic 要求 leadingCoeff = 1（精确），不是 unit。所以 h̄ 不一定 monic in ZMod(m²)[x]。

**修正**：不需要 h̄ monic in ZMod(m²)。用 `divByMonic` **in ZMod m[x]**（π(h̄) = map_m h IS monic），
然后提升结果。

具体：设 `s̃·e = s̃·(f̄ - ḡ·h̄)`。对 map_π(s̃·e) 做 divByMonic (by map_m h)：
`map_π(s̃·e) = q̄·map_m h + σ̄`，deg(σ̄) < deg(map_m h)。

但 map_π(e) = 0（Step 2），所以 map_π(s̃·e) = map_π(s̃)·map_π(e) = s̄·0 = 0。
故 `q̄ = 0, σ̄ = 0`。这又回到了平凡情况。

**问题**：在 ZMod m[x] 层 e 映射为 0，无法在该层做有意义的 divmod。

**正确方法**：在 ZMod(m²)[x] 中做 divByMonic。需要 h̄ 在 ZMod(m²)[x] 中 "essentially monic"。

**关键观察**：h̄ 的 leading coefficient 是 `1 + δ₀`（unit in ZMod(m²)）。
`modByMonic` 需要 Monic，但我们可以用 unit 归一化：
设 `h̄_norm = (1+δ₀)⁻¹ · h̄`。则 h̄_norm monic in ZMod(m²)[x]。
`modByMonic h̄_norm` 给出 `s̃·e = q·h̄_norm + σ`，deg(σ) < deg(h̄_norm) = deg(h̄)。

但这引入了 `(1+δ₀)⁻¹`，增加了证明复杂度。

### 2.3 构造中的度数控制（divByMonic）

**问题**：g̃ = ḡ + t̃·e, h̃ = h̄ + s̃·e 不控制度数。`natDegree(s̃·e)` 可能 > `natDegree(h̄)`。

**修正**：用 divByMonic 将修正项的度数限制在 deg(h̄) 以下。

**Step D1**：h̄ 在 ZMod(m²)[x] 中 "essentially monic"。

`π(leadingCoeff(h̄)) = leadingCoeff(map_m h) = 1`（hh_monic）。
设 u := leadingCoeff(h̄)。π(u) = 1 ≠ 0 → u ≠ 0。
u - 1 ∈ ker(π) → (u-1)² = 0 → u-1 nilpotent → `IsNilpotent.isUnit_one_add` → u unit。✓

设 h̄' := C(u⁻¹) * h̄。h̄' monic in ZMod(m²)[x]
（leadingCoeff(h̄') = u⁻¹ · u = 1）。natDegree(h̄') = natDegree(h̄)。✓

**Step D2**：divByMonic 分解。

`modByMonic_add_div` 对任意交换环 + monic 除数成立（Mathlib）。
对 s̃·e 做 divByMonic by h̄'：
```
s̃ · e = q · h̄' + σ    (modByMonic_add_div)
deg(σ) < deg(h̄') = deg(h̄)
```

**关键引理**：σ 和 q 的系数仍在 ker(π) 中。

证明：divByMonic 的算法逐步从 s̃·e 中减去 `coeff / leadingCoeff(h̄')` 乘以 h̄' 的位移。
- s̃·e 的系数 ∈ ker(π)（因 e 系数 ∈ ker(π)，ker(π) 是理想）。
- 每步商系数 = (当前 remainder 的 leading coeff) / 1 = (ker(π) element) / 1 = ker(π) element。
- 减去的项 = ker(π) element · h̄'，其系数 = ker(π) · ZMod(m²) ⊆ ker(π)。
- 归纳：每步 remainder 的系数仍 ∈ ker(π)。
故 q 系数 ∈ ker(π)，σ 系数 ∈ ker(π)。✓

（Lean 路径：对 `Polynomial.divModByMonic` 做归纳，或用 `modByMonic_eq_sub_mul_div` + 理想性质。
 若 Lean 中此引理不易直接证，可用替代方案：σ = s̃·e - q·h̄'，s̃·e 系数 ∈ ker(π)，
 q·h̄' 系数 = ∑ qᵢ·h̄'ⱼ。若 q 系数 ∈ ker(π)... 这又需要先证 q 系数 ∈ ker(π)。
 **实际简化**：不证 q, σ 各自在 ker(π) 中。直接证 σ ∈ ker(π)[x] 如下：
 `π(σ) = π(s̃·e) - π(q)·π(h̄') = π(s̃)·π(e) - π(q)·π(h̄') = s̄·0 - π(q)·π(h̄') = -π(q)·π(h̄')`。
 但 `π(s̃·e) = 0`（因 π(e) = 0），且 `π(q·h̄') = π(q)·π(h̄')`。
 `π(σ) = π(s̃·e) - π(q·h̄') = 0 - π(q)·π(h̄')`... 这不一定是 0。

 **更直接**：`π(s̃·e) = 0`。`π(q·h̄' + σ) = π(q)·π(h̄') + π(σ) = 0`。
 所以 `π(σ) = -π(q)·π(h̄')`。但 π(h̄') = π(u⁻¹)·π(h̄) = u⁻¹_mod_m · map_m h。这不是 0。
 所以 π(σ) 不一定是 0！

 **问题**：divByMonic 不保证余式系数在 ker(π) 中。）

**重新审视**：上面的论证有误。让我用更基本的方法。

**替代方案（直接构造，避免 divByMonic 的 ker 保持问题）**：

我们不需要 σ, q 各自在 ker(π) 中。我们需要的是：
1. σ·ḡ + τ·h̄ = e（Bézout 分解）
2. τ·σ = 0（保证乘积正确）
3. deg(σ) < deg(h̄)（度数控制）

**直接利用 ker(π)² = 0 得到 τ·σ = 0**：
只需 τ 或 σ 的系数在 ker(π) 中即可（不需要两者都在）。
实际上需要 τ·σ 的每个系数 ∈ ker(π)² = 0，即 ∑ τᵢ·σⱼ = 0。
若 τ 系数 ∈ ker(π) 且 σ 系数任意：τᵢ·σⱼ 不一定在 ker(π)² 中。
若 τ 系数 ∈ ker(π) 且 σ 系数 ∈ ker(π)：τᵢ·σⱼ ∈ ker(π)² = 0 ✓。

所以确实**需要两者都在 ker(π)**。

**最终方案**：在原始构造 g̃ = ḡ + t̃·e, h̃ = h̄ + s̃·e 的基础上，做一次 "度数修正"。

设 h̃₀ = h̄ + s̃·e（原始构造，系数在 ker(π) 中，但度数可能过大）。
设 h̃₀ = q₀·h̄' + r₀（divByMonic by h̄'，deg(r₀) < deg(h̄)）。

定义 h̃ = r₀ + C(u)·X^{deg(h̄')} ·... 不行，这改变了 h̃ 的结构。

**根本问题**：在 ZMod(m²)[x] 中，当 ker(π)² = 0 时，h̄ 的"高次修正项"（来自 s̃·e）可能使度数增大，但这些高次项在 mod m 下为 0。

**正确的度数论证**：即使 natDegree(h̃) > natDegree(h̄)，只要 `π(leadingCoeff(h̃)) = 0`（高次项在 ker(π) 中），就有 `natDegree(map_π h̃) = natDegree(map_m h)` < natDegree(h̃)。

问题是 HenselCorrect 比较 `natDegree(map_p gᵢ)` 和 `natDegree(map_{p^k} h'ᵢ)`。后者 = `natDegree(h̃)`。
若 natDegree(h̃) > natDegree(map_m h)，则 natDegree(map_{p^k} h'ᵢ) ≠ natDegree(map_p gᵢ)。

**终极修正**：选择 Z[x] lift g', h' 使得 natDegree 正确。

`Polynomial.map_surjective` 给出**某个** lift g' with map_{m²} g' = g̃。
但 natDegree(g') 可能 > natDegree(g̃)（Z[x] lift 的高次系数可以是 m² 的倍数）。
反之，natDegree(g') 也可能 < natDegree(g̃)（不太可能，map 不增加度数）。
实际上 natDegree(map φ g') ≤ natDegree(g')，且 map φ g' = g̃，所以 natDegree(g̃) ≤ natDegree(g')。

**但我们可以选择 g' 使 natDegree(g') = natDegree(g̃)**：
g̃ 有度数 d = natDegree(g̃)。选择 g' 的系数为 g̃ 的系数的 canonical lift to ℤ
（即 0 ≤ g'.coeff(i) < m²），高于 d 次的系数设为 0。
则 natDegree(g') = d = natDegree(g̃)。且 map_{m²} g' = g̃。✓

同理选择 h' with natDegree(h') = natDegree(h̃)。

那么 `natDegree(map_{m²} h') = natDegree(h̃)`（因为 map_{m²} h' = h̃，natDegree 保持）。

但 HenselCorrect 要 `natDegree(map_{p^k} h')` = ... hmm, 我们需要的是 map_{p^k} 的度数，不是 map_{m²}。

当通过 castHom 投射 ZMod(m²) → ZMod(p^k)：
`natDegree(map_π' h̃) ≤ natDegree(h̃)` (natDegree_map_le)。
也 `natDegree(map_π' h̃) ≥ natDegree(map_p h̃) = natDegree(map_m h)`
（castHom from p^k to p 的 natDegree_map_le 给出 `natDegree(map_p(map_π' h̃)) ≤ natDegree(map_π' h̃)`，
 即 `natDegree(map_m h) ≤ natDegree(map_π' h̃)`）。

所以 `natDegree(map_m h) ≤ natDegree(map_{p^k} h') ≤ natDegree(h̃)`。

若 natDegree(h̃) = natDegree(map_m h)：则一切夹逼得等式。✓
若 natDegree(h̃) > natDegree(map_m h)：无法保证等式。✗

**结论**：度数保持**要求 natDegree(h̃) = natDegree(map_m h)**，这要求构造中控制度数。

---

### 2.3 最终方案：构造中控制度数

**将 h 侧的修正限制为 deg < deg(h̄)**。

重新定义构造。设 h̄' = C(u⁻¹) · h̄（monic in ZMod(m²)[x]，u = leadingCoeff(h̄)）。

对 Bézout 分解 `s̃·e·ḡ + t̃·e·h̄ = e` 中的 `s̃·e` 做 modByMonic by h̄'：
```
s̃ · e = q · h̄' + σ     (modByMonic_add_div，deg(σ) < deg(h̄'))
```

定义：
```
h̃ := h̄ + σ           (deg(σ) < deg(h̄) → natDegree(h̃) = natDegree(h̄))
g̃ := ḡ + t̃·e + q·C(u⁻¹)·ḡ
```

**验证 g̃·h̃ = f̄**：

先验证 Bézout 分解在新分割下成立：
```
σ·ḡ + (t̃·e + q·C(u⁻¹)·ḡ)·h̄
= σ·ḡ + t̃·e·h̄ + q·C(u⁻¹)·ḡ·h̄
= σ·ḡ + t̃·e·h̄ + q·h̄'·ḡ          [因 C(u⁻¹)·h̄ = h̄']
= (σ + q·h̄')·ḡ + t̃·e·h̄
= s̃·e·ḡ + t̃·e·h̄                   [因 s̃·e = q·h̄' + σ]
= e·(s̃·ḡ + t̃·h̄)
= e·(1 + δ)
= e                                    [因 e·δ = 0]
```
✓

即：设 τ := t̃·e + q·C(u⁻¹)·ḡ。则 `σ·ḡ + τ·h̄ = e`。

现在验证 `σ 系数 ∈ ker(π)` 和 `τ 系数 ∈ ker(π)`：

**σ 系数 ∈ ker(π)**：
`σ = s̃·e - q·h̄'`。`π(σ) = π(s̃·e) - π(q)·π(h̄') = s̄·0 - π(q)·π(h̄')`。
但 `π(s̃·e) = 0`（因 π(e) = 0），所以 `π(σ) = -π(q)·π(h̄')`。
又 `0 = π(s̃·e) = π(q·h̄' + σ) = π(q)·π(h̄') + π(σ)`。
故 `π(σ) = -π(q)·π(h̄')`，且 `π(q)·π(h̄') + π(σ) = 0`。
这不给出 π(σ) = 0。所以 **σ 系数不一定在 ker(π) 中**。

**问题仍然存在**。divByMonic 后的余式不保证在 ker(π) 中。

**根本问题重新审视**：

在原始构造（h̃ = h̄ + s̃·e）中：
- s̃·e 系数 ∈ ker(π) ✓（因 e 系数 ∈ ker(π)，ker(π) 是理想）
- 所以 h̃ - h̄ 系数 ∈ ker(π) ✓
- 乘积正确（ker² = 0）✓
- **但度数不受控** ✗

若我们做 modByMonic 分割 s̃·e = q·h̄' + σ，那么 σ 和 q 不一定在 ker(π) 中。

**关键观察**：设 s̃·e 的次数为 d + k（d = deg(h̄')，k ≥ 0）。divByMonic 的第一步：
商系数 = coeff_{d+k}(s̃·e) / 1 = coeff_{d+k}(s̃·e) ∈ ker(π)（s̃·e 系数 ∈ ker(π)）。
减去 coeff_{d+k}(s̃·e) · X^k · h̄'：
新 remainder = s̃·e - coeff_{d+k}(s̃·e) · X^k · h̄'。

`coeff_{d+k}(s̃·e) · X^k · h̄'` 的系数：coeff_{d+k}(s̃·e) ∈ ker(π)，乘以 h̄' 的系数。
ker(π) · ZMod(m²) = ker(π)（理想性质）。所以每个系数 ∈ ker(π)。✓

新 remainder 的系数 = (s̃·e 的系数) - (ker(π) element) = (ker(π) element) - (ker(π) element) = ker(π) element。✓

归纳：divByMonic 的每一步，remainder 的系数都 ∈ ker(π)。✓
因此 σ 系数 ∈ ker(π) ✓，q 系数 ∈ ker(π) ✓。

**之前的反例分析有误**：`π(σ) = -π(q)·π(h̄')` 中，若 σ, q 系数 ∈ ker(π)，则 π(σ) = 0, π(q) = 0，等式变为 0 = 0。✓

**divByMonic 保持 ker(π) 系数的完整证明**：

引理：设 a ∈ Polynomial(ZMod(m²))，a 的所有系数 ∈ ker(π)。
设 b monic in Polynomial(ZMod(m²))。
则 a modByMonic b 和 a divByMonic b 的系数都 ∈ ker(π)。

证明（对 natDegree(a) 归纳）：
- 若 natDegree(a) < natDegree(b)：q = 0（系数 = 0 ∈ ker(π)），r = a（系数 ∈ ker(π)）。✓
- 若 natDegree(a) ≥ natDegree(b)：
  一步 divByMonic：q₀ = leadingCoeff(a)（∈ ker(π)）。
  a' = a - C(q₀) · X^{deg(a)-deg(b)} · b。
  C(q₀) · X^k · b 的系数 = q₀ · bⱼ ∈ ker(π)（q₀ ∈ ker(π)，ker(π) 是理想）。
  a' 的系数 = a 的系数 - ker(π) element = ker(π) element。
  natDegree(a') < natDegree(a)（b monic → leading term 消去）。
  归纳假设适用于 a'。✓

（Lean 路径：可能需要对 `Polynomial.divModByMonicAux` 做归纳。
 或者直接证 `map π (a modByMonic b) = (map π a) modByMonic (map π b) = 0 modByMonic (map π b) = 0`，
 因此 `a modByMonic b` 系数 ∈ ker(π)。
 **这是最简洁的路径**：`Polynomial.map_modByMonic` 给出 `map π (a modByMonic b) = (map π a) modByMonic (map π b)`。
 `map π a = 0`（因 a = s̃·e，π(e) = 0）→ `map π (a modByMonic b) = 0 modByMonic (map π b) = 0`。
 故 σ = a modByMonic b 的系数 ∈ ker(π)。✓ 同理 q。）

### 2.4 hensel_step 完整证明总结

1. 设 e = f̄ - ḡ·h̄ in ZMod(m²)[x]。e 系数 ∈ ker(π)。e² = 0。
2. Bézout lift：s̃·ḡ + t̃·h̄ = 1 + δ，δ ∈ ker(π)[x]。e·δ = 0。
   故 `s̃·e·ḡ + t̃·e·h̄ = e`（精确）。
3. h̄' = C(u⁻¹)·h̄ monic（u = leadingCoeff(h̄), u unit in ZMod(m²)）。
4. divByMonic：s̃·e = q·h̄' + σ，deg(σ) < deg(h̄')。
   σ, q 系数 ∈ ker(π)（by `map_modByMonic` 或归纳论证）。
5. 定义 τ = t̃·e + q·C(u⁻¹)·ḡ。τ 系数 ∈ ker(π)（t̃·e ∈ ker(π)[x]，q·anything ∈ ker(π)[x]）。
6. 验证 `σ·ḡ + τ·h̄ = e`（§2.3 的推导）。
7. 定义 g̃ = ḡ + τ，h̃ = h̄ + σ。
8. **乘积**：g̃·h̃ = ḡ·h̄ + (σ·ḡ + τ·h̄) + τ·σ = ḡ·h̄ + e + 0 = f̄。
   （τ·σ = 0 因 τ, σ 系数 ∈ ker(π)，ker(π)² = 0。）✓
9. **mod m 保持**：π(g̃) = π(ḡ) + π(τ) = π(ḡ) + 0 = map_m g。
   π(h̃) = π(h̄) + π(σ) = π(h̄) + 0 = map_m h。✓
10. **度数**：deg(σ) < deg(h̄) → leadingCoeff(h̃) = leadingCoeff(h̄) = u。
    π(u) = 1 ≠ 0 → `natDegree_map_of_leadingCoeff_ne_zero` →
    natDegree(map_π h̃) = natDegree(h̃) → natDegree(map_m h) = natDegree(h̃)。✓
11. **lift**：∃ g', h' ∈ Z[x]，map_{m²} g' = g̃, map_{m²} h' = h̃。
    natDegree(map_{m²} h') = natDegree(h̃) = natDegree(map_m h)。✓

---

## 3. IsCoprime 传播

**引理**：`IsCoprime (map_m g) (map_m h)` → `IsCoprime (map_{m²} g') (map_{m²} h')`
（其中 g', h' 是 hensel_step 的输出）。

**证明**：
map_π(map_{m²} g') = map_m g'= map_m g（H2）。
map_π(map_{m²} h') = map_m h' = map_m h（H3）。

从 `IsCoprime (map_m g) (map_m h)` 取 Bézout：∃ ŝ, t̂，`ŝ·map_m g + t̂·map_m h = 1`。
Lift to Z[x]：ŝ₀, t̂₀。
`map_{m²} ŝ₀ · map_{m²} g' + map_{m²} t̂₀ · map_{m²} h'` maps under π to
`ŝ·map_m g + t̂·map_m h = 1`。
故 = `1 + δ'`，δ' ∈ ker(π)[x]。

`δ'² = 0`（ker(π)² = 0 in ZMod(m²)，推广到多项式环）。
`1 + δ'` 的 nilpotence：δ' nilpotent of order 2 in Polynomial(ZMod(m²))。
（`IsNilpotent δ'` 因 δ'² = 0。）
`IsNilpotent.isUnit_one_add`：`1 + δ'` unit in Polynomial(ZMod(m²))。✓

设 v = (1+δ')⁻¹。`(v·map_{m²} ŝ₀)·map_{m²} g' + (v·map_{m²} t̂₀)·map_{m²} h' = 1`。
故 `IsCoprime (map_{m²} g') (map_{m²} h')` ✓。

---

## 4. hensel_two_factor（迭代）

### 4.1 定理陈述

```lean
theorem hensel_two_factor
    (p : ℕ) (hp : Nat.Prime p) (k : ℕ) (hk : 0 < k)
    (f g h : Polynomial ℤ)
    (hprod : map_p f = map_p g * map_p h)
    (hcop : IsCoprime (map_p g) (map_p h))
    (hh_monic : Monic (map_p h))
    : ∃ g' h' : Polynomial ℤ,
        map_{p^k} f = map_{p^k} g' * map_{p^k} h'
        ∧ map_p g' = map_p g
        ∧ map_p h' = map_p h
        ∧ natDegree (map_{p^k} h') = natDegree (map_p h)
```

### 4.2 证明

**对 n 归纳**（二次 doubling），目标 p^{2^n}：

**n = 0**：p^1 = p。取 g' = g, h' = h。✓

**n → n+1**：
IH：∃ g_n, h_n with
- `map_{p^{2^n}} f = map_{p^{2^n}} g_n * map_{p^{2^n}} h_n`
- `map_p g_n = map_p g`，`map_p h_n = map_p h`
- `natDegree(map_{p^{2^n}} h_n) = natDegree(map_p h)`
- `IsCoprime (map_{p^{2^n}} g_n) (map_{p^{2^n}} h_n)`（从 §3 传播）

应用 hensel_step with m = p^{2^n}：
- hprod ✓
- hcop ✓
- hh_monic：`map_{p^{2^n}} h_n` monic?
  `natDegree(map_{p^{2^n}} h_n) = natDegree(map_p h)` and `map_p h` monic。
  Need: `map_{p^{2^n}} h_n` monic in ZMod(p^{2^n})[x]。
  `leadingCoeff(map_{p^{2^n}} h_n)` at position `natDegree(map_p h)`。
  Hmm, actually `Monic(map_{p^{2^n}} h_n)` would require leadingCoeff = 1 in ZMod(p^{2^n})。
  由 `map_p h_n = map_p h` monic，leadingCoeff(map_p h_n) = 1。
  `castHom(leadingCoeff(map_{p^{2^n}} h_n)) = leadingCoeff(map_p h_n) = 1`
  （当 natDegree 在 castHom 下保持——由 IH 的度数条件给出）。
  故 `leadingCoeff(map_{p^{2^n}} h_n) = 1 + δ`（mod p^{2^n}，not necessarily monic）。

**问题**：hensel_step 需要 `Monic (map_m h)` 但我们只有 "leadingCoeff maps to 1 under castHom"。

**修正**：放宽 hensel_step 的条件：不需要 `Monic (map_m h)`，只需要
`π(leadingCoeff(map_{m²} h)) ≠ 0`（或等价地，`Monic (map_m h)` 的弱化版本）。

实际上度数保持论证（§2.4）只需 `π(leadingCoeff(h̃)) ≠ 0`，不需要精确 = 1。
`leadingCoeff(h̃) = leadingCoeff(h̄ + s̃·e)`。
当 `natDegree(s̃·e) < natDegree(h̄)` 时，`leadingCoeff(h̃) = leadingCoeff(h̄)`。
但我们不知道 `natDegree(s̃·e) < natDegree(h̄)`。

**更好的方法**：直接在迭代中跟踪 `map_p h' monic`（始终 mod p monic），
而不是跟踪 `map_{p^{2^n}} h'` monic。

因为 `map_p h' = map_p h`（mod p 保持），且 `map_p h` monic，
所以 h' 始终满足 "mod p monic"。

对于度数保持：我们需要在 map_{p^{2^n}} 层的度数。

**关键简化**：hensel_step 只输出 (H4) `natDegree(map_{m²} h') = natDegree(map_m h)`。
迭代时链式得到：
- `natDegree(map_{p^2} h₁) = natDegree(map_p h)` ← hensel_step(m=p)
- `natDegree(map_{p^4} h₂) = natDegree(map_{p^2} h₁) = natDegree(map_p h)` ← hensel_step(m=p²)
- ...
- `natDegree(map_{p^{2^n}} h_n) = natDegree(map_p h)`

每步的 hensel_step 需要 `Monic (map_m h_{prev})`。但 `map_m h_{prev}` 可能不 monic。

**真正需要的**：hensel_step 的 hh_monic 条件改为：
`π(leadingCoeff(map_{m²} h)) ≠ 0`
即 `leadingCoeff(map_m h) ≠ 0`（不需要 = 1）。

在 `natDegree_map_of_leadingCoeff_ne_zero` 中只需 `φ(leadingCoeff p) ≠ 0`，不需要 = 1。
故度数论证不需要 Monic，只需 leading coeff 非零 under π。

**修正 hensel_step**：将 `hh_monic` 替换为 `leadingCoeff (map_m h) ≠ 0`。

迭代时：
- `leadingCoeff(map_p h) ≠ 0`（因 map_p h monic → leadingCoeff = 1 ≠ 0）✓
- `leadingCoeff(map_{p²} h₁) ≠ 0`？
  由 `natDegree(map_{p²} h₁) = natDegree(map_p h)`，
  且 `castHom(leadingCoeff(map_{p²} h₁)) = leadingCoeff(map_p h₁) = leadingCoeff(map_p h) = 1 ≠ 0`。
  若 `leadingCoeff(map_{p²} h₁) = 0`：castHom(0) = 0 ≠ 1。矛盾。✓

故 `leadingCoeff(map_{p^{2^n}} h_n) ≠ 0` 在每步归纳中保持。✓

### 4.3 投射到 p^k

doubling 给出 mod p^{2^n} 的结果。对于目标 p^k（2^n ≥ k）：

`ZMod.castHom (p^k ∣ p^{2^n})` 给出投射 π_k : ZMod(p^{2^n}) → ZMod(p^k)。
`map π_k (map_{p^{2^n}} f) = map_{p^k} f`（by `Polynomial.map_map`）。
`map π_k (map_{p^{2^n}} g' * map_{p^{2^n}} h') = map_{p^k} g' * map_{p^k} h'`（by `map_mul`）。
结合 product identity → `map_{p^k} f = map_{p^k} g' * map_{p^k} h'` ✓。

mod p 保持不变（因 p | p^{2^n} 和 p | p^k，castHom 组合一致）。✓

度数：`natDegree(map_{p^k} h') ≤ natDegree(map_{p^{2^n}} h') = natDegree(map_p h)`。
但也 `natDegree(map_{p^k} h') ≥ natDegree(map_p h')`（castHom from p^k to p 的 natDegree_map_le）。
`natDegree(map_p h') = natDegree(map_p h)`。
故 `natDegree(map_p h) ≤ natDegree(map_{p^k} h') ≤ natDegree(map_p h)`。等式 ✓。

---

## 5. hensel_multifactor

### 5.1 定理陈述

```lean
theorem hensel_multifactor
    (p : ℕ) (hp : Nat.Prime p) (k : ℕ) (hk : 0 < k)
    (f : Polynomial ℤ)
    (factors : List (Polynomial (ZMod p)))
    (hprod : map_p f = factors.prod)
    (hcop : factors.Pairwise IsCoprime)
    (hmonic : ∀ g ∈ factors, Monic g)
    : ∃ lifted : List (Polynomial (ZMod (p^k))),
        map_{p^k} f = lifted.prod
        ∧ factors.length = lifted.length
        ∧ List.Forall₂ (fun g h => natDegree g = natDegree h) factors lifted
```

### 5.2 证明（对 factors 长度归纳）

**Base [] (length = 0)**：
`map_p f = [].prod = 1`。
lifted = []。`map_{p^k} f = [].prod = 1`？只有 f ≡ 1 (mod p) 时。
实际上 `map_p f = 1` 意味着 f 的所有系数 ≡ 0 (mod p) 除了常数项 ≡ 1。
`map_{p^k} f` 的常数项 ≡ f 的常数项 (mod p^k)，不一定 = 1。

**但 HenselCorrect spec 要求精确等式 `map_{p^k} f = lifted.prod`**。
当 factors = []：prod = 1。需 map_{p^k} f = 1。不一定成立。

**修正**：Base case 不需要处理（Hensel 提升的 input 来自 Zp 因式分解，factors 非空）。
或：添加前置条件 `factors.length ≥ 1`。

**实际上 factors = [] 意味着 f 是 unit mod p**。这种情况下 Hensel 提升无意义。可以安全跳过。

**Base [g] (length = 1)**：
factors = [g]。`map_p f = g`。f 本身就是唯一因子。
lifted = [map_{p^k} f]。
- `map_{p^k} f = [map_{p^k} f].prod` ✓
- length 1 = 1 ✓
- `natDegree g = natDegree (map_{p^k} f)`？
  g = map_p f。natDegree g = natDegree(map_p f)。
  natDegree(map_{p^k} f)：由 g monic → map_p f monic → p ∤ leadingCoeff f。
  `natDegree_map_of_leadingCoeff_ne_zero`：p ∤ lc(f) → natDegree(map_p f) = natDegree(f)。
  p^k 也 ∤ lc(f)（因 p | p^k，p ∤ lc(f) → p^k ∤ lc(f)）。
  故 natDegree(map_{p^k} f) = natDegree(f) = natDegree(map_p f) = natDegree(g) ✓。

**Step (g :: rest)**：
factors = g :: rest（length ≥ 2）。
`map_p f = g · rest.prod`。

`IsCoprime g rest.prod`（from `factors.Pairwise IsCoprime` + `IsCoprime.prod_right`/`prod_left_iff`）。
`rest.prod` monic（from `hmonic` + `Monic.prod`... `monic_list_prod`）。

应用 hensel_two_factor with g 和 h (where map_p h = rest.prod)：
需 Z[x] lift h with `map_p h = rest.prod`（`Polynomial.map_surjective`）。

hensel_two_factor 给出 g', h' with：
- `map_{p^k} f = map_{p^k} g' · map_{p^k} h'`
- `map_p g' = map_p g = g`... wait, g ∈ ZMod p[x] 而 g' ∈ Z[x]。

**修正类型**：hensel_two_factor 的输入 g, h ∈ Z[x]，输出也是 Z[x]。
factors 是 ZMod p[x]。需要 Z[x] lifts。

设 g̃ ∈ Z[x] with map_p g̃ = g（from map_surjective）。
设 h̃ ∈ Z[x] with map_p h̃ = rest.prod。
`map_p f = map_p g̃ · map_p h̃ = g · rest.prod = factors.prod`。
`IsCoprime (map_p g̃) (map_p h̃) = IsCoprime g rest.prod` ✓。
`Monic (map_p h̃) = Monic rest.prod` ✓。

hensel_two_factor → ∃ g', h'：
- `map_{p^k} f = map_{p^k} g' · map_{p^k} h'`
- `map_p g' = g`，`map_p h' = rest.prod`
- `natDegree(map_{p^k} h') = natDegree(rest.prod)`

归纳假设对 h' 和 rest 适用：
- `map_p h' = rest.prod` ✓
- `rest.Pairwise IsCoprime`（from factors.Pairwise 的子列表）✓
- `∀ r ∈ rest, Monic r`（from hmonic）✓

IH → ∃ lifted_rest：
- `map_{p^k} h' = lifted_rest.prod`
- `rest.length = lifted_rest.length`
- `Forall₂ natDeg_eq rest lifted_rest`

组合：lifted = [map_{p^k} g'] ++ lifted_rest。
- `map_{p^k} f = map_{p^k} g' · map_{p^k} h' = map_{p^k} g' · lifted_rest.prod = lifted.prod` ✓
- length = 1 + rest.length = factors.length ✓
- `natDegree g = natDegree (map_{p^k} g')`？
  From product：`natDegree(map_{p^k} f) = natDegree(map_{p^k} g') + natDegree(map_{p^k} h')`
  （需要乘积度数加法——要求两者 leading coeff 乘积非零，由 monic/non-vanishing 条件保证）。
  `natDegree(map_p f) = natDegree(g) + natDegree(rest.prod)`（ZMod p 是域，度数加法成立）。
  `natDegree(map_{p^k} f) = natDegree(map_p f)`（lc(f) 非零 mod p → 非零 mod p^k）。
  `natDegree(map_{p^k} h') = natDegree(rest.prod)`（hensel_two_factor 给出）。
  故 `natDegree(map_{p^k} g') = natDegree(map_{p^k} f) - natDegree(map_{p^k} h') = natDegree(g)` ✓。

  **但**：`natDegree(map_{p^k} f) = natDegree(map_{p^k} g') + natDegree(map_{p^k} h')` 需要两者乘积的 leading coeff 非零。如上 §2.3 论证：h' 的 leading coeff 在某种意义上非零。

  具体：从 `natDegree(map_{p^k} h') = natDegree(rest.prod)` 和 `map_p h' = rest.prod` monic：
  castHom(leadingCoeff(map_{p^k} h')) = leadingCoeff(map_p h') = 1 ≠ 0。
  故 leadingCoeff(map_{p^k} h') ≠ 0（若 = 0 则 castHom = 0 ≠ 1）。
  同理从 natDegree 等式推出 leadingCoeff(map_{p^k} g') ≠ 0。

  积的度数加法：`natDegree(a*b) = natDegree(a) + natDegree(b)` 当
  `leadingCoeff(a) * leadingCoeff(b) ≠ 0`。
  由上：两者非零，且其中 leadingCoeff(h') 的 castHom = 1（unit under castHom），
  但在 ZMod(p^k) 中是否 non-zero-divisor？

  leadingCoeff(map_{p^k} h') = c，castHom(c) = 1。在 ZMod(p^k) 中 c = 1 + p·d。
  c 是 unit（1 + nilpotent）。故 c · anything = 0 → anything = 0。
  故 leadingCoeff(g') · c ≠ 0（因 leadingCoeff(g') ≠ 0 且 c unit）。
  `natDegree(map_{p^k} g' · map_{p^k} h') = natDegree(map_{p^k} g') + natDegree(map_{p^k} h')` ✓。

---

## 6. 辅助引理清单

| 引理 | Mathlib | 需自证 |
|------|---------|--------|
| `Polynomial.map_surjective` | ✅ | — |
| `IsNilpotent.isUnit_one_add` | ✅ | — |
| `ZMod.castHom` + `castHom_comp` | ✅ | — |
| `IsCoprime.prod_left_iff` | ✅ | — |
| `ZMod.natCast_pow_eq_zero_of_le` | ✅ | — |
| `natDegree_map_of_leadingCoeff_ne_zero` | ✅ | — |
| `Polynomial.map_map` | ✅ | — |
| ker(π)² = 0 in ZMod(m²) | 组合已有引理 | ~3 行 |
| poly_coeff_in_ker → poly² = 0 | 从 ker² = 0 推 | ~5 行 |
| `Polynomial.map_modByMonic` | modByMonic 与 map 交换 | 确认 Mathlib |
| divByMonic 保持 ker(π) 系数 | via `map_modByMonic` + π(a)=0 | ~5 行 |
| u = leadingCoeff(h̄) is unit | 1 + nilpotent is unit | ~3 行 |

---

## 7. Lean 翻译方案（ℤ[x] 方法）

### 7.1 总体策略

**放弃在 ZMod(m²)[x] 内部工作**。改为在 ℤ[x] 中构造 g', h'，只在最后用 `Polynomial.map` 投射。

原因：ZMod(m²) 中 "π(x) = 0 → ∃ a, x = m*a" 的提取在 Lean 中非常困难（需要 ZMod 结构论）。
在 ℤ[x] 中直接乘以 C(m) 避免了这个问题。

### 7.2 Mathlib API 速查（已确认存在）

| 需要 | 正确 Lean 名 | 文件 |
|------|-------------|------|
| `(a : ZMod b) = 0 ↔ b ∣ a`（ℤ版） | `ZMod.intCast_zmod_eq_zero_iff_dvd` | ZMod/Basic.lean:505 |
| `(a : ZMod b) = 0 ↔ b ∣ a`（ℕ版） | `ZMod.natCast_eq_zero_iff` | ZMod/Basic.lean:511 |
| `castHom h R i = cast i` | `ZMod.castHom_apply`（= rfl） | ZMod/Basic.lean:334 |
| `cast((k:ℤ) : ZMod n) = k` in target | `ZMod.cast_intCast` | ZMod/Basic.lean:354 |
| `coeff (map f p) n = f (coeff p n)` | `Polynomial.coeff_map`（@[simp]） | Coeff.lean:79 |
| `(p - q).map f = p.map f - q.map f` | `Polynomial.map_sub`（protected） | Eval/Defs.lean:724 |
| `map f (map g p) = map (f.comp g) p` | `Polynomial.map_map` | Eval/Defs.lean |
| `φ surj → map φ surj` | `Polynomial.map_surjective` | Coeff.lean:109 |
| `(p %ₘ q).map f = (p.map f) %ₘ (q.map f)` | `Polynomial.map_modByMonic` | Div.lean:382 |
| `ZMod.val a` cast 回 ZMod = a | `ZMod.natCast_zmod_val` | ZMod/Basic.lean:204 |
| `ℤ → ZMod n` 满射 | `ZMod.intCast_surjective` | ZMod/Basic.lean |
| `IsNilpotent a → IsUnit (1+a)` | `IsNilpotent.isUnit_one_add` | Nilpotent/Basic.lean:71 |

### 7.3 hensel_step 的 ℤ[x] 构造

**输入**：f, g, h : ℤ[x], map_m f = map_m g * map_m h, IsCoprime, Monic(map_m h)。

**Step A**：提取 ℤ[x] 中的误差多项式。
```
map_m (f - g * h) = 0
→ ∀ i, (m : ℤ) ∣ (f - g * h).coeff i        [by intCast_zmod_eq_zero_iff_dvd + coeff_map]
→ ∃ e_int : ℤ[x], f - g * h = C (↑m) * e_int  [逐系数构造]
```

**Lean 构造 e_int**：
```lean
-- 方法 1：Finsupp.mapRange
let e_int := Polynomial.ofFinsupp
  ((f - g * h).toFinsupp.mapRange (· / (m : ℤ)) (by simp))
-- 验证：C(m) * e_int = f - g * h
-- 每个系数：m * ((f-g*h).coeff i / m) = (f-g*h).coeff i  [Int.ediv_mul_cancel' (hdvd_i)]
```

```lean
-- 方法 2（更简洁）：直接构造，避免 Finsupp 操作
-- 由 map_m (f - g*h) = 0，存在性由 Polynomial 的 coeff 整除性保证
-- 可用 Polynomial.ext + Int.ediv_mul_cancel'
```

**Step B**：Bézout 提升到 ℤ[x]。
```
IsCoprime (map_m g) (map_m h)
→ ∃ s_bar t_bar, s_bar * map_m g + t_bar * map_m h = 1   [IsCoprime.exists 或解构]
→ ∃ s t : ℤ[x], map_m s = s_bar, map_m t = t_bar          [Polynomial.map_surjective]
→ map_m (s * g + t * h - 1) = 0
→ ∃ w : ℤ[x], s * g + t * h - 1 = C(↑m) * w               [同 Step A]
即 s * g + t * h = 1 + C(↑m) * w
```

**Step C**：在 ZMod m[x] 中做 divByMonic + Bézout 分解。
```
map_m h monic（hh_monic）
(map_m s * map_m e_int) %ₘ (map_m h) = σ_bar   [modByMonic in ZMod m[x]]
deg(σ_bar) < deg(map_m h)

τ_bar = map_m t * map_m e_int + ((map_m s * map_m e_int) /ₘ (map_m h)) * map_m g

验证：σ_bar * map_m g + τ_bar * map_m h = map_m e_int
  [展开 + Bézout + modByMonic_add_div]
```

**Step D**：提升 σ_bar, τ_bar 到 ℤ[x]（有界度数）。
```
∃ σ_int : ℤ[x], map_m σ_int = σ_bar ∧ natDegree σ_int ≤ natDegree σ_bar
∃ τ_int : ℤ[x], map_m τ_int = τ_bar
```

**σ_int 的有界度数 lift**：
```lean
-- 用 ZMod.val 做 canonical lift：σ_int.coeff(i) = ZMod.val(σ_bar.coeff(i))
-- natDegree(σ_int) ≤ natDegree(σ_bar)  [所有高次系数 = 0]
-- map_m σ_int = σ_bar  [ZMod.natCast_zmod_val]
```
注：τ_int 不需要有界度数（只有 σ_int 需要，用于 h' 的度数控制）。

**Step E**：定义 g', h'。
```
g' := g + C(↑m) * τ_int
h' := h + C(↑m) * σ_int
```

**Step F**：验证 (H1) 乘积 mod m²。
```
g' * h' = (g + C(m)*τ)(h + C(m)*σ)
        = g*h + C(m)*(τ*h + g*σ) + C(m)²*τ*σ
        = g*h + C(m)*(τ*h + g*σ) + C(m²)*τ*σ    [C(m)² = C(m²)]
```
map_{m²}(g'*h') = map_{m²}(g*h) + map_{m²}(C(m)*(τ*h + g*σ)) + map_{m²}(C(m²)*τ*σ)

- `map_{m²}(C(m²)*τ*σ) = 0`：因 `(m² : ZMod(m²)) = 0` → `map_{m²}(C(m²)) = C(0) = 0`。
  Lean: `simp [Polynomial.map_C, ZMod.natCast_self]`

- `map_{m²}(C(m)*(τ*h + g*σ)) = map_{m²}(C(m)*e_int)`：
  因 `τ*h + g*σ ≡ e_int (mod m)`（from Step C: map_m(σ*g + τ*h) = map_m(e_int)）。
  即 `map_m(τ*h + g*σ - e_int) = 0` → `∃ d, τ*h + g*σ - e_int = C(m)*d`。
  → `C(m)*(τ*h + g*σ) = C(m)*e_int + C(m²)*d`。
  → `map_{m²}(C(m)*(τ*h + g*σ)) = map_{m²}(C(m)*e_int) + 0 = map_{m²}(C(m)*e_int)`。

- `map_{m²}(g*h + C(m)*e_int) = map_{m²}(g*h + (f - g*h)) = map_{m²}(f)`。

故 `map_{m²}(g'*h') = map_{m²}(f)` ✓。

Lean：全部是 `ring` + `simp [map_C, map_mul, map_add, ZMod.natCast_self]`。

**Step G**：验证 (H2, H3) mod m 保持。
```
map_m(g') = map_m(g) + map_m(C(m)*τ)
map_m(C(m)) = C((m : ZMod m)) = C(0) = 0
→ map_m(C(m)*τ) = 0
→ map_m(g') = map_m(g)  ✓
```
Lean: `simp [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_C, ZMod.natCast_self]`

**Step H**：验证 (H4) 度数保持。
```
h' = h + C(m) * σ_int
leadingCoeff(h) mod p 非零（from monic chain）
natDegree(σ_int) ≤ natDegree(σ_bar) < natDegree(map_m h) = natDegree(h)
  [后者因 p ∤ leadingCoeff(h) → natDegree_map_of_leadingCoeff_ne_zero]
C(m)*σ_int 的 natDegree ≤ natDegree(σ_int) < natDegree(h)
  [natDegree_C_mul_le]
h' = h + (lower degree term) → natDegree(h') = natDegree(h)
  [Polynomial.natDegree_add_of_natDegree_lt_of_ne_zero]
```
进而：
```
map_{m²}(h') 的 natDegree = natDegree(h')  [natDegree_map_of_leadingCoeff_ne_zero, 因 m² ∤ lc(h')]
= natDegree(h)
= natDegree(map_m h)  [natDegree_map_of_leadingCoeff_ne_zero, 因 m ∤ lc(h)]
```

### 7.4 避免 `set` 的策略

**不使用 `set`**。所有中间量用 `have` 定义或直接内联。
`set` 创建的别名阻碍 `rw`（已知 CLAUDE.md 陷阱）。

**使用 `show` 明确目标类型**，使 `rw` 的模式匹配更可预测。

### 7.5 形式化步骤

| Step | 内容 | 估计行数 |
|------|------|---------|
| 辅助引理 | `div_coeff_poly`（逐系数整除 → ∃ e, f = C(m)*e） | ~20 |
| hensel_step | Steps A-H 组合 | ~60 |
| IsCoprime 传播 | 1+δ unit（ZMod(m²) 方法仍适用） | ~15 |
| hensel_two_factor | 对 n 归纳（doubling） | ~25 |
| hensel_multifactor | 对 factors 长度归纳 | ~30 |
| hensel_correct | 满足 HenselCorrect spec | ~10 |
| **总计** | | **~160 行** |

### 总计：~160 行，预期 0 sorry
