# Hensel 构造性修正：lc-baking + 度数保持

> 状态：nl-proof v1
> 目标：修改 hensel_step 使其输出 (H4) 度数保持，匹配 C++ __hensel_step

---

## 0. 当前差距

当前 `hensel_step` 构造 `h' = h + C(m)·(s·e_int)`，`s·e_int` 可能有任意高度数 → h' 度数不受控。

C++ 的 `__hensel_step`（line 431-453）：
```cpp
se = s * e;
(q_se, r_se) = divmod(se, h, mod m);  // r_se = se mod h (mod m), deg(r_se) < deg(h)
tau = t*e + q_se*g;
g_new = g + m * tau;
h_new = h + m * r_se;   // deg(r_se) < deg(h) → deg(h_new) = deg(h)
```

关键：C++ 对 `s*e` 做 divmod by `h`，取余式 `r_se` 作为 h 的修正项。`deg(r_se) < deg(h)` 保证度数不增。

---

## 1. 修改方案

### 1.1 新增假设

`hensel_step` 新增 `Monic (map_m h)` 假设。这保证在 ZMod m[x] 中可以做 `modByMonic`。

### 1.2 新构造

在 ZMod m[x] 中做 Bézout + divmod：

```
-- 在 ZMod m[x] 中
s_bar := map_m s
e_bar := map_m e_int    -- 注：map_m e_int = 0（因 e_int = (f-g*h)/m）
```

等等—— `map_m e_int = 0` 因为 `f - g*h = C(m)·e_int` 所以 `map_m(f - g*h) = 0` 但 `map_m(C(m)·e_int) = C(0)·map_m(e_int) = 0`，这不给出 `map_m(e_int)` 的信息。实际上 `C(m)·e_int = m·e_int` 的每个系数是 m 的倍数乘以 e_int 的系数... 不对，`C(m)·e_int` 的第 k 个系数 = m · (e_int 的第 k 个系数)。`map_m` 给零只因为 m ≡ 0 (mod m)。e_int 本身的系数不需要被 m 整除。

所以 `map_m(e_int)` 不一定是零——e_int 是任意 Z[x] 多项式。

**核心观察**：C++ 在 `ZMod m[x]` 中做 divmod `s*e mod h`。在我们的 ℤ[x] 方法中：
- s, e_int ∈ ℤ[x]
- 在 ZMod m[x] 中：σ_bar = (map_m s · map_m e_int) modByMonic (map_m h)
- deg(σ_bar) < deg(map_m h)
- 提升 σ_bar 回 ℤ[x]：∃ σ_int with map_m σ_int = σ_bar

然后 h' = h + C(m)·σ_int（而不是 h + C(m)·s·e_int）。

### 1.3 问题：lift σ_bar 的度数控制

`Polynomial.map_surjective` 给出**某个** σ_int with map_m σ_int = σ_bar。但 σ_int 的 natDegree 可能任意大（高次系数可以是 m 的倍数，不影响 mod m 像）。

**解决**：不需要 σ_int 在 ℤ[x] 中度数受控。只需要 `map_{m²} h'` 的度数受控。

`h' = h + C(m)·σ_int`。
`map_{m²} h' = map_{m²} h + C(m_m2)·map_{m²} σ_int`。

C(m_m2) = C(m : ZMod(m²))。m 在 ZMod(m²) 中不是零（m < m²），但 m 也不是 unit。

`C(m)·map_{m²} σ_int` 的度数 ≤ deg(map_{m²} σ_int)。
`map_{m²} σ_int` 的度数可能很大...

**这又回到了同样的问题**。即使在 ZMod m 中 σ_bar 度数受控，lift 到 ℤ 或 ZMod(m²) 后度数可能增大。

### 1.4 正确方案：在乘积验证中利用度数

实际上，对于匹配 C++，我们不需要在 hensel_step 中证度数保持。

**关键洞察**：HenselCorrect 比较的是 ZMod(p^k) 多项式的度数。`hensel_two_factor` 输出的 g', h' ∈ ℤ[x]，它们的 `map_{p^k}` 像的度数是什么？

由 Hensel 唯一性（hensel_unique）：map_{p^k} g' = map_{p^k} (唯一的 Hensel lift)。
唯一 Hensel lift 的 ZMod(p^k) 像有正确的度数（因为真因子 gⱼ 的 map_{p^k} 有正确度数）。
所以 natDegree(map_{p^k} g') = natDegree(map_{p^k} gⱼ 的调整) = natDegree(map_p gⱼ)。

**这个论证不需要修改 hensel_step！只需要 Hensel 唯一性 + 真因子存在。**

### 1.5 结论

**不修改 hensel_step 的构造**。度数保持通过 Hensel 唯一性推导：

```lean
theorem hensel_degree_preservation
    (p k : ℕ) (hp : Nat.Prime p) (hk : 0 < k)
    (f g h g' h' : Polynomial ℤ)
    (hprod : ...)  -- hensel_two_factor 给的 3 个性质
    (hg_monic : Monic (map_p g))  -- map_p g monic
    (hh_monic : Monic (map_p h))  -- map_p h monic
    : natDegree (map_{p^k} h') = natDegree (map_p h)
```

证明：
1. map_p h' = map_p h（hensel_two_factor H3）
2. map_p h monic → 系数位置 natDegree(map_p h) 处 h' 的系数 ≡ 1 (mod p)
3. → p ∤ 该系数 → p^k ∤ 该系数
4. → natDegree(map_{p^k} h') ≥ natDegree(map_p h)
5. castHom 不增度数 → natDegree(map_p h') ≤ natDegree(map_{p^k} h')
6. 但 natDegree(map_p h') = natDegree(map_p h)
7. 所以 natDegree(map_{p^k} h') ≥ natDegree(map_p h)

还需要上界... 这需要 h' 在位置 > natDegree(map_p h) 的系数都被 p^k 整除。
这等价于 natDegree(map_{p^k} h') ≤ natDegree(map_p h)。

但 h' = h + C(m)·(s·e_int)，位置 > natDegree(h) 的系数 = m·(s·e_int 的对应系数)。
这些系数被 m 整除但不一定被 m² = p^{2k} 整除... 所以 map_{m²} h' 可能度数更大。

**但 map_{p^k} h' 通过 castHom 从 map_{m²} h' 投射得到**... natDegree 可能变小也可能不变。

实际上：如果 map_{p^k} h' 的某个高次系数不为零（在 ZMod(p^k) 中），那就度数更大。
而 h' 的高次系数 = m·(s·e_int 在该位置的系数)。m = p^j 对某个 j。
如果 j < k：这个系数在 ZMod(p^k) 中不为零 → 度数增大 → 度数保持失败。

**结论**：不改构造确实不能保证度数保持。之前的分析是对的。

---

## 2. 最终方案

**必须修改构造**来控制度数。但不修改现有 `hensel_step`（它的 3 个输出性质仍然有用）。而是**新增一个带度数保持的版本**：

```lean
theorem hensel_step_with_degree
    (m : ℕ) (hm : 1 < m)
    (f g h : Polynomial ℤ)
    (hprod : map_m f = map_m g * map_m h)
    (hcop : IsCoprime (map_m g) (map_m h))
    (hh_monic : Monic (map_m h))  -- 新增
    : ∃ g' h' : Polynomial ℤ,
        map_{m²} f = map_{m²} g' * map_{m²} h'
        ∧ map_m g' = map_m g
        ∧ map_m h' = map_m h
        ∧ natDegree (map_{m²} h') = natDegree (map_m h)  -- 新增 H4
```

**构造**（匹配 C++）：

Step A-B 同现有 hensel_step。

Step C（新）：在 ZMod m[x] 中做 modByMonic。
```
σ_bar := (map_m s · map_m e_int) %ₘ (map_m h)    -- deg(σ_bar) < deg(map_m h)
q_bar := (map_m s · map_m e_int) /ₘ (map_m h)
τ_bar := map_m t · map_m e_int + q_bar · map_m g
```

Step D：提升回 ℤ[x]。
```
∃ σ_int with map_m σ_int = σ_bar    (map_surjective)
∃ τ_int with map_m τ_int = τ_bar    (map_surjective)
```

Step E：定义 g', h'。
```
g' := g + C(m) · τ_int
h' := h + C(m) · σ_int
```

**乘积验证**（同 nl-proof phase3-t34-hensel.md）：
g'·h' - f = C(m)·(σ_int·g + τ_int·h - e_int) + C(m²)·σ_int·τ_int

map_m(σ_int·g + τ_int·h) = σ_bar·map_m(g) + τ_bar·map_m(h) = map_m(e_int)
（Bézout 分解）。所以 map_m(σ_int·g + τ_int·h - e_int) = 0。
→ ∃ d, σ_int·g + τ_int·h - e_int = C(m)·d。
→ C(m)·(σ_int·g + τ_int·h - e_int) = C(m²)·d。
→ g'·h' - f = C(m²)·(d + σ_int·τ_int)。
→ map_{m²}(g'·h' - f) = 0。✓

**度数保持**：
map_m(h') = map_m(h) + C(0)·σ_bar = map_m(h)。✓
（因 C(m : ZMod m) = C(0) = 0。）

对于 map_{m²} h' 的度数：
h' = h + C(m)·σ_int。
在 ZMod(m²) 中：map_{m²}(h') = map_{m²}(h) + C(m_m2)·map_{m²}(σ_int)。

deg(C(m_m2)·map_{m²}(σ_int)) ≤ deg(map_{m²} σ_int) ≤ deg(σ_int)。

**关键**：σ_int 是 σ_bar 的 ℤ[x] 提升。σ_bar = (map_m s · map_m e_int) %ₘ (map_m h)。
deg(σ_bar) < deg(map_m h)。

如果 σ_int 的度数 = deg(σ_bar)（精确提升）：
deg(C(m_m2)·map_{m²}(σ_int)) ≤ deg(σ_bar) < deg(map_m h) = deg(map_{m²} h)
（后者因 h monic mod m → p ∤ lc(h) → natDegree_map 保持）。

所以 map_{m²}(h') = map_{m²}(h) + 低次项 → natDegree(map_{m²} h') = natDegree(map_{m²} h) = natDegree(map_m h)。✓

**但**：`Polynomial.map_surjective` 不保证 deg(σ_int) = deg(σ_bar)。
如何得到精确度数的提升？

**方案**：不用 map_surjective。构造 σ_int 为 σ_bar 的 canonical lift（用 ZMod.val 逐系数提升）：

```lean
-- σ_int.coeff i = ZMod.val (σ_bar.coeff i) for i < deg(σ_bar), 0 otherwise
```

这保证 deg(σ_int) ≤ deg(σ_bar) < deg(map_m h)。✓

**Lean 构造**（~20 行）：
```lean
let σ_int := ∑ i in Finset.range (σ_bar.natDegree + 1),
    C (ZMod.val (σ_bar.coeff i) : ℤ) * X ^ i
```

验证 map_m σ_int = σ_bar（用 ZMod.natCast_zmod_val）。

---

## 3. 形式化估计

| 改动 | 行数 |
|------|------|
| canonical lift 构造 + 验证 | ~30 |
| hensel_step_with_degree（含 modByMonic + 度数） | ~80 |
| hensel_two_factor_with_degree（迭代版） | ~30 |
| **总计** | **~140** |

**不修改现有 hensel_step**（保持向后兼容）。新增 hensel_step_with_degree。
