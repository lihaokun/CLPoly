# Hensel 构造性修正：lc-baking + 度数保持

> 状态：nl-proof v2（审核修正：新增 hh_deg 假设 + 度数保持完整推导链）
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

`hensel_step_with_degree` 新增两个假设：
1. `Monic (map_m h)` — 保证在 ZMod m[x] 中可以做 `modByMonic`
2. `natDegree h = natDegree (map_m h)` — h 没有多余高次项（lc 不被 m 整除）

假设 2 的合理性：
- 初始 h 用 canonical lift（ZMod.val 逐系数），natDegree 精确
- 迭代中 h' = h + C(m)·σ_int（deg(σ_int) < natDegree(h)）→ natDegree(h') = natDegree(h) 不变
- 所以假设 2 在每步迭代中自动保持

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
    (hh_monic : Monic (map_m h))           -- 新增：modByMonic 需要
    (hh_deg : h.natDegree = (map_m h).natDegree)  -- 新增：h 无多余高次项
    : ∃ g' h' : Polynomial ℤ,
        map_{m²} f = map_{m²} g' * map_{m²} h'
        ∧ map_m g' = map_m g
        ∧ map_m h' = map_m h
        ∧ (map_{m²} h').natDegree = (map_m h).natDegree  -- 新增 H4
        ∧ h'.natDegree = h.natDegree  -- 新增 H5：迭代时保持 hh_deg
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

**mod m 保持（H3）**：
map_m(h') = map_m(h) + C(0)·map_m(σ_int) = map_m(h)。✓
（因 C(m : ZMod m) = C(0) = 0。）

**度数保持（H4 + H5）**：

h' = h + C(m)·σ_int。

**Step D1**：σ_int 的度数控制。

使用 canonical lift（不用 map_surjective）：
```lean
let σ_int := ∑ i in Finset.range (σ_bar.natDegree + 1),
    C (ZMod.val (σ_bar.coeff i) : ℤ) * X ^ i
```

性质：
- map_m σ_int = σ_bar（每个系数：`(ZMod.val c : ℤ) mod m = c`，by `ZMod.natCast_zmod_val`）
- natDegree(σ_int) ≤ natDegree(σ_bar)（构造只在 range(natDegree+1) 内）
- natDegree(σ_bar) < natDegree(map_m h)（modByMonic 保证）

所以 **natDegree(σ_int) < natDegree(map_m h) = natDegree(h)**（后者由假设 hh_deg）。

**Step D2**：h' 的 natDegree。

h' = h + C(m)·σ_int。natDegree(C(m)·σ_int) ≤ natDegree(σ_int) < natDegree(h)。
加一个严格更低度数的多项式不改变 natDegree 和 leading coefficient：
**natDegree(h') = natDegree(h)**。✓ （H5）
**lc(h') = lc(h)**（未被修改）。

**Step D3**：map_{m²}(h') 的 natDegree。

lc(h') = lc(h)。由 hh_monic：map_m(h) monic → lc(map_m h) = 1 → (lc(h) : ZMod m) = 1。
故 m ∤ lc(h)。进而 m² ∤ lc(h)（因 m | m² 且 m ∤ lc(h)）。
natDegree_map_of_leadingCoeff_ne_zero：(lc(h') : ZMod(m²)) ≠ 0 →
**natDegree(map_{m²} h') = natDegree(h') = natDegree(h) = natDegree(map_m h)**。✓ （H4）

**Lean 路径**：
- canonical lift：`Finset.sum` + `C` + `X^i`（~20 行）
- map_m σ_int = σ_bar：`ZMod.natCast_zmod_val` + `Polynomial.ext`（~10 行）
- natDegree(h') = natDegree(h)：`natDegree_add_of_natDegree_lt` + `natDegree_C_mul_le`（~5 行）
- natDegree(map_{m²} h')：`natDegree_map_of_leadingCoeff_ne_zero`（~10 行）

---

## 3. 形式化估计

| 改动 | 行数 |
|------|------|
| canonical lift 构造 + 验证 | ~30 |
| hensel_step_with_degree（含 modByMonic + 度数） | ~80 |
| hensel_two_factor_with_degree（迭代版） | ~30 |
| **总计** | **~140** |

**不修改现有 hensel_step**（保持向后兼容）。新增 hensel_step_with_degree。
