# FactorZp.lean 自然语言证明

> 状态：全部已形式化，0 sorry

---

## 辅助引理

### 1. poly_unit_eq_C

**定理**：域上多项式环的单位是非零常数，即 u 是单位 → ∃ a ≠ 0, ↑u = C a。

**证明**：
1. Mathlib 提供 `natDegree_coe_units u`：单位的 natDegree = 0
2. `eq_C_of_natDegree_eq_zero`：natDegree = 0 → 多项式 = C (coeff 0)
3. `coeff_coe_units_zero_ne_zero u`：单位的常数项 ≠ 0
4. 组合即得。

### 2. eq_of_associated_monic

**定理**：首一多项式的 Associated 即相等，即 Monic a → Monic b → Associated a b → a = b。

**证明**：
1. Associated a b 给出单位 u 使得 a * ↑u = b
2. 由 `natDegree_coe_units`，↑u 的 natDegree = 0，所以 ↑u = C c
3. 比较两边 leadingCoeff：lc(a) * c = lc(b)，即 1 * c = 1，得 c = 1
4. 代入：a * C 1 = a * 1 = a = b

**关键依赖**：`natDegree_coe_units`、`eq_C_of_natDegree_eq_zero`、`leadingCoeff_mul`（需要 NoZeroDivisors，ZMod p 是域所以满足）

### 3. monic_list_prod

**定理**：首一多项式列表的积是首一。

**证明**：列表归纳。
- nil：空积 = 1，Monic 1 成立
- cons a t：积 = a * t.prod，由 Monic.mul 和归纳假设得证

### 4. list_prod_pow

**定理**：(∏ aᵢ)^n = ∏(aᵢ^n)，在交换幺半群中成立。

**证明**：列表归纳。
- nil：1^n = 1 = 空积
- cons a t：(a * t.prod)^n = a^n * (t.prod)^n = a^n * ∏(tᵢ^n)（由 mul_pow + 归纳假设）

### 5. list_prod_flatMap

**定理**：(l.flatMap f).prod = (l.map (fun x => (f x).prod)).prod

**证明**：列表归纳。
- nil：[].flatMap f = []，积 = 1；[].map ... = []，积 = 1
- cons a t：(a :: t).flatMap f = f(a) ++ t.flatMap f
  - 积 = f(a).prod * (t.flatMap f).prod（List.prod_append）
  - = f(a).prod * (t.map ...).prod（归纳假设）
  - = ((a :: t).map ...).prod（List.prod_cons）

---

## ddf_edf_combine

**定理**：对单个首一无平方多项式 g，DDF + EDF 组合给出完整不可约分解。

**证明**：

**第一部分（乘积相等）**：
1. 对每个 (gd, d) ∈ ddf(g)：
   - gd 是首一的（DDFCorrect.5）
   - gd 是无平方的（Squarefree.squarefree_of_dvd，因 gd ∣ g 且 g 无平方）
   - EDF 给出 EDFCorrect gd d (edf gd d)
   - Associated gd (edf gd d).prod，两边都首一 → 由 eq_of_associated_monic 得 gd = (edf gd d).prod
2. g 与 (ddf g).map(Prod.fst).prod 是 Associated 的（DDFCorrect.4），两边首一 → 相等
3. 逐项替换 Prod.fst → (edf ...).prod，得 g = flatMap 的积

**第二部分（每个因子不可约且首一）**：
- q ∈ flatMap → ∃ (gd, d) ∈ ddf(g) 使得 q ∈ edf gd d
- EDFCorrect 保证 q 不可约且首一

---

## 主定理 factor_Zp_correct

**定理**：假设 SQF、DDF、EDF 各自正确，则组合结果满足 FactorZpCorrect。

**证明**：

1. **SQF 展开**：f ≈ ∏ gᵢ^eᵢ，每个 gᵢ 首一无平方，eᵢ ≥ 1
2. **对每个 SQF 分量**：用 ddf_edf_combine 得到 gᵢ = ∏ qₖ（不可约首一因子）
3. **构造 factors**：对每个 (g, e) ∈ sqf(f)，展开为 [(q, e) | q ∈ DDF+EDF(g)]
4. **key 证明**：factors.map(·^·).prod = sqf(f).map(·^·).prod
   - 展开 factors 定义，map 穿过 flatMap（List.map_flatMap + List.map_map）
   - 外层 list_prod_flatMap 把 flatMap 积转为 map-prod 积
   - 逐项：内层积 = gᵢ^eᵢ（由 component_ok 得 gᵢ = ∏qₖ，再 list_prod_pow 得 (∏qₖ)^eᵢ = ∏(qₖ^eᵢ)）
   - **注意**：必须用 conv_rhs 定向重写，否则 rw 会替换两边所有 ge.1
5. **Associated → 单位提取**：sqf_assoc + key → Associated f factors_prod → ∃ u, f * u = prod
6. **Step 6 单位代数**：f = f·1 = f·(u·u⁻¹) = (f·u)·u⁻¹ = prod·C(a) = C(a)·prod
7. **条件二**：每个 (q, e) ∈ factors，q 来自 EDF → Irreducible + Monic，e 来自 SQF → ≥ 1
