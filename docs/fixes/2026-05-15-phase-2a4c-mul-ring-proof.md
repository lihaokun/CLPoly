# 修正方案：Phase 2A.4c — SparsePolyZp 乘法 ring 公理证明 + canonical injection

> **状态**：草稿，待用户审核
> 对应 workflow.md §3.1 新功能开发流程
> 分支：`feature/formal-proofs`
> 任务 ID：#75

---

## 第一部分：背景与现状

### Phase 2A.4 拆分
| 子阶段 | 范围 | 状态 |
|--------|------|------|
| 2A.4a | Bridge 类型（Zp.toZMod / WellFormed / toPoly）+ Empty/Zero 证明 + listSum helpers | ✅ 已完成 |
| 2A.4b | `toPoly_add / _neg / _sub` 同态 + `WellFormed_arr.add/sub/neg` 保持 + 加法 ring 公理 + `Zp.toZMod_mul` (单项乘法 helper) | ✅ 已完成 |
| **2A.4c** | **toPoly_mul / WellFormed_arr.mul / mul ring 公理 / toPoly_inj_canonical** | ❌ **未做（本 doc 修复）** |

### 当前文件状态
- `proof/lean/CLPoly/Math/Univariate.lean:558-567` 注释列出 6 项余下工作，**全部未实施**
- `proof/lean/CLPoly/Impl/Types.lean` 零 sorry（含 `toPoly_ne_zero_of_nonempty`）
- TODO.md L313 误标 "[x] 2026-05-07 完成" — 实际仅 add/neg/sub 子集完成，本 doc 完工后才能 [x]

### 现有可复用 infrastructure
- `Zp.toZMod : Nat → Zp → ZMod p`（L40）
- `Zp.toZMod_add : (val toZMod 是同态，前提 2p ≤ UInt64.size)`（L50）
- `Zp.toZMod_neg`（L260）
- `Zp.toZMod_mul`（L283，前提 a.val * b.val < UInt64.size）
- `listSum p : List (UMonomial × Zp) → Polynomial (ZMod p)`（L84-89）
- `listSum_cons / _nil / _append`（L87-94）
- `listSum_mergeAdd`（L147，核心加法）
- `SparsePolyZp.AllReduced p xs := ∀ x ∈ xs, Zp.Reduced p x.snd`（L143）
- `SparsePolyZp.WellFormed_arr p f := AllReduced p f.toList`（L299）
- `SparsePolyZp.toPoly_add / _neg / _sub`（L302-355）
- `SparsePolyZp.WellFormed_arr.add / .neg / .sub`（L471-509）
- 加法 ring 公理（add_comm/_assoc 等，通过 toPoly bridge 已证）

### 待证目标
**Mul 同态链 + canonical injection + mul ring 公理**：

1. `listSum_filterMap_scale` — 标量乘 helper（filterMap 内层）
2. `toPoly_scaleByMonomial` — 单 monomial × poly 同态
3. `WellFormed_arr.scaleByMonomial` — 闭合性
4. `toPoly_mul` — 全乘法同态（核心，最复杂）
5. `WellFormed_arr.mul` — 闭合性
6. `toPoly_inj_canonical` — canonical 唯一对应（独立大块）
7. Ring 公理：`mul_comm / mul_assoc / left_distrib / right_distrib` via toPoly bridge

---

## 第二部分：定义回顾

```lean
-- Model.lean:560
def scaleByMonomial (m : UMonomial) (c : Zp) (f : SparsePolyZp) : SparsePolyZp :=
  if c.val = 0 then #[]
  else f.filterMap (fun (mf, cf) =>
    let prod := c * cf
    if prod.val = 0 then none
    else some (UMonomial.mk (m.deg + mf.deg), prod))

-- Model.lean:569
def mulImpl (f g : SparsePolyZp) : SparsePolyZp :=
  f.foldl (fun acc (mf, cf) => addImpl acc (scaleByMonomial mf cf g)) #[]

-- Model.lean instance L571+:
instance : HMul SparsePolyZp SparsePolyZp SparsePolyZp where
  hMul := SparsePolyZp.mulImpl

-- Math/Univariate.lean:84
def listSum (p : Nat) : List (UMonomial × Zp) → Polynomial (ZMod p)
  | [] => 0
  | (m, c) :: rest => Polynomial.monomial m.deg (Zp.toZMod p c) + listSum p rest
```

---

## 第三部分：nl-proof（按依赖序）

### §3.1 Lemma 1 — `listSum_filterMap_scale`

**陈述**：
```
theorem listSum_filterMap_scale (p : Nat) (h_p2 : p * p ≤ UInt64.size)
    (m : UMonomial) (c : Zp) (hc : Zp.Reduced p c)
    (xs : List (UMonomial × Zp)) (hxs : SparsePolyZp.AllReduced p xs) :
    listSum p (xs.filterMap (fun (mf, cf) =>
      let prod := c * cf
      if prod.val = 0 then none
      else some (UMonomial.mk (m.deg + mf.deg), prod)))
    = Polynomial.C (Zp.toZMod p c) * Polynomial.X^m.deg * listSum p xs
```

**前提**：`p * p ≤ UInt64.size` 防 Mul Zp 溢出（来自 `Zp.toZMod_mul` 前提 `a.val * b.val < UInt64.size`）。

**证明**（归纳 xs）：
- **Base** xs = []：
  - LHS = `listSum p [] = 0`
  - RHS = `C(toZMod c) * X^m.deg * 0 = 0` ✓
- **Step** xs = (mf, cf) :: rest，IH 假设对 rest 成立：
  - 设 `prod := c * cf`
  - **Case 1: `prod.val = 0`**（filterMap drop）：
    - LHS = `listSum p (filterMap _ rest)` = IH 应用 = `C(toZMod c) * X^m.deg * listSum p rest`
    - RHS = `C(toZMod c) * X^m.deg * (monomial mf.deg (toZMod cf) + listSum p rest)`
    - 差值 = `C(toZMod c) * X^m.deg * monomial mf.deg (toZMod cf)`
    - 由 `Zp.toZMod_mul`：`toZMod (c*cf) = toZMod c * toZMod cf`
    - prod.val = 0 → `toZMod prod = 0` → `toZMod c * toZMod cf = 0`
    - `monomial mf.deg (toZMod c * toZMod cf) = monomial mf.deg 0 = 0`
    - `C(toZMod c) * X^m.deg * monomial mf.deg (toZMod cf) = monomial (m.deg + mf.deg) (toZMod c * toZMod cf) = 0`
    - 故差值 = 0 ✓
  - **Case 2: `prod.val ≠ 0`**（filterMap keep）：
    - LHS = `monomial (m.deg + mf.deg) (toZMod prod) + listSum p (filterMap _ rest)`
    - 由 IH：`listSum p (filterMap _ rest) = C(toZMod c) * X^m.deg * listSum p rest`
    - 由 `Zp.toZMod_mul`：`toZMod prod = toZMod c * toZMod cf`
    - `monomial (m.deg + mf.deg) (toZMod c * toZMod cf)`
      `= C(toZMod c) * X^m.deg * monomial mf.deg (toZMod cf)`（Mathlib `monomial_mul_monomial` + `C_mul_X_pow_eq_monomial`）
    - 故 LHS = `C(toZMod c) * X^m.deg * (monomial mf.deg (toZMod cf) + listSum p rest)` = RHS ✓

**Mathlib 依赖**：
- `Polynomial.C_mul_X_pow_eq_monomial : C a * X^n = monomial n a`
- `Polynomial.monomial_mul_monomial : monomial m a * monomial n b = monomial (m+n) (a*b)`
- `Polynomial.monomial_zero_right : monomial n 0 = 0`
- 分配律：`a * (b + c) = a * b + a * c`（`mul_add`）

**预期行数**：~25-30 行

---

### §3.2 Lemma 2 — `toPoly_scaleByMonomial`

**陈述**：
```
theorem SparsePolyZp.toPoly_scaleByMonomial (p : Nat) (h_p2 : p * p ≤ UInt64.size)
    (m : UMonomial) (c : Zp) (hc : Zp.Reduced p c)
    (f : SparsePolyZp) (hf : SparsePolyZp.WellFormed_arr p f) :
    SparsePolyZp.toPoly p (SparsePolyZp.scaleByMonomial m c f)
    = Polynomial.monomial m.deg (Zp.toZMod p c) * SparsePolyZp.toPoly p f
```

**证明**（case split on c.val）：
- **Case 1**: `c.val = 0`：
  - scaleByMonomial 返 `#[]` → toPoly = 0
  - `c.val = 0 → toZMod c = 0`（unfold + `ZMod.natCast_eq_zero_iff`，用 `c.val < p`）
  - RHS = `monomial m.deg 0 * toPoly f = 0 * _ = 0` ✓
- **Case 2**: `c.val ≠ 0`：
  - scaleByMonomial 返 `f.filterMap ...`，转 List
  - 套 Lemma 1：`listSum p (filterMap_body m c f.toList) = C(toZMod c) * X^m.deg * listSum p f.toList`
  - 用 `Polynomial.C_mul_X_pow_eq_monomial` 改写为 `monomial m.deg (toZMod c) * toPoly p f`

**预期行数**：~20 行

---

### §3.3 Lemma 3 — `WellFormed_arr.scaleByMonomial`

**陈述**：
```
theorem SparsePolyZp.WellFormed_arr.scaleByMonomial (p : Nat)
    (m : UMonomial) (c : Zp) (hc : c.prime.toNat = p)
    (f : SparsePolyZp) (hf : SparsePolyZp.WellFormed_arr p f) :
    SparsePolyZp.WellFormed_arr p (SparsePolyZp.scaleByMonomial m c f)
```

**证明**：
- **Case 1**: `c.val = 0` → 返 #[] → WellFormed_arr 平凡成立
- **Case 2**: filterMap 输出每项 `(m+mf, c*cf)`
  - `(c*cf).prime = c.prime`（看 `Mul Zp` 实例：`⟨..., a.prime⟩`）
  - `c.prime.toNat = p` by hypothesis ✓

但注意 `WellFormed_arr` 用的是 `Zp.Reduced`（除 prime 还需 `val < prime`）。`(c*cf).val = c.val * cf.val % c.prime`，自动 `< c.prime = p`。✓

**预期行数**：~15 行

---

### §3.4 Lemma 4 — `toPoly_mul` (核心)

**陈述**：
```
theorem SparsePolyZp.toPoly_mul (p : Nat)
    (h_2p : 2 * p ≤ UInt64.size) (h_p2 : p * p ≤ UInt64.size)
    (f g : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g) :
    SparsePolyZp.toPoly p (f * g)
    = SparsePolyZp.toPoly p f * SparsePolyZp.toPoly p g
```

**证明**：
- 展开 `*` = `mulImpl` = `f.foldl (...) #[]`
- 转 `Array.foldl → List.foldl on f.toList`（Mathlib `Array.foldl_eq_foldl_toList`）
- 定义辅助 `mulPartial xs g acc := xs.foldl (fun acc (mf, cf) => addImpl acc (scaleByMonomial mf cf g)) acc`
- 不变量：对 `processed ++ rest`（其中 acc = mulPartial processed g #[]，rest 未处理），有：
  ```
  toPoly p (mulPartial rest g acc) 
    = toPoly p acc + listSum p rest * toPoly p g
  ```
  当 acc = #[]，rest = f.toList 时：左侧 = toPoly p (f * g)，右侧 = 0 + toPoly p f * toPoly p g ✓
- 归纳 rest：
  - **Base** rest = []：foldl 不变 → toPoly acc + 0 = toPoly acc ✓
  - **Step** rest = (mf, cf) :: rest'：
    - `mulPartial (cons :: rest') g acc = mulPartial rest' g (addImpl acc (scaleByMonomial mf cf g))`
    - 设 `acc' = addImpl acc (scaleByMonomial mf cf g)`
    - 需 `acc'` 仍 WellFormed_arr（由 WellFormed_arr.add + Lemma 3）
    - 由 IH for rest' with acc'：
      ```
      toPoly p (mulPartial rest' g acc')
        = toPoly p acc' + listSum p rest' * toPoly p g
      ```
    - `toPoly p acc' = toPoly p acc + toPoly p (scaleByMonomial mf cf g)`（toPoly_add）
                  `= toPoly p acc + monomial mf.deg (toZMod cf) * toPoly p g`（Lemma 2）
    - listSum p ((mf, cf) :: rest') = monomial mf.deg (toZMod cf) + listSum p rest' (listSum_cons)
    - 故：
      ```
      RHS_after = (toPoly acc + monomial mf.deg (toZMod cf) * toPoly g) 
                  + listSum rest' * toPoly g
                = toPoly acc + (monomial _ + listSum rest') * toPoly g
                = toPoly acc + listSum (cons :: rest') * toPoly g
                = LHS for cons :: rest' ✓
      ```

**Mathlib 依赖**：
- `Array.foldl_eq_foldl_toList`
- `add_mul`（分配律）
- `mul_add`

**预期行数**：~50-60 行（最复杂）

---

### §3.5 Lemma 5 — `WellFormed_arr.mul`

**陈述**：
```
theorem SparsePolyZp.WellFormed_arr.mul (p : Nat) (f g : SparsePolyZp)
    (hf : SparsePolyZp.WellFormed_arr p f) (hg : SparsePolyZp.WellFormed_arr p g) :
    SparsePolyZp.WellFormed_arr p (f * g)
```

**证明**：
- 同 Lemma 4 的 foldl 不变量结构
- 不变量：`WellFormed_arr p acc`
- Base：`#[]` 平凡满足
- Step：用 `WellFormed_arr.add` + `WellFormed_arr.scaleByMonomial`(Lemma 3)
- 注：每个 f 中的项 `(mf, cf)` 都满足 `cf.prime.toNat = p`（由 hf），故 Lemma 3 前提满足

**预期行数**：~15 行

---

### §3.6 Lemma 6 — `toPoly_inj_canonical`

**Canonical 定义**：
```
def SparsePolyZp.Canonical (p : Nat) (f : SparsePolyZp) : Prop :=
  SparsePolyZp.WellFormed_arr p f ∧
  -- 排序：按 deg 严格降序
  List.Chain' (fun a b => a.fst.deg > b.fst.deg) f.toList ∧
  -- 无零系数
  ∀ x ∈ f.toList, x.snd.val ≠ 0
```

**陈述**：
```
theorem SparsePolyZp.toPoly_inj_canonical (p : Nat) [Fact (Nat.Prime p)]
    (f g : SparsePolyZp)
    (hf : SparsePolyZp.Canonical p f) (hg : SparsePolyZp.Canonical p g)
    (heq : SparsePolyZp.toPoly p f = SparsePolyZp.toPoly p g) :
    f = g
```

**证明**：归纳 f.toList.length。
- **Base** `f.toList = []`：
  - toPoly f = 0
  - 假设 g.toList ≠ []，设 (mg₀, cg₀) :: rest = g.toList
  - toPoly g = monomial mg₀.deg (toZMod cg₀) + listSum p rest
  - toPoly g 的 leading coeff（在 Polynomial 中）≠ 0：
    - cg₀.val ≠ 0 (Canonical) + cg₀.val < p (WellFormed) → `toZMod cg₀ : ZMod p` 非零（用 `ZMod.natCast_zmod_eq_zero_iff_dvd` + `Nat.Prime` + val < p → ¬(p ∣ val)）
    - rest 所有 deg < mg₀.deg（Chain' 性质递推），故 listSum p rest 在 X^mg₀.deg 系数 = 0
    - toPoly g 在 X^mg₀.deg 系数 = toZMod cg₀ ≠ 0
  - 故 toPoly g ≠ 0，与 heq 矛盾
  - 故 g.toList = []，f = g（空数组等）
- **Step** `f.toList = (mf₀, cf₀) :: rest_f`，类似 g：
  - 同理 g.toList ≠ []，设 = (mg₀, cg₀) :: rest_g
  - 比较两端 leading degree：
    - deg(toPoly f) = mf₀.deg（同上推理）
    - deg(toPoly g) = mg₀.deg
    - toPoly f = toPoly g → mf₀.deg = mg₀.deg
  - 比较 leading coef：
    - toPoly f.coeff mf₀.deg = toZMod cf₀
    - toPoly g.coeff mg₀.deg = toZMod cg₀
    - 由 mf₀.deg = mg₀.deg 和 heq → `toZMod cf₀ = toZMod cg₀`
    - toZMod injective on Reduced：`cf₀.val < p, cg₀.val < p` → `cf₀.val = cg₀.val`
    - cf₀.prime = cg₀.prime = p（WellFormed）→ `cf₀ = cg₀`
  - 故 head 相等：`(mf₀, cf₀) = (mg₀, cg₀)`
  - 余下：`monomial(_, toZMod cf₀) + listSum rest_f = monomial(_, toZMod cg₀) + listSum rest_g`
    → `listSum rest_f = listSum rest_g`
  - rest_f, rest_g 仍 Canonical（Chain' 头尾分解）
  - 由 IH（长度严格递减）：rest_f.toArray = rest_g.toArray
  - 故 f = g

**Mathlib 依赖**：
- `Polynomial.degree / leadingCoeff / coeff_monomial`
- `Polynomial.coeff_add`
- `ZMod.natCast_zmod_eq_zero_iff_dvd : (n : ZMod p) = 0 ↔ p ∣ n`
- `Nat.Prime.dvd_iff_eq`
- `List.Chain'` 性质

**预期行数**：~100-120 行（最难，多 lemma 拆分）

---

### §3.6.5 — Canonical preservation 子定理 nl（2026-05-15 补）

`§3.7 Ring 公理` 假定 `Canonical.add / .mul` 保持性，本节提供详细 nl-proof。

#### Lemma 7.1 — `chain_R_head_to_mem` (generic helper)

**陈述**：若 `R` 传递 + `xs` 满足 `IsChain R` + 对 `xs.head?` 有 `∀ y ∈ xs.head?, R x y`，
则 `∀ a ∈ xs, R x a`。

**证明**：归纳 xs。
- nil：vacuous
- cons z zs：a = z 时由 hypothesis 直接；a ∈ zs 时由 IH on (R x zs.head?) — 需 `R x z` + `R z zs.head` → `R x zs.head` (R 传递)

#### Lemma 7.2 — `List.IsChain_filterMap_of_mono` (generic helper)

**陈述**：若 `R, S` 传递 + 单调 `f : α → Option β`（即 `R a₁ a₂ + f a_i = some b_i → S b₁ b₂`），
则 `xs.filterMap f` 满足 `IsChain S`。

**证明**：归纳 xs。
- nil：filterMap = []，trivial
- cons x rest：case split on `f x`：
  - `f x = none`：filterMap 跳过，直接用 IH
  - `f x = some bx`：用 IH 得 `rest.filterMap f` 是 IsChain S；需证 head 关系
    - 设 y 是 `rest.filterMap f` 的首项，要 `S bx y`
    - 由 `List.mem_filterMap` 拆解：∃ a ∈ rest, f a = some y
    - 由 Lemma 7.1（R 传递）+ `R x rest.head?`：`R x a`
    - 由 mono：`S bx y` ✓

#### Lemma 7.3 — `scaleByMonomial_Chain'`

**陈述**：`Canonical p f → IsChain (>) ((scaleByMonomial m c f).toList.map deg)`
等价：`Canonical p f → scaleByMonomial m c f` 的 deg 列表严格降序。

**证明**：
- case `c.val = 0`：返 #[]，trivial Chain'
- case `c.val ≠ 0`：scaleByMonomial = filterMap-over-orig
  - 套 Lemma 7.2 with：
    - R = `(>)`on orig 的 deg, S = `(>)` on new 的 deg
    - R 传递：Nat 上 `>` 当然传递
    - mono：原 `mf₁.deg > mf₂.deg` → 新 `(m.deg + mf₁.deg) > (m.deg + mf₂.deg)` (加法保单调)

#### Lemma 7.4 — `mergeAdd_Chain'`

**陈述**：若 xs 和 ys 都满足 IsChain (>) 且 AllReduced (Reduced)，
则 `mergeAdd xs ys` 满足 IsChain (>)。

**证明**：归纳 `mergeAdd.induct` 的 6 case：
- case1 (xs=[])：mergeAdd = ys，由前提 ys Chain'
- case2 (ys=[])：对称
- case3 (xs head > ys head)：合并 = xs.head :: mergeAdd xs.rest ys
  - 用 IH on xs.rest, ys 得 mergeAdd 结果 Chain'
  - 需要 head 关系：xs.head.deg 大于 mergeAdd 结果 head 的 deg
  - mergeAdd 结果 head 是 max(xs.rest.head, ys.head)
  - xs.head.deg > xs.rest.head.deg (xs 自身 Chain')
  - xs.head.deg > ys.head.deg (case 假设)
  - 故 xs.head.deg > max(...) ✓
- case4 (xs.head < ys.head)：对称
- case5 (deg 相等, sum.val = 0)：返 mergeAdd xs.rest ys.rest，IH 直接
- case6 (deg 相等, sum.val ≠ 0)：合并 = (xs.head.deg, sum) :: mergeAdd xs.rest ys.rest
  - 用 IH on xs.rest, ys.rest
  - head 关系：xs.head.deg 大于 mergeAdd xs.rest ys.rest 的 head deg
  - 同 case3 推理

预期行数：~80（与 listSum_mergeAdd 同复杂度）

#### Lemma 7.5 — `Canonical.add`

**陈述**：`Canonical p f ∧ Canonical p g → Canonical p (f + g)`

**证明**：拆 3 子句：
- WellFormed：已证 (Stage 2A.4b)
- Chain'：用 Lemma 7.4 mergeAdd_Chain'
- no-zero：mergeAdd 内置 filter val = 0（case5 删除）

#### Lemma 7.6 — `Canonical.mul`

**陈述**：`Canonical p f ∧ Canonical p g + h_p2 → Canonical p (f * g)`

**证明**：f * g = foldl over (addImpl acc (scaleByMonomial mf cf g))
- foldl 不变量：每步 acc Canonical
- 初值 `#[]` Canonical（vacuously）
- 每步：addImpl (Canonical, Canonical) → Canonical (Lemma 7.5 + Lemma 7.3)

预期行数：~25

### §3.7 Ring 公理

**陈述模式**：
```
theorem SparsePolyZp.mul_comm_via_toPoly (p : Nat) [Fact (Nat.Prime p)]
    (h_2p : 2 * p ≤ UInt64.size) (h_p2 : p * p ≤ UInt64.size)
    (f g : SparsePolyZp)
    (hf : SparsePolyZp.Canonical p f) (hg : SparsePolyZp.Canonical p g)
    (hfg_canon : SparsePolyZp.Canonical p (f * g))
    (hgf_canon : SparsePolyZp.Canonical p (g * f)) :
    f * g = g * f
```

**证明**：
- `toPoly p (f * g) = toPoly p f * toPoly p g`（Lemma 4）
- `= toPoly p g * toPoly p f`（Mathlib `mul_comm` on Polynomial）
- `= toPoly p (g * f)`（Lemma 4 反向）
- 由 `toPoly_inj_canonical` 配合 `hfg_canon, hgf_canon`：`f * g = g * f` ✓

**等价模式**：`mul_assoc`, `left_distrib`, `right_distrib`。

**关于 Canonical 假设**：需另证 `mulImpl` 保持 Canonical（不仅 WellFormed）。`mulImpl` 经过 `addImpl` 由 mergeAdd 自动产生 Canonical 结构（降序、无重复 deg、无零）。这需要 Lemma `Canonical.mul`。

**预期行数**：4 个 ring 公理 × ~15 行 + Canonical.mul ~30 行 = ~90 行

---

## 第三-bis 部分：实施状态（2026-05-15 更新）

### 已完成
- ✅ **Lemma 1-3** Stage 1 (commit `5da715b`) — 135 行
- ✅ **Lemma 4-5** Stage 2 (commit `a445719`) — 106 行
- ✅ **Ring 公理 (toPoly-level)** Stage 3a (commit `ac84154`) — 70 行
- ✅ **Lemma 6 toPoly_inj_canonical 核心** Stage 3b core (commit `1a983d2`) — 243 行
  - 配套 helpers: listSum_coeff_zero_of_all_lt, Zp.toZMod_inj_of_reduced,
    chain_gt_all_after_head, listSum_coeff_at_head, listSum_inj_canonical

**当前 lake build 0 sorry，3071 jobs 全过**。

### 未完成（nl-proof 已补足 §3.6.5；Lean 形式化为 Stage 3b 续 / 独立 Phase 2A.4d）
- ❌ Lemma 7.1-7.2 generic chain helpers (~30 行)
- ❌ Lemma 7.3 scaleByMonomial_Chain' (~15 行)
- ❌ Lemma 7.4 mergeAdd_Chain' (~80 行；与 listSum_mergeAdd 同复杂度)
- ❌ Lemma 7.5-7.6 Canonical.add / .mul (~30 行)
- ❌ §3.7 Ring 公理 (4 个，~80 行)

预计剩余形式化工作量：~250 行 Lean 代码 + 多轮 Mathlib API 探索（已知卡点：`List.IsChain` deprecated form, `List.mem_filterMap`, `head?` 处理等）

## 第四部分：实施分阶段

按依赖序，分 3 个 commit：

### Stage 1: scaleByMonomial 同态 + WellFormed
- Lemma 1 (`listSum_filterMap_scale`)
- Lemma 2 (`toPoly_scaleByMonomial`)
- Lemma 3 (`WellFormed_arr.scaleByMonomial`)
- 预期：~70 行新代码，1-2 小时

### Stage 2: toPoly_mul + WellFormed_arr.mul
- Lemma 4 (`toPoly_mul`)（核心）
- Lemma 5 (`WellFormed_arr.mul`)
- 预期：~75 行新代码，2-3 小时

### Stage 3: Canonical injection + ring 公理
- `SparsePolyZp.Canonical` 定义
- Lemma 6 (`toPoly_inj_canonical`)
- `Canonical.add / .mul` 保持性
- Ring 公理 4 个（mul_comm / mul_assoc / left_distrib / right_distrib）
- 预期：~150 行新代码，3-4 小时

---

## 第五部分：测试计划

### 编译验证
- 每 Stage 结束 `lake build` 全过（3071 jobs，加上新增 1-2 jobs）
- 0 sorry 验证：`grep "sorry" proof/lean/CLPoly/Math/Univariate.lean`

### 不变性验证
- 现有 137 B2B 测试不退化
- `__mul_zp_poly` B2B 测试已 PASS（4 case），新增 Lemma 仅形式化已有行为，不应改变 `#eval` 输出

### Mathlib API 调研
- 实施前每 Stage 先 grep 确认所有引用 Mathlib lemma 名 & 签名正确（参 CLAUDE.md §"nl-proof 审核标准"）

---

## 第六部分：风险与未决

### R1: Mathlib `monomial_mul_monomial` API 是否完全匹配
- 用 `Polynomial.monomial_mul_monomial` 还是 `Polynomial.monomial_mul`？
- 待 Stage 1 实施前 grep 验证

### R2: `toPoly_inj_canonical` 的归纳变量
- 用 `f.toList.length` 还是 `f.size`？
- 用 `Nat.strongRecOn` 还是直接 List.rec？
- 待 Stage 3 实施前判断

### R3: Canonical 定义是否需要扩展现有 WellFormed
- `WellFormed_arr` 仅含 AllReduced
- Canonical 额外加 Chain' (严格降序) + no-zero-val
- 是否要把 Canonical 提升为 default predicate？
- 决策：保持 WellFormed_arr 作为 Lemma 1-5 的前提，Canonical 仅 Lemma 6 + ring 公理用

### R4: ring 公理"实用 form"问题
- 若 f, g 不是 Canonical，f * g = g * f 在 Array 层可能不成立（中间 acc 顺序不同）
- 解决：先用 `normalization`（已有 def L113）规范化，再比较 — 这会让 ring 公理表述变成"normalize ∘ mul 满足 ring 公理"
- 待 Stage 3 决定最终陈述形式

### R5: 不打算证的
- `WellFormed_arr` 在 mergeAdd 输出上的 Chain' 性质（即 mergeAdd 保持降序）—— 是 Canonical.add 的子引理，可能本身需 ~30 行
- 范围内若发现独立工作量很大，移到 follow-up

---

## 第七部分：执行步骤（待批准）

1. **[等批准]** 用户审核本 doc
2. **Stage 1 实施**：Lemma 1-3 + 编译 + 单元 #eval 检查（若有）+ commit
3. **Stage 2 实施**：Lemma 4-5 + commit
4. **Stage 3 实施**：Canonical 定义 + Lemma 6 + ring 公理 + commit
5. 三个 Stage 间用户都可中断 / 审视；每 Stage 完成后 lake build 全过 + commit + push
6. 完成后：
   - 更新 TODO.md L313 状态（移除"2026-05-07 完成"误标，改成"2026-05-15 完成"）
   - 更新 task #75 status → completed
   - 更新 Univariate.lean 末尾注释（移除"余下工作"列表）

---

## 第八部分：依赖图（一图速查）

```
                    Zp.toZMod_mul (已证)
                          │
                    Lemma 1 (listSum_filterMap_scale)
                    /        \
        Lemma 2 (toPoly_scaleByMonomial)    Lemma 3 (WellFormed.scaleByMonomial)
                    \        /
                     \      /
                      ↓    ↓
                    Lemma 4 (toPoly_mul)  ← + toPoly_add (已证) + addImpl
                          │
                    Lemma 5 (WellFormed.mul) ← + WellFormed.add (已证)
                          │
                          ↓
                ┌────────────────┐
                │ Canonical 定义  │
                └────────────────┘
                          │
                    Lemma 6 (toPoly_inj_canonical)
                    │
                    │      Canonical.add (新)
                    │      Canonical.mul (新)
                    ↓
                Ring 公理（mul_comm/_assoc/left_distrib/right_distrib）
```
