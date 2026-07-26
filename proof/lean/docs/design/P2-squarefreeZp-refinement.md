# P2-a：SquarefreeZp 精化设计文档（`__squarefree_Zp_ir_refines` Branch A + B）

日期：2026-07-26
目标：清 `CLPoly/Refinement/SquarefreeZp.lean` 的 2 个 admit：
- **Admit 1**（Branch A，line 1969，导数=0 分支）—— 指数 UInt64 转换。
- **Admit 2**（Branch B，line 1986，导数≠0 即 Yun 分支）—— 整个 Yun 循环对应。

定理签名（line 1594）：
```
__squarefree_Zp_ir_refines p f (hwf_f) (hred_f) (hp_size : 2*p ≤ 2^64)
  (h_no_overflow : ∀ x ∈ f.toList, x.2.val.toNat * x.1.deg < 2^64)
  (h_deg_bound   : ∀ x ∈ f.toList, x.1.deg < 2^64)
  (h_val_nonzero_f : ∀ x ∈ f.toList, x.2.val.toNat ≠ 0)
  : toPolyList (__squarefree_Zp_ir_safe p f) p = sqfZp (toPoly p f)
```
强归纳于 `natDegree`。已完成：Branch A 的 95%（导数映射、g_1/g_2 构造、递归 IH h_ih_g2、
循环正确性 h_loop_lemma、单位乘不变 h_sqfz_smul、分解同构 h_sqfz_g）。Branch B 未开始。

---

## Part 1：Branch A（Admit 1，line 1969）

### 1.1 目标命题
```
h_toPolyList_specific :
  toPolyList ((__squarefree_Zp_ir_safe p g_2).map (fun (g_h, e) => (g_h, e * p_1))) p
    = (toPolyList (__squarefree_Zp_ir_safe p g_2) p).map (fun (h, e) => (h, e * p))
```
`simp [Function.comp, hp_1_eq_p]` 后归结为**逐元素**：`(e * UInt64.ofNat p).toNat = e.toNat * p`
（`e : UInt64` 为 C++ 数组 `__squarefree_Zp_ir_safe p g_2` 中每个因子的重数；p 侧 `e` 已转 ℕ）。

### 1.2 UInt64 转换引理（✓ 已验证可证）
```lean
lemma uint64_mul_ofNat_toNat (e : UInt64) (b : ℕ) (hb : b < 2^64) (h : e.toNat * b < 2^64) :
    (e * UInt64.ofNat b).toNat = e.toNat * b := by
  rw [UInt64.toNat_mul, UInt64.toNat_ofNat_of_lt' hb, Nat.mod_eq_of_lt h]
```
`p < 2^64` 由 `hp_size : 2*p ≤ 2^64` 得（p < 2^63）。**唯一缺口**：`e.toNat * p < 2^64`（指数界）。

### 1.3 指数界（关键新引理 —— Branch A 的真正难点）

**需证**：`∀ e ∈ (__squarefree_Zp_ir_safe p g_2 的重数), e.toNat * p < 2^64`。

**推导链**（每步引理）：
1. `e.toNat ≤ natDegree(toPoly p g_2)` —— **新结构引理**（见下），核心。
2. `natDegree(toPoly p g_2) = natDegree(g_1) = natDegree(contract p (toPoly p g))`
   —— 由已有 `g_2 = makeMonic(g_1)`（natDegree 不变，h_makeMonic 单位乘）+ `g_1 = extract_pth_root`。
3. `natDegree(contract p F) ≤ natDegree(F) / p` —— 已有 `natDegree_contract_le`（line 477-486）。
   故 `p · natDegree(contract p F) ≤ natDegree(F)`（p ≥ 1）。
4. `natDegree(toPoly p g) = natDegree(g) < 2^64` —— 由 `h_deg_bound`（首项度 < 2^64）+
   `natDegree = 首项.deg`（稀疏降序表示）。
5. 合：`e.toNat · p ≤ natDegree(toPoly g_2) · p ≤ natDegree(contract p (toPoly g)) · p
        ≤ natDegree(toPoly g) < 2^64`。

**新结构引理 `sqfZp_exponent_le_natDegree`**：
```lean
lemma sqfZp_exponent_le_natDegree (f : Polynomial (ZMod p)) :
    ∀ (he : (h, e) ∈ sqfZp f), e ≤ f.natDegree
```
证：`sqf_correct f` ⇒ `f ~ ∏ (hᵢ^eᵢ)`（乘积结构），每 `hᵢ` 非常数（deg ≥ 1）⇒
`Σ eᵢ·deg(hᵢ) = deg(f)` ⇒ `eᵢ ≤ Σ eⱼ·deg(hⱼ) = deg(f)`。
Mathlib 路径：`SquarefreeDecomp` 的乘积字段 + `natDegree_prod` + `Finset.single_le_sum`。

**从 sqfZp 指数界到 C++ 原始 UInt64 e**：C++ 数组 `__squarefree_Zp_ir_safe p g_2` 的 `e.toNat`
= `sqfZp (toPoly p g_2)` 的对应指数（由 `h_ih_g2 : toPolyList (…g_2) p = sqfZp(toPoly g_2)`
逐元素对应；toPolyList 保重数的 toNat）。故需一条「toPolyList 保指数 toNat」+ `h_ih_g2` 的
成员对应引理，或直接对 C++ 输出证结构不变量 `∀ e ∈ __squarefree_Zp_ir output, e.toNat ≤ natDegree`
（对 C++ 递归归纳，复用 h_loop_lemma）。

> **实现注记**：1.3 是 Branch A 的实际工作量所在（非 1-2 小时）。两条路径二选一：
> (a) 证 sqfZp 指数界 + toPolyList 成员/指数对应（偏数学，复用 sqf_correct）；
> (b) 对 `__squarefree_Zp_ir` 输出证「指数 ≤ 输入 natDegree」结构不变量（偏 C++ 归纳，复用 h_loop_lemma）。
> 建议先试 (a)（数学侧引理更可复用，且 h_ih_g2 已把两侧接上）。

### 1.4 组装
```
h_toPolyList_specific:
  unfold toPolyList; simp [Function.comp, hp_1_eq_p]
  refine List.map_congr_left ?_  -- 或对 Array/List 逐元素
  intro (g_h, e) hmem
  have hbound : e.toNat * p < 2^64 := <1.3 的指数界 applied to this element>
  rw [uint64_mul_ofNat_toNat e p (p < 2^64) hbound]
```

---

## Part 2：Branch B（Admit 2，line 1986）—— Yun 循环对应

### 2.1 结构对应
| C++ (`_loop___squarefree_Zp_1_ir_def`, Corpus 182-204) | L2 (`yunLoop`, Algo 62-82) |
|---|---|
| 参数 `(i, w, c, result)` | `(w, c, i, acc, hc)` |
| `if deg(w) > 0 then yun_step else return result` | `if deg(w) = 0 then acc else 递归` |
| `y := gcd(w, c)`；`w/y` 度 > 0 时 append `(w/y, i)` | 同 |
| 递归 `(i+1, y, c/y, result')` | `yunLoop y (c/y) (i+1) acc'` |

`w = normalize(f/c)`，`c = normalize(gcd(f, f'))` 初值。

### 2.2 nl-proof 骨架
1. **L1↔L2 逐步对应**：C++ yun 循环第 i 步的 `(w_i, c_i)`（稀疏表示）经 toPoly 映射到 L2 `yunLoop`
   第 i 步的 `(toPoly w_i, toPoly c_i)`。需：`gcd_l1 ↔ GCDMonoid.gcd`（经 Model.gcd 精化——注意
   Model.gcdAux 仍 partial，需先解决或用 divmodAux 的 WF 版；见风险）、`normalize` 对应、
   `divmod` 对应（Model.divmod 已 WF ✓）。
2. **提取因子对应**：C++ append `(w/y, i)` ↔ L2 `(w/y, i)`，逐元素 toPoly 一致。
3. **循环不变量**：i 计数、acc 累积、终止（w.natDegree 递减）逐步保持 —— 复用 L2 已证的
   `yunLoop_extracts_factor` / `yunLoop_c_natDegree_le` / `yunLoop_crem_dvd_c`。
4. **第二轮 p-th root**（`deg(c_rem) > 0`）：c_rem 导数=0（已证 `derivative_of_yun_remainder_eq_zero`，
   Algo 728-827），递归进入 Branch A 型缩放，指数再乘 p。复用 Part 1 的指数缩放。
5. **safe wrapper + partial def**：C++ 内循环是 partial def；精化只需对 terminating 输入
   （safe wrapper 保证 deg > 0 时进入）的**输出等式**，用 L2 yunLoop 的 WF 结构 + Model divmod/gcd
   的 WF 版做归纳。

### 2.3 关键风险
- **Model.gcdAux 仍 partial def**：Yun 用 gcd。C++ yun 循环调 `polynomial_GCD`。若精化需对 gcd
  归纳，得先把 Model.gcdAux 做成 WF（复用本会话 divmodAux 的度数递减技术 + Euclid natDegree 递减）。
  这是 Branch B 的前置依赖，**可能需要单独一步**。
- **normalize/Associated 摩擦**（CLAUDE.md 已点名）：gcd 返回非 monic，normalize 引入 Associated，
  toPoly 等式需 Associated 桥接。这是 Branch B 最大时间消耗源。

### 2.4 分步计划（Branch B 建议拆成子步，逐个 commit）
1. （前置）Model.gcdAux → WF（若归纳需要）。
2. gcd_l1 ↔ GCDMonoid.gcd 精化 + normalize 对应。
3. yun 循环逐步对应引理（归纳）。
4. 第二轮 p-th root 缩放（复用 Part 1）。
5. 组装 Branch B。

---

## Part 3：5 点审核

1. **数学正确性**：Branch A 指数缩放 e→e·p 对应 f=g^p 重数升 p 倍 ✓；指数界经 contract 度数减半链 ✓。
   Branch B Yun 循环对应经 gcd/divmod 精化 + 已证 L2 不变量 ✓。
2. **无跳步**：Part 1.3 指数界的每步给了引理（natDegree_contract_le / sqf_correct / h_deg_bound）；
   Part 2 每步复用已证 L2 yunLoop 引理。唯一"显然"是 toPolyList 成员对应 —— 需展开（实现注记标出）。
3. **Lean 可行**：UInt64 转换引理已编译验证 ✓；natDegree_contract_le 已存在；sqf_correct 已存在。
   Branch B 的 gcd 精化依赖 Model.gcd WF 化（风险已标）。
4. **工程**：Branch A 无终止性问题（等式证明）。Branch B 的 partial def 由 safe wrapper 规避，
   精化只证输出等式；Model.gcd WF 化是前置。
5. **边界**：常数多项式（natDegree=0）由 safe wrapper 处理 ✓；deg(c_rem)=0 时无第二轮 ✓；
   e=0 不出现（重数 ≥ 1）✓。

---

## 建议实现顺序
1. **先 Branch A**（Part 1）：UInt64 引理（已备）+ 指数界结构引理（路径 a）。较 contained，先落地。
2. **再 Branch B**（Part 2）：先评估 Model.gcd WF 化是否必需，再逐子步。这是全项目最复杂精化之一。

---

## Part 4：实现进度（2026-07-27 更新）

### 已验证落地
- ✅ **`yunLoop_factors_natDegree_pos`**（Algorithm/SquarefreeZp.lean，已编译）：
  `(∀ pr ∈ acc, 0 < pr.1.natDegree) → ∀ pr ∈ (yunLoop w c i acc hc).1, 0 < pr.1.natDegree`。
  强归纳 on measure `w.natDegree + c.natDegree`（复用 yunLoop_c_natDegree_le 的测量递减）。
  是 Branch B 因子非常数的基石。

### 已验证策略（待落地）
- **UInt64 指数转换**（scratch 验证 ✓）：
  ```lean
  lemma uint64_mul_ofNat_toNat (e : UInt64) (b : ℕ) (hb : b < 2^64) (h : e.toNat * b < 2^64) :
      (e * UInt64.ofNat b).toNat = e.toNat * b := by
    rw [UInt64.toNat_mul, UInt64.toNat_ofNat_of_lt' hb, Nat.mod_eq_of_lt h]
  ```
- **指数界 `sqfZp_exponent_le_natDegree`**（简洁证法，避开求和）：
  `pr ∈ sqfZp f ⇒ pr.1^pr.2 ∣ f`（List.dvd_prod + sqf_correct 的 Associated.symm.dvd）
  → `natDegree_le_of_dvd` + `natDegree_pow` → `pr.2 ≤ pr.2 * pr.1.natDegree ≤ f.natDegree`
  （需 `sqfZp_factor_natDegree_pos`：0 < pr.1.natDegree）。
- **contract-lt 事实**（复用 sqfZp decreasing_by，line 573-616）：
  Branch A `(contract p f).natDegree < f.natDegree`；Branch B `(contract p c_rem).natDegree < f.natDegree`。

### 剩余步骤（下轮）
1. **`sqfZp_factor_natDegree_pos`**：`∀ pr ∈ sqfZp f, 0 < pr.1.natDegree`。
   难点：sqfZp 多 `let`（c/w/c_rem/yun_output）结构使 functional induction (`sqfZp.induct`) 的
   case3/case4 binder（~13 个含 let）难匹配，strong induction 又难干净访问中间 let。
   **建议路径**：functional induction + 每 case 用 `rw [sqfZp, dif_neg ‹_›, dif_pos/neg ‹_›]` 展开 +
   `‹_›` 匿名访问 case 假设 + `rename_i` 取 IH；Branch B 两处用 `yunLoop_factors_natDegree_pos _ _ 1 [] _ (by simp)`。
   case3/case4 的内层 `if 0 < c_rem.natDegree` 需 `simp only [‹0 < _›]` 或 split 归约。
2. **`sqfZp_exponent_le_natDegree`**：按上述简洁证法（sqf_correct + factor_pos）。
3. **line 1969 组装**：`unfold toPolyList; simp [Function.comp, hp_1_eq_p]` → `List.map_congr_left` 逐元素 →
   对每个 `(g_h, e) ∈ __squarefree_Zp_ir_safe p g_2`，经 `h_ih_g2 : toPolyList(…) = sqfZp(toPoly g_2)`
   得 `e.toNat ≤ natDegree(toPoly g_2)`，再经 contract 度数链 + h_deg_bound 得 `e.toNat * p < 2^64`，
   套 uint64_mul_ofNat_toNat 收尾。

### Branch B（Part 2）—— 仍待 Model.gcdAux WF 化前置 + 全部实现（下轮大工程）。

---

## Part 5：Branch A 完成（2026-07-27）✅

**line 1969 admit 清零。** 落地引理链：
- `yunLoop_factors_natDegree_pos`、`sqfZp_factor_natDegree_pos`、`sqfZp_exponent_le_natDegree`（Algorithm）
- `listSum_natDegree_lt`（Refinement）：natDegree(toPoly g) < 2^64
- 组装：h_ih_g2 接数组↔sqfZp → 指数界 e≤deg(g_2) → 度数链 deg(g_2)·p=deg(g)<2^64 →
  UInt64.toNat_mul 收尾。

**SquarefreeZp 剩余**：Branch B（Yun 循环，line 2015）+ get_deg_toPoly（line 54）。
**下一步**：Branch B —— 按 Part 2/2.4，先评估 Model.gcdAux WF 化前置。
