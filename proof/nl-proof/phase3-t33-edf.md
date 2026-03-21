# T3.3: EDF 算法建模与正确性

> 状态：nl-proof v2（审核修正版）
> 对应 C++：`polynomial_factorize_zp.hh:294-354` `__edf_Zp`

---

## 0. 数学背景

设 `f ∈ F_p[x]` monic squarefree，所有不可约因子度 = d，`p` 奇素数。

**EDF 三分性（已证，T2.4）**：对任意不可约 q，deg(q) | d，对任意 a ∈ F_p[x]：
```
设 m = (p^d - 1) / 2。
q | a  ∨  q | (a^m - 1)  ∨  q | (a^m + 1)
```

**关键推论**：设 f = q₁ · q₂ · ... · qₖ（两两不同的 monic 不可约，deg = d）。
对任意 a，每个 qⱼ 恰好整除 a、a^m - 1、a^m + 1 之一。
因此 `gcd(a^m - 1, f)` = ∏{qⱼ : qⱼ | a^m - 1}。
若 ∃ qⱼ | (a^m - 1) 且 ∃ qₖ ∤ (a^m - 1)，则 gcd 是 f 的真因子（非平凡分裂）。

---

## 1. 算法模型

### 1.1 EDF 递归分解

**设计决策**：将 edfSplit 逻辑 **inline** 到 edf 主函数中（不拆分为独立函数），
使 `if` 条件的度数 bound 直接出现在 `decreasing_by` 的 context 中，避免跨函数传递证明。

```lean
noncomputable def edf
    (f : Polynomial (ZMod p)) (d : ℕ)
    (splits : List (Polynomial (ZMod p)))
    : List (Polynomial (ZMod p)) :=
  if _hf : f.natDegree ≤ d then [f]
  else
    match splits with
    | [] => [f]
    | a :: rest =>
      let m := (p ^ d - 1) / 2
      let g := normalize (EuclideanDomain.gcd (a ^ m - 1) f)
      if hg : 0 < g.natDegree ∧ g.natDegree < f.natDegree then
        let h := normalize (f /ₘ g)
        edf g d rest ++ edf h d rest
      else
        edf f d rest
termination_by (f.natDegree, splits.length)
```

**终止性**（`decreasing_by` 细节）：

- **`else` 分支**（gcd 不分裂）：`(f.natDegree, rest.length) <_lex (f.natDegree, (a::rest).length)`。
  第一分量相等，第二分量 `rest.length < (a::rest).length`。✓

- **`edf g d rest`**：`(g.natDegree, rest.length) <_lex (f.natDegree, (a::rest).length)`。
  第一分量 `g.natDegree < f.natDegree`，直接从 `hg.2` 获得（`hg` 在 context 中）。✓

- **`edf h d rest`**：`(h.natDegree, rest.length) <_lex (f.natDegree, (a::rest).length)`。
  需证 `h.natDegree < f.natDegree`。推导：
  - `g | f`：`normalize_dvd_iff.mpr (EuclideanDomain.gcd_dvd_right ...)`
  - `g` monic：`Polynomial.monic_normalize hgcd_ne`（gcd ≠ 0 因 context 有 `¬f.natDegree ≤ d` → f 非常数 → f ≠ 0）
  - `f = g * (f /ₘ g)`：`modByMonic_add_div` + mod = 0
  - `natDegree(f) = natDegree(g) + natDegree(f /ₘ g)`：`Polynomial.natDegree_mul`（两因子非零）
  - `natDegree(h) = natDegree(f /ₘ g)`：`natDegree_normalize_eq`
  - `natDegree(h) = natDegree(f) - natDegree(g) < natDegree(f)`（因 `hg.1 : 0 < g.natDegree`）。✓

  此推导约 ~8 行 tactic proof。关键：`hg` 在 `if` 分支 context 中直接可用。

**其他设计决策**：
- 用 `splits : List` 代替 `while(true)` + 随机生成器，彻底去除随机性
- 移除 `hp_odd` 参数：p 的奇偶性只影响 m 的计算，不影响正确性证明；在 L2 中 m = (p^d-1)/2 对任意 p 有定义
- 正确性是结构性的：定理 A、B 对任意 splits 无条件成立

---

## 2. 正确性证明

### 2.1 前置条件

输入 f 满足：
- (P1) f monic
- (P2) Squarefree f
- (P3) ∀ 不可约 q, q | f → q.natDegree = d
- (P4) d ≥ 1
- (P5) p 奇素数
- (P6) 0 < f.natDegree（f 非平凡，对应 C++ `if (get_deg(f) <= 0) return;`）

### 2.2 分裂步骤的性质（inline 在 edf 证明中）

在 edf 的 `if hg : ...` 分支中，g 和 h 满足：
1. `Associated f (g * h)`
2. `0 < g.natDegree < f.natDegree`
3. `0 < h.natDegree < f.natDegree`
4. `Squarefree g` 且 `Squarefree h`（继承自 f）
5. g 和 h 的所有不可约因子度 = d（继承自 f）
6. `Monic g` 且 `Monic h`

**证明**：

**(1) 乘积还原**：`Associated f (g * h)`。

逐步推导：
- `gcd(a^m - 1, f) | f`（gcd 整除第二参数）。
- `normalize(gcd) | f`（`normalize_dvd_iff`）。即 `g | f`。
- g monic（`normalize` 保证，需 `gcd ≠ 0`；gcd ≠ 0 因 f ≠ 0 from P1）。
- `f = g * (f /ₘ g)`（`modByMonic_add_div` + `g | f` → mod = 0）。
- `h = normalize(f /ₘ g)`。`normalize_associated` 给 `Associated h (f /ₘ g)`。
- `Associated.mul_left g (Associated h (f /ₘ g)).symm` 给 `Associated (g * h) (g * (f /ₘ g))`。
- `g * (f /ₘ g) = f`。故 `Associated (g * h) f`，即 `Associated f (g * h)`。✓

**(2) 度数**：`0 < g.natDegree < f.natDegree` 且 `0 < h.natDegree < f.natDegree`。

- 由 `if` 条件直接得 `0 < g.natDegree` 和 `g.natDegree < f.natDegree`。
- `f = g * (f /ₘ g)`，f ≠ 0，g ≠ 0（deg > 0），(f /ₘ g) ≠ 0（`divByMonic_ne_zero_of_ne_zero`）。
- 整环中 `natDegree(f) = natDegree(g) + natDegree(f /ₘ g)`（`Polynomial.natDegree_mul`）。
- `h ~ (f /ₘ g)` → `natDegree(h) = natDegree(f /ₘ g)`（`natDegree_normalize_eq`）。
- 故 `natDegree(h) = natDegree(f) - natDegree(g)`。
- `0 < g.natDegree` → `h.natDegree = f.natDegree - g.natDegree < f.natDegree`。✓
- `g.natDegree < f.natDegree` → `h.natDegree = f.natDegree - g.natDegree > 0`。✓

**(3) Squarefree 继承**：

对 g：`g | f`（步骤 1）且 `Squarefree f`（P2）→ `Squarefree g`（`Squarefree.squarefree_of_dvd`）。✓

对 h：`(f /ₘ g) | f`（因 `f = g * (f /ₘ g)`，dvd_mul_left 给 `(f /ₘ g) | f`）。
`Squarefree f` → `Squarefree (f /ₘ g)`（同上）。
`Associated h (f /ₘ g)` → `Squarefree h`（Squarefree 在 Associated 下保持：
若 d² | h，则 d² | (f /ₘ g)（Associated.dvd），IsUnit d（Squarefree），故 Squarefree h）。✓

**(4) 等度继承**：

对 g：若 q 不可约且 q | g，则 q | f（dvd_trans q | g | f）。由 (P3)：q.natDegree = d。✓
对 h：若 q 不可约且 q | h，则 q | (f /ₘ g)（Associated.dvd），q | f（dvd_trans）。同理 q.natDegree = d。✓

**(5) Monic 继承**：g = normalize(...) monic。h = normalize(...) monic。✓

### 2.3 edf 正确性（主定理）

**定理结构**：无条件证明两个性质，EDFCorrect 作为推论。

---

**定理 A (edf_associated)**：对任意 splits，若 f 满足 (P1)：
```
Associated f (edf f d hp_odd splits).prod
```

**定理 B (edf_preserves)**：对任意 splits，若 f 满足 (P1)-(P6)：
```
∀ q ∈ edf f d hp_odd splits,
  Monic q ∧ Squarefree q ∧ 0 < q.natDegree ∧
  (∀ r, Irreducible r → r ∣ q → r.natDegree = d)
```

**推论 (edf_correct)**：若 f 满足 (P1)-(P6) 且 edf 的所有输出元素 natDegree ≤ d（即 splits 充分），则 `EDFCorrect f d (edf f d hp_odd splits)`。

推论的证明：由定理 A 得乘积还原。由定理 B 得每个元素 monic + squarefree + 所有不可约因子度 = d。由子引理（下面）+ natDegree ≤ d 得每个元素 Irreducible 且 natDegree = d。✓

---

**子引理 (irreducible_of_le_deg)**：
```
f monic, Squarefree f, 0 < f.natDegree, f.natDegree ≤ d,
∀ 不可约 q, q | f → q.natDegree = d
→ Irreducible f ∧ f.natDegree = d
```

**证明**：

**Step 1**：`f.natDegree = d`。
- f 非 unit：`0 < f.natDegree` → `natDegree_eq_zero_of_isUnit` 的逆否 → `¬IsUnit f`。
- f ≠ 0：monic → 非零（`Monic.ne_zero`）。
- `WfDvdMonoid.exists_irreducible_factor (¬IsUnit f) (f ≠ 0)`：∃ 不可约 q₁, q₁ | f。
- q₁.natDegree = d（前置条件 P3）。
- `Polynomial.natDegree_le_of_dvd (q₁ | f) (f ≠ 0)`：q₁.natDegree ≤ f.natDegree。
- 故 d ≤ f.natDegree。结合 f.natDegree ≤ d 得 f.natDegree = d。✓

**Step 2**：`Irreducible f`。

使用 `irreducible_iff`：`Irreducible f ↔ ¬IsUnit f ∧ ∀ a b, f = a * b → IsUnit a ∨ IsUnit b`。
`¬IsUnit f`：Step 1 已证。✓

需证：`∀ a b, f = a * b → IsUnit a ∨ IsUnit b`。
设 f = a * b。反证：假设 `¬IsUnit a ∧ ¬IsUnit b`。

- a ≠ 0, b ≠ 0：f = a * b ≠ 0（f monic → f ≠ 0），整环中 a ≠ 0 且 b ≠ 0。
- ¬IsUnit a, a ≠ 0 → natDegree(a) ≥ 1：
  `natDegree = 0` 时多项式 = C(c)，c ≠ 0 → IsUnit（域上非零常数是 unit）。逆否：¬IsUnit → natDegree ≥ 1。
  （Lean 路径：`Polynomial.eq_C_of_natDegree_eq_zero` + `Polynomial.isUnit_C` + `IsUnit.mk0`）
- 同理 natDegree(b) ≥ 1。
- a 有不可约因子 q₁：`WfDvdMonoid.exists_irreducible_factor (¬IsUnit a) (a ≠ 0)`。
- q₁ | a | f → q₁ | f → q₁.natDegree = d（P3）。
- `Polynomial.natDegree_le_of_dvd (q₁ | a) (a ≠ 0)`：natDegree(a) ≥ d。
- 同理 natDegree(b) ≥ d。
- `Polynomial.natDegree_mul (a ≠ 0) (b ≠ 0)`：natDegree(f) = natDegree(a) + natDegree(b)。
- d = natDegree(f) = natDegree(a) + natDegree(b) ≥ d + d = 2d。
- d ≥ 1（P4）→ 2d > d → d ≥ 2d 矛盾。✓

故 ¬(¬IsUnit a ∧ ¬IsUnit b)，即 IsUnit a ∨ IsUnit b。✓
（Lean 路径：`by_contra h; push_neg at h; obtain ⟨ha, hb⟩ := h; ...`）

---

**定理 A 的证明**（对 `(f.natDegree, splits.length)` 的 well-founded 归纳）：

**Base case**（`f.natDegree ≤ d`）：
result = `[f]`。`[f].prod = f`。`Associated f f`。✓

**Case `splits = []`**：
result = `[f]`。同上。✓

**Case `none`**（edfSplit 返回 none）：
result = `edf f d hp_odd rest`。
归纳假设（splits.length 减 1，同一 lex 序）：`Associated f result.prod`。✓

**Case `some (g, h)`**（edfSplit 返回分裂）：
result = `edf g d hp_odd rest ++ edf h d hp_odd rest`。

edfSplit_correct (1)：`Associated f (g * h)`。
g monic（normalize）→ 归纳假设适用于 g（natDegree 严格减小）：
  `Associated g (edf g ...).prod`。
h monic（normalize）→ 归纳假设适用于 h（natDegree 严格减小）：
  `Associated h (edf h ...).prod`。
`List.prod_append`：`(edf_g ++ edf_h).prod = (edf_g).prod * (edf_h).prod`。
链：`f ~ g * h ~ (edf_g).prod * (edf_h).prod = result.prod`。
（用 `Associated.mul` 组合两个 Associated）。✓

---

**定理 B 的证明**（同一归纳结构）：

**Base case**（`f.natDegree ≤ d`）：
result = `[f]`。唯一元素是 f。
- Monic f：(P1)。✓
- Squarefree f：(P2)。✓
- 0 < f.natDegree：(P6)。✓
- ∀ r irred, r | f → r.natDegree = d：(P3)。✓

**Case `splits = []`**：
result = `[f]`。同上。✓

**Case `none`**：
result = edf f d hp_odd rest。归纳假设直接适用（f 不变，splits 缩短）。✓

**Case `some (g, h)`**：
result = edf g ... ++ edf h ...。

edfSplit_correct 给出 g, h 满足：
- Monic：(5)。✓
- Squarefree：(3)。✓
- 0 < natDegree：(2)。✓
- 等度继承：(4)。✓

即 g, h 满足 (P1)-(P6)（P4, P5 不变）。

归纳假设适用于 g 和 h（natDegree 严格减小）。
∀ q ∈ result：q ∈ edf g ... 或 q ∈ edf h ...。
两边由归纳假设得满足所有性质。✓

---

## 3. 辅助引理清单

| 引理 | 内容 | 依赖 | 状态 |
|------|------|------|------|
| `edf_trichotomy` | q \| a ∨ q \| (a^m-1) ∨ q \| (a^m+1) | T2.4 | ✅ 已证（不用于正确性证明） |
| `irreducible_of_le_deg` | f monic sqfree, 0 < deg(f) ≤ d, ∀ irred factor deg=d → Irreducible f ∧ deg=d | 度数论证 | 需新证 ~15 行 |
| 分裂性质 (1)-(6) | inline 在 edf 证明中 | normalize, gcd, modByMonic | inline ~20 行 |
| `Squarefree.squarefree_of_dvd` | g \| f, Squarefree f → Squarefree g | Mathlib | ✅ |
| `List.prod_append` | (l₁ ++ l₂).prod = l₁.prod * l₂.prod | Mathlib | ✅ |
| `Associated.mul` | Associated a b → Associated c d → Associated (a*c) (b*d) | Mathlib | ✅ |

---

## 4. 形式化策略

### Step 1: 函数定义 + 终止性（~40 行）
- `edf` 定义（edfSplit 逻辑 inline）
- `termination_by (f.natDegree, splits.length)`
- `decreasing_by`：3 个分支各需证明 lex 递减
  - `else`（不分裂）：`Nat.lt_of_lt_of_le` on list length — ~1 行
  - `edf g`：`hg.2` 直接给 `g.natDegree < f.natDegree` — ~1 行
  - `edf h`：度数算术 ~8 行（modByMonic_add_div + natDegree_mul + natDegree_normalize_eq + omega）
- **预期 sorry: 0**

### Step 2: irreducible_of_le_deg（~20 行）
- Step 1（natDegree = d）：exists_irreducible_factor + natDegree_le_of_dvd
- Step 2（Irreducible）：正向构造 `irreducible_iff.mpr`，内部反证 `by_contra` + 度数矛盾
- **预期 sorry: 0**

### Step 3: edf_correct 主定理（~60 行）
- `edf_associated` + `edf_preserves` 合并为单定理（同一归纳）
- 分裂性质 (1)-(6) inline 在 `some` 分支中证明（~20 行）
- Associated 链：`Associated.mul` + `List.prod_append`
- `edf_correct` 推论：+ `irreducible_of_le_deg` 得 EDFCorrect
- **预期 sorry: 0**

### 总计：~120 行，预期 0 sorry

**注**：整个 EDF 证明**不使用 `edf_trichotomy`**。三分性保证存在有效随机元素（概率 ≥ 1/2），但正确性证明是结构性的（"若分裂发生则正确"），不依赖三分性。三分性仅用于概率分析（不形式化）。

---

## 5. 与 C++ 的对应

| C++ (L294-354) | Lean 模型 | 说明 |
|----------------|-----------|------|
| `deg(f) == d → return {f}` | base case `f.natDegree ≤ d` | Lean 用 ≤ 更通用 |
| `while(true) { r = random; ... }` | `splits : List` 参数 | 去随机化 |
| `g_pow = r^{(p^d-1)/2} mod f` | `a ^ m`（纯数学） | L2 不建模 powmod |
| `g = gcd(g_pow - 1, f)` | `EuclideanDomain.gcd (a^m - 1) f` | inline 在 edf 中 |
| `deg(g) > 0 ∧ deg(g) < deg(f)` | `if hg : ...`（携带证明） | hg 在 decreasing_by 可用 |
| 递归 `__edf_Zp(g)` + `__edf_Zp(f/g)` | `edf g ... ++ edf h ...` | 直接对应 |
