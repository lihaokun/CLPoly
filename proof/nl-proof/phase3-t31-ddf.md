# T3.1: DDF 算法建模与正确性

> 状态：待形式化（nl-proof v2，审核修正版）
> 对应 C++：`polynomial_factorize_zp.hh:247-292` `__ddf_Zp`

---

## 1. 算法模型

### C++ 算法（简化）

```
输入: 首一无平方 f ∈ F_p[x]
输出: [(g_d, d)] 其中 g_d = ∏{f 的 d 次不可约因子}

h ← x,  f* ← f
for d = 1, 2, ...:
    if deg(f*) < 2d: break
    h ← h^p mod f*
    g_d ← monic_gcd(h - x, f*)
    if deg(g_d) > 0:
        output (g_d, d)
        f* ← f* /ₘ g_d
        h ← h mod f*
if deg(f*) > 0: output (f*, deg(f*))
```

### Lean 4 模型

```lean
noncomputable def ddfLoop
    (h f_star : Polynomial (ZMod p)) (d : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    : List (Polynomial (ZMod p) × ℕ) :=
  if f_star.natDegree < 2 * d then
    if 0 < f_star.natDegree then acc ++ [(f_star, f_star.natDegree)]
    else acc
  else
    let h' := (h ^ p) %ₘ f_star
    let gd := normalize (EuclideanDomain.gcd (h' - X) f_star)
    if 0 < gd.natDegree then
      let f_new := f_star /ₘ gd
      ddfLoop (h' %ₘ f_new) f_new (d + 1) (acc ++ [(gd, d)])
    else
      ddfLoop h' f_star (d + 1) acc
termination_by f_star.natDegree + 1 - 2 * d

noncomputable def ddf (f : Polynomial (ZMod p)) :
    List (Polynomial (ZMod p) × ℕ) :=
  ddfLoop X f 1 []
```

**建模说明**：
- `%ₘ` = `modByMonic`，`/ₘ` = `divByMonic`（f_star、gd 始终 monic，保证操作正确）
- `normalize` 对应 C++ `__upoly_make_monic`（`EuclideanDomain.gcd` 不保证 monic）
- 终止度量：`f_star.natDegree + 1 - 2 * d`（Nat 截断减法）

---

## 2. 终止性

**度量**：`μ = f_star.natDegree + 1 - 2 * d`

只在 `¬(f_star.natDegree < 2 * d)` 即 `f_star.natDegree ≥ 2 * d` 时递归，此时 `μ ≥ 1`。

### 2.1 无分裂情况（gd 度 = 0）

f_star 不变，d → d+1。
新 μ = `f_star.natDegree + 1 - 2*(d+1)` = μ - 2。

在 ℕ 截断算术中：μ ≥ 1，若 μ ≥ 2 则 μ - 2 < μ；若 μ = 1 则新 μ = 0 < 1 = μ。✓

**decreasing_by 路径**：`omega`（给定 `f_star.natDegree ≥ 2 * d`）。

### 2.2 分裂情况（gd 度 > 0）

f_new = f_star /ₘ gd。需要在函数体内（无不变量假设）建立度数关系。

**`decreasing_by` 推导链**：

1. `gd_raw := EuclideanDomain.gcd (h' - X) f_star`
2. `EuclideanDomain.gcd_dvd_right (h' - X) f_star` → `gd_raw ∣ f_star`
3. `gd = normalize gd_raw` → `Associated gd gd_raw` → `gd ∣ f_star`
   （`dvd_of_associated` + 传递性）
4. `0 < gd.natDegree` → `gd ≠ 0` → `Monic gd`（`normalize` 对非零多项式产生 monic）
5. `Monic gd ∧ gd ∣ f_star` → `f_star %ₘ gd = 0`（`dvd_iff_modByMonic_eq_zero`）
   → `f_star = gd * (f_star /ₘ gd)`（由 `modByMonic_add_div`）
6. `natDegree (gd * f_new) = natDegree gd + natDegree f_new`（`natDegree_mul`，域上无零因子）
   → `f_new.natDegree = f_star.natDegree - gd.natDegree`
7. `gd.natDegree ≥ 1` → `f_new.natDegree ≤ f_star.natDegree - 1`
8. 故新 μ = `f_new.natDegree + 1 - 2*(d+1)` ≤ `(f_star.natDegree - 1) + 1 - 2*d - 2`
   = `f_star.natDegree - 2 - 2*d` < `f_star.natDegree + 1 - 2*d` = μ（ℕ 中，μ ≥ 1 保证）

**decreasing_by 路径**：步骤 1-7 提供 `f_new.natDegree ≤ f_star.natDegree - 1`，然后 `omega`。

**注意**：步骤 6 中 `natDegree_mul` 要求域上无零因子。`Polynomial (ZMod p)` 是域上多项式环，`ZMod p` 是整环（p 素数 → 域），故 `Polynomial (ZMod p)` 无零因子。✓

---

## 3. 循环不变量

在 `ddfLoop h f_star d acc` 的每次调用入口处，假设以下不变量成立：

| 标号 | 不变量 | 含义 |
|------|--------|------|
| P0 | `d ≥ 1` | d 始终 ≥ 1（初始 d=1，递增 d+1 ≥ 2） |
| P1 | `f = f_star * (acc.map Prod.fst).prod` | 乘积分解（精确等式，均 monic） |
| P2 | `f_star ∣ (h - X ^ (p ^ (d-1)))` | h 是 X^{p^{d-1}} 模 f* 的代表元 |
| P3 | `Monic f_star` | f* 始终首一 |
| P4 | `Squarefree f_star` | f* 始终无平方 |
| P5 | `∀ q, Irreducible q → q ∣ f_star → q.natDegree ≥ d` | f* 的因子度 ≥ d |
| P6 | 对每个 `(g, d') ∈ acc`：`g ∣ f`、`Monic g`、`∀ q irred, q ∣ g → deg(q) = d'` | acc 中条目正确 |

### 初始调用满足不变量

`ddfLoop X f 1 []`：
- P0: d = 1 ≥ 1 ✓
- P1: `f = f * [].prod = f * 1 = f` ✓
- P2: `f ∣ (X - X^{p^0}) = f ∣ (X - X) = f ∣ 0` ✓（任何多项式整除 0）
- P3: 给定 `Monic f` ✓
- P4: 给定 `Squarefree f` ✓
- P5: 不可约多项式度 ≥ 1 = d（`Irreducible.natDegree_pos`） ✓
- P6: acc = []，空真 ✓

---

## 4. 分裂步正确性

以下在 else 分支（`f_star.natDegree ≥ 2 * d`）内论证。

### 4.1 h' ≡ X^{p^d} (mod f_star)

由 P2: `f_star ∣ (h - X^{p^{d-1}})`。

**Step A**: `f_star ∣ (h^p - X^{p^d})`

选择 Frobenius 路径（char p 下最直接）：
1. `f_star ∣ (h - X^{p^{d-1}})` （P2）
2. `dvd_pow_self`：若 `a ∣ b` 则 `a ∣ b^n`。故 `f_star ∣ (h - X^{p^{d-1}})^p`
3. Frobenius：在 char p 环中 `(a - b)^p = a^p - b^p`
   Mathlib: `sub_pow_char (ZMod p) h (X^{p^{d-1}})` 给出
   `(h - X^{p^{d-1}})^p = h^p - (X^{p^{d-1}})^p = h^p - X^{p^d}`
   （最后一步用 `pow_mul_comm` 或 `pow_pow`：`(X^{p^{d-1}})^p = X^{p^{d-1}·p} = X^{p^d}`）
4. 由 2 和 3：`f_star ∣ (h^p - X^{p^d})` ✓

**Mathlib 路径**：
- `sub_pow_char R a b : (a - b) ^ ringChar R = a ^ ringChar R - b ^ ringChar R`（或 `sub_pow_char_of_commute` + `Commute.all`）
- `dvd_pow (h : a ∣ b) (n : ℕ) (hn : n ≠ 0) : a ∣ b ^ n`
  - **注意**：`dvd_pow` 签名是 `a ∣ b → n ≠ 0 → a ∣ b ^ n`。不是 `dvd_pow_self`。
  - 实际上，我们需要的是：`a ∣ b → a ∣ b ^ p`。用 `dvd_pow h (Nat.Prime.ne_zero hp)`。

**备选路径**（不依赖 char p）：
- `sub_dvd_pow_sub_pow`（或 `Commute.sub_dvd_pow_sub_pow`）：`(a-b) ∣ (a^n - b^n)` 在交换环中成立
- 若 Mathlib 有此引理，则无需 Frobenius，一步完成：
  `(h - X^{p^{d-1}}) ∣ (h^p - X^{p^d})` + P2 传递性

**Step B**: `f_star ∣ (h^p - h')`

`h' = h^p %ₘ f_star`。由 `modByMonic_add_div`（P3 保证 Monic f_star）：
`h^p %ₘ f_star + f_star * (h^p /ₘ f_star) = h^p`
即 `h^p - h' = f_star * (h^p /ₘ f_star)`，故 `f_star ∣ (h^p - h')`。 ✓

**Step C**: 合并

`h' - X^{p^d} = (h' - h^p) + (h^p - X^{p^d})`

f_star 整除两个加数（Step B 给 `f_star ∣ (h^p - h')` 即 `f_star ∣ -(h' - h^p)` 即 `f_star ∣ (h' - h^p)`，Step A 给 `f_star ∣ (h^p - X^{p^d})`），故 `f_star ∣ (h' - X^{p^d})`。 ✓

### 4.2 gd 的不可约因子刻画

**引理（核心桥接）**：设 P0-P5 成立。对任意不可约 q：

```
q ∣ gd  ↔  q ∣ f_star ∧ q.natDegree = d
```

**证明**：

设 `gd_raw = EuclideanDomain.gcd (h' - X) f_star`，`gd = normalize gd_raw`。
`Associated gd gd_raw`（`normalize_associated`）。

**Claim**：对不可约 q，`q ∣ gd ↔ q ∣ gd_raw`。
证明：`Associated a b` → (`q ∣ a ↔ q ∣ b`)。
详细：`a = b * u`（u 单位）。若 `q ∣ a = b * u`，q 不可约故 prime（UFD），
`q ∣ b` 或 `q ∣ u`。`q` 不可约不是单位，不能整除单位 `u`，故 `q ∣ b`。反向类似。 ✓

**Step 1**: `q ∣ gcd_raw ↔ q ∣ (h'-X) ∧ q ∣ f_star`

正向（两个传递性）：
- `gcd_dvd_left : gcd_raw ∣ (h' - X)` → `q ∣ gcd_raw ∣ (h'-X)` → `q ∣ (h'-X)`
- `gcd_dvd_right : gcd_raw ∣ f_star` → `q ∣ gcd_raw ∣ f_star` → `q ∣ f_star`

反向（GCD 的通用性质）：
- `q ∣ (h'-X)` 且 `q ∣ f_star` → `q ∣ gcd_raw`（`EuclideanDomain.dvd_gcd`）✓

**Step 2**: `q ∣ (h'-X) ∧ q ∣ f_star  ↔  q ∣ (X^{p^d}-X) ∧ q ∣ f_star`

由 §4.1，`f_star ∣ (h' - X^{p^d})`。

正向：给定 `q ∣ (h'-X)` 且 `q ∣ f_star`。
1. `q ∣ f_star` 且 `f_star ∣ (h' - X^{p^d})` → `q ∣ (h' - X^{p^d})`（传递性）
2. `(h'-X) - (h'-X^{p^d}) = X^{p^d} - X`
3. `q ∣ (h'-X)` 且 `q ∣ (h'-X^{p^d})` → `q ∣ ((h'-X) - (h'-X^{p^d}))` → `q ∣ (X^{p^d} - X)` ✓

反向：给定 `q ∣ (X^{p^d}-X)` 且 `q ∣ f_star`。
1. `q ∣ f_star` 且 `f_star ∣ (h' - X^{p^d})` → `q ∣ (h' - X^{p^d})`（传递性）
2. `(h'-X^{p^d}) + (X^{p^d}-X) = h' - X`
3. `q ∣ (h'-X^{p^d})` 且 `q ∣ (X^{p^d}-X)` → `q ∣ ((h'-X^{p^d}) + (X^{p^d}-X))` → `q ∣ (h'-X)` ✓

**Step 3**: `q ∣ (X^{p^d}-X)  ↔  q.natDegree ∣ d`

直接由 T2.2 + T2.2'：
- T2.2 `irreducible_dvd_X_pow_sub_X`：Irreducible q ∧ q.natDegree ∣ d → q ∣ X^{p^d}-X
- T2.2' `dvd_X_pow_sub_X_imp_natDegree_dvd`：Irreducible q ∧ q ∣ X^{p^d}-X → q.natDegree ∣ d

注：T2.2 不要求 Monic q。 ✓

**Step 4**: `q ∣ f_star ∧ q.natDegree ∣ d  ↔  q ∣ f_star ∧ q.natDegree = d`

正向：
- `q ∣ f_star` → `q.natDegree ≥ d`（P5）
- `q.natDegree ∣ d`：存在 k 使 `d = q.natDegree * k`。
  由 P0（`d ≥ 1`），`k ≥ 1`（若 k = 0 则 d = 0 矛盾）。
  故 `d = q.natDegree * k ≥ q.natDegree * 1 = q.natDegree`，即 `q.natDegree ≤ d`。
- 合并 `≥` 和 `≤`：`q.natDegree = d` ✓

反向：`q.natDegree = d → q.natDegree ∣ d`（`dvd_refl`）✓

**综合 Step 1-4**：`q ∣ gd  ↔  q ∣ f_star ∧ q.natDegree = d` ✓

### 4.3 不变量保持（分裂情况）

设 `f_new = f_star /ₘ gd`，`h_new = h' %ₘ f_new`，`acc_new = acc ++ [(gd, d)]`。

**前提事实**（在分裂分支内可推导）：
- F1: `gd ∣ f_star`（§2.2 推导链步骤 2-3）
- F2: `Monic gd`（§2.2 推导链步骤 4）
- F3: `f_star = gd * f_new`（F1 + F2 + `modByMonic_add_div` + `dvd_iff_modByMonic_eq_zero`）
- F4: `Monic f_new`（F3 + `leadingCoeff_mul` + P3：`1 = LC(f_star) = LC(gd) * LC(f_new) = 1 * LC(f_new)`）

**P0**: d+1 ≥ 1 + 1 = 2 ≥ 1 ✓

**P1**: `f = f_new * (acc_new.map Prod.fst).prod`
- `(acc_new.map Prod.fst).prod = (acc.map Prod.fst).prod * gd`（`List.prod_append`）
- `f_new * acc_prod_new = f_new * (acc_prod * gd) = (f_new * gd) * acc_prod`（结合律+交换律）
- `f_new * gd = f_star`（由 F3：`f_star = gd * f_new`，交换律）
- `= f_star * acc_prod = f`（P1 旧值）✓

**P2**: `f_new ∣ (h_new - X^{p^d})`
1. §4.1 给出 `f_star ∣ (h' - X^{p^d})`
2. F3: `f_star = gd * f_new` → `f_new ∣ f_star` → `f_new ∣ (h' - X^{p^d})`（传递性）
3. `h_new = h' %ₘ f_new`，由 `modByMonic_add_div`（F4 保证 Monic f_new）：
   `f_new ∣ (h' - h_new)`
4. `h_new - X^{p^d} = (h_new - h') + (h' - X^{p^d})`
   f_new 整除两个加数（3 给 `f_new ∣ (h' - h_new)` 即 `f_new ∣ -(h_new - h')`，2 给第二项），
   故 `f_new ∣ (h_new - X^{p^d})` ✓

**P3**: `Monic f_new` — 由 F4 ✓

**P4**: `Squarefree f_new`
- F3: `f_new ∣ f_star`（f_star = gd * f_new → f_new | f_star 通过 dvd_mul_left）
- `Squarefree.squarefree_of_dvd (dvd_mul_left f_new gd) P4` ✓

**P5**: `∀ q irred, q ∣ f_new → q.natDegree ≥ d+1`

假设不可约 q 满足 `q ∣ f_new`。

Step a: `q ∣ f_new` → `q ∣ f_star`（F3 传递性）→ `q.natDegree ≥ d`（P5 旧值）

Step b: 排除 `q.natDegree = d`。假设 `q.natDegree = d`，推矛盾：
- `q ∣ f_star ∧ q.natDegree = d` → `q ∣ gd`（§4.2 反方向）
- 但 `q ∣ f_new` 且 `q ∣ gd`
- **需要：`IsCoprime f_new gd`**

  **`squarefree_coprime_of_mul` 的完整论证**：

  假设存在不可约 r 满足 `r ∣ f_new` 且 `r ∣ gd`。
  则 `r ∣ f_new` → `f_new = r * a`，`r ∣ gd` → `gd = r * b`。
  `f_star = gd * f_new = (r * b) * (r * a) = r * r * (b * a)` → `r * r ∣ f_star`。
  即 `r ^ 2 ∣ f_star`（用 `sq` = `r * r`）。
  但 `Squarefree f_star`：`∀ c, c * c ∣ f_star → IsUnit c`（定义）。
  故 `IsUnit r`。但 r 不可约 → `¬ IsUnit r`（`Irreducible.not_unit`）。矛盾。

  故不存在同时整除 f_new 和 gd 的不可约元。
  在 UFD 中，"无公共不可约因子" 等价于 `IsCoprime`。 ✓

- `IsCoprime f_new gd` → `q ∣ f_new ∧ q ∣ gd → IsUnit q`（`IsCoprime.isUnit_of_dvd`）
- `IsUnit q` 与 `Irreducible q` 矛盾

Step c: 由 a 和 b，`q.natDegree ≥ d` 且 `q.natDegree ≠ d`，故 `q.natDegree ≥ d + 1` ✓

**P6**: acc_new 正确
- 旧 acc 条目不变，满足 P6 ✓
- 新条目 (gd, d)：
  - `gd ∣ f`：`gd ∣ f_star`（F1）且 `f_star ∣ f`（P1：`f = f_star * acc_prod` → `dvd_mul_right`）→ 传递性 ✓
  - `Monic gd`：F2 ✓
  - `∀ q irred, q ∣ gd → q.natDegree = d`：§4.2 正方向（`q ∣ gd → q ∣ f_star ∧ deg(q) = d`）的第二分量 ✓

### 4.4 不变量保持（无分裂情况）

`gd.natDegree = 0`：gd_raw = gcd(h'-X, f_star) ≠ 0（因 gcd_dvd_right 且 f_star ≠ 0（P3 Monic → 非零）），
故 gd = normalize gd_raw ≠ 0。natDegree = 0 + 非零 → gd 是非零常数 → gcd(h'-X, f_star) ~ 常数 ~ 1。

递归调用 `ddfLoop h' f_star (d+1) acc`。

**P0**: d+1 ≥ 2 ≥ 1 ✓

**P1-P4**: f_star 和 acc 不变 ✓

**P5**: `∀ q irred, q ∣ f_star → q.natDegree ≥ d+1`

由 §4.2：`q ∣ gd ↔ q ∣ f_star ∧ q.natDegree = d`。
gd.natDegree = 0 → 不存在度 ≥ 1 的不可约 q 满足 `q ∣ gd`
（若存在则 q.natDegree ≥ 1 → gd.natDegree ≥ q.natDegree ≥ 1，矛盾）。

但 §4.2 的反方向说 `q ∣ f_star ∧ q.natDegree = d → q ∣ gd`。
对位否定：`¬(q ∣ gd) → ¬(q ∣ f_star ∧ q.natDegree = d)`。
即：对任何不可约 q，若 `q ∣ f_star`，则 `q.natDegree ≠ d`。

旧 P5 给 `q.natDegree ≥ d`。排除 `= d`，得 `q.natDegree ≥ d + 1` ✓

**P6**: 不变 ✓

**P2**: `f_star ∣ (h' - X^{p^d})`（已在 §4.1 证明，此处 d 即为新 d+1 的 (d+1)-1 = d）✓

---

## 5. 终止情况正确性

### 5.1 deg(f_star) ≥ 1 且 deg(f_star) < 2d（输出 f_star）

result = acc ++ [(f_star, deg(f_star))]

**f_star 是不可约的**：

1. P3: `Monic f_star`（特别地 f_star ≠ 0）
2. P4: `Squarefree f_star`
3. P5: 所有不可约因子 q 满足 `q.natDegree ≥ d`
4. 在 UFD 中，f_star 有不可约分解 `f_star ~ q₁ * q₂ * ... * qₖ`（各 qᵢ 不可约）
5. Squarefree → 各 qᵢ 两两不 associated（否则 qᵢ² ∣ f_star）
6. 若 k ≥ 2：`deg(f_star) = deg(q₁) + ... + deg(qₖ) ≥ k * d ≥ 2d`，与 `deg(f_star) < 2d` 矛盾
7. `deg(f_star) ≥ 1` → f_star 不是单位 → k ≥ 1
8. 由 6 和 7：k = 1，即 `f_star ~ q₁`（Associated）
9. `Associated` 保持 `Irreducible`（`Irreducible.associated`）→ `Irreducible f_star` ✓

**验证 DDFCorrect 五条**：

**条件 1**：`∀ pr ∈ result, pr.1 ∣ f`
- acc 条目：P6 给出 `g ∣ f` ✓
- f_star：P1 → `f = f_star * acc_prod` → `dvd_mul_right f_star acc_prod` ✓

**条件 2**：`∀ pr ∈ result, ∀ q irred, q ∣ pr.1 → q.natDegree = pr.2`
- acc 条目：P6 ✓
- (f_star, deg(f_star))：
  f_star 不可约（上述），q 不可约且 `q ∣ f_star`。
  UFD 中两个不可约元，一个整除另一个 → 它们 Associated：
  `f_star = q * r`，f_star 不可约 → `IsUnit r` 或 `IsUnit q`。
  q 不可约 → `¬ IsUnit q`，故 `IsUnit r`，即 `Associated f_star q`。
  `Associated` → `natDegree` 相等 → `q.natDegree = f_star.natDegree` ✓

**条件 3**：`∀ q irred, q ∣ f → ∃ pr ∈ result, q ∣ pr.1`
- P1: `f = f_star * acc_prod`
- q 不可约 → q prime（`Irreducible.prime`，UFD 中）
- `q ∣ f = f_star * acc_prod` → `q ∣ f_star` 或 `q ∣ acc_prod`（`Prime.dvd_mul`）
- 情况 1：`q ∣ f_star` → (f_star, deg(f_star)) ∈ result，且 `q ∣ f_star` ✓
- 情况 2：`q ∣ acc_prod = g₁ * g₂ * ... * gₖ`（其中 gᵢ = acc[i].fst）
  对列表归纳使用 `Prime.dvd_mul`：
  `q ∣ g₁ * (g₂ * ... * gₖ)` → `q ∣ g₁` 或 `q ∣ g₂ * ... * gₖ`
  最终得 `∃ i, q ∣ gᵢ`，对应某个 acc 条目 ✓

**条件 4**：`Associated f (result.map Prod.fst).prod`
- `result.map Prod.fst = (acc.map Prod.fst) ++ [f_star]`
- `List.prod_append`：乘积 = `acc_prod * f_star`
- 交换律：`= f_star * acc_prod = f`（P1 精确等式）
- `f = result_prod` → `Associated f result_prod`（`Associated.refl`）✓

**条件 5**：`∀ pr ∈ result, Monic pr.1`
- acc 条目：P6 ✓
- f_star：P3 ✓

### 5.2 deg(f_star) = 0（不输出 f_star）

result = acc

- P3 Monic + natDegree = 0 → `f_star = 1`
  （Monic 意味 `leadingCoeff = 1`，natDegree = 0 意味多项式是常数，故 `f_star = C 1 = 1`）
- P1: `f = 1 * acc_prod = acc_prod`

**条件 1**：P6 给出每个 `g ∣ f` ✓

**条件 2**：P6 给出不可约因子度数正确 ✓

**条件 3**：`q irred ∧ q ∣ f = q ∣ acc_prod`。
q prime → 归纳 `Prime.dvd_mul` 得 `∃ gᵢ ∈ acc, q ∣ gᵢ` ✓

**条件 4**：`Associated f acc_prod` 由 `f = acc_prod`（精确等式）✓

**条件 5**：P6 ✓

---

## 6. 顶层定理

```lean
theorem ddf_correct
    (f : Polynomial (ZMod p)) (hm : Monic f) (hsq : Squarefree f) :
    DDFCorrect f (ddf f)
```

由 `ddf f = ddfLoop X f 1 []`，初始不变量 P0-P6 满足（§3），由 §4-§5 递归论证得 DDFCorrect。

**实际证明结构**：对 `ddfLoop` 使用 `ddfLoop.induct` 归纳原理（Lean 4 自动生成），
在每个分支验证不变量保持（§4）或终止条件正确性（§5）。

```lean
-- 草案签名：中间引理（归纳携带不变量）
theorem ddfLoop_correct
    (f h f_star : Polynomial (ZMod p)) (d : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    (hd : d ≥ 1)                                          -- P0
    (hprod : f = f_star * (acc.map Prod.fst).prod)         -- P1
    (hcong : f_star ∣ (h - X ^ (p ^ (d - 1))))            -- P2
    (hmonic : Monic f_star)                                -- P3
    (hsqf : Squarefree f_star)                             -- P4
    (hfac : ∀ q, Irreducible q → q ∣ f_star → q.natDegree ≥ d)  -- P5
    (hacc : ∀ pr ∈ acc, pr.1 ∣ f ∧ Monic pr.1 ∧
            (∀ q, Irreducible q → q ∣ pr.1 → q.natDegree = pr.2)) -- P6
    : DDFCorrect f (ddfLoop h f_star d acc) := by
  -- 对 ddfLoop 使用 functional induction
  ...
```

---

## 7. 辅助引理清单

| 引理 | 内容 | 预计难度 | Mathlib 状态 |
|------|------|---------|-------------|
| `sub_pow_char_apply` | char p 下 `(a-b)^p = a^p - b^p` | 低 | `sub_pow_char` 已有 |
| `dvd_pow_of_dvd` | `a ∣ b → a ∣ b ^ n`（n ≥ 1） | 低 | `dvd_pow` 已有 |
| `irred_dvd_associated` | `Irreducible q → (q ∣ normalize a ↔ q ∣ a)` | 低 | 从 `Associated` 推 |
| `monic_normalize_of_ne_zero` | `a ≠ 0 → Monic (normalize a)` | 低 | `normalize_monic`? |
| `divByMonic_mul_cancel` | `Monic g → g ∣ f → g * (f /ₘ g) = f` | 低 | 从 `modByMonic_add_div` + `dvd_iff_modByMonic_eq_zero` 推 |
| `monic_divByMonic_of_dvd` | `Monic f → Monic g → g ∣ f → Monic (f /ₘ g)` | 低 | `leadingCoeff_mul` 推 |
| `natDegree_divByMonic_of_dvd` | `Monic g → g ∣ f → g ≠ 0 → deg(f /ₘ g) = deg(f) - deg(g)` | 低 | `natDegree_mul` + F3 |
| `squarefree_mul_imp_coprime` | `Squarefree (a * b) → IsCoprime a b` | 中 | 需自证（q² 矛盾） |
| `prime_dvd_list_prod` | `Prime q → q ∣ l.prod → ∃ a ∈ l, q ∣ a` | 低 | `Prime.dvd_list_prod`? 或归纳 `Prime.dvd_mul` |

---

## 8. 预计行数与风险

- 函数定义 + 终止性（含 decreasing_by）：~40 行
- 辅助引理：~80 行
- 主定理（归纳证明）：~190 行
- **总计：~310 行**

**主要风险**：
1. `decreasing_by` 分裂情况：需在函数体内推导 `gd ∣ f_star` + `Monic gd` + `natDegree` 计算，约 10-15 行
2. `ddfLoop.induct` 使用：Lean 4 自动生成的归纳原理可能需要适配（参数顺序、分支结构）
3. `squarefree_mul_imp_coprime`：约 10-15 行，核心是构造 `q² ∣ a*b` 的 dvd 证据
4. 列表乘积操作：`List.prod_append`、`List.map_append` 等 API 拼接

---

## 9. Mathlib API 依赖（待确认）

| 需要 | 可能的 Mathlib 名 | 确认状态 |
|------|-------------------|---------|
| char p Frobenius | `sub_pow_char` | 大概率已有 |
| a ∣ b → a ∣ b^n | `dvd_pow` | 已有 |
| monic 除法后度数 | `Polynomial.natDegree_mul` + 推导 | 已有 natDegree_mul |
| normalize associated | `normalize_associated` | 大概率已有 |
| dvd ↔ modByMonic = 0 | `Polynomial.dvd_iff_modByMonic_eq_zero` | Phase 0 已确认 |
| modByMonic + div = orig | `Polynomial.modByMonic_add_div` | Phase 0 已确认 |
| squarefree 因子也 squarefree | `Squarefree.squarefree_of_dvd` | Phase 2 已使用 |
| irreducible → prime (UFD) | `Irreducible.prime` | 大概率已有 |
| prime 整除列表积 | 归纳 `Prime.dvd_mul` | 需自行归纳或查找 |
| Associated 保持 Irreducible | `Irreducible.associated` 或反向 | 大概率已有 |
