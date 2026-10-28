# Phase 2 完成：T2.4 EDF 三分性 + T2.5 半幂根计数

> 日期：2026-03-17
> 分支：`feature/formal-proofs`

---

## 做了什么

在 `FiniteFieldFact.lean` 中完成 T2.4 和 T2.5 的 Lean 4 形式化，Phase 2（L3 数学基石）全部 **0 sorry**。

### 新增辅助引理：`all_pow_eq_of_dvd_X_pow_sub_X`

从 T2.2' 的 Step A+B 提取为独立引理：若 `g` 不可约且 `g | X^{p^d} - X`，则 AdjoinRoot g 中所有元素满足 `x^{p^d} = x`。

- 构造 `iterateFrobenius` 为 `AlgHom`，由 `AdjoinRoot.algHom_ext`（在 root 上一致）得 φ = id
- T2.2' 的 Step A+B 简化为单行调用
- T2.4 直接复用

### T2.4: `edf_trichotomy`（~40 行）

| 步骤 | 方法 |
|------|------|
| `p^d = 2m+1` | `Nat.odd_iff` + `Nat.pow_mod` |
| `ā(ā^m-1)(ā^m+1) = 0` | `ring` 恒等式 + `all_pow_eq_of_dvd_X_pow_sub_X` |
| 三分 | `mul_eq_zero.mp` 两次 |
| 翻译回整除 | `AdjoinRoot.mk_eq_zero` + `map_sub/map_add/map_pow/map_one` |

### T2.5: `card_pow_half_eq_one`（~95 行）

| 步骤 | 方法 |
|------|------|
| `q ≥ 3` | `FiniteField.card` → `ringChar ≥ 3` → `Nat.pow_le_pow_right` |
| `1 ≠ -1` | `nth_rewrite` 精确重写 + `CharP.cast_eq_zero_iff` → `ringChar ∣ 2` 矛盾 |
| `\|A\| ≤ m` | A ⊆ `(X^m-1).roots.toFinset` + `card_roots'` + `natDegree_sub_eq_left_of_natDegree_lt` |
| `\|B\| ≤ m` | `C_neg, C_1` 将 `X^m+1` 改写为 `X^m - C(-1)` 后同上 |
| 覆盖性 | `pow_card_sub_one_eq_one` (Fermat) → `(a^m)²=1` → `mul_eq_zero` |
| 夹逼 | `\|A\|+\|B\| ≥ 2m`, `\|A\| ≤ m`, `\|B\| ≤ m` → `omega` |

## 关键技术点

1. **`FiniteField.card` 返回 `ℕ+`**：`obtain ⟨npn, _, hcard⟩` 中 `npn : ℕ+`，使用 `(npn : ℕ)` 显式转换和 `npn.pos` 获取 `0 < n`。

2. **`ring` 不识别 `C (-1)`**：多项式环中 `C (-1)` 对 `ring` 不透明。解法：先 `rw [C_neg, C_1]` 将 `C (-1)` 化为 `-(1 : Polynomial F)`，再 `ring`。

3. **`pow_card_sub_one_eq_one` 签名**：第一个显式参数是元素 `a`，非证明。正确调用：`pow_card_sub_one_eq_one a ha_ne`。

4. **`nth_rewrite` 精确控制**：证明 `1 ≠ -1` 时，`rw [h]` 会替换所有 `1` 为 `-1`，导致 `1+1=0` 变成 `-1+(-1)=0`。用 `nth_rewrite 1 [h]` 仅替换第一个。

## Phase 2 验收状态

| 定理 | 内容 | 行数 | 状态 |
|------|------|------|------|
| T2.1 | X^{p^d}-X 的 Separable 性 | ~8 | ✓ |
| T2.2 | 不可约多项式整除 X^{p^d}-X | ~18 | ✓ |
| 辅助 | `all_pow_eq_of_dvd_X_pow_sub_X` | ~26 | ✓ 本次提取 |
| T2.2' | 逆命题 | ~50 | ✓ 本次简化 |
| T2.3 | gcd 不可约因子刻画 | ~12 | ✓ |
| T2.4 | EDF 三分性 | ~40 | ✓ 本次新增 |
| T2.5 | 半幂根计数 | ~95 | ✓ 本次新增 |
| **合计** | | **~351** | **0 sorry** |

## 度量
- 耗时：~6 小时（nl-proof 草稿+审核 ~1h、T2.4 形式化 ~1h、T2.5 形式化 ~3.5h、辅助引理提取+T2.2'重构 ~0.5h）
- 迭代：T2.4 1 轮（首次编译仅 unused variable 警告）；T2.5 4 轮（`ℕ+` 类型、`C(-1)` ring 不透明、`pow_card_sub_one_eq_one` 签名、`omega` 需 `q≥3`）
- Lean 新增/修改行数：+215 行 / -31 行（净增 184 行）
- 对应 C++ 行数：N/A（纯数学定理，无直接 C++ 对应；为 DDF/EDF 算法正确性提供数学基础）
- 放弃的方案：(1) T2.5 `hB_le` 最初直接用 `ring` 处理 `X^m+1 = X^m - C(-1)`——`ring` 不识别 `C`，改用 `rw [C_neg, C_1]; ring`；(2) T2.5 `hone_ne_neg` 最初用 `rw [h]; ring`——`rw` 替换所有 `1` 导致目标变形，改用 `nth_rewrite 1 [h]`

## 涉及文件

| 文件 | 操作 |
|------|------|
| `proof/lean/CLPoly/Math/FiniteFieldFact.lean` | 新增 T2.4、T2.5；提取辅助引理；简化 T2.2' |
| `proof/nl-proof/phase2-t24-edf-trichotomy.md` | 状态更新：待形式化 → 已形式化 |
| `proof/nl-proof/phase2-t25-half-pow-roots.md` | 状态更新：待形式化 → 已形式化 |
