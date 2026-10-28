# Phase 3 T3.1: DDF 算法建模与正确性证明完成

> 日期：2026-03-19
> 分支：`feature/formal-proofs`

---

## 做了什么

在 `CLPoly/Algorithm/DDF.lean` 中完成 DDF（Distinct Degree Factorization）算法的 L2 模型定义与完整正确性证明，**435 行 Lean，0 sorry**。这是 Phase 3（L2 算法模型）的第一个定理，也是首次对 CLPoly 的递归算法进行形式化验证。

### 新增文件
- `proof/lean/CLPoly/Algorithm/DDF.lean` — L2 DDF 算法模型 + 正确性
- `proof/nl-proof/phase3-t31-ddf.md` — 自然语言证明草稿 v2（审核修正版）

### 核心内容

#### 1. `ddfLoop` 递归函数（~25 行定义 + ~20 行终止性）
- 对应 C++ `__ddf_Zp`（polynomial_factorize_zp.hh:247-292）
- 终止度量：`f_star.natDegree + 1 - 2 * d`（Nat 截断减法）
- `decreasing_by`：split case 用 `natDegree_divByMonic_lt`（自定义引理），no-split case 用 `omega`

#### 2. 辅助引理（~120 行）
| 引理 | 内容 | 行数 |
|------|------|------|
| `natDegree_divByMonic_lt` | Monic g 整除 f 且 deg(g)>0 → deg(f/ₘg) < deg(f) | ~15 |
| `monic_divByMonic_mul_eq` | Monic g ∣ f → g * (f /ₘ g) = f | ~5 |
| `prime_dvd_list_prod` | Prime q ∣ l.prod → ∃ a ∈ l, q ∣ a | ~10 |
| `h_cong_step` | h-congruence 递推（Frobenius 替代路径） | ~25 |
| `gd_irred_characterization` | q ∣ gd ↔ q ∣ f* ∧ deg(q) = d（核心桥接） | ~55 |

#### 3. `ddfLoop_correct` 主定理（~175 行）
Functional induction（`ddfLoop.induct`）+ 携带 7 个不变量（P0-P6）：

| Case | 描述 | 关键技术 |
|------|------|---------|
| 1 | Termination, deg>0 | `WfDvdMonoid.exists_irreducible_factor` + 度数计数 |
| 2 | Termination, deg=0 | `eq_one_of_monic_natDegree_zero` |
| 3 | Split (gd.deg>0) | 全部 7 个不变量保持，Squarefree 对角论证 |
| 4 | No-split | `gd_irred_characterization` 反向 + 度数矛盾 |

#### 4. `ddf_correct` 顶层定理（~15 行）
```
theorem ddf_correct (f : Polynomial (ZMod p)) (hm : Monic f) (hsq : Squarefree f) :
    DDFCorrect f (ddf f)
```

## 为什么做

Phase 3 路线图：T3.1 DDF → T3.2 SQF → T3.3 EDF → T3.4 Zp pipeline assembly。DDF 是 Zp 因式分解管线的核心子过程，Phase 1 骨架证明（`FactorZp.lean`）以 `DDFCorrect` 为假设。本步骤将假设替换为具体算法的机器检查证明。

## 关键决策

1. **`sub_dvd_pow_sub_pow` 替代 Frobenius**：nl-proof 原方案用 `sub_pow_char`（char p 下 Frobenius），实际形式化时发现 `sub_pow_char` 的参数传递复杂（需要 `CharP (Polynomial (ZMod p)) p` 实例）。改用 `Commute.sub_dvd_pow_sub_pow`（`(a-b) ∣ (a^n - b^n)` 在交换环中成立），一步完成，不依赖特征。

2. **Functional induction（`ddfLoop.induct`）**：Lean 4 自动生成的归纳原理，匹配递归函数的分支结构。需要 `revert` 不变量后再 induction，使 IH 携带所有不变量的量化。`rename_i` 需要精确计数 inaccessible names（含 `let` 绑定：Case 3 有 10 个，Case 4 有 9 个）。

3. **Squarefree 直接使用定义**：P5 split case 需要 `IsCoprime f_new gd`，nl-proof 规划了 `squarefree_mul_coprime` 辅助引理。实际形式化时直接用 `Squarefree` 的定义（`q * q ∣ f → IsUnit q`）+ `mul_dvd_mul` 构造 `q² ∣ f_star`，避免了额外引理。

## 遇到的问题

1. **`rename_i` 计数**：`ddfLoop.induct` 的 Case 3（split）比 Case 4（no-split）多 1 个 inaccessible name（`f_new` let 绑定）。调试需要逐步尝试不同数量的名称。
2. **`hF3 ▸` 过度替换**：`▸` 替换目标中所有出现的 `fv`（包括 `fv /ₘ gdv` 中的），导致类型错误。解决：改用显式 dvd witness `⟨gdv, hF3.trans (mul_comm _ _)⟩`。
3. **`dv + 1 - 1 ≠ dv`（Nat 截断）**：IH 的 h-congruence 条件中出现 `p ^ (dv + 1 - 1)`，需要手动 `rw [show dv + 1 - 1 = dv by omega]`。

## 涉及的文件

| 文件 | 变更 |
|------|------|
| `proof/lean/CLPoly/Algorithm/DDF.lean` | **新增** 435 行 |
| `proof/lean/CLPoly.lean` | 添加 `import CLPoly.Algorithm.DDF` |
| `proof/nl-proof/phase3-t31-ddf.md` | **新增** nl-proof 草稿 v2 |

## 度量

- 耗时：~6 小时（含 nl-proof 草稿 + 审核 + 形式化 + 调试）
- 迭代：~30 轮编译-修复循环
- Lean 新增行数：435 行
- 对应 C++ 行数：~45 行（`__ddf_Zp` 函数 + `__upoly_powmod`/`__upoly_subtract_x` 辅助）
- 放大系数：~9.7x（435/45）
- 放弃的方案：
  - `sub_pow_char` Frobenius 路径（参数传递困难，改用 `sub_dvd_pow_sub_pow`）
  - `squarefree_mul_coprime` 独立引理（直接用 Squarefree 定义更简洁）
  - `simp only [ddfLoop, ...]` 展开递归（无限循环，改用 `unfold ddfLoop` + 手动 `dif_pos/dif_neg`）
