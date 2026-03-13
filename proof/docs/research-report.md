# Lean 4 形式化验证：深度调研报告

基于 `implementation-roadmap.md` 的实施路线，对关键技术问题进行系统调研。

---

## 1. Mathlib API 精确评估

### 1.1 `Polynomial (ZMod p)` 的 typeclass 实例

| 实例 | 可用？ | 来源 | 备注 |
|------|-------|------|------|
| `CommRing` | ✅ | `Mathlib.Data.ZMod.Basic` | 对任意 `n` 均有 |
| `Field`（`p` 素数） | ✅ | `ZMod.instField` | 需 `Fact (Nat.Prime p)` |
| `EuclideanDomain`（多项式） | ✅ | `Mathlib.Algebra.Polynomial.FieldDivision` | 仅当系数环是域时 |
| `GCDMonoid`（多项式） | ✅ | 继承自 `EuclideanDomain` | `Polynomial.gcd` 可用 |
| `DecidableEq (ZMod p)` | ✅ | `ZMod.decidableEq` | 多项式的 `DecidableEq` 随之可用 |
| `Fintype (ZMod p)` | ✅ | `ZMod.instFintype` | `Fintype.card (ZMod p) = p` |

**结论**：Phase 1-2 所需的 `Polynomial (ZMod p)` 基础设施完备。`gcd` 可计算，`divByMonic`/`modByMonic` 可用。

### 1.2 有限域核心定理

#### 定理 2.1 的 Mathlib 原料

Mathlib 中**没有**直接等价于 "x^{p^d} - x = ∏{首一不可约 g : deg g | d}" 的定理。但有以下原料：

| Mathlib 定理 | 精确内容 | 可用性 |
|-------------|---------|-------|
| `FiniteField.roots_X_pow_card_sub_X` | `(X^q - X).roots = Finset.univ.val` where `q = #K` | ✅ 给出根的完整刻画 |
| `FiniteField.splits_X_pow_nat_card_sub_X` | `X^(Nat.card K) - X` splits over `K` | ✅ 分裂性 |
| `GaloisField` 定义 | `GaloisField p n` = `X^(p^n) - X` 在 `ZMod p` 上的分裂域 | ✅ |
| `FiniteField.nonempty_algHom_of_finrank_dvd` | `[K:F] ∣ [L:F]` ⟹ ∃ 嵌入 `K →ₐ[F] L` | ✅ 域嵌入 |
| Frobenius 自同构 (`frobeniusAlgEquiv`) | `F_{p^d}` 上的 Frobenius 是域自同构 | ✅ |
| `FiniteField.pow_card` | `a ^ #K = a`（Fermat 小定理推广） | ✅ |

**桥接策略**（~100-150 行）：
1. `roots_X_pow_card_sub_X` → 根 = 全部 $\mathbb{F}_{p^d}$ 元素
2. 不可约 $g$, $\deg g = k \mid d$ → $g$ 的根在 $\mathbb{F}_{p^k} \hookrightarrow \mathbb{F}_{p^d}$ 中（用 `nonempty_algHom_of_finrank_dvd`）
3. → $g \mid (x^{p^d} - x)$
4. $x^{p^d}-x$ 无平方（导数 $= -1 \neq 0$）+ UFD → 精确乘积分解

**缺失的关键步骤**：从"根"视角到"不可约因子乘积"视角的桥接。需要证明：在 UFD 中，无平方多项式 = 其不可约因子的无重复乘积。这在 Mathlib 中有 `UniqueFactorizationMonoid` 相关工具，但需要手动组装。

#### Squarefree ↔ 导数互素

| 定理 | 方向 | 条件 |
|------|------|------|
| `Polynomial.Separable.squarefree` | `Separable f → Squarefree f` | 任意 `CommSemiring` |
| `PerfectField.separable_iff_squarefree` | `Squarefree f ↔ Separable f` | 完美域（包括 `ZMod p`） |
| `Polynomial.Separable` 定义 | `IsCoprime f (derivative f)` | — |

**结论**：`Squarefree f ↔ IsCoprime f f'` 在 `ZMod p` 上完全可用。

### 1.3 `ZMod (p^k)` 的可用性（Phase 3 Hensel 提升）

| 问题 | 答案 | 影响 |
|------|------|------|
| `CommRing` 实例？ | ✅ 有 | 多项式环可用 |
| `IsDomain`？ | ❌ 否（k ≥ 2 有零因子） | 不能用域特有的定理 |
| `Polynomial.divByMonic` 可用？ | ✅ **仅需 `Ring`，不需要 `IsDomain`** | **关键发现**：首一除法在非整环上也工作 |
| `IsUnit a ↔ Nat.Coprime a.val (p^k)`？ | ✅ `ZMod.coe_int_isUnit_iff_isCoprime` | 可构造模逆元 |
| `ZMod.unitOfCoprime`？ | ✅ 有 | 从互素构造单位 |
| `ZMod.castHom`（层间投影）？ | ✅ 有 | `ZMod (p^k) → ZMod (p^j)` 当 `j ∣ k` |
| `Polynomial.map` 配合 `castHom`？ | ✅ 有 | 多项式系数约化 |
| 对称模表示 `sym_mod`？ | ❌ **Mathlib 无** | 需自定义 ~30 行 |

**关键发现**：`Polynomial.divByMonic` 只要求 `Ring` 实例，不要求 `IsDomain`。这意味着 Hensel 提升中的首一多项式除法在 `ZMod (p^k)` 上**直接可用**，无需额外证明。对应半形式化证明引理 1.1 的核心内容。

**sym_mod 定义**（需自建）：
```lean
def sym_mod (a : ZMod m) (hm : 0 < m) : ℤ :=
  let v := (ZMod.val a : ℤ)
  if v ≤ m / 2 then v else v - m
```

### 1.4 总结：Mathlib 缺口修正

| 蓝图中的需求 | 原预估 | 修正后预估 | 原因 |
|-------------|--------|-----------|------|
| Thm 2.1 桥接 | ~100 行 | **~100-150 行** | 根→因子乘积桥接比预期稍复杂，需要 UFD 工具 |
| Cor 2.2 | ~50 行 | ~50 行 | 不变 |
| 引理 1.1（首一除法） | ~50 行 | **~10 行** | `divByMonic` 已内置，仅需包装 |
| Hensel 唯一性 | ~200-300 行 | **~300-500 行**（见 §2） | 根据 Isabelle 经验调高 |
| Mignotte 界 | ~100 行 | ~100 行 | 不变 |
| sym_mod | 未计划 | **~30 行** | 新增 |

### 1.5 UFD 分解与不可约因子 API

Mathlib 提供了完整的 UFD 分解设施，**不需要构建任何基础理论**。

#### `normalizedFactors`

```lean
def UniqueFactorizationMonoid.normalizedFactors (a : α) : Multiset α
```

返回 `a` 的全部不可约因子的 `Multiset`，每个因子已归一化。对 `Polynomial (ZMod p)`，"归一化"即**首一**（`normalize f = f * C (leadingCoeff f)⁻¹`）。

#### 核心定理

| 定理 | 签名 | 用途 |
|------|------|------|
| `prod_normalizedFactors` | `Associated (normalizedFactors a).prod a` | 乘积还原（差一个单位） |
| `irreducible_of_mem_normalizedFactors` | `p ∈ normalizedFactors a → Irreducible p` | 每个因子不可约 |
| `squarefree_iff_nodup_normalizedFactors` | `Squarefree x ↔ (normalizedFactors x).Nodup` | 无平方 ↔ 无重复 |
| `dvd_iff_normalizedFactors_le_normalizedFactors` | `x ∣ y ↔ normalizedFactors x ≤ normalizedFactors y` | 整除 ↔ 多重集包含 |
| `normalizedFactors_mul` | `normalizedFactors (x * y) = normalizedFactors x + normalizedFactors y` | 乘积因子分解 |

#### 过滤与乘积

可以直接写：
```lean
(normalizedFactors f).filter (fun g => g.natDeg = d) |>.prod
```
`Multiset.filter` + `Multiset.prod` 组合可用。

#### 实际使用建议

虽然 `normalizedFactors` 可用，但 DDF 正确性定理**不必依赖它**。用纯整除性谓词更简洁：

```lean
-- "gd 恰为 f 中全部 d 次不可约因子之积" 等价于：
gd ∣ f ∧ Monic gd ∧ Squarefree gd
∧ (∀ q, Irreducible q → q ∣ gd → q.natDeg = d)
∧ (∀ q, Irreducible q → q ∣ f → q.natDeg = d → q ∣ gd)
```

两种表述等价，后者避免了操作 `Multiset` 的证明开销。

### 1.6 `Polynomial.gcd` 归一化行为

对 `Polynomial K`（`K` 为域），`gcd` 返回**已归一化（首一）**的结果：

```lean
-- Polynomial K 有 NormalizedGCDMonoid 实例
-- 所以 normalize (gcd f g) = gcd f g
-- 即 gcd 直接返回首一多项式
```

这意味着 DDF 中 `gcd(h-x, f*)` 的结果可以直接用于 `divByMonic`，无需额外归一化。

### 1.7 `Monic` 与 `divByMonic` API

```lean
-- Monic 是 Prop，不是子类型
def Polynomial.Monic (p : Polynomial R) : Prop := leadingCoeff p = 1

-- divByMonic 不需要 Monic 证明作为参数（unchecked）
def Polynomial.divByMonic (p q : Polynomial R) : Polynomial R
-- 若 q 非首一则返回 0

-- 关键恒等式（仅当 Monic q 时成立）
theorem modByMonic_add_div (p : Polynomial R) {q} (hq : Monic q) :
    p %ₘ q + q * (p /ₘ q) = p
```

`divByMonic` 是 unchecked 的——不需要传入 `Monic` 证明，但行为仅在首一时正确。

### 1.8 `content` / `IsPrimitive` / Gauss 引理（Z[x]）

Mathlib 完整支持：

```lean
-- 内容
def Polynomial.content (p : Polynomial ℤ) : ℤ  -- GCD of all coefficients

-- 本原
def Polynomial.IsPrimitive (p : Polynomial ℤ) : Prop  -- content = 1

-- 本原部分
def Polynomial.primPart (p : Polynomial ℤ) : Polynomial ℤ  -- p / content(p)

-- 分解
theorem Polynomial.eq_content_mul_primPart (p : Polynomial ℤ) :
    p = C (content p) * primPart p

-- Gauss 引理（本原 × 本原 = 本原）
theorem Polynomial.IsPrimitive.mul (hp : IsPrimitive p) (hq : IsPrimitive q) :
    IsPrimitive (p * q)

-- 内容可乘性
theorem Polynomial.content_mul (p q : Polynomial ℤ) :
    (p * q).content = p.content * q.content
```

**结论**：Phase 3 需要的 `pp(sym_m(P)) = g_j` 中的 `pp`（本原部分）对应 `Polynomial.primPart`，Gauss 引理已形式化。无需自建。

---

## 2. Hensel 引理形式化：深度分析

### 2.1 现有形式化成果

| 系统 | 内容 | 类型 | 状态 |
|------|------|------|------|
| **Lean 4 / Mathlib** | `hensels_lemma` in `Mathlib.NumberTheory.Padics.Hensel` | **标量根提升**（不是多项式因子提升） | ✅ 已有但不适用 |
| **Isabelle / AFP** | Berlekamp-Zassenhaus entry | **完整多项式因子提升**（二进制 + 二次 + 多因子 + 平衡树） | ✅ 完整 |
| **Coq / SSReflect** | Martin-Dorel et al. | 单变量 + 双变量 Hensel | 部分 |

### 2.2 Mathlib 的 `hensels_lemma` — 不可直接用

Mathlib 的 Hensel 引理是关于 **p-adic 整数上的标量根提升**：

> 给定 $F \in \mathbb{Z}_p[x]$ 和近似根 $a$ 满足 $\|F(a)\| < \|F'(a)\|^2$，存在唯一精确根 $z$ 使得 $F(z) = 0$。

这是 Newton 迭代的收敛性证明，与我们需要的**多项式因子提升**（$f \equiv g \cdot h \pmod{p}$ 提升到 $f \equiv G \cdot H \pmod{p^k}$）本质不同。

**结论**：不能复用 Mathlib 的 `hensels_lemma`，需从零证明多项式版本。

### 2.3 Isabelle AFP 的经验（核心参考）

Divasón、Thiemann 等人的 Isabelle 形式化是**最重要的参考**：

**覆盖范围**：
- Berlekamp 算法（$\mathbb{F}_p$ 上因式分解）
- 二进制 Hensel 提升（线性：$p^i \to p^{i+1}$）
- 二次 Hensel 提升（$p^i \to p^{2i}$）
- 多因子提升（$n$-ary，使用因子树）
- 平衡因子树（度数均衡，优化效率）
- Mignotte 界 + Graeffe 变换
- 因子重组

**关键技术决策**：
- 类型问题：$\mathbb{F}_p$ 和 $\mathbb{Z}/p^k\mathbb{Z}$ 的模数 $p$、$k$ 是运行时值，无法编码在类型中
- 解决方案：使用**整数表示**（在 $\mathbb{Z}$ 上运算，显式携带模数参数），而非类型级表示
- 性能：因式分解 500 次多项式在秒级完成，仅比 Mathematica 慢 2.5 倍

**对 CLPoly 的启示**：
- Lean 4 有更好的依赖类型，`ZMod n` 的 `n` 可以是变量，可能比 Isabelle 处理得更优雅
- 但 Hensel 提升中 $k$ 逐步增加，$\mathbb{Z}/p^{k}\mathbb{Z}$ 的类型在每一步都不同，仍需仔细处理
- **可能的替代方案**：在 $\mathbb{Z}$ 上做算术，用 `%` 操作表示模约化，避免类型级模数变化

### 2.4 多项式 Hensel 引理的技术挑战

| 挑战 | 描述 | 应对 |
|------|------|------|
| **非整环** | $\mathbb{Z}/p^k\mathbb{Z}$ 有零因子，唯一性证明更难 | 利用首一约束 + `divByMonic` 的环级可用性 |
| **Bézout 系数度界** | $sG + tH \equiv 1 \pmod{p}$ 中 $\deg s < \deg H$, $\deg t < \deg G$ | 需显式形式化度数界 |
| **唯一性需要度约束** | 没有度约束时，提升不唯一 | 形式化时显式传递度不变量 |
| **多因子版本** | 从二因子递推到多因子，需要因子树结构 | 可推迟到 Phase 3 后期 |
| **类型变化** | 每步提升 $k \to 2k$，`ZMod (p^k)` 类型变化 | 选择：(a) 在 `ℤ` 上工作 + 显式模约化，或 (b) 使用 `ZMod` 但需要 `cast` |

### 2.5 工作量估算修正

基于 Isabelle 经验（de Bruijn factor ~17，即非形式化到形式化的膨胀率约 17 倍）：

| 组件 | 半形式化行数 | 预估 Lean 4 行数 | 备注 |
|------|------------|-----------------|------|
| 二因子 Hensel 唯一性（Thm 2.1） | ~15 行 | **300-500 行** | 核心定理，需精细处理度约束 |
| 多因子 Hensel 提升（Thm 2.2） | ~20 行 | **200-400 行** | 从二因子归纳 |
| Bézout 系数提升 | ~10 行 | **150-250 行** | 辅助，但不可跳过 |
| 因子树结构 | 新增 | **200-300 行** | CLPoly 使用二叉树提升 |
| **Hensel 子系统总计** | ~45 行 | **~1000-1500 行** | |

**与路线图对比**：原预估 200-300 行明显偏低。修正为 **1000-1500 行**，时间从 1-2 周调整为 **2-3 周**。

### 2.6 替代方案评估

| 方案 | 描述 | 优劣 |
|------|------|------|
| **A: ZMod 类型级** | 在 `Polynomial (ZMod (p^k))` 上做所有运算 | 类型安全但每步类型变化，需大量 `cast` |
| **B: ℤ 上 + 显式模** | 在 `Polynomial ℤ` 上运算，所有等式模 $m$ 表述 | Isabelle 的做法，灵活但丢失类型约束 |
| **C: axiom 占位** | 先 `axiom hensel_uniqueness`，完成上层证明，后补 | 快速推进，风险是占位可能永远不补 |
| **D: p-adic 路径** | 在 $\mathbb{Z}_p$ 上证明，投影到 $\mathbb{Z}/p^k\mathbb{Z}$ | 理论优雅但失去算法结构，不利于 1:1 对应 |

**推荐**：Phase 3 初期用方案 C 快速推进，同时并行用方案 A 或 B 逐步填充。

---

## 3. Lean 4 递归与终止性机制

### 3.1 `termination_by` — DDF 循环建模的核心工具

DDF 的循环结构需要翻译为递归 + 终止性证明。关键问题：`f_star.natDeg - 2 * d` 可以作为终止度量吗？

**答案：可以。**

```lean
def ddf_loop (h f_star : ZpPoly) (d : ℕ) (acc : List (ZpPoly × ℕ)) : List (ZpPoly × ℕ) :=
  if f_star.natDeg < 2 * d then ...  -- 终止
  else
    ...
    if 0 < gd.natDeg then
      ddf_loop h'' f_star' (d + 1) (acc ++ [(gd, d)])  -- f_star' 度数减小
    else
      ddf_loop h' f_star (d + 1) acc  -- d 增大，2*d 趋近 f_star.natDeg
termination_by f_star.natDeg + 1 - 2 * d
```

**关键点**：
- 自然数减法 `n - m` 在 `n < m` 时返回 0，天然有界
- 两个分支都使度量递减：
  - 非平凡 gcd：`f_star'.natDeg < f_star.natDeg`，所以 `f_star'.natDeg + 1 - 2*(d+1) < f_star.natDeg + 1 - 2*d`
  - 平凡 gcd：`d+1 > d`，所以 `f_star.natDeg + 1 - 2*(d+1) < f_star.natDeg + 1 - 2*d`
- `omega` tactic 通常能自动证明这类线性算术

### 3.2 `decreasing_by` — 何时需要

```lean
termination_by f_star.natDeg + 1 - 2 * d
decreasing_by
  all_goals simp_wf
  · -- 分支 1：非平凡 gcd
    -- 需证明 f_star'.natDeg + 1 - 2 * (d+1) < f_star.natDeg + 1 - 2 * d
    omega  -- 通常能自动解决
  · -- 分支 2：平凡 gcd
    omega
```

**经验法则**：先试 `decreasing_by all_goals simp_wf; omega`。只有 `omega` 失败时才写手动证明。

### 3.3 `partial def` — EDF 的特殊处理

EDF（Cantor-Zassenhaus）是概率算法，理论上可能不终止（虽然概率为 0）。

**关键限制**：`partial def` 的函数**无法证明任何性质**。

**替代方案**：

| 方案 | 实现 | 可证明性 |
|------|------|---------|
| `partial def` | 最简单，直接写 | ❌ 不可证明 |
| 燃料参数 `fuel : ℕ` | `if fuel = 0 then none else ...` | ✅ 可证明（返回 `Option`） |
| 外部终止性证明 | `def edf ... (h : terminates ...) : ...` | ✅ 完全可证明 |
| 概率终止文档论证 | `partial def` + 外部文档 | ❌ 形式化层面不可证 |

**推荐**：
- **Phase 2 初版**：用 `partial def`，附文档论证期望迭代次数 ≤ 2
- **Phase 2 加强版**：改为燃料参数，证明"给定足够燃料，输出正确"
- `factor_Zp_correct` 顶层定理中，EDF 的终止性作为 `sorry` 或显式假设

### 3.4 `let rec` 与嵌套递归

**可以在 `def` 内部使用 `let rec`**：

```lean
def ddf_Zp (f : ZpPoly) : List (ZpPoly × ℕ) :=
  let rec loop (h f_star : ZpPoly) (d : ℕ) (acc : List (ZpPoly × ℕ)) :=
    ...
  termination_by f_star.natDeg + 1 - 2 * d
  loop X f 1 []
```

`termination_by` 直接附在 `let rec` 上，Lean 4.6+ 完全支持。

### 3.5 函数归纳原理 — 证明利器

Lean 4 为每个递归函数自动生成 `.induct` 归纳原理：

```lean
-- 自动生成
ddf_loop.induct : ∀ (motive : ZpPoly → ZpPoly → ℕ → List ... → Prop),
  (∀ h f_star d acc, f_star.natDeg < 2 * d → motive h f_star d acc) →  -- 终止分支
  (∀ h f_star d acc, ¬(f_star.natDeg < 2 * d) → ... →                   -- 非平凡 gcd
    motive h'' f_star' (d+1) ... → motive h f_star d acc) →
  (∀ h f_star d acc, ¬(f_star.natDeg < 2 * d) → ... →                   -- 平凡 gcd
    motive h' f_star (d+1) acc → motive h f_star d acc) →
  ∀ h f_star d acc, motive h f_star d acc
```

**用法**：
```lean
theorem ddf_loop_invariant ... :=
  ddf_loop.induct (motive := fun h f_star d acc => ...)
    (fun h f_star d acc hlt => ...)    -- 终止分支
    (fun h f_star d acc _ _ ih => ...) -- 非平凡 gcd，ih 是归纳假设
    (fun h f_star d acc _ _ ih => ...) -- 平凡 gcd
```

**优势**：证明结构自动匹配函数的递归结构，不需要手动重复 `if-then-else` 分析。

### 3.6 各函数的终止性策略汇总

| 函数 | 递归结构 | 终止度量 | 方案 |
|------|---------|---------|------|
| `upoly_powmod` | while (e > 0), e >>= 1 | `e` | `termination_by e` + `omega` |
| `ddf_Zp` | for d = 1, ..., break when deg < 2d | `f_star.natDeg + 1 - 2 * d` | `termination_by` + `omega` |
| `edf_Zp` | while(true), 随机分裂 | 概率终止 | `partial def` 或 燃料 |
| `squarefree_Zp` | while (deg w > 0) + 递归（p-th root） | `f.natDeg`（外层递归减小） | 嵌套：外层 `termination_by natDeg`，内层 `termination_by` 计数 |
| `factor_Zp` | 嵌套 for 循环 | 非递归（遍历列表） | 无需终止性证明（List.map / for-in） |
| `hensel_lift` | while (m ≤ target), m ← m² | `target - m`（粗略）或 $\lceil\log_2(target/m)\rceil$ | `termination_by` |
| `zassenhaus_recombine` | for s = 1..n/2, 内层组合遍历 | 活跃集 `T` 大小单调递减 | `termination_by T.card` |

---

## 4. 现有形式化项目全景

### 4.1 Isabelle / AFP — 最完整的参考

**Berlekamp-Zassenhaus 形式化**（Divasón, Joosten, Thiemann, Yamada）：
- **论文**：CPP 2017 + Journal of Automated Reasoning 2019
- **AFP 入口**：`Berlekamp_Zassenhaus`
- **覆盖**：Berlekamp + Hensel（二进制/二次/多因子/平衡树）+ 因子重组 + Mignotte 界
- **性能**：500 次多项式秒级因式分解，比 Mathematica 慢 2.5 倍
- **关键技术**：用整数表示绕过类型级模数问题

**LLL 基约化 + LLL 因式分解**（Bottesch, Divasón, Haslbeck, Thiemann 等）：
- **论文**：ITP 2018 + Journal of Automated Reasoning 2020
- **AFP 入口**：`LLL_Basis_Reduction` + `LLL_Factorization`
- **规模**：**~14,811 行 Isabelle 代码**，de Bruijn factor ~17
- **开发时间**：~23 人月
- **重要发现**：形式化过程中**发现并修正了教科书中的错误**
- **van Hoeij 尚未形式化**，但作者指出可基于已有 LLL 工作构建

### 4.2 Coq

- **Martin-Dorel et al.**：Hensel 引理的 Coq/SSReflect 形式化（单变量 + 双变量）
- **Kirkels**：不可约性证书（Magma 计算 + Coq 验证的混合方法）
- **CoqEAL**：精化式多项式计算库（MathComp 多项式到 seq 表示的精化）

### 4.3 Lean 4

- **Hensel 引理**：仅 p-adic 标量版（Robert Y. Lewis, CPP 2019），已集成 Mathlib
- **Gröbner 基**：正在进行（arXiv 2602.12772, 2025），覆盖 Buchberger 算法
- **多项式因式分解**：**尚无任何形式化**

### 4.4 尚无形式化的领域

| 算法 | 状态 | 对 CLPoly 的影响 |
|------|------|-----------------|
| DDF（按度分组） | ❌ 任何证明器中均无 | Phase 1 目标，**首个形式化** |
| EDF（Cantor-Zassenhaus） | ❌ 任何证明器中均无 | Phase 2 目标 |
| van Hoeij 算法 | ❌ 任何证明器中均无 | Phase 3-4 目标 |
| 多变量因式分解 | ❌ 任何证明器中均无 | 超出当前范围 |

**意义**：CLPoly 的 Phase 1（DDF 形式化）将是**全球首个** DDF 算法的机器检查形式化。

### 4.5 从 Isabelle 经验中学到的

| 教训 | 内容 | 对 CLPoly 的影响 |
|------|------|-----------------|
| **de Bruijn factor ~17** | 1 页非形式化 → ~17 页形式化代码 | 半形式化证明 ~30 页 → 预估 Lean **~500 页 / ~15,000 行** |
| **类型级模数是大坑** | Isabelle 被迫用整数表示 | Lean 4 的 `ZMod n`（n 可为变量）可能更好，但 Hensel 提升中 k 变化仍需处理 |
| **平衡因子树很重要** | 效率差异显著 | CLPoly 已使用二叉树，建模时需如实反映 |
| **发现教科书错误** | LLL 形式化过程中发现 | 形式化本身有价值，可能发现 CLPoly 代码中的隐藏问题 |
| **开发时间 ~23 人月** | LLL 子系统 | 整个管线可能需要 **6-12 人月** |

---

## 5. 路线图修正建议

### 5.1 工作量修正

| Phase | 原预估 | 修正预估 | 原因 |
|-------|--------|---------|------|
| Phase 0 | 1 周 | **1 周** | 不变 |
| Phase 1 (DDF) | 2-3 周 | **3-4 周** | Thm 2.1 桥接可能更复杂；函数归纳证明有学习曲线 |
| Phase 2 (Zp 管线) | 3-4 周 | **3-4 周** | EDF 用 `partial def` 可控 |
| Phase 3 (Hensel + Z[x]) | 4-6 周 | **6-10 周** | Hensel 引理从 200-300 行修正为 1000-1500 行 |
| Phase 4 (LLL) | 可选 | 可选（若做，**8-12 周**） | 参考 Isabelle 23 人月 |

### 5.2 关键风险更新

| 风险 | 严重性 | 缓解 |
|------|--------|------|
| **Hensel 引理工作量远超预期** | 🔴 高 | 方案 C（axiom 占位）+ 并行推进上层 |
| **类型级 vs 整数级模数选择** | 🟡 中 | Phase 0 做小规模实验对比两种方案 |
| **`partial def` 限制 EDF 可证明性** | 🟡 中 | 燃料参数替代；顶层定理中显式假设 EDF 终止 |
| **Lean 4 + Mathlib 编译时间** | 🟡 中 | `sorry` 占位加速迭代；分模块编译 |
| **Thm 2.1 根→因子乘积桥接** | 🟢 低 | UFD 工具链在 Mathlib 中较成熟 |

### 5.3 Phase 0 新增任务

基于调研结果，Phase 0 应增加以下验证实验：

1. **Hensel 类型实验**：在 `ZMod (p^k)` 上做一个 toy Hensel 提升（二因子、一步），验证 `divByMonic` 和 `castHom` 的实际使用体验
2. **Thm 2.1 桥接实验**：尝试从 `roots_X_pow_card_sub_X` 推导出因子乘积形式的前 2 步
3. **终止性实验**：写一个类似 DDF 的 toy 递归函数，验证 `termination_by` + `decreasing_by omega` 的工作情况
4. **编译时间基准**：导入完整 Mathlib 后测量 `lake build` 时间

---

## 6. 参考文献

### 6.1 Isabelle 形式化

- Divasón, Joosten, Thiemann, Yamada. "A Verified Implementation of the Berlekamp–Zassenhaus Factorization Algorithm." *J. Automated Reasoning*, 2019. [Springer](https://link.springer.com/article/10.1007/s10817-019-09526-y)
- Divasón et al. "A Formalization of the Berlekamp-Zassenhaus Factorization Algorithm." *CPP 2017*. [PDF](http://cl-informatik.uibk.ac.at/users/thiemann/paper/CPP17_Berlekamp_Zassenhaus.pdf)
- Bottesch, Divasón, Haslbeck, Thiemann, Joosten, Yamada. "Formalizing the LLL Basis Reduction Algorithm and the LLL Factorization Algorithm in Isabelle/HOL." *J. Automated Reasoning* 64(5), 2020. [Springer](https://link.springer.com/article/10.1007/s10817-020-09552-1)
- AFP: [Berlekamp_Zassenhaus](https://www.isa-afp.org/entries/Berlekamp_Zassenhaus.html), [LLL_Basis_Reduction](https://www.isa-afp.org/entries/LLL_Basis_Reduction.html), [LLL_Factorization](https://www.isa-afp.org/entries/LLL_Factorization.html)

### 6.2 Lean 4 / Mathlib

- Lewis. "A Formal Proof of Hensel's Lemma over the p-adic Integers." *CPP 2019*. [PDF](https://robertylewis.com/padics/padics.pdf)
- Lean 4 递归文档：[Well-Founded Recursion](https://lean-lang.org/doc/reference/latest/Recursive-Definitions/Well-Founded-Recursion/)
- Lean 4 函数归纳：[Blog post](https://lean-lang.org/blog/2024-5-17-functional-induction)
- Gröbner 基形式化：[arXiv 2602.12772](https://arxiv.org/abs/2602.12772)

### 6.3 Coq

- Martin-Dorel et al. "Formalization of Hensel's lemma in Coq." [HAL](https://ens-lyon.hal.science/ensl-00560449)
- Kirkels. "Irreducibility Certificates for Polynomials with Integer Coefficients." [MSc thesis](https://www.math.ru.nl/~bosma/students/kirkels/mscthesis.pdf)
- CoqEAL: [GitHub](https://github.com/rocq-community/coqeal)

### 6.4 Mathlib 关键模块

| 模块 | 用途 |
|------|------|
| `Mathlib.Data.ZMod.Basic` | `ZMod n` 基本定义和实例 |
| `Mathlib.Data.ZMod.Coprime` | 单位 ↔ 互素 |
| `Mathlib.Algebra.Polynomial.Div` | `divByMonic` / `modByMonic`（Ring 级） |
| `Mathlib.Algebra.Polynomial.FieldDivision` | `EuclideanDomain` 实例 |
| `Mathlib.FieldTheory.Finite.Basic` | Frobenius, `pow_card` |
| `Mathlib.FieldTheory.Finite.GaloisField` | `GaloisField p n`, 分裂域 |
| `Mathlib.FieldTheory.Separable` | `Separable ↔ IsCoprime f f'` |
| `Mathlib.Algebra.Squarefree.Basic` | `Squarefree` 定义和性质 |
| `Mathlib.NumberTheory.Padics.Hensel` | p-adic Hensel（标量，不适用） |
