# 因式分解模块 Lean 4 形式化验证：实施路线（v2）

> v2: 自顶向下验证 — 先锁定接口规约，再逐层填充

## 0. 前置背景

本文档是 `lean-verification-blueprint.md`（蓝图）的可执行实施方案。

**输入文档**：
- `proof/docs/lean-verification-blueprint.md` — 验证架构（L3/L2/L1 三层 + 文件结构 + 证明义务）
- `proof/docs/pre-verification-architecture.md` — 方法论（seL4 精化链 + 分支二路径）
- `docs/design/factorization/formal-proof-ddf-edf.md` — DDF/EDF 半形式化证明（待机械化）
- `docs/design/factorization/formal-proof-univar-factorization.md` — Z[x] 因式分解半形式化证明（待机械化）

**C++ 源文件**（建模对象）：
- `clpoly/polynomial_factorize_zp.hh` — Zp[x] 因式分解管线
- `clpoly/polynomial_factorize_univar.hh` — Z[x] 因式分解管线

**当前状态**：Phase 0 完成，`lean/` 项目可编译，Mathlib API 验证通过。

---

## 核心原则：自顶向下验证

```
传统（自底向上）：       本项目（自顶向下）：
  powmod 证明              顶层正确性骨架（各子过程 = axiom）
    ↓                        ↓
  subtractx 证明           锁定每个子过程的接口规约
    ↓                        ↓
  DDF 循环证明             逐个填充子过程证明
    ↓                        ↓
  顶层定理                 消除所有 sorry/axiom
```

**为什么自顶向下**：
1. **先锁定接口** — 每个 `axiom` 的类型签名就是子过程必须满足的 spec，后续填充时目标明确
2. **早期发现设计错误** — 如果接口规约无法组合成顶层定理，说明分解方案有问题，越早发现越好
3. **可并行** — 接口锁定后，各子过程的证明彼此独立，可以并行推进
4. **渐进信心** — 每消除一个 `sorry` 就增加一份确定性，而不是全部完成才能验证

---

## 1. Phase 0：环境搭建 ✅ 已完成

- Lean 4 + Mathlib 项目可编译
- E1-E4 实验全部通过
- Mathlib Gap Report 完成
- 详见 `docs/devlog/2026-03-14-phase0-lean4-experiments.md`

---

## 2. Phase 1：顶层正确性骨架（1-2 周）

### 2.1 目标

写出 `factor_Zp_correct` 和 `factor_ZZ_correct` 的**完整定理陈述**，假设所有子过程正确，证明它们组合起来给出完整因式分解。

**Phase 1 结束时**：两个顶层定理通过编译，内部全是 `sorry`，但**定理陈述本身是精确的**。

### 2.2 工作流

```
T1.1 定义接口规约
  ↓
T1.2 Zp[x] 顶层骨架
  ↓
T1.3 Z[x] 顶层骨架
  ↓
T1.4 评估 + 锁定接口
```

### 2.3 任务详述

#### T1.1 定义所有子过程的接口规约

**文件**：`CLPoly/Spec.lean`

为每个子过程定义一个 `Prop`，描述它**做了什么**，不描述**怎么做**：

```lean
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.FieldTheory.Finite.Basic
import Mathlib.Algebra.Squarefree.Basic

open Polynomial

variable {p : ℕ} [Fact (Nat.Prime p)]

-- ============================================================
-- 子过程接口规约：每个 Prop 定义 "做了什么"
-- ============================================================

/-- SQF 规约：输入多项式 → 输出无平方因子分解 -/
def SquarefreeDecomp (f : Polynomial (ZMod p))
    (result : List (Polynomial (ZMod p) × ℕ)) : Prop :=
  -- 1. 乘积还原（到 associate）
  Associated f (∏ pr in result, pr.1 ^ pr.2)
  -- 2. 每个因子无平方且首一
  ∧ (∀ pr ∈ result, Squarefree pr.1 ∧ Monic pr.1)
  -- 3. 重数 ≥ 1
  ∧ (∀ pr ∈ result, pr.2 ≥ 1)
  -- 4. 因子两两互素
  ∧ (∀ pr₁ ∈ result, ∀ pr₂ ∈ result, pr₁ ≠ pr₂ → IsCoprime pr₁.1 pr₂.1)

/-- DDF 规约：输入首一无平方多项式 → 输出按度分组的不可约因子积 -/
def DDFCorrect (f : Polynomial (ZMod p))
    (result : List (Polynomial (ZMod p) × ℕ)) : Prop :=
  -- 1. 每个 (gd, d) 中 gd 整除 f
  (∀ pr ∈ result, pr.1 ∣ f)
  -- 2. gd 的每个不可约因子度数 = d
  ∧ (∀ pr ∈ result, ∀ q, Irreducible q → q ∣ pr.1 → q.natDeg = pr.2)
  -- 3. 不遗漏：f 的每个不可约因子都在某个 gd 中
  ∧ (∀ q, Irreducible q → q ∣ f → ∃ pr ∈ result, q ∣ pr.1)
  -- 4. 乘积还原
  ∧ Associated f (∏ pr in result, pr.1)

/-- EDF 规约：输入等度多项式 → 输出不可约因子列表 -/
def EDFCorrect (g : Polynomial (ZMod p)) (d : ℕ)
    (result : List (Polynomial (ZMod p))) : Prop :=
  -- 1. 乘积还原
  Associated g (∏ q in result, q)
  -- 2. 每个因子不可约、首一、度 = d
  ∧ (∀ q ∈ result, Irreducible q ∧ Monic q ∧ q.natDeg = d)

/-- Hensel 提升规约：模 p 因子 → 模 p^k 因子 -/
def HenselCorrect
    (f : Polynomial ℤ) (k : ℕ)
    (factors_mod_p : List (Polynomial (ZMod p)))
    (factors_mod_pk : List (Polynomial (ZMod (p ^ k)))) : Prop :=
  -- 1. 模 p^k 下乘积 ≡ f
  sorry -- 待精化
  -- 2. 模 p 下与原始因子一致
  -- 3. 度数保持

/-- 因子重组规约：Hensel 因子 → 真正的 Z[x] 因子 -/
def RecombineCorrect
    (f : Polynomial ℤ)
    (result : List (Polynomial ℤ)) : Prop :=
  -- 1. 乘积还原
  f = ∏ g in result, g
  -- 2. 每个因子不可约
  ∧ (∀ g ∈ result, Irreducible g)
```

#### T1.2 Zp[x] 顶层正确性骨架

**文件**：`CLPoly/Pipeline/FactorZp.lean`

```lean
/-- Zp[x] 因式分解的顶层正确性：
    假设 SQF、DDF、EDF 各自正确，则组合结果是完整因式分解 -/
theorem factor_Zp_correct
    (f : Polynomial (ZMod p)) (hf : f ≠ 0)
    -- 假设各子过程存在且正确
    (sqf : Polynomial (ZMod p) → List (Polynomial (ZMod p) × ℕ))
    (hsqf : ∀ g, g ≠ 0 → SquarefreeDecomp g (sqf g))
    (ddf : Polynomial (ZMod p) → List (Polynomial (ZMod p) × ℕ))
    (hddf : ∀ g, Monic g → Squarefree g → DDFCorrect g (ddf g))
    (edf : Polynomial (ZMod p) → ℕ → List (Polynomial (ZMod p)))
    (hedf : ∀ g d, Monic g → Squarefree g →
            (∀ q, Irreducible q → q ∣ g → q.natDeg = d) →
            EDFCorrect g d (edf g d))
    : -- 结论：组合后得到完整不可约分解
      ∃ (lc : ZMod p) (factors : List (Polynomial (ZMod p) × ℕ)),
        -- 乘积还原
        f = C lc * ∏ pr in factors, pr.1 ^ pr.2
        -- 每个因子不可约、首一
        ∧ (∀ pr ∈ factors, Irreducible pr.1 ∧ Monic pr.1 ∧ pr.2 ≥ 1) := by
  sorry -- Phase 1 目标：填充此证明
```

#### T1.3 Z[x] 顶层正确性骨架

**文件**：`CLPoly/Pipeline/FactorZZ.lean`

```lean
/-- Z[x] 因式分解的顶层正确性：
    假设 Zp[x] 因式分解、Hensel 提升、因子重组各自正确，
    则组合结果是完整不可约分解 -/
theorem factor_ZZ_correct
    (f : Polynomial ℤ) (hf : Squarefree f) (hprim : IsPrimitive f)
    -- 假设子过程
    (factor_zp_ok : ...)
    (hensel_ok : ...)
    (recombine_ok : ...)
    : ∃ factors : List (Polynomial ℤ),
        f = ∏ g in factors, g
        ∧ (∀ g ∈ factors, Irreducible g) := by
  sorry
```

#### T1.4 评估与接口锁定

- 检查所有 `Prop` 的类型签名是否自洽（能否组合成顶层定理）
- 检查子过程之间的输入/输出类型是否匹配
- 如果发现接口不匹配，**在此阶段修正**，代价极低

### 2.4 Phase 1 验收标准

- `Spec.lean` 中所有规约 Prop 定义完成，编译通过
- `FactorZp.lean` 顶层定理陈述完成（可有 `sorry`）
- `FactorZZ.lean` 顶层定理陈述完成（可有 `sorry`）
- 各子过程的接口规约**已锁定**，后续仅填充不修改签名

---

## 3. Phase 2：数学基石 L3（2-3 周）

### 3.1 目标

证明 DDF/EDF 正确性所需的**数学定理**，不涉及算法实现细节。

### 3.2 工作流

```
T2.1 X^{p^d} - X 的 Separable 性
  ↓
T2.2 Thm 2.1: 不可约 g, deg g | d → g ∣ (X^{p^d} - X)
  ↓
T2.3 Cor 2.2: gcd(X^{p^d} - X, f) 的不可约因子刻画
  ↓
T2.4 CRT 分解（EDF 数学基础）
  ↓
T2.5 分裂概率下界（EDF 概率分析）
```

### 3.3 任务详述

#### T2.1 X^{p^d} - X 的 Separable 性 (~20-50 行)

**文件**：`CLPoly/Math/FiniteFieldFact.lean`

```lean
theorem X_pow_sub_X_separable (d : ℕ) (hd : 0 < d) :
    Separable (X ^ (p ^ d) - X : Polynomial (ZMod p))
```

证明：导数 = -1（char p 下 p^d = 0），gcd(f, -1) = 1。

#### T2.2 Thm 2.1: 不可约因子整除 X^{p^d} - X (~100-150 行)

**文件**：`CLPoly/Math/FiniteFieldFact.lean`

```lean
/-- 不可约多项式 g, deg g = k, k | d → g | (X^{p^d} - X) -/
theorem irreducible_dvd_X_pow_sub_X
    (g : Polynomial (ZMod p)) (hg : Irreducible g) (hm : Monic g)
    (k d : ℕ) (hk : g.natDeg = k) (hdvd : k ∣ d) :
    g ∣ (X ^ (p ^ d) - X : Polynomial (ZMod p))
```

桥接路径（已在 E2 确认可行）：
1. `GaloisField.finrank` → finrank = k
2. `nonempty_algHom_iff_finrank_dvd` → k | d → F_{p^k} ↪ F_{p^d}
3. g 的根在 F_{p^k} 中 → 在 F_{p^d} 中 → 是 X^{p^d}-X 的根
4. 域上多项式由根决定 → g | (X^{p^d}-X)

#### T2.3 Cor 2.2: gcd 的不可约因子刻画 (~50 行)

```lean
/-- gcd(X^{p^d} - X, f) 恰为 f 中所有度整除 d 的不可约因子之积 -/
theorem gcd_X_pow_sub_X_factors
    (f : Polynomial (ZMod p)) (hm : Monic f) (hsq : Squarefree f) (d : ℕ) :
    ∀ q, Irreducible q → (q ∣ gcd (X ^ (p ^ d) - X) f ↔ q ∣ f ∧ q.natDeg ∣ d)
```

#### T2.4 EDF 的 CRT 分解 (~80 行)

等度多项式在 F_q[x] 上的中国剩余定理分解。

#### T2.5 分裂概率下界 (~80 行)

Cantor-Zassenhaus 的 gcd 非平凡概率 ≥ 1 - 2^{1-k}。

### 3.4 Phase 2 验收标准

- 上述 5 个数学定理全部通过，无 `sorry`
- 这些定理直接对接 Phase 1 中 `Spec.lean` 的规约

---

## 4. Phase 3：算法模型 L2 — Zp[x] 管线（3-4 周）

### 4.1 目标

逐个实现 SQF、DDF、EDF 算法模型，证明它们满足 Phase 1 定义的规约，**消除 `factor_Zp_correct` 中的所有 `sorry`**。

### 4.2 工作流

```
T3.1 DDF 算法建模 + 正确性     ← 最关键，概念验证
  ↓
T3.2 SQF 算法建模 + 正确性
  ↓
T3.3 EDF 算法建模 + 正确性     ← partial def，终止性文档论证
  ↓
T3.4 Zp 管线组装：消除 factor_Zp_correct 的 sorry
```

### 4.3 任务详述

#### T3.1 DDF 建模与正确性 (~300 行)

**文件**：`CLPoly/Algorithm/DDF.lean`

对应 C++ `__ddf_Zp` 的算法逻辑，操作 Mathlib 多项式类型，建模为顶层递归函数：

```lean
def ddfLoop (h f_star : Polynomial (ZMod p)) (d : ℕ)
    (acc : List (Polynomial (ZMod p) × ℕ))
    : List (Polynomial (ZMod p) × ℕ) :=
  if f_star.natDeg < 2 * d then
    if 0 < f_star.natDeg then acc ++ [(f_star, f_star.natDeg)] else acc
  else
    let h' := (h ^ p) %ₘ f_star           -- powmod
    let gd := gcd (h' - X) f_star          -- gcd 步
    if 0 < gd.natDeg then
      ddfLoop (h' %ₘ (f_star /ₘ gd)) (f_star /ₘ gd) (d + 1) (acc ++ [(gd, d)])
    else
      ddfLoop h' f_star (d + 1) acc
termination_by f_star.natDeg + 1 - 2 * d
```

证明义务：
- 循环不变量：h ≡ X^{p^d} (mod f_star)
- 分裂正确性：gd = gcd(X^{p^d}-X, f_star)（利用 T2.3）
- 提前终止：deg(f_star) < 2d → 至多一个不可约因子
- 顶层 `DDFCorrect f (ddfLoop X f 1 [])` 成立

#### T3.2 SQF 建模与正确性 (~200 行)

**文件**：`CLPoly/Algorithm/SquarefreeZp.lean`

对应 C++ `__squarefree_Zp` 的算法逻辑。递归结构（char p 下提取 p-th root）。

证明 `SquarefreeDecomp f (squarefree_Zp f)` 成立。

#### T3.3 EDF 建模与正确性 (~150 行)

**文件**：`CLPoly/Algorithm/EDF.lean`

对应 C++ `__edf_Zp` 的算法逻辑。用 `partial def`（概率终止）。

证明义务是**条件化的**：假设随机选取找到非平凡因子，则分裂正确。

#### T3.4 Zp 管线组装 (~50 行)

**文件**：`CLPoly/Pipeline/FactorZp.lean`

将 T3.1-T3.3 的具体实现代入 Phase 1 的骨架，消除 `sorry`：

```lean
-- 现在不再需要假设子过程正确，直接用已证明的版本
theorem factor_Zp_correct (f : Polynomial (ZMod p)) (hf : f ≠ 0) :
    ... := by
  -- 用 squarefree_Zp_correct, ddf_correct, edf_correct 组合
  ...
```

### 4.4 每函数 6 步工作流

每个 C++ 函数严格执行：

| 步骤 | 操作 | 产出 |
|------|------|------|
| Step 1: 提取 | C++ 源码添加标注 | `CLPoly/Annotated/` |
| Step 2: L2 建模 | 提取算法逻辑 → Lean 4 `def`（Mathlib 类型） | `CLPoly/Algorithm/` |
| Step 3: 陈述 | 写 `theorem ... := by sorry` | 同上 |
| Step 4: 证明 | 替换 `sorry`，利用 Phase 2 的 L3 定理 | 同上 |
| Step 5: 审查 | 对照 Lean 算法模型 ↔ C++ 算法逻辑 | `CLPoly/Review/` |

**算法建模原则**：Lean 函数捕获 C++ 的算法逻辑（循环结构、分支条件），但操作 Mathlib 数学类型，抽象掉实现细节。1:1 控制流对应在 Phase 5（L1 实现模型）中处理。

### 4.5 Phase 3 验收标准

- `factor_Zp_correct` 通过 `#check`，无 `sorry`（EDF 终止性除外）
- DDF 的 `.induct` 归纳原理可用于循环不变量证明
- `lake build` 全部通过

---

## 5. Phase 4：算法模型 L2 — Z[x] 管线（4-6 周）

### 5.1 目标

完成 Hensel 提升和因子重组，**消除 `factor_ZZ_correct` 中的所有 `sorry`**。

### 5.2 工作流

```
T4.1 L3: Hensel 唯一性 (Thm 2.1)                 T4.3 L3: Mignotte 界
  ↓                                                   ↓
T4.2 L3: 多因子 Hensel 提升 (Thm 2.2)                ↓
  ↓                                                   ↓
T4.4 L2: __hensel_lift 建模+证明                      ↓
  ↓                                                   ↓
T4.5 L2: 因子恢复 + 剪枝 + CLD (Thm 6.1-6.5)  ←────┘
  ↓
T4.6 L2: __zassenhaus_recombine 建模+证明
  ↓
T4.7 Z[x] 管线组装：消除 factor_ZZ_correct 的 sorry
```

### 5.3 关键任务

#### T4.1 Hensel 唯一性（最难的单个任务，~200-300 行）

```lean
theorem hensel_uniqueness
    (F : Polynomial (ZMod (p ^ k)))
    (A B : Polynomial (ZMod p))
    (hc : IsCoprime A B)
    (hF : ...) :
    ∃! (Am Bm : Polynomial (ZMod (p ^ k))), ...
```

**风险缓解**：若过于困难，先 `axiom` 占位，优先完成上层证明。

#### T4.2 多因子 Hensel 提升 (~300 行)

构造性地给出提升步骤 + 不变量 I1-I4。

#### T4.5 修复方案形式化（核心价值）

对应 `formal-proof-univar-factorization.md` §6 的 Thm 6.1-6.5。这是 CLPoly bug 修复的形式化保证。

### 5.4 Phase 4 验收标准

- `factor_ZZ_correct` 通过 `#check`
- Hensel 引理：完成或以 `axiom` 占位
- Thm 6.1-6.5（修复方案）全部无 `sorry`

---

## 6. Phase 5：L1 实现模型 + LLL

### 6.1 L1 实现模型

1:1 对应 C++ 实现：U64/I64 补码语义、Vec 越界安全、Move 语义。控制流必须与 C++ 一一对应，操作对象是 C++ 语义模型类型。

**依赖 L2 完成**：L1 精化证明证明 L1 行为与 L2 算法模型一致，因此排在 L2 之后。

#### B3 Bug 形式化证明

```lean
/-- 旧代码在 p > 2^63 时 (int64_t)(p-1) 溢出导致错误 -/
theorem old_neg_one_wrong (p : U64) (hp : p.val > 2 ^ 63) :
    old_neg_one p ≠ ⟨p.val - 1, ...⟩

/-- 新代码 (p-1) % p 对任意 p ≥ 2 正确 -/
theorem new_neg_one_correct (p : U64) (hp : 2 ≤ p.val) :
    new_neg_one p = ⟨p.val - 1, ...⟩
```

### 6.2 LLL 验证（高难度，视经验决定）

- L3: LLL 格基约化正确性 (~500 行)
- L2: van Hoeij 重组因子绑定

---

## 7. 文件结构

```
lean/
├── lakefile.lean
├── CLPoly.lean
└── CLPoly/
    ├── Spec.lean                        -- 子过程接口规约（Phase 1 核心）
    ├── Pipeline/                        -- 顶层骨架（Phase 1）
    │   ├── FactorZp.lean               -- Zp[x] 顶层正确性
    │   └── FactorZZ.lean               -- Z[x] 顶层正确性
    ├── Math/                            -- L3 数学基石（Phase 2）
    │   ├── FiniteFieldFact.lean        -- Thm 2.1, Cor 2.2, Separable
    │   ├── HenselLemma.lean            -- Hensel 唯一性 + 提升
    │   ├── MignotteBound.lean          -- Mignotte 界
    │   └── CantorZassenhaus.lean       -- EDF 概率分析
    ├── Algorithm/                       -- L2 算法模型（Phase 3-4）
    │   ├── DDF.lean                    -- ddf_Zp 建模 + 证明
    │   ├── SquarefreeZp.lean           -- squarefree_Zp
    │   ├── EDF.lean                    -- edf_Zp
    │   ├── HenselLift.lean             -- __hensel_lift
    │   ├── Recombine.lean              -- __zassenhaus_recombine
    │   └── LLLFactorize.lean           -- __lll_factorize
    ├── Impl/                            -- L1 实现模型（Phase 5）
    │   ├── UInt64.lean
    │   ├── Vector.lean
    │   └── Ownership.lean
    ├── Experiment/                      -- Phase 0 实验（保留参考）
    │   ├── E1_ZpPolyAPI.lean
    │   ├── E2_TheoremBridge.lean
    │   ├── E3_ZModPkDiv.lean
    │   └── E4_Termination.lean
    └── Annotated/                       -- Step 1 标注版 C++
```

---

## 8. 里程碑

| 里程碑 | 内容 | 验收标准 |
|--------|------|---------|
| **M0** Phase 0 | 环境 + API 验证 | ✅ 已完成 |
| **M1** Phase 1 | 顶层骨架 + 接口锁定 | 定理陈述编译通过，接口自洽 |
| **M2** Phase 2 | L3 数学基石 | 5 个数学定理无 `sorry` |
| **M3** Phase 3 | L2 Zp[x] 管线 | `factor_Zp_correct` 无 `sorry` |
| **M4** Phase 4 | L2 Z[x] 管线 | `factor_ZZ_correct` 无 `sorry`（Hensel 可 `axiom`） |
| **M5** Phase 5 | L1 实现模型 + 补全 | L1 精化证明完成，消除所有 `sorry` 和 `axiom` |

---

## 9. 风险与缓解

| 风险 | 级别 | 缓解策略 |
|------|------|---------|
| 接口规约无法组合成顶层定理 | 🟡 中 | **Phase 1 就会暴露**，修正代价极低 |
| Thm 2.1 桥接困难 | 🟢 低 | Phase 0 已确认所有 API 可用 |
| Hensel 提升工作量超预期 | 🟡 中 | `axiom` 占位，优先上层 |
| EDF 概率终止性 | 🟢 低 | `partial def` + 文档论证 |
| 编译时间 | 🟢 低 | 增量编译 ~1s，预编译缓存已下载 |

---

## 10. Mathlib 可用性速查

Phase 0 实验已确认的 API（详见 `mathlib-gap-report.md`）：

| 类别 | 关键 API | 状态 |
|------|---------|------|
| 基础结构 | `Field (ZMod p)`, `EuclideanDomain`, `GCDMonoid`, `DecidableEq` | ✅ |
| 多项式运算 | `divByMonic`, `modByMonic`, `modByMonic_add_div`, `Monic.map` | ✅ |
| 有限域 | `roots_X_pow_card_sub_X`, `splits_X_pow_card_sub_X`, `nonempty_algHom_iff_finrank_dvd` | ✅ |
| 因子谓词 | `Irreducible`, `Squarefree`, `Separable` | ✅ |
| ZMod p^k | `castHom`, `valMinAbs`, `divByMonic`(非整环) | ✅ |
| 终止性 | `termination_by` + `omega` 自动, `.induct` 自动生成 | ✅ |
