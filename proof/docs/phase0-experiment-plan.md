# Phase 0：环境搭建与关键假设验证

## 目标

在投入 Phase 1（DDF 端到端验证）之前，通过 5 个小规模实验确认关键技术假设。每个实验独立、可在数小时内完成，产出是一个可编译的 `.lean` 文件 + 结论。

**Phase 0 结束时的决策**：确认或修正 `implementation-roadmap.md` 中的技术方案，特别是：
- Thm 2.1 的证明路径是否可行
- `ZMod (p^k)` 上的多项式除法是否如预期工作
- DDF 循环的终止性证明是否顺畅
- 编译时间是否可接受

---

## E0：项目初始化

### 任务

创建 `lean/` 目录，配置 Lean 4 + Mathlib 依赖，确认 `lake build` 通过。

### 步骤

1. 创建 `lean/lakefile.lean`：

```lean
import Lake
open Lake DSL

package CLPoly where
  leanOptions := #[⟨`autoImplicit, false⟩]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "stable"

@[default_target]
lean_lib CLPoly where
  srcDir := "."
```

2. 创建 `lean/lean-toolchain`：对齐 Mathlib 的 Lean 版本
3. 创建 `lean/CLPoly.lean`（空主入口）
4. 运行 `lake update && lake build`

### 验收

`lake build` 成功，无错误。记录首次编译时间。

---

## E1：Mathlib 多项式 API 可用性

### 假设

> `Polynomial (ZMod p)` 在 `p` 为素数时拥有 `EuclideanDomain`、`GCDMonoid`、`DecidableEq` 实例，`gcd`、`divByMonic`、`modByMonic` 均可正常使用。

### 实验

`lean/CLPoly/Experiment/E1_ZpPolyAPI.lean`：

```lean
import Mathlib.RingTheory.Polynomial.Basic
import Mathlib.Algebra.Polynomial.FieldDivision
import Mathlib.Data.ZMod.Basic
import Mathlib.FieldTheory.Finite.Basic

open Polynomial

-- 1. ZMod p 是域
example : Field (ZMod 7) := inferInstance

-- 2. Polynomial (ZMod 7) 是 EuclideanDomain
example : EuclideanDomain (Polynomial (ZMod 7)) := inferInstance

-- 3. GCDMonoid
example : GCDMonoid (Polynomial (ZMod 7)) := inferInstance

-- 4. DecidableEq
example : DecidableEq (Polynomial (ZMod 7)) := inferInstance

-- 5. gcd 可用
#check (gcd : Polynomial (ZMod 7) → Polynomial (ZMod 7) → Polynomial (ZMod 7))

-- 6. 构造具体多项式并计算
noncomputable def f7 : Polynomial (ZMod 7) := X ^ 3 + C 2 * X + C 1
noncomputable def g7 : Polynomial (ZMod 7) := X ^ 2 + C 3

-- 7. divByMonic / modByMonic
noncomputable def q7 := f7 /ₘ (X ^ 2 + C 3 : Polynomial (ZMod 7))
noncomputable def r7 := f7 %ₘ (X ^ 2 + C 3 : Polynomial (ZMod 7))

-- 8. 除法恒等式
example : (X ^ 2 + C 3 : Polynomial (ZMod 7)) * q7 + r7 = f7 := by
  simp [q7, r7, modByMonic_add_div]

-- 9. Squarefree 与 Separable
#check @Polynomial.Separable.squarefree (ZMod 7) _
```

### 关注点

- `noncomputable` 是否是必须的？如果所有定义都需要 `noncomputable`，是否影响后续证明？
- `inferInstance` 的求解速度——是否有超时风险？
- `gcd` 的实际行为：是否返回首一多项式？

### 产出

确认/否认假设，记录遇到的 API 陷阱。

---

## E2：定理 2.1 桥接路径验证

### 假设

> 可以从 Mathlib 的 `roots_X_pow_card_sub_X`（根的刻画）桥接到不可约因子乘积形式，工作量 ~100-150 行。

### 实验

`lean/CLPoly/Experiment/E2_TheoremBridge.lean`：

不要求完成完整证明，只要求走通前 2-3 步，评估可行性。

```lean
import Mathlib.FieldTheory.Finite.Basic
import Mathlib.FieldTheory.Finite.GaloisField
import Mathlib.Algebra.Squarefree.UniqueFactorizationDomain

open Polynomial FiniteField

variable (p : ℕ) [hp : Fact (Nat.Prime p)]

-- 步骤 1：确认 roots 定理的精确签名
#check @roots_X_pow_card_sub_X

-- 步骤 2：X^{p^d} - X 无平方
-- 导数 = p^d · X^{p^d-1} - 1 = -1（特征 p 下 p^d = 0）
-- 所以 gcd(X^{p^d}-X, -1) = 1，即 Separable
example (d : ℕ) (hd : 0 < d) :
    Separable (X ^ (p ^ d) - X : Polynomial (ZMod p)) := by
  sorry -- 尝试证明

-- 步骤 3：UFD 中无平方 → 不可约因子无重复
-- 在 Mathlib 中找到：Squarefree 与 UniqueFactorizationMonoid 的关系
#check UniqueFactorizationMonoid
#check Squarefree
-- 寻找：squarefree_iff_nodup_normalizedFactors 或类似

-- 步骤 4：不可约 g, deg g = k | d → g 的根在 F_{p^d} 中
-- 需要：GaloisField 嵌入 + 根属于分裂域
#check GaloisField
#check GaloisField.instFintype  -- card = p^n
-- 寻找 F_{p^k} ↪ F_{p^d} 当 k | d
#check @FiniteField.algHomEquivOfCardEq  -- 或类似
```

### 关注点

- 步骤 2（Separable 证明）能否自动化？`simp` + `ring` 是否足够？
- 步骤 3 在 Mathlib 中的具体 API 是什么？`squarefree_iff_nodup_normalizedFactors`？
- 步骤 4 域嵌入的确切定理名称
- 是否需要处理 `noncomputable` / `Decidable` 问题

### 产出

- 确认桥接路径的每一步在 Mathlib 中有对应工具
- 识别缺失的胶水引理
- 修正工作量预估

---

## E3：`ZMod (p^k)` 上的多项式除法

### 假设

> `Polynomial.divByMonic` 在 `ZMod (p^k)`（非整环）上对首一除数工作正常。`castHom` 可用于层间投影。

### 实验

`lean/CLPoly/Experiment/E3_ZModPkDiv.lean`：

```lean
import Mathlib.Algebra.Polynomial.Div
import Mathlib.Data.ZMod.Basic

open Polynomial

-- 1. ZMod (p^k) 是 CommRing（不是 IsDomain）
example : CommRing (ZMod (7 ^ 3)) := inferInstance
-- 确认不是 IsDomain
-- example : IsDomain (ZMod (7 ^ 3)) := inferInstance  -- 应该失败

-- 2. 首一多项式除法
noncomputable def f343 : Polynomial (ZMod 343) := X ^ 3 + C 2 * X + C 1
noncomputable def g343 : Polynomial (ZMod 343) := X ^ 2 + C 3  -- 首一

-- divByMonic 应该工作
noncomputable def q343 := f343 /ₘ g343
noncomputable def r343 := f343 %ₘ g343

-- 3. 恒等式
example : g343 * q343 + r343 = f343 := by
  simp [q343, r343, g343, f343, modByMonic_add_div]

-- 4. castHom 层间投影
#check @ZMod.castHom
-- ZMod 343 → ZMod 7（343 = 7^3，7 | 343）
noncomputable def proj : ZMod 343 →+* ZMod 7 :=
  ZMod.castHom (show 7 ∣ 343 by norm_num) (ZMod 7)

-- 5. Polynomial.map 配合 castHom
noncomputable def f7_from_343 := Polynomial.map proj f343
-- 验证：map 后与直接构造一致？

-- 6. 单位与模逆元
-- 3 和 343 互素，所以 3 是 ZMod 343 的单位
example : IsUnit (3 : ZMod 343) := by
  sorry -- 尝试 ZMod.unitOfCoprime 或 decide

-- 7. sym_mod 自定义
def symMod (n : ℕ) (a : ZMod n) : ℤ :=
  let v := (ZMod.val a : ℤ)
  if v ≤ (n : ℤ) / 2 then v else v - n

-- 测试 symMod 的行为
```

### 关注点

- `divByMonic` 在非整环上的**实际行为**是否与理论一致
- `castHom` 需要的 `dvd` 证明如何提供（`norm_num`？`decide`？）
- `IsUnit` 的构造是否方便
- `Polynomial.map` 的系数映射是否如预期
- `ZMod.val` 在大模数下的行为

### 产出

确认 Hensel 提升所需的 `ZMod (p^k)` 基础设施是否可用，识别需要自建的部分。

---

## E4：递归函数终止性与函数归纳

### 假设

> DDF 循环可以用 `let rec loop ... termination_by f_star.natDeg + 1 - 2 * d` 建模，`omega` 可自动证明递减性。Lean 4 生成的 `.induct` 归纳原理可用于证明循环不变量。

### 实验

`lean/CLPoly/Experiment/E4_Termination.lean`：

```lean
-- Toy DDF：模拟 DDF 的递归结构，不涉及多项式
-- 用自然数模拟 deg(f_star)

/-- 模拟 DDF 循环：
    - n 模拟 deg(f_star)
    - d 从 1 开始递增
    - 每步要么 n 减小（提取了非平凡 gcd），要么 d 增大
    - 当 n < 2*d 时终止 -/
def toyDDF (n : ℕ) : List (ℕ × ℕ) :=
  let rec loop (n d : ℕ) (acc : List (ℕ × ℕ)) : List (ℕ × ℕ) :=
    if n < 2 * d then
      if 0 < n then acc ++ [(n, n)] else acc
    else
      -- 模拟：50% 概率提取一个度为 d 的因子
      if n % 3 = 0 then  -- 模拟非平凡 gcd
        let n' := n - d   -- f_star 度数减小
        loop n' (d + 1) (acc ++ [(d, d)])
      else
        loop n (d + 1) acc
  termination_by n + 1 - 2 * d
  -- decreasing_by all_goals simp_wf; omega  -- 是否需要？
  loop n 1 []

-- 测试
#eval toyDDF 10  -- 应该返回某个列表

-- 函数归纳原理
#check toyDDF.loop.induct  -- 自动生成？

-- 用归纳原理证明一个简单性质
theorem toyDDF_loop_sum (n d : ℕ) (acc : List (ℕ × ℕ)) :
    -- 输出列表的第二分量之和 ≤ n
    True := by  -- 先用 True 占位，看 .induct 的签名
  trivial
```

### 关注点

- `termination_by n + 1 - 2 * d` 是否被接受（自然数减法的 well-founded）
- 是否需要 `decreasing_by`，还是 `omega` 自动搞定
- `.induct` 归纳原理的精确签名——有多少个分支？参数是什么？
- `#eval` 是否能工作（computability）

### 进阶实验

如果基础实验通过，用真正的 `Polynomial (ZMod p)` 替换自然数：

```lean
-- 用真正的多项式类型测试
noncomputable def toyDDF_poly (f : Polynomial (ZMod 7)) : List (Polynomial (ZMod 7) × ℕ) :=
  let rec loop (h f_star : Polynomial (ZMod 7)) (d : ℕ)
      (acc : List (Polynomial (ZMod 7) × ℕ)) :=
    if f_star.natDeg < 2 * d then acc
    else
      let h' := (h ^ 7) %ₘ f_star  -- 模拟 powmod
      let gd := gcd (h' - X) f_star
      if 0 < gd.natDeg then
        loop (h' %ₘ (f_star /ₘ gd)) (f_star /ₘ gd) (d + 1) (acc ++ [(gd, d)])
      else
        loop h' f_star (d + 1) acc
  termination_by f_star.natDeg + 1 - 2 * d
  loop X f 1 []
```

- `noncomputable` 是否阻止 `#eval`？
- `natDeg` 在 `termination_by` 中是否工作？
- `gcd` 返回的多项式是否首一？（影响 `/ₘ` 的行为）

### 产出

确认 DDF 建模的技术可行性。如果 `termination_by` + `omega` 不够，确定替代方案（手动 `decreasing_by`、`fuel`、`WellFounded.fix`）。

---

## E5：编译时间与开发体验

### 实验

1. 测量 `lake build` 首次编译时间（含 Mathlib 下载 + 编译/缓存获取）
2. 测量增量编译时间（修改一个 `.lean` 文件后重新 build）
3. 测量单文件 `import Mathlib` 的处理时间
4. 测试 VS Code + lean4 插件的响应速度

### 关注点

- Mathlib 是否提供预编译 `.olean` 缓存（`lake exe cache get`）
- 增量编译是否可接受（< 30 秒为优，< 2 分钟可接受）
- `sorry` 占位是否能加速编译（跳过证明 elaboration）

### 产出

确定开发节奏：是否需要拆分模块、是否需要 CI 离线构建。

---

## 执行顺序

```
E0 (项目初始化)
 ↓
E5 (编译时间) ← 立即可测
 ↓
E1 (ZpPoly API) + E4 (终止性)  ← 并行
 ↓
E2 (Thm 2.1 桥接) + E3 (ZMod p^k 除法)  ← 并行，依赖 E1
 ↓
撰写 mathlib-gap-report.md（Phase 0 交付物）
```

预估总用时：**3-5 天**。

---

## 决策矩阵

| 实验 | 如果成功 | 如果失败 |
|------|---------|---------|
| E1 | 按路线图继续 | 调查具体缺失的实例，可能需要手动构造 |
| E2 | Thm 2.1 桥接预估 ~150 行不变 | 若关键步骤缺失，考虑从零证明 Thm 2.1（~300 行），或先 `axiom` 占位 |
| E3 | Phase 3 Hensel 方案用 `ZMod (p^k)` 类型级 | 改为方案 B（ℤ 上 + 显式模约化），参考 Isabelle 做法 |
| E4 | DDF 建模方案确认 | 若 `termination_by` 失败 → 用 `fuel` 参数；若 `.induct` 不好用 → 手动归纳 |
| E5 | 正常开发 | 若首次编译 > 30 分钟 → 配置 olean 缓存；若增量 > 5 分钟 → 拆分模块减少依赖 |
