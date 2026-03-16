# 框架组合证明：填补 Phase 1 设计缺陷

> 日期：2026-03-16
> 分支：`feature/formal-proofs`

---

## 1. 问题发现

Phase 1 原始设计的验收标准是：

> - Spec.lean 中所有规约 Prop 定义完成，编译通过
> - FactorZp.lean 顶层定理陈述完成（可有 sorry）
> - FactorZZ.lean 顶层定理陈述完成（可有 sorry）
> - 各子过程的接口规约已锁定

T1.4 接口审查只做了**类型级一致性检查**（输出类型匹配下游输入类型）。这是必要的，但不充分。

**遗漏**：没有验证规约在逻辑层面能否组合。也就是说——即使假设每个子过程都满足规约，"组合起来能达成顶层目标"这一命题本身没有证明。

这个组合证明被原计划归入 Phase 3（与算法实现混在一起），但它实际上：
- **不依赖算法实现**（纯粹是规约间的逻辑推导）
- **不依赖 Phase 2 的数学定理**（只需标准代数事实）
- **如果失败，说明接口设计有根本缺陷**，Phase 2-4 的工作全部白费

这是一个应该在 Phase 1 完成的验证步骤，却被推迟到了 Phase 3。

## 2. 修正：框架组合证明

在 `FactorZp.lean` 中，将 `factor_Zp_correct` 的单个 `sorry` 替换为结构化证明。

### 证明架构

分两层：

**第一层：`ddf_edf_combine`**（完全证明，无 sorry）

对单个首一无平方多项式 g，证明 DDF + EDF 组合给出完整不可约分解：

```
DDFCorrect g (ddf g)
  → 对每个 (gd, d)：gd 首一（DDFCorrect.5）
  → gd 无平方（Squarefree.squarefree_of_dvd，因 gd ∣ g 且 g 无平方）
  → EDFCorrect gd d (edf gd d)
  → gd = (edf gd d).prod（首一 Associated → 相等）
  → g = flatMap 所有不可约因子的乘积
```

关键步骤：首一 Associated 消除了 DDF/EDF 引入的单位差异，使 Associated 退化为相等。

**第二层：主定理**

```
SQF: f ≈ ∏ gᵢ^eᵢ
  → 对每个 gᵢ：由第一层得 gᵢ = ∏ q_k（不可约因子）
  → f ≈ ∏∏ q_k^eᵢ
  → f = C(lc) * ∏∏ q_k^eᵢ（提取单位为非零常数）
```

### sorry 分析

| sorry | 内容 | 类别 |
|-------|------|------|
| `poly_unit_eq_C` | 域上多项式环的单位是非零常数 | 标准代数 |
| `eq_of_associated_monic` | 首一 Associated → 相等 | 标准代数 |
| `monic_list_prod` | 首一列表的积是首一 | 归纳法 |
| `list_prod_pow` | (∏aᵢ)^n = ∏(aᵢ^n) | 交换幺半群 |
| `list_prod_flatMap` | flatMap 的积 = map-prod 的积 | 列表归纳 |
| 主定理 `key` | factors 积 = SQF 积（组合以上引理） | 列表操作 |
| 主定理 Step 6 | 单位代数推导 f = C(lc) * prod | 代数操作 |

已完全证明（无 sorry）的部分：
- `ddf_edf_combine` 全部逻辑（DDF+EDF 组合正确性）
- 主定理第二条：每个因子 Irreducible + Monic + 重数 ≥ 1

**结论：所有 sorry 都是标准数学事实和列表操作，不涉及接口设计问题。框架可组合性已验证。**

## 3. 启示：Phase 1 验收标准的修正

原始标准只有"编译通过 + 类型匹配"。应增加：

> **框架组合证明**：顶层定理的 sorry 替换为结构化证明，剩余 sorry 仅限于标准数学引理和列表操作，不含接口设计相关的 sorry。

这等价于要求：**接口锁定 = 类型一致性 + 逻辑可组合性**。

类型一致性（T1.4）回答的是"能拼在一起吗？"，
组合证明回答的是"拼在一起能得到我们想要的吗？"

后者更重要——T1.4 发现的 DDFCorrect 缺少 `Monic` 正是在为组合证明铺路。没有 `Monic`，`eq_of_associated_monic` 无法应用，DDF→EDF 的 Associated 无法退化为相等，整个乘积链断裂。

## 4. 技术笔记

- Lean 4 中 `List.bind` 已不存在，使用 `List.flatMap`
- `List.map_congr` → `List.map_congr_left`
- `List.mem_bind` → `List.mem_flatMap`
- `let` 绑定在证明 goal 中会被 `rw` 穿透——避免在引理声明中使用 `let`，或用 `exact ... .trans ...` 代替 `rw`
- `Squarefree.squarefree_of_dvd : x ∣ y → Squarefree y → Squarefree x` 已在 Mathlib 中

## 5. 涉及文件

| 文件 | 操作 |
|------|------|
| `proof/lean/CLPoly/Pipeline/FactorZp.lean` | 重写：5 个辅助引理 + ddf_edf_combine + 主定理结构化证明 |
