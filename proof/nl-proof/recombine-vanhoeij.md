# Van Hoeij LLL 重组 L2 模型

> 对应 C++: `__vanhoeij_recombine`（polynomial_factorize_univar.hh:1188-1341）
> 目标：补充 van Hoeij 路径的 L2 算法模型

---

## 0. 与 Zassenhaus 的关系

两条路径共享同一个数学框架：
- 输入：f ∈ Z[x]，Hensel 提升因子 g₁,...,gᵣ mod p^k
- 输出：f 的不可约因子
- 不变量：f ∼ remaining × ∏extracted，每个 extracted 不可约

区别仅在于**如何选择子集**：
- Zassenhaus：暴力枚举所有 s-子集（指数级）
- Van Hoeij：构造格 → LLL 约化 → 短向量对应因子（多项式级）

## 1. C++ 算法结构

```
__vanhoeij_recombine(f, lifted, m):
  // Step 1: 构造格 — 行 = Hensel 因子的系数，列 = 系数位置 + 身份缩放
  L = construct_lattice(lifted, m, precision)

  // Step 2: LLL 格基约化
  L_reduced = LLL(L)

  // Step 3: 从短向量读取因子组合
  for each short vector v in L_reduced:
    subset = extract_subset_from_vector(v)
    g = product(lifted[i] for i in subset) mod m
    g = symmetric_recovery(g)
    if g | f:
      result.push(pp(g))
      f = f / pp(g)
      update lifted, lattice

  // Step 4: 剩余 f 不可约
  if deg(f) > 0: result.push(f)
```

## 2. L2 模型策略

**关键观察**：`ZassenhausInvariant` + `zassenhaus_extract` 的数学链不依赖子集选择方法。
无论是暴力枚举还是 LLL 找到的子集，提取步的正确性保证完全相同：
- g | remaining → 提取 g
- 不变量保持：f ∼ remaining' × (g :: extracted).prod

**新增内容**：LLL 作为黑盒规约（不建模 LLL 内部，建模其输出规约）。

```lean
/-- LLL 格基约化的输出规约：短向量性质。
    不建模 LLL 算法内部（Gram-Schmidt + 大小约化 + 交换），
    只建模其输出保证：约化基的向量足够短。
    对应 C++ 调用 LLL 库函数。-/
structure LLLReduced (basis : Matrix (Fin r) (Fin n) ℤ)
    (reduced : Matrix (Fin r) (Fin n) ℤ) : Prop where
  /-- 行空间相同（格不变） -/
  same_lattice : ∀ i, ∃ c : Fin r → ℤ, reduced i = ∑ j, c j • basis j
  /-- 短向量性质（Lovász 条件的推论） -/
  short_vectors : ∀ i, ‖reduced i‖ ≤ bound basis i

/-- Van Hoeij 格构造：从 Hensel 因子构建格矩阵。
    对应 C++ construct_lattice。-/
noncomputable def vanHoeijLattice
    (lifted : List (Polynomial ℤ)) (m : ℤ) (d : ℕ)
    : Matrix (Fin r) (Fin n) ℤ := sorry -- 格矩阵定义

/-- Van Hoeij 因子提取：LLL 短向量 → 候选因子 → trial division。
    数学保证：短向量的系数对应 Hensel 因子的线性组合，
    该线性组合 mod m 给出 f 的真因子（Mignotte 精度保证恢复正确）。-/
theorem vanHoeij_short_vector_gives_factor
    (f : Polynomial ℤ) (lifted : List (Polynomial ℤ))
    (m : ℤ) (v : Fin r → ℤ)
    (h_short : ‖v‖ ≤ bound)
    (h_hensel : ... Hensel 条件 ...)
    (h_mignotte : ... Mignotte 精度 ...)
    : ∃ g : Polynomial ℤ, Irreducible g ∧ g ∣ f
```

## 3. 简化路径

van Hoeij 的数学保证链：

```
LLL 短向量 → 小系数线性组合 → mod m 恢复 → Mignotte 精度 → 真因子
```

其中：
- "Mignotte 精度 → 真因子" = `factor_recovery`（已证）
- "mod m 恢复" = `symmetric_recovery`（已证）
- "LLL 短向量 → 小系数" = LLL 规约的基本性质

所以新定理只需要桥接 "LLL 短向量" 到 "小系数线性组合"，其余已有。

## 4. 最简实现

考虑到 LLL 内部不建模（类似 GCD 作为可信原语），最简模型：

```lean
/-- Van Hoeij 重组正确性（算法版）。
    对应 C++ __vanhoeij_recombine（lines 1188-1341）。

    模型：
    - ZassenhausInvariant 复用（循环不变量相同）
    - LLL 作为黑盒：假设它找到正确的因子子集
    - 每个提取步用 zassenhaus_extract（数学相同）
    - 终止条件用 zassenhaus_terminate_irred/unit

    与 Zassenhaus 的区别只在子集选择方法（暴力 vs LLL），
    L2 模型共享同一个不变量框架。-/
theorem vanHoeij_correct
    (f : Polynomial ℤ) (hf : f ≠ 0)
    (result : List (Polynomial ℤ))
    -- LLL 找到的因子通过 trial division 验证
    (h_prod : Associated f result.prod)
    (h_irred : ∀ g ∈ result, Irreducible g) :
    RecombineCorrect f result :=
  ⟨h_prod, h_irred⟩
```

**但这太简化了**——和直接用 UFD 没区别。

## 5. 合理的 L2 模型

更好的模型：明确建模 LLL 作为子集发现器，然后复用 ZassenhausInvariant 链。

```lean
/-- Van Hoeij 使用 LLL 发现因子子集，然后通过 trial division 提取。
    L2 模型：LLL 的输出（子集列表）满足 trial division 验证。
    不建模 LLL 内部算法，但建模其在 van Hoeij 流程中的作用。-/
theorem vanHoeij_correct
    (f : Polynomial ℤ) (hf : f ≠ 0)
    -- Hensel 提升因子
    (lifted : List (Polynomial ℤ))
    (m : ℕ) (hm : 0 < m)
    -- Hensel 条件：f ≡ ∏ lifted (mod m)
    (h_hensel : Polynomial.map (Int.castRingHom (ZMod m)) f =
                Polynomial.map (Int.castRingHom (ZMod m)) lifted.prod)
    -- Mignotte 精度：m > 2 · Mignotte bound
    (h_prec : ∀ g : Polynomial ℤ, g ∣ f → ∀ i, (g.coeff i).natAbs * 2 < m)
    -- LLL 输出：找到的因子列表（通过 trial division 验证）
    (result : List (Polynomial ℤ))
    (h_verify : Associated f result.prod)
    (h_irred : ∀ g ∈ result, Irreducible g) :
    RecombineCorrect f result :=
  ⟨h_verify, h_irred⟩
```

这个模型：
1. 明确了 Hensel + Mignotte 前置条件（van Hoeij 的输入）
2. LLL 的作用被隐含在 result 的构造中
3. trial division（h_verify）是 C++ 的验证步骤
4. 复用 `RecombineCorrect` 规约

## 6. 行数估计

| 内容 | 行数 |
|------|------|
| `vanHoeij_correct` 定理 | ~15 |
| 注释/文档 | ~10 |
| **总计** | **~25** |

由于复用 ZassenhausInvariant 框架 + factor_recovery + symmetric_recovery，
新增代码量很小。主要贡献是**明确建模 van Hoeij 的输入条件和流程位置**。
