# P2 多变量因式分解优化：架构文档

> 阶段：**待调研**
> 目标：将 CLPoly 多变量 Z[x₁,...,xₙ] 因式分解性能对齐 Maple 2019（MTSHL）
> 目标文件：`clpoly/polynomial_factorize_wang.hh`
> 对应调研：`docs/research/vanhoeij-lll-research.md` §5.2（初步参考材料）
> 专项调研：`docs/research/mtshl-research.md`（待创建）

---

## P2 优化范围

| 组件 | 目标函数 | 当前实现 | 目标实现 | 状态 |
|------|---------|---------|---------|------|
| **P2a** 稀疏 Diophantine 求解 | `__multivar_diophantine` | 稠密递归 + Bézout 链 | Zippel 稀疏插值嵌入（MTSHL） | **待调研** |
| **P2b** Wang 控制流调整 | `__wang_core` / `__multivar_hensel_lift` | 稠密控制流 | 稀疏探针求值点管理 | **待调研，依赖 P2a** |

> P2 与 P1 完全独立（不同文件、不同 Hensel 方向），可在 P1 稳定后按序进行。

---

## 1. 当前多变量因式分解管线

```
__factor_multivar(f ∈ Z[x₁,...,xₙ])
  → squarefreefactorize(f)
  → __wang_core(g)                           ← 每个无平方因子
        ├── __select_eval_point              选求值点 α
        ├── factorize(g(x₁, α))             单变量分解（P1 目标）
        ├── __wang_leading_coeff            LC 分配
        └── __multivar_hensel_lift          多变量 Hensel 提升
                └── __hensel_lift_one_var   逐变量 × 逐阶 Taylor 循环
                                └── __multivar_diophantine  ← P2a 目标（每阶调用一次，line 475）
```

**关键接口**：
```cpp
// __multivar_diophantine（polynomial_factorize_wang.hh:475）
// 求解：h = Σ δᵢ · ∏_{j≠i} G[j]，返回 {δ₁, ..., δᵣ}
std::vector<Poly>
__multivar_diophantine(h, G[], bezout_s[], bezout_denom, v_factors,
                       main_var, eval_vars, pa, depth)
```

---

## 2. 当前瓶颈分析

### 2.1 稠密性代价

`__multivar_diophantine` 采用**递归稠密求解**：

- **递归结构**：对 n 个变量递归 n 层，基本情形用 Bézout 链求解单变量 Diophantine
- **误差修正**（lines 560-630）：每层计算 Ĝᵢ = ∏_{j≠i} G[j]，是全稠密多项式乘法
- **调用频次**：`__hensel_lift_one_var` 对每个提升变量 xₖ、每阶 j = 1..deg(f, xₖ) 调用一次
- **总调用次数**：O(Σ deg(f, xₖ)) 次，每次处理 k 个变量的稠密多项式

设多项式有 n 个变量、各变量度数为 d：

| 量 | 当前代价（稠密，估算） | MTSHL 目标（稀疏，估算） |
|----|-------------------|-----------------------|
| 每次 Diophantine 调用（k 变量） | O(r · d^k) 次多项式运算 | O(T · k · d)（T = 解的非零项数） |
| 总 Hensel 提升（n 变量，逐变量累计） | O(r · d^n)（末项主导） | O(r · n · d · T) |

> **估算说明**：以上为单精度多项式运算次数的粗估，各变量度数均简化为 d。精确复杂度分析见 Monagan-Tuncer (2020) §3。对稀疏多项式（T << d^n）差距可达 10–100x。

### 2.2 稠密瓶颈的代码位置

```
__multivar_diophantine（line 475）
    递归情形（line 547-）:
      1. h_bar = assign(h, {xk: αk})       ← 求值（O(T_h) 项）
      2. 递归调用                           ← O(r · d^{k-1}) 代价
      3. 误差 e = h - Σ δᵢ·Ĝᵢ             ← 稠密乘法（瓶颈）
         Ĝᵢ = ∏_{j≠i} G[j]              ← 对每个 i 一次全乘法 O(d^k)
      4. Taylor 展开误差，逐阶求解         ← 递归调用

__hensel_lift_one_var（line 698）
    for j = 1..dk:                          ← 对每个变量的每个度调用
        ej = __taylor_coeff(e, xk, αk, j)  ← O(T_e) 项
        deltas = __multivar_diophantine(ej, ...)  ← O(r·d^{n-1})
```

---

## 3. 目标算法：MTSHL

### 3.1 核心思想（Monagan-Tuncer 2016）

**MTSHL（Monagan-Tuncer Sparse Hensel Lifting）**：在求解 Diophantine 方程时，不计算全稠密解，而是用 **Zippel 稀疏插值**逐项确定解的非零单项式。

```
稠密路径（当前）：
  解 h = Σ δᵢ · Ĝᵢ  → 对全稠密多项式 h 做 r 次 pseudo-division → 得稠密 δᵢ

稀疏路径（MTSHL）：
  1. 在随机点 β 处探针求值 h(β), Ĝᵢ(β) → 得到 δᵢ(β) 的候选值（线性方程组）
  2. 再取若干 β 点 → Zippel 稀疏插值重建 δᵢ 的非零项
  3. 只对实际出现的单项式做精确验证
  → 代价 = O(T_δᵢ · n · d)，T_δᵢ = 解的项数
```

**关键性质（Monagan-Tuncer 2020 证明）**：对稀疏输入多项式（T 个非零项），MTSHL 运行时间为多项式时间 O(r · n · d · T²)，而稠密算法最坏为指数时间（在 T << d^n 时）。

### 3.2 与当前代码的改动关系

| 函数 | 当前 | MTSHL 改动 |
|------|------|-----------|
| `__multivar_diophantine` | 递归稠密 pseudo-division | 改为稀疏 Zippel 插值基本情形 |
| `__multivar_hensel_lift` | 顺序调用 `__hensel_lift_one_var` | 可能需要调整求值点管理（插值需额外点） |
| `__wang_core` | 稠密控制流 | 可能需要"项数估计"来决定何时用稀疏路径 |
| `__hensel_lift_one_var` | 不变（外层循环结构保留） | 接口可保持兼容 |

**最小改动原则**：P2a 的核心改动在 `__multivar_diophantine` 内部。`__hensel_lift_one_var` 的接口（调用签名）可以保持不变，只替换内部求解策略。

---

## 4. 模块划分（待细化，依赖调研）

当前可识别的目标模块：

```
P2a: __multivar_diophantine（改造为稀疏插值版本）
    ├── [新] __dioph_probe_eval(h, G_hats, beta[])
    │         → 在若干随机点求值，构造线性方程组
    │         → 返回 δᵢ(β_j) 的候选值
    ├── [新] __zippel_interpolate(samples, vars, deg_bound)
    │         → 从点值重建稀疏多项式
    │         → 返回 δᵢ 的单项式结构
    └── [改] 基本情形：仅在 δᵢ 结构确认后做精确除法验证

P2b: 控制流调整（待调研后确定）
    └── 求值点管理：MTSHL 需要比 Wang 更多的随机点
```

> **注**：以上模块划分是基于 Monagan-Tuncer 2016 摘要的推断，具体设计需读原文后修订。

---

## 5. 与 P1 的关系

| 维度 | P1（单变量） | P2（多变量） |
|------|------------|------------|
| 目标文件 | `polynomial_factorize_univar.hh` | `polynomial_factorize_wang.hh` |
| Hensel 方向 | p-adic（模精度提升） | x-adic（变量展开） |
| 核心改动 | `__factor_recombine` / `__hensel_lift` | `__multivar_diophantine` |
| 共享代码 | 无 | 无 |
| 依赖关系 | 无 | 无 |

P1 和 P2 完全独立，可并行调研，按性价比顺序实现。P1a（van Hoeij）是当前最大瓶颈，建议优先。

---

## 6. 参考文献

| 文献 | 角色 |
|------|------|
| Monagan & Tuncer (2016). *CASC 2016*, LNCS 9890:171–186 | 核心算法（稀疏 Hensel + Zippel 插值）|
| Monagan & Tuncer (2018). *ICMS 2018*, LNCS 10931:378–387 | 高性能实现（Maple 实际采用版本）|
| Monagan & Tuncer (2020). *J. Symbolic Comput.* 99:189–230 | 复杂度证明（稀疏情形多项式时间）|
| Monagan & Tuncer (2019). *ISSAC 2019* | Maple 2019 算法全貌（单/多变量管线综述）|
| 调研报告 §5.2 | Maple P2 路线摘录（一手资料待补充）|

---

## 7. 下一步

1. 读取 Monagan-Tuncer (2016) 及 (2018) 原文，确认 Zippel 插值嵌入细节
2. 调研 FLINT `fmpz_mpoly_factor` 中稀疏多变量因式分解的实现策略（FLINT 不一定用 MTSHL，但可作参照）
3. 撰写 `docs/research/mtshl-research.md`（对应 P1 的 `vanhoeij-lll-research.md`）
4. 基于调研结果细化 P2a 模块划分（补充完整的功能/接口规约，参照 P1a M1-M5 格式）
