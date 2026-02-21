# 调研报告：van Hoeij LLL 因子重组

> 调研日期：2026-02-21
> 目标：替换 `__factor_recombine` 中的 Zassenhaus O(2^r) 子集枚举

---

## 1. 现有基础设施

### 1.1 当前单变量因式分解管线

```
factorize(f ∈ Z[x])
  → squarefreefactorize(f)
  → 对每个无平方因子 g:
      → __select_prime(g)          §8.3  选择模因子数最少的素数 p（试 30-50 个）
      → __hensel_lift(g, factors, p)  §6.6  二次 Hensel 提升（mod p^k，k 使 m > 2·lc·B_Mig）
      → __factor_recombine(g, lifted, m)  §7.5  ← 目标替换
```

### 1.2 `__factor_recombine` 现状（polynomial_factorize.hh:1922）

**算法**：Zassenhaus 子集枚举
- 按子集大小 s = 1, 2, ... 枚举所有 C(r, s) 个子集
- 剪枝 1：首项系数必整除 lc(f*)²（mod m 后检查）
- 剪枝 2：常数项乘积必整除 f* 常数项（mod m 后检查）
- 完整验证：计算子集乘积，试除 f*

**复杂度**：O(2^r) 最坏情形，r = 模因子数

**接口**：
```cpp
std::vector<upolynomial_<ZZ>>
__factor_recombine(
    const upolynomial_<ZZ>& f,
    const std::vector<upolynomial_<ZZ>>& lifted,  // Hensel 提升后因子 mod m
    const ZZ& m)
```

**已有辅助函数（可复用）**：
- `__symmetric_mod(a, m)`：对称模（输出在 (-m/2, m/2]）
- `__upoly_symmetric_mod(f, m)`：多项式各项对称模
- `__mignotte_bound(f)`：Mignotte 界 B（系数绝对值上界）
- `__upoly_primitive(f)`：本原化（提内容 + lc > 0）
- `pair_vec_div`：精确整数试除

---

## 2. van Hoeij 算法调研

### 2.1 核心思想（van Hoeij 2002）

**问题**：给定 r 个提升因子 h_1,...,h_r（mod p^k），每个真因子 g 对应 {0,1}^r 的 indicator 向量 v（v_i = 1 iff h_i | g mod p^k）。目标：找出所有这样的 v。

**关键洞察**：

设 f 的真因子 g 对应 S ⊆ {1,...,r}，则：
- g 的系数（对称模后）可从 ∑_{i∈S} "h_i 的信息" 推导
- 这产生线性约束：v 必须满足若干整数线性方程

**格框架**：

定义格 L ⊆ Z^r，其中包含所有合法的 indicator 向量 v（真因子子集对应的 {0,1}^r 向量）。

初始格：L_0 = Z^r（平凡，无约束）

逐步加约束：每次喂入一个"数据行"（Newton 迹或 CLD 系数），格缩小。

最终：L = span(v_1, ..., v_k)，其中 v_i 对应真因子 g_i。

**LLL 规约**：

格中合法向量范数小（O(√r)）。LLL 找到短向量后：
- 短向量的 {0,1} 分量近似等于某个真因子的 indicator
- 验证：用 indicator 指定的子集乘积试除 f，确认因子关系

### 2.2 格矩阵的具体构造

#### 方案 A：Newton 迹（van Hoeij 2002 原版）

**Newton 幂和**（Power sums）：
```
p_j(h_i) = Σ_{α: h_i(α)≡0 (mod p^k)} α^j  (mod p^k)
```

可通过 Newton 恒等式从 h_i 的系数 O(n) 时间递推计算：
```
p_j + e_1·p_{j-1} + ... + j·e_j = 0   (j ≤ deg h_i)
p_j + e_1·p_{j-1} + ... + e_n·p_{j-n} = 0  (j > deg h_i)
```
其中 e_k 是 h_i 的初等对称多项式（即系数±1）。

**矩阵构造**（r 个提升因子，J 列 Newton 迹）：

```
M = [p^k · I_r | 0_{r×J}  ]  ← 格的初始生成元（p^k · Z^r 部分）
    [T          | W · I_J  ]  ← 每行对应一个 Newton 迹约束
```

其中：
- T[j, i] = p_j(h_i) mod p^k（对称约化），j = 1,...,J，i = 1,...,r
- W = 一个大权重（使格向量的迹部分较小），通常 W = p^k 或者类似量级
- 矩阵尺寸：(r + J) × (r + J)

真因子 g 对应的 indicator 向量 v = (v_1,...,v_r,0,...,0) 在格中：
- 前 r 个分量：v_i ∈ {0,1}
- 后 J 个分量：T · v（mod 某个大数），约为 1/W 级别（应为短向量）

LLL 后，短向量（‖b‖² 小）的前 r 个分量给出 indicator。

#### 方案 B：CLD（Coefficient Logarithmic Derivative，Belabas-van Hoeij-Klueners-Steel 2009）

**思路**：用多项式系数（而非代数根的幂和）描述约束。

对真因子 g = ∏_{i∈S} h_i，有：
```
f · g'/g = Σ_{i∈S} f · h_i'/h_i   (mod p^k)
```

即 f 乘以 g 的对数导数等于各 h_i 对数导数的 v_i 加权和。

取若干系数（低次端和高次端），得到数据行：
- 每行：`[coeff_j(f · h_1'/h_1), ..., coeff_j(f · h_r'/h_r)]`（mod p^k 后对称约化）
- 拼入格矩阵（类似方案 A）

**优势**：
- 系数界更紧（CLD 给出更准确的 Mignotte 型界）
- 可避免"根"的概念，纯多项式运算

**FLINT 的实现**：
- 函数 `_fmpz_poly_factor_CLD_mat()` 计算此数据矩阵
- 取低次端（`lo_n` 列）和高次端（`hi_n` 列）各若干列，排除中间（界太大）
- 每轮取 30-50 列（`num_coeffs` 参数），每次迭代翻倍

### 2.3 LLL 规约

**目标**：对格矩阵做 LLL，找出所有范数约为 √r 的短向量。

**FLINT 的 LLL 变体**：`fmpz_lll_wrapper_with_removal_knapsack`

**"knapsack removal" 优化**：
- LLL 规约后，若某个格向量 b 满足 ‖b‖² ≤ B（阈值），则其前 r 个分量解读为 indicator v
- 验证 v 是否为合法因子（试除）。若合法，从格中移除对应行，矩阵缩小
- 反复进行直到矩阵无短向量

**参数**：
- δ = 3/4（标准 Lovász 条件）
- B ≈ (r + 1) · 2^(2 · U_exp)，其中 U_exp = ⌊log₂(r)⌋（FLINT 默认）

### 2.4 渐进喂入（Hart-van Hoeij-Novocin 2011）

**思想**：不一次性构建大矩阵，而是逐步加列再做 LLL。

```
初始：M = p^k · I_r  (r×r)

循环：
  1. 计算 num_cols 列 CLD/Newton 数据
  2. 将这些列拼入 M（M 扩展为 (r+J) × (r+J+num_cols)）
  3. 做 LLL with removal
  4. 若找到短向量 → 验证并提取因子，缩小矩阵
  5. 未找到 → 增大 num_cols（翻倍），精度翻倍（Hensel 再提升），重试
  6. 直到矩阵秩降到 1（所有因子找到）
```

**好处**：
- 早期就能发现易识别的因子（如次数小的因子）
- 避免一次性构建昂贵的大矩阵
- FLINT 初始 num_cols = 30-50，上界 = (deg_f + 1) / 2

---

## 3. LLL 实现方案对比

| 方案 | 优点 | 缺点 | 适用 CLPoly？ |
|------|------|------|--------------|
| **fplll 外部库** | 工业级，支持 BKZ、HLLL | 新增 C++ 依赖（~100KB） | 可以，但侵入性强 |
| **自写整数 LLL（ZZ）** | 无新依赖，精确，与 ZZ 直接集成 | 实现成本高（需 Gram-Schmidt over Q） | 推荐初版 |
| **浮点 LLL（double）** | 实现简单 | 精度问题（大矩阵可能失败） | 不推荐 |

**结论**：初版自写整数 LLL（基于 ZZ），矩阵规模通常 r ≤ 50，足够精确。

**参考实现**：
- NTL `src/LLL.cpp`（`LLL_XD` 算法：以 `xdouble` 扩精度浮点做 Gram-Schmidt 近似，失败时回退精确有理数）
- FLINT `src/fmpz_lll/lll.c`（纯整数 LLL，ZZ 系数 Gram-Schmidt）
- Griffiths & Steinfeld (1996)：整数 LLL 无浮点参考

---

## 4. 实现方案推荐

### 4.1 算法选择

**第一版**：Newton 迹 + 渐进喂入 + 自写整数 LLL

理由：
- Newton 迹比 CLD 概念更简单，便于正确实现和验证
- 渐进喂入是关键优化（避免 Zassenhaus 指数爆炸）
- 自写整数 LLL 消除外部依赖

**第二版**（可选）：将 Newton 迹替换为 CLD，界更紧，性能更优。

### 4.2 与现有代码的接口

`__factor_recombine` 的签名保持不变，内部替换实现：
```cpp
std::vector<upolynomial_<ZZ>>
__factor_recombine(
    const upolynomial_<ZZ>& f,
    const std::vector<upolynomial_<ZZ>>& lifted,
    const ZZ& m)
```

辅助函数（新增，局限在 `polynomial_factorize.hh` 内部）：
- `__newton_power_sum(h, j, m)` → p_j(h) mod m（对称）
- `__lll_reduce(M)` → 原地 LLL 规约（ZZ 矩阵）
- `__vanhoeij_recombine(f, lifted, m)` → 核心循环（渐进喂入 + LLL）

### 4.3 回退策略

当 r ≤ 10（或 deg_f ≤ 20）时，仍使用现有 Zassenhaus（已有剪枝，性能可接受）。
当 r > 10 时，使用 van Hoeij。
阈值可通过 benchmark 微调。

### 4.4 测试策略

1. 回归：现有 256 个 test_factorize 用例全部通过
2. 性能对比：bench_comparative 中 "70 factors" 用例，目标 < 5 ms（对比 FLINT 3.6 ms）
3. 压力测试：stress_factorize 4 个测试用例全部通过（release 模式）

---

## 5. 外部参考文献

| 文献 | 内容 |
|------|------|
| van Hoeij 2002, J. Number Theory 95:167-189 | 算法原版，格构造，Newton 迹方法 |
| Hart-van Hoeij-Novocin 2011, ISSAC | 渐进喂入（gradual feeding）优化 |
| Belabas-van Hoeij-Klueners-Steel 2009 | CLD 替代 Newton 迹，界更紧 |
| Klueners, "The van Hoeij Algorithm" (SpringerLink) | 算法综述，清晰的教学版伪代码 |
| FLINT `src/fmpz_poly_factor/factor_van_hoeij.c` | 参考实现（CLD + 渐进喂入 + knapsack LLL） |
| FLINT `src/fmpz_poly_factor/CLD_mat.c` | CLD 矩阵构造参考 |
| FLINT `src/fmpz_lll/` | 整数 LLL 参考实现 |

---

## 6. 结论

**推荐方案**：van Hoeij + Newton 迹 + 渐进喂入 + 自写整数 LLL。

实现范围：
- 新增 `__lll_reduce` 模块（~150-200 行 ZZ 矩阵 LLL）
- 新增 `__newton_power_sum` 函数（~30 行，Newton 恒等式递推）
- 新增 `__vanhoeij_recombine`（~150 行，主循环）
- `__factor_recombine` 路由逻辑（小 r 回退 Zassenhaus，大 r 用 van Hoeij）

预期改善：单变量 70 因子 80 ms → < 5 ms（~16x），与 FLINT 差距缩小到 <2x。
