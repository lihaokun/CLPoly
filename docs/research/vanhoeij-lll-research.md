# 调研报告：van Hoeij LLL 因子重组

> 调研日期：2026-02-21
> 目标：替换 `__factor_recombine` 中的 Zassenhaus O(2^r) 子集枚举
> **一手资料状态**：§2.2–§2.4 已基于 FLINT 源码（`factor_van_hoeij.c`、`CLD_mat.c`）核实并修正；van Hoeij 2002 原文 PDF 未能直接获取，算法逻辑以 FLINT 实现为准。

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

### 1.2 `__factor_recombine` 现状（`polynomial_factorize_univar.hh`）

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
p_j - e_1·p_{j-1} + e_2·p_{j-2} - ... + (-1)^{j-1} e_{j-1}·p_1 + (-1)^j j·e_j = 0  (j ≤ d)
p_j - e_1·p_{j-1} + e_2·p_{j-2} - ... + (-1)^d e_d·p_{j-d} = 0                      (j > d)
```
其中 e_k 是根的第 k 个初等对称多项式，由 Vieta 公式从首一因子 h_i(x) = x^d + c_{d-1}x^{d-1} + ⋯ + c_0 直接读出：
```
e_k = (-1)^k · c_{d-k}    (k = 1,...,d)
```
即 e_1 = -c_{d-1}，e_2 = c_{d-2}，e_3 = -c_{d-3}，以此类推。

**矩阵构造**（r 个提升因子，J 列数据约束）：

> 注：以下基于 FLINT `factor_van_hoeij.c` 实现，Newton 迹仅为原始算法（§方案 A）的理论描述；FLINT/Maple 实际使用 CLD 数据（§方案 B），矩阵结构相同。

```
M = [2^U_exp · I_r | 0_{r×J}  ]  ← 格的初始生成元（预缩放单位矩阵）
    [T              | I_J      ]  ← 每行对应一列数据约束
```

其中：
- `U_exp = FLINT_BIT_COUNT(max(r, 20))`（即 max(r,20) 的二进制位数）
- T[j, i] = 第 j 列数据对第 i 个因子的约束值（对称约化至 (-p^k/2, p^k/2]）
- 下半矩阵 I_J 权重为 1（不额外缩放）
- 矩阵尺寸：(r + J) × (r + J)

**矩阵预缩放的作用**：`2^U_exp` 不是 Hensel 模数 p^k，而是一个精度预缩放因子。
- 使 `2^U_exp · I_r` 的基向量范数大于真 indicator 向量（范数 ~√r），令 LLL 能区分二者
- 与短向量阈值 B = (r+1) · 2^(2·U_exp) 配合使用（见 §2.3）

**短向量分析**：真因子 g 对应 indicator v = (v_1,...,v_r) ∈ {0,1}^r，对应格中向量范数 ~√r，远小于 `2^U_exp` 量级的非合法组合，LLL 能区分并找到它。

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

> **记法注**：`f · h_i'/h_i` 与 `(f/h_i) · h_i'` 在数学上等价（因 Hensel 提升保证 h_i | f mod m，精确可除）。后者（先精确除法再乘导数）是实现时更自然的写法，见 §4.1 及架构文档 M1。

**优势**：
- 系数界更紧（CLD 给出更准确的 Mignotte 型界）
- 可避免"根"的概念，纯多项式运算

**FLINT 的实现**（`CLD_mat.c`，一手来源）：
- 函数 `_fmpz_poly_factor_CLD_mat()` 构造 (r+1)×2k 数据矩阵：前 r 行为各因子数据，末行为各列 CLD 系数界
- 取低次端（低次系数列）和高次端（高次系数列），排除中间部分（界过大、信噪比低）
- 列过滤：若某列界 × √N > p^a / 2^(1.5r)，则跳过该列
- 每轮初始列数见 §2.4

### 2.3 LLL 规约

**目标**：对格矩阵做 LLL，找出所有范数² ≤ B 的短向量（对应真因子 indicator）。

**FLINT 的实现**（`factor_van_hoeij.c`，一手来源）：

```c
// 调用方式：
num_rows = fmpz_lll_wrapper_with_removal_knapsack(M, NULL, B, fl);
// B = (r + 1) * 2^(2 * U_exp)，其中 U_exp = FLINT_BIT_COUNT(max(r, 20))
```

- **Knapsack removal 内嵌于 LLL**：不是后处理步骤，而是 LLL 规约过程中一旦发现范数² ≤ B 的向量即触发移除
- LLL 返回后矩阵自动缩减至 num_rows 行（已移除短向量对应行）
- 因子提取：通过 `fmpz_mat_col_partition(U)` 对 LLL 的幺模变换矩阵 U 做列分组（partition），每个分组对应一个候选因子的提升因子子集——**不是直接读取向量的 {0,1} 分量**

**`fmpz_mat_col_partition` 算法（FLINT 文档 + 上下文推断）**：

FLINT 文档定义：
```
fmpz_mat_col_partition(part[], M, short_circuit):
  返回 M 的不重复列数 p
  part[j] ∈ [0, p)，使得 part[i] == part[j] ⟺ M 的第 i 列完全相等
  （short_circuit = true 时：若 p > 行数则提前返回 0）
```

在 van Hoeij 上下文中的推断语义（基于算法逻辑，具体调用方式需读 FLINT `factor_van_hoeij.c` 确认）：

设 `U_short` 为 LLL 的幺模变换矩阵 U 中对应已移除短向量的行，限制到前 r 列（第 j 列对应模因子 h_j）：

- `U_short[:,j]` 记录了模因子 h_j 在所有已找到 indicator 向量中的参与情况
- **若 `U_short[:,i] == U_short[:,j]`**：h_i 和 h_j 在所有已找到的真因子中 indicator 模式完全相同——它们永远共同出现或共同不出现 → 属于同一真因子的提升因子子集
- 每个等价类（`part[k] == c`）对应一个候选真因子，其提升因子子集 = `{h_k : part[k] == c}`

> **CLPoly `__extract_candidates` 实现指引**：按 `part[]` 对前 r 列分组，每组子集乘积即一个候选因子，按 §2.5 方式验证。**具体 U 矩阵的构造和取行方式**待实现时读 FLINT `src/fmpz_poly_factor/factor_van_hoeij.c` `_fmpz_poly_factor_van_hoeij_partial` 函数确认（约 40 行内）。

**参数**（均来自 FLINT 源码）：
- δ = 3/4（标准 Lovász 条件，`fmpz_lll_context_init_default`）
- B = (r + 1) · 2^(2 · U_exp)，U_exp = FLINT_BIT_COUNT(max(r, 20))

### 2.4 渐进喂入（Hart-van Hoeij-Novocin 2011）

**思想**：不一次性构建大矩阵，而是逐步加列再做 LLL。

**FLINT 实际流程**（`factor_van_hoeij.c`，一手来源）：

```
初始：M = 2^U_exp · I_r  (r×r)，U_exp = FLINT_BIT_COUNT(max(r, 20))

外层循环（Hensel 精度倍增）：
  初始 num_coeffs:
    若 3*r > N+1（因子数多于多项式结构允许量）：num_coeffs = 50（r>200）或 30
    否则：num_coeffs = 10
  num_coeffs = min(num_coeffs, (N+1)/2)

  内层循环（列数翻倍直到收敛）：
    1. 计算 CLD 数据矩阵（num_coeffs 列，过滤噪声列）
    2. 按内螺旋顺序逐列喂入 M：
         奇数步取左端列，偶数步取右端列（从两端向中间交替）
         每列喂入后若触发 LLL：运行 fmpz_lll_wrapper_with_removal_knapsack
         若 LLL 找到短向量 → 验证候选因子（模检验 + 精确试除）
         若分解完成 → 返回结果
    3. num_coeffs *= 2，重复直到 num_coeffs 不变（上界收敛）

  若内层循环未完成分解：Hensel 精度加倍（a *= 2），重复外层循环
```

> **CLPoly 简化注**：FLINT 在每列喂入后即可触发 LLL（细粒度控制），且外层有精度倍增安全网（`a *= 2`）。CLPoly 初版架构做两处简化：① 喂入全部 J_target 列后统一运行一次 LLL；② 若 J_target 超出 J_max = (N+1)/2 则回退 Zassenhaus（替代 FLINT 的精度倍增）。正确性等价，性能略低，见架构文档 M2-M3 节。

**关键参数来源（FLINT 实测）**：
- 初始 num_coeffs：10 或 30/50（依条件 `3r > N+1`）
- 列上界：(N+1)/2，其中 N = deg(f)
- 列喂入顺序：内螺旋（左右交替），不是顺序喂入，目的是尽早暴露因子结构
- 内层循环终止：num_coeffs 达到上界后不再增长（循环收敛）

### 2.5 Hensel 精度与首项系数处理

**Hensel 精度**：van Hoeij 对精度要求与 Zassenhaus 相同——需提升至 m > 2·lc(f)·B_Mig（Mignotte 界）。精度只需在进入 LLL 循环前一次性确定，渐进喂入过程中增加的是"约束列数"（Newton 迹/CLD 列），而非精度。

> 注：§2.4 循环第 5 步"精度翻倍"是指当格约简无法找到短向量时（极少发生），回退增大精度的安全网；正常情况下不需要。**对 P1a 而言**，CLPoly 现有 `__hensel_lift` 的二叉树提升架构无需修改，提升一次到完整精度即可。**P1b** 将修改此架构实现"提升-LLL 交织"，精度可早停，见 §5.4。

**首项系数处理**：Hensel 提升后，每个 h_i 在 Z/m 中的首项系数满足 ∏ lc(h_i) ≡ lc(f) (mod m)。确认 indicator v 对应候选因子时，构造方式与 Zassenhaus 相同：
```
g_trial = lc(f) · ∏_{v_i = 1} h_i   (mod m)
```
对称约化后，先试整除 lc(f)（提取内容），再用 `pair_vec_div` 试除 f，完整验证。实现上可直接复用 `__factor_recombine` 中现有的候选构造与验证代码，只替换"枚举子集"这一步为"读取 LLL 短向量"。

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
- Cohen, H. (1993). *A Course in Computational Algebraic Number Theory*. §2.6（整数 LLL，基于有理 Gram-Schmidt，无浮点；含完整伪代码）

**CLPoly `__lll_reduce` 实现选型**：使用 **QQ 有理 Gram-Schmidt**（Cohen §2.6 算法 2.6.3）：

- **选型理由**：
  - r ≤ 50，QQ 分数大小可控（μ_{ij} 的分子分母在 GMP 处理范围内）
  - CLPoly 已有 QQ 基础设施，无额外依赖，可直接按 Cohen §2.6 伪代码照写
  - 精确算术消除浮点失效风险（浮点 LLL 对近病态格可能给出错误结果，需要复杂回退逻辑）
  - Cohen §2.6 含完整伪代码 + Lovász 条件 (δ = 3/4)，实现目标明确

- **数据结构**（`__lll_reduce` 内部）：
  - `mu[i][j]`（i > j）：`QQ` 型 Gram-Schmidt 系数
  - `B[i]`：`QQ` 型 Gram-Schmidt 平方范数（B[i] = ‖b*_i‖²）
  - `M`：`std::vector<std::vector<ZZ>>` 整数格矩阵（原地 LLL 规约）
  - `U`：`std::vector<std::vector<ZZ>>` 幺模变换矩阵，与 M 同步行操作，LLL 结束后 M_new = U · M_old

- **幺模变换 U 的用途**：LLL 结束后，U 的各行对应格规约后的各基向量是原始基向量的哪种整数线性组合。短向量行（norm² ≤ B 的行）在 U 中对应的列信息，用于 `__extract_candidates` 的因子子集恢复（见 §2.3）。

- **不推荐**：浮点 Gram-Schmidt（double/xdouble）在 r > 30 时可能失效，FLINT 纯整数路线实现更复杂，对 r ≤ 50 性能优势不明显；两者均不适合 CLPoly 初版。

---

## 4. 实现方案推荐

### 4.1 算法选择

**方案**：CLD + 渐进喂入（内螺旋列喂入）+ 自写整数 LLL

理由：
- CLD 约束界不随列索引指数增长，收敛所需列数更少（FLINT 实测 30-50 列即可）
- Maple/FLINT 均采用 CLD，CLPoly 与参考实现路线一致
- CLD 计算用 `C_i = (f/h_i)·h_i'`（精确除法），无需多项式幂级数求逆，实现并不复杂
- 渐进喂入减少每轮 LLL 的矩阵规模，加速收敛（van Hoeij 整体算法替换了 Zassenhaus 枚举，渐进喂入是 van Hoeij 内部的效率优化）
- 自写整数 LLL 消除外部依赖

> Newton 迹（van Hoeij 2002 原版）已被 CLD 取代，不再作为实现选项。

### 4.2 与现有代码的接口

`__factor_recombine` 的签名保持不变，内部替换实现：
```cpp
std::vector<upolynomial_<ZZ>>
__factor_recombine(
    const upolynomial_<ZZ>& f,
    const std::vector<upolynomial_<ZZ>>& lifted,
    const ZZ& m)
```

辅助函数（新增，局限在 `polynomial_factorize_univar.hh` 内部）：
- `__cld_polys(f, active, m)` → C_i = (f/h_i)·h_i' for all i（CLD 多项式）
- `__build_cld_matrix(M, cld, J_cur, J_target, m)` → 内螺旋喂列，扩展格矩阵
- `__lll_reduce(M, U, B)` → 原地 LLL 规约，返回幺模变换 U + 短向量行下标
- `__extract_candidates(short_rows, U, r)` → 从 U 提取候选因子子集
- `__vanhoeij_recombine(f, lifted, m)` → 核心循环（CLD + 渐进喂入 + LLL）

### 4.3 回退策略

当 r ≤ 10 时，仍使用现有 Zassenhaus（已有剪枝，性能可接受）。
当 r > 10 时，使用 van Hoeij。
阈值（ZASSENHAUS_THRESHOLD = 10）为编译期常量，可通过 benchmark 微调。

### 4.4 测试策略

1. 回归：现有 256 个 test_factorize 用例全部通过
2. 性能对比：bench_comparative 中 "70 factors" 用例，目标 < 5 ms（对比 FLINT 3.6 ms）
3. 压力测试：stress_factorize 4 个测试用例全部通过（release 模式）

---

## 5. Maple 算法路线参考

> CLPoly 以 Maple 为核心算法标杆。以下记录 Maple 因式分解的演进路线，供选型和优先级决策参考。

### 5.1 单变量 Z[x] 因式分解（P1 直接参考）

Maple 的单变量算法分两个演进阶段：

**阶段一：van Hoeij LLL 重组（替换 Zassenhaus 枚举）**

| 文献 | 内容 |
|------|------|
| van Hoeij (2002). *J. Number Theory* 95:167–189 | **算法原版**：格构造 + Newton 迹（Power sums）+ knapsack 规约 |
| Belabas, van Hoeij, Klüners, Steel (2009). *J. Théor. Nombres Bordeaux* 21(1):15–39 | **CLD 改进**：用 `f·h_i'/h_i` 的系数替代 Newton 迹；界更紧；**FLINT/Maple 当前均采用 CLD** |
| Klüners (2010). In: *The LLL Algorithm* (Springer), pp. 283–291 | 教学综述，含完整伪代码 |

**阶段二：线性 Hensel 提升 + 早期因子检测（Maple 2019 新增）**

| 文献 | 内容 |
|------|------|
| Monagan (2019). *ISSAC 2019*, pp. 275–282 | **线性 Hensel 提升**：将提升复杂度从 O(n²d²) 降至 O(nd²)（n = deg(f)，d = 目标模数 bit 数）；在提升过程中即可做早期因子试除，无需等到提升完成 |
| Monagan & Tuncer (2019). *Polynomial Factorization in Maple 2019* | 描述 Maple 2019 单/多变量算法全貌：单变量用线性提升 + van Hoeij LLL；多变量用 MTSHL（亦见 §5.2） |

> **对 CLPoly P1 的含义**：CLPoly 当前 `__hensel_lift` 用二叉树二次提升。Maple 2019 将 LLL 重组与线性 Hensel 提升**交织使用**（提升-检测互嵌）。P1a（LLL 重组）是收益最大的改动；P1b（线性 Hensel）进一步消除精度过提升。两者均属 P1 单变量优化范畴，定位分析见 §5.4。

### 5.2 多变量 Z[x₁,...,xₙ] 因式分解（P2 参考）

Maple 2019 起默认使用稀疏 Hensel 提升（MTSHL 算法），对应 P2 优化：

| 文献 | 内容 |
|------|------|
| Monagan & Tuncer (2016). *CASC 2016*, LNCS 9890:171–186 | **核心思路**：在 Hensel 提升内嵌入 Zippel 稀疏插值，解多变量 Diophantine 方程 |
| Monagan & Tuncer (2018). *ICMS 2018*, LNCS 10931:378–387 | **高性能实现**：稀疏 Hensel 提升的工程优化，Maple 实际采用 |
| Monagan & Tuncer (2020). *J. Symbolic Comput.* 99:189–230 | **复杂度分析**：稀疏情形多项式时间界（CLPoly P2 设计的理论依据） |
| Monagan & Tuncer (2019)（同 §5.1） | 多变量 MTSHL 管线描述，与单变量线性 Hensel 合并发布 |

### 5.3 Hensel 提升本身的改进（备查）

| 文献 | 内容 |
|------|------|
| Monagan (2019). *ISSAC 2019*, pp. 275–282 | **线性 Hensel 提升**：Z[x] 提升由 O(n²d²) 降至 O(nd²)；与 van Hoeij LLL 交织（P1b 直接参考） |
| Monagan (2022). *ISSAC 2022* | **n 个因子线性提升**：进一步推广至任意多因子场景 |

### 5.4 P1b 定位分析：单变量线性 Hensel 属于 P1 还是 P2？

> **结论**：P1b 属于 P1，不应与 P2 多变量 Hensel 合并。

#### 5.4.1 CLPoly 中存在两个正交的 Hensel 提升

"Hensel 提升"在 CLPoly 中有两套互不相关的实现，解决两个完全不同的问题：

| 维度 | **P1b 目标** | **P2 目标** |
|------|-------------|------------|
| 函数 | `__hensel_lift` (polynomial_factorize_univar.hh:495) | `__hensel_lift_one_var` (polynomial_factorize_wang.hh:698) |
| 提升方向 | **p-adic**：模精度 `p → p² → p⁴ → ...` | **x-adic**：变量展开 `Z[x₁..xₖ₋₁] → Z[x₁..xₖ]` |
| 当前策略 | 二叉树二次提升（每步精度平方） | Taylor 逐阶展开（每步增加 xₖ 的一个次项） |
| 目标策略 | Monagan 2019：线性精度提升 + 早期因子检测 | Monagan-Tuncer MTSHL：Zippel 稀疏插值嵌入 Diophantine 求解 |
| 数据结构 | `__hensel_node` 二叉树（s/t/g/h） | `G[]` 多项式数组 + `bezout_s[]` Bézout 链 |
| 调用链 | `__factor_squarefree_primitive_ZZ` → `__hensel_lift` → `__factor_recombine` | `__wang_core` → `__multivar_hensel_lift` → `__hensel_lift_one_var` |
| 所在文件 | polynomial_factorize_univar.hh | polynomial_factorize_wang.hh |
| 共享代码 | **无** | **无** |

**核心区分**：p-adic 提升处理的是"同一个因子从小模数升到大模数"；x-adic 提升处理的是"从低变量多项式扩展到多变量多项式"。算法方向正交，代码完全不重叠。

#### 5.4.2 P1b 与 P1a 的协同关系（是设计约束，不是理由合并）

Maple 2019 将 P1b 和 P1a **交织使用**，形成"提升-检测互嵌"模式：

```
线性 Hensel 过程（Maple 2019）：
  步骤 1: h_i mod p¹ → h_i mod p²     （精度指数 +1）
  步骤 2: h_i mod p² → h_i mod p³     （精度指数 +1）
  ...
  每 k 步：以当前精度构建 CLD 矩阵 → 运行 LLL
    → 若 LLL 找到完整短向量（因子分组收敛）：
        立刻停止提升，试除验证
        → 精度只需提升到"LLL 能收敛"，而非完整 Mignotte 界
    → 若 LLL 尚未收敛：继续提升
```

这与现有 CLPoly 的顺序架构（先完整提升到 Mignotte 精度，再做 LLL）不同。这种交织意味着：

- P1b 需要在提升循环内部调用 P1a 的 LLL 模块（`__vanhoeij_recombine` 或其子函数）
- P1a 实现完成后，P1b 才能定义"何时触发 LLL 检测"的接口
- **实施顺序**：必须先完成 P1a，再在其基础上实现 P1b 的交织循环

这是 P1a 和 P1b 都在 P1 的充分理由——它们共享同一个单变量控制流，而非因为"都叫 Hensel"。

#### 5.4.3 P1b 与 P2 无法合并的原因

P2（MTSHL 稀疏多变量提升）需要修改：
- `__hensel_lift_one_var`（变量展开提升主体）
- `__multivar_diophantine`（多变量 Diophantine，嵌入 Zippel 插值）
- `__wang_core`（Wang 算法核心控制流）

P1b（线性单变量 Hensel）需要修改：
- `__hensel_lift`（删除二叉树二次提升，改线性逐步）
- `__factor_squarefree_primitive_ZZ`（调整提升-重组接口为交织模式）

**完全无交集**。若将 P1b 推入 P2：
1. 延误单变量性能改善（P1b 是独立的单变量优化）
2. P2 的稀疏插值工作无法简化 P1b 的任何实现
3. 两套改动不共享测试路径（单变量 vs 多变量测试用例不同）
4. P1b 的"提升-LLL 交织"依赖 P1a 完成，不依赖 P2

#### 5.4.4 预期收益评估

| 优化 | 主要收益 | 量级估算 |
|------|---------|---------|
| **P1a** 仅 van Hoeij LLL | 消除 O(2^r) 枚举 → 多项式时间 | ~10–100x（r=30 时效果最显著） |
| **P1a + P1b** 线性 Hensel 交织 | 精度过提升减少 + 早期因子检测 | 在 P1a 基础上额外 2–5x |
| **P2** MTSHL 稀疏多变量 | 多变量 Hensel 稀疏化（稠密 → 稀疏项数） | 多变量场景 5–50x（依稀疏度） |

P1b 对单变量的额外贡献相对 P1a 较小，但实现成本也较低（P1b 只需改 `__hensel_lift`，约 150 行），且为后续 P1a/P1b 交织奠定基础。

---

## 6. 其他参考文献

| 文献 | 内容 |
|------|------|
| Hart, van Hoeij, Novocin (2011). *ISSAC 2011* | 渐进喂入（gradual feeding）优化 |
| Cohen, H. (1993). *A Course in Computational Algebraic Number Theory*, §2.6 | 整数 LLL 完整伪代码（无浮点） |
| FLINT `src/fmpz_poly_factor/factor_van_hoeij.c` | 参考实现（CLD + 渐进喂入 + knapsack LLL） |
| FLINT `src/fmpz_poly_factor/CLD_mat.c` | CLD 矩阵构造参考 |
| FLINT `src/fmpz_lll/` | 整数 LLL 参考实现 |

---

## 7. 结论

**推荐方案**：van Hoeij + CLD + 渐进喂入（内螺旋）+ 自写整数 LLL。

实现范围：
- 新增 `__cld_polys`（~30 行，精确除法 + 导数 + 乘法）
- 新增 `__build_cld_matrix`（~60 行，内螺旋列选取 + 矩阵扩展）
- 新增 `__lll_reduce`（~150-200 行，ZZ 矩阵 LLL + 幺模变换记录）
- 新增 `__extract_candidates`（~30 行，U 矩阵列分组）
- 新增 `__vanhoeij_recombine`（~100 行，主循环）
- `__factor_recombine` 路由逻辑（r ≤ 10 回退 Zassenhaus，r > 10 用 van Hoeij）

预期改善：单变量 70 因子 80 ms → < 5 ms（~16x），与 FLINT 差距缩小到 <2x。
