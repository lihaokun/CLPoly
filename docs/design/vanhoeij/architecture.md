# P1 单变量因式分解优化：架构文档

> 阶段：架构（§2.2）
> 目标：将 CLPoly 单变量 Z[x] 因式分解性能对齐 Maple/FLINT
> 目标文件：`clpoly/polynomial_factorize_univar.hh`

## P1 组成

| 组件 | 目标函数 | 当前实现 | 目标实现 | 状态 |
|------|---------|---------|---------|------|
| **P1a** van Hoeij LLL 重组 | `__factor_recombine` | Zassenhaus O(2^r) | CLD + LLL（Belabas 2009） | **架构完成** |
| **P1b** 线性 Hensel 提升 | `__linear_hensel_lift_with_lll` | 二叉树二次提升 → Mignotte 精度 | **n 因子直接**线性提升 + van Hoeij LLL 交织 + 早期因子检测 | **架构完成** |

> P1a 和 P1b 初版实现相互独立，可按序进行。P1a 是当前主要瓶颈（70因子 80ms 的直接原因），优先实现。P1b 详细设计在 P1a 完成后展开（完整形态下两者交织）。

---

## Part A：van Hoeij LLL 重组（P1a）

> 对应调研：`docs/research/vanhoeij-lll-research.md`
> 数据约束方案：**CLD（Coefficient Logarithmic Derivative，Belabas 2009）**

---

## 1. 核心流程（P1a）

```
__factor_recombine(f, lifted[0..r-1], m)
    │
    ├── r ≤ ZASSENHAUS_THRESHOLD (默认 10)
    │       └── 现有 Zassenhaus 子集枚举（保留不变）
    │
    └── r > ZASSENHAUS_THRESHOLD
            └── __vanhoeij_recombine(f, lifted, m)
                    │
                    ├── 初始化：
                    │     f_star = f，active = {0..r-1}
                    │     U_exp = bit_length(max(r, 20))
                    │     B = (r+1) · 2^(2·U_exp)
                    │     M = 2^U_exp · I_r  (初始格，r×r)
                    │     N = deg(f)
                    │     J_max = (N+1) / 2
                    │     J₀ = (3r > N+1) ? 30 : 10，然后 J₀ = min(J₀, J_max)
                    │     J_cur = 0，J_target = J₀
                    │
                    └── 主循环（直到 |active| ≤ 1）：
                            │
                            ├── [M1] __cld_polys(f_star, active因子, m)
                            │         → C_i = (f_star / h_i) · h_i'  for i ∈ active
                            │
                            ├── [M2] __build_cld_matrix(M, cld, J_cur, J_target, m)
                            │         → 内螺旋顺序逐列喂入 M，扩展至 (r+J_cur+J_new)×(r+J_cur+J_new)
                            │         → J_new = 0 时（全列被过滤）：
                            │             倍增 J_target 后重入外层循环
                            │             若 J_target > J_max 则 fallback Zassenhaus
                            │
                            ├── [M3] __lll_reduce(M, U, B)
                            │         → M 原地规约，U 记录幺模变换
                            │         → 返回短向量行下标（‖row‖² ≤ B）
                            │
                            ├── [M4] __extract_candidates(short_rows, U, r)
                            │         → 对每个短向量行：从 U 的对应行读出 active 因子分组
                            │         → 返回候选子集列表
                            │
                            ├── 对每个候选子集 S：
                            │         → g_trial = lc(f_star) · ∏_{i∈S} h_i  (mod m)
                            │         → 对称约化 + 本原化 → pp_g
                            │         → pair_vec_div(f_star, pp_g) → 试除
                            │         → 成功：记录因子，f_star /= pp_g
                            │                   缩减 active（r' = r - |S|）
                            │                   重算 U_exp' = bit_length(max(r', 20))
                            │                        B' = (r'+1) · 2^(2·U_exp')
                            │                   重建 M = 2^U_exp' · I_{r'}（r'×r'）
                            │                   重置 J_cur = 0，J_target = J₀
                            │         → 失败：继续
                            │
                            └── 本轮无新因子：J_target *= 2
                                （若 J_target > J_max 则 fallback Zassenhaus）
```

---

## 2. 模块划分

### M1 `__cld_polys`

**功能规约**

```
函数名称：__cld_polys

功能描述：对所有活跃提升因子计算 CLD 多项式
          C_i = (f / h_i) · h_i'  (在 Z/m[x] 中精确计算，m = p^k)
          即 f 乘以 h_i 的对数导数 h_i'/h_i。

前置条件（Requires）：
  - f 本原，lc(f) > 0
  - 每个 h_i 首一，h_i | f (mod m)（Hensel 提升保证）
  - f 的系数在 ZZ 中（非 mod m）；h_i 的系数在 (-m/2, m/2]

后置条件（Ensures）：
  - 返回 r 个多项式 C_i（系数 ∈ Z，对 m 对称约化至 (-m/2, m/2]）
  - 对真因子 g = ∏_{i∈S} h_i，满足：
      Σ_{i∈S} C_i ≡ f · g'/g  (mod m)   ← CLD 可加性
  - deg(C_i) ≤ deg(f) + deg(h_i) - 1

计算步骤（无需多项式求逆）：
  1. q_i = f / h_i  (Z/m[x] 精确除法，h_i | f 在 mod m 下精确整除)
  2. h_i' = 形式导数（系数 × 次数，mod m）
  3. C_i = q_i · h_i'  (mod m)，结果对称约化

副作用：无。
```

**接口规约**

```
接口：__vanhoeij_recombine → __cld_polys

输入数据：
  - f_star: upolynomial_<ZZ>（当前待分解多项式，ZZ 系数）
  - active_factors: vector<upolynomial_<ZZ>>（mod m，对称约化）
  - m: ZZ（Hensel 模数）

输出数据：
  - cld: vector<upolynomial_<ZZ>>（r 个 CLD 多项式，mod m 对称约化）

协议约定：
  - 调用方保证每个 h_i 精确整除 f_star (mod m)
  - 被调用方保证 cld[i] 的系数在 (-m/2, m/2]
```

---

### M2 `__build_cld_matrix`

**功能规约**

```
函数名称：__build_cld_matrix

功能描述：从 CLD 多项式中选取系数列，按内螺旋顺序逐列喂入格矩阵 M，
          扩展 M 直到达到 J_target 新列（或 CLD 可用列耗尽）。
          不运行 LLL——LLL 由调用方在本函数返回后执行。

前置条件（Requires）：
  - M 当前为 (r + J_cur) × (r + J_cur) 整数矩阵
  - cld 包含 r 个 CLD 多项式
  - J_target > 0

后置条件（Ensures）：
  - M 扩展为 (r + J_cur + J_new) × (r + J_cur + J_new)，0 ≤ J_new ≤ J_target
  - 新增的每一列对应一个系数索引 idx：
      M 的第 (r+J_cur+j) 行 = [cld[0][idx], ..., cld[r-1][idx], 0,...,1,...,0]
      （第 r+J_cur+j 列为 1，其余数据位为 0）
  - 列选取顺序：内螺旋（从两端交替向中间），跳过 CLD 界超过阈值的列
  - 若所有可用列均被过滤（J_new = 0），调用方倍增 J_target 后重试；
    若 J_target 超过 J_max 则对剩余因子回退 Zassenhaus

列过滤条件（来自 FLINT CLD_mat.c）：
  - bound_j = max_i |cld[i][idx_j]| · √(deg(f)+1)
  - 若 log2(bound_j) > log2(m / 2^(1.5·r))，跳过该列

副作用：原地修改 M（追加行列）。
```

**接口规约**

```
接口：__vanhoeij_recombine → __build_cld_matrix

输入数据：
  - M: vector<vector<ZZ>>（格矩阵，原地扩展）
  - cld: vector<upolynomial_<ZZ>>（r 个 CLD 多项式，r = cld.size() 内部推导）
  - J_cur: int（已喂入列数）
  - J_target: int（目标新增列数）
  - m: ZZ

输出数据：
  - J_new: int（实际新增列数，≤ J_target；为 0 时表示全列被过滤）

协议约定：
  - 调用方保证 M 已初始化为 2^U_exp · I_r 或之前的扩展状态
  - 被调用方不运行 LLL；LLL（M3）由调用方在本函数返回后调用
  - 若 J_new = 0（全列被过滤），调用方倍增 J_target 重试；若 J_target > J_max 则
    fallback Zassenhaus（CLPoly 无在线 Hensel 精度扩展接口，不在本次范围内）
```

---

### M3 `__lll_reduce`

**功能规约**

```
函数名称：__lll_reduce

功能描述：对整数矩阵 M 做 LLL 基规约（δ = 3/4），同时记录幺模变换矩阵 U，
          返回规约后范数² ≤ B 的行下标。

前置条件（Requires）：
  - M 非空，行数 = 列数 = n（方阵），满秩
  - B ≥ 1

后置条件（Ensures）：
  - M 原地规约为满足 Lovász 条件（δ = 3/4）的 LLL 基
  - U 为 n×n 整数幺模矩阵，满足 M_new = U · M_old
    （U 记录从旧基到新基的坐标变换）
  - 返回 short_rows：所有满足 ‖M[i]‖² ≤ B 的行下标 i（按范数升序）

不变式：
  - 格（行向量张成的 Z-模）在规约过程中保持不变

副作用：原地修改 M；写入 U。
```

**接口规约**

```
接口：__vanhoeij_recombine → __lll_reduce

输入数据：
  - M: vector<vector<ZZ>>（格基，原地修改）
  - B: ZZ（短向量范数²阈值）

输出数据：
  - U: vector<vector<ZZ>>（幺模变换矩阵，n×n）
  - short_rows: vector<int>（范数² ≤ B 的行下标，按范数升序）

协议约定：
  - 调用方保证 M 满秩
  - 被调用方保证 M 规约后仍张成同一格，且 det(U) = ±1

实现规约：
  - Gram-Schmidt 方案：QQ 有理 Gram-Schmidt（Cohen §2.6 算法 2.6.3）
  - 内部数据结构：
      mu[i][j]  : QQ  （Gram-Schmidt 系数 μ_{ij}，i > j，按需更新）
      B_gs[i]   : QQ  （Gram-Schmidt 平方范数 B_i = ‖b*_i‖²）
      M         : vector<vector<ZZ>>  （整数格矩阵，原地规约）
      U         : vector<vector<ZZ>>  （幺模变换矩阵，初始为 I_n）
  - 每次行操作（大小规约 / Lovász 条件交换）同步施加于 M 和 U
  - δ = 3/4（标准 Lovász 条件，参数写死，不对外暴露）
```

---

### M4 `__extract_candidates`

**功能规约**

```
函数名称：__extract_candidates

功能描述：从 LLL 幺模变换矩阵 U 的短向量行中，提取候选因子对应的
          active 因子子集（分组）。

前置条件（Requires）：
  - U 为 LLL 规约后的幺模变换矩阵（n×n，n = r + J）
  - r = 当前 active 因子数

后置条件（Ensures）：
  - 若 short_rows 为空，返回空列表（正常情形：LLL 本轮无收敛）
  - 否则返回 candidates：每项为一个 vector<int>，
    包含 [0, r) 中的 active 因子下标，表示候选因子子集
  - 算法（对应 `fmpz_mat_col_partition` 语义，来自 FLINT 文档）：
      设 U_short = U 中 short_rows 行、前 r 列构成的子矩阵
      对 U_short 的列做相等性分组：
        part[j] = c  ⟺  U_short[:,j] 与第 c 类的列向量完全相等
      等价类 c 对应的候选子集 = { j : part[j] == c }
    （具体 U_short 的取行方式待实现时读 FLINT `factor_van_hoeij.c`
      `_fmpz_poly_factor_van_hoeij_partial` 函数确认，约 40 行）

副作用：无。
```

---

### M5 `__vanhoeij_recombine`

**功能规约**

```
函数名称：__vanhoeij_recombine

功能描述：用 CLD + van Hoeij LLL 方法从 Hensel 提升因子中重组出
          f ∈ Z[x] 的所有不可约因子。

前置条件（Requires）：
  - f 本原，deg(f) ≥ 2，lc(f) > 0
  - lifted 是 f 的 Hensel 提升因子（mod m）
  - m > 2·lc(f)·B_Mig(f)（Mignotte 精度充足）
  - |lifted| = r ≥ 2

后置条件（Ensures）：
  - 返回 g_1,...,g_k ∈ Z[x]，满足 f = g_1 · ... · g_k，各 g_i 本原且 lc > 0
  - 若 f 不可约，返回 {f}

不变式（循环中）：
  - f_star · ∏(result) = f
  - ∏_{i ∈ active} h_i ≡ lc(f_star)^{|active|-1} · f_star  (mod m)

副作用：无（不修改 lifted）。
```

---

### 路由层 `__factor_recombine`（修改现有函数）

```
修改内容：
  - r = |lifted|
  - r ≤ ZASSENHAUS_THRESHOLD：执行现有 Zassenhaus 逻辑（不变）
  - r > ZASSENHAUS_THRESHOLD：调用 __vanhoeij_recombine

常量：ZASSENHAUS_THRESHOLD = 10（编译期常量，可调）
```

---

## 3. 模块依赖与调用关系

```
__factor_recombine
    ├── [r ≤ 10] 现有 Zassenhaus（不变）
    └── [r > 10] __vanhoeij_recombine (M5)
                    ├── __cld_polys (M1)
                    │       └── 复用：pair_vec_div（Z/m[x] 精确除法），__symmetric_mod
                    ├── __build_cld_matrix (M2)
                    │       └── 复用：__symmetric_mod
                    ├── __lll_reduce (M3)       ← 新实现（Cohen §2.6）
                    ├── __extract_candidates (M4)
                    └── 因子验证
                            └── 复用：__upoly_symmetric_mod, __upoly_primitive,
                                       pair_vec_div（ZZ 试除）
```

新增代码全部位于 `polynomial_factorize_univar.hh`，无新文件，无新依赖。

---

## 4. 关键设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 数据约束方案 | **CLD**（非 Newton 迹） | 界更紧（不随列索引指数增长）；有多项式时间证明（Belabas 2009）；Maple/FLINT 均采用 |
| CLD 计算方式 | `C_i = (f/h_i) · h_i'`（精确除法） | 避免多项式幂级数求逆；h_i \| f 保证精确可除 |
| 矩阵初始化 | `2^U_exp · I_r`，U_exp = bit_length(max(r,20)) | 来自 FLINT 实现；与阈值 B = (r+1)·2^(2·U_exp) 配合，保证 LLL 能区分短向量 |
| 列喂入顺序 | 内螺旋（两端交替向中间） | 来自 FLINT；尽早暴露两端因子结构，比顺序喂入收敛更快 |
| 因子提取 | 幺模变换矩阵 U 列分组 | 比直接读{0,1}更稳健；短向量不一定有整{0,1}分量 |
| LLL 实现 | 自写整数 LLL（Cohen §2.6，δ=3/4）| 无外部依赖；r ≤ 100，精度足够 |
| 初始列数 J₀ | 30（3r > N+1 时）或 10 | 来自 FLINT 实测参数 |
| 列上界 J_max | (deg_f + 1) / 2 | 来自 FLINT；超过此范围信息无增益 |
| Zassenhaus 阈值 | r ≤ 10 | 小 r 已有剪枝，Zassenhaus 足够快；阈值可调 |
| Fallback 策略 | J_target > J_max 时回退 Zassenhaus | 确保正确性；理论上 m 充足时 van Hoeij 必然终止 |
| 矩阵表示 | `vector<vector<ZZ>>` | 实现简单；r ≤ 100，无性能瓶颈 |
| 代码位置 | 全部在 polynomial_factorize_univar.hh | 与 __factor_recombine 同文件，无跨文件接口问题 |

---

## 5. 不在本次范围内（P1a）

- Newton 迹方案（已由 CLD 替代，理论描述保留在调研报告 §2.2 方案 A）
- 浮点加速 LLL（v1 用纯整数）
- ZZ 矩阵专用内存优化

---

## Part B：线性 Hensel 提升（P1b）

> 阶段：**架构完成**
> 对应调研文档：`docs/research/linear-hensel-research.md`
> 目标文件：`clpoly/polynomial_factorize_univar.hh`

---

## 6. 核心流程（P1b）

P1b 将现有的"提升到 Mignotte 精度 → 重组"两阶段流程，改为"**n 因子直接**线性步提升 + LLL 检测"的**交织循环**，实现早期因子检测（early termination）。不使用二叉树结构，所有 r 个因子直接在一个扁平循环中同步提升。

```
__factor_squarefree_primitive_ZZ
    │
    │  [P1b 修改前]
    │  __hensel_lift(f, sel.factors, sel.prime)         ← 二次提升到 Mignotte 精度
    │  __factor_recombine(f, lifted, modulus)            ← P1a LLL 重组
    │
    └── [P1b 修改后]
        __linear_hensel_lift_with_lll(f, sel.factors, sel.prime)  [M3]
            │
            ├── 预处理：
            │     lc_f = lc(f)；对 factors[0] 乘以 lc_f mod p（首项系数处理）
            │     s[0..r-1] = __linear_bezout_chain(factors_adj, p)  [M0]  ← 一次性预计算
            │     h[0..r-1] = Zp_to_ZZ(factors_adj[i])                    ← 初始提升基
            │     a₀ = __heuristic_starting_precision(f, r, p)  [M4]
            │     a_mig = ⌈log_p(2·lc_f·B_Mig(f))⌉
            │     p_a = p^1
            │
            ├── 初始线性提升（a₀ 步）：
            │     for step = 1 to a₀:
            │         __hensel_step_linear_nfactor(h, s, f, p, p_a)   [M1]
            │         p_a *= p
            │
            ├── 初始化 van Hoeij LLL 状态（复用 P1a 逻辑）：
            │     active = {0..r-1}，f_star = f，result = []
            │     U_exp = bit_length(max(r, 20))
            │     B = (r+1) · 2^(2·U_exp)
            │     M = 2^U_exp · I_r（r×r）
            │     J_cur = 0，J_target = J₀
            │
            └── 主交织循环（直到 |active| ≤ 1）：
                    │
                    ├── [P1a M1] __cld_polys(f_star, active_h_sym, p_a)
                    ├── [P1a M2] __build_cld_matrix(M, cld, J_cur, J_target, p_a)
                    │               → J_new = 0 时倍增 J_target，重入循环顶
                    ├── [P1a M3] __lll_reduce(M, U, B) → short_rows
                    ├── [P1a M4] __extract_candidates(short_rows, U, r') → candidates
                    │
                    ├── 因子验证（同 P1a M5 逻辑）：
                    │     对每个候选子集 S：
                    │         g_trial = lc(f_star) · ∏_{i∈S} h[i]  (mod p_a)，对称约化，本原化
                    │         pair_vec_div(f_star, g_trial) → 试除
                    │         成功：result.push(g_trial)，f_star /= g_trial
                    │               缩减 active，重置 M/J_cur/J_target  ← 早期因子发现
                    │
                    └── 本轮无进展：
                            若 J_target 可增大（< J_max）：J_target *= 2
                            否则（列已耗尽）：                    ← 增量精度提升
                                for step = 1 to k（默认 k = 5）：
                                    __hensel_step_linear_nfactor(h, s, f, p, p_a)  [M1]
                                    p_a *= p；a_cur += 1
                                若 a_cur > a_mig：  ← 安全网（理论上不应触发）
                                    对剩余 active 调用 __zassenhaus_recombine
                                    return result + zassenhaus_result
```

---

## 7. 模块划分（P1b）

### M0 `__linear_bezout_chain`

**功能规约**

```
函数名称：__linear_bezout_chain

功能描述：对 r 个两两互素的模因子 h_1,...,h_r ∈ Z/p[x]，
          计算 Bézout 系数链 s_1,...,s_r 满足：
              Σᵢ sᵢ · ∏_{j≠i} hⱼ ≡ 1  (mod p)
          且 deg sᵢ < deg hᵢ。
          此系数链在整个线性提升过程中保持不变（固定 mod p）。

前置条件（Requires）：
  - h_1,...,h_r ∈ Z/p[x]，两两互素，首一
  - r ≥ 2

后置条件（Ensures）：
  - 返回 s_1,...,s_r ∈ Z/p[x]（以 upolynomial_<ZZ> 形式，系数 ∈ [0,p)）
  - Σᵢ sᵢ · ∏_{j≠i} hⱼ ≡ 1  (mod p)
  - deg sᵢ < deg hᵢ，∀i

计算方法：
  1. P = ∏ᵢ hᵢ  (mod p)                     ← 一次全乘积
  2. 对每个 i：
       Pᵢ = P / hᵢ  (mod p)                 ← Z/p[x] 精确除法
       extGCD(Pᵢ, hᵢ) → aᵢ 使 aᵢ·Pᵢ ≡ 1  (mod hᵢ, mod p)
       sᵢ = aᵢ mod hᵢ（确保 deg sᵢ < deg hᵢ）

副作用：无。
```

**接口规约**

```
接口：__linear_hensel_lift_with_lll → __linear_bezout_chain

输入数据：
  - factors: const vector<upolynomial_<Zp>>&（lc_f 调整后的 Zp 因子）
  - p: uint32_t

输出数据：
  - vector<upolynomial_<ZZ>>（r 个 Bézout 系数，系数 ∈ [0,p)，值语义与 Zp 一致）

协议约定：
  - 调用方保证因子两两互素（由 __select_prime 保证）
  - 被调用方只调用一次（M3 预处理阶段）；输出在整个提升过程中不变
```

---

### M1 `__hensel_step_linear_nfactor`

**功能规约**

```
函数名称：__hensel_step_linear_nfactor

功能描述：对 r 个提升因子执行一步 n 因子直接线性 Hensel 提升。
          利用预计算的 Bézout 链（固定 mod p），将
              ∏ hᵢ ≡ f  (mod p^a)
          同步提升至
              ∏ hᵢ' ≡ f  (mod p^(a+1))
          无需二叉树结构，所有因子在同一层处理。

前置条件（Requires）：
  - f ∈ Z[x]（ZZ 系数，完整无截断）
  - h[i] ∈ Z[x]，系数 mod p^a，∈ [0, p^a)，∀i
  - s[i] ∈ Z[x]，系数 ∈ [0, p)（来自 M0，固定不变），∀i
  - Σᵢ s[i] · ∏_{j≠i} h[j] ≡ 1 (mod p)（M0 后置保证）
  - ∏ h[i] ≡ f (mod p^a)
  - p_a = p^a

后置条件（Ensures）：
  - 每个 h[i] 原地更新为 h'[i]，系数 mod p^(a+1)，∈ [0, p^(a+1))
  - ∏ h'[i] ≡ f (mod p^(a+1))
  - s[i] 全部不变

算法：
  e = (f - ∏ h[i]) / p_a  mod p           ← 误差，在 Z/p 中
  for i = 1..r:
      σ[i] = s[i] · e  mod h[i]  (mod p)  ← 独立 r 个 mod-p 除法
      h'[i] = h[i] + p_a · σ[i]  mod p^(a+1)

副作用：原地修改 h[0..r-1]；不修改 s[0..r-1]。
```

**正确性（一步验证）**

```
∏ h'[i] ≡ ∏ h[i] + p_a · Σᵢ σ[i] · ∏_{j≠i} h[j]  (mod p^(a+1))
        ≡ (f - p_a·e·p_a⁰) + p_a · e · Σᵢ s[i]·∏_{j≠i}h[j]  (mod p^(a+1))
        ≡ f - p_a·e + p_a·e  = f  (mod p^(a+1))   ✓
（交叉项 p^(2a) 在 2a ≥ a+1 ∀a≥1 时消去）
（Diophantine 解：Σᵢ σ[i]·∏_{j≠i}h[j] ≡ e (mod p) 由 M0 Bézout 链保证）
```

**接口规约**

```
接口：__linear_hensel_lift_with_lll → __hensel_step_linear_nfactor

输入数据：
  - h: vector<upolynomial_<ZZ>>&（r 个提升因子，原地更新）
  - s: const vector<upolynomial_<ZZ>>&（Bézout 链，只读）
  - f: const upolynomial_<ZZ>&
  - p: uint32_t，p_a: const ZZ&

输出数据：无（原地修改 h[0..r-1]）

协议约定：
  - 调用方在每次调用后将 p_a 乘以 p（调用方维护精度）
  - 被调用方不修改 s
```

---

### M3 `__linear_hensel_lift_with_lll`

**功能规约**

```
函数名称：__linear_hensel_lift_with_lll

功能描述：P1b 主控函数。对 f ∈ Z[x] 的模因子列表，通过线性 Hensel 提升
          与 van Hoeij LLL（P1a M1-M4）交织，实现早期因子检测。
          一旦 LLL 找到所有因子立即返回，无需提升到完整 Mignotte 精度。

前置条件（Requires）：
  - f 本原，lc(f) > 0，deg(f) ≥ 2
  - factors 是 f mod p 的不可约因子（互素，首一）
  - |factors| ≥ 2

后置条件（Ensures）：
  - 返回 f 的所有不可约因子（ZZ 系数，本原，lc > 0）
  - 结果等价于原 __hensel_lift + __factor_recombine 的组合，但通常更早完成

副作用：无（不修改输入）。
```

**接口规约**

```
接口：__factor_squarefree_primitive_ZZ → __linear_hensel_lift_with_lll

输入数据：
  - f: const upolynomial_<ZZ>&
  - factors: const vector<upolynomial_<Zp>>&（mod p 模因子）
  - p: uint32_t

输出数据：
  - vector<upolynomial_<ZZ>>（不可约因子列表）

协议约定：
  - 调用方保证 factors 是 f mod p 的无平方因式分解
  - 被调用方负责首项系数处理（lc_f 分配），调用方无需预处理
  - 被调用方保证即使早期终止，结果因子的乘积等于 f
```

---

### M4 `__heuristic_starting_precision`

**功能规约**

```
函数名称：__heuristic_starting_precision

功能描述：估计线性 Hensel 提升的初始精度 a₀（步数，使首次 LLL 检测有合理的
          收敛概率），采用 FLINT 启发式公式，取其与 Mignotte 精度的较小值。

前置条件（Requires）：
  - f ∈ Z[x]（用于计算 f->length）
  - r ≥ 2（模因子数）
  - p（素数）

后置条件（Ensures）：
  - 返回 a₀：整数，1 ≤ a₀ ≤ a_mig
  - 公式：a_heuristic = ⌈(2.5·r + min_b) · ln(2) / ln(p) + ln(N+1) / (2·ln(p))⌉
    其中 min_b ≈ log₂(p)（首个模因子最小系数的 bit 数），N = f->length - 1
  - a₀ = min(a_mig, a_heuristic)

副作用：无。
```

---

### 接口修改：`__factor_squarefree_primitive_ZZ`

```
修改内容：
  [P1b 前]
    auto [lifted, modulus] = __hensel_lift(f, sel.factors, sel.prime);
    return __factor_recombine(f, lifted, modulus);

  [P1b 后]
    return __linear_hensel_lift_with_lll(f, sel.factors, sel.prime);

原 __hensel_lift 和 __factor_recombine 保留（M3 内部安全网的 fallback 路径可调用）。
```

---

## 8. 模块依赖与调用关系（P1b）

```
__factor_squarefree_primitive_ZZ
    └── __linear_hensel_lift_with_lll (M3)
            ├── __linear_bezout_chain (M0)        ← 预处理，调用一次
            │       └── 复用：polynomial_GCD（Z/p 扩展 GCD），Zp 多项式除法
            ├── __heuristic_starting_precision (M4)
            │       └── 无外部依赖（纯计算）
            ├── __hensel_step_linear_nfactor (M1)  ← 每提升步调用
            │       └── 复用：__upoly_mod_coeff, __upoly_divmod_mod,
            │                  ZZ::fdiv_q, ZZ::fdiv_r
            ├── [P1a M1] __cld_polys              ← 直接复用 P1a 模块
            ├── [P1a M2] __build_cld_matrix        ← 直接复用 P1a 模块
            ├── [P1a M3] __lll_reduce              ← 直接复用 P1a 模块
            ├── [P1a M4] __extract_candidates      ← 直接复用 P1a 模块
            ├── 因子验证逻辑（内联于 M3 主循环）
            │       └── 复用：__upoly_symmetric_mod, __upoly_primitive, pair_vec_div
            └── [安全网] __zassenhaus_recombine    ← fallback，现有代码，不修改
```

新增代码全部位于 `polynomial_factorize_univar.hh`，无新文件，无新依赖。
不再依赖 `__hensel_node` 结构和 `__hensel_tree_build`（二叉树完全废弃于 P1b 提升路径）。

---

## 9. 关键设计决策（P1b）

| 决策 | 选择 | 理由 |
|------|------|------|
| Bézout 系数精度 | **固定 mod p**（不随 a 更新） | 线性提升核心：每步运算在 Z/p 上（单精度），消除 Bézout 更新代价；详见调研报告 §2.3 |
| 提升结构 | **n 因子直接方法**（无二叉树） | 结构更简单（平铺循环）；误差一次计算；与 Monagan 2022 路线一致；Bézout 链预计算代价 O(r·n²/p)，仅一次 |
| 外层精度管理 | **倍增触发**（每次 LLL 未收敛时 J_target *= 2，列耗尽后增量 k 步） | FLINT 实证：精度倍增优于单步线性；k=5 平衡检测频率与提升代价 |
| 初始精度 a₀ | **FLINT 启发式**（调研报告 §3.2），取与 Mignotte 精度的较小值 | 避免从 Mignotte 精度全量提升；uni-70 实测约节省 40% 提升步数 |
| LLL 检测时机 | 每次列耗尽后的增量提升 k 步后 | 与 P1a 主循环逻辑保持一致，代码复用最大化 |
| P1a 模块复用 | M3 主循环**直接内联调用** P1a M1-M4 | P1b 吸收 P1a M5（`__vanhoeij_recombine`）的主循环逻辑；P1a M5 在 M3 完成后不再独立调用 |
| 安全网 | a > a_mig 时 fallback Zassenhaus | 保证正确性；理论上 Mignotte 精度充足时不应触发 |
| 代码位置 | 全部在 `polynomial_factorize_univar.hh` | 与 P1a 一致，无跨文件接口 |

---

## 10. 不在本次范围内（P1b）

- **[TODO v2] Monagan 2022 矩阵优化**：将全部 D 步提升的 Bézout 修正序列组织为多项式矩阵乘法，把 n 因子直接方法的总代价从 O(rn²D) 降至 O(rnD²/r) = O(nD²)（"立方代价"）。与早期终止叠加后收益有限，留作后续增强。
- 浮点辅助精度估计（初版保守使用 FLINT 启发式整数公式）
- FLINT 精确的 `_fmpz_poly_hensel_continue_lift` 实现（初版用简单循环）

---

## P2（多变量优化）

P2 是独立的多变量因式分解性能优化，目标为 `polynomial_factorize_wang.hh`。

详见 [`docs/design/mtshl/architecture.md`](../mtshl/architecture.md)。
