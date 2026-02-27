# M6：MTSHL-d 论文对齐 & 优化

> 创建日期：2026-02-27
> 前置：M5 完成（`__mtshl_lift` 已接入 `__wang_core`，全量测试通过）
> 关联文档：`detailed-design.md`（M1–M5）、`architecture.md`
> 论文：CASC 2016、CASC 2018、ICMS 2018（`docs/research/papers/`）

---

## 背景

M5 集成测试中发现两个 bug 并修复。修复后与论文全文对照，确认了 3 处
实现与设计文档/论文的差异。三处差异均不影响正确性（全量测试通过），
但偏离了论文的最优设计。M6 目标：逐项对齐论文，消除技术债。

---

## 差异 1：Taylor 循环终止条件

### 现状

```cpp
// __mtshl_step_j, wang.hh L859-862
int64_t D_j = degree(aj, xj);
for (int k = 1; k <= D_j; ++k)
```

### 论文

**CASC 2018 Algorithm 4, Line 3**（多因子，r ≥ 2）：
```
for k = 1, 2, 3, ... while error ≠ 0 and ∑ᵢ deg(fj,i, xj) < deg(aj, xj) do
```

双终止条件：
1. `error = 0` → 提前退出（分解已精确完成）
2. `∑ deg(fj,i, xj) ≥ deg(aj, xj)` → 度数上界（保证终止）

**ICMS 2018 Algorithm 2, Line 2**（2 因子特化）：
```
for i = 1, 2, ... while df + dg < da do
```

### 分析

| 属性 | 当前实现 `k <= D_j` | 论文双条件 |
|------|-------------------|-----------|
| 正确性 | 正确 | 正确 |
| 终止保证 | 有（固定上界） | 有（度数上界） |
| 提前退出 | 无（但 `if (ck.empty()) continue` 使空迭代代价极低） | 有（`error = 0` 时立即退出） |
| 等价性 | `k <= D_j` 是度数上界的充分条件：当 k = D_j 时 ∑deg ≥ D_j | 完全匹配 |

**M5 bug 回顾**：原代码 `while (!error.empty())` 只有条件 1（无度数上界），
r ≥ 3 时交叉项导致 error 永不为零 → 死循环。修复为 `k <= D_j` 即条件 2 的等价形式。

### M6 改动

**优先级**：低（性能收益可忽略，空迭代代价 ≈ 0）

可选改进：在循环体内加提前退出，完全匹配论文双条件：

```cpp
int64_t D_j = degree(aj, xj);
for (int k = 1; k <= D_j; ++k)
{
    // ... 现有代码 ...

    // 论文条件 1：提前退出
    if (error.empty()) break;
}
```

---

## 差异 2：Mignotte 界检查

### 现状

```cpp
// __wang_core, wang.hh L2161-2164
// MTSHL prime: 使用最大 31-bit 素数
// 不做 Mignotte 界检查——依靠 trial division 验证因子正确性。
static constexpr uint32_t mtshl_p = 2147483629;
```

设计文档 M5 规定了 `__mtshl_mignotte_check`（已删除），要求 `p > 2·B(f_scaled)`。

### 论文

**CASC 2018 不要求 Mignotte 检查**。

论文使用 p-adic 提升（Algorithm 5）：选单个机器素数 p，迭代提升到 p^l > 2B。
不需要 `p > 2B`——多次迭代自然处理大系数。

### 分析

| 方案 | Mignotte 检查 | 系数恢复 | 大系数处理 |
|------|-------------|---------|-----------|
| 设计文档（已废弃） | `p > 2B` 检查 | 单次对称约化 | p 不满足时 continue |
| 当前实现 | 无 | 单次对称约化 | trial division 兜底 |
| 论文 Alg.5 | 无 | p-adic 多次迭代 | 自然处理 |

**M5 bug 回顾**：LC 缩放后 `f_scaled` 系数可达 10⁷+，Mignotte 界 `2^(d+1)·√N·‖f‖∞`
经常超过 31-bit 素数上限 (~2×10⁹) → 几乎所有求值点被拒绝 → 死循环。

**当前策略的正确性**：
- 若 p 不足以恢复真实系数 → 对称约化给出错误系数 → trial division 拒绝
  → `__wang_core` 换求值点重试
- 理论风险：错误因子恰好整除原多项式（概率极低，实测未发生）

**根本解法**：实现论文的 p-adic 提升（差异 3 的一部分），消除对单次对称约化的依赖。

### M6 改动

**优先级**：中（随差异 3 一并解决）

1. 实现 p-adic 提升后，Mignotte 检查不再需要
2. 更新设计文档 M5 节：删除 `__mtshl_mignotte_check` 相关内容
3. 更新错误处理策略汇总表：删除 `__mtshl_mignotte_check` 行

---

## 差异 3：素数选择 & 系数恢复

### 现状

```cpp
// __wang_core, wang.hh L2165
static constexpr uint32_t mtshl_p = 2147483629;  // 最大 31-bit 素数
```

- 固定 31-bit 素数，独立于 `__wang_core` 的原始素数
- 系数恢复：单次对称约化 `mods(coeff, p)`
- `scaled_factors` 是精确 Z[x1] 多项式 → 约化到任意素数均有效

### 论文

**CASC 2018 Algorithm 5**：
```
1. 选 p 为 (m+1)-bit 随机机器素数（64 位机器上 m=62 → 63-bit 素数）
2. 在 Zp 中完整运行 MTSHL-d（Algorithm 4）→ 得 f1,...,fr mod p
3. p-adic 提升：
   for k = 1, 2, ..., l（l = ⌈log_p(2B)⌉）:
     dk = (a - ∏ fi) / p^k       ← 差商，Zp 运算
     SparseInt 解 MDP → σk,i      ← 强支撑假设 (Theorem 1)
     fi += σk,i · p^k
4. 最终 ∏ fi = a over Z
```

**Theorem 1（p-adic 支撑保持）**：
```
Pr(Supp{uj} ⊆ Supp{uj-1} ∀ 1 ≤ j ≤ l) > 1 − tl/((π(2^(m+1)) − π(2^m)) − l)
```
对 m=62, l=5, t=500：概率 > 1 − 10⁻¹⁵（几乎确定）。

### 分析

| 属性 | 当前实现 | 论文 Algorithm 5 |
|------|---------|-----------------|
| 素数大小 | 31-bit (uint32_t) | 63-bit (uint64_t) |
| 素数选择 | 固定 2147483629 | 随机 (m+1)-bit |
| 系数恢复 | 单次 mods(c, p)，p/2 ≈ 10⁹ | p-adic 迭代，上限 p^l |
| 可恢复系数范围 | |c| < 10⁹ | |c| < p^l/2 ≈ 10^(19l) |
| Theorem 1 保证 | 固定素数无概率保证 | 随机素数有高概率保证 |
| 失败后果 | trial division 拒绝，换求值点 | 极低概率失败 |

**核心限制**：31-bit 素数 + 单次约化只能恢复绝对值 < 10⁹ 的系数。
对 LC 缩放后的大系数多项式（如 Swinnerton-Dyer），系数可达 10¹⁰+，
此时当前实现必然失败（靠 trial division 兜底不断重试，直到换到系数更小的求值点）。

63-bit 素数 + p-adic 提升从根本上解决此问题。

### M6 改动

**优先级**：高（63-bit 素数是 M6 核心任务）

分两步实施：

#### M6a：63-bit 素数支持（Zp64 类）

扩展 `Zp` 类支持 64-bit 模数，或新增 `Zp64` 类。

**方案对比**：

| 方案 | 改动范围 | 优势 | 劣势 |
|------|---------|------|------|
| A: 扩展 `Zp` 为 `uint64_t` | 全局（所有 Zp 用户） | 统一类型，无模板分支 | 单变量因式分解不需要 64-bit，乘法开销增加 |
| B: 新增 `Zp64` 类 | 仅 MTSHL 模块 | 不影响现有代码 | 需在 MTSHL 中参数化类型或硬编码 Zp64 |
| C: 模板化 `Zp<W>` | 全局（模板化） | 灵活，编译期选择 | 改动大，需重构所有 Zp 使用点 |

**建议方案 A**：将 `Zp` 内部从 `uint32_t` 扩展到 `uint64_t`。

关键改动点：
1. **`Zp` 存储**：`int _number` → `int64_t _number`，`int _p` → `uint64_t _p`
2. **乘法**：`a * b % p` 需 128-bit 中间值 → 用 `__int128` 或 `_mulmod`
3. **逆元**：扩展 GCD 改为 64-bit 版本
4. **接口兼容**：`number()` 返回类型从 `int` → `int64_t`；所有 `(int, int)` 构造需兼容

影响评估：
- `Zp` 在 MTSHL、`__select_prime`、GCD 模块中使用
- 单变量因式分解的 `__berlekamp_small_p` 等仅用小素数（p < 10⁶），uint64_t 无性能差异
- 乘法热路径需基准测试确认（`__int128` 在 x86-64 上编译为 `mul` + `div`，约 4ns）

#### M6b：p-adic 提升（Algorithm 5 外循环）

在 63-bit 素数基础上，为 `__mtshl_lift` 增加 p-adic 提升外循环。

**改动函数**：`__mtshl_lift`

```cpp
// 当前：单次 Zp 计算 + 对称约化
auto mv_factors = __mtshl_lift(f_scaled, scaled_factors, lc_targets, eval, x1, p);

// M6b 后：单次 Zp 计算 + p-adic 迭代
// __mtshl_lift 内部：
// 1. Zp 中运行 __mtshl_step_j (j=2..n) → F[i] mod p
// 2. error = (f_scaled_zz - ∏ sym_mod(F[i], p)) / p
// 3. while error ≠ 0:
//      c = error mod p                      ← Zp 运算
//      SparseInt 解 MDP → σk,i              ← 复用 __mtshl_sparse_int
//      F[i] += σk,i · p^k                   ← 大整数加法
//      error = (error - Σ σk,i·bi) / p
// 4. return sym_mod(F[i], p^l)
```

注意：p-adic 提升的 SparseInt 使用**强支撑假设**（Theorem 1）：
`Supp(σk,i) ⊆ Supp(σk-1,i)`，即每一步的支撑只会缩小。
63-bit 随机素数下概率 > 1 − 10⁻¹⁵。

**提升界**：
```
B = Mignotte_bound(f_scaled)
l = ⌈log_p(2B)⌉
```
对 63-bit 素数 p ≈ 2⁶³，大多数情况 l = 1（单次约化即可），
极端大系数情况 l = 2–3（如 Swinnerton-Dyer deg-128：系数约 10³⁸ → l ≈ 2）。

---

## 改动汇总

### M6a：63-bit 素数（核心前置）

| 文件 | 改动 |
|------|------|
| `clpoly/number.hh` | `Zp` 类：uint32→uint64，乘法用 `__int128`，逆元 64-bit |
| `clpoly/polynomial_factorize_wang.hh` | `mtshl_p` 改为 63-bit 素数 |
| 全局 | `Zp` 接口返回类型 int→int64_t，调用点适配 |

验收：
- 全量测试通过（`bash test/run_all_tests.sh`）
- `make bench-all` 无性能退化（基准测试对比）

### M6b：p-adic 提升（可选，M6a 后）

| 文件 | 改动 |
|------|------|
| `clpoly/polynomial_factorize_wang.hh` | `__mtshl_lift` 增加 p-adic 外循环 |
| `clpoly/polynomial_factorize_wang.hh` | 新增 `__mtshl_padic_step`（单步 p-adic 修正） |
| `clpoly/polynomial_factorize_wang.hh` | `__wang_core`：删除 `mtshl_p` 固定值，改为随机 63-bit 素数 |
| `test/test_mtshl_m6.cc` | p-adic 提升单元测试（大系数恢复验证） |

验收：
- 大系数测试用例通过（Swinnerton-Dyer、SymPy f_1 等）
- crosscheck 通过
- 性能对比：bench-all 与 M5 基线持平或更优

### 差异 1 修复（可选，低优先级）

| 文件 | 改动 |
|------|------|
| `clpoly/polynomial_factorize_wang.hh` | `__mtshl_step_j` Taylor 循环加 `if (error.empty()) break;` |

---

## 设计文档同步更新

M6 实施后需同步更新以下设计文档内容：

1. **`detailed-design.md` M5 节**：
   - 删除 `__mtshl_mignotte_check` 相关内容（改动 3、Mignotte 界计算）
   - 更新 `__wang_core` 调用示例（63-bit 素数）
   - 更新错误处理策略汇总表（删除 mignotte 行）

2. **`detailed-design.md` 新增 M6 节引用**：指向本文档

3. **`architecture.md`**：更新素数策略描述

---

## 附录：M5 后其他优化方向（参考）

以下优化不在 M6 范围内，记录供后续参考。

### HenselLift1 + 双变量归约（ICMS 2018）

r=2 时的 j-th step 优化：将多变量问题分层为双变量 HenselLift1 调用。
对 r=2、低维（j=3,4）情形约快 2–3x。对 bivar-70（r=70）无显著影响。

详见 ICMS 2018 Algorithm 2。实现思路：
1. 新增 `__hensellift1(aj_bivar, f0, g0, αj, p)`
2. 修改 `__mtshl_step_j`（j≥3, r=2）：重构为 HenselLift1 + 插值路径
3. 并行化

### Bernardin 增量误差更新（ICMS 2018 §2）

Taylor 循环中避免每步重算完整乘积 `∏F[i]`，改用已存储的 Taylor 系数卷积：
- r=2：`ck = a^(k)(αj)/k! - Σ σi·τ_{k-i}`（ICMS 2018 原文，仅适用于 2 因子）
- r>2：逐对卷积 `Σ_{k₁+...+kᵣ=k} ∏σ_{l,kₗ}`（泛化，论文中无描述）

当前 `__mtshl_step_j` 每步全量重算 `error = aj - ∏F[i]`（与 CASC 2018 一致）。
Bernardin 消除这个冗余：只计算第 k 阶系数，不重算全部阶。

**状态**：待实验论证。需要：
1. 构造 r=2、高变量数（n≥6）、高次数（d≥40）的 benchmark
2. 对比全量重算 vs 卷积的实际耗时
3. 评估 r>2 泛化的收益（σk 逐步收缩，卷积项快速变少）
4. 决定是否值得增加代码复杂度

### ~~双变量基底（Bivariate Base）~~ [已废弃]

> **勘误**：经核实论文原文（MC 2019, ICMS 2018 等），Monagan & Tuncer 全部论文
> 均使用**单变量基底 Zp[x₁]**，不存在"双变量基底"算法。
> CLPoly 的单变量基底与 Maple MTSHL 完全一致，此优化方向已取消。
