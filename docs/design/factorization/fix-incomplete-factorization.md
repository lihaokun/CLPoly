# 不完全因式分解 Bug：根因分析与修复方案

## 1. 现象

### 1.1 复现用例

**Case 1**: `f = 10707840x^12 - 4358016x^11 - ... + 4287360x`

| | 因子数 | deg6 部分 |
|---|---|---|
| CLPoly | 4 | `2288x^6+1216x^5-559x^4-428x^3+1819x^2-1276`（作为不可约返回） |
| FLINT | 5 | `(44x^3-13x^2+29)(52x^3+43x^2-44)` |

乘积验证均正确，问题是 CLPoly 未能拆分 deg6。

### 1.2 频率

压力测试 580 个用例中 5 个失败（~0.9%），均为同一模式：两个不可约因子被合并为一个可约多项式返回。

### 1.3 有趣的对照

**单独重新分解 deg6 因子时，CLPoly 正确拆分为两个 deg3！** 这排除了算法本身的理论缺陷，指向特定的实现 Bug。

## 2. 因式分解流水线追踪

### 2.1 整体流程

```
factorize(f)
 ├─ 提取 content = 3, f_prim = f/3 (lc = 3569280)
 ├─ squarefreefactorize → f_prim 本身无平方 (deg 12)
 └─ __factor_squarefree_primitive_ZZ(f_prim)
     ├─ __select_prime → p=17, r=8
     └─ __lll_factorize(f_prim, 8 mod-p factors, p=17)
         ├─ Phase 1: Hensel lift to a_h=7, m=17^8≈7×10^9
         │   └─ __vanhoeij_recombine → 4 factors (BUG)
         └─ Phase 2: triggered (4 < 8 && 7 < 14)
             ├─ Hensel lift to a_mig=14, m=17^16≈4.9×10^19
             └─ __vanhoeij_recombine → 4 factors (SAME BUG)
```

### 2.2 素数选择

| 素数 | 状态 | r | 模因子度数分布 |
|------|------|---|---------------|
| 2,3,5 | skip (lc ≡ 0) | — | — |
| 7 | skip (非无平方) | — | — |
| 11,13 | skip (lc ≡ 0) | — | — |
| **17** | **选中 (tried #1, 最小 r)** | **8** | [1,1,1,1,1,2,2,3] |
| 19 | tried #2 | 9 | [1,1,1,1,1,1,1,2,3] |
| 23 | tried #3 | 9 | [1,1,1,1,1,1,1,2,3] |

### 2.3 真因子的模 17 分解

| 真因子 | mod 17 分解 | 模因子数 |
|--------|-----------|---------|
| x | 1^1 | 1 |
| 65x²-61x-80 | 2^1 (不可约二次) | 1 |
| 24x³+9x+14 | 3^1 (不可约三次) | 1 |
| 44x³-13x²+29 | 1^1 + 2^1 (一次+二次) | 2 |
| 52x³+43x²-44 | 1^3 (三个一次) | 3 |

**总模因子数 = 1+1+1+2+3 = 8 = r**。每个真因子至少贡献一个模因子，理论上足以拆分。

### 2.4 Hensel 提升的因子

```
lifted[0] deg=1 lc=3569280   ← lc(f_prim) baked in
lifted[1] deg=1 lc=1          ← monic
lifted[2] deg=1 lc=1
lifted[3] deg=1 lc=1
lifted[4] deg=1 lc=1
lifted[5] deg=2 lc=1
lifted[6] deg=2 lc=1
lifted[7] deg=3 lc=1
```

## 3. Bug 根因分析

### 3.1 van Hoeij 第一轮（J_target=0 单因子预提取）

8 个单因子候选中：
- `cand{4}` (deg1) → **x** ✓
- `cand{5}` (deg2) → **65x²-61x-80** ✓
- `cand{7}` (deg3) → **24x³+9x+14** ✓
- 其余 5 个 FAIL（正确，它们是 h3/h4 的部分因子）

提取后：`active = {0,1,2,3,6}`（5 个模因子），`f_star = h3×h4 = 2288x^6+...`

**关键状态变化**：
- `lc(f_star)` 从 3569280 变为 **2288**
- 但 `lifted[0]` 仍然携带 `lc = 3569280`（原始 lc baked in，未更新）

### 3.2 van Hoeij 第二轮（5 个模因子分组）

**J_target=0**: 5 个单因子候选全部 FAIL（正确——h3 和 h4 各需要多个模因子组合）。

**J_target=6**: 添加 6 个 CLD 列，11×11 LLL 规约。**结果：仍然是 5 个单因子组**。LLL 未能正确分组。

**原因**：CLD 计算使用 `f_star / active_lifted[k]` mod m。但 `active_lifted[0] = lifted[0]` 具有 `lc = 3569280`，而 `f_star` 的 `lc = 2288`。虽然模除法在 Z_m[x] 上是合法的，但 CLD 系数中混入了 `3569280/2288 = 1560` 的缩放因子，破坏了 CLD 的等价类结构。

### 3.3 Zassenhaus 回退——**Bug 爆发点**

LLL 失败后（J_target=12 > J_max=6），回退到 Zassenhaus：
```cpp
auto zass = __zassenhaus_recombine(f_star, active_lifted, m);
```

Zassenhaus 对所有 s=1,2 的子集进行剪枝检查。**实测：所有 15 个子集全部被错误剪枝！**

#### Bug 1: lc 剪枝（含 index 0 的子集）

```
lc_sq = 2288² = 5,234,944
lc_prod = 3,569,280   (lifted[0] 的 lc)
5,234,944 % 3,569,280 = 1,665,664 ≠ 0  → 错误拒绝
```

**原因**：代码假设 `lifted[0]` 携带 `lc(f_star)` 但实际携带 `lc(f_original)`。

```cpp
// 错误逻辑：
bool subset_has_lc = false;
for (size_t i : S_idx)
    if (i == 0) { subset_has_lc = true; break; }
ZZ lc_mult = subset_has_lc ? ZZ(1) : lc_fstar;
```

当 Zassenhaus 被 van Hoeij 回退调用时：
- `lifted` = `active_lifted`（active 子集的因子）
- `active_lifted[0]` = 原始 `lifted[0]`，其 lc = `lc(f_original)` ≠ `lc(f_star)`
- 代码用 `lc_mult = 1` 是基于 "index 0 已有 lc baked in" 的假设，但此 lc 是错误的

#### Bug 2: 常数项剪枝（不含 index 0 的子集）

```
fstar_const = 2288 × (-1276) = -2,919,488
c_prod = 2288 × const(lifted[k1]) × const(lifted[k2]) mod m  （巨大的模值）
```

对于正确的 h3 因子对 {k, 4}：
- `const(monic h3) ≡ 29/44 (mod m) ≈ 大整数`
- `c_prod = 2288 × (29 × 44⁻¹ mod m) ≈ O(m)` — 远大于 `fstar_const`
- `-2,919,488 % c_prod ≠ 0` → 错误拒绝

**原因**：模因子的常数项是以 `lc(f_original)` 为基准的模表示，而剪枝用的是 `lc(f_star)`。两者之间有 `lc(f_original)/lc(f_star) = 3569280/2288 = 1560` 倍的缩放差异，导致常数项是大模值而非小整数。

### 3.4 为何 Phase 2 同样失败

Phase 2 将 Hensel 精度从 33 bits 提升到 66 bits，但 **Bug 与精度无关**——它源于 lc 不匹配。Phase 2 的 van Hoeij + Zassenhaus 回退展现完全相同的失败模式。

### 3.5 为何单独重新分解 deg6 时成功

当 `factorize(2288x^6+...)` 被调用时：
1. `__select_prime` 为 deg6 选择 p=19（不同的素数）
2. r=3：`[1, 2, 3]`（DDF 分出 deg1+deg2+deg3）
3. `lifted[0]` 的 lc = 2288 = `lc(f)` ← **匹配！**
4. J_target=0 提取 deg3 因子 (h3)
5. 剩余 2 个模因子，`f_star = 52x³+43x²-44`，`lc(f_star) = 52`

此时 `lifted[0]` 也被消耗（因为 r=3 中 deg3 的因子是 lifted[2]），active 重新构建后不再包含 lc-baked 因子，因此 Zassenhaus 回退时使用正确的 `lc_mult = lc(f_star) = 52`。

### 3.6 根因总结

**核心 Bug**：van Hoeij 提取部分真因子后，`f_star` 的 lc 改变，但 `lifted[0]` 的 baked-in lc 仍为原始值。这导致：

1. **Zassenhaus lc 剪枝**：含 index 0 的子集被错误拒绝（lc_mult = 1 对应错误的 lc）
2. **Zassenhaus 常数项剪枝**：模因子常数项以原始 lc 为基准，与新 `f_star` 的 lc 不匹配，导致常数项是大模值
3. **Van Hoeij CLD**：CLD 系数中混入缩放因子，破坏 LLL 分组
4. **Van Hoeij 试除**：`g_trial = lc(f_star) × lifted[0]` 产生双 lc，但 pp() 可以恢复（此处影响较小）

**触发条件**：
- `lc(f)` 有多个素因子（使 lc 在提取因子后改变）
- 至少一个真因子对应多个模因子（需要组合重组）
- `lifted[0]`（lc-baked 因子）恰好对应需要组合重组的因子之一

## 4. 修复方案

### 4.1 方案选择

| 方案 | 描述 | 优点 | 缺点 |
|------|------|------|------|
| A: 归一化 lifted[0] | 提取因子后，将 lifted[0] 除以其 lc 再乘以 lc(f_star)（mod m） | 精确修复 | 需要跟踪 lc-baked index |
| **B: 消除 lc 特殊处理** | 不再将 lc bake into lifted[0]；所有 lifted 因子保持 monic，始终用 lc(f_star) 做前缀 | 最简洁，消除根因 | 需修改 Hensel 提升入口 |
| C: 修复剪枝 | 保留当前 lc 分配，修复剪枝条件以正确处理 lc 不匹配 | 最小改动 | 治标不治本 |

**推荐方案 B**：从根源消除问题。

### 4.2 方案 B 详细设计

#### 修改点 1: `__hensel_lift`（line 562-568）

**当前**：将 `lc(f)` bake into `factors_adj[0]`
```cpp
std::vector<upolynomial_<Zp>> factors_adj = factors;
Zp lc_mod_p(f.front().second, p);
for (auto& term : factors_adj[0])
    term.second *= lc_mod_p;
factors_adj[0].normalization();
```

**修改**：不 bake lc，所有因子保持 monic。Hensel 提升的目标函数改为 `f/lc(f)` (monic)。
```cpp
// 所有因子保持 monic，不 bake lc
// Hensel 目标：f_monic = f / lc(f)
upolynomial_<ZZ> f_monic = f;
ZZ lc_f = f.front().second;
ZZ lc_inv;
// 注：不能在 ZZ 上精确除——改在 Hensel step 中处理
```

**但有一个问题**：标准 Hensel 提升要求 `f ≡ lc(f) × H₁ × ... × Hᵣ (mod p)`，其中 Hᵢ monic。如果不 bake lc，则 `∏ Hᵢ = f/lc(f) (mod p)`，这不是整数除法。

**更好的方案 B'**：保留 lc bake-in，但在 van Hoeij 提取因子后**重新归一化** `active_lifted`。

#### 修改点 1 (B'): `__vanhoeij_recombine` 提取因子后的重建

在 `found_any` 分支中，提取因子并重建 active 后，如果原始 `lifted[0]` 仍在 active 中，将其归一化为 monic，并记录 lc 已被移除：

```cpp
if (found_any)
{
    // 移除已消耗因子
    for (int j = (int)active.size() - 1; j >= 0; --j)
        if (consumed[j]) active.erase(active.begin() + j);

    // 归一化：如果 lc-baked factor (original index 0) 仍在 active 中，
    // 将其变为 monic（mod m），消除 lc 不匹配
    // ... (见下方完整代码)

    // 重建格矩阵参数
    // ...
}
```

#### 修改点 2: 追踪 lc-baked index

在 `__vanhoeij_recombine` 开头添加 `lc_baked_idx` 变量：

```cpp
int lc_baked_idx = 0;  // lifted[0] 初始携带 lc(f)
```

提取因子后，如果 `lc_baked_idx` 对应的因子被消耗：
```cpp
lc_baked_idx = -1;  // 无 lc-baked 因子，全部 monic
```

如果 `lc_baked_idx` 仍在 active 中：
```cpp
// 将 lifted[lc_baked_idx] 变为 monic: 除以其 lc (mod m)
ZZ lc_inv;
bool ok = ZZ::invert(lc_inv, active_lifted[local_idx].front().second, m);
assert(ok);
for (auto& term : active_lifted[local_idx])
    term.second = (term.second * lc_inv) % m;
lc_baked_idx = -1;  // 现在所有因子都是 monic
```

#### 修改点 3: 统一试除逻辑

所有因子均为 monic 后，试除始终使用 `lc(f_star)` 前缀：
```cpp
g_trial = lc(f_star) × ∏ active_lifted[k]
pp_g = pp(g_trial)
```

这与现有 van Hoeij 代码一致（它已经始终用 `lc(f_star)` 前缀）。

#### 修改点 4: 修复 Zassenhaus 的 lc 特殊处理

当 Zassenhaus 被 van Hoeij 回退调用时，所有 active_lifted 已 monic 化，不再需要 index 0 特殊处理。

但 Zassenhaus 也被 `__factor_recombine` 直接调用（r ≤ 10 时），此时 `lifted[0]` 确实有 lc baked in。因此需要区分两种调用场景。

**最简方案**：在 van Hoeij 回退 Zassenhaus 前，保证 active_lifted 全部 monic，然后 Zassenhaus 始终使用 `lc(f_star)` 前缀（移除 index 0 特殊处理）。

对于 `__factor_recombine` 直接调用的场景：也将 lifted[0] 变 monic（在入口处处理）。这样 Zassenhaus 代码统一为 "所有因子 monic + lc_mult = lc(f)" 的简单模式。

### 4.3 改动清单

| 文件 | 函数 | 改动 |
|------|------|------|
| `polynomial_factorize_univar.hh` | `__hensel_lift` | lifted[0] 提取后立即 monic 化 |
| `polynomial_factorize_univar.hh` | `__vanhoeij_recombine` | 始终使用 lc(f_star) 前缀（已是） |
| `polynomial_factorize_univar.hh` | `__zassenhaus_recombine` | 移除 index 0 特殊处理，始终 lc_mult = lc_fstar |
| `polynomial_factorize_univar.hh` | `__factor_recombine` | 入口处 monic 化 lifted[0] |

### 4.4 正确性论证

设 f ∈ Z[x] 本原无平方，f ≡ lc(f) × H₁ × ... × Hᵣ (mod m)，Hᵢ monic。

若 monic 化 lifted[0]（除以 lc(f)），则：
∏ᵢ lifted[i] = ∏ᵢ Hᵢ ≡ f / lc(f) (mod m)

对任何真因子 h 对应的子集 S：
lc(f) × ∏_{i∈S} Hᵢ ≡ lc(f) × h / lc(h) (mod m)

pp(lc(f) × h / lc(h)) = pp((lc(f)/lc(h)) × h) = h（因为 h 本原，lc(f)/lc(h) 是整数常数）

当 m > 2 × |lc(f)| × B_mig(f) 时，symmetric mod 正确恢复系数。✓

提取因子后 f_star 的处理：
- f_star = pp(f / h)，lc(f_star) 可能 ≠ lc(f)
- 剩余 monic 因子的乘积 ≡ f_star / lc(f_star) (mod m)
- 试除：lc(f_star) × ∏ monic → pp 得到真因子 ✓
- 剪枝：lc_prod = lc(f_star)（因为所有因子 lc=1），lc_sq % lc_prod = 0 ✓
- 常数项：c_prod = lc(f_star) × ∏ const(monic_Hᵢ)，这是小整数（Mignotte 保证） ✓

## 5. 补充：van Hoeij LLL 未能分组的分析

即使修复了 Zassenhaus 剪枝 Bug，LLL 也未能正确分组 5 个模因子为 {2, 3}。原因：

1. **CLD 缩放污染**：`active_lifted[0]` 的 lc=3569280 引入 1560 倍缩放，CLD 系数不在统一基准上
2. **J_max 限制**：N=12（原始 deg），J_max=6。修复后 J_max 应基于 deg(f_star)=6 重新计算

LLL 分组失败不是致命问题（有 Zassenhaus 回退），但修复 CLD 缩放和 J_max 计算可以提升 LLL 效率。

## 6. 测试要求

1. **Case 1 & Case 2 回归**：修复后 Case 1/2 应返回 5 个因子
2. **全量测试**：`bash test/run_all_tests.sh`
3. **交叉验证**：`make crosscheck`
4. **压力测试**：重新运行 `test_factorize_stress` 验证 cross_fail = 0
5. **性能回归**：`make bench-clpoly` 确保无性能退化
