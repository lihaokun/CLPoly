# P1b 启发式精度提前终止：细化设计文档

> 阶段：细化（workflow.md §2.3）
> 基础文档：`docs/design/vanhoeij/architecture.md` Part B
> 目标文件：`clpoly/polynomial_factorize_univar.hh`
> 方案：**两阶段精度提前终止**（先试 a_h，再 fallback a_mig）

---

## 0. 旧设计失败原因分析

旧设计（"n 因子直接线性步交织循环"）在实现后产生性能回退，根本原因有两点：

### 原因 1：在不充足精度下运行 LLL → O(n⁴) 迭代

旧设计每提升 k=5 步就运行一次 LLL。当 Hensel 精度 a < a_h 时，CLD 多项式的系数界远大于提升精度，LLL 格基无法快速收敛——每次 LLL 调用需 O(n⁴) 次 Gram-Schmidt 迭代，而非 O(n)。结果：运行了数十次低效 LLL，累积代价远超直接提升到 a_h 再运行一次 LLL 的代价。

代码中已有注释（§8.3 前）：
> "Van Hoeij LLL 需要在满精度（a_mig 或 a_h）下运行才能高效收敛；低精度下运行 LLL 会导致格基规约需 O(n⁴) 次迭代而非 O(n)，性能急剧下降。"

### 原因 2：旧 M3 初始化 J_target = J₀（跳过对角 LLL）

旧设计在 `__linear_hensel_lift_with_lll` 中初始化 `J_target = J0`，跳过了对角 LLL（即先以纯 2^U_exp·I_r 对角矩阵运行一次 LLL 的步骤）。

当前 `__vanhoeij_recombine` 已正确实现对角 LLL（第 1162 行 `J_target = 0`），旧设计重新实现主循环时丢失了这一优化。

### 正确结论

**LLL 必须在充足精度（a_h 或 a_mig）下运行**，且每个精度级别只运行一次。不应交织多次低精度 LLL。

---

## 1. 正确设计概述

用最少的代码变更实现"启发式精度提前终止"：

```
__lll_factorize(f, factors, p)
    │
    ├── 计算 a_h = __heuristic_starting_precision(f, r, p)   ← 已存在 §6.7
    │
    ├── [阶段 1] __hensel_lift(f, factors, p, a_h)  → lifted_h, m_h
    │       __vanhoeij_recombine(f, lifted_h, m_h, &zass_triggered)
    │       if !zass_triggered → return（快速路径，大多数情况）
    │
    └── [阶段 2] __hensel_lift(f, factors, p)        → lifted_mig, m_mig（全 Mignotte）
                __vanhoeij_recombine(f, lifted_mig, m_mig)  → return
```

**阶段 1 快速路径** 适用条件：
- `a_heuristic < a_mig`（启发式精度更小，才有节省）
- LLL 在 a_h 精度下收敛（不触发 Zassenhaus fallback）

典型用例 `(x-1)...(x-70)`（p=71, r=70）：
- a_h = 31，a_mig = 66
- 因子 x-i 的系数最大为 70，而 71^31 ≈ 10^56 >> 70，LLL 在 a_h=31 必然收敛
- 阶段 1 即返回，Hensel 步数从 7 步（p^64 ≥ Mignotte）降至 5 步（p^32 ≥ p^31）

**阶段 2 fallback** 仅在 a_h 精度不足时触发（罕见）。

---

## 2. 代码变更

所有变更均在 `clpoly/polynomial_factorize_univar.hh`，共 2 处修改：

### 2.1 `__hensel_lift`：增加 `a_target` 参数

**位置**：§6.6，line 494

**签名变更**：

```cpp
// 改前：
inline std::pair<std::vector<upolynomial_<ZZ>>, ZZ>
__hensel_lift(
    const upolynomial_<ZZ>& f,
    const std::vector<upolynomial_<Zp>>& factors,
    uint32_t p)

// 改后：
inline std::pair<std::vector<upolynomial_<ZZ>>, ZZ>
__hensel_lift(
    const upolynomial_<ZZ>& f,
    const std::vector<upolynomial_<Zp>>& factors,
    uint32_t p,
    int a_target = 0)   // 0 = 全 Mignotte；> 0 = 停在 p^a >= p^a_target
```

**循环停止条件变更**：

```cpp
// 改前：
ZZ B = __mignotte_bound(f);
ZZ lc_f = f.front().second;
if (lc_f < ZZ(0)) lc_f = -lc_f;
ZZ target = ZZ(2) * lc_f * B;
// ...
ZZ m(p);
while (m <= target)
{
    __hensel_lift_recursive(nodes, 0, f, m);
    m = m * m;
}

// 改后：
ZZ target;
if (a_target == 0) {
    ZZ B = __mignotte_bound(f);
    ZZ lc_f_abs = f.front().second;
    if (lc_f_abs < ZZ(0)) lc_f_abs = -lc_f_abs;
    target = ZZ(2) * lc_f_abs * B;  // Mignotte 停止条件
} else {
    target = ZZ(1);
    for (int i = 0; i < a_target; ++i)
        target *= ZZ(p);              // 目标 = p^a_target
    target -= ZZ(1);                 // while(m <= p^a_target-1) ↔ while(m < p^a_target)
}
// ... (lc_f、因子调整代码保持不变，但 lc_f 已在 target 分支内计算)
// 二次提升循环不变：
ZZ m(p);
while (m <= target)
{
    __hensel_lift_recursive(nodes, 0, f, m);
    m = m * m;
}
```

> **注**：`a_target > 0` 时停止条件为 `m <= p^a_target - 1`，即 `m < p^a_target`。循环退出后 m ≥ p^a_target（通常为某个 p^(2^k)，二次提升不能精确命中 a_target）。这是正确的：LLL 只需 m > Mignotte 界，精度略超 p^a_target 不影响正确性。

### 2.2 `__vanhoeij_recombine`：无变更

`__vanhoeij_recombine` 签名保持不变。

> **实现说明**：初始设计考虑增加 `bool* zassenhaus_triggered` 输出参数来检测 Phase 1 失败，
> 但最终实现采用了更简洁的方案（见 2.3），无需修改本函数。

### 2.3 `__lll_factorize`：两阶段逻辑

**位置**：§8.3，line 1428

**完整替换**：

```cpp
// 改前：
inline std::vector<upolynomial_<ZZ>>
__lll_factorize(
    const upolynomial_<ZZ>&              f,
    const std::vector<upolynomial_<Zp>>& factors,
    uint32_t                             p)
{
    auto [lifted, m] = __hensel_lift(f, factors, p);
    return __vanhoeij_recombine(f, lifted, m);
}

// 改后：
inline std::vector<upolynomial_<ZZ>>
__lll_factorize(
    const upolynomial_<ZZ>&              f,
    const std::vector<upolynomial_<Zp>>& factors,
    uint32_t                             p)
{
    int r   = (int)factors.size();
    int a_h = __heuristic_starting_precision(f, r, p);

    // Phase 1：提升到启发式精度 a_h
    auto [lifted_h, m_h] = __hensel_lift(f, factors, p, a_h);
    auto result = __vanhoeij_recombine(f, lifted_h, m_h);

    // Phase 1 返回单个因子且 r > 1，说明精度不足（Mignotte 条件未满足）。
    // Phase 2：提升到完整 Mignotte 精度，保证正确性。
    if ((int)result.size() == 1 && r > 1) {
        auto [lifted_mig, m_mig] = __hensel_lift(f, factors, p);
        return __vanhoeij_recombine(f, lifted_mig, m_mig);
    }
    return result;
}
```

**Phase 2 触发条件说明**：

原设计用 `zassenhaus_triggered` 检测 LLL 是否触发了 Zassenhaus fallback。
但实测发现更可靠的触发条件是 `result.size() == 1 && r > 1`：

- `result.size() == 1` 说明所有子集试除均失败，f 被误判为不可约
- `r > 1` 确认存在多个模因子，f 应当可分解
- **根本原因**：二次提升会"超调"——以 `a_target=k` 停止时，实际精度为某个 `p^(2^j) > p^k`，
  但该精度仍可能小于真正的 Mignotte 界（`a_mig`）。例如：`lc_f=10000, r=2, a_h=3, a_mig=5`
  时，二次提升到 `p^4 < Mignotte`，Zassenhaus 全部失败返回 `{f}`（而非触发 zassenhaus_triggered）。
- `zassenhaus_triggered` 方案下，LLL→Zassenhaus→{f} 也会设置标志，但等价于检查 `result.size()==1`，
  更直接简洁。

---

## 3. 不变部分

| 模块 | 变更 |
|------|------|
| `__heuristic_starting_precision` | **不变**（§6.7，已实现） |
| `__vanhoeij_recombine` 主循环逻辑 | **不变**（J_target=0 对角 LLL 等已正确实现） |
| `__hensel_lift` 二次提升结构 | **不变**（仅添加停止条件参数） |
| `__hensel_tree_build` / `__hensel_step` | **不变**（二叉树结构完整保留） |
| `__zassenhaus_recombine` | **不变** |
| `__factor_recombine` | **不变**（现有调用兼容默认参数） |
| `__factor_squarefree_primitive_ZZ` | **不变** |

**不实现**：`__linear_bezout_chain`、`__hensel_step_linear_nfactor`、`__linear_hensel_lift_with_lll`（旧 M0/M1/M3）。这些是"完整线性交织"路线的组件，已确认其会导致性能回退，不在本次范围内。

---

## 4. 预期性能

| 用例 | 当前 (a_mig) | P1b 后 (a_h) | Hensel 步数 | 预期 total |
|------|------------|-------------|-----------|----------|
| uni-70 | 16ms | ~8ms | 5 vs 7 步 | ~13ms |
| W(15) | <1ms | <1ms | a_h=a_mig | 不变 |
| W(20) | ~1ms | ~1ms | a_h=a_mig | 不变 |
| deg25-5fac | ~2ms | ~2ms | a_h≈a_mig | 不变 |

> uni-70 中 a_h=31 < a_mig=66，Hensel 从 7 步（到 p^64）减至 5 步（到 p^32）。
> 其他用例 a_h ≥ a_mig（启发式无收益），行为与当前完全一致。

---

## 5. 测试策略

### T1. 回归测试

```
bash test/run_all_tests.sh
```

全部通过（277 个测试）。a_target=0 路径（默认参数）与改前完全一致。

### T2. 性能验证

```bash
g++ -O3 -DNDEBUG -std=c++20 -I. test/profile_factorize.cc lib/libclpoly.a -lgmpxx -lgmp -o /tmp/pf && /tmp/pf
```

确认 `uni 70 fac` 的 `total` 列时间降低（实测 ~13ms，原约 21ms）。
注：profile 中 `hensel` 列仍测量全 Mignotte 提升（非 P1b 路径），该列无意义；请看 `total` 列。

### T3. 压力测试

```
make stress
```

确认 uni-70 因子分解结果正确性（尤其 a_h 路径不触发 Zassenhaus 时的因子数 = 70）。

### T4. 交叉验证（可选）

```
make crosscheck
```

确认 FLINT/NTL 交叉验证无新失败。

---

## 6. 关键设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 提升方式 | 保留**二次提升**（非线性逐步） | 二次提升总步数更少（O(log a) vs O(a)）；线性交织方案已确认回退 |
| 阶段 2 重新提升 | **从头提升**（非从 a_h 树继续） | 实现简单；Phase 2 罕见，重新提升代价可接受；避免维护中间 Hensel 树 |
| J_target 初始化 | 由 `__vanhoeij_recombine` 内部管理（J_target=0） | 对角 LLL 已在现有代码中正确实现，不在调用方重新实现 |
| 参数类型 | `a_target: int`（步数，非 ZZ 精度值） | 与 `__heuristic_starting_precision` 返回类型一致；避免 ZZ 大整数 |
| Phase 2 触发条件 | `result.size()==1 && r>1` | 直接检测失败结果（精度不足 → 全部试除失败 → 返回 {f}）；比 `zassenhaus_triggered` 更简洁且覆盖相同场景 |
| 阶段 2 后不再检查 | 不需要 | m_mig > 2·lc(f)·B_Mig 保证 LLL 必然收敛（理论保证）；Zassenhaus fallback 仍作为安全网保留 |
