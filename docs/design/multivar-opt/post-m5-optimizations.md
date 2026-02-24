# MTSHL-d 后续优化方向（M5 之后）

> 创建日期：2026-02-25
> 关联文档：`architecture.md`（MTSHL-d 架构，M1–M5 范围）
> 依据：Maple 2019 完整框架 vs 当前架构的三项差距（见下）

---

## 背景

架构文档（M1–M5）实现了 Maple 2019 的 MTSHL-d 核心算法
（CASC 2016 + CASC 2018 + JSC 2020），正确性完整覆盖。

Maple 2019 的完整框架还包含三项额外特性，当前均不在 M1–M5 范围内，
记录如下，供后续性能调优参考。

---

## 优化 1：HenselLift1 + 双变量归约（ICMS 2018）

**论文**：Monagan & Tuncer, ICMS 2018, LNCS 10931:378–387

**Maple 2019 的做法**（r=2 情形）：

```
提升 xj（j≥3）时，不将 x2,...,xj-1 同等处理，而是分层：
  1. 将 x3,...,xj-1 稀疏求值到 s 个双变量像 aj(x1,x2,βk,xj)
  2. 对每个双变量像，将 x2 稠密求值到 deg(aj,x2)+1 个点
  3. 对每个单变量像，调用 HenselLift1（可并行）
  4. 稠密插值恢复 x2 系数
  5. 稀疏插值恢复 x3,...,xj-1 系数
```

**HenselLift1**（Algorithm 2, ICMS 2018）：
- 输入：Zp[x1,xj] 双变量问题（已求值到单变量）
- 输出：fj, gj ∈ Zp[x1,xj]
- 优化：σk 在 x1 = 0,...,d1 处的值缓存复用，Σ 代价从 O(id1²) → O(d1²)
- 复杂度：O(d1²·dj + d1·dj²)
- 可并行：每个 HenselLift1 实例独立，用 Cilk spawn 并行

**当前架构的做法**（CASC 2016/CASC 2018 风格）：
- x2,...,xj-1 均用 θ-array 几何点一次性求值到单变量
- 逐点解 Zp[x1] MDP，再 Vandermonde 恢复所有辅助变量

**差距影响**：
- 正确性：无影响
- 性能：对 r=2、低维（j=3,4）情形，ICMS 2018 实现约快 2–3x
- 对主要基准（bivar-70，r=70）：瓶颈在多因子 MDP，HenselLift1 无显著影响

**实现思路**（M6 参考）：
1. 新增 `__hensellift1(aj_bivar, f0, g0, αj, p)` → Zp[x1,xj] 双变量提升
2. 修改 `__mtshl_step_j`（j≥3）：将 j=3 情形重构为 HenselLift1 + x2 稠密插值路径
3. 并行化：`cilk_spawn` 或 C++17 `std::for_each` 并行调用

---

## 优化 2：双变量基底（Bivariate Base）

**出处**：multivar-mtshl-research.md §7.4；MC 2019 隐含；Maple 2019 实现

**问题背景**：

当前架构使用单变量基底：
```
factorize(f|_{x2=α2,...,xn=αn})  → r 个 Zp[x1] 因子
```

Maple 2019 使用双变量基底：
```
factorize(f|_{x3=α3,...,xn=αn})  → r 个 Zp[x1,x2] 因子（保留 x2）
```

**双变量基底的优势**：
- Zp[x1,x2] 因式分解比 Zp[x1] 更能保留因子结构
- 好的求值点下，每个 Zp[x1,x2] 因子以高概率对应唯一真因子（s=1）
- s=1 时 Zassenhaus 重组代价降为 O(r)（无指数爆炸）

**当前影响**：
- 对"所有因子在单变量层面分离"的基准（如 bivar-70）：无影响，
  `__select_eval_point` 已保证单变量层面的因子分离性
- 潜在问题：某些多项式单变量层面 s>1，Zassenhaus 重组指数退化；
  双变量基底可消除此退化

**实现思路**（M6 参考）：
1. 修改 `__wang_core`：将求值从 `f|_{x2=α2,...}` 改为 `f|_{x3=α3,...}` 保留 x2
2. 改 `factorize` 基底为双变量因式分解（调用现有 `__factor_bivar` 或新增）
3. 改 `__wang_leading_coeff` 以处理双变量 LC 分配
4. 改 Zassenhaus 重组以处理双变量因子的重组（取代单变量 Zassenhaus）

注：此优化改动范围较大（涉及 `__wang_core` 主流程），建议在 M5 验收且
    基准测试显示重组步骤是瓶颈时再启动。

---

## 优化 3：p-adic 大系数提升（CASC 2018 Algorithm 5）

**论文**：Monagan & Tuncer, CASC 2018 §4, Algorithm 5

**问题背景**：

当前架构（决策 1）使用单机器素数 p（uint32_t），要求 `p > 2·Mignotte_bound(f_scaled)`。
若 Mignotte 界超过约 2×10⁹（uint32_t 极限），则无合适的机器素数，当前直接报告失败。

**CASC 2018 Algorithm 5 方案**：
```
1. 选 63-bit 机器素数 p（uint64_t），在 Zp 中完整运行 MTSHL-d
   → 得 f1,...,fr mod p（机器精度因子）
2. for k = 1, 2, ..., l（l = ⌈log_p(2B)⌉）:
   dk = (a - Π fi) / p^k       ← p-adic 差商，Zp 运算
   SparseInt 解 Σ σk,i·bi = dk  ← 强支撑假设：Supp(σk,i) ⊆ Supp(σk-1,i) w.h.p.
   fi += σk,i·p^k
3. 最终 Π fi = a over Z         ← p-adic 恢复
```

**性能优势**（CASC 2018 §4）：
- 99%+ 运算在 64-bit 模运算下完成，无大数开销
- 对大系数多项式（Swinnerton-Dyer 等），比传统 p^l 方案快 3–10x

**实现条件**：
- 需将 `Zp` 类从 uint32_t 扩展到 uint64_t（Zp 乘法用 `__int128`）
- `__mtshl_lift` 增加 p-adic 外循环（l 步修正）
- `SparseInt` 保持不变（仍在 Zp 中工作）

**实现思路**（M6 参考）：
1. 新增 `Zp64` 类（uint64_t 系数，乘法用 `__int128`）
2. 修改 `__mtshl_lift`：增加 p-adic 提升外循环
3. 修改 `__wang_core`：选 p 时优先选 63-bit 机器素数

---

## 优先级建议

| 优化 | 影响范围 | 实现复杂度 | 建议时机 |
|------|---------|-----------|---------|
| HenselLift1（优化 1） | r=2 情形，低维性能 | 中（新增双变量提升单元）| M5 后，如 r=2 用例性能不达标 |
| p-adic（优化 3） | 大系数多项式 | 中（新增 Zp64 + 外循环） | M5 后，如 Mignotte 界溢出报告多 |
| 双变量基底（优化 2） | 重组步骤（边界情形） | 大（改动 wang_core 主流程）| M5 后，如 stress 测试中重组成为瓶颈 |
