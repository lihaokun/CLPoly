# P2a 架构文档：素数试验次数优化

> 状态：待确认
> 调研依据：`docs/research/factorize-profiling.md` §3.3

---

## 1. 核心流程

当前 `__select_prime` 流程：

```
for tried = 0..max_tries(30-50):
    p = next_prime()
    if lc(f) mod p == 0: continue
    if deg(f mod p) != deg(f): continue
    if GCD(fp, fp') != 1: continue   ← 通过则 tried++
    DDF(fp) + EDF(fp)                ← 每次都执行完整分解
    if r_p < best_count: update best
return best
```

目标流程（改动极小）：

```
for tried = 0..MAX_TRIES(3):         ← 只改这一个常数
    p = next_prime()
    if lc(f) mod p == 0: continue
    if deg(f mod p) != deg(f): continue
    if GCD(fp, fp') != 1: continue   ← 通过则 tried++
    DDF(fp) + EDF(fp)
    if r_p < best_count: update best
return best
```

## 2. 模块划分

改动范围：`polynomial_factorize_univar.hh §8.1`，仅修改 `__select_prime` 中一个常量。

```
受影响模块       改动内容
─────────────────────────────
__select_prime   max_tries: 30 → 3
__lll_factorize  无改动（调用方）
```

## 3. 关键设计决策

### 决策 1：试验次数选 3（参照 FLINT）

FLINT 使用 `for (i = 0; i < 3; i++)`，经过工程实践验证。3 次的意义：
- 有足够概率避开"坏"素数（lc 整除、次数下降）
- 能在少数候选中选出 r 较小的那个
- 相较于 1 次：多 2 次 DDF+EDF 开销，但 r 略优（对 LLL 影响极小）

### 决策 2：保留 min-r 策略（不改为"取第一个"）

虽然 LLL 不需要最小 r，保留 3 次中取最优是零成本的（3 次都要做 DDF+EDF）。去掉 min-r 不减少开销，可以保留。

### 决策 3：`max_tries` 扩展逻辑一并移除

当前代码：
```cpp
size_t max_tries = 30;
if (best_count > (size_t)(deg_f / 2) && tried == 30)
    max_tries = 50;  // 扩展到 50
```

此逻辑（"r 大时多试"）是 Zassenhaus 时代的遗产。切换到 3 次后，这个分支直接删除。

## 4. 接口规约

```
接口：__select_prime → __lll_factorize

输入数据：upolynomial_<ZZ> f（无平方本原多项式）
输出数据：__prime_selection_result { prime, factors, irreducible }
          — 与当前完全相同，接口不变
协议约定：
  - 调用方保证：f 无平方、本原、deg ≥ 2
  - 被调用方保证：返回一个使 f mod prime 无平方的 prime
```

## 5. 预期收益与风险

| 指标 | 当前 | 改后 | 变化 |
|------|------|------|------|
| DDF+EDF 调用次数 | 30-50 | 3 | **-90%** |
| sel_prime 时间 | ~4ms（W15） | ~0.4ms | **-90%** |
| 总时间 | ~4.5ms（W15） | ~0.9ms | **-80%** |
| r 的质量 | 全局最优 | 3 次中最优 | 略降，对 LLL 无实质影响 |

**风险**：无。r 略大不影响 LLL 正确性，FLINT 在生产中已验证 3 次足够。
