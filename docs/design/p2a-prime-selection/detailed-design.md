# P2a 细化文档：素数试验次数优化

> 状态：已实现
> 依据：`docs/design/p2a-prime-selection/architecture.md`

---

## 1. 修改范围

单一函数：`__select_prime`，位于 `clpoly/polynomial_factorize_univar.hh §8.1`，第 1323-1403 行。

---

## 2. 函数级改动

### `__select_prime`

**签名**（不变）：
```cpp
inline __prime_selection_result __select_prime(const upolynomial_<ZZ>& f)
```

**改动 1**：`max_tries` 初始值 30 → 3

```cpp
// 改前：
size_t max_tries = 30;

// 改后：
size_t max_tries = 3;
```

**改动 2**：删除 50 次扩展逻辑（Zassenhaus 时代遗产，LLL 路径无需）

```cpp
// 删除以下两行：
if (best_count > (size_t)(deg_f / 2) && tried == 30)
    max_tries = 50;
```

**改动后完整函数体逻辑**（伪代码，供审查）：

```
for idx = 0, tried = 0; tried < 3; ++idx:
    p = boost::math::prime(idx)
    if lc(f) mod p == 0: continue
    fp = polynomial_mod(f, p)
    if deg(fp) != deg(f): continue
    if deg(GCD(fp, fp')) > 0: continue
    ++tried
    __upoly_make_monic(fp)
    ddf = __ddf_Zp(fp)
    irr_factors = flatten(EDF each ddf group)
    if size == 1: return immediately (irreducible)
    if size < best_count: update best
return best
```

**不变**：函数签名、返回类型 `__prime_selection_result`、内部数据结构、min-r 选择逻辑。

---

## 3. 复用点

无新函数，无新文件。仅修改已有常量和条件分支。

---

## 4. 测试要求

1. **功能测试**：运行全量测试 `bash test/run_all_tests.sh`，确认 278/278 通过
2. **性能验证**：运行 `_build/release/bin/profile_factorize`，确认 `sel_prime` 时间降至原来的约 1/10
3. **压力测试**：`make stress`，确认 uni-70 用例正确性
