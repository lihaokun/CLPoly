# GCD 优化细化设计

> 状态：已确认（2026-03-01）
> 架构依据：`docs/design/gcd-opt/architecture.md`（已确认 2026-02-28）

---

## 实现顺序

```
M0.1 → M0.2 → M0.4 → M0.3 → M0.5 → P0
```

每个模块独立提交，逐步验证。

---

## M0.1: ZZ 接口扩展

**文件**：`clpoly/number/ZZ.hh`

### 1. GMP 64 位 helper（ZZ 类定义之前，namespace 内）

`mpz_set_si` / `mpz_get_si` / `mpz_set_ui` / `mpz_get_ui` 在 LLP64（Windows）上 `long` 为 32 位，会截断 `int64_t` / `uint64_t`。提供 4 个 helper，ZZ 内部及后续 `next_prime_64` 等统一使用：

```cpp
// ---- GMP 64-bit helpers ----

inline void _mpz_set_s64(mpz_ptr z, int64_t v) {
    if constexpr (sizeof(long) == 8) {
        mpz_set_si(z, v);
    } else {
        uint64_t abs_v = (v >= 0) ? static_cast<uint64_t>(v)
                                  : static_cast<uint64_t>(-(v + 1)) + 1;
        mpz_import(z, 1, 1, sizeof(uint64_t), 0, 0, &abs_v);
        if (v < 0) mpz_neg(z, z);
    }
}

inline int64_t _mpz_get_s64(mpz_srcptr z) {
    if constexpr (sizeof(long) == 8) {
        return mpz_get_si(z);
    } else {
        uint64_t abs_val = 0;
        mpz_t tmp;
        mpz_init(tmp);
        mpz_abs(tmp, z);
        mpz_export(&abs_val, NULL, 1, sizeof(uint64_t), 0, 0, tmp);
        mpz_clear(tmp);
        // 避免 -INT64_MIN 的有符号溢出 UB，与 _mpz_set_s64 保持一致
        if (mpz_sgn(z) >= 0)
            return static_cast<int64_t>(abs_val);
        else
            return static_cast<int64_t>(-(abs_val - 1) - 1);
    }
}

inline void _mpz_set_u64(mpz_ptr z, uint64_t v) {
    if constexpr (sizeof(unsigned long) == 8) {
        mpz_set_ui(z, v);
    } else {
        mpz_import(z, 1, 1, sizeof(uint64_t), 0, 0, &v);
    }
}

inline uint64_t _mpz_get_u64(mpz_srcptr z) {
    if constexpr (sizeof(unsigned long) == 8) {
        return mpz_get_ui(z);
    } else {
        uint64_t result = 0;  // mpz 为 0 时 export 不写入
        mpz_export(&result, NULL, 1, sizeof(uint64_t), 0, 0, z);
        return result;
    }
}
```

**ZZ 内部替换**：详见 §3。

### 2. ZZ 接口扩展

在 `public:` 区域的 state query 段（line 240 `// ---- state query ----` 之后）新增：

```cpp
// ---- accessors for internal representation ----
bool       is_small()  const { return _is_small(); }
int64_t    get_val()   const { return _val; }      // 仅 is_small()==true 时有意义
mpz_srcptr get_mpz_v() const { return _mpz; }      // is_small()==true 时返回 nullptr
```

**命名说明**：`get_mpz_v()` 而非 `get_mpz()`，因为 GMP 宏 `mpz_t` 展开后可能与 `get_mpz` 产生歧义。若实测无冲突可改为 `get_mpz()`。

在构造函数区域（line 150 `ZZ(const std::string&...)` 之后）新增：

```cpp
ZZ(mpz_srcptr z) : _val(0), _mpz(nullptr) {
    if (_fits_si(z)) {
        _val = _mpz_get_s64(z);
    } else {
        _mpz = _mpz_new();
        mpz_set(_mpz, z);
    }
}
```

### 3. ZZ 内部 LLP64 整体梳理

ZZ.hh 内部所有 `mpz_set_si` / `mpz_get_si` / `mpz_init_set_si` 在 LLP64 上均存在 `int64_t` → `long`（32 位）截断风险。需**全量替换**为 `_mpz_set_s64` / `_mpz_get_s64`：

| 模式 | 替换为 | 数量 | 说明 |
|------|--------|------|------|
| `mpz_set_si(z, v)` | `_mpz_set_s64(z, v)` | ~11 | 含 `_promote()`、算术、`INT64_MIN/MAX` |
| `mpz_init_set_si(z, v)` | `mpz_init(z); _mpz_set_s64(z, v)` | ~29 | 算术运算中的临时变量提升 |
| `mpz_get_si(z)` | `_mpz_get_s64(z)` | ~5 | 含 `_demote_if_small()`、`get_si()`、`operator==` |

**特别注意**：`_fits_si()` 内部 line 66 `mpz_set_si(max_val, INT64_MAX)` 和 line 76 `mpz_set_si(min_val, INT64_MIN)` 也必须替换为 `_mpz_set_s64`。`_fits_si` 判断的是 `int64_t` 范围，但它自身用 `mpz_set_si` 设置上下界——LLP64 上被截断为 32 位，导致 `_fits_si` 判断范围错误，进而所有依赖 `_fits_si` 的逻辑（`_demote_if_small`、`operator==`、构造函数等）全部失效。修复后 line 79-81 的 string fallback 可一并删除（`_mpz_set_s64` 已处理 LLP64）。

**测试**：

```cpp
// test/test_zz_accessors.cc
void test_is_small() {
    ZZ a(42);
    assert(a.is_small() && a.get_val() == 42);
    assert(a.get_mpz_v() == nullptr);

    ZZ b("99999999999999999999");  // > INT64_MAX
    assert(!b.is_small());
    assert(b.get_mpz_v() != nullptr);
    assert(mpz_cmp_si(b.get_mpz_v(), 0) > 0);
}

void test_construct_from_mpz() {
    mpz_t tmp;
    mpz_init_set_si(tmp, 123);
    ZZ a(static_cast<mpz_srcptr>(tmp));
    assert(a.is_small() && a == 123);
    mpz_clear(tmp);

    mpz_t big;
    mpz_init_set_str(big, "99999999999999999999", 10);
    ZZ b(static_cast<mpz_srcptr>(big));
    assert(!b.is_small() && b == ZZ("99999999999999999999"));
    mpz_clear(big);
}

void test_mpz_s64_helpers() {
    // 正常值
    mpz_t z;
    mpz_init(z);
    _mpz_set_s64(z, 42);
    assert(_mpz_get_s64(z) == 42);

    // 负值
    _mpz_set_s64(z, -100);
    assert(_mpz_get_s64(z) == -100);

    // INT64_MAX
    _mpz_set_s64(z, INT64_MAX);
    assert(_mpz_get_s64(z) == INT64_MAX);

    // INT64_MIN（边界：避免 -INT64_MIN 溢出）
    _mpz_set_s64(z, INT64_MIN);
    assert(_mpz_get_s64(z) == INT64_MIN);

    // 零
    _mpz_set_s64(z, 0);
    assert(_mpz_get_s64(z) == 0);
    mpz_clear(z);
}

void test_mpz_u64_helpers() {
    mpz_t z;
    mpz_init(z);

    _mpz_set_u64(z, 0);
    assert(_mpz_get_u64(z) == 0);

    _mpz_set_u64(z, UINT64_MAX);
    assert(_mpz_get_u64(z) == UINT64_MAX);

    _mpz_set_u64(z, UINT64_C(18446744073709551557));  // 2^64 - 59
    assert(_mpz_get_u64(z) == UINT64_C(18446744073709551557));
    mpz_clear(z);
}
```

---

## M0.2: 素数生成组件

**文件**：`clpoly/number/ZZ.hh`，在 ZZ 类定义之后、namespace 结尾之前。

```cpp
// ---- prime generation (GMP-based) ----

/// 返回最小素数 p > n。n >= 2。
/// 若 p 超出 uint64_t 范围抛出 std::overflow_error。
inline uint64_t next_prime_64(uint64_t n) {
    mpz_t z;
    mpz_init(z);
    _mpz_set_u64(z, n);
    mpz_nextprime(z, z);
    if (mpz_sizeinbase(z, 2) > 64) {
        mpz_clear(z);
        throw std::overflow_error("next_prime_64: result exceeds uint64_t");
    }
    uint64_t result = _mpz_get_u64(z);
    mpz_clear(z);
    return result;
}

/// 返回最大素数 p < n。n >= 3。
/// 需要 GMP >= 6.3.0（mpz_prevprime）。
inline uint64_t prev_prime_64(uint64_t n) {
    if (n <= 2)
        throw std::domain_error("prev_prime_64: no prime less than 2");
    mpz_t z;
    mpz_init(z);
    _mpz_set_u64(z, n);
    mpz_prevprime(z, z);
    uint64_t result = _mpz_get_u64(z);
    mpz_clear(z);
    return result;
}

/// 返回最小素数 p > n（任意精度版本）。
inline ZZ next_prime(const ZZ& n) {
    mpz_t z;
    mpz_init(z);
    if (n.is_small()) {
        _mpz_set_s64(z, n.get_val());
    } else {
        mpz_set(z, n.get_mpz_v());
    }
    mpz_nextprime(z, z);
    ZZ result(static_cast<mpz_srcptr>(z));
    mpz_clear(z);
    return result;
}
```

**注意**：`_mpz_set_u64` / `_mpz_get_u64` / `_mpz_set_s64` / `_mpz_get_s64` 已在 M0.1 中定义，此处直接使用。

**测试**：

```cpp
void test_next_prime_64() {
    assert(next_prime_64(2) == 3);
    assert(next_prime_64(10) == 11);
    assert(next_prime_64(100) == 101);
    // 大素数附近
    uint64_t p = next_prime_64(UINT64_C(1) << 62);
    assert(p > (UINT64_C(1) << 62));
}

void test_prev_prime_64() {
    assert(prev_prime_64(3) == 2);
    assert(prev_prime_64(12) == 11);
    // 验证 2^64 - 59 的前一个素数
    uint64_t p = UINT64_C(18446744073709551557);  // 2^64 - 59
    uint64_t pp = prev_prime_64(p);
    assert(pp < p);
    // prev_prime_64(p) 再 next 回来应该 == p
    assert(next_prime_64(pp) == p);
}

void test_next_prime_zz() {
    ZZ big("99999999999999999999");
    ZZ p = next_prime(big);
    assert(p > big);
}
```

---

## M0.4: Zp 类大素数安全修复

**文件**：`clpoly/number.hh`

### 修复 1-2：构造函数 `(int64_t)p` 溢出

**当前代码**（line 142-147）：

```cpp
Zp(int64_t i, uint64_t p) : _p(p)
{
    __precompute(p, _ninv, _norm);
    int64_t r = i % (int64_t)p;           // ← 溢出
    _i = (uint64_t)(r >= 0 ? r : r + (int64_t)p);  // ← 溢出
}
```

**改为**：

```cpp
Zp(int64_t i, uint64_t p) : _p(p)
{
    __precompute(p, _ninv, _norm);
    uint64_t abs_i = (i >= 0) ? (uint64_t)i : -(uint64_t)i;
    uint64_t r = abs_i % p;
    _i = (r == 0 || i > 0) ? r : p - r;
}
```

**同样修改 `operator=`**（line 161-167）：

```cpp
constexpr Zp& operator=(int64_t i)
{
    if (this->_p == 0) { assert(i == 0); this->_i = 0; return *this; }
    uint64_t abs_i = (i >= 0) ? (uint64_t)i : -(uint64_t)i;
    uint64_t r = abs_i % this->_p;
    this->_i = (r == 0 || i > 0) ? r : this->_p - r;
    return *this;
}
```

### 修复 3-4：加法溢出

**当前代码**（line 190-195）：

```cpp
friend inline Zp operator+(Zp op1, const Zp& op2)
{
    assert(op1._p == op2._p);
    uint64_t r = op1._i + op2._i;        // ← p > 2^63 时溢出
    op1._i = (r >= op1._p) ? r - op1._p : r;
    return op1;
}
```

**改为**：

```cpp
friend inline Zp operator+(Zp op1, const Zp& op2)
{
    assert(op1._p == op2._p);
    const uint64_t neg = op1._p - op1._i;  // p - a，不溢出（a < p）
    if (neg > op2._i)
        op1._i = op1._i + op2._i;          // a + b < p，安全
    else
        op1._i = op2._i - neg;             // = a + b - p
    return op1;
}
```

**同样修改 `operator+=`**（line 197-203）：

```cpp
inline Zp& operator+=(const Zp& op2)
{
    assert(this->_p == op2._p);
    const uint64_t neg = this->_p - this->_i;
    if (neg > op2._i)
        this->_i = this->_i + op2._i;
    else
        this->_i = op2._i - neg;
    return *this;
}
```

### 修复 5-6：放宽 assert

**`__precompute`**（line 109）：

```cpp
// 当前：assert(p >= 2 && p < (1ULL << 63));
// 改为：
assert(p >= 2);
```

**`set_prime`**（line 136）：

```cpp
// 当前：assert(p >= 2 && p < (1ULL << 63));
// 改为：
assert(p >= 2);
```

### Barrett 乘法验证

`__nmod_mul` 使用 `unsigned __int128` 全程无符号运算，对任何 `p < 2^64` 均安全。`__preinvert_limb` 要求 `pn = p << norm` 最高位为 1：
- `p > 2^63` 时 `norm = 0`，`pn = p`，最高位为 1 ✓
- `p ≤ 2^63` 时 `norm = clzll(p) ≥ 0`，`pn = p << norm`，最高位为 1 ✓

无需修改。

### 测试

```cpp
void test_zp_large_prime() {
    // 2^64 - 59，最大 64-bit 素数
    uint64_t p = UINT64_C(18446744073709551557);

    // 构造函数
    Zp a(-1LL, p);
    assert(a.number() == p - 1);

    Zp b(INT64_MIN, p);
    // INT64_MIN mod p = p - (|INT64_MIN| mod p)
    uint64_t expected = p - ((uint64_t)INT64_MAX + 1) % p;
    if (expected == p) expected = 0;
    assert(b.number() == expected);

    // 加法
    Zp x(p - 1, p), y(2, p);
    Zp sum = x + y;
    assert(sum.number() == 1);  // (p-1) + 2 = p + 1 ≡ 1

    // 乘法
    Zp m1(p - 1, p), m2(p - 1, p);
    Zp prod = m1 * m2;
    assert(prod.number() == 1);  // (-1) * (-1) = 1

    // 减法（已确认安全，但验证一下）
    Zp s1(1, p), s2(p - 1, p);
    Zp diff = s1 - s2;
    assert(diff.number() == 2);  // 1 - (p-1) = 2 - p ≡ 2
}
```

---

## M0.3: boost::math::prime 完全迁移

**依赖**：M0.2 完成。

### 步骤 1：删除 include

```
polynomial_gcd.hh  line 10:  删除 #include <boost/math/special_functions/prime.hpp>
polynomial_factorize_univar.hh line 17:  删除同上
```

### 步骤 2：逐点替换

**polynomial_gcd.cc**（单变量 GCD）：

| 行号 | 当前 | 改为 |
|------|------|------|
| 43-44 | `uint32_t p_index=0; uint64_t prime=boost::math::prime(p_index);` | `uint64_t prime = UINT64_C(18446744073709551557);  // 2^64 - 59` |
| 57-61 | `while (...) { if (++p_index >= 9999) break; prime=boost::math::prime(p_index); }` | `while (...) { prime = prev_prime_64(prime); }` |
| 62 | `if (p_index >= 9999) break;` | 删除（prev_prime_64 无下限问题，实际最多用 ~7 个素数） |
| 73-77 | 同上模式 | 同上 |
| 194-195 | `if (++p_index >= 9999) break; prime=boost::math::prime(p_index);` | `prime = prev_prime_64(prime);` |

**polynomial_gcd.hh**（多变量 GCD）：

| 行号 | 当前 | 改为 |
|------|------|------|
| 310-318 | `tmp_x/log(tmp_x)` 查表找 > degree 的素数 | `uint64_t prime = UINT64_C(18446744073709551557);  // 2^64 - 59` |
| 312-313 | `p_index` 初始化 + 9999 检查 | 删除 |
| 331-335 | lc 检查循环内 `++p_index; prime=boost::math::prime(p_index);` | `prime = prev_prime_64(prime);` |
| 336 | `if (p_index >= 9999) break;` | 删除 |
| 347-351 | 同上 | 同上 |
| 463-464 | 同上 | 同上 |

**polynomial_factorize_univar.hh**（`__select_prime`）：

| 行号 | 当前 | 改为 |
|------|------|------|
| 1397 | `boost::math::prime((unsigned)idx)` | `prime = next_prime_64(prime)` |

当前代码用 `idx` 做随机访问（`boost::math::prime(idx)` 返回第 idx 个素数）。迁移后改为顺序推进：

```cpp
// 当前结构：
for (size_t idx = start_idx; tried < max_tries; ++idx) {
    uint64_t p = boost::math::prime((unsigned)idx);
    // ... 检查 p 是否可用 ...
}

// 改为：
uint64_t p = 2;
for (size_t tried = 0; tried < max_tries; ) {
    // ... 检查 p 是否可用 ...
    p = next_prime_64(p);  // 推进到下一个素数
}
```

**注意**：因式分解的 `__select_prime` 仍从**小素数**开始（p=2），不改为大素数。因式分解需要小素数做 Zp 分解 + DDF/EDF，大素数无优势（因子数不变，且 Zp 算术成本更高）。

**polynomial_factorize_wang.hh**（Wang 因式分解中的素数选择）：

| 行号 | 当前 | 改为 |
|------|------|------|
| 2024 | `boost::math::prime(idx)` | `prime = next_prime_64(prime)` |

同上，从小素数顺序推进。重构方式与 `__select_prime` 相同。

### 步骤 3：清理

- 删除所有 `p_index` / `uint32_t p_index` 变量声明
- 删除所有 `>= 9999` 分支
- 删除 `PRIME_TABLE_SIZE` 相关常量（如存在）
- 确认 Makefile 不再链接 boost_math（当前仅 header-only，无需改 Makefile）

### 验证

```bash
make test/test_gcd && _build/debug/bin/test_gcd
make test/test_factorize && _build/debug/bin/test_factorize
bash test/run_all_tests.sh
make crosscheck
```

---

## M0.5: MTSHL 素数选择修复

**依赖**：M0.2（`prev_prime_64`）+ M0.4（Zp 支持 `p > 2^63`）

**文件**：`clpoly/polynomial_factorize_wang.hh`

### 改动 1：删除硬编码，动态选择素数

**当前**（line 2272）：
```cpp
static constexpr uint64_t mtshl_p = 9223372036854775783ULL;
```

**改为**：per-variable 缓存素数，避免多轮重复检查：

```cpp
// 在 for(;;) 循环外声明
std::vector<uint64_t> var_mtshl_p(main_vars.size(), 0);  // 0 = 未初始化

// ... 进入 for(;;) → vi 循环后：
for (size_t vi = 0; vi < main_vars.size(); ++vi)
{
    if (var_dead[vi]) continue;
    variable x1 = main_vars[vi];

    // ---- 首次访问此变量时，选择安全的 mtshl_p ----
    if (var_mtshl_p[vi] == 0) {
        var_mtshl_p[vi] = UINT64_C(18446744073709551557);  // 2^64 - 59
        auto L = leadcoeff(g, x1);
        auto all_div = [&](uint64_t p) {
            for (const auto& term : L)
                if (term.second.fdiv_ui(p) != 0) return false;
            return true;
        };
        while (all_div(var_mtshl_p[vi]))
            var_mtshl_p[vi] = prev_prime_64(var_mtshl_p[vi]);
    }
    uint64_t mtshl_p = var_mtshl_p[vi];

    int batch_end = var_skip[vi] + BATCH_SIZE;
    for (int skip = var_skip[vi]; skip < batch_end; ++skip)
    {
        auto eval = __select_eval_point(g, x1, skip, mtshl_p);
        // ...
```

**说明**：`var_mtshl_p[vi] == 0` 表示首次访问，执行 content(lc) 检查并缓存结果；后续轮次直接复用。lc(g, x1) 在整个 `for(;;)` 循环中不变，缓存安全。不同主变量可能有不同的 lc，所以 per-variable 独立缓存。

### 改动 2：`__select_eval_point` 增加 mod p 检查

**当前**（line 1247-1258）：

```cpp
// 条件 (b): lc 非零
ZZ delta;
if (is_number(L))
    delta = L.front().second;
else
{
    auto L_eval = assign(L, alpha);
    if (L_eval.empty()) continue;
    if (!is_number(L_eval)) continue;
    delta = L_eval.front().second;
}
if (delta == 0) continue;
```

**在 `if (delta == 0) continue;` 之后追加**：

```cpp
// 条件 (b'): lc(f)(α) 不被 MTSHL 素数整除
if (lc_coprime_mod != 0 && delta.fdiv_ui(lc_coprime_mod) == 0) continue;
```

**函数签名变更**：

```cpp
// 当前：
__select_eval_point(const polynomial_<ZZ, lex_<var_order>>& f,
                    const variable& main_var, int skip = 0)

// 改为：
__select_eval_point(const polynomial_<ZZ, lex_<var_order>>& f,
                    const variable& main_var, int skip = 0,
                    uint64_t lc_coprime_mod = 0)
```

`lc_coprime_mod = 0` 表示不检查（兼容单变量因式分解等不需要 MTSHL 的调用方）。

### 改动 3：调用处传入 mtshl_p

**当前**（line 2245）：
```cpp
auto eval = __select_eval_point(g, x1, skip);
```

**改为**：
```cpp
auto eval = __select_eval_point(g, x1, skip, mtshl_p);
```

### 改动 4：删除原硬编码行

删除 line 2272 的 `static constexpr uint64_t mtshl_p = ...;`。

### 测试

```cpp
// 构造一个 lc 为 mtshl_p 倍数的多项式
void test_mtshl_bad_prime() {
    // lc(f) = 2 * (2^64 - 59)，这要求 ZZ 能表示 > 2^64 的整数
    ZZ big_lc = ZZ("36893488147419103114");  // 2 * (2^64-59)
    // 构造 f = big_lc * x^2 + x + 1 的多变量版本
    // 验证 MTSHL 素数选择跳过 2^64-59，选择下一个素数
    // 验证因式分解结果正确
}
```

实际验证通过 `make crosscheck`（FLINT 对照测试）。

---

## P0: 大素数 + content 修复

**依赖**：M0.1-M0.4 完成（需要 Zp 支持 64-bit 素数）。

**说明**：P0.1（素数起始位置）和 P0.3（移除 p_index）的代码改动已在 M0.3 中完成。P0 的独立改动仅为 P0.2（content 重复计算修复）。以下列出 P0 的完整逻辑，供验证 M0.3 改动的正确性。

### P0.1: 素数起始位置（已在 M0.3 中实施）

**polynomial_gcd.cc**（line 43-44）：

```cpp
// 当前：
std::uint32_t p_index=0;
std::uint64_t prime=boost::math::prime(p_index);

// 改为（M0.3 已完成此替换）：
uint64_t prime = UINT64_C(18446744073709551557);  // 2^64 - 59
```

**polynomial_gcd.hh**（line 310-318）：

```cpp
// 当前：
std::uint32_t tmp_x=std::max(degree(F),degree(G));
if (tmp_x<2) tmp_x=2;
std::uint32_t p_index=tmp_x/std::log(tmp_x);
if (p_index >= 9999) p_index = 9998;
std::uint64_t prime=boost::math::prime(p_index);
while (prime <tmp_x)
{
    if (++p_index >= 9999) break;
    prime=boost::math::prime(p_index);
}

// 改为：
uint64_t prime = UINT64_C(18446744073709551557);  // 2^64 - 59
```

CRT 循环内取下一个素数（M0.3 替换表已覆盖）：

```cpp
// 当前：
if (++p_index >= 9999) break;
prime = boost::math::prime(p_index);

// 改为：
prime = prev_prime_64(prime);
```

**说明**：当前代码选素数 > degree 是为了保证 Zp 下不退化。`2^64 - 59` 远大于任何实际 degree，此条件自动满足。素数从最大向前取（`prev_prime_64`），与 MTSHL 策略统一。`p_index` 和 9999 上限检查一并删除（P0.3）。

### P0.2: content 重复计算修复

**polynomial_gcd.hh line 299**：

```cpp
// 当前：
polynomial_<ZZ,lex_<var_order>> cont_gcd=polynomial_GCD(cont(F),cont(G));

// 改为：
polynomial_<ZZ,lex_<var_order>> cont_gcd=polynomial_GCD(F_cont,G_cont);
```

复用 line 297-298 已计算的 `F_cont` 和 `G_cont`。

### P0.3: 移除 p_index 和 9999 限制（已在 M0.3 中实施）

验证所有 `p_index` 变量和 `>= 9999` 检查已删除。

### CRT 循环中 Zp 构造的安全性

CRT 合并代码（`polynomial_gcd.cc` line 107, `polynomial_gcd.hh` line 375）：

```cpp
Zp tmp_inv(Pout_prime, prime);
```

这里 `Pout_prime` 是 ZZ（CRT 累积乘积），`prime` 是 uint64_t。调用的是 `Zp(const ZZ& i, uint64_t p)` 构造函数（line 149），内部使用 `fdiv_ui` 做 unsigned 取模，不受 `(int64_t)p` 问题影响。**无需修改**。

### 验证

```bash
# 1. 全量测试
bash test/run_all_tests.sh

# 2. FLINT 对照
make crosscheck

# 3. 性能基准
make bench-clpoly
# 对比 benchmarks/2026-02-23.txt（优化前基准）
# 预期：GCD 相关 benchmark 提速 20-30x
```

### 预期效果

以 Benchmark `gcd_univar_400bit` 为例：

| 指标 | 优化前 | 优化后 |
|------|--------|--------|
| 起始素数 | 2 | 2^64 - 59 |
| 每素数贡献 | ~1-17 bits | ~64 bits |
| CRT 循环次数 | ~200 | ~7 |
| 总时间 | ~Xms | ~X/30 ms |

---

## P1: GCDHEU（启发式 GCD）

**依赖**：无（独立于 P0，仅依赖已有 ZZ 接口）

**理论依据**：Parisse 2002（arXiv:cs/0206032v1）Theorem 1 + FLINT `fmpz_poly/gcd_heuristic.c`

### 算法概述

GCDHEU 在 ξ = 2^pack_bits 处求值多项式为大整数，用 GMP `mpz_gcd` 求整数 GCD，再通过 **对称 ξ-adic 重构** 还原为多项式并验证。

```
Input:  F, G ∈ Z[x]  （已提取 content，primitive）
Output: gcd(F, G)  或  FAIL（回退到模算法）

1. bits_f = max_{coeff c in F} sizeinbase(|c|, 2)
   bits_g = max_{coeff c in G} sizeinbase(|c|, 2)
   if bits_f + bits_g >= 128:  return FAIL

2. pack_bits = max(min(bits_f, bits_g) + 6, max(bits_f, bits_g) + 1)
   ξ = 2^pack_bits

3. f_val = F(ξ)   // Horner 求值: shift-and-add（有符号）
   g_val = G(ξ)

4. g_int = gcd(|f_val|, |g_val|)    // GMP 大整数 GCD（始终 ≥ 0）

5. H = symmetric_ξ_adic(g_int)      // 对称 ξ-adic 重构 → 多项式

6. H = pp(H)                        // 提取原始部分

7. 验证: H | F  且  H | G ?
   - 是 → return H * cont_gcd
   - 否 → return FAIL
```

**正确性（Parisse 2002 Theorem 1）**：

> 设 P, Q ∈ Z[X₁,...,Xₖ]，z 为整数，|z| ≥ 2·min(|P|, |Q|) + 2（|P| 为 P 的最大系数绝对值）。若 gcd(P(z), Q(z)) 的 z-adic 对称重构的原始部分 G 整除 P 和 Q，则 G = gcd(P, Q)。

**下界验证**：`pack_bits = min_bits + 6` ⇒ `ξ = 2^(min_bits+6) = 64·2^min_bits`。
最大系数 ≤ `2^min_bits - 1`。`ξ ≥ 64·2^min_bits ≫ 2·(2^min_bits - 1) + 2`。✓

**适用条件**：单变量 Z[x]，`bits_f + bits_g < 128`。Benchmark 中系数 ∈ [-20, 20]（5 bits），总是满足。

### 关键数学原理：为什么 Horner 求值而非 offset 编码

~~原设计使用 offset 编码（系数 + 2^(pack_bits-1) 使非负），这是**错误的**~~：

```
offset 编码: a_int = Σ (c_i + offset)·ξ^i = F(ξ) + offset·Σξ^i ≠ F(ξ)
```

Parisse 定理要求 `gcd(F(z), G(z))`，不是 `gcd(F(z)+extra, G(z)+extra)`。offset 项破坏了求值的正确性。

**正确做法**：直接计算 `F(ξ)`（Horner 方法，允许中间结果为负）。FLINT 用二补码 + 借位链在 limb 层面实现同一语义；CLPoly 用 mpz 有符号算术更简洁。

### 文件：`clpoly/polynomial_gcd.cc`

#### 改动 1：新增 `__gcdheu` 函数

在 `polynomial_GCD` 函数之前，新增静态辅助函数：

```cpp
// ---- GCDHEU: 启发式 GCD (Parisse 2002 Theorem 1) ----
// 返回 true = 成功（result 为 primitive GCD），false = 失败（回退到模算法）
static bool __gcdheu(
    upolynomial_<ZZ>& result,
    const upolynomial_<ZZ>& F,
    const upolynomial_<ZZ>& G)
{
    // F, G 必须已 primitive、非空、size > 1

    // Step 1: 计算最大系数位长
    size_t bits_f = 0, bits_g = 0;
    for (auto& term : F) {
        size_t b = term.second.sizeinbase(2);  // sizeinbase 已取绝对值
        if (b > bits_f) bits_f = b;
    }
    for (auto& term : G) {
        size_t b = term.second.sizeinbase(2);
        if (b > bits_g) bits_g = b;
    }

    // 阈值检查：FLINT 用 2*FLINT_BITS = 128
    if (bits_f + bits_g >= 128) return false;

    // Step 2: pack_bits (Parisse bound 要求 +3，FLINT 用 +6 做安全裕量)
    size_t mn = std::min(bits_f, bits_g);
    size_t mx = std::max(bits_f, bits_g);
    size_t pack_bits = std::max(mn + 6, mx + 1);

    // Step 3: 稠密化 — 将 F, G 转为 dense 系数数组 [c_0, c_1, ..., c_deg]
    int64_t deg_f = get_deg(F);
    int64_t deg_g = get_deg(G);

    std::vector<ZZ> fc(deg_f + 1), gc(deg_g + 1);
    for (auto& term : F) fc[term.first.deg()] = term.second;
    for (auto& term : G) gc[term.first.deg()] = term.second;

    // Step 4: Horner 求值 F(ξ), G(ξ)  其中 ξ = 2^pack_bits
    // 乘以 ξ 即左移 pack_bits 位，加系数直接 mpz_add（允许负数）
    mpz_t f_val, g_val, g_int, tmp;
    mpz_inits(f_val, g_val, g_int, tmp, NULL);

    // Pack F: f_val = c[deg] * ξ^deg + ... + c[0]
    // Horner: f_val = (...((c[deg]) * ξ + c[deg-1]) * ξ + ...) * ξ + c[0]
    mpz_set_ui(f_val, 0);
    for (int64_t i = deg_f; i >= 0; --i) {
        mpz_mul_2exp(f_val, f_val, pack_bits);       // f_val <<= pack_bits (×ξ)
        if (fc[i].is_small()) {
            int64_t v = fc[i].get_val();
            if (v >= 0) mpz_add_ui(f_val, f_val, (uint64_t)v);
            else        mpz_sub_ui(f_val, f_val, (uint64_t)(-v));
        } else {
            mpz_add(f_val, f_val, fc[i].get_mpz_v());  // 有符号加法
        }
    }

    // Pack G
    mpz_set_ui(g_val, 0);
    for (int64_t i = deg_g; i >= 0; --i) {
        mpz_mul_2exp(g_val, g_val, pack_bits);
        if (gc[i].is_small()) {
            int64_t v = gc[i].get_val();
            if (v >= 0) mpz_add_ui(g_val, g_val, (uint64_t)v);
            else        mpz_sub_ui(g_val, g_val, (uint64_t)(-v));
        } else {
            mpz_add(g_val, g_val, gc[i].get_mpz_v());
        }
    }

    // Step 5: 整数 GCD（取绝对值确保非负）
    mpz_abs(f_val, f_val);
    mpz_abs(g_val, g_val);
    mpz_gcd(g_int, f_val, g_val);

    // Step 6: 对称 ξ-adic 重构 g_int → 多项式 H
    //   重复: r = g_int mod ξ (对称), h[i] = r, g_int = (g_int - r) / ξ
    mpz_t mask, half_xi, coeff_mpz;
    mpz_inits(mask, half_xi, coeff_mpz, NULL);
    mpz_set_ui(mask, 1);
    mpz_mul_2exp(mask, mask, pack_bits);    // mask = ξ = 2^pack_bits
    mpz_set(half_xi, mask);
    mpz_tdiv_q_ui(half_xi, half_xi, 2);    // half_xi = ξ/2
    mpz_sub_ui(mask, mask, 1);             // mask = ξ - 1 (用于提取低位)

    upolynomial_<ZZ> H;
    int64_t deg_bound = std::min(deg_f, deg_g);
    for (int64_t i = 0; i <= deg_bound && mpz_sgn(g_int) != 0; ++i) {
        // 提取低 pack_bits 位: r = g_int & (ξ-1), 0 ≤ r < ξ
        mpz_and(coeff_mpz, g_int, mask);
        // 右移（先做，因为 and 已保存低位）
        mpz_tdiv_q_2exp(g_int, g_int, pack_bits);
        // 对称表示: 若 r > ξ/2，则 r -= ξ 并向高位进位 +1
        if (mpz_cmp(coeff_mpz, half_xi) > 0) {
            mpz_sub(coeff_mpz, coeff_mpz, mask);
            mpz_sub_ui(coeff_mpz, coeff_mpz, 1);   // r -= ξ (= r - (mask+1))
            mpz_add_ui(g_int, g_int, 1);            // 进位
        }

        if (mpz_sgn(coeff_mpz) != 0)
            H.push_back({umonomial(i), ZZ(coeff_mpz)});
    }

    mpz_clears(f_val, g_val, g_int, tmp, mask, half_xi, coeff_mpz, NULL);

    if (H.empty()) return false;

    // 降幂排列（upolynomial 约定高次在前）
    std::reverse(H.data().begin(), H.data().end());
    H.normalization();

    // Step 7: pp(H) — 提取原始部分，确保首项正
    ZZ h_cont = cont(H);
    if (h_cont != 1 && h_cont != -1)
        for (auto& term : H) term.second /= h_cont;
    if (!H.empty() && H.front().second < 0)
        for (auto& term : H) term.second = -term.second;

    // Step 8: 验证 H | F 且 H | G（多项式精确除法）
    upolynomial_<ZZ> Q, R;
    pair_vec_div(Q.data(), R.data(), F.data(), H.data(), F.comp());
    if (!R.empty()) return false;

    pair_vec_div(Q.data(), R.data(), G.data(), H.data(), F.comp());
    if (!R.empty()) return false;

    result = std::move(H);
    return true;
}
```

**设计决策**：

- **Horner 求值（而非 offset 编码）**：直接计算 `F(2^pack_bits)` 为有符号大整数。乘以 ξ = 2^pack_bits 仅需位移（`mpz_mul_2exp`），加系数用 `mpz_add/sub`。结果可能为负，用 `mpz_abs` 取绝对值后求 GCD。FLINT 通过二补码借位链在 limb 层面实现同一语义，但需要复杂的 `fmpz_bit_pack/unpack` 机制；CLPoly 用 mpz 有符号算术更简洁，开销仅为每系数一次 `mpz_add_ui`/`mpz_sub_ui`。
- **对称 ξ-adic 重构**：从低位逐窗口提取系数，`r > ξ/2` 时做 `r -= ξ` 并向高位进位 +1。这与 FLINT 的 `fmpz_bit_unpack` 中的 borrow 传播等价。
- **稠密数组**：GCDHEU 仅在小系数（< 64 bits）情况下调用，多项式通常稠密。稀疏情况下稠密数组浪费少量内存但简化了索引逻辑。
- **验证用多项式精确除法**：FLINT 先用整数整除检查（`flint_mpn_divides`），再用位长检查或乘法验证。CLPoly 直接做多项式除法（`pair_vec_div`），更简洁。对于 GCDHEU 的目标场景（小系数），验证开销远小于 GCD 本身。

#### 改动 2：在 `polynomial_GCD` 中插入快速路径

**位置**：`polynomial_gcd.cc` line 42（content 提取之后，CRT 循环之前）

```cpp
        for (auto &i:G)
            i.second/=g_cont;

        // GCDHEU 快速路径
        {
            upolynomial_<ZZ> heu_result;
            if (__gcdheu(heu_result, F, G)) {
                // 恢复 content
                if (cont_gcd != 1) {
                    for (auto& term : heu_result) term.second *= cont_gcd;
                }
                return heu_result;
            }
        }

        uint64_t prime = UINT64_C(18446744073709551557);  // 2^64 - 59
        // ... 现有 CRT 循环
```

**说明**：

- GCDHEU 在 content 提取后调用，输入是 primitive 多项式
- 成功时乘回 `cont_gcd` 并直接返回，跳过整个 CRT 循环
- 失败时无开销地进入模算法（`__gcdheu` 内部的位长检查几乎免费）
- 无阈值检查包装：`__gcdheu` 内部已有 `bits_f + bits_g >= 128` 检查

### ZZ 接口需求

| 需求 | 已有接口 | 说明 |
|------|---------|------|
| 判断小/大模式 | `is_small()` (M0.1 新增) | ✓ |
| 读取小值 | `get_val()` (M0.1 新增) | ✓ |
| 读取 mpz | `get_mpz()` (M0.1 新增) | ✓ |
| 系数位长 | `sizeinbase(2)` （内部已取绝对值）| ✓ |
| ZZ(mpz_srcptr) 构造 | `ZZ(mpz_srcptr)` (M0.1 新增) | ✓ |

**所有接口已在 M0.1 中实现，无需新增。**

### 与 FLINT 实现的差异

| 方面 | FLINT | CLPoly (P1) | 说明 |
|------|-------|-------------|------|
| Pack 方式 | 二补码 + 借位链（limb 级操作） | Horner 求值（mpz 有符号算术） | 数学等价，CLPoly 更简洁 |
| Unpack 方式 | `fmpz_bit_unpack` + sign bit + borrow | 对称 mod + 进位 | 数学等价 |
| 验证 | 整数整除 → 位长检查 → 乘法验证 | 多项式精确除法 | CLPoly 更简洁 |
| content 处理 | packed limb 原地除 content | unpack 后 `cont()` + 除法 | CLPoly 更简洁 |
| 特殊情况 | len < 6 用 subresultant | 无（直接尝试 GCDHEU） | 可后续添加 |
| 首项符号 | pack 前 negate 使首项正 | GCD 后 `mpz_abs` | 等价 |

### 测试

```cpp
// 1. 小系数已知 GCD
void test_gcdheu_small_coeff() {
    // f = (x+1)(x+2) = x^2 + 3x + 2
    // g = (x+1)(x-3) = x^2 - 2x - 3
    // gcd = x + 1
    // 系数位长: 2 bits, 总计 4 < 128 → GCDHEU 应触发
}

// 2. 互素多项式
void test_gcdheu_coprime() {
    // gcd(x^2+1, x^2-1) = 1
    // GCDHEU 成功返回 1
}

// 3. 大系数回退
void test_gcdheu_fallback() {
    // 系数 > 2^64 → bits_f + bits_g > 128 → GCDHEU 跳过
    // 验证结果仍正确（由模算法计算）
}

// 4. 随机测试
void test_gcdheu_random() {
    // 随机小系数多项式，验证 gcd 结果与已有实现一致
}
```

实际验证通过 `bash test/run_all_tests.sh` + `make crosscheck`。

### 预期效果

| Benchmark | 当前 (P0) | 预期 (P1) | 说明 |
|-----------|----------|----------|------|
| gcd deg30+common15 | ~0.1ms | ~0.02ms | 小系数，GCDHEU 命中 |
| gcd deg80+common40 | ~0.5ms | ~0.1ms | 小系数，GCDHEU 命中 |
| gcd deg200+common100 | ~1.9ms | ~0.5ms | 小系数，GCDHEU 命中 |
| gcd deg8+common4 (multivar) | ~0.5ms | ~0.3ms | 间接受益（递归 GCD 调用） |
| gcd deg15+common8 (multivar) | ~1.6ms | ~0.8ms | 间接受益 |

---

## 变更摘要

| 文件 | 改动类型 | 涉及模块 |
|------|---------|---------|
| `clpoly/number/ZZ.hh` | GMP 64-bit helper + 新增接口 + LLP64 整体梳理(~45处) + 素数函数 | M0.1, M0.2 |
| `clpoly/number.hh` | Zp bug fix（6 处） | M0.4 |
| `clpoly/polynomial_gcd.cc` | 素数策略 + 删 boost + GCDHEU 快速路径 | M0.3, P0, P1 |
| `clpoly/polynomial_gcd.hh` | 素数策略 + content fix + 删 boost | M0.3, P0 |
| `clpoly/polynomial_factorize_univar.hh` | 删 boost include + 素数推进 | M0.3 |
| `clpoly/polynomial_factorize_wang.hh` | MTSHL 素数 + eval point 检查 | M0.5 |
