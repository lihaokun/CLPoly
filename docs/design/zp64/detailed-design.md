# M6a 细化设计文档：Zp 类 64-bit 扩展

> 阶段：细化（workflow.md §2.3）
> 前置文档：`docs/design/zp64/architecture.md`
> 目标文件：主要 `clpoly/number.hh` + `clpoly/number/ZZ.hh`，下游机械适配
> 实施里程碑：Step 1 → Step 4（顺序推进）

---

## Step 1：`number/ZZ.hh` — `fdiv_ui` 参数扩展

**文件**：`clpoly/number/ZZ.hh` L598-607

```cpp
// 旧（L599）：
uint64_t fdiv_ui(uint32_t d) const {

// 新：
uint64_t fdiv_ui(uint64_t d) const {
```

内部逻辑不变：小整数路径 `_val % d`（uint64 除法自动适配），大整数路径 `mpz_fdiv_ui(_mpz, d)`（GMP 接受 `unsigned long` = 64-bit）。

验收：编译通过，test_number 通过。

---

## Step 2：`number.hh` — `inv_prime` + `Zp` 类重写（核心）

**文件**：`clpoly/number.hh` L72-320

### 2a. `inv_prime`（L72-84）

```cpp
// 旧：
inline uint64_t inv_prime(uint64_t _i, uint32_t _p)
{
    assert(_p!=0 && _i!=0);
    uint64_t a=_p,b=_i,c;
    uint64_t s1=0,s2=1,s3;
    while (c=(a%b))
    {
        s3=(s1+_p-(s2*(a/b))%_p)%_p;
        a=b;b=c;
        s1=s2;s2=s3;
    }
    return s2;
}

// 新：
inline uint64_t inv_prime(uint64_t _i, uint64_t _p)
{
    assert(_p != 0 && _i != 0);
    uint64_t a = _p, b = _i, c;
    uint64_t s1 = 0, s2 = 1, s3;
    while ((c = a % b))
    {
        uint64_t q = a / b;
        uint64_t sq_mod = (unsigned __int128)s2 * q % _p;
        s3 = s1 >= sq_mod ? s1 - sq_mod : _p - sq_mod + s1;
        a = b; b = c;
        s1 = s2; s2 = s3;
    }
    return s2;
}
```

改动点：
1. `uint32_t _p` → `uint64_t _p`
2. `s2*(a/b)` → `(unsigned __int128)s2 * q`（防 64-bit 溢出）
3. `s1+_p-sq_mod` → 条件减法（防 `s1+_p` 溢出 uint64）

### 2b. `Zp` 类数据成员（L85-92）

```cpp
// 旧：
uint32_t _i;     // 值 ∈ [0, p)
uint32_t _p;     // 素数 p（0 = 未初始化哨兵）
uint64_t _ninv;  // Barrett 常数：UINT64_MAX / p（0 = 未初始化）
inline static uint32_t _s_prime = 0;

// 新：
uint64_t _i;     // 值 ∈ [0, p)
uint64_t _p;     // 素数 p（0 = 未初始化哨兵）
uint64_t _ninv;  // FLINT preinvert(p << _norm)
uint32_t _norm;  // __builtin_clzll(p)（p=0 时为 0）
inline static uint64_t _s_prime = 0;
```

### 2c. 预计算函数（L94-111）

替换 `__barrett_ninv` 和 `__barrett_reduce`：

```cpp
// 新增：FLINT 归一化逆元（pn 必须最高位为 1）
static uint64_t __preinvert_limb(uint64_t pn)
{
    assert(pn >> 63);
    unsigned __int128 num = ((unsigned __int128)(~pn) << 64) | ~(uint64_t)0;
    return (uint64_t)(num / pn);
}

// 替换 __barrett_ninv：预计算 ninv 和 norm
static void __precompute(uint64_t p, uint64_t& ninv, uint32_t& norm)
{
    if (p == 0) { ninv = 0; norm = 0; return; }
    assert(p >= 2 && p < (1ULL << 63));
    norm = (uint32_t)__builtin_clzll(p);
    ninv = __preinvert_limb(p << norm);
}

// 替换 __barrett_reduce：FLINT 归一化 Barrett 乘法归约
uint64_t __nmod_mul(uint64_t a, uint64_t b) const
{
    uint64_t pn = _p << _norm;
    uint64_t a_shifted = a << _norm;
    unsigned __int128 prod = (unsigned __int128)a_shifted * b;
    uint64_t hi = (uint64_t)(prod >> 64);
    uint64_t lo = (uint64_t)prod;
    unsigned __int128 qm = (unsigned __int128)hi * _ninv;
    uint64_t q1 = (uint64_t)(qm >> 64);
    uint64_t q0 = (uint64_t)qm;
    q0 += lo;
    q1 += hi + (q0 < lo ? 1 : 0);
    uint64_t r = lo - (q1 + 1) * pn;
    if (r > q0) r += pn;
    if (r >= pn) r -= pn;
    return r >> _norm;
}
```

### 2d. 类级别素数管理（L115-116）

```cpp
// 旧：
static void     set_prime(uint32_t p) { assert(p >= 2 && p < (1u << 31)); _s_prime = p; }
static uint32_t cur_prime()           { return _s_prime; }

// 新：
static void     set_prime(uint64_t p) { assert(p >= 2 && p < (1ULL << 63)); _s_prime = p; }
static uint64_t cur_prime()           { return _s_prime; }
```

### 2e. 构造函数（L118-133）

```cpp
// 新：
Zp() : _i(0), _p(_s_prime) { __precompute(_p, _ninv, _norm); }

explicit Zp(uint64_t p) : _i(0), _p(p) { __precompute(p, _ninv, _norm); }

Zp(uint64_t i, uint64_t p) : _p(p) { __precompute(p, _ninv, _norm); _i = i % p; }

// 删除 Zp(uint32_t i, uint32_t p) — 被 Zp(uint64_t, uint64_t) 覆盖

Zp(int64_t i, uint64_t p) : _p(p)
{
    __precompute(p, _ninv, _norm);
    int64_t r = i % (int64_t)p;
    _i = (uint64_t)(r >= 0 ? r : r + (int64_t)p);
}

Zp(int i, uint64_t p) : _p(p)
{
    __precompute(p, _ninv, _norm);
    int64_t r = (int64_t)i % (int64_t)p;
    _i = (uint64_t)(r >= 0 ? r : r + (int64_t)p);
}

Zp(const ZZ& i, uint64_t p) : _p(p)
{
    __precompute(p, _ninv, _norm);
    _i = i.fdiv_ui(p);
}
```

### 2f. `inv()` 方法（L135-140）

```cpp
// 旧：
inline Zp inv() const
{
    Zp new_op(this->_p);
    new_op._i = (uint32_t)inv_prime(this->_i, this->_p);
    return new_op;
}

// 新：
inline Zp inv() const
{
    Zp new_op(this->_p);
    new_op._i = inv_prime(this->_i, this->_p);
    return new_op;
}
```

### 2g. 赋值运算符（L141-153）

```cpp
// 新 operator=(int64_t)：
constexpr Zp& operator=(int64_t i)
{
    if (this->_p == 0) { assert(i == 0); this->_i = 0; return *this; }
    int64_t r = i % (int64_t)this->_p;
    this->_i = (uint64_t)(r >= 0 ? r : r + (int64_t)this->_p);
    return *this;
}

// 新 operator=(const ZZ&)：
inline Zp& operator=(const ZZ& i)
{
    if (this->_p == 0) { assert(!i); this->_i = 0; return *this; }
    this->_i = i.fdiv_ui(this->_p);
    return *this;
}
```

### 2h. 访问器（L154-161）

```cpp
// 旧 → 新：
explicit constexpr operator std::uint64_t() const { return this->_i; }  // 不变
explicit constexpr operator bool() const { return this->_i != 0; }      // 不变
constexpr uint64_t  prime() const { return this->_p; }                  // uint32→uint64
constexpr uint64_t  number() const { return this->_i; }                 // uint32→uint64
constexpr uint64_t& number()       { return this->_i; }                 // uint32→uint64
constexpr void normalization() { assert(this->_p); this->_i %= this->_p; }  // 不变
void prime(uint64_t p) { this->_p = p; __precompute(p, _ninv, _norm); } // setter 同步 ninv+norm
```

### 2i. 取负（L163-167）

```cpp
// 新：
Zp operator-() const
{
    Zp r(*this);
    r._i = (_i == 0) ? 0 : _p - _i;
    return r;
}
```

注意：不能用 `Zp(p - _i, _p)` 构造（会触发 `__precompute`）。
用拷贝构造后修改 `_i` 避免重复预计算。

### 2j. 加减法（L169-194）

所有 `uint32_t r` → `uint64_t r`，去掉 `(uint32_t)` 强转：

```cpp
friend inline Zp operator+(Zp op1, const Zp& op2)
{
    assert(op1._p == op2._p);
    uint64_t r = op1._i + op2._i;
    op1._i = (r >= op1._p) ? r - op1._p : r;
    return op1;
}
// += / - / -= 同理
```

### 2k. 乘法（L195-211）

所有 `__barrett_reduce(...)` → `__nmod_mul(...)`：

```cpp
friend inline Zp operator*(Zp op1, const Zp& op2)
{
    assert(op1._p == op2._p);
    op1._i = op1.__nmod_mul(op1._i, op2._i);
    return op1;
}

friend inline Zp operator*(Zp op1, int64_t op2)
{
    op1._i = op1.__nmod_mul(op1._i, Zp(op2, op1._p)._i);
    return op1;
}

friend inline Zp operator*(int64_t op1, const Zp& op2) { return op2 * op1; }

inline Zp& operator*=(const Zp& op2)
{
    assert(this->_p == op2._p);
    this->_i = this->__nmod_mul(this->_i, op2._i);
    return *this;
}
```

### 2l. 除法（L213-224）

```cpp
friend inline Zp operator/(Zp op1, const Zp& op2)
{
    assert(op1._p == op2._p);
    op1._i = op1.__nmod_mul(op1._i, inv_prime(op2._i, op1._p));
    return op1;
}

inline Zp& operator/=(const Zp& op2)
{
    assert(this->_p == op2._p);
    this->_i = this->__nmod_mul(this->_i, inv_prime(op2._i, this->_p));
    return *this;
}
```

### 2m. 比较运算符（L225-256）

`operator>=` 需修改强转：

```cpp
// 旧 L255：
return op1._i >= (uint32_t)op2;

// 新：
return op1._i >= (uint64_t)op2;
```

其余比较运算符（`==`、`!=`）无需改动。

### 2n. 类外函数（L263-320）

`hash_value`、`zore_check`、`addmul`、`submul`、`pow`、`__div`：

这些函数通过 `.number()`、`.prime()` 和运算符间接使用 Zp，返回类型自动跟随。
**无需改动**（`addmul`/`submul` 中 `op.prime(op1.prime())` 的 setter 已在 §2h 中适配）。

验收：编译通过，test_number 中已有 Zp 测试通过 + 新增 63-bit 素数测试。

---

## Step 3：下游签名适配

### 3a. `upolynomial.hh`（L104）

```cpp
// 旧：
inline upolynomial_<Zp> polynomial_mod(const upolynomial_<ZZ>& p, uint32_t prime)

// 新：
inline upolynomial_<Zp> polynomial_mod(const upolynomial_<ZZ>& p, uint64_t prime)
```

### 3b. `polynomial_factorize_zp.hh`（L20）

```cpp
// 旧：
inline Zp __make_zp(int64_t val, uint32_t p) { return Zp(val, p); }

// 新：
inline Zp __make_zp(int64_t val, uint64_t p) { return Zp(val, p); }
```

### 3c. `polynomial_factorize_univar.hh`

4 处 `uint32_t p` → `uint64_t p`（L338, L395, L542, L1397）。
这些是函数参数或局部变量，值来自 `boost::math::prime()`（返回 uint32，隐式转 uint64）。

### 3d. `polynomial_factorize_wang.hh`

~15 处 `uint32_t p` → `uint64_t p`（所有函数参数和局部变量）。

关键位置：
- L455 `__mtshl_sparse_int` 参数
- L772 `__mtshl_step_j` 参数
- L939, L960, L976 辅助函数参数
- L1003 `__mtshl_lift` 参数

`mtshl_p` 更新（L2165）：
```cpp
// 旧：
static constexpr uint32_t mtshl_p = 2147483629;

// 新：
static constexpr uint64_t mtshl_p = 9223372036854775783ULL;
```

随机数分布（L485）：
```cpp
// 旧：
std::uniform_int_distribution<uint32_t> dist(1, p - 1);

// 新：
std::uniform_int_distribution<uint64_t> dist(1, p - 1);
```

### 3e. `polynomial_gcd.hh` + `polynomial_gcd.cc`

~6 处 `uint32_t prime` → `uint64_t prime`（局部变量）。
值来自 `boost::math::prime()`，隐式转换无问题。

验收：编译通过，全量测试通过。

---

## Step 4：测试

### 4a. 已有测试适配

`test/test_number.cc`、`test/test_mtshl_*.cc`、`test/test_factorize_zp.cc`：
所有 `uint32_t p = ...` → `uint64_t p = ...`。

### 4b. 新增 63-bit 素数测试

在 `test/test_number.cc` 中新增：

```cpp
CLPOLY_TEST("Zp 63-bit prime basic arithmetic");
{
    uint64_t p = 9223372036854775783ULL;  // 最大 < 2^63 的素数
    Zp a(12345678901234567ULL, p);
    Zp b(98765432109876543ULL, p);

    // 加法
    Zp c = a + b;
    CLPOLY_ASSERT_EQ(c, Zp((unsigned __int128)12345678901234567ULL + 98765432109876543ULL % p, p));

    // 乘法
    Zp d = a * b;
    uint64_t expected = (unsigned __int128)a.number() * b.number() % p;
    CLPOLY_ASSERT_EQ(d.number(), expected);

    // 逆元
    Zp e = a.inv();
    CLPOLY_ASSERT_EQ((a * e).number(), (uint64_t)1);

    // 负值构造
    Zp f(-1, p);
    CLPOLY_ASSERT_EQ(f.number(), p - 1);

    // ZZ 构造
    Zp g(ZZ("123456789012345678901234567890"), p);
    uint64_t expected_g = ZZ("123456789012345678901234567890").fdiv_ui(p);
    CLPOLY_ASSERT_EQ(g.number(), expected_g);
}

CLPOLY_TEST("Zp 63-bit prime pow and division");
{
    uint64_t p = 9223372036854775783ULL;
    Zp a(42, p);

    // Fermat 小定理：a^(p-1) = 1
    CLPOLY_ASSERT_EQ(pow(a, (int64_t)(p - 1)).number(), (uint64_t)1);

    // 除法
    Zp b(7, p);
    Zp c = a / b;
    CLPOLY_ASSERT_EQ((c * b).number(), a.number());
}
```

### 4c. 验收

```bash
make test/test_number && _build/debug/bin/test_number       # Step 1-2
bash test/run_all_tests.sh                                   # Step 3-4 全量
make bench-all                                               # 性能无退化
```

---

## 复用点汇总

| 复用对象 | 位置 | 说明 |
|---------|------|------|
| `__builtin_clzll` | GCC/Clang 内建 | 计算 norm |
| `unsigned __int128` | GCC/Clang 扩展 | 128-bit 乘法 |
| `mpz_fdiv_ui` | GMP | ZZ → Zp 转换（已支持 unsigned long = 64-bit） |
| `boost::math::prime` | Boost | 素数表（返回 uint32，隐式转 uint64） |
| FLINT `NMOD_MUL_PRENORM` 算法 | FLINT nmod.h | `__nmod_mul` 的参考实现 |
