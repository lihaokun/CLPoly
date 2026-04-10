# M6a 架构设计：Zp 类 64-bit 扩展

> 阶段：架构（workflow.md §2.2）
> 前置文档：`docs/research/zp64-research.md`（调研报告）
> 方案选择：方案 A — 统一 FLINT 归一化 Barrett

---

## 1. 核心流程

```
                    构造/初始化
                         │
                ┌────────┴────────┐
                │  __precompute   │  norm = clz(p)
                │                 │  pn = p << norm
                │                 │  ninv = preinvert(pn)
                └────────┬────────┘
                         │
          ┌──────────────┼──────────────┐
          │              │              │
     ┌────┴────┐   ┌────┴────┐   ┌────┴────┐
     │  加/减  │   │   乘法  │   │  逆/除  │
     │ uint64  │   │ Barrett │   │  EEA    │
     │ 条件减  │   │ 归一化  │   │ __int128│
     └─────────┘   └─────────┘   └─────────┘
```

所有素数（2 ≤ p < 2⁶³）走**同一条代码路径**，无分支。

---

## 2. 模块划分

### 2.1 数据成员

```cpp
class Zp {
private:
    uint64_t _i;      // 值 ∈ [0, p)
    uint64_t _p;      // 素数 p（0 = 未初始化哨兵）
    uint64_t _ninv;   // FLINT preinvert(p << _norm)
    uint32_t _norm;   // __builtin_clzll(p)（p=0 时为 0）

    inline static uint64_t _s_prime = 0;
};
```

对象大小：8 + 8 + 8 + 4 = 28，padding 到 32 bytes。

**约束**：`p ∈ {0} ∪ [2, 2⁶³)`。

- p = 0：未初始化哨兵（`addmul`/`submul` 延迟初始化语义保持）
- p < 2⁶³：保证加法 `a + b < 2p < 2⁶⁴` 不溢出 uint64

### 2.2 预计算模块

```
函数：__precompute(p)

功能描述：为素数 p 计算 Barrett 归约所需的归一化参数。

前置条件：
  - p = 0（哨兵，返回全零）或 2 ≤ p < 2^63

后置条件：
  - norm = __builtin_clzll(p)（p 的前导零位数）
  - pn = p << norm（归一化素数，最高位为 1）
  - ninv = floor(2^128 / pn) - 2^64
  - 由于 pn ∈ [2^63, 2^64)，ninv ∈ [0, 2^64)，恰好放入 uint64
```

**`__preinvert_limb(pn)` 计算**（参考 FLINT `n_preinvert_limb_prenorm`）：

```cpp
// 前置：pn 最高位为 1（即 pn ≥ 2^63）
// 后置：返回 floor((2^128 - pn·2^64 - 1) / pn) = floor(2^128/pn) - 2^64
static uint64_t __preinvert_limb(uint64_t pn)
{
    assert(pn >> 63);  // 最高位为 1
    unsigned __int128 num = ((unsigned __int128)(~pn) << 64) | ~(uint64_t)0;
    return (uint64_t)(num / pn);
}
```

推导：`~pn = 2^64 - 1 - pn`，`num = (~pn)·2^64 + (2^64-1) = 2^128 - pn·2^64 - 1`。

### 2.3 乘法归约模块

```
函数：__nmod_mul(a, b)

功能描述：计算 a·b mod p（归一化 Barrett 归约）。

前置条件：
  - 0 ≤ a, b < p < 2^63
  - _ninv, _norm 由 __precompute 正确初始化

后置条件：
  - 返回值 = a·b mod p，∈ [0, p)

不变式：
  - 中间乘积 a_shifted·b < pn·p < 2^64·2^63 = 2^127，放入 __int128
```

**算法**（对应 FLINT `NMOD_MUL_PRENORM` + `NMOD_RED2`）：

```cpp
uint64_t __nmod_mul(uint64_t a, uint64_t b) const
{
    uint64_t pn = _p << _norm;

    // Step 1: 归一化左操作数
    uint64_t a_shifted = a << _norm;  // a < p < 2^63 → a_shifted < pn ≤ 2^64-1

    // Step 2: 128-bit 乘法
    unsigned __int128 prod = (unsigned __int128)a_shifted * b;
    uint64_t hi = (uint64_t)(prod >> 64);
    uint64_t lo = (uint64_t)prod;

    // Step 3: Barrett 估计商
    unsigned __int128 qm = (unsigned __int128)hi * _ninv;
    uint64_t q1 = (uint64_t)(qm >> 64);
    uint64_t q0 = (uint64_t)qm;

    // Step 4: 加 hi:lo（128-bit 加法）
    q0 += lo;
    q1 += hi + (q0 < lo ? 1 : 0);  // 进位

    // Step 5: 近似余数
    uint64_t r = lo - (q1 + 1) * pn;

    // Step 6: 至多两次修正
    if (r > q0) r += pn;    // 修正下溢
    if (r >= pn) r -= pn;   // 修正超范围

    // Step 7: 反归一化
    return r >> _norm;
}
```

**性能**：2 次 `__int128` 乘法 + 若干 64-bit 运算 ≈ 15 cycles（x86-64）。

### 2.4 加减法模块

```
前置条件：a, b < p < 2^63
不变式：a + b < 2p < 2^64（不溢出 uint64）
```

逻辑与现有相同，仅类型从 uint32 → uint64：

```
加法：r = a + b; return (r >= p) ? r - p : r
减法：return (a >= b) ? a - b : p - b + a
取负：return (a == 0) ? 0 : p - a
```

`p - b + a`：p < 2⁶³, a < 2⁶³ → p + a < 2⁶⁴，不溢出。

### 2.5 逆元模块

```
函数：inv_prime(i, p)

功能描述：计算 i 的模 p 逆元（扩展 Euclidean）。

前置条件：0 < i < p，p ≥ 2

后置条件：返回值 x 满足 i·x ≡ 1 (mod p)，x ∈ [0, p)
```

与现有算法相同，关键修复：中间乘法 `s2 * q` 用 `__int128` 防溢出。

```cpp
inline uint64_t inv_prime(uint64_t _i, uint64_t _p)
{
    uint64_t a = _p, b = _i, c;
    uint64_t s1 = 0, s2 = 1, s3;
    while ((c = a % b))
    {
        uint64_t q = a / b;
        // s2 * q 可能溢出 uint64（s2, q 均可达 p-1 ≈ 2^63）
        uint64_t sq_mod = (unsigned __int128)s2 * q % _p;
        s3 = s1 >= sq_mod ? s1 - sq_mod : _p - sq_mod + s1;
        a = b; b = c;
        s1 = s2; s2 = s3;
    }
    return s2;
}
```

注意：旧代码 `(s1 + _p - sq_mod) % _p` 对 64-bit p 有 `s1 + _p` 溢出风险，
改为条件减法避免溢出。

---

## 3. 接口规约

### 3.1 Zp 对外接口变更

| 接口 | 旧签名 | 新签名 |
|------|--------|--------|
| `number()` | `uint32_t` | `uint64_t` |
| `number()` (mutable) | `uint32_t&` | `uint64_t&` |
| `prime()` | `uint32_t` | `uint64_t` |
| `prime(p)` setter | `void(uint32_t)` | `void(uint64_t)` |
| `set_prime(p)` | `void(uint32_t)` | `void(uint64_t)` |
| `cur_prime()` | `uint32_t` | `uint64_t` |
| 构造函数 | `Zp(*, uint32_t p)` | `Zp(*, uint64_t p)` |

### 3.2 调用方责任

所有接收 `Zp::prime()` 或 `Zp::number()` 的局部变量：`uint32_t` → `uint64_t`。

### 3.3 依赖模块变更

```
接口：ZZ → Zp 构造

输入数据：ZZ 值 + uint64_t 素数
输出数据：Zp 值 ∈ [0, p)

协议约定：
  - ZZ::fdiv_ui(uint64_t d) 必须支持 64-bit 除数
  - 调用方保证 p ∈ [2, 2^63)
```

---

## 4. 改动范围

### 4.1 核心改动（number 模块）

| 文件 | 函数/类 | 改动 |
|------|--------|------|
| `number/ZZ.hh` | `fdiv_ui` | 参数 `uint32_t d` → `uint64_t d` |
| `number.hh` | `inv_prime` | 参数 `uint32_t _p` → `uint64_t _p`；中间乘法用 `__int128` |
| `number.hh` | `Zp` 类 | 数据成员重构 + Barrett 归约重写 + 所有构造/访问器改签名 |

### 4.2 下游适配（机械类型替换）

| 文件 | 改动点数 | 说明 |
|------|---------|------|
| `upolynomial.hh` | 1 | `polynomial_mod` 签名 |
| `polynomial_factorize_zp.hh` | 1 | `__make_zp` 签名 |
| `polynomial_factorize_univar.hh` | ~4 | `uint32_t p` 参数 → `uint64_t` |
| `polynomial_factorize_wang.hh` | ~15 | 参数 + `mtshl_p` 改 63-bit 素数 + `uniform_int_distribution<uint64_t>` |
| `polynomial_gcd.hh` | ~3 | `uint32_t prime` 局部变量 |
| `polynomial_gcd.cc` | ~3 | `uint32_t prime` 局部变量 |

### 4.3 `mtshl_p` 更新

```cpp
// 旧：
static constexpr uint32_t mtshl_p = 2147483629;        // 最大 < 2^31 的素数

// 新：
static constexpr uint64_t mtshl_p = 9223372036854775783ULL;  // 最大 < 2^63 的素数
```

### 4.4 测试

| 文件 | 改动 |
|------|------|
| `test/test_number.cc` | `uint32_t p` → `uint64_t`；**新增** 63-bit 素数算术测试 |
| `test/test_mtshl_*.cc` | `uint32_t p` → `uint64_t` |
| `test/test_factorize_zp.cc` | `uint32_t p` → `uint64_t` |

---

## 5. 设计决策

### D1：p < 2⁶³ 约束

加法 `a + b` 需 `2p < 2⁶⁴`（不溢出 uint64），故限制 `p < 2⁶³`。
CASC 2018 Algorithm 5 指定 63-bit 素数，此约束与论文一致。

若需支持 p ∈ [2⁶³, 2⁶⁴)，加法需进位检测（参考 FLINT `nmod_add`），
但增加复杂度且无实际需求，不做。

### D2：统一路径 vs 双路径

选择统一路径（方案 A）：所有素数用 FLINT 归一化 Barrett。

小素数（p < 2³¹）理论慢 ~1.5x（15 vs 10 cycles），但：
- 单变量因式分解瓶颈是 `powmod` 的 O(d²·log p) 算法复杂度
- 5 cycles 差异在总时间中占比极小
- 统一路径消除分支，简化维护

### D3：`_norm` 字段导致对象增大

从 24 bytes → 32 bytes（+8 bytes padding）。
Zp 作为多项式系数类型，每个单项式存一个 Zp。
对中等规模多项式（~1000 项），额外内存 ~8KB，可接受。

### D4：构造函数简化

删除 `Zp(uint32_t i, uint32_t p)` 构造函数（被 `Zp(uint64_t i, uint64_t p)` 覆盖，
`uint32_t` → `uint64_t` 隐式转换）。

保留 `Zp(int i, uint64_t p)` 和 `Zp(int64_t i, uint64_t p)` 处理负值。
