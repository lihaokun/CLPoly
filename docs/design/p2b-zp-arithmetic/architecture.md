# P2b 架构文档：Zp 标量算术优化

> 状态：待确认
> 调研依据：`docs/research/factorize-profiling.md` §3.4
> 范围限定：**只改 `Zp` 类**（`clpoly/number.hh`），不改 `upolynomial_<Zp>` 结构

---

## 1. 问题定位

`Zp` 类当前所有模运算均使用 `%` 运算符（x86 `div` 指令，约 20-40 cycles）：

```cpp
// 乘法（number.hh:152-158）
op1._i *= op2._i;
op1._i %= op1._p;   // ← div 指令

// 加法（number.hh:131-136）
op1._i += op2._i;
op1._i %= op1._p;   // ← div 指令（加法后也用 div！）
```

目标：替换为无除法的快速归约。

---

## 2. 核心流程

### 当前 `Zp` 结构（16 bytes，存在设计缺陷）

```
┌─────────────────┐
│ uint64_t _i     │  8 bytes  — 设计问题：值始终 < p < 2^32，
│                 │             却用 uint64_t 存储；且被兼作中间
│                 │             计算累加器（number() += x*y），
│                 │             导致类型职责混乱
│ uint32_t _p     │  4 bytes  — 素数 p
│ (4 bytes 填充)  │
└─────────────────┘
```

**设计缺陷**：`_i` 同时承担"归约后的值存储"和"未归约的中间计算"两个职责，导致不得不用 `uint64_t`。`number()` 的可变引用被外部代码（`polynomial_gcd.hh/cc`）直接用于 `number() += a*b` 式的 in-place 乘加，跳过了 `Zp` 的算术接口。

### 新 `Zp` 结构（16 bytes，设计正确）

```
┌─────────────────┐
│ uint32_t _i     │  4 bytes  — 值 ∈ [0, p)，类型与语义一致
│ uint32_t _p     │  4 bytes  — 素数 p
│ uint64_t _ninv  │  8 bytes  — Barrett 常数 floor(2^64 / p)
└─────────────────┘
```

总大小：16 bytes（比旧设计+Barrett 的 24 bytes 更小，且无填充）。

**原则**：中间计算全部使用局部 `uint64_t` 临时变量，`_i` 只存归约后的结果。外部代码不再通过 `number()` 做 in-place 计算，必须使用 `Zp` 算术运算符。

### 两种操作的优化策略

**加法**：a + b < 2p，用一次条件减法，无需除法：
```cpp
uint32_t r = _i + op2._i;
_i = (r >= _p) ? r - _p : r;   // 约 2-3 cycles（vs 25 cycles for %）
```

**乘法**：a, b < p < 2^32，乘积 < 2^64，用局部 `uint64_t` + Barrett 归约：
```cpp
uint64_t product = (uint64_t)_i * op2._i;
uint64_t q = (unsigned __int128)(product) * _ninv >> 64;
uint64_t r = product - q * _p;
_i = (uint32_t)((r >= _p) ? r - _p : r);   // 约 6-8 cycles（vs 25 cycles for %）
```

---

## 3. 模块划分

```
受影响文件                            改动内容
───────────────────────────────────────────────────────────────────
clpoly/number.hh                     Zp 类完全重写：
                                       - _i: uint64_t → uint32_t
                                       - 新增字段 _ninv (uint64_t)

                                       - 新增静态方法 __barrett_ninv(p)
                                       - 重写所有构造函数
                                       - 重写 operator+/-/*//= 等
                                       - number() 返回类型: uint64_t& → uint32_t&
                                       - number() const 返回: uint64_t → uint32_t

clpoly/polynomial_gcd.hh/.cc         修复 number() 滥用：
                                       将 op.number() += a*b; op.number() %= p;
                                       替换为 op = op + Zp(a,p) * Zp(b,p); 等正确调用

不改动文件
───────────────────────────────────────────────────────────────────
clpoly/upolynomial_.hh               不变（自动受益于更快的 Zp）
clpoly/polynomial_factorize_zp.hh    不变（DDF/EDF 自动受益）
clpoly/polynomial_factorize_univar.hh 不变（__select_prime 自动受益）
```

---

## 4. 关键设计决策

### 决策 0：`_i` 改为 `uint32_t`，消除 dual-role 设计缺陷

旧设计将 `_i` 既作为值存储又作为中间累加器，强迫使用 `uint64_t`。
新设计严格分离：`_i` 只存归约后的值（< p < 2^32，`uint32_t` 足够），所有中间计算用函数内局部 `uint64_t` 临时变量。`polynomial_gcd.hh/cc` 中绕过 `Zp` 接口的 `number()` 直接操作模式一并修复为正确的 `Zp` 运算。

### 决策 1：Barrett 常数存在 `Zp` 对象里（per-element），而非全局 context

**备选方案 A（选定）**：存在每个 `Zp` 对象里（per-element）
- Pro：API 完全不变，可同时存在不同素数的 `Zp` 对象（CLPoly 内部多处用到）
- Con：无（新结构 uint32_t _i + uint32_t _p + uint64_t _ninv = 16 bytes，与旧结构相同，且无填充浪费）

**备选方案 B**：全局/线程局部 context（NTL 风格，`Zp::init(p)` 设全局素数）
- Pro：元素只存值（8 bytes）
- Con：需要修改所有调用点；CLPoly 中多素数并发使用的场景（如素数选择循环中 `Zp` 对象跨越不同 p）会有全局状态冲突风险
- **结论：不选**，API 改动太大，风险高

### 决策 2：加法用条件减法（而非 Barrett）

加法结果 ∈ [0, 2p)，一次条件减即可归约到 [0, p)，无需 Barrett。
`if (_i >= _p) _i -= _p;` ≈ 2-3 cycles（远优于 `% _p` 的 20-40 cycles）。

### 决策 3：Barrett 实现使用 `unsigned __int128`

GCC/Clang 在 x86-64 上将 `unsigned __int128` 乘法编译为 `mulq` 指令（3 cycles），比调用 GMP 快得多。不引入新依赖。

### 决策 4：`_ninv` 的计算公式

对于 `uint32_t p`（p > 1），预计算：

```cpp
// ninv = UINT64_MAX / p（等价于 floor(2^64/p)，无需 __int128）
// 等价性：两常数差 δ≤1，product < p²≤2^64（严格），q 差值 < product/2^64 < 1，整数相等
static uint64_t __barrett_ninv(uint32_t p) {
    assert(p >= 2);
    return UINT64_MAX / p;
}
```

Barrett 归约正确性：
- `q = floor(product * ninv / 2^64)` 满足 `q ≤ floor(product / p) ≤ q + 1`
- 一次条件减 `if (r >= p) r -= p` 保证结果在 `[0, p)`

---

## 5. 接口规约

```
模块：Zp

功能描述：模 p 的整数算术，p 为 uint32_t 素数

前置条件（Requires）：
  - 构造时 p > 1（实际使用中均为素数；uint32_t 足够，因式分解所用素数 < 10^6 << 2^32）
  - 二元运算符要求两操作数有相同的 _p（assert 保留）

后置条件（Ensures）：
  - 所有算术操作返回值 ∈ [0, p)，且存储在 uint32_t _i 中

变更点（相对旧 API）：
  - number() 返回 uint32_t / uint32_t&（旧为 uint64_t / uint64_t&）
  - 元素大小：16 bytes（旧为 16 bytes，但现在无填充浪费，结构更紧凑）
  - polynomial_gcd.hh/cc 中 number() 直接操作改为 Zp 运算符

副作用：无（纯值语义）
```

---

## 6. 预期收益与风险

| 操作 | 当前（cycles） | 优化后（cycles） | 倍率 |
|------|---------------|-----------------|------|
| `Zp operator+` | ~25（div） | ~3（条件减） | **8x** |
| `Zp operator*` | ~25（div） | ~8（Barrett） | **3x** |
| `Zp operator/` | ~50（inv+div） | ~35（inv+Barrett） | ~1.4x |
| DDF per trial（W15） | ~0.133ms | ~0.04ms | **~3x** |
| sel_prime 总时间（30 trials） | ~4.0ms | ~1.2ms | **~3x** |

**P2b 单独实施**：总时间从 4.5ms 降到约 1.5ms（W15）。
**P2a + P2b 叠加**（3 次试验 + 快速算术）：预期降到 ~0.15ms，vs FLINT 0.63ms（剩余差距来自 upolynomial 数据结构，暂不优化）。

**风险**：
- `unsigned __int128` 在 GCC/Clang x86-64 上有效；如需移植到 MSVC 或 ARM 需改用 `__umulh`
- `polynomial_gcd.hh/cc` 中 `number()` 滥用点需逐一修复并回归测试
- Barrett 正确性边界：`p < 2^32`，CLPoly 始终满足（`uint32_t` prime）
