# P3-mul 细化文档：dense_upoly_zp 快速算术

> 状态：待确认
> 架构依据：`docs/design/p3-mul/architecture.md`

---

## M1: Karatsuba 乘法

### M1.1 `_classical_mul` — schoolbook 基例（惰性累加）

```
函数签名：
  void _classical_mul(uint64_t* C,
                      const uint64_t* A, size_t len_a,
                      const uint64_t* B, size_t len_b) const;

功能描述：计算 C[0..len_a+len_b-2] = A × B mod _p（schoolbook，点积形式）。
         使用 3-word 惰性累加，每输出系数仅 1 次模归约。
         与 FLINT _NMOD_VEC_DOT3 同等效率。

前置条件（Requires）：
  - len_a > 0, len_b > 0
  - C 预分配长度 >= len_a + len_b - 1
  - C 与 A, B 无重叠
  （C 无需预清零——每个 C[k] 由归约直接赋值）

后置条件（Ensures）：
  - C[k] = sum_{i+j=k} A[i]·B[j] mod _p，对 k = 0..len_a+len_b-2

副作用：写入 C[0..len_a+len_b-2]
```

**实现**：

```cpp
void _classical_mul(uint64_t* C,
                    const uint64_t* A, size_t len_a,
                    const uint64_t* B, size_t len_b) const
{
    size_t len_c = len_a + len_b - 1;
    for (size_t k = 0; k < len_c; ++k) {
        word3 acc = {0, 0, 0};
        size_t j_min = (k >= len_b) ? k - len_b + 1 : 0;
        size_t j_max = (k < len_a) ? k : len_a - 1;
        for (size_t j = j_min; j <= j_max; ++j) {
            uint64_t p1, p0;
            _umul128(p1, p0, A[j], B[k - j]);   // 1 次 mulq
            _add_carry3(acc, p1, p0);             // adcq 链
        }
        C[k] = _lll_mod_preinv(acc.hi, acc.mid, acc.lo,
                                _p, _ninv, _norm);  // 每系数 1 次归约
    }
}
```

**与旧设计的差异**：

| 指标 | 旧（逐元素 Barrett） | 新（惰性累加） |
|------|---------------------|---------------|
| 内循环乘法 | `nmod_mul`：2 次 mulq + 归约 | `_umul128`：1 次 mulq |
| 内循环加法 | `nmod_add`：分支比较 | `_add_carry3`：无分支 adcq |
| 模归约次数 | len\_a × len\_b | len\_c = len\_a + len\_b - 1 |
| C 预清零 | 需要 | 不需要 |

**复用**：`word3`、`_umul128`、`_add_carry3`、`_lll_mod_preinv` 全部复用 M2.0 辅助函数，零新增代码。

**参考**：FLINT `_nmod_poly_mul_classical` 对 64-bit 素数使用 `_NMOD_VEC_DOT3` 宏，内部结构完全相同（`umul_ppmm` + `add_sssaaaaaa` + `NMOD_RED3`）。

---

### M1.2 `_kar_mul` — Karatsuba 递归核心

```
函数签名：
  void _kar_mul(uint64_t* C,
                const uint64_t* A, const uint64_t* B,
                size_t n, uint64_t* scratch) const;

功能描述：计算 C[0..2n-2] = A[0..n-1] × B[0..n-1] mod _p。
         n < KARATSUBA_THRESHOLD 时调用 _classical_mul。

前置条件（Requires）：
  - n > 0
  - A, B 长度均为 n（调用方已零填充）
  - C 预分配长度 >= 2n-1
  - scratch 预分配长度 >= 6n
  - C 与 A, B, scratch 无重叠

后置条件（Ensures）：
  - C[k] = sum_{i+j=k} A[i]·B[j] mod _p，对 k = 0..2n-2

副作用：写入 C[0..2n-2]；scratch 内容未定义
```

**实现**：

```cpp
void _kar_mul(uint64_t* C,
              const uint64_t* A, const uint64_t* B,
              size_t n, uint64_t* scratch) const
{
    if (n < KARATSUBA_THRESHOLD) {
        _classical_mul(C, A, n, B, n);   // 惰性累加，无需预清零
        return;
    }

    size_t m = n / 2;
    size_t h = n - m;   // h = m 或 m+1

    // --- scratch 布局 ---
    // t1[h] | t2[h] | P0[2m-1] | P1[2h-1] | rec_scratch[...]
    uint64_t* t1  = scratch;
    uint64_t* t2  = t1 + h;
    uint64_t* sP0 = t2 + h;
    uint64_t* sP1 = sP0 + (2 * m - 1);
    uint64_t* rec = sP1 + (2 * h - 1);

    // --- t1 = A_lo + A_hi, t2 = B_lo + B_hi ---
    // A_lo 有 m 项, A_hi 有 h 项; h >= m
    for (size_t i = 0; i < m; ++i) {
        t1[i] = nmod_add(A[i], A[m + i]);
        t2[i] = nmod_add(B[i], B[m + i]);
    }
    if (h > m) {   // n 为奇数
        t1[m] = A[m + m];   // A_lo[m] = 0, 所以 t1[m] = A_hi[m]
        t2[m] = B[m + m];
    }

    // --- P0 = A_lo × B_lo → sP0[0..2m-2] ---
    _kar_mul(sP0, A, B, m, rec);

    // --- P1 = t1 × t2 → sP1[0..2h-2] ---
    _kar_mul(sP1, t1, t2, h, rec);

    // --- P2 = A_hi × B_hi → C[2m..2n-2] ---
    _kar_mul(C + 2 * m, A + m, B + m, h, rec);

    // --- P1 -= P0 ---
    for (size_t i = 0; i < 2 * m - 1; ++i)
        sP1[i] = nmod_sub(sP1[i], sP0[i]);

    // --- P1 -= P2 ---
    for (size_t i = 0; i < 2 * h - 1; ++i)
        sP1[i] = nmod_sub(sP1[i], C[2 * m + i]);

    // --- 组装 C ---
    // C[0..2m-2] = P0
    std::memcpy(C, sP0, (2 * m - 1) * sizeof(uint64_t));

    // C[2m-1] = 0（间隙）
    C[2 * m - 1] = 0;

    // C[m..m+2h-2] += P1（交叉项）
    for (size_t i = 0; i < 2 * h - 1; ++i)
        C[m + i] = nmod_add(C[m + i], sP1[i]);
}
```

**调用关系**：`_kar_mul` → `_kar_mul`（递归）/ `_classical_mul`（基例）

**关键实现细节**：

1. **无需预清零**：`_classical_mul` 使用惰性累加（点积形式），直接赋值每个 C[k]，无需 C 预清零。递归路径的组装步骤（memcpy P0 + 间隙清零 + P1 加法）覆盖 C 的所有位置。因此递归调用前也无需 memset。
2. **rec 共享**：P0、P1、P2 的递归 scratch 共享同一块 `rec` 区域。安全性：P0 完成后结果在 sP0，rec 可被 P1 覆盖；P1 完成后结果在 sP1，rec 可被 P2 覆盖。
3. **n 为奇数时**：h = m+1，t1/t2 最后一项直接取 A\_hi/B\_hi 的最后一项（因为 A\_lo 只有 m 项，缺失位用 0 填充）。

---

### M1.3 修改 `mul` — 添加 Karatsuba 分派

```
函数签名：（不变）
  static void mul(dense_upoly_zp& C,
                  const dense_upoly_zp& A, const dense_upoly_zp& B);

修改内容：在 aliasing 保护和空输入检查之后，根据 max(len_a, len_b)
         决定 schoolbook 或 Karatsuba 路径。
```

**实现**：

```cpp
static void mul(dense_upoly_zp& C,
                const dense_upoly_zp& A, const dense_upoly_zp& B)
{
    assert(A._p == B._p);
    // aliasing 保护（不变）
    if (&C == &A || &C == &B) {
        dense_upoly_zp tmp;
        mul(tmp, A, B);
        C = std::move(tmp);
        return;
    }

    C._p = A._p; C._ninv = A._ninv; C._norm = A._norm;
    if (A.empty() || B.empty()) { C._coeffs.clear(); return; }

    size_t len_a = A._coeffs.size(), len_b = B._coeffs.size();
    size_t n = std::max(len_a, len_b);

    if (n < KARATSUBA_THRESHOLD) {
        // schoolbook
        C._coeffs.resize(len_a + len_b - 1);   // 无需清零，_classical_mul 直接赋值
        C._classical_mul(C._coeffs.data(),
                         A._coeffs.data(), len_a,
                         B._coeffs.data(), len_b);
        C.__strip();
    } else {
        // Karatsuba: 零填充到等长 n
        std::vector<uint64_t> a_pad(n, 0), b_pad(n, 0);
        std::memcpy(a_pad.data(), A._coeffs.data(), len_a * sizeof(uint64_t));
        std::memcpy(b_pad.data(), B._coeffs.data(), len_b * sizeof(uint64_t));

        C._coeffs.resize(2 * n - 1);
        std::vector<uint64_t> scratch(6 * n);
        C._kar_mul(C._coeffs.data(), a_pad.data(), b_pad.data(),
                   n, scratch.data());
        C.__strip();
    }
}
```

**复用**：aliasing 保护、空输入检查、`_p/_ninv/_norm` 复制均复用现有代码。

**与现有代码的差异**：现有 `mul`（第 148-175 行）的内循环逻辑移入 `_classical_mul`，`mul` 本身变为分派器。

---

### M1.4 新增常量

```cpp
// dense_upoly_zp 类 private 区域
static constexpr size_t KARATSUBA_THRESHOLD = 16;
```

后续 A/B 调优范围：12-24。

---

## M2: 惰性归约 divrem

### M2.0 辅助类型与内联函数

**`word3` 结构体**（类 private）：

```cpp
struct word3 { uint64_t lo, mid, hi; };
```

**`_umul128`** — 128-bit 无符号乘法：

```
函数签名：
  static void _umul128(uint64_t& hi, uint64_t& lo, uint64_t a, uint64_t b);

功能描述：计算 a × b 的 128-bit 乘积，hi:lo = a * b。
         替代 FLINT 宏 umul_ppmm。

前置条件：无
后置条件：(hi << 64) | lo == a * b（数学精确）
副作用：写入 hi, lo
```

**实现**：

```cpp
static void _umul128(uint64_t& hi, uint64_t& lo, uint64_t a, uint64_t b)
{
    unsigned __int128 prod = (unsigned __int128)a * b;
    lo = (uint64_t)prod;
    hi = (uint64_t)(prod >> 64);
}
```

编译器生成：单条 `mulq` 指令，与 FLINT `umul_ppmm` 完全等价。

---

**`_add_carry3`** — 3-word 进位加法：

```
函数签名：
  static void _add_carry3(word3& s, uint64_t b1, uint64_t b0);

功能描述：将 128-bit 值 (b1:b0) 加到 192-bit 累加器 s 上。
         替代 FLINT 宏 add_sssaaaaaa。

前置条件：s.hi + carry <= 2^64-1（由溢出分析保证，见架构文档 §3）
后置条件：s_new = s_old + (b1 << 64) + b0（数学精确）
副作用：修改 s
```

**实现**（平台分派 inline asm + 可移植 fallback）：

```cpp
static void _add_carry3(word3& s, uint64_t b1, uint64_t b0)
{
#if defined(__x86_64__) || defined(_M_X64)
    __asm__(
        "addq %[b0], %[lo]\n\t"
        "adcq %[b1], %[mid]\n\t"
        "adcq $0, %[hi]"
        : [lo] "+&r"(s.lo), [mid] "+&r"(s.mid), [hi] "+r"(s.hi)
        : [b0] "rme"(b0), [b1] "rme"(b1)
        : "cc"
    );
#elif defined(__aarch64__) || defined(_M_ARM64)
    __asm__(
        "adds %[lo], %[lo], %[b0]\n\t"
        "adcs %[mid], %[mid], %[b1]\n\t"
        "adc  %[hi], %[hi], xzr"
        : [lo] "+&r"(s.lo), [mid] "+&r"(s.mid), [hi] "+r"(s.hi)
        : [b0] "r"(b0), [b1] "r"(b1)
        : "cc"
    );
#else
    // __int128 fallback（可移植，编译器通常生成 ~16 条指令）
    unsigned __int128 sum = (unsigned __int128)s.lo + b0;
    s.lo = (uint64_t)sum;
    sum = (unsigned __int128)s.mid + b1 + (uint64_t)(sum >> 64);
    s.mid = (uint64_t)sum;
    s.hi += (uint64_t)(sum >> 64);
#endif
}
```

**平台 codegen**：

| 平台 | 实现 | 核心指令 | 指令数 |
|------|------|---------|--------|
| x86_64 | inline asm | `addq` + `adcq` + `adcq` | 3 |
| AArch64 | inline asm | `adds` + `adcs` + `adc` | 3 |
| 其他 | `__int128` fallback | 编译器生成 | ~16 |

与 FLINT `add_sssaaaaaa` 完全等价。

**Earlyclobber 约束**：`lo` 和 `mid` 使用 `"+&r"`（earlyclobber），防止编译器将输入 `b1` 分配到与 `lo` 相同的寄存器——`addq` 先写 `lo`，若 `b1` 与 `lo` 同寄存器则后续 `adcq` 读到污染值。`hi` 最后写入，不需要 `&`。参考 FLINT `add_sssaaaaaa` 同样使用 `=&r`。

**为什么不用纯 C 比较链**：C 语言无法表达 CPU 进位标志（carry flag），编译器无法将 `(s.lo < old_lo)` 等溢出检测优化为 `adcq`。实测 C 比较链生成 ~18 条指令（含 `setc`/`sete`/`andl`/`orl`），在 deg199 内循环中比 inline asm 慢 68%（0.039ms vs 0.023ms）。

---

### M2.1 `_lll_mod_preinv` — 3-word 模归约

```
函数签名：
  static uint64_t _lll_mod_preinv(uint64_t hi, uint64_t mid, uint64_t lo,
                                   uint64_t p, uint64_t pinv, uint32_t norm);

功能描述：将非负 192-bit 值 (hi·2^128 + mid·2^64 + lo) 归约为 r ∈ [0, p)。

前置条件（Requires）：
  - p >= 2
  - pinv = __preinvert_limb(p << norm)
  - norm = __builtin_clzll(p)
  - hi < p（Barrett 归一化要求 h_shifted < pn = p << norm）
  - (hi, mid, lo) 表示非负 192-bit 值

后置条件（Ensures）：
  - 返回值 r = (hi·2^128 + mid·2^64 + lo) mod p
  - r ∈ [0, p)

副作用：无
```

**实现**（参考 FLINT `n_lll_mod_preinv`）：

```cpp
static uint64_t _lll_mod_preinv(uint64_t hi, uint64_t mid, uint64_t lo,
                                 uint64_t p, uint64_t pinv, uint32_t norm)
{
    // 第一步：归约 (hi, mid) → r1 ∈ [0, p)
    // 等价于 (hi·2^64 + mid) mod p
    uint64_t r1;
    {
        uint64_t pn = p << norm;
        uint64_t h_shifted = hi << norm;
        // 补上 mid 的高 norm 位
        if (norm > 0) h_shifted |= (mid >> (64 - norm));
        uint64_t m_shifted = mid << norm;

        unsigned __int128 qm = (unsigned __int128)h_shifted * pinv;
        uint64_t q1 = (uint64_t)(qm >> 64);
        uint64_t q0 = (uint64_t)qm;
        q0 += m_shifted;
        q1 += h_shifted + (q0 < m_shifted ? 1 : 0);
        uint64_t r = m_shifted - (q1 + 1) * pn;
        if (r > q0) r += pn;
        if (r >= pn) r -= pn;
        r1 = r >> norm;
    }

    // 第二步：归约 (r1, lo) → r2 ∈ [0, p)
    // 等价于 (r1·2^64 + lo) mod p
    // 复用现有 nmod_mul 的 Barrett 结构，但输入是 (r1, lo) 而非 (a*b)
    {
        uint64_t pn = p << norm;
        uint64_t h = r1 << norm;
        if (norm > 0) h |= (lo >> (64 - norm));
        uint64_t l = lo << norm;

        unsigned __int128 qm = (unsigned __int128)h * pinv;
        uint64_t q1 = (uint64_t)(qm >> 64);
        uint64_t q0 = (uint64_t)qm;
        q0 += l;
        q1 += h + (q0 < l ? 1 : 0);
        uint64_t r = l - (q1 + 1) * pn;
        if (r > q0) r += pn;
        if (r >= pn) r -= pn;
        return r >> norm;
    }
}
```

**复用**：每步归约的 Barrett 逻辑与现有 `nmod_mul` 第 72-86 行结构相同，区别在于输入是 (hi, lo) 两个独立 word 而非乘积 a*b。

---

### M2.2 重写 `divrem` — 惰性归约版

```
函数签名：（不变）
  static void divrem(dense_upoly_zp& Q, dense_upoly_zp& R,
                     const dense_upoly_zp& A, const dense_upoly_zp& B);

修改内容：替换内层循环为 3-word 累加器 + negmod + add。
         外层框架（aliasing 保护、deg 检查、Q 初始化）不变。
```

**实现**：

```cpp
static void divrem(dense_upoly_zp& Q, dense_upoly_zp& R,
                   const dense_upoly_zp& A, const dense_upoly_zp& B)
{
    assert(!B.empty() && A._p == B._p);
    // aliasing 保护（不变）
    if (&Q == &A || &Q == &B || &R == &A || &R == &B) {
        dense_upoly_zp tq, tr;
        divrem(tq, tr, A, B);
        Q = std::move(tq);
        R = std::move(tr);
        return;
    }

    Q._p = A._p; Q._ninv = A._ninv; Q._norm = A._norm;
    R._p = A._p; R._ninv = A._ninv; R._norm = A._norm;

    if (A.deg() < B.deg()) {
        Q._coeffs.clear();
        R._coeffs = A._coeffs;
        return;
    }

    int64_t d = B.deg();
    int64_t q_len = A.deg() - d + 1;
    Q._coeffs.assign(q_len, 0);
    uint64_t inv_lc = R.nmod_inv(B.lead());

    // --- 3-word 累加器 ---
    size_t r_len = A._coeffs.size();
    std::vector<word3> R3(r_len);
    for (size_t i = 0; i < r_len; ++i)
        R3[i] = {A._coeffs[i], 0, 0};

    uint64_t p = A._p;
    uint64_t pinv = A._ninv;
    uint32_t norm = A._norm;

    // --- 外层循环 ---
    for (int64_t i = q_len - 1; i >= 0; --i) {
        // 归约首项
        uint64_t r = _lll_mod_preinv(R3[i + d].hi, R3[i + d].mid,
                                      R3[i + d].lo, p, pinv, norm);
        uint64_t q_i = R.nmod_mul(r, inv_lc);
        Q._coeffs[i] = q_i;

        if (q_i != 0) {
            uint64_t c = p - q_i;   // negmod

            // 内层循环：逐系数 3-word 加法
            for (int64_t j = 0; j <= d; ++j) {
                uint64_t p1, p0;
                _umul128(p1, p0, c, B._coeffs[j]);
                _add_carry3(R3[i + j], p1, p0);
            }
        }
    }

    // --- 最终归约：仅余式部分 ---
    R._coeffs.resize(d > 0 ? d : 0);
    for (int64_t i = 0; i < d; ++i)
        R._coeffs[i] = _lll_mod_preinv(R3[i].hi, R3[i].mid, R3[i].lo,
                                        p, pinv, norm);
    R.__strip();
    Q.__strip();
}
```

**复用**：
- aliasing 保护：现有 divrem 的 aliasing 逻辑不变
- deg 检查 + Q 初始化：现有逻辑不变
- `__strip()` / `nmod_mul` / `nmod_inv`：现有方法
- Barrett 预计算参数 `_ninv`, `_norm`：现有成员

**关键实现细节**：

1. **`_umul128` + `_add_carry3` 封装**：将 FLINT 宏 `umul_ppmm` / `add_sssaaaaaa` 的等价逻辑封装为类的 private static inline 函数。类型安全，可调试。`_umul128` 使用 `__int128` 编译器自动生成 `mulq`；`_add_carry3` 使用平台分派 inline asm 保证生成最优的 `addq` + `adcq` 链（详见 M2.0）。内循环仅 3 行，可读性好。

2. **`word3` 结构体**：类 private 结构体，比 `uint64_t R3[3*r_len]` 可读性更好，编译器生成相同代码。

3. **R3 大小为 `r_len = deg(A)+1`**：涵盖所有被读写的位置。位置 i+d 的最大值 = (q_len-1)+d = deg(A)，在范围内。

---

## 错误处理策略

| 场景 | 处理 | 位置 |
|------|------|------|
| A._p ≠ B._p | `assert` 失败（debug 构建）| mul, divrem 入口 |
| B 为空（除零）| `assert` 失败 | divrem 入口 |
| A, B 均为空 | `mul` 返回空多项式，`divrem` 不会到达此分支 | mul 空输入检查 |
| aliasing (C == A) | 临时对象递归调用 | mul/divrem 入口 |
| n = 0 (Karatsuba) | 不可能：由 `mul` 空输入检查保证 | — |
| n = 1 (Karatsuba) | 走 `_classical_mul` 基例 | _kar_mul |
| deg(A) < deg(B) | Q = 0, R = A | divrem 度数检查 |

所有错误条件均由 `assert` 保护（与现有代码风格一致）。不引入异常或错误码。

---

## 实施步骤

### 步骤 1：M2.0 共享辅助函数

1a. 添加 `word3`、`_umul128`、`_add_carry3`、`_lll_mod_preinv`（M1、M2 共用）
1b. 编译验证

### 步骤 2：M1 Karatsuba mul

2a. 添加 `KARATSUBA_THRESHOLD`、`_classical_mul`（惰性累加，依赖步骤 1 的辅助函数）
2b. 添加 `_kar_mul`
2c. 修改 `mul` 为分派逻辑
2d. 运行 `bash test/run_all_tests.sh` 验证正确性
2e. `make bench-clpoly` 对比 mul 性能

### 步骤 3：M2 惰性归约 divrem

3a. 重写 `divrem` 内循环（依赖步骤 1 的辅助函数）
3b. 运行 `bash test/run_all_tests.sh` 验证正确性
3c. `make bench-clpoly` 对比 divrem/gcd 性能

### 步骤 4：全量验证

4a. `make crosscheck` — FLINT 交叉验证
4b. `make bench-comparative` — FLINT 对比性能
4c. 记录 benchmark 结果
