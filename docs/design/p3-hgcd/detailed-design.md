# P3-HGCD 细化设计

> 状态：待确认
> 依据：`architecture.md`（架构）、`correctness-proof.md`（正确性推导）
> 实现文件：`clpoly/dense_upoly_zp.hh`

---

## 1. 现有代码结构

`dense_upoly_zp` 当前 private 区域包含：

```
word3, _umul128, _add_carry3, _lll_mod_preinv    // 3-word 惰性归约
KARATSUBA_THRESHOLD = 16
_classical_mul(C, A, len_a, B, len_b)            // schoolbook（点积形式）
_kar_mul(C, A, B, n, scratch)                     // Karatsuba（等长，需 scratch）
nmod_add, nmod_sub                                // 模加减（private）
__strip, __precompute, __preinvert_limb           // 辅助
```

public 区域：

```
nmod_mul, nmod_inv                                // 模乘、模逆
构造、查询、to_upoly、scalar_mul
mul(C, A, B)                                      // 对象级乘法
divrem(Q, R, A, B)                                // 对象级除法
gcd(G, A, B)                                      // 对象级 GCD（当前纯 Euclid）
```

新增代码全部在 `dense_upoly_zp` 类内部，按模块顺序添加。

---

## 2. 常量与辅助类型

```cpp
// === 新增常量 ===
static constexpr size_t GCD_CUTOFF  = 340;   // gcd: Euclid vs HGCD
static constexpr size_t HGCD_CUTOFF = 100;   // HGCD: 迭代 vs 递归

// === 新增辅助结构 ===
struct hgcd_mat {
    uint64_t* poly[4];   // [0]=M00, [1]=M01, [2]=M10, [3]=M11
    size_t len[4];        // 各多项式的 length（非 degree）
};
```

`hgcd_mat` 的 4 个指针指向工作区内的预分配空间，不独立分配。每个多项式槽最大 `(n+1)/2` 元素。

---

## 3. M0: 基础操作

### M0.1: `_poly_normalise`

```cpp
// 返回去除高位零后的 length
static size_t _poly_normalise(const uint64_t* A, size_t len)
{
    while (len > 0 && A[len - 1] == 0)
        --len;
    return len;
}
```

### M0.2: `_poly_add` / `_poly_sub`

```cpp
// 逐系数模加，C 预分配 max(len_a, len_b) 元素
// 支持 aliasing：C 可与 A 或 B 重叠
// 返回 normalize 后的 length
size_t _poly_add(uint64_t* C,
                 const uint64_t* A, size_t len_a,
                 const uint64_t* B, size_t len_b) const
{
    size_t i, min_len = std::min(len_a, len_b);
    size_t max_len = std::max(len_a, len_b);

    for (i = 0; i < min_len; ++i)
        C[i] = nmod_add(A[i], B[i]);

    // 较长的部分直接拷贝（若 C 已与较长者重叠则跳过）
    if (len_a > len_b) {
        if (C != A)
            std::memcpy(C + min_len, A + min_len,
                        (len_a - min_len) * sizeof(uint64_t));
    } else if (len_b > len_a) {
        if (C != B)
            std::memcpy(C + min_len, B + min_len,
                        (len_b - min_len) * sizeof(uint64_t));
    }

    return _poly_normalise(C, max_len);
}

// 逐系数模减，C = A - B
// 支持 aliasing：C 可与 A 或 B 重叠
size_t _poly_sub(uint64_t* C,
                 const uint64_t* A, size_t len_a,
                 const uint64_t* B, size_t len_b) const
{
    size_t i, min_len = std::min(len_a, len_b);
    size_t max_len = std::max(len_a, len_b);

    for (i = 0; i < min_len; ++i)
        C[i] = nmod_sub(A[i], B[i]);

    if (len_a > len_b) {
        // C[min..len_a-1] = A[min..len_a-1] - 0
        if (C != A)
            std::memcpy(C + min_len, A + min_len,
                        (len_a - min_len) * sizeof(uint64_t));
    } else if (len_b > len_a) {
        // C[min..len_b-1] = 0 - B[min..len_b-1] = p - B[i]
        for (i = min_len; i < len_b; ++i)
            C[i] = (B[i] == 0) ? 0 : _p - B[i];
    }

    return _poly_normalise(C, max_len);
}
```

### M0.3: `_poly_divrem`（raw API，惰性 3-word）

复用现有 `divrem` 的 3-word 累加器逻辑，改为 raw pointer 接口。

```cpp
// 完整除法：
//   len_a >= len_b 时：Q 预分配 len_a - len_b + 1，R 预分配 len_b - 1
//   len_a < len_b 时：Q 不写入（len_q=0），R 预分配 len_a
// W3: 预分配 len_a 个 word3 的工作区（3-word 累加器，由调用方提供）
// 返回 {len_q, len_r}，均已 normalize
// 前置：len_b > 0
// 对标 FLINT __divrem 宏：len_a < len_b 时退化为 Q=0, R=A
std::pair<size_t, size_t>
_poly_divrem(uint64_t* Q, uint64_t* R,
             const uint64_t* A, size_t len_a,
             const uint64_t* B, size_t len_b,
             word3* W3) const
{
    assert(len_b > 0);

    // FLINT 兼容：len_a < len_b 时 Q=0, R=A（不执行除法）
    // HGCD 递归重构后 lena2 < lenb2 理论上可能发生（Zp 首项抵消）
    if (len_a < len_b) {
        std::memcpy(R, A, len_a * sizeof(uint64_t));
        return {0, len_a};
    }

    size_t d = len_b - 1;       // deg(B)
    size_t q_len = len_a - d;   // len_a - len_b + 1
    uint64_t inv_lc = nmod_inv(B[d]);

    // 3-word 累加器（使用调用方提供的工作区，零动态分配）
    for (size_t i = 0; i < len_a; ++i)
        W3[i] = {A[i], 0, 0};

    for (size_t ii = q_len; ii > 0; --ii) {
        size_t i = ii - 1;      // 从 q_len-1 downto 0
        // 归约首项
        uint64_t r = _lll_mod_preinv(W3[i + d].hi, W3[i + d].mid,
                                      W3[i + d].lo, _p, _ninv, _norm);
        uint64_t q_i = nmod_mul(r, inv_lc);
        Q[i] = q_i;

        if (q_i != 0) {
            uint64_t c = _p - q_i;   // negmod
            for (size_t j = 0; j <= d; ++j) {
                uint64_t p1, p0;
                _umul128(p1, p0, c, B[j]);
                _add_carry3(W3[i + j], p1, p0);
            }
        }
    }

    // 最终归约余式部分
    for (size_t i = 0; i < d; ++i)
        R[i] = _lll_mod_preinv(W3[i].hi, W3[i].mid, W3[i].lo,
                                _p, _ninv, _norm);

    size_t len_q = _poly_normalise(Q, q_len);
    size_t len_r = _poly_normalise(R, d);
    return {len_q, len_r};
}
```

**零动态分配**：`word3* W3` 由调用方从预分配工作区提供。HGCD 中所有 `_poly_divrem` 调用共享同一递归层的 W3 空间（同层的多次 divrem 不并发），零 heap 操作。

对象级 `divrem` 仍使用内部 `std::vector<word3>` 分配（非 HGCD 热路径，保持简洁）。

### M0.4: `_mul`（不等长 raw API）

```cpp
// 不等长乘法
// C 预分配 2*len_a - 1 元素（Karatsuba 零填充后写满 2n-1，
//   其中 C[len_a+len_b-1 .. 2*len_a-2] 数学上为零但会被写入）
// schoolbook 分支只写 len_a + len_b - 1 元素
// scratch 预分配 7*len_a（b_pad[len_a] + kar_scratch[6*len_a]）
// 前置：len_a >= len_b > 0
void _mul(uint64_t* C,
          const uint64_t* A, size_t len_a,
          const uint64_t* B, size_t len_b,
          uint64_t* scratch) const
{
    assert(len_a >= len_b && len_b > 0);

    if (len_b < KARATSUBA_THRESHOLD) {
        // schoolbook（短操作数 → O(len_a * len_b)，写 len_a+len_b-1 元素）
        _classical_mul(C, A, len_a, B, len_b);
    } else {
        // 零填充到等长后 Karatsuba
        size_t n = len_a;
        uint64_t* b_pad = scratch;
        uint64_t* kar_scratch = scratch + n;
        std::memcpy(b_pad, B, len_b * sizeof(uint64_t));
        std::memset(b_pad + len_b, 0, (n - len_b) * sizeof(uint64_t));

        // _kar_mul 写入 C[0 .. 2n-2]（共 2n-1 元素）
        // 其中 C[len_a+len_b-1 .. 2n-2] 数学上为零（B 高位全零），
        // 但 Karatsuba 递归的中间加减可能写入非零后再消回零，
        // 所以调用方必须确保 C 至少 2n-1 = 2*len_a-1 元素
        _kar_mul(C, A, b_pad, n, kar_scratch);
    }
}
```

**空间需求**：
- **C 空间**：`2*len_a - 1` 元素（Karatsuba 分支）。schoolbook 分支只需 `len_a + len_b - 1`，但统一按 `2*len_a - 1` 分配以避免分支判断。
- **scratch**：`7*len_a`（`b_pad[len_a]` + `kar_scratch[6*len_a]`）。schoolbook 分支不使用 scratch。

**安全性**：HGCD 中所有 `_mul` 调用点的 C 空间均已验证足够（详见正确性验证记录）——`_mat_row_update` 的 *pT 为 n 大小，`_mat_mul_entry` 的 C 为父层 half ≈ 2*子层 half。

**不等长策略**：当 `len_a > len_b >= KARATSUBA_THRESHOLD` 时，零填充 B 到 len_a 长度后做 Karatsuba。虽然这浪费了一些（因为 B 高位全零），但 Karatsuba 仍然是 O(n^1.585)，且实现简单。后续可用分段 Karatsuba 优化不等长场景，但第一版不需要。

### M0.5: 对象级 `add` / `sub`

```cpp
static void add(dense_upoly_zp& C,
                const dense_upoly_zp& A, const dense_upoly_zp& B)
{
    assert(A._p == B._p);
    // aliasing 保护
    if (&C == &A || &C == &B) {
        dense_upoly_zp tmp;
        add(tmp, A, B);
        C = std::move(tmp);
        return;
    }
    C._p = A._p; C._ninv = A._ninv; C._norm = A._norm;
    size_t max_len = std::max(A._coeffs.size(), B._coeffs.size());
    C._coeffs.resize(max_len);
    size_t len = C._poly_add(C._coeffs.data(),
                              A._coeffs.data(), A._coeffs.size(),
                              B._coeffs.data(), B._coeffs.size());
    C._coeffs.resize(len);
}

static void sub(dense_upoly_zp& C,
                const dense_upoly_zp& A, const dense_upoly_zp& B)
{
    assert(A._p == B._p);
    if (&C == &A || &C == &B) {
        dense_upoly_zp tmp;
        sub(tmp, A, B);
        C = std::move(tmp);
        return;
    }
    C._p = A._p; C._ninv = A._ninv; C._norm = A._norm;
    size_t max_len = std::max(A._coeffs.size(), B._coeffs.size());
    C._coeffs.resize(max_len);
    size_t len = C._poly_sub(C._coeffs.data(),
                              A._coeffs.data(), A._coeffs.size(),
                              B._coeffs.data(), B._coeffs.size());
    C._coeffs.resize(len);
}
```

---

## 4. M1: HGCD 迭代 base case

### 4.1 函数签名

```cpp
// 迭代 Euclid + 矩阵追踪
// 输入：a[len_a], b[len_b]，前置 len_a > len_b > 0
// 输出：M（矩阵），*pA/*pB（缩减后多项式的指针），len_A/len_B（长度）
// 缓冲区：
//   pA, pB, pT: 3 个多项式缓冲区（各 len_a 大小），通过指针旋转复用
//   M.poly[0..3], pt: 5 个矩阵项缓冲区（各 (len_a+1)/2 大小）
//   Q: 商缓冲区（(len_a+1)/2 大小）
//   W3: word3 累加器（len_a 个 word3，_poly_divrem 使用）
//   scratch: _mul 所需工作区
// 返回：sgn（行列式符号）
int _hgcd_iter(hgcd_mat& M,
               uint64_t** pA, size_t& len_A,
               uint64_t** pB, size_t& len_B,
               const uint64_t* a, size_t len_a,
               const uint64_t* b, size_t len_b,
               uint64_t* Q, word3* W3,
               uint64_t** pT, uint64_t** pt,
               uint64_t* scratch) const
```

### 4.2 完整实现

```cpp
int _hgcd_iter(hgcd_mat& M,
               uint64_t** pA, size_t& len_A,
               uint64_t** pB, size_t& len_B,
               const uint64_t* a, size_t len_a,
               const uint64_t* b, size_t len_b,
               uint64_t* Q, word3* W3,
               uint64_t** pT, uint64_t** pt,
               uint64_t* scratch) const
{
    const size_t m = len_a / 2;
    int sgn = 1;

    // 初始化矩阵为 identity
    _mat_one(M);

    // 初始化 A = a, B = b（拷贝到缓冲区）
    // 注意：A ← a 必须在 B ← b 之前（支持 {*pB, a} aliasing）
    std::memcpy(*pA, a, len_a * sizeof(uint64_t));
    len_A = len_a;
    std::memcpy(*pB, b, len_b * sizeof(uint64_t));
    len_B = len_b;

    while (len_B >= m + 1) {
        // divrem: *pT = *pA mod *pB, Q = *pA div *pB
        size_t len_Q, len_T;
        std::tie(len_Q, len_T) = _poly_divrem(Q, *pT, *pA, len_A, *pB, len_B, W3);

        // 指针旋转：A_new = B_old, B_new = T (= A_old mod B_old)
        // T_new = A_old（可被下次覆盖）
        std::swap(*pB, *pT);       // *pB 现在指向余式
        std::swap(*pA, *pT);       // *pA 现在指向旧 *pB
        len_A = len_B;
        len_B = len_T;

        // 矩阵更新：M_new = M_old * [[Q, 1], [1, 0]]
        // 行 2：M[2]_new = M[3]_old + Q * M[2]_old, M[3]_new = M[2]_old
        if (len_Q > 0 && M.len[2] > 0) {
            size_t len_prod;
            if (len_Q >= M.len[2])
                _mul(*pT, Q, len_Q, M.poly[2], M.len[2], scratch);
            else
                _mul(*pT, M.poly[2], M.len[2], Q, len_Q, scratch);
            len_prod = len_Q + M.len[2] - 1;

            size_t lent = _poly_add(*pt, M.poly[3], M.len[3], *pT, len_prod);
            std::swap(M.poly[3], M.poly[2]);
            std::swap(M.len[3], M.len[2]);
            std::swap(M.poly[2], *pt);
            M.len[2] = lent;
        } else if (len_Q == 0) {
            // Q = 0 → M[2]_new = M[3], M[3]_new = M[2]
            std::swap(M.poly[3], M.poly[2]);
            std::swap(M.len[3], M.len[2]);
        } else {
            // M[2] = 0 → M[2]_new = M[3], M[3]_new = 0
            std::swap(M.poly[3], M.poly[2]);
            std::swap(M.len[3], M.len[2]);
            // M[2] 现在是旧 M[3]，M[3] 现在是 0（旧 M[2]）
        }

        // 行 1：M[0]_new = M[1]_old + Q * M[0]_old, M[1]_new = M[0]_old
        if (len_Q > 0 && M.len[0] > 0) {
            size_t len_prod;
            if (len_Q >= M.len[0])
                _mul(*pT, Q, len_Q, M.poly[0], M.len[0], scratch);
            else
                _mul(*pT, M.poly[0], M.len[0], Q, len_Q, scratch);
            len_prod = len_Q + M.len[0] - 1;

            size_t lent = _poly_add(*pt, M.poly[1], M.len[1], *pT, len_prod);
            std::swap(M.poly[1], M.poly[0]);
            std::swap(M.len[1], M.len[0]);
            std::swap(M.poly[0], *pt);
            M.len[0] = lent;
        } else if (len_Q == 0) {
            std::swap(M.poly[1], M.poly[0]);
            std::swap(M.len[1], M.len[0]);
        } else {
            std::swap(M.poly[1], M.poly[0]);
            std::swap(M.len[1], M.len[0]);
        }

        sgn = -sgn;
    }

    return sgn;
}
```

**简化说明**：上面 `len_Q == 0` 和 `M.len[i] == 0` 的分支处理是为了正确性——当商为零或矩阵项为零时，`_mul` 的前置条件 `len_b > 0` 不满足。但在实际 HGCD 中，`len_Q >= 1`（divrem 的商至少是常数），且 M 的对角线项从 identity 开始永远非零。off-diagonal 项 M[1], M[3] 初始为零，第一步后 M[1] = 1、M[3] = 1，之后永远非零。所以只有第一步可能触发 `M.len[2] == 0` 或 `M.len[1] == 0` 的分支。

**进一步简化**：实际可以检测到 FLINT 的做法是把 `__mul` 宏内部就处理了 `lenA == 0 || lenB == 0` 的情况（输出 lenC = 0）。我们也这样做——让 `_mul` 内部处理零长度输入即可，不需要在 _hgcd_iter 里分支：

```cpp
// _mul 修改：支持 len_b == 0（输出 C 不写入，返回 length 0）
// 调用方需检查返回值
```

但这会让 `_mul` 签名变复杂。更简洁的做法是用一个 `_mat_update` 辅助函数封装行更新逻辑。

**最终方案**：

```cpp
// 矩阵行更新辅助：M[i0] = M[i1] + Q * M[i0], M[i1] = old M[i0]
// 使用 T, t 作为临时缓冲区，通过 swap 避免拷贝
void _mat_row_update(hgcd_mat& M, size_t i0, size_t i1,
                     const uint64_t* Q, size_t len_Q,
                     uint64_t** pT, size_t& len_T,
                     uint64_t** pt,
                     uint64_t* scratch) const
{
    if (len_Q != 0 && M.len[i0] != 0) {
        if (len_Q >= M.len[i0])
            _mul(*pT, Q, len_Q, M.poly[i0], M.len[i0], scratch);
        else
            _mul(*pT, M.poly[i0], M.len[i0], Q, len_Q, scratch);
        len_T = len_Q + M.len[i0] - 1;

        size_t lent = _poly_add(*pt, M.poly[i1], M.len[i1], *pT, len_T);
        std::swap(M.poly[i1], M.poly[i0]);
        std::swap(M.len[i1], M.len[i0]);
        std::swap(M.poly[i0], *pt);
        M.len[i0] = lent;
    } else {
        // Q == 0 或 M[i0] == 0：M_new[i0] = M_old[i1], M_new[i1] = M_old[i0]
        std::swap(M.poly[i1], M.poly[i0]);
        std::swap(M.len[i1], M.len[i0]);
    }
}
```

然后 `_hgcd_iter` 的矩阵更新部分简化为：

```cpp
        size_t len_T_dummy;
        _mat_row_update(M, 2, 3, Q, len_Q, pT, len_T_dummy, pt, scratch);
        _mat_row_update(M, 0, 1, Q, len_Q, pT, len_T_dummy, pt, scratch);
```

注意：两行更新都用 `*pT` 和 `*pt` 作为临时空间——第一行更新结束后 `*pT` 和 `*pt` 的内容可能已变（swap 了指针），但第二行更新的 `_mul` 输出和 `_poly_add` 输出会覆盖它们，所以没有问题。

---

## 5. M3: 矩阵乘法

M3 在 M2 之前实现，因为 M2 依赖 M3。

### 5.1 `_mat_one`

```cpp
static void _mat_one(hgcd_mat& M)
{
    M.poly[0][0] = 1; M.len[0] = 1;    // M00 = 1
    M.len[1] = 0;                        // M01 = 0
    M.len[2] = 0;                        // M10 = 0
    M.poly[3][0] = 1; M.len[3] = 1;    // M11 = 1
}
```

### 5.2 `_mat_mul`（classical）

```cpp
// C = A * B（classical 2x2 矩阵乘法）
// C[0] = A[0]*B[0] + A[1]*B[2]
// C[1] = A[0]*B[1] + A[1]*B[3]
// C[2] = A[2]*B[0] + A[3]*B[2]
// C[3] = A[2]*B[1] + A[3]*B[3]
// T: 临时多项式空间（足够容纳最长乘积）
// 不支持 aliasing（C != A, C != B）
void _mat_mul(hgcd_mat& C, const hgcd_mat& A, const hgcd_mat& B,
              uint64_t* T, uint64_t* scratch) const
{
    size_t lenT;

    // C[0] = A[0]*B[0] + A[1]*B[2]
    _mat_mul_entry(C.poly[0], C.len[0],
                   A.poly[0], A.len[0], B.poly[0], B.len[0],
                   A.poly[1], A.len[1], B.poly[2], B.len[2],
                   T, scratch);

    // C[1] = A[0]*B[1] + A[1]*B[3]
    _mat_mul_entry(C.poly[1], C.len[1],
                   A.poly[0], A.len[0], B.poly[1], B.len[1],
                   A.poly[1], A.len[1], B.poly[3], B.len[3],
                   T, scratch);

    // C[2] = A[2]*B[0] + A[3]*B[2]
    _mat_mul_entry(C.poly[2], C.len[2],
                   A.poly[2], A.len[2], B.poly[0], B.len[0],
                   A.poly[3], A.len[3], B.poly[2], B.len[2],
                   T, scratch);

    // C[3] = A[2]*B[1] + A[3]*B[3]
    _mat_mul_entry(C.poly[3], C.len[3],
                   A.poly[2], A.len[2], B.poly[1], B.len[1],
                   A.poly[3], A.len[3], B.poly[3], B.len[3],
                   T, scratch);
}

// 辅助：C_entry = P*Q + R*S
void _mat_mul_entry(uint64_t* C, size_t& lenC,
                    const uint64_t* P, size_t lenP,
                    const uint64_t* Q, size_t lenQ,
                    const uint64_t* R, size_t lenR,
                    const uint64_t* S, size_t lenS,
                    uint64_t* T, uint64_t* scratch) const
{
    // PQ
    size_t lenPQ = 0;
    if (lenP > 0 && lenQ > 0) {
        if (lenP >= lenQ)
            _mul(C, P, lenP, Q, lenQ, scratch);
        else
            _mul(C, Q, lenQ, P, lenP, scratch);
        lenPQ = lenP + lenQ - 1;
    }

    // RS
    size_t lenRS = 0;
    if (lenR > 0 && lenS > 0) {
        if (lenR >= lenS)
            _mul(T, R, lenR, S, lenS, scratch);
        else
            _mul(T, S, lenS, R, lenR, scratch);
        lenRS = lenR + lenS - 1;
    }

    // C = PQ + RS
    if (lenPQ > 0 && lenRS > 0)
        lenC = _poly_add(C, C, lenPQ, T, lenRS);
    else if (lenPQ > 0)
        lenC = lenPQ;
    else if (lenRS > 0) {
        std::memcpy(C, T, lenRS * sizeof(uint64_t));
        lenC = lenRS;
    } else
        lenC = 0;
}
```

---

## 6. M2: HGCD 递归

### 6.1 函数签名

```cpp
// 递归 HGCD
// 输入：a[len_a], b[len_b]，前置 len_a > len_b > 0
// 输出：A[len_A], B[len_B]（拷贝到调用方提供的缓冲区）
//       M（矩阵，若 compute_M 为 true）
// W: 工作区（每层消耗 9*len_a + 10*(len_a+1)/2 元素，含 word3 累加器）
// scratch: _mul 所需工作区
// 支持 aliasing：{B, a} 可重叠（见 architecture.md Decision 7）
// 返回 sgn
int _hgcd_recursive(hgcd_mat& M, bool compute_M,
                    uint64_t* A, size_t& len_A,
                    uint64_t* B, size_t& len_B,
                    const uint64_t* a, size_t len_a,
                    const uint64_t* b, size_t len_b,
                    uint64_t* W, uint64_t* scratch) const
```

### 6.2 完整实现

```cpp
int _hgcd_recursive(hgcd_mat& M, bool compute_M,
                    uint64_t* A, size_t& len_A,
                    uint64_t* B, size_t& len_B,
                    const uint64_t* a, size_t len_a,
                    const uint64_t* b, size_t len_b,
                    uint64_t* W, uint64_t* scratch) const  // W 含 word3 空间
{
    const size_t m = len_a / 2;

    // --- Base case ---
    if (len_b < m + 1) {
        if (compute_M)
            _mat_one(M);
        // A ← a 必须在 B ← b 之前（aliasing 安全）
        std::memcpy(A, a, len_a * sizeof(uint64_t));
        len_A = len_a;
        std::memcpy(B, b, len_b * sizeof(uint64_t));
        len_B = len_b;
        return 1;   // sgn = +1 (identity)
    }

    // --- 工作区切片 ---
    // 布局：a2[n] b2[n] a3[n] b3[n] q[(n+1)/2] d[n] T0[n] T1[(n+1)/2]
    //        R[0..3] 各 (n+1)/2    S[0..3] 各 (n+1)/2
    //        W3[3n]（word3 累加器，复用为 n 个 word3 = 3n 个 uint64_t）
    size_t n = len_a;
    size_t half = (n + 1) / 2;
    uint64_t* a2 = W;
    uint64_t* b2 = a2 + n;
    uint64_t* a3 = b2 + n;
    uint64_t* b3 = a3 + n;
    uint64_t* q  = b3 + n;
    uint64_t* d  = q + half;
    uint64_t* T0 = d + n;
    uint64_t* T1 = T0 + n;

    hgcd_mat R;
    R.poly[0] = T1 + half;
    R.poly[1] = R.poly[0] + half;
    R.poly[2] = R.poly[1] + half;
    R.poly[3] = R.poly[2] + half;

    hgcd_mat S;
    S.poly[0] = R.poly[3] + half;
    S.poly[1] = S.poly[0] + half;
    S.poly[2] = S.poly[1] + half;
    S.poly[3] = S.poly[2] + half;

    // word3 累加器空间：n 个 word3 = 3n 个 uint64_t
    // 同层的多次 _poly_divrem 调用共享此空间（不并发）
    word3* W3 = reinterpret_cast<word3*>(S.poly[3] + half);

    uint64_t* W_next = W + 6 * n + 10 * half + 3 * n;

    // --- 第一次递归：对高位 a0, b0 ---
    const uint64_t* a0 = a + m;
    size_t len_a0 = len_a - m;
    const uint64_t* b0 = b + m;
    size_t len_b0 = (len_b >= m) ? len_b - m : 0;

    size_t len_a3, len_b3;
    int sgnR;

    if (len_a0 < HGCD_CUTOFF) {
        // 迭代 base case
        uint64_t* pA_buf = a3;           // 复用 a3 空间作为 A 缓冲区
        uint64_t* pB_buf = b3;           // 复用 b3 空间作为 B 缓冲区
        uint64_t* pT_buf = T0;           // T 缓冲区
        uint64_t* pt_buf = T1;           // t 缓冲区

        uint64_t* pA = pA_buf;
        uint64_t* pB = pB_buf;
        uint64_t* pT = pT_buf;
        uint64_t* pt = pt_buf;

        sgnR = _hgcd_iter(R,
                          &pA, len_a3,
                          &pB, len_b3,
                          a0, len_a0, b0, len_b0,
                          q, W3, &pT, &pt, scratch);

        // 迭代后 pA, pB 可能指向 a3, b3, T0 中的任一个（指针旋转）
        // 需要确保数据回到固定位置 a3, b3
        // 注意交叉情况：pA=T0, pB=a3 时必须先保存 B 再覆盖 a3
        if (pA != a3 && pB == a3) {
            // 交叉：先拷 B（从 a3），再覆盖 a3
            std::memcpy(b3, pB, len_b3 * sizeof(uint64_t));
            std::memcpy(a3, pA, len_a3 * sizeof(uint64_t));
        } else {
            if (pA != a3)
                std::memcpy(a3, pA, len_a3 * sizeof(uint64_t));
            if (pB != b3)
                std::memcpy(b3, pB, len_b3 * sizeof(uint64_t));
        }
    } else {
        // 递归
        sgnR = _hgcd_recursive(R, true,
                                a3, len_a3,
                                b3, len_b3,
                                a0, len_a0, b0, len_b0,
                                W_next, scratch);
    }

    // --- 重构 a2, b2 ---
    // a_lo = a mod x^m (truncate), b_lo = b mod x^m
    const uint64_t* a_lo = a;
    size_t len_a_lo = std::min(len_a, m);
    const uint64_t* b_lo = b;
    size_t len_b_lo = std::min(len_b, m);

    // b2_lo：
    //   sgnR > 0: b2_lo = R[0]*b_lo - R[2]*a_lo
    //   sgnR < 0: b2_lo = R[2]*a_lo - R[0]*b_lo
    size_t lenb2;
    {
        // b2 = R[2] * a_lo
        size_t len_t1 = 0;
        if (R.len[2] > 0 && len_a_lo > 0) {
            if (R.len[2] >= len_a_lo)
                _mul(b2, R.poly[2], R.len[2], a_lo, len_a_lo, scratch);
            else
                _mul(b2, a_lo, len_a_lo, R.poly[2], R.len[2], scratch);
            len_t1 = R.len[2] + len_a_lo - 1;
        }

        // T0 = R[0] * b_lo
        size_t len_t2 = 0;
        if (R.len[0] > 0 && len_b_lo > 0) {
            if (R.len[0] >= len_b_lo)
                _mul(T0, R.poly[0], R.len[0], b_lo, len_b_lo, scratch);
            else
                _mul(T0, b_lo, len_b_lo, R.poly[0], R.len[0], scratch);
            len_t2 = R.len[0] + len_b_lo - 1;
        }

        if (sgnR < 0)
            lenb2 = _poly_sub(b2, b2, len_t1, T0, len_t2);     // R2*a_lo - R0*b_lo
        else
            lenb2 = _poly_sub(b2, T0, len_t2, b2, len_t1);     // R0*b_lo - R2*a_lo
    }

    // 零填充 b2[lenb2 .. m+lenb3-1] = 0，然后 b2[m:] += b3
    {
        if (lenb2 < m + len_b3)
            std::memset(b2 + lenb2, 0, (m + len_b3 - lenb2) * sizeof(uint64_t));

        // b2[m:] += b3
        uint64_t* b2_hi = b2 + m;
        size_t len_b2_hi = (lenb2 >= m) ? lenb2 - m : 0;
        size_t new_len_b2_hi = _poly_add(b2_hi, b2_hi, len_b2_hi, b3, len_b3);
        lenb2 = std::max(m + len_b3, lenb2);
        // normalize（可能高位加法后有 carry 或 cancellation）
        lenb2 = _poly_normalise(b2, lenb2);
    }

    // a2_lo：
    //   sgnR > 0: a2_lo = R[3]*a_lo - R[1]*b_lo
    //   sgnR < 0: a2_lo = R[1]*b_lo - R[3]*a_lo
    size_t lena2;
    {
        size_t len_t1 = 0;
        if (R.len[3] > 0 && len_a_lo > 0) {
            if (R.len[3] >= len_a_lo)
                _mul(a2, R.poly[3], R.len[3], a_lo, len_a_lo, scratch);
            else
                _mul(a2, a_lo, len_a_lo, R.poly[3], R.len[3], scratch);
            len_t1 = R.len[3] + len_a_lo - 1;
        }

        size_t len_t2 = 0;
        if (R.len[1] > 0 && len_b_lo > 0) {
            if (R.len[1] >= len_b_lo)
                _mul(T0, R.poly[1], R.len[1], b_lo, len_b_lo, scratch);
            else
                _mul(T0, b_lo, len_b_lo, R.poly[1], R.len[1], scratch);
            len_t2 = R.len[1] + len_b_lo - 1;
        }

        if (sgnR < 0)
            lena2 = _poly_sub(a2, T0, len_t2, a2, len_t1);     // R1*b_lo - R3*a_lo
        else
            lena2 = _poly_sub(a2, a2, len_t1, T0, len_t2);     // R3*a_lo - R1*b_lo
    }

    // 零填充 + a2[m:] += a3
    {
        if (lena2 < m + len_a3)
            std::memset(a2 + lena2, 0, (m + len_a3 - lena2) * sizeof(uint64_t));

        uint64_t* a2_hi = a2 + m;
        size_t len_a2_hi = (lena2 >= m) ? lena2 - m : 0;
        _poly_add(a2_hi, a2_hi, len_a2_hi, a3, len_a3);
        lena2 = std::max(m + len_a3, lena2);
        lena2 = _poly_normalise(a2, lena2);
    }

    // --- 提前终止 ---
    if (lenb2 < m + 1) {
        std::memcpy(A, a2, lena2 * sizeof(uint64_t));
        len_A = lena2;
        std::memcpy(B, b2, lenb2 * sizeof(uint64_t));
        len_B = lenb2;

        if (compute_M) {
            // M = R
            for (int i = 0; i < 4; ++i) {
                std::memcpy(M.poly[i], R.poly[i], R.len[i] * sizeof(uint64_t));
                M.len[i] = R.len[i];
            }
        }
        return sgnR;
    }

    // --- 中间 divrem ---
    size_t lenq, lend;
    std::tie(lenq, lend) = _poly_divrem(q, d, a2, lena2, b2, lenb2, W3);

    // --- 第二次递归：对高位 c0, d0 ---
    size_t k = 2 * m - lenb2 + 1;
    const uint64_t* c0 = b2 + k;
    size_t lenc0 = (lenb2 >= k) ? lenb2 - k : 0;
    const uint64_t* d0 = d + k;
    size_t lend0 = (lend >= k) ? lend - k : 0;

    size_t len_a3b, len_b3b;
    int sgnS;

    if (lenc0 < HGCD_CUTOFF) {
        uint64_t* pA_buf = a3;
        uint64_t* pB_buf = b3;
        uint64_t* pT_buf = T0;
        uint64_t* pt_buf = T1;

        uint64_t* pA = pA_buf;
        uint64_t* pB = pB_buf;
        uint64_t* pT = pT_buf;
        uint64_t* pt = pt_buf;

        sgnS = _hgcd_iter(S,
                          &pA, len_a3b,
                          &pB, len_b3b,
                          c0, lenc0, d0, lend0,
                          a2, W3, &pT, &pt, scratch);
        // 注：FLINT 在此处复用 a2 作为 Q 缓冲区——我们用 a2 因为此时 a2 已不需要

        // 交叉保护（同第一次迭代）
        if (pA != a3 && pB == a3) {
            std::memcpy(b3, pB, len_b3b * sizeof(uint64_t));
            std::memcpy(a3, pA, len_a3b * sizeof(uint64_t));
        } else {
            if (pA != a3)
                std::memcpy(a3, pA, len_a3b * sizeof(uint64_t));
            if (pB != b3)
                std::memcpy(b3, pB, len_b3b * sizeof(uint64_t));
        }
    } else {
        sgnS = _hgcd_recursive(S, true,
                                a3, len_a3b,
                                b3, len_b3b,
                                c0, lenc0, d0, lend0,
                                W_next, scratch);
    }

    // --- 第二次重构：输出到 A, B ---
    // b2_lo_trunc = b2 mod x^k, d_lo_trunc = d mod x^k
    const uint64_t* s2 = b2;
    size_t lens2 = std::min(lenb2, k);
    const uint64_t* t2 = d;
    size_t lent2 = std::min(lend, k);

    // B：
    //   sgnS > 0: B_lo = S[0]*t2 - S[2]*s2
    //   sgnS < 0: B_lo = S[2]*s2 - S[0]*t2
    {
        size_t len_t1 = 0;
        if (S.len[2] > 0 && lens2 > 0) {
            if (S.len[2] >= lens2)
                _mul(B, S.poly[2], S.len[2], s2, lens2, scratch);
            else
                _mul(B, s2, lens2, S.poly[2], S.len[2], scratch);
            len_t1 = S.len[2] + lens2 - 1;
        }

        size_t len_t2 = 0;
        if (S.len[0] > 0 && lent2 > 0) {
            if (S.len[0] >= lent2)
                _mul(T0, S.poly[0], S.len[0], t2, lent2, scratch);
            else
                _mul(T0, t2, lent2, S.poly[0], S.len[0], scratch);
            len_t2 = S.len[0] + lent2 - 1;
        }

        if (sgnS < 0)
            len_B = _poly_sub(B, B, len_t1, T0, len_t2);
        else
            len_B = _poly_sub(B, T0, len_t2, B, len_t1);

        if (len_B < k + len_b3b)
            std::memset(B + len_B, 0, (k + len_b3b - len_B) * sizeof(uint64_t));

        uint64_t* B_hi = B + k;
        size_t len_B_hi = (len_B >= k) ? len_B - k : 0;
        _poly_add(B_hi, B_hi, len_B_hi, b3, len_b3b);
        len_B = std::max(k + len_b3b, len_B);
        len_B = _poly_normalise(B, len_B);
    }

    // A：
    //   sgnS > 0: A_lo = S[3]*s2 - S[1]*t2
    //   sgnS < 0: A_lo = S[1]*t2 - S[3]*s2
    {
        size_t len_t1 = 0;
        if (S.len[3] > 0 && lens2 > 0) {
            if (S.len[3] >= lens2)
                _mul(A, S.poly[3], S.len[3], s2, lens2, scratch);
            else
                _mul(A, s2, lens2, S.poly[3], S.len[3], scratch);
            len_t1 = S.len[3] + lens2 - 1;
        }

        size_t len_t2 = 0;
        if (S.len[1] > 0 && lent2 > 0) {
            if (S.len[1] >= lent2)
                _mul(T0, S.poly[1], S.len[1], t2, lent2, scratch);
            else
                _mul(T0, t2, lent2, S.poly[1], S.len[1], scratch);
            len_t2 = S.len[1] + lent2 - 1;
        }

        if (sgnS < 0)
            len_A = _poly_sub(A, T0, len_t2, A, len_t1);
        else
            len_A = _poly_sub(A, A, len_t1, T0, len_t2);

        if (len_A < k + len_a3b)
            std::memset(A + len_A, 0, (k + len_a3b - len_A) * sizeof(uint64_t));

        uint64_t* A_hi = A + k;
        size_t len_A_hi = (len_A >= k) ? len_A - k : 0;
        _poly_add(A_hi, A_hi, len_A_hi, a3, len_a3b);
        len_A = std::max(k + len_a3b, len_A);
        len_A = _poly_normalise(A, len_A);
    }

    // --- 矩阵合并：M = R * ([[q,1],[1,0]] * S) ---
    if (compute_M) {
        // S_modified = [[q,1],[1,0]] * S
        // = [[S[2] + q*S[0],  S[3] + q*S[1]],
        //    [S[0],           S[1]          ]]
        // 实现：swap(S[0],S[2]); swap(S[1],S[3]);
        //       S[0] += q * S[2];  S[1] += q * S[3]
        std::swap(S.poly[0], S.poly[2]);
        std::swap(S.len[0], S.len[2]);
        std::swap(S.poly[1], S.poly[3]);
        std::swap(S.len[1], S.len[3]);

        if (lenq > 0 && S.len[2] > 0) {
            if (lenq >= S.len[2])
                _mul(T0, q, lenq, S.poly[2], S.len[2], scratch);
            else
                _mul(T0, S.poly[2], S.len[2], q, lenq, scratch);
            size_t lenT = lenq + S.len[2] - 1;
            S.len[0] = _poly_add(S.poly[0], S.poly[0], S.len[0], T0, lenT);
        }

        if (lenq > 0 && S.len[3] > 0) {
            if (lenq >= S.len[3])
                _mul(T0, q, lenq, S.poly[3], S.len[3], scratch);
            else
                _mul(T0, S.poly[3], S.len[3], q, lenq, scratch);
            size_t lenT = lenq + S.len[3] - 1;
            S.len[1] = _poly_add(S.poly[1], S.poly[1], S.len[1], T0, lenT);
        }

        // M = R * S_modified
        // 使用 a2, b2 作为矩阵乘法的临时空间（此时 a2, b2 已不需要）
        _mat_mul(M, R, S, a2, scratch);
    }

    return -(sgnR * sgnS);
}
```

### 6.3 对照 FLINT 验证清单

| FLINT 行号 | 操作 | 对应位置 |
|------------|------|---------|
| 390-401 | base case: mat_one + set(A,a) + set(B,b) | base case 段 |
| 416-434 | 工作区切片 | 工作区切片段 |
| 436-437 | a0 = a >> m, b0 = b >> m | `a0 = a + m` |
| 447-455 | 第一次递归/迭代 | 第一次递归段 |
| 464-465 | s = truncate(a, m), t = truncate(b, m) | `a_lo, b_lo` |
| 467-481 | 重构 b2 | b2 重构段 |
| 483-496 | 重构 a2 | a2 重构段 |
| 498-512 | 提前终止 | 提前终止段 |
| 516 | k = 2*m - lenb2 + 1 | 同 |
| 552-554 | divrem + attach_shift | 中间 divrem 段 |
| 588-595 | 第二次递归/迭代 | 第二次递归段 |
| 604-633 | 第二次重构 → A, B | 第二次重构段 |
| 637-644 | S_modified + mat_mul | 矩阵合并段 |
| 647 | sgn = -(sgnR * sgnS) | 返回值 |

---

## 7. M4: GCD 分派与 HGCD-based GCD

### 7.1 `_gcd_euclid`（raw API）

从现有 Euclid 循环提取的 raw 版本，供 `_gcd_hgcd` fallback 调用。

```cpp
// Euclid GCD（raw API）
// G 预分配 len_a 元素
// 前置：len_a >= len_b > 0
void _gcd_euclid(uint64_t* G, size_t& len_G,
                 const uint64_t* A, size_t len_a,
                 const uint64_t* B, size_t len_b) const
{
    // 需要 3 个临时缓冲区 + Q + word3 累加器
    std::vector<uint64_t> a_buf(A, A + len_a);
    std::vector<uint64_t> b_buf(B, B + len_b);
    std::vector<uint64_t> q_buf(len_a);
    std::vector<uint64_t> r_buf(len_a);
    std::vector<word3> w3(len_a);

    uint64_t* a = a_buf.data();
    size_t la = len_a;
    uint64_t* b = b_buf.data();
    size_t lb = len_b;

    while (lb > 0) {
        auto [lq, lr] = _poly_divrem(q_buf.data(), r_buf.data(), a, la, b, lb, w3.data());
        // a ← b, b ← r
        std::swap(a_buf, b_buf);
        a = a_buf.data();
        la = lb;
        b = b_buf.data();

        // r_buf 的数据需要到 b
        std::memcpy(b_buf.data(), r_buf.data(), lr * sizeof(uint64_t));
        lb = lr;
    }

    std::memcpy(G, a, la * sizeof(uint64_t));
    len_G = la;
}
```

注：此实现内部分配临时缓冲区。因为 `_gcd_euclid` 只在 `len < GCD_CUTOFF = 340` 时调用，分配开销可忽略。后续可优化为使用调用方提供的工作区。

### 7.2 `_gcd_hgcd`

```cpp
// HGCD-based GCD
// 前置：len_a >= len_b >= GCD_CUTOFF
void _gcd_hgcd(uint64_t* G, size_t& len_G,
               const uint64_t* A, size_t len_a,
               const uint64_t* B, size_t len_b) const
{
    size_t n = len_a;

    // 工作区分配
    // GCD 外层：J[n] + R[n] + Q[n] + W3_gcd[3n]
    // HGCD 递归：28*n + 16*(ceil_log2(n)+1)（含每层 word3 空间）
    // scratch: 7*n（_mul 所需）
    size_t log2n = 0;
    { size_t tmp = n; while (tmp > 1) { tmp = (tmp + 1) / 2; ++log2n; } }

    size_t hgcd_ws = 28 * n + 16 * (log2n + 1);
    size_t scratch_size = 7 * n + 1;
    size_t half = (n + 1) / 2;
    size_t mat_size = 4 * half;  // M_dummy 的 4 个行向量，compute_M=false 不写入
    size_t total = 6 * n + hgcd_ws + scratch_size + mat_size;

    std::vector<uint64_t> workspace(total, 0);
    uint64_t* J = workspace.data();
    uint64_t* R_buf = J + n;
    uint64_t* Q_buf = R_buf + n;
    word3* W3_gcd = reinterpret_cast<word3*>(Q_buf + n);  // 3n 个 uint64_t
    uint64_t* W = Q_buf + n + 3 * n;
    uint64_t* scratch = W + hgcd_ws;
    uint64_t* mat_buf = scratch + scratch_size;

    // M_dummy：compute_M = false 时矩阵空间不写入，但参数要有效
    hgcd_mat M_dummy;
    for (int i = 0; i < 4; ++i) {
        M_dummy.poly[i] = mat_buf + i * half;
        M_dummy.len[i] = 0;
    }

    // --- 初始 divrem ---
    size_t len_Q, len_R;
    std::tie(len_Q, len_R) = _poly_divrem(Q_buf, R_buf, A, len_a, B, len_b, W3_gcd);

    if (len_R == 0) {
        // B 整除 A
        std::memcpy(G, B, len_b * sizeof(uint64_t));
        len_G = len_b;
        return;
    }

    // --- 首次 HGCD：缩减 (B, R) → (G, J) ---
    size_t len_g, len_j;
    _hgcd_recursive(M_dummy, false,
                    G, len_g,
                    J, len_j,
                    B, len_b,
                    R_buf, len_R,
                    W, scratch);

    // --- 主循环 ---
    while (len_j != 0) {
        // divrem: R = G mod J
        // _poly_divrem 已处理 len_g < len_j 的情况（Q=0, R=G）
        std::tie(len_Q, len_R) = _poly_divrem(Q_buf, R_buf, G, len_g, J, len_j, W3_gcd);

        if (len_R == 0) {
            // J 整除 G
            std::memcpy(G, J, len_j * sizeof(uint64_t));
            len_G = len_j;
            return;
        }

        if (len_j < GCD_CUTOFF) {
            // fallback Euclid
            _gcd_euclid(G, len_G, J, len_j, R_buf, len_R);
            return;
        }

        // HGCD：缩减 (J, R) → (G, J)
        // 注意 aliasing：输出 J 与输入 a (= J) 重叠（Decision 7）
        _hgcd_recursive(M_dummy, false,
                        G, len_g,
                        J, len_j,
                        J, len_j,
                        R_buf, len_R,
                        W, scratch);
    }

    // J == 0：G 即为 GCD
    len_G = len_g;
}
```

### 7.3 修改 `gcd` 分派

```cpp
static void gcd(dense_upoly_zp& G,
                const dense_upoly_zp& A, const dense_upoly_zp& B)
{
    dense_upoly_zp a = A, b = B;
    if (a.deg() < b.deg()) std::swap(a, b);

    if (b.empty()) {
        G = std::move(a);
        return;
    }

    size_t len_b = b._coeffs.size();

    if (len_b < GCD_CUTOFF) {
        // 现有 Euclid
        dense_upoly_zp q, r;
        while (!b.empty()) {
            divrem(q, r, a, b);
            a = std::move(b);
            b = std::move(r);
        }
        G = std::move(a);
    } else {
        // HGCD-based GCD
        G._p = a._p; G._ninv = a._ninv; G._norm = a._norm;
        G._coeffs.resize(a._coeffs.size());
        size_t len_g;
        G._gcd_hgcd(G._coeffs.data(), len_g,
                     a._coeffs.data(), a._coeffs.size(),
                     b._coeffs.data(), b._coeffs.size());
        G._coeffs.resize(len_g);
    }
}
```

---

## 8. 工作区总预算

每层递归消耗（含 word3 累加器）：

```
a2[n] + b2[n] + a3[n] + b3[n] + q[half] + d[n] + T0[n] + T1[half]
+ R[0..3]*half + S[0..3]*half + W3[3n]
= 9n + 10*half
```

总递归工作区：

```
W_hgcd = sum_{i=0}^{log n} (9*(n/2^i) + 10*(n/2^i + 1)/2)
       ≈ 28n + 16*(ceil_log2(n) + 1)
```

| 组件 | 大小（uint64_t 个数） | 说明 |
|------|------|------|
| GCD 外层：J, R, Q | 3n | 临时多项式 |
| GCD 外层：W3_gcd | 3n | GCD 外层 divrem 的 word3 空间 |
| HGCD 递归工作区 W | 28n + 16*log2(n) | 逐层切片（含 word3） |
| _mul scratch | 7n | Karatsuba 零填充 + scratch |
| 矩阵 M_dummy | 4 * (n+1)/2 | GCD 不需要矩阵但参数要有效 |
| **总计** | ~43n | 对 n=10000 约 3.4MB |

**零动态分配**：整个 HGCD GCD 过程只在 `_gcd_hgcd` 入口做一次 `std::vector<uint64_t>` 分配，内部所有操作（divrem、乘法、矩阵运算）均使用工作区指针切片。

---

## 9. _hgcd_recursive 中 _hgcd_iter 的缓冲区复用

迭代 base case 需要 5 个缓冲区（A, B, T 各 `len_a0` 大小，再加 `M.poly[0..3]` 和 `t` 各 `half` 大小）。这些全部从工作区 W 中已分配的空间复用：

| 迭代缓冲区 | 复用的 W 空间 | 大小 |
|------------|-------------|------|
| pA (多项式 A) | a3 | n |
| pB (多项式 B) | b3 | n |
| pT (多项式 T) | T0 | n |
| pt (矩阵临时) | T1 | half |
| Q (商) | q | half |
| R.poly[0..3] | R.poly[0..3] | 各 half |
| scratch | 外部传入 | 7n |

迭代结束后 `pA` 和 `pB` 可能已被 swap 到 `a3`/`b3`/`T0` 中的任意位置。如果 `pA != a3`，需要拷贝回 `a3`（重构步需要 `a3` 在固定位置）。

---

## 10. 实施步骤

### Step 1: M0 基础操作
新增到 `dense_upoly_zp` private 区域：
- `_poly_normalise`（static）
- `_poly_add`、`_poly_sub`
- `_poly_divrem`
- `_mul`

新增到 public 区域：
- `add`、`sub`

### Step 2: M3 矩阵操作
新增到 private 区域：
- `hgcd_mat` 结构
- `_mat_one`（static）
- `_mat_mul_entry`
- `_mat_mul`

### Step 3: M1 迭代 base case
新增到 private 区域：
- `_mat_row_update`
- `_hgcd_iter`

### Step 4: M2 HGCD 递归
新增到 private 区域：
- `_hgcd_recursive`

### Step 5: M4 GCD 分派
新增到 private 区域：
- `_gcd_euclid`
- `_gcd_hgcd`

修改 public `gcd`：加分派逻辑。

### Step 6: 测试
1. M0 单元测试：add/sub/divrem/_mul 正确性
2. M1-M3 集成测试：HGCD 不变量 `M * (A,B) = (a,b)` + det(M) = sgn
3. M4 全量回归：`bash test/run_all_tests.sh` + `make crosscheck`
4. 性能：`make bench-clpoly` + `make bench-comparative`

---

## 11. 代码行数估算

| 模块 | 新增行数 |
|------|---------|
| 常量 + hgcd_mat | ~10 |
| M0: add/sub/divrem/_mul + 对象级 | ~130 |
| M3: mat_one/mat_mul_entry/mat_mul | ~60 |
| M1: mat_row_update + hgcd_iter | ~70 |
| M2: hgcd_recursive | ~200 |
| M4: gcd_euclid + gcd_hgcd + gcd 修改 | ~80 |
| **总计** | ~550 行 |
