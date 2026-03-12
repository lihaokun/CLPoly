/**
 * @file dense_upoly_zp.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 稠密单变量 Zp 多项式类型，用于性能关键路径
 *
 * 设计动机：upolynomial_<Zp> 使用稀疏表示 (40 字节/项)，
 * 而 Zp GCD 中的中间多项式实际稠密。此类型使用连续 uint64_t
 * 数组 (8 字节/系数)，内存减少 5x，cache 友好。
 *
 * 参考实现：FLINT nmod_poly（相同的稠密 + Barrett 设计）
 */
#ifndef CLPOLY_DENSE_UPOLY_ZP_HH
#define CLPOLY_DENSE_UPOLY_ZP_HH

#include <vector>
#include <cstdint>
#include <cassert>
#include <algorithm>
#include <cstring>
#include <clpoly/upolynomial.hh>

namespace clpoly {

class dense_upoly_zp {
private:
    std::vector<uint64_t> _coeffs;  // 升幂 [c_0, ..., c_deg], 保证无 leading zero
    uint64_t _p;
    uint64_t _ninv;
    uint32_t _norm;

    // FLINT 归一化逆元（pn 必须最高位为 1）
    static uint64_t __preinvert_limb(uint64_t pn)
    {
        assert(pn >> 63);
        unsigned __int128 num = ((unsigned __int128)(~pn) << 64) | ~(uint64_t)0;
        return (uint64_t)(num / pn);
    }

    // 预计算 ninv 和 norm
    void __precompute()
    {
        if (_p == 0) { _ninv = 0; _norm = 0; return; }
        assert(_p >= 2);
        _norm = (uint32_t)__builtin_clzll(_p);
        _ninv = __preinvert_limb(_p << _norm);
    }

    // 溢出安全模加：64-bit 素数下 a+b 可能溢出，先算 p-a 再比较
    uint64_t nmod_add(uint64_t a, uint64_t b) const
    {
        uint64_t neg = _p - a;      // p - a ∈ [1, p], 不溢出
        return (neg > b) ? a + b    // a + b < p
                         : b - neg; // = a + b - p
    }

    // 溢出安全模减
    uint64_t nmod_sub(uint64_t a, uint64_t b) const
    {
        return (a >= b) ? a - b : _p - b + a;
    }

    // 去除高次零系数
    void __strip()
    {
        while (!_coeffs.empty() && _coeffs.back() == 0)
            _coeffs.pop_back();
    }

    // --- M2.0 辅助类型与函数（_classical_mul 和 divrem 共用）---

    struct word3 { uint64_t lo, mid, hi; };

    static void _umul128(uint64_t& hi, uint64_t& lo, uint64_t a, uint64_t b)
    {
        unsigned __int128 prod = (unsigned __int128)a * b;
        lo = (uint64_t)prod;
        hi = (uint64_t)(prod >> 64);
    }

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
        unsigned __int128 sum = (unsigned __int128)s.lo + b0;
        s.lo = (uint64_t)sum;
        sum = (unsigned __int128)s.mid + b1 + (uint64_t)(sum >> 64);
        s.mid = (uint64_t)sum;
        s.hi += (uint64_t)(sum >> 64);
#endif
    }

    static uint64_t _lll_mod_preinv(uint64_t hi, uint64_t mid, uint64_t lo,
                                     uint64_t p, uint64_t pinv, uint32_t norm)
    {
        // 第一步：归约 (hi, mid) → r1 ∈ [0, p)
        uint64_t r1;
        {
            uint64_t pn = p << norm;
            uint64_t h_shifted = hi << norm;
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

    // --- M1 Karatsuba 乘法 ---

    static constexpr size_t KARATSUBA_THRESHOLD = 16;

    // schoolbook 基例（惰性累加，点积形式）
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
                _umul128(p1, p0, A[j], B[k - j]);
                _add_carry3(acc, p1, p0);
            }
            C[k] = _lll_mod_preinv(acc.hi, acc.mid, acc.lo,
                                    _p, _ninv, _norm);
        }
    }

    // Karatsuba 递归核心
    void _kar_mul(uint64_t* C,
                  const uint64_t* A, const uint64_t* B,
                  size_t n, uint64_t* scratch) const
    {
        if (n < KARATSUBA_THRESHOLD) {
            _classical_mul(C, A, n, B, n);
            return;
        }

        size_t m = n / 2;
        size_t h = n - m;

        // scratch 布局：t1[h] | t2[h] | P0[2m-1] | P1[2h-1] | rec[...]
        uint64_t* t1  = scratch;
        uint64_t* t2  = t1 + h;
        uint64_t* sP0 = t2 + h;
        uint64_t* sP1 = sP0 + (2 * m - 1);
        uint64_t* rec = sP1 + (2 * h - 1);

        // t1 = A_lo + A_hi, t2 = B_lo + B_hi
        for (size_t i = 0; i < m; ++i) {
            t1[i] = nmod_add(A[i], A[m + i]);
            t2[i] = nmod_add(B[i], B[m + i]);
        }
        if (h > m) {
            t1[m] = A[m + m];
            t2[m] = B[m + m];
        }

        // P0 = A_lo × B_lo → sP0
        _kar_mul(sP0, A, B, m, rec);

        // P1 = t1 × t2 → sP1
        _kar_mul(sP1, t1, t2, h, rec);

        // P2 = A_hi × B_hi → C[2m..2n-2]
        _kar_mul(C + 2 * m, A + m, B + m, h, rec);

        // P1 -= P0
        for (size_t i = 0; i < 2 * m - 1; ++i)
            sP1[i] = nmod_sub(sP1[i], sP0[i]);

        // P1 -= P2
        for (size_t i = 0; i < 2 * h - 1; ++i)
            sP1[i] = nmod_sub(sP1[i], C[2 * m + i]);

        // 组装 C
        std::memcpy(C, sP0, (2 * m - 1) * sizeof(uint64_t));
        C[2 * m - 1] = 0;
        for (size_t i = 0; i < 2 * h - 1; ++i)
            C[m + i] = nmod_add(C[m + i], sP1[i]);
    }

    // --- P3-HGCD: 常量与辅助类型 ---

    static constexpr size_t GCD_CUTOFF  = 340;   // gcd: Euclid vs HGCD
    static constexpr size_t HGCD_CUTOFF = 100;   // HGCD: 迭代 vs 递归

    struct hgcd_mat {
        uint64_t* poly[4];   // [0]=M00, [1]=M01, [2]=M10, [3]=M11
        size_t len[4];        // 各多项式的 length（非 degree）
    };

private:
    // --- M0: 基础操作（raw pointer API）---

    // 返回去除高位零后的 length
    static size_t _poly_normalise(const uint64_t* A, size_t len)
    {
        while (len > 0 && A[len - 1] == 0)
            --len;
        return len;
    }

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

    // 稠密长除法（惰性 3-word 归约，raw API）
    // len_a >= len_b 时：Q 预分配 len_a - len_b + 1，R 预分配 len_b - 1
    // len_a < len_b 时：Q 不写入（len_q=0），R 预分配 len_a
    // W3: 预分配 len_a 个 word3 的工作区
    // 返回 {len_q, len_r}，均已 normalize
    std::pair<size_t, size_t>
    _poly_divrem(uint64_t* Q, uint64_t* R,
                 const uint64_t* A, size_t len_a,
                 const uint64_t* B, size_t len_b,
                 word3* W3) const
    {
        assert(len_b > 0);

        // FLINT 兼容：len_a < len_b 时 Q=0, R=A
        if (len_a < len_b) {
            std::memcpy(R, A, len_a * sizeof(uint64_t));
            return {0, len_a};
        }

        size_t d = len_b - 1;       // deg(B)
        size_t q_len = len_a - d;   // len_a - len_b + 1
        uint64_t inv_lc = nmod_inv(B[d]);

        // 初始化 3-word 累加器
        for (size_t i = 0; i < len_a; ++i)
            W3[i] = {A[i], 0, 0};

        for (size_t ii = q_len; ii > 0; --ii) {
            size_t i = ii - 1;
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

    // 不等长乘法（raw API）
    // C 预分配 2*len_a - 1 元素
    // scratch 预分配 7*len_a（b_pad[len_a] + kar_scratch[6*len_a]）
    // 前置：len_a >= len_b > 0
    void _mul(uint64_t* C,
              const uint64_t* A, size_t len_a,
              const uint64_t* B, size_t len_b,
              uint64_t* scratch) const
    {
        assert(len_a >= len_b && len_b > 0);

        if (len_b < KARATSUBA_THRESHOLD) {
            // schoolbook：写 len_a + len_b - 1 元素
            _classical_mul(C, A, len_a, B, len_b);
        } else {
            // 零填充到等长后 Karatsuba
            size_t n = len_a;
            uint64_t* b_pad = scratch;
            uint64_t* kar_scratch = scratch + n;
            std::memcpy(b_pad, B, len_b * sizeof(uint64_t));
            std::memset(b_pad + len_b, 0, (n - len_b) * sizeof(uint64_t));
            _kar_mul(C, A, b_pad, n, kar_scratch);
        }
    }

    // --- M3: 矩阵操作 ---

    // 矩阵初始化为单位矩阵
    static void _mat_one(hgcd_mat& M)
    {
        M.poly[0][0] = 1; M.len[0] = 1;    // M00 = 1
        M.len[1] = 0;                        // M01 = 0
        M.len[2] = 0;                        // M10 = 0
        M.poly[3][0] = 1; M.len[3] = 1;    // M11 = 1
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

    // C = A * B（2x2 矩阵乘法）
    // T: 临时多项式空间（足够容纳最长乘积）
    // 不支持 aliasing（C != A, C != B）
    void _mat_mul(hgcd_mat& C, const hgcd_mat& A, const hgcd_mat& B,
                  uint64_t* T, uint64_t* scratch) const
    {
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
            // Q == 0 或 M[i0] == 0
            std::swap(M.poly[i1], M.poly[i0]);
            std::swap(M.len[i1], M.len[i0]);
        }
    }

    // --- M1: HGCD 迭代 base case ---

    // 迭代 Euclid + 矩阵追踪
    // 输入：a[len_a], b[len_b]，前置 len_a > len_b > 0
    // 输出：M（矩阵），*pA/*pB（缩减后多项式的指针），len_A/len_B（长度）
    // 缓冲区：pA, pB, pT（各 len_a）, pt（half）, Q（half）, W3, scratch
    // 返回 sgn（行列式符号）
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
            std::swap(*pB, *pT);
            std::swap(*pA, *pT);
            len_A = len_B;
            len_B = len_T;

            // 矩阵更新
            size_t len_T_dummy;
            _mat_row_update(M, 2, 3, Q, len_Q, pT, len_T_dummy, pt, scratch);
            _mat_row_update(M, 0, 1, Q, len_Q, pT, len_T_dummy, pt, scratch);

            sgn = -sgn;
        }

        return sgn;
    }

    // --- M2: HGCD 递归 ---

    // 递归 HGCD
    // 输入：a[len_a], b[len_b]，前置 len_a > len_b > 0
    // 输出：A[len_A], B[len_B]（拷贝到调用方提供的缓冲区）
    //       M（矩阵，若 compute_M 为 true）
    // W: 工作区（每层消耗 9n + 10*half + 3n 元素）
    // scratch: _mul 所需工作区
    // 支持 aliasing：{B, a} 可重叠
    // 返回 sgn
    int _hgcd_recursive(hgcd_mat& M, bool compute_M,
                        uint64_t* A, size_t& len_A,
                        uint64_t* B, size_t& len_B,
                        const uint64_t* a, size_t len_a,
                        const uint64_t* b, size_t len_b,
                        uint64_t* W, uint64_t* scratch) const
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
            uint64_t* pA = a3;
            uint64_t* pB = b3;
            uint64_t* pT = T0;
            uint64_t* pt = T1;

            // 保存原始 R.poly 指针（_hgcd_iter 内部 _mat_row_update 会
            // 将 pt(=T1) 与 R.poly[i] 交换，导致 R.poly 指向 T0/T1。
            // 第二次 _hgcd_iter 复用 T0/T1 时会破坏 R 数据。）
            uint64_t* R_orig[4] = {R.poly[0], R.poly[1], R.poly[2], R.poly[3]};

            sgnR = _hgcd_iter(R,
                              &pA, len_a3,
                              &pB, len_b3,
                              a0, len_a0, b0, len_b0,
                              q, W3, &pT, &pt, scratch);

            // 稳定化 R：矩阵项可能被 _mat_row_update 旋转到 T0/T1，
            // 需要拷贝回原始缓冲区，避免第二次 _hgcd_iter 复用 T0/T1 时破坏数据。
            // 分两阶段：先全部暂存到 a2（与 R.poly 不重叠），再拷回原位。
            {
                uint64_t* stage = a2;  // a2 此时可安全用作暂存
                size_t off = 0;
                for (int i = 0; i < 4; ++i) {
                    std::memcpy(stage + off, R.poly[i], R.len[i] * sizeof(uint64_t));
                    off += R.len[i];
                }
                off = 0;
                for (int i = 0; i < 4; ++i) {
                    std::memcpy(R_orig[i], stage + off, R.len[i] * sizeof(uint64_t));
                    R.poly[i] = R_orig[i];
                    off += R.len[i];
                }
            }

            // 交叉保护：pA=T0, pB=a3 时必须先保存 B 再覆盖 a3
            if (pA != a3 && pB == a3) {
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

        // --- 重构 b2 ---
        const uint64_t* a_lo = a;
        size_t len_a_lo = std::min(len_a, m);
        const uint64_t* b_lo = b;
        size_t len_b_lo = std::min(len_b, m);

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
                lenb2 = _poly_sub(b2, b2, len_t1, T0, len_t2);
            else
                lenb2 = _poly_sub(b2, T0, len_t2, b2, len_t1);
        }

        // 零填充 + b2[m:] += b3
        {
            if (lenb2 < m + len_b3)
                std::memset(b2 + lenb2, 0, (m + len_b3 - lenb2) * sizeof(uint64_t));

            uint64_t* b2_hi = b2 + m;
            size_t len_b2_hi = (lenb2 >= m) ? lenb2 - m : 0;
            _poly_add(b2_hi, b2_hi, len_b2_hi, b3, len_b3);
            lenb2 = std::max(m + len_b3, lenb2);
            lenb2 = _poly_normalise(b2, lenb2);
        }

        // --- 重构 a2 ---
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
                lena2 = _poly_sub(a2, T0, len_t2, a2, len_t1);
            else
                lena2 = _poly_sub(a2, a2, len_t1, T0, len_t2);
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
            uint64_t* pA = a3;
            uint64_t* pB = b3;
            uint64_t* pT = T0;
            uint64_t* pt = T1;

            // 保存原始 S.poly 指针（同 R 的稳定化理由）
            uint64_t* S_orig[4] = {S.poly[0], S.poly[1], S.poly[2], S.poly[3]};

            sgnS = _hgcd_iter(S,
                              &pA, len_a3b,
                              &pB, len_b3b,
                              c0, lenc0, d0, lend0,
                              a2, W3, &pT, &pt, scratch);
            // 注：复用 a2 作为 Q 缓冲区（此时 a2 已不需要）

            // 稳定化 S（同 R 的理由）
            {
                uint64_t* stage = a2;
                size_t off = 0;
                for (int i = 0; i < 4; ++i) {
                    std::memcpy(stage + off, S.poly[i], S.len[i] * sizeof(uint64_t));
                    off += S.len[i];
                }
                off = 0;
                for (int i = 0; i < 4; ++i) {
                    std::memcpy(S_orig[i], stage + off, S.len[i] * sizeof(uint64_t));
                    S.poly[i] = S_orig[i];
                    off += S.len[i];
                }
            }

            // 交叉保护
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
        const uint64_t* s2 = b2;
        size_t lens2 = std::min(lenb2, k);
        const uint64_t* t2 = d;
        size_t lent2 = std::min(lend, k);

        // B
        {
            size_t len_t1 = 0;
            if (S.len[2] > 0 && lens2 > 0) {
                if (S.len[2] >= lens2)
                    _mul(B, S.poly[2], S.len[2], s2, lens2, scratch);
                else
                    _mul(B, s2, lens2, S.poly[2], S.len[2], scratch);
                len_t1 = S.len[2] + lens2 - 1;
            }

            size_t len_t2b = 0;
            if (S.len[0] > 0 && lent2 > 0) {
                if (S.len[0] >= lent2)
                    _mul(T0, S.poly[0], S.len[0], t2, lent2, scratch);
                else
                    _mul(T0, t2, lent2, S.poly[0], S.len[0], scratch);
                len_t2b = S.len[0] + lent2 - 1;
            }

            if (sgnS < 0)
                len_B = _poly_sub(B, B, len_t1, T0, len_t2b);
            else
                len_B = _poly_sub(B, T0, len_t2b, B, len_t1);

            if (len_B < k + len_b3b)
                std::memset(B + len_B, 0, (k + len_b3b - len_B) * sizeof(uint64_t));

            uint64_t* B_hi = B + k;
            size_t len_B_hi = (len_B >= k) ? len_B - k : 0;
            _poly_add(B_hi, B_hi, len_B_hi, b3, len_b3b);
            len_B = std::max(k + len_b3b, len_B);
            len_B = _poly_normalise(B, len_B);
        }

        // A
        {
            size_t len_t1 = 0;
            if (S.len[3] > 0 && lens2 > 0) {
                if (S.len[3] >= lens2)
                    _mul(A, S.poly[3], S.len[3], s2, lens2, scratch);
                else
                    _mul(A, s2, lens2, S.poly[3], S.len[3], scratch);
                len_t1 = S.len[3] + lens2 - 1;
            }

            size_t len_t2b = 0;
            if (S.len[1] > 0 && lent2 > 0) {
                if (S.len[1] >= lent2)
                    _mul(T0, S.poly[1], S.len[1], t2, lent2, scratch);
                else
                    _mul(T0, t2, lent2, S.poly[1], S.len[1], scratch);
                len_t2b = S.len[1] + lent2 - 1;
            }

            if (sgnS < 0)
                len_A = _poly_sub(A, T0, len_t2b, A, len_t1);
            else
                len_A = _poly_sub(A, A, len_t1, T0, len_t2b);

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

            // M = R * S_modified（a2 作为临时空间）
            _mat_mul(M, R, S, a2, scratch);
        }

        return -(sgnR * sgnS);
    }

    // --- M4: GCD 分派 ---

    // Euclid GCD（raw API）
    // G 预分配 len_a 元素
    // 前置：len_a >= len_b > 0
    void _gcd_euclid(uint64_t* G, size_t& len_G,
                     const uint64_t* A, size_t len_a,
                     const uint64_t* B, size_t len_b) const
    {
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
            std::swap(a_buf, b_buf);
            a = a_buf.data();
            la = lb;
            b = b_buf.data();
            std::memcpy(b_buf.data(), r_buf.data(), lr * sizeof(uint64_t));
            lb = lr;
        }

        std::memcpy(G, a, la * sizeof(uint64_t));
        len_G = la;
    }

    // HGCD-based GCD
    // 前置：len_a >= len_b >= GCD_CUTOFF
    void _gcd_hgcd(uint64_t* G, size_t& len_G,
                   const uint64_t* A, size_t len_a,
                   const uint64_t* B, size_t len_b) const
    {
        size_t n = len_a;

        // 工作区分配
        size_t log2n = 0;
        { size_t tmp = n; while (tmp > 1) { tmp = (tmp + 1) / 2; ++log2n; } }

        size_t hgcd_ws = 28 * n + 16 * (log2n + 1);
        size_t scratch_size = 7 * n + 1;
        size_t half = (n + 1) / 2;
        size_t mat_size = 4 * half;
        size_t total = 6 * n + hgcd_ws + scratch_size + mat_size;

        std::vector<uint64_t> workspace(total, 0);
        uint64_t* J = workspace.data();
        uint64_t* R_buf = J + n;
        uint64_t* Q_buf = R_buf + n;
        word3* W3_gcd = reinterpret_cast<word3*>(Q_buf + n);  // 3n 个 uint64_t
        uint64_t* W = Q_buf + n + 3 * n;
        uint64_t* scratch = W + hgcd_ws;
        uint64_t* mat_buf = scratch + scratch_size;

        hgcd_mat M_dummy;
        for (int i = 0; i < 4; ++i) {
            M_dummy.poly[i] = mat_buf + i * half;
            M_dummy.len[i] = 0;
        }

        // --- 初始 divrem ---
        size_t len_Q, len_R;
        std::tie(len_Q, len_R) = _poly_divrem(Q_buf, R_buf, A, len_a, B, len_b, W3_gcd);

        if (len_R == 0) {
            std::memcpy(G, B, len_b * sizeof(uint64_t));
            len_G = len_b;
            return;
        }

        // --- 首次 HGCD ---
        size_t len_g, len_j;
        _hgcd_recursive(M_dummy, false,
                        G, len_g,
                        J, len_j,
                        B, len_b,
                        R_buf, len_R,
                        W, scratch);

        // --- 主循环 ---
        while (len_j != 0) {
            std::tie(len_Q, len_R) = _poly_divrem(Q_buf, R_buf, G, len_g, J, len_j, W3_gcd);

            if (len_R == 0) {
                std::memcpy(G, J, len_j * sizeof(uint64_t));
                len_G = len_j;
                return;
            }

            if (len_j < GCD_CUTOFF) {
                _gcd_euclid(G, len_G, J, len_j, R_buf, len_R);
                return;
            }

            // HGCD：输出 J 与输入 a (= J) aliasing
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

public:
    // FLINT 归一化 Barrett 模乘（同 Zp::__nmod_mul）
    uint64_t nmod_mul(uint64_t a, uint64_t b) const
    {
        assert(_p != 0);
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

    // 模逆元（复用全局 inv_prime）
    uint64_t nmod_inv(uint64_t a) const
    {
        return inv_prime(a, _p);
    }

    // --- 构造 ---

    dense_upoly_zp() : _p(0), _ninv(0), _norm(0) {}

    explicit dense_upoly_zp(uint64_t p) : _p(p) { __precompute(); }

    // 从稀疏表示构造（稀疏 → 稠密）
    dense_upoly_zp(const upolynomial_<Zp>& sparse, uint64_t p) : _p(p)
    {
        __precompute();
        if (sparse.empty()) return;
        int64_t d = sparse.front().first.deg();  // 稀疏降幂，首项最高次
        _coeffs.resize(d + 1, 0);
        for (auto& term : sparse)
            _coeffs[term.first.deg()] = term.second.number();
    }

    // --- 查询 ---

    bool empty() const { return _coeffs.empty(); }
    int64_t deg() const { return (int64_t)_coeffs.size() - 1; }  // 空 → -1
    uint64_t lead() const { return _coeffs.empty() ? 0 : _coeffs.back(); }
    uint64_t operator[](size_t i) const { return i < _coeffs.size() ? _coeffs[i] : 0; }
    uint64_t prime() const { return _p; }
    size_t size() const { return _coeffs.size(); }

    // --- 转换：稠密 → 稀疏 ---

    upolynomial_<Zp> to_upoly() const
    {
        upolynomial_<Zp> result;
        result.reserve(_coeffs.size());
        // 降幂输出（upolynomial 约定高次在前）
        for (int64_t i = (int64_t)_coeffs.size() - 1; i >= 0; --i) {
            if (_coeffs[i] != 0)
                result.push_back({umonomial(i), Zp(_coeffs[i], _p)});
        }
        return result;
    }

    // --- 标量乘法 ---

    void scalar_mul(uint64_t c)
    {
        if (c == 0) { _coeffs.clear(); return; }
        if (c == 1) return;
        for (auto& coeff : _coeffs)
            coeff = nmod_mul(coeff, c);
        // 不需要 __strip()：c != 0 且 lead != 0 → 乘积 != 0（素数域）
    }

    // --- 乘法（Karatsuba / schoolbook 自动分派）---

    static void mul(dense_upoly_zp& C,
                    const dense_upoly_zp& A, const dense_upoly_zp& B)
    {
        assert(A._p == B._p);
        // aliasing 保护
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
            C._coeffs.resize(len_a + len_b - 1);
            C._classical_mul(C._coeffs.data(),
                             A._coeffs.data(), len_a,
                             B._coeffs.data(), len_b);
            C.__strip();
        } else {
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

    // --- 稠密长除法（惰性 3-word 归约）---

    static void divrem(dense_upoly_zp& Q, dense_upoly_zp& R,
                       const dense_upoly_zp& A, const dense_upoly_zp& B)
    {
        assert(!B.empty() && A._p == B._p);
        // aliasing 保护
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

        // 3-word 累加器
        size_t r_len = A._coeffs.size();
        std::vector<word3> R3(r_len);
        for (size_t i = 0; i < r_len; ++i)
            R3[i] = {A._coeffs[i], 0, 0};

        uint64_t p = A._p;
        uint64_t pinv = A._ninv;
        uint32_t norm = A._norm;

        for (int64_t i = q_len - 1; i >= 0; --i) {
            // 归约首项
            uint64_t r = _lll_mod_preinv(R3[i + d].hi, R3[i + d].mid,
                                          R3[i + d].lo, p, pinv, norm);
            uint64_t q_i = R.nmod_mul(r, inv_lc);
            Q._coeffs[i] = q_i;

            if (q_i != 0) {
                uint64_t c = p - q_i;   // negmod
                for (int64_t j = 0; j <= d; ++j) {
                    uint64_t p1, p0;
                    _umul128(p1, p0, c, B._coeffs[j]);
                    _add_carry3(R3[i + j], p1, p0);
                }
            }
        }

        // 最终归约：仅余式部分
        R._coeffs.resize(d > 0 ? d : 0);
        for (int64_t i = 0; i < d; ++i)
            R._coeffs[i] = _lll_mod_preinv(R3[i].hi, R3[i].mid, R3[i].lo,
                                            p, pinv, norm);
        R.__strip();
        Q.__strip();
    }

    // --- Euclid GCD ---

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
            // Euclid
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

    // --- 对象级加减法 ---

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
};

} // namespace clpoly
#endif
