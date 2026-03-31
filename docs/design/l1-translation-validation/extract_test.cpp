// 最小测试：从 CLPoly 提取 __nmod_mul 供 AST 分析
#include <cstdint>
#include <cassert>

// 简化版 Barrett 模乘（从 number.hh 提取核心逻辑）
uint64_t nmod_mul(uint64_t a, uint64_t b, uint64_t p, uint64_t ninv, uint32_t norm) {
    assert(p != 0);
    uint64_t pn = p << norm;
    uint64_t a_shifted = a << norm;
    unsigned __int128 prod = (unsigned __int128)a_shifted * b;
    uint64_t hi = (uint64_t)(prod >> 64);
    uint64_t lo = (uint64_t)prod;
    unsigned __int128 qm = (unsigned __int128)hi * ninv;
    uint64_t q1 = (uint64_t)(qm >> 64);
    uint64_t q0 = (uint64_t)qm;
    q0 += lo;
    q1 += hi + (q0 < lo ? 1 : 0);
    uint64_t r = lo - (q1 + 1) * pn;
    if (r > q0) r += pn;
    if (r >= pn) r -= pn;
    return r >> norm;
}

// 简化版多项式 make_monic（测试循环翻译）
struct Term { uint64_t deg; uint64_t coeff; };

void upoly_make_monic(Term* poly, int n, uint64_t p, uint64_t ninv, uint32_t norm) {
    if (n == 0) return;
    uint64_t lc = poly[0].coeff;
    if (lc == 1) return;
    uint64_t lc_inv = 0; // 简化：假设已计算
    for (int i = 0; i < n; i++) {
        poly[i].coeff = nmod_mul(poly[i].coeff, lc_inv, p, ninv, norm);
    }
}
