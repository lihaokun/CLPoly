// 翻译器测试：模拟 CLPoly 核心函数
#include <cstdint>
#include <cassert>

// Barrett reduction 模乘（真实 CLPoly 逻辑）
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
