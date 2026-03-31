// Phase 2 测试：模拟 CLPoly 的多项式类型和操作
#include <cstdint>
#include <vector>
#include <utility>
#include <cassert>

// 简化版 Zp 系数
struct Zp {
    uint64_t val;
    uint64_t p;
    uint64_t prime() const { return p; }
};

// 简化版单项式
struct umonomial {
    uint64_t d;
    umonomial() : d(0) {}
    explicit umonomial(uint64_t deg) : d(deg) {}
    uint64_t deg() const { return d; }
};

// 简化版稀疏多项式 = vector<pair<umonomial, Zp>>
using SparsePolyZp = std::vector<std::pair<umonomial, Zp>>;

uint64_t get_deg(const SparsePolyZp& f) {
    if (f.empty()) return 0;
    return f.front().first.deg();
}

// 模拟 CLPoly 的 p 次根提取
SparsePolyZp extract_pth_root(const SparsePolyZp& f) {
    uint64_t p = f.front().second.prime();
    SparsePolyZp g;
    g.reserve(f.size());
    for (const auto& term : f) {
        assert(term.first.deg() % p == 0);
        g.push_back(std::make_pair(umonomial(term.first.deg() / p), term.second));
    }
    return g;
}

// 模拟 CLPoly 的 make_monic
void upoly_make_monic(SparsePolyZp& f, uint64_t p) {
    if (f.empty()) return;
    uint64_t lc = f.front().second.val;
    if (lc == 1) return;
    // 简化：假设 lc_inv 已知
    uint64_t lc_inv = 1; // placeholder
    for (auto& term : f) {
        term.second.val = (term.second.val * lc_inv) % p;
    }
}

// 模拟 Yun while 循环核心
void yun_loop(SparsePolyZp& w, SparsePolyZp& c,
              std::vector<std::pair<SparsePolyZp, uint64_t>>& result,
              uint64_t p) {
    uint64_t i = 1;
    while (!w.empty() && get_deg(w) > 0) {
        // 简化：y = gcd(w, c)，z = w / y
        SparsePolyZp y = w; // placeholder
        SparsePolyZp z = w; // placeholder

        if (!z.empty() && get_deg(z) > 0) {
            upoly_make_monic(z, p);
            result.push_back({std::move(z), i});
        }

        c = std::move(y);
        w = c; // placeholder
        ++i;
    }
}
