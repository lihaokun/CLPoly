#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

// Helper: 构造 upolynomial_<ZZ> from {(deg, coeff), ...} 降序
upolynomial_<ZZ> make_upoly_zz(std::initializer_list<std::pair<int64_t, int64_t>> terms)
{
    upolynomial_<ZZ> poly;
    for (auto& t : terms)
        if (t.second != 0)
            poly.push_back(std::make_pair(umonomial(t.first), ZZ(t.second)));
    return poly;
}

// Helper: 构造 upolynomial_<Zp>
upolynomial_<Zp> make_upoly_zp(std::initializer_list<std::pair<int64_t, uint64_t>> terms, uint32_t p)
{
    upolynomial_<Zp> poly;
    for (auto& t : terms)
        if (t.second % p != 0)
            poly.push_back(std::make_pair(umonomial(t.first), Zp(t.second, p)));
    return poly;
}

// Helper: 清除 ZZ 多项式中的零项
void clean_zeros(upolynomial_<ZZ>& f)
{
    auto it = f.data().begin();
    auto out = it;
    for (; it != f.data().end(); ++it)
        if (it->second) { if (out != it) *out = std::move(*it); ++out; }
    f.data().erase(out, f.data().end());
}

// Helper: ZZ 多项式系数 mod m 并清除零项
void mod_coeff(upolynomial_<ZZ>& f, const ZZ& m)
{
    for (auto& term : f.data())
        ZZ::fdiv_r(term.second, term.second, m);
    clean_zeros(f);
}

// Helper: ZZ 多项式乘积 mod m
upolynomial_<ZZ> product_mod(const std::vector<upolynomial_<ZZ>>& factors, const ZZ& m)
{
    upolynomial_<ZZ> result;
    result.push_back(std::make_pair(umonomial(0), ZZ(1)));
    for (auto& f : factors)
    {
        result = result * f;
        result.normalization();
        mod_coeff(result, m);
    }
    return result;
}

int main() {

    // ========================================
    // §4.2 对称模测试
    // ========================================

    CLPOLY_TEST("__symmetric_mod");
    {
        // 7 mod 10 → 7 > 5 → -3
        CLPOLY_ASSERT_EQ(__symmetric_mod(ZZ(7), ZZ(10)), ZZ(-3));
        // 3 mod 10 → 3 <= 5 → 3
        CLPOLY_ASSERT_EQ(__symmetric_mod(ZZ(3), ZZ(10)), ZZ(3));
        // -7 mod 10 → fdiv_r = 3 → 3
        CLPOLY_ASSERT_EQ(__symmetric_mod(ZZ(-7), ZZ(10)), ZZ(3));
        // 5 mod 10 → stays 5 ((-m/2, m/2])
        CLPOLY_ASSERT_EQ(__symmetric_mod(ZZ(5), ZZ(10)), ZZ(5));
        // 0 mod 7 → 0
        CLPOLY_ASSERT_EQ(__symmetric_mod(ZZ(0), ZZ(7)), ZZ(0));
        // 15 mod 7 → 1
        CLPOLY_ASSERT_EQ(__symmetric_mod(ZZ(15), ZZ(7)), ZZ(1));
        // -1 mod 7 → fdiv_r=6 → 6 > 3 → -1
        CLPOLY_ASSERT_EQ(__symmetric_mod(ZZ(-1), ZZ(7)), ZZ(-1));
    }

    CLPOLY_TEST("__upoly_symmetric_mod");
    {
        // f = 8x^2 + 6x + 11, mod 7
        auto f = make_upoly_zz({{2, 8}, {1, 6}, {0, 11}});
        auto r = __upoly_symmetric_mod(f, ZZ(7));
        // 8 mod 7 = 1, 6 mod 7 = 6 → -1, 11 mod 7 = 4 → -3
        auto expected = make_upoly_zz({{2, 1}, {1, -1}, {0, -3}});
        CLPOLY_ASSERT_EQ(r, expected);
    }

    // ========================================
    // §4.3 范数测试
    // ========================================

    CLPOLY_TEST("__upoly_norm_l2_sq");
    {
        // f = 3x^2 + 4x + 5 → 9+16+25=50
        auto f = make_upoly_zz({{2, 3}, {1, 4}, {0, 5}});
        CLPOLY_ASSERT_EQ(__upoly_norm_l2_sq(f), ZZ(50));
    }

    // ========================================
    // §4.4 Mignotte 界测试
    // ========================================

    CLPOLY_TEST("__mignotte_bound");
    {
        // f = x^2 + 5x + 6 = (x+2)(x+3)
        auto f = make_upoly_zz({{2, 1}, {1, 5}, {0, 6}});
        ZZ B = __mignotte_bound(f);
        // 因子系数最大是 3, bound 应 >= 3
        CLPOLY_ASSERT(B >= ZZ(3));
    }

    // ========================================
    // §4.1.4 模除法测试
    // ========================================

    CLPOLY_TEST("__upoly_divmod_mod");
    {
        // f = x^3 + 2x^2 + 3x + 4, g = x + 1, mod 7
        auto f = make_upoly_zz({{3, 1}, {2, 2}, {1, 3}, {0, 4}});
        auto g = make_upoly_zz({{1, 1}, {0, 1}});
        upolynomial_<ZZ> q, r;
        __upoly_divmod_mod(q, r, f, g, ZZ(7));

        // 验证 f ≡ q*g + r (mod 7)
        auto check = q * g + r;
        check.normalization();
        mod_coeff(check, ZZ(7));
        auto f_mod = f;
        mod_coeff(f_mod, ZZ(7));
        CLPOLY_ASSERT_EQ(check, f_mod);
        // deg(r) < deg(g) = 1
        CLPOLY_ASSERT(r.empty() || get_deg(r) < get_deg(g));
    }

    CLPOLY_TEST("__upoly_divmod_mod (larger)");
    {
        // f = 6x^3 + 5x^2 + 4x + 3, g = 2x^2 + 3x + 1, mod 11
        auto f = make_upoly_zz({{3, 6}, {2, 5}, {1, 4}, {0, 3}});
        auto g = make_upoly_zz({{2, 2}, {1, 3}, {0, 1}});
        upolynomial_<ZZ> q, r;
        __upoly_divmod_mod(q, r, f, g, ZZ(11));

        auto check = q * g + r;
        check.normalization();
        mod_coeff(check, ZZ(11));
        auto f_mod = f;
        mod_coeff(f_mod, ZZ(11));
        CLPOLY_ASSERT_EQ(check, f_mod);
        CLPOLY_ASSERT(r.empty() || get_deg(r) < get_deg(g));
    }

    // ========================================
    // §6.4 Hensel 树构建测试
    // ========================================

    CLPOLY_TEST("__hensel_tree_build (2 factors)");
    {
        uint32_t p = 5;
        auto f1 = make_upoly_zp({{1, 1}, {0, 1}}, p);  // x+1
        auto f2 = make_upoly_zp({{1, 1}, {0, 2}}, p);  // x+2
        std::vector<upolynomial_<Zp>> factors = {f1, f2};

        auto nodes = __hensel_tree_build(factors, p);
        CLPOLY_ASSERT(nodes.size() >= 1);

        // 验证 s*g + t*h ≡ 1 (mod p)
        auto sg = nodes[0].s * nodes[0].g;
        auto th = nodes[0].t * nodes[0].h;
        auto sum = sg + th;
        sum.normalization();
        mod_coeff(sum, ZZ(p));
        auto one = make_upoly_zz({{0, 1}});
        CLPOLY_ASSERT_EQ(sum, one);
    }

    CLPOLY_TEST("__hensel_tree_build (3 factors)");
    {
        uint32_t p = 7;
        auto f1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto f2 = make_upoly_zp({{1, 1}, {0, 2}}, p);
        auto f3 = make_upoly_zp({{1, 1}, {0, 3}}, p);
        std::vector<upolynomial_<Zp>> factors = {f1, f2, f3};

        auto nodes = __hensel_tree_build(factors, p);
        CLPOLY_ASSERT(nodes.size() >= 2);

        // 验证根节点 s*g + t*h ≡ 1 (mod p)
        auto sg = nodes[0].s * nodes[0].g;
        auto th = nodes[0].t * nodes[0].h;
        auto sum = sg + th;
        sum.normalization();
        mod_coeff(sum, ZZ(p));
        auto one = make_upoly_zz({{0, 1}});
        CLPOLY_ASSERT_EQ(sum, one);
    }

    // ========================================
    // §6.5 单步 Hensel 提升测试
    // ========================================

    CLPOLY_TEST("__hensel_step (x+1)(x+2)");
    {
        auto f = make_upoly_zz({{2, 1}, {1, 3}, {0, 2}});
        uint32_t p = 5;
        auto f1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto f2 = make_upoly_zp({{1, 1}, {0, 2}}, p);
        std::vector<upolynomial_<Zp>> factors = {f1, f2};
        auto nodes = __hensel_tree_build(factors, p);

        __hensel_step(nodes[0], f, ZZ(p));

        // f ≡ g*h (mod 25)
        auto gh = nodes[0].g * nodes[0].h;
        gh.normalization();
        mod_coeff(gh, ZZ(25));
        auto f_mod = f;
        mod_coeff(f_mod, ZZ(25));
        CLPOLY_ASSERT_EQ(gh, f_mod);

        // s*g + t*h ≡ 1 (mod 25)
        auto sg = nodes[0].s * nodes[0].g;
        auto th = nodes[0].t * nodes[0].h;
        auto sum = sg + th;
        sum.normalization();
        mod_coeff(sum, ZZ(25));
        auto one = make_upoly_zz({{0, 1}});
        CLPOLY_ASSERT_EQ(sum, one);
    }

    CLPOLY_TEST("__hensel_step (x+3)(x-5)");
    {
        // f = x^2 - 2x - 15
        auto f = make_upoly_zz({{2, 1}, {1, -2}, {0, -15}});
        uint32_t p = 7;
        // mod 7: (x+3)(x+2)
        auto f1 = make_upoly_zp({{1, 1}, {0, 3}}, p);
        auto f2 = make_upoly_zp({{1, 1}, {0, 2}}, p);
        std::vector<upolynomial_<Zp>> factors = {f1, f2};
        auto nodes = __hensel_tree_build(factors, p);

        __hensel_step(nodes[0], f, ZZ(p));

        // f ≡ g*h (mod 49)
        auto gh = nodes[0].g * nodes[0].h;
        gh.normalization();
        mod_coeff(gh, ZZ(49));
        auto f_mod = f;
        mod_coeff(f_mod, ZZ(49));
        CLPOLY_ASSERT_EQ(gh, f_mod);
    }

    // ========================================
    // §6.6 完整 Hensel 提升测试
    // ========================================

    CLPOLY_TEST("__hensel_lift (x+1)(x+2)");
    {
        auto f = make_upoly_zz({{2, 1}, {1, 3}, {0, 2}});
        uint32_t p = 5;
        auto f1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto f2 = make_upoly_zp({{1, 1}, {0, 2}}, p);
        std::vector<upolynomial_<Zp>> factors = {f1, f2};

        auto [lifted, modulus] = __hensel_lift(f, factors, p);
        CLPOLY_ASSERT_EQ(lifted.size(), (size_t)2);

        // ∏ lifted ≡ f (mod modulus)
        auto prod = product_mod(lifted, modulus);
        auto f_mod = f;
        mod_coeff(f_mod, modulus);
        CLPOLY_ASSERT_EQ(prod, f_mod);

        // 对称 mod 后恢复原因子
        auto g1 = __upoly_symmetric_mod(lifted[0], modulus);
        auto g2 = __upoly_symmetric_mod(lifted[1], modulus);
        auto product_check = g1 * g2;
        product_check.normalization();
        CLPOLY_ASSERT_EQ(product_check, f);
    }

    CLPOLY_TEST("__hensel_lift (x+1)(x-2)(x+3)");
    {
        // f = x^3 + 2x^2 - 5x - 6
        auto f = make_upoly_zz({{3, 1}, {2, 2}, {1, -5}, {0, -6}});
        uint32_t p = 7;
        auto f1 = make_upoly_zp({{1, 1}, {0, 1}}, p);  // x+1
        auto f2 = make_upoly_zp({{1, 1}, {0, 5}}, p);  // x-2 mod 7
        auto f3 = make_upoly_zp({{1, 1}, {0, 3}}, p);  // x+3
        std::vector<upolynomial_<Zp>> factors = {f1, f2, f3};

        auto [lifted, modulus] = __hensel_lift(f, factors, p);
        CLPOLY_ASSERT_EQ(lifted.size(), (size_t)3);

        // ∏ lifted ≡ f (mod modulus)
        auto prod = product_mod(lifted, modulus);
        auto f_mod = f;
        mod_coeff(f_mod, modulus);
        CLPOLY_ASSERT_EQ(prod, f_mod);

        // 对称 mod 后恢复
        upolynomial_<ZZ> total;
        total.push_back(std::make_pair(umonomial(0), ZZ(1)));
        for (auto& li : lifted)
        {
            auto gi = __upoly_symmetric_mod(li, modulus);
            total = total * gi;
            total.normalization();
        }
        CLPOLY_ASSERT_EQ(total, f);
    }

    CLPOLY_TEST("__hensel_lift (2x+3)(x+1)");
    {
        // f = 2x^2 + 5x + 3
        auto f = make_upoly_zz({{2, 2}, {1, 5}, {0, 3}});
        uint32_t p = 7;
        // mod 7 首一: (x+5)(x+1)
        auto f1 = make_upoly_zp({{1, 1}, {0, 5}}, p);
        auto f2 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        std::vector<upolynomial_<Zp>> factors = {f1, f2};

        auto [lifted, modulus] = __hensel_lift(f, factors, p);
        CLPOLY_ASSERT_EQ(lifted.size(), (size_t)2);

        // ∏ lifted ≡ f (mod modulus)
        auto prod = product_mod(lifted, modulus);
        auto f_mod = f;
        mod_coeff(f_mod, modulus);
        CLPOLY_ASSERT_EQ(prod, f_mod);
    }

    CLPOLY_TEST("__hensel_lift 4 factors");
    {
        // f = (x+1)(x-1)(x+2)(x-3) = x^4 - x^3 - 7x^2 + x + 6
        auto f = make_upoly_zz({{4, 1}, {3, -1}, {2, -7}, {1, 1}, {0, 6}});
        uint32_t p = 11;
        auto f1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto f2 = make_upoly_zp({{1, 1}, {0, 10}}, p); // x-1
        auto f3 = make_upoly_zp({{1, 1}, {0, 2}}, p);
        auto f4 = make_upoly_zp({{1, 1}, {0, 8}}, p);  // x-3
        std::vector<upolynomial_<Zp>> factors = {f1, f2, f3, f4};

        auto [lifted, modulus] = __hensel_lift(f, factors, p);
        CLPOLY_ASSERT_EQ(lifted.size(), (size_t)4);

        // 对称 mod 后恢复
        upolynomial_<ZZ> total;
        total.push_back(std::make_pair(umonomial(0), ZZ(1)));
        for (auto& li : lifted)
        {
            auto gi = __upoly_symmetric_mod(li, modulus);
            total = total * gi;
            total.normalization();
        }
        CLPOLY_ASSERT_EQ(total, f);
    }

    CLPOLY_TEST("__hensel_lift big coefficients");
    {
        // f = (x+100)(x-200) = x^2 - 100x - 20000
        auto f = make_upoly_zz({{2, 1}, {1, -100}, {0, -20000}});
        uint32_t p = 11;
        // mod 11: -100 mod 11 = 10, -20000 mod 11 = 9
        // f ≡ x^2 + 10x + 9 = (x+1)(x+9) mod 11
        auto g1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto g2 = make_upoly_zp({{1, 1}, {0, 9}}, p);
        std::vector<upolynomial_<Zp>> factors = {g1, g2};

        auto [lifted, modulus] = __hensel_lift(f, factors, p);
        CLPOLY_ASSERT_EQ(lifted.size(), (size_t)2);

        // 对称 mod 后恢复
        auto sym1 = __upoly_symmetric_mod(lifted[0], modulus);
        auto sym2 = __upoly_symmetric_mod(lifted[1], modulus);
        auto prod = sym1 * sym2;
        prod.normalization();
        CLPOLY_ASSERT_EQ(prod, f);
    }

    // ========================================
    // ZZ 新增功能测试
    // ========================================

    CLPOLY_TEST("ZZ::fdiv_r");
    {
        ZZ r;
        ZZ::fdiv_r(r, ZZ(17), ZZ(5));
        CLPOLY_ASSERT_EQ(r, ZZ(2));
        ZZ::fdiv_r(r, ZZ(-17), ZZ(5));
        CLPOLY_ASSERT_EQ(r, ZZ(3));
        ZZ::fdiv_r(r, ZZ(0), ZZ(5));
        CLPOLY_ASSERT_EQ(r, ZZ(0));
        ZZ::fdiv_r(r, ZZ(15), ZZ(5));
        CLPOLY_ASSERT_EQ(r, ZZ(0));
    }

    CLPOLY_TEST("ZZ::invert");
    {
        ZZ result;
        // 3^{-1} mod 7 = 5
        bool ok = ZZ::invert(result, ZZ(3), ZZ(7));
        CLPOLY_ASSERT(ok);
        CLPOLY_ASSERT_EQ(result, ZZ(5));

        // 2^{-1} mod 11 = 6
        ok = ZZ::invert(result, ZZ(2), ZZ(11));
        CLPOLY_ASSERT(ok);
        CLPOLY_ASSERT_EQ(result, ZZ(6));

        // 大模数: 3^{-1} mod 49, 验证 3*result mod 49 = 1
        ok = ZZ::invert(result, ZZ(3), ZZ(49));
        CLPOLY_ASSERT(ok);
        ZZ check = ZZ(3) * result;
        ZZ::fdiv_r(check, check, ZZ(49));
        CLPOLY_ASSERT_EQ(check, ZZ(1));
    }

    return clpoly_test::test_summary();
}
