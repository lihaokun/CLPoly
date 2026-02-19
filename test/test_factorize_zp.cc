#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

// Helper: 构造 upolynomial_<Zp> from {(deg, coeff), ...} 降序
upolynomial_<Zp> make_upoly_zp(std::initializer_list<std::pair<int64_t, uint64_t>> terms, uint32_t p)
{
    upolynomial_<Zp> poly;
    for (auto& t : terms)
        if (t.second % p != 0)
            poly.push_back(std::make_pair(umonomial(t.first), Zp(t.second, p)));
    return poly;
}



// Helper: 验证 product of factors == f in Zp[x]
bool verify_factorization_Zp(
    const upolynomial_<Zp>& f,
    Zp lc,
    const std::vector<std::pair<upolynomial_<Zp>, uint64_t>>& factors)
{
    uint32_t p = f.front().second.prime();
    upolynomial_<Zp> product({std::make_pair(umonomial(0), lc)});
    for (auto& fac : factors)
    {
        auto fi_pow = pow(fac.first, (int64_t)fac.second);
        product = product * fi_pow;
        product.normalization();
    }
    return product == f;
}

// Helper: 验证 Zp 多项式在 Zp 上是否不可约
bool is_irreducible_Zp(const upolynomial_<Zp>& f)
{
    if (get_deg(f) <= 1) return true;
    auto f_copy = f;
    __upoly_make_monic(f_copy);
    auto ddf = __ddf_Zp(f_copy);
    return ddf.size() == 1 && get_deg(ddf[0].first) == get_deg(f);
}

int main() {

    // ========================================
    // 辅助函数测试
    // ========================================

    CLPOLY_TEST("__upoly_make_monic");
    {
        uint32_t p = 7;
        auto f = make_upoly_zp({{2, 3}, {1, 2}, {0, 5}}, p);
        Zp lc = __upoly_make_monic(f);
        CLPOLY_ASSERT_EQ(lc.number(), (uint64_t)3);
        CLPOLY_ASSERT_EQ(f.front().second.number(), (uint64_t)1);
        auto rebuilt = f;
        for (auto& t : rebuilt) t.second *= lc;
        auto orig = make_upoly_zp({{2, 3}, {1, 2}, {0, 5}}, p);
        CLPOLY_ASSERT_EQ(rebuilt, orig);
    }

    CLPOLY_TEST("__upoly_mod");
    {
        uint32_t p = 7;
        auto f = make_upoly_zp({{3, 1}, {1, 2}, {0, 1}}, p);
        auto g = make_upoly_zp({{2, 1}, {0, 1}}, p);
        auto r = __upoly_mod(f, g);
        auto expected = make_upoly_zp({{1, 1}, {0, 1}}, p);
        CLPOLY_ASSERT_EQ(r, expected);
    }

    CLPOLY_TEST("__upoly_divmod");
    {
        uint32_t p = 5;
        auto f = make_upoly_zp({{3, 1}, {0, 1}}, p);
        auto g = make_upoly_zp({{1, 1}, {0, 1}}, p);
        upolynomial_<Zp> q, r;
        __upoly_divmod(q, r, f, g);
        auto expected_q = make_upoly_zp({{2, 1}, {1, 4}, {0, 1}}, p);
        CLPOLY_ASSERT(r.empty());
        CLPOLY_ASSERT_EQ(q, expected_q);
    }

    CLPOLY_TEST("__upoly_gcd_Zp");
    {
        uint32_t p = 7;
        auto a = make_upoly_zp({{2, 1}, {0, 6}}, p);
        auto b = make_upoly_zp({{1, 1}, {0, 6}}, p);
        auto g = polynomial_GCD(a, b);
        auto expected = make_upoly_zp({{1, 1}, {0, 6}}, p);
        CLPOLY_ASSERT_EQ(g, expected);
    }

    CLPOLY_TEST("__upoly_gcd_extended");
    {
        uint32_t p = 7;
        auto a = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto b = make_upoly_zp({{1, 1}, {0, 2}}, p);
        upolynomial_<Zp> s, t;
        auto g = polynomial_GCD(a, b, s, t);
        auto sa = s * a;
        auto tb = t * b;
        auto sum = sa + tb;
        sum.normalization();
        CLPOLY_ASSERT_EQ(sum, g);
    }

    CLPOLY_TEST("__upoly_gcd_extended_2");
    {
        uint32_t p = 13;
        auto a = make_upoly_zp({{2, 1}, {0, 1}}, p);
        auto b = make_upoly_zp({{1, 1}, {0, 3}}, p);
        upolynomial_<Zp> s, t;
        auto g = polynomial_GCD(a, b, s, t);
        auto sa = s * a;
        auto tb = t * b;
        auto sum = sa + tb;
        sum.normalization();
        CLPOLY_ASSERT_EQ(sum, g);
    }

    CLPOLY_TEST("__upoly_powmod");
    {
        uint32_t p = 5;
        auto base = make_upoly_zp({{1, 1}}, p);
        auto mod = make_upoly_zp({{2, 1}, {0, 1}}, p);
        auto result = __upoly_powmod(base, ZZ(5), mod);
        auto expected = make_upoly_zp({{1, 1}}, p);
        CLPOLY_ASSERT_EQ(result, expected);
    }

    CLPOLY_TEST("__upoly_powmod_large");
    {
        uint32_t p = 7;
        auto base = make_upoly_zp({{1, 1}}, p);
        auto mod = make_upoly_zp({{3, 1}, {1, 1}, {0, 1}}, p);
        auto result = __upoly_powmod(base, ZZ(7), mod);
        CLPOLY_ASSERT(get_deg(result) < 3);
        auto expected = make_upoly_zp({{2, 2}, {0, 6}}, p);
        CLPOLY_ASSERT_EQ(result, expected);
    }

    // ========================================
    // 无平方分解测试
    // ========================================

    CLPOLY_TEST("__squarefree_Zp_already_squarefree");
    {
        uint32_t p = 7;
        auto f = make_upoly_zp({{2, 1}, {1, 1}, {0, 1}}, p);
        __upoly_make_monic(f);
        auto sqf = __squarefree_Zp(f);
        CLPOLY_ASSERT_EQ((int)sqf.size(), 1);
        CLPOLY_ASSERT_EQ(sqf[0].second, (uint64_t)1);
    }

    CLPOLY_TEST("__squarefree_Zp_with_square");
    {
        uint32_t p = 7;
        auto xp1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto xp2 = make_upoly_zp({{1, 1}, {0, 2}}, p);
        auto f = xp1 * xp1 * xp2;
        f.normalization();
        __upoly_make_monic(f);
        auto sqf = __squarefree_Zp(f);
        // 验证乘积
        uint32_t pp = p;
        upolynomial_<Zp> product({std::make_pair(umonomial(0), Zp((int64_t)1, pp))});
        for (auto& si_ei : sqf)
        {
            product = product * pow(si_ei.first, (int64_t)si_ei.second);
            product.normalization();
        }
        CLPOLY_ASSERT_EQ(product, f);
        CLPOLY_ASSERT(sqf.size() >= 2);
    }

    CLPOLY_TEST("__squarefree_Zp_pth_power");
    {
        uint32_t p = 3;
        auto xp1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto f = pow(xp1, 3);
        __upoly_make_monic(f);
        auto sqf = __squarefree_Zp(f);
        CLPOLY_ASSERT_EQ((int)sqf.size(), 1);
        CLPOLY_ASSERT_EQ(sqf[0].second, (uint64_t)3);
        CLPOLY_ASSERT_EQ(get_deg(sqf[0].first), (int64_t)1);
    }

    // ========================================
    // DDF 测试
    // ========================================

    CLPOLY_TEST("__ddf_Zp_all_linear");
    {
        uint32_t p = 7;
        auto f1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto f2 = make_upoly_zp({{1, 1}, {0, 2}}, p);
        auto f3 = make_upoly_zp({{1, 1}, {0, 3}}, p);
        auto f = f1 * f2 * f3;
        f.normalization();
        __upoly_make_monic(f);
        auto ddf = __ddf_Zp(f);
        CLPOLY_ASSERT_EQ((int)ddf.size(), 1);
        CLPOLY_ASSERT_EQ(ddf[0].second, (uint64_t)1);
        CLPOLY_ASSERT_EQ(get_deg(ddf[0].first), (int64_t)3);
    }

    CLPOLY_TEST("__ddf_Zp_irreducible_deg2");
    {
        uint32_t p = 5;
        auto f = make_upoly_zp({{2, 1}, {0, 2}}, p);
        auto ddf = __ddf_Zp(f);
        CLPOLY_ASSERT_EQ((int)ddf.size(), 1);
        CLPOLY_ASSERT_EQ(ddf[0].second, (uint64_t)2);
    }

    CLPOLY_TEST("__ddf_Zp_mixed_degrees");
    {
        uint32_t p = 5;
        auto f1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto f2 = make_upoly_zp({{2, 1}, {0, 2}}, p);
        auto f = f1 * f2;
        f.normalization();
        __upoly_make_monic(f);
        auto ddf = __ddf_Zp(f);
        CLPOLY_ASSERT_EQ((int)ddf.size(), 2);
        CLPOLY_ASSERT_EQ(ddf[0].second, (uint64_t)1);
        CLPOLY_ASSERT_EQ(ddf[1].second, (uint64_t)2);
        CLPOLY_ASSERT_EQ(get_deg(ddf[0].first), (int64_t)1);
        CLPOLY_ASSERT_EQ(get_deg(ddf[1].first), (int64_t)2);
    }

    // ========================================
    // EDF 测试
    // ========================================

    CLPOLY_TEST("__edf_Zp_split_linear");
    {
        uint32_t p = 7;
        auto f1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto f2 = make_upoly_zp({{1, 1}, {0, 2}}, p);
        auto f3 = make_upoly_zp({{1, 1}, {0, 3}}, p);
        auto f = f1 * f2 * f3;
        f.normalization();
        __upoly_make_monic(f);

        std::mt19937 rng(42);
        std::vector<upolynomial_<Zp>> factors;
        __edf_Zp(factors, f, 1, rng);

        CLPOLY_ASSERT_EQ((int)factors.size(), 3);
        upolynomial_<Zp> product({std::make_pair(umonomial(0), Zp((int64_t)1, p))});
        for (auto& fi : factors)
        {
            product = product * fi;
            product.normalization();
        }
        CLPOLY_ASSERT_EQ(product, f);
    }

    CLPOLY_TEST("__edf_Zp_split_deg2");
    {
        uint32_t p = 5;
        auto f1 = make_upoly_zp({{2, 1}, {0, 2}}, p);
        auto f2 = make_upoly_zp({{2, 1}, {0, 3}}, p);
        auto f = f1 * f2;
        f.normalization();
        __upoly_make_monic(f);

        std::mt19937 rng(42);
        std::vector<upolynomial_<Zp>> factors;
        __edf_Zp(factors, f, 2, rng);

        CLPOLY_ASSERT_EQ((int)factors.size(), 2);
        upolynomial_<Zp> product({std::make_pair(umonomial(0), Zp((int64_t)1, p))});
        for (auto& fi : factors)
        {
            product = product * fi;
            product.normalization();
        }
        CLPOLY_ASSERT_EQ(product, f);
    }

    // ========================================
    // __factor_Zp 完整分解测试
    // ========================================

    CLPOLY_TEST("__factor_Zp_irreducible");
    {
        uint32_t p = 5;
        auto f = make_upoly_zp({{2, 1}, {0, 2}}, p);
        auto result = __factor_Zp(f);
        CLPOLY_ASSERT_EQ(result.first.number(), (uint64_t)1);
        CLPOLY_ASSERT_EQ((int)result.second.size(), 1);
        CLPOLY_ASSERT_EQ(result.second[0].second, (uint64_t)1);
        CLPOLY_ASSERT_EQ(result.second[0].first, f);
    }

    CLPOLY_TEST("__factor_Zp_linear_factors");
    {
        uint32_t p = 7;
        auto f1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto f2 = make_upoly_zp({{1, 1}, {0, 2}}, p);
        auto f3 = make_upoly_zp({{1, 1}, {0, 4}}, p);
        auto f = f1 * f2 * f3;
        f.normalization();
        // 乘以 3
        for (auto& t : f) t.second *= Zp((int64_t)3, p);
        f.normalization();

        auto result = __factor_Zp(f);
        CLPOLY_ASSERT_EQ(result.first.number(), (uint64_t)3);
        CLPOLY_ASSERT_EQ((int)result.second.size(), 3);
        CLPOLY_ASSERT(verify_factorization_Zp(f, result.first, result.second));
    }

    CLPOLY_TEST("__factor_Zp_with_multiplicity");
    {
        uint32_t p = 7;
        auto xp1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto xp2 = make_upoly_zp({{1, 1}, {0, 2}}, p);
        auto f = pow(xp1, 3) * xp2;
        f.normalization();

        auto result = __factor_Zp(f);
        CLPOLY_ASSERT(verify_factorization_Zp(f, result.first, result.second));
        CLPOLY_ASSERT_EQ((int)result.second.size(), 2);
    }

    CLPOLY_TEST("__factor_Zp_mixed");
    {
        uint32_t p = 5;
        auto f1 = make_upoly_zp({{1, 1}, {0, 1}}, p);
        auto f2 = make_upoly_zp({{2, 1}, {0, 2}}, p);
        auto f = f1 * f2;
        f.normalization();

        auto result = __factor_Zp(f);
        CLPOLY_ASSERT_EQ((int)result.second.size(), 2);
        CLPOLY_ASSERT(verify_factorization_Zp(f, result.first, result.second));
        CLPOLY_ASSERT_EQ(get_deg(result.second[0].first), (int64_t)1);
        CLPOLY_ASSERT_EQ(get_deg(result.second[1].first), (int64_t)2);
    }

    CLPOLY_TEST("__factor_Zp_constant");
    {
        uint32_t p = 5;
        auto f = make_upoly_zp({{0, 3}}, p);
        auto result = __factor_Zp(f);
        CLPOLY_ASSERT_EQ(result.first.number(), (uint64_t)3);
        CLPOLY_ASSERT((int)result.second.size() == 0);
    }

    CLPOLY_TEST("__factor_Zp_complete_split");
    {
        uint32_t p = 5;
        // x^5 - x = x^5 + 4x in Z5
        auto f = make_upoly_zp({{5, 1}, {1, 4}}, p);
        auto result = __factor_Zp(f);
        CLPOLY_ASSERT_EQ((int)result.second.size(), 5);
        CLPOLY_ASSERT(verify_factorization_Zp(f, result.first, result.second));
        for (auto& fac : result.second)
        {
            CLPOLY_ASSERT_EQ(get_deg(fac.first), (int64_t)1);
            CLPOLY_ASSERT_EQ(fac.second, (uint64_t)1);
        }
    }

    CLPOLY_TEST("__factor_Zp_larger_prime");
    {
        uint32_t p = 101;
        auto f1 = make_upoly_zp({{1, 1}, {0, 50}}, p);
        auto f2 = make_upoly_zp({{1, 1}, {0, 99}}, p);
        auto f3 = make_upoly_zp({{2, 1}, {1, 1}, {0, 1}}, p);
        auto f = f1 * f2 * f3;
        f.normalization();

        auto result = __factor_Zp(f);
        CLPOLY_ASSERT(verify_factorization_Zp(f, result.first, result.second));
        for (auto& fac : result.second)
            CLPOLY_ASSERT(is_irreducible_Zp(fac.first));
    }

    return clpoly_test::test_summary();
}
