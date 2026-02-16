#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

// Helper: 构造 upolynomial_<ZZ>
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

// 验证因子分解: ∏ factors = f (精确)
bool verify_ZZ_factorization(const upolynomial_<ZZ>& f,
                              const std::vector<upolynomial_<ZZ>>& factors)
{
    upolynomial_<ZZ> product;
    product.push_back(std::make_pair(umonomial(0), ZZ(1)));
    for (auto& fi : factors)
    {
        product = product * fi;
        product.normalization();
    }
    return product == f;
}

// 端到端测试辅助: 从 f 出发, 做 Zp 分解 → Hensel 提升 → 因子重组
std::vector<upolynomial_<ZZ>> full_factor_pipeline(
    const upolynomial_<ZZ>& f, uint32_t p)
{
    // 1. 转换到 Zp 并分解
    auto f_zp = polynomial_mod(f, p);
    auto [lc_zp, factors_zp] = __factor_Zp(f_zp);

    // 提取首一因子
    std::vector<upolynomial_<Zp>> monic_factors;
    for (auto& fi_ei : factors_zp)
        monic_factors.push_back(fi_ei.first);

    if (monic_factors.size() <= 1)
        return {f};  // 不可约

    // 2. Hensel 提升
    auto [lifted, modulus] = __hensel_lift(f, monic_factors, p);

    // 3. 因子重组
    return __factor_recombine(f, lifted, modulus);
}

int main() {

    // ========================================
    // 辅助函数测试
    // ========================================

    CLPOLY_TEST("__upoly_norm_l1");
    {
        // f = 3x^2 - 4x + 5 → |3|+|-4|+|5| = 12
        auto f = make_upoly_zz({{2, 3}, {1, -4}, {0, 5}});
        CLPOLY_ASSERT_EQ(__upoly_norm_l1(f), ZZ(12));
    }

    CLPOLY_TEST("__upoly_primitive");
    {
        // f = 6x^2 + 9x + 3 → cont=3, pp = 2x^2 + 3x + 1
        auto f = make_upoly_zz({{2, 6}, {1, 9}, {0, 3}});
        auto [c, pp] = __upoly_primitive(std::move(f));
        CLPOLY_ASSERT_EQ(c, ZZ(3));
        auto expected = make_upoly_zz({{2, 2}, {1, 3}, {0, 1}});
        CLPOLY_ASSERT_EQ(pp, expected);

        // 负首项: f = -4x + 2 → cont=-2, pp = 2x - 1
        auto f2 = make_upoly_zz({{1, -4}, {0, 2}});
        auto [c2, pp2] = __upoly_primitive(std::move(f2));
        CLPOLY_ASSERT_EQ(c2, ZZ(-2));
        auto expected2 = make_upoly_zz({{1, 2}, {0, -1}});
        CLPOLY_ASSERT_EQ(pp2, expected2);
    }

    // ========================================
    // 简单因子重组测试
    // ========================================

    CLPOLY_TEST("recombine (x+1)(x+2)");
    {
        // f = x^2 + 3x + 2 = (x+1)(x+2)
        auto f = make_upoly_zz({{2, 1}, {1, 3}, {0, 2}});
        auto result = full_factor_pipeline(f, 5);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine (x+1)(x-2)(x+3)");
    {
        // f = x^3 + 2x^2 - 5x - 6
        auto f = make_upoly_zz({{3, 1}, {2, 2}, {1, -5}, {0, -6}});
        auto result = full_factor_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)3);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine irreducible (x^2+1)");
    {
        // f = x^2 + 1, 在 Z[x] 上不可约
        // mod 5: x^2+1 = (x+2)(x+3) → 提升后重组应发现不可分
        auto f = make_upoly_zz({{2, 1}, {0, 1}});
        auto result = full_factor_pipeline(f, 5);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(result[0], f);
    }

    CLPOLY_TEST("recombine (x+1)(x-1)(x+2)(x-3)");
    {
        // f = x^4 - x^3 - 7x^2 + x + 6
        auto f = make_upoly_zz({{4, 1}, {3, -1}, {2, -7}, {1, 1}, {0, 6}});
        auto result = full_factor_pipeline(f, 11);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)4);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine (x+100)(x-200)");
    {
        // f = x^2 - 100x - 20000
        auto f = make_upoly_zz({{2, 1}, {1, -100}, {0, -20000}});
        auto result = full_factor_pipeline(f, 11);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine (2x+3)(x+1)");
    {
        // f = 2x^2 + 5x + 3
        auto f = make_upoly_zz({{2, 2}, {1, 5}, {0, 3}});
        auto result = full_factor_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine (x^2+x+1)(x^2-x+1)");
    {
        // f = x^4 + x^2 + 1
        // = (x^2+x+1)(x^2-x+1)
        auto f = make_upoly_zz({{4, 1}, {2, 1}, {0, 1}});
        auto result = full_factor_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine (x^2+1)(x^2+2)");
    {
        // f = x^4 + 3x^2 + 2 = (x^2+1)(x^2+2)
        auto f = make_upoly_zz({{4, 1}, {2, 3}, {0, 2}});
        // 需要找一个素数使得两个二次因子 mod p 都能分裂
        // mod 5: x^2+1 = (x+2)(x+3), x^2+2 = (x+1*? no)
        // mod 13: x^2+1 = (x+5)(x+8), x^2+2 = (x+4)(x+9)
        auto result = full_factor_pipeline(f, 13);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    // 需要重组的情况: mod p 分裂得更多, 但 Z 上因子更少
    CLPOLY_TEST("recombine x^4+1 (irreducible over Z, splits mod p)");
    {
        // x^4+1 在 Z[x] 上不可约
        // 但 mod 任何奇素数 p, 它都分裂
        // mod 5: x^4+1 = (x^2+2x+4)(x^2+3x+4) mod 5? 验证:
        // (x^2+2x+4)(x^2+3x+4) = x^4 + 3x^3 + 4x^2 + 2x^3 + 6x^2 + 8x + 4x^2 + 12x + 16
        // = x^4 + 5x^3 + 14x^2 + 20x + 16 ≡ x^4 + 0 + 4x^2 + 0 + 1 mod 5
        // 那是 x^4+4x^2+1 ≠ x^4+1, 不对
        // 实际: mod 5, x^4+1 的根: x^4 ≡ -1 ≡ 4 mod 5
        // 2^4 = 16 ≡ 1, 3^4 = 81 ≡ 1, 不是 4. 没有根
        // 但可能分裂为两个二次: 用 CZ 分解
        // 让我换个素数. mod 17: 2^4=16≡-1, 所以 x=2 是根
        // x^4+1 mod 17 有根 2, 也有根 -2=15, 也有 2^{-1}=9 (因为 9*2=18≡1), 15^{-1}=8
        // 所以 4 个根: 2, 8, 9, 15
        // x^4+1 = (x-2)(x-8)(x-9)(x-15) mod 17
        auto f = make_upoly_zz({{4, 1}, {0, 1}});
        auto result = full_factor_pipeline(f, 17);
        // x^4+1 在 Z 上不可约 (但注意: 实际上 x^4+1 = (x^2+√2 x+1)(x^2-√2 x+1)
        // 在 Z 上不可约!)
        CLPOLY_ASSERT_EQ(result.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(result[0], f);
    }

    CLPOLY_TEST("recombine (3x^2+2x+1)(x+5)");
    {
        // f = 3x^3 + 17x^2 + 11x + 5
        auto factor1 = make_upoly_zz({{2, 3}, {1, 2}, {0, 1}});
        auto factor2 = make_upoly_zz({{1, 1}, {0, 5}});
        auto f = factor1 * factor2;
        f.normalization();
        auto result = full_factor_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine Swinnerton-Dyer (x^2-2)(x^2-3)+(...)");
    {
        // f = (x^2-2)(x^2-3) = x^4 - 5x^2 + 6
        auto f = make_upoly_zz({{4, 1}, {2, -5}, {0, 6}});
        auto result = full_factor_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));

        // 检查因子是 x^2-2 和 x^2-3
        auto expected1 = make_upoly_zz({{2, 1}, {0, -2}});
        auto expected2 = make_upoly_zz({{2, 1}, {0, -3}});
        CLPOLY_ASSERT(
            (result[0] == expected1 && result[1] == expected2) ||
            (result[0] == expected2 && result[1] == expected1));
    }

    return clpoly_test::test_summary();
}
