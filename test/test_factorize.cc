#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include <set>

using namespace clpoly;

// Helper: 验证 factorization 重组回原多项式
template<class Poly>
bool verify_factorization(const Poly& f, const factorization<Poly>& fac)
{
    Poly product(f.comp_ptr());
    product.push_back({typename Poly::monomial_type(f.comp_ptr()), fac.content});
    for (auto& [fi, ei] : fac.factors)
        product = product * pow(fi, (int64_t)ei);
    product.normalization();
    return product == f;
}

int main()
{
    variable x("x");
    polynomial_ZZ f, g;

    // ========================================
    // 零 / 常数 / 线性 边界情况
    // ========================================
    CLPOLY_TEST("factorize_zero");
    {
        polynomial_ZZ zero;
        auto fac = factorize(zero);
        CLPOLY_ASSERT_EQ(fac.content, ZZ(0));
        CLPOLY_ASSERT(fac.factors.empty());
    }

    CLPOLY_TEST("factorize_constant");
    {
        polynomial_ZZ c = polynomial_ZZ(ZZ(42));
        auto fac = factorize(c);
        CLPOLY_ASSERT_EQ(fac.content, ZZ(42));
        CLPOLY_ASSERT(fac.factors.empty());
    }

    CLPOLY_TEST("factorize_linear");
    {
        f = 3*x + 5;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
    }

    // ========================================
    // 不可约多项式
    // ========================================
    CLPOLY_TEST("factorize_irreducible_x2p1");
    {
        f = pow(x,2) + 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.content, ZZ(1));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)1);
    }

    CLPOLY_TEST("factorize_irreducible_cyclotomic6");
    {
        // Φ₆ = x² - x + 1
        f = pow(x,2) - x + 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
    }

    CLPOLY_TEST("factorize_irreducible_x4p1");
    {
        // x^4+1 is irreducible over Z
        f = pow(x,4) + 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
    }

    CLPOLY_TEST("factorize_irreducible_x3m2");
    {
        f = pow(x,3) - 2;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
    }

    // ========================================
    // 完全分裂
    // ========================================
    CLPOLY_TEST("factorize_x2m1");
    {
        f = pow(x,2) - 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)2);
    }

    CLPOLY_TEST("factorize_x4m1");
    {
        // (x-1)(x+1)(x^2+1)
        f = pow(x,4) - 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)3);
    }

    CLPOLY_TEST("factorize_4_linear_factors");
    {
        // (x-1)(x-2)(x-3)(x-4)
        f = (x - 1) * (x - 2) * (x - 3) * (x - 4);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)4);
        for (auto& [fi, ei] : fac.factors)
            CLPOLY_ASSERT_EQ(ei, (uint64_t)1);
    }

    CLPOLY_TEST("factorize_x6m1");
    {
        // (x-1)(x+1)(x^2-x+1)(x^2+x+1)
        f = pow(x,6) - 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)4);
    }

    // ========================================
    // 含重因子
    // ========================================
    CLPOLY_TEST("factorize_with_multiplicity");
    {
        // (x+1)^2 * (x-2)^3
        f = pow(x + 1, 2) * pow(x - 2, 3);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)2);
        // 检查 multiplicities
        std::set<uint64_t> mults;
        for (auto& [fi, ei] : fac.factors)
            mults.insert(ei);
        CLPOLY_ASSERT(mults.count(2));
        CLPOLY_ASSERT(mults.count(3));
    }

    // ========================================
    // 非首一多项式
    // ========================================
    CLPOLY_TEST("factorize_nonmonic");
    {
        // 6x^2 + 5x + 1 = (2x+1)(3x+1), content = 1
        f = 6*pow(x,2) + 5*x + 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)2);
    }

    CLPOLY_TEST("factorize_with_content");
    {
        // 12x^3 - 4x^2 + 6x - 2 = 2(2x^2+1)(3x-1)
        f = 12*pow(x,3) - 4*pow(x,2) + 6*x - 2;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)2);
    }

    // ========================================
    // x^8-1
    // ========================================
    CLPOLY_TEST("factorize_x8m1");
    {
        f = pow(x,8) - 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        // (x-1)(x+1)(x^2+1)(x^4+1)
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)4);
    }

    // ========================================
    // Swinnerton-Dyer x^4 - 10x^2 + 1 (irreducible)
    // ========================================
    CLPOLY_TEST("factorize_swinnerton_dyer");
    {
        f = pow(x,4) - 10*pow(x,2) + 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
    }

    // ========================================
    // 负首项系数
    // ========================================
    CLPOLY_TEST("factorize_negative_leading");
    {
        f = -(pow(x,2) - 1);  // -x^2 + 1
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
    }

    // ========================================
    // QQ[x] 分解
    // ========================================
    CLPOLY_TEST("factorize_QQ_simple");
    {
        // (1/2)(x^2 - 1) = (1/2)(x-1)(x+1)
        polynomial_ZZ fz = pow(x,2) - 1;
        polynomial_QQ fq;
        poly_convert(fz, fq);
        fq = fq * QQ(1,2);
        auto fac = factorize(fq);
        CLPOLY_ASSERT(verify_factorization(fq, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)2);
        // 每个因子应是首一的
        for (auto& [fi, ei] : fac.factors)
            CLPOLY_ASSERT_EQ(fi.front().second, QQ(1));
    }

    CLPOLY_TEST("factorize_QQ_irreducible");
    {
        polynomial_ZZ fz = pow(x,2) + 1;
        polynomial_QQ fq;
        poly_convert(fz, fq);
        auto fac = factorize(fq);
        CLPOLY_ASSERT(verify_factorization(fq, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
    }

    // ========================================
    // grlex ordering (通用 comp 包装)
    // ========================================
    CLPOLY_TEST("factorize_grlex");
    {
        polynomial_ZZ fg;
        fg = pow(x,2) - 1;
        auto fac = factorize(fg);
        CLPOLY_ASSERT(verify_factorization(fg, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)2);
    }

    return clpoly_test::test_summary();
}
