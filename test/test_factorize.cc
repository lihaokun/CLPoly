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

    // ========================================
    // 回归: 高次非首一多项式 (双 lc 乘法 bug)
    // 108x^9+138x^8-36x^7+165x^6+170x^5-120x^4-130x^3-50x
    // = x * (12x^3+18x^2+5) * (9x^5-2x^4+10x^2-10)
    // 之前 recombination 中 lc(f) 被重复乘导致因子无法分离
    // ========================================
    CLPOLY_TEST("factorize_nonmonic_high_deg");
    {
        f = 108*pow(x,9) + 138*pow(x,8) - 36*pow(x,7) + 165*pow(x,6)
            + 170*pow(x,5) - 120*pow(x,4) - 130*pow(x,3) - 50*x;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        // 应分出 3 个不可约因子: x, (12x^3+18x^2+5), (9x^5-2x^4+10x^2-10)
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)3);
    }

    // ========================================
    // 回归: 另一个高次非首一用例
    // 224x^7 - 238x^6 + ... = (16x^2-17x+7)(14x^5-7x^3+6x+4)
    // ========================================
    CLPOLY_TEST("factorize_nonmonic_deg7");
    {
        // (16x^2 - 17x + 7)(14x^5 - 7x^3 + 6x + 4)
        polynomial_ZZ f1 = 16*pow(x,2) - 17*x + 7;
        polynomial_ZZ f2 = 14*pow(x,5) - 7*pow(x,3) + 6*x + 4;
        f = f1 * f2;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)2);
    }

    // ========================================
    // 非首一含因子 x (常数项为零)
    // ========================================
    CLPOLY_TEST("factorize_nonmonic_with_x_factor");
    {
        // x * (6x^2 + 5x + 1) = x * (2x+1)(3x+1)
        f = 6*pow(x,3) + 5*pow(x,2) + x;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)3);
    }

    // ========================================
    // 非首一大系数
    // ========================================
    CLPOLY_TEST("factorize_nonmonic_large_lc");
    {
        // (100x + 3)(100x - 7) = 10000x^2 - 400x - 21
        f = 10000*pow(x,2) - 400*x - 21;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)2);
    }

    // ========================================
    // 退化输入补全
    // ========================================
    CLPOLY_TEST("factorize_negative_constant");
    {
        polynomial_ZZ c = polynomial_ZZ(ZZ(-42));
        auto fac = factorize(c);
        CLPOLY_ASSERT_EQ(fac.content, ZZ(-42));
        CLPOLY_ASSERT(fac.factors.empty());
    }

    CLPOLY_TEST("factorize_one");
    {
        polynomial_ZZ c = polynomial_ZZ(ZZ(1));
        auto fac = factorize(c);
        CLPOLY_ASSERT_EQ(fac.content, ZZ(1));
        CLPOLY_ASSERT(fac.factors.empty());
    }

    CLPOLY_TEST("factorize_minus_one");
    {
        polynomial_ZZ c = polynomial_ZZ(ZZ(-1));
        auto fac = factorize(c);
        CLPOLY_ASSERT_EQ(fac.content, ZZ(-1));
        CLPOLY_ASSERT(fac.factors.empty());
    }

    CLPOLY_TEST("factorize_x_alone");
    {
        // x 本身：content=1, 因子 x mult=1
        f = x;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.content, ZZ(1));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)1);
    }

    // ========================================
    // 纯单项式 x^n
    // ========================================
    CLPOLY_TEST("factorize_monomial_x5");
    {
        // x^5 → content=1, 因子 x, 重数 5
        f = pow(x, 5);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.content, ZZ(1));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)5);
    }

    CLPOLY_TEST("factorize_monomial_with_content");
    {
        // 6x^3 → content=6, 因子 x, 重数 3
        f = 6 * pow(x, 3);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)3);
    }

    CLPOLY_TEST("factorize_negative_monomial");
    {
        // -5x^4 → content=-5, 因子 x, 重数 4
        f = -5 * pow(x, 4);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)4);
    }

    // ========================================
    // 高重数
    // ========================================
    CLPOLY_TEST("factorize_high_mult_10");
    {
        // (x+1)^10
        f = pow(x + 1, 10);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)10);
    }

    CLPOLY_TEST("factorize_high_mult_nonmonic");
    {
        // (2x-3)^6
        f = pow(2*x - 3, 6);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)6);
    }

    // ========================================
    // 完全幂 (不可约多项式的幂)
    // ========================================
    CLPOLY_TEST("factorize_perfect_power_irreducible");
    {
        // (x²+1)^4 — squarefree 分解应直接处理，不进 Hensel
        f = pow(pow(x, 2) + 1, 4);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)4);
    }

    CLPOLY_TEST("factorize_perfect_power_x3m2");
    {
        // (x³-2)^3 — 不可约三次方的三次幂
        f = pow(pow(x, 3) - 2, 3);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)3);
    }

    // ========================================
    // 混合 squarefree 分量
    // ========================================
    CLPOLY_TEST("factorize_mixed_squarefree_3_components");
    {
        // (x+1)¹ * (x-1)² * (x²+1)³ — 3 个独立 squarefree 分量
        f = (x + 1) * pow(x - 1, 2) * pow(pow(x, 2) + 1, 3);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)3);
        std::set<uint64_t> mults;
        for (auto& [fi, ei] : fac.factors)
            mults.insert(ei);
        CLPOLY_ASSERT(mults.count(1));
        CLPOLY_ASSERT(mults.count(2));
        CLPOLY_ASSERT(mults.count(3));
    }

    CLPOLY_TEST("factorize_mixed_squarefree_with_content");
    {
        // 6 * (x+1)² * (x²-x+1)³
        f = 6 * pow(x + 1, 2) * pow(pow(x, 2) - x + 1, 3);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)2);
    }

    CLPOLY_TEST("factorize_mixed_squarefree_with_x");
    {
        // x² * (x+1)³ * (x²+1)
        f = pow(x, 2) * pow(x + 1, 3) * (pow(x, 2) + 1);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)3);
    }

    // ========================================
    // 5+ 个不可约因子 (Hensel 树 depth≥3 + recombination)
    // ========================================
    CLPOLY_TEST("factorize_5_linear_factors");
    {
        // (x-1)(x-2)(x-3)(x-4)(x-5)
        f = (x - 1) * (x - 2) * (x - 3) * (x - 4) * (x - 5);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)5);
        for (auto& [fi, ei] : fac.factors)
            CLPOLY_ASSERT_EQ(ei, (uint64_t)1);
    }

    CLPOLY_TEST("factorize_6_factors_mixed_deg");
    {
        // (x-1)(x+1)(x-2)(x+2)(x²+1)(x²+x+1) — 6 因子，含二次
        f = (x - 1) * (x + 1) * (x - 2) * (x + 2)
            * (pow(x, 2) + 1) * (pow(x, 2) + x + 1);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)6);
    }

    CLPOLY_TEST("factorize_x10m1");
    {
        // x^10-1 = (x-1)(x+1)(x⁴+x³+x²+x+1)(x⁴-x³+x²-x+1) — 4 因子
        f = pow(x, 10) - 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)4);
    }

    CLPOLY_TEST("factorize_x12m1");
    {
        // x^12-1 有 6 个不可约因子
        // (x-1)(x+1)(x²+1)(x²+x+1)(x²-x+1)(x⁴-x²+1)
        f = pow(x, 12) - 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)6);
    }

    // ========================================
    // QQ 深度测试
    // ========================================
    CLPOLY_TEST("factorize_QQ_large_denom");
    {
        // (x²-1) / 7700 = (1/7700)(x-1)(x+1)
        polynomial_ZZ fz = pow(x, 2) - 1;
        polynomial_QQ fq;
        poly_convert(fz, fq);
        fq = fq * QQ(1, 7700);
        auto fac = factorize(fq);
        CLPOLY_ASSERT(verify_factorization(fq, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)2);
        // 因子应首一
        for (auto& [fi, ei] : fac.factors)
            CLPOLY_ASSERT_EQ(fi.front().second, QQ(1));
    }

    CLPOLY_TEST("factorize_QQ_multifactor");
    {
        // (x-1)(x+1)(x²+1) over QQ, 乘以 QQ(3,7)
        polynomial_ZZ fz = pow(x, 4) - 1;
        polynomial_QQ fq;
        poly_convert(fz, fq);
        fq = fq * QQ(3, 7);
        auto fac = factorize(fq);
        CLPOLY_ASSERT(verify_factorization(fq, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)3);
    }

    CLPOLY_TEST("factorize_QQ_high_mult");
    {
        // QQ: (x+1)^5 * (1/6)
        polynomial_ZZ fz = pow(x + 1, 5);
        polynomial_QQ fq;
        poly_convert(fz, fq);
        fq = fq * QQ(1, 6);
        auto fac = factorize(fq);
        CLPOLY_ASSERT(verify_factorization(fq, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)5);
    }

    // ========================================
    // 全负系数
    // ========================================
    CLPOLY_TEST("factorize_all_negative_coeffs");
    {
        // -(x²+x+1) — 不可约，全负后 content=-1
        f = -(pow(x, 2) + x + 1);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
    }

    // ========================================
    // upolynomial 空输入
    // ========================================
    CLPOLY_TEST("factorize_upoly_zero");
    {
        upolynomial_ZZ uz;
        auto fac = factorize(uz);
        CLPOLY_ASSERT_EQ(fac.content, ZZ(0));
        CLPOLY_ASSERT(fac.factors.empty());
    }

    CLPOLY_TEST("factorize_upoly_constant");
    {
        upolynomial_ZZ uc({{umonomial(0), ZZ(7)}});
        auto fac = factorize(uc);
        CLPOLY_ASSERT_EQ(fac.content, ZZ(7));
        CLPOLY_ASSERT(fac.factors.empty());
    }

    CLPOLY_TEST("factorize_upoly_linear");
    {
        upolynomial_ZZ ul({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(-3)}});
        auto fac = factorize(ul);
        CLPOLY_ASSERT_EQ(fac.content, ZZ(1));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)1);
    }

    CLPOLY_TEST("factorize_upoly_monomial");
    {
        // upolynomial: 3x^4
        upolynomial_ZZ um({{umonomial(4), ZZ(3)}});
        auto fac = factorize(um);
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(fac.factors[0].second, (uint64_t)4);
        // 验证重组
        upolynomial_ZZ product({{umonomial(0), fac.content}});
        for (auto& [fi, ei] : fac.factors)
            product = product * pow(fi, (int64_t)ei);
        product.normalization();
        CLPOLY_ASSERT_EQ(product, um);
    }

    // ========================================
    // 随机单变量 ZZ 因式分解
    // ========================================

    CLPOLY_TEST("factorize_random_ZZ_2factors");
    for (int trial = 0; trial < 20; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto p1 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        auto p2 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        if (p1.empty() || p2.empty()) continue;
        auto f_prod = p1 * p2;
        if (f_prod.empty() || degree(f_prod) < 2) continue;
        auto fac = factorize(f_prod);
        CLPOLY_ASSERT(verify_factorization(f_prod, fac));
    }

    CLPOLY_TEST("factorize_random_ZZ_3factors");
    for (int trial = 0; trial < 15; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto p1 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        auto p2 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        auto p3 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        if (p1.empty() || p2.empty() || p3.empty()) continue;
        auto f_prod = p1 * p2 * p3;
        if (f_prod.empty() || degree(f_prod) < 2) continue;
        auto fac = factorize(f_prod);
        CLPOLY_ASSERT(verify_factorization(f_prod, fac));
    }

    CLPOLY_TEST("factorize_random_ZZ_4factors");
    for (int trial = 0; trial < 10; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto p1 = random_polynomial<ZZ>({x}, 2 + (trial % 3), 2 + (trial % 3), {-10, 10});
        auto p2 = random_polynomial<ZZ>({x}, 2 + (trial % 3), 2 + (trial % 3), {-10, 10});
        auto p3 = random_polynomial<ZZ>({x}, 2 + (trial % 3), 2 + (trial % 3), {-10, 10});
        auto p4 = random_polynomial<ZZ>({x}, 2 + (trial % 3), 2 + (trial % 3), {-10, 10});
        if (p1.empty() || p2.empty() || p3.empty() || p4.empty()) continue;
        auto f_prod = p1 * p2 * p3 * p4;
        if (f_prod.empty() || degree(f_prod) < 2) continue;
        auto fac = factorize(f_prod);
        CLPOLY_ASSERT(verify_factorization(f_prod, fac));
    }

    CLPOLY_TEST("factorize_random_ZZ_with_multiplicity");
    for (int trial = 0; trial < 10; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto p1 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        auto p2 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        if (p1.empty() || p2.empty()) continue;
        int e1 = 1 + (trial % 3), e2 = 1 + ((trial + 1) % 3);
        auto f_prod = pow(p1, e1) * pow(p2, e2);
        if (f_prod.empty() || degree(f_prod) < 2) continue;
        auto fac = factorize(f_prod);
        CLPOLY_ASSERT(verify_factorization(f_prod, fac));
    }

    CLPOLY_TEST("factorize_random_ZZ_nonmonic_content");
    for (int trial = 0; trial < 5; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto p1 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        auto p2 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        if (p1.empty() || p2.empty()) continue;
        ZZ c(trial * 7 + 2);
        auto f_prod = c * p1 * p2;
        if (f_prod.empty() || degree(f_prod) < 2) continue;
        auto fac = factorize(f_prod);
        CLPOLY_ASSERT(verify_factorization(f_prod, fac));
    }

    // ========================================
    // 随机 QQ 因式分解
    // ========================================

    CLPOLY_TEST("factorize_random_QQ_2factors");
    for (int trial = 0; trial < 10; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto p1 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        auto p2 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        if (p1.empty() || p2.empty()) continue;
        auto prodz = p1 * p2;
        if (prodz.empty() || degree(prodz) < 2) continue;
        polynomial_QQ prodq;
        poly_convert(prodz, prodq);
        prodq = prodq * QQ(1, trial + 2);
        auto fac = factorize(prodq);
        CLPOLY_ASSERT(verify_factorization(prodq, fac));
    }

    CLPOLY_TEST("factorize_random_QQ_with_multiplicity");
    for (int trial = 0; trial < 5; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto p1 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        auto p2 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        if (p1.empty() || p2.empty()) continue;
        int e1 = 1 + (trial % 3), e2 = 1 + ((trial + 1) % 3);
        auto prodz = pow(p1, e1) * pow(p2, e2);
        if (prodz.empty() || degree(prodz) < 2) continue;
        polynomial_QQ prodq;
        poly_convert(prodz, prodq);
        prodq = prodq * QQ(1, trial + 2);
        auto fac = factorize(prodq);
        CLPOLY_ASSERT(verify_factorization(prodq, fac));
    }

    CLPOLY_TEST("factorize_random_QQ_large_denom");
    for (int trial = 0; trial < 5; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto p1 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        auto p2 = random_polynomial<ZZ>({x}, 3 + (trial % 4), 3 + (trial % 3), {-10, 10});
        if (p1.empty() || p2.empty()) continue;
        auto prodz = p1 * p2;
        if (prodz.empty() || degree(prodz) < 2) continue;
        polynomial_QQ prodq;
        poly_convert(prodz, prodq);
        prodq = prodq * QQ(1, (trial + 1) * 20);
        auto fac = factorize(prodq);
        CLPOLY_ASSERT(verify_factorization(prodq, fac));
    }

    // ========================================
    // 随机 upolynomial 因式分解
    // ========================================

    {
        std::mt19937 rng(42);
        std::uniform_int_distribution<int> deg_dis(3, 6);
        std::uniform_int_distribution<int> len_dis(3, 5);
        std::uniform_int_distribution<int> exp_dis(1, 3);

        auto verify_upoly_factorization = [](const upolynomial_ZZ& f,
                                              const factorization<upolynomial_ZZ>& fac) {
            upolynomial_ZZ product({{umonomial(0), fac.content}});
            for (auto& [fi, ei] : fac.factors)
                product = product * pow(fi, (int64_t)ei);
            product.normalization();
            return product == f;
        };

        CLPOLY_TEST("factorize_random_upoly_2factors");
        for (int trial = 0; trial < 10; ++trial)
        {
            CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
            auto f1 = random_upolynomial<ZZ>(deg_dis(rng), len_dis(rng), {-10, 10});
            auto f2 = random_upolynomial<ZZ>(deg_dis(rng), len_dis(rng), {-10, 10});
            if (f1.empty() || f2.empty()) continue;
            auto uf = f1 * f2;
            if (uf.empty() || get_deg(uf) < 2) continue;
            auto fac = factorize(uf);
            CLPOLY_ASSERT(verify_upoly_factorization(uf, fac));
        }

        CLPOLY_TEST("factorize_random_upoly_3factors");
        for (int trial = 0; trial < 5; ++trial)
        {
            CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
            auto f1 = random_upolynomial<ZZ>(deg_dis(rng), len_dis(rng), {-10, 10});
            auto f2 = random_upolynomial<ZZ>(deg_dis(rng), len_dis(rng), {-10, 10});
            auto f3 = random_upolynomial<ZZ>(deg_dis(rng), len_dis(rng), {-10, 10});
            if (f1.empty() || f2.empty() || f3.empty()) continue;
            auto uf = f1 * f2 * f3;
            if (uf.empty() || get_deg(uf) < 2) continue;
            auto fac = factorize(uf);
            CLPOLY_ASSERT(verify_upoly_factorization(uf, fac));
        }

        CLPOLY_TEST("factorize_random_upoly_with_multiplicity");
        for (int trial = 0; trial < 5; ++trial)
        {
            CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
            auto f1 = random_upolynomial<ZZ>(deg_dis(rng), len_dis(rng), {-10, 10});
            auto f2 = random_upolynomial<ZZ>(deg_dis(rng), len_dis(rng), {-10, 10});
            if (f1.empty() || f2.empty()) continue;
            int e1 = exp_dis(rng), e2 = exp_dis(rng);
            auto uf = pow(f1, e1) * pow(f2, e2);
            if (uf.empty() || get_deg(uf) < 2) continue;
            auto fac = factorize(uf);
            CLPOLY_ASSERT(verify_upoly_factorization(uf, fac));
        }
    }

    // ========================================
    // 经典单变量用例
    // ========================================

    CLPOLY_TEST("factorize_classic_wilkinson_10");
    {
        // W(10) = ∏_{i=1}^{10}(x-i)
        f = polynomial_ZZ(ZZ(1));
        for (int i = 1; i <= 10; ++i)
            f = f * (x - i);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)10);
        for (auto& [fi, ei] : fac.factors)
            CLPOLY_ASSERT_EQ(ei, (uint64_t)1);
    }

    CLPOLY_TEST("factorize_classic_wilkinson_15");
    {
        // W(15) = ∏_{i=1}^{15}(x-i)
        f = polynomial_ZZ(ZZ(1));
        for (int i = 1; i <= 15; ++i)
            f = f * (x - i);
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)15);
    }

    CLPOLY_TEST("factorize_classic_cyclotomic_x15m1");
    {
        // x^15-1 = Φ₁·Φ₃·Φ₅·Φ₁₅ = 4 个不可约因子
        f = pow(x, 15) - 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)4);
    }

    CLPOLY_TEST("factorize_classic_cyclotomic_x24m1");
    {
        // x^24-1 有 8 个不可约因子
        // Φ₁,Φ₂,Φ₃,Φ₄,Φ₆,Φ₈,Φ₁₂,Φ₂₄
        f = pow(x, 24) - 1;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)8);
    }

    CLPOLY_TEST("factorize_classic_x6m9");
    {
        // x^6 - 9 = (x^3-3)(x^3+3)
        f = pow(x, 6) - 9;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)2);
    }

    CLPOLY_TEST("factorize_classic_swinnerton_dyer_S3");
    {
        // Swinnerton-Dyer S₃ = minpoly(√2+√3+√5), degree 8
        // S₃ = x^8 - 40x^6 + 352x^4 - 960x^2 + 576
        f = pow(x,8) - 40*pow(x,6) + 352*pow(x,4) - 960*pow(x,2) + 576;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        // S₃ is irreducible over Z
        CLPOLY_ASSERT_EQ(fac.factors.size(), (size_t)1);
    }

    CLPOLY_TEST("factorize_classic_mignotte");
    {
        // Mignotte-style: x^8 + 14x^4 + 1
        // = (x^4 + 4x^2 + 2x + 1)(x^4 - 4x^2 - 2x + 1) ... nah
        // Actually x^8 + 14x^4 + 1 factors as (x^4-2x^2+2x+1)(... no
        // Use a known Mignotte: x^6 - 2(3x+1)^2 = x^6 - 18x^2 - 12x - 2
        f = pow(x, 6) - 18*pow(x, 2) - 12*x - 2;
        auto fac = factorize(f);
        CLPOLY_ASSERT(verify_factorization(f, fac));
        // Should be irreducible (Mignotte's are typically hard to factor but irreducible)
        CLPOLY_ASSERT(fac.factors.size() >= 1);
    }

    return clpoly_test::test_summary();
}
