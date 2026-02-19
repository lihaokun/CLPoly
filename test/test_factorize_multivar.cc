/**
 * @file test_factorize_multivar.cc
 * @brief Phase 4 端到端测试: 多变量 factorize()
 */
#include "clpoly_test.hh"
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <iostream>

using namespace clpoly;
using PolyZZ = polynomial_<ZZ, lex>;

PolyZZ make_lex(const polynomial_ZZ& p)
{
    PolyZZ result;
    poly_convert(p, result);
    return result;
}

// 验证: content * ∏ fᵢ^eᵢ == f
void verify_factorization(const PolyZZ& f, const factorization<PolyZZ>& fac, const char* name)
{
    PolyZZ prod;
    if (fac.content != ZZ(0))
    {
        basic_monomial<lex> m0;
        prod.push_back(std::make_pair(m0, fac.content));
    }

    for (auto& [fi, ei] : fac.factors)
    {
        PolyZZ fi_pow = fi;
        for (uint64_t e = 1; e < ei; ++e)
        {
            fi_pow = fi_pow * fi;
            fi_pow.normalization();
        }
        prod = prod * fi_pow;
        prod.normalization();
    }

    CLPOLY_ASSERT_EQ(prod, f);
}

int main()
{
    variable x("x"), y("y"), z("z"), w("w");

    // ============================================================
    // 基本二变量
    // ============================================================

    CLPOLY_TEST("factorize: x^2 - y^2 = (x+y)(x-y)");
    {
        auto f = make_lex(pow(x,2) - pow(y,2));
        auto fac = factorize(f);
        verify_factorization(f, fac, "x^2-y^2");
        CLPOLY_ASSERT(fac.factors.size() == 2);
    }

    CLPOLY_TEST("factorize: x^2 + 2xy + y^2 = (x+y)^2");
    {
        auto f = make_lex(pow(x,2) + 2*x*y + pow(y,2));
        auto fac = factorize(f);
        verify_factorization(f, fac, "x^2+2xy+y^2");
        // (x+y)^2: 一个因子, 重数 2
        CLPOLY_ASSERT(fac.factors.size() == 1);
        CLPOLY_ASSERT(fac.factors[0].second == 2);
    }

    CLPOLY_TEST("factorize: (x+y)(x-y+1)");
    {
        auto f = make_lex(pow(x,2) + x - pow(y,2) + y);
        auto fac = factorize(f);
        verify_factorization(f, fac, "(x+y)(x-y+1)");
        CLPOLY_ASSERT(fac.factors.size() == 2);
    }

    CLPOLY_TEST("factorize: x^3 + y^3 = (x+y)(x^2-xy+y^2)");
    {
        auto f = make_lex(pow(x,3) + pow(y,3));
        auto fac = factorize(f);
        verify_factorization(f, fac, "x^3+y^3");
        CLPOLY_ASSERT(fac.factors.size() == 2);
    }

    // ============================================================
    // 非平凡 lc
    // ============================================================

    CLPOLY_TEST("factorize: nontrivial lc, ((y+1)x+y)(x+1)");
    {
        // f = (y+1)x^2 + (2y+1)x + y
        auto f = make_lex((y + polynomial_ZZ(1))*pow(x,2) + (2*y + polynomial_ZZ(1))*x + y);
        auto fac = factorize(f);
        verify_factorization(f, fac, "nontrivial lc");
        CLPOLY_ASSERT(fac.factors.size() == 2);
    }

    CLPOLY_TEST("factorize: ((y-1)x+1)((y+1)x+1)");
    {
        // f = (y^2-1)x^2 + 2yx + 1
        auto f = make_lex((pow(y,2) - polynomial_ZZ(1))*pow(x,2) + 2*y*x + polynomial_ZZ(1));
        auto fac = factorize(f);
        verify_factorization(f, fac, "two lc factors");
        CLPOLY_ASSERT(fac.factors.size() == 2);
    }

    // ============================================================
    // 三因子
    // ============================================================

    CLPOLY_TEST("factorize: (x+y)(x-y)(x+1)");
    {
        auto f = make_lex(pow(x,3) + pow(x,2) - x*pow(y,2) - pow(y,2));
        auto fac = factorize(f);
        verify_factorization(f, fac, "three factors");
        CLPOLY_ASSERT(fac.factors.size() == 3);
    }

    // ============================================================
    // 三变量
    // ============================================================

    CLPOLY_TEST("factorize: (x+y+z)(x-y+z)");
    {
        auto f = make_lex(pow(x,2) + 2*x*z + pow(z,2) - pow(y,2));
        auto fac = factorize(f);
        verify_factorization(f, fac, "trivariate");
        CLPOLY_ASSERT(fac.factors.size() == 2);
    }

    CLPOLY_TEST("factorize: (x+y)(x+z)");
    {
        auto f = make_lex(pow(x,2) + x*z + x*y + y*z);
        auto fac = factorize(f);
        verify_factorization(f, fac, "(x+y)(x+z)");
        CLPOLY_ASSERT(fac.factors.size() == 2);
    }

    // ============================================================
    // 含 content
    // ============================================================

    CLPOLY_TEST("factorize: 2*(x^2 - y^2)");
    {
        auto f = make_lex(2*pow(x,2) - 2*pow(y,2));
        auto fac = factorize(f);
        verify_factorization(f, fac, "2*(x^2-y^2)");
        CLPOLY_ASSERT(fac.content == ZZ(2));
        CLPOLY_ASSERT(fac.factors.size() == 2);
    }

    CLPOLY_TEST("factorize: (y+1)*(x^2-y^2)");
    {
        // f = (y+1)(x+y)(x-y) = (y+1)(x^2-y^2)
        auto f = make_lex((y + polynomial_ZZ(1))*(pow(x,2) - pow(y,2)));
        auto fac = factorize(f);
        verify_factorization(f, fac, "(y+1)(x^2-y^2)");
        CLPOLY_ASSERT(fac.factors.size() == 3);
    }

    // ============================================================
    // 含重因子
    // ============================================================

    CLPOLY_TEST("factorize: (x+y)^2*(x-y)");
    {
        auto f = make_lex(pow(x,3) + pow(x,2)*y - x*pow(y,2) - pow(y,3));
        auto fac = factorize(f);
        verify_factorization(f, fac, "(x+y)^2*(x-y)");
        // 应有两个不同因子, 其中一个重数 2
        bool has_mult_2 = false;
        for (auto& [fi, ei] : fac.factors)
            if (ei == 2) has_mult_2 = true;
        CLPOLY_ASSERT(has_mult_2);
    }

    // ============================================================
    // 不可约
    // ============================================================

    CLPOLY_TEST("factorize: irreducible x^2 + y^2 + 1");
    {
        auto f = make_lex(pow(x,2) + pow(y,2) + polynomial_ZZ(1));
        auto fac = factorize(f);
        verify_factorization(f, fac, "irreducible");
        CLPOLY_ASSERT(fac.factors.size() == 1);
        CLPOLY_ASSERT(fac.factors[0].second == 1);
    }

    // ============================================================
    // 负首项系数
    // ============================================================

    CLPOLY_TEST("factorize: -(x^2 - y^2)");
    {
        auto f = make_lex(-pow(x,2) + pow(y,2));
        auto fac = factorize(f);
        verify_factorization(f, fac, "negative lc");
        CLPOLY_ASSERT(fac.factors.size() == 2);
    }

    // ============================================================
    // grlex 输入 (通过通用 comp dispatch)
    // ============================================================

    CLPOLY_TEST("factorize: grlex input x^2 - y^2");
    {
        polynomial_ZZ f = pow(x,2) - pow(y,2);
        auto fac = factorize(f);
        // 验证 content * ∏ fi^ei == f
        polynomial_ZZ prod;
        basic_monomial<grlex> m0;
        prod.push_back({m0, fac.content});
        for (auto& [fi, ei] : fac.factors)
        {
            polynomial_ZZ fi_pow = fi;
            for (uint64_t e = 1; e < ei; ++e)
            {
                fi_pow = fi_pow * fi;
                fi_pow.normalization();
            }
            prod = prod * fi_pow;
            prod.normalization();
        }
        CLPOLY_ASSERT_EQ(prod, f);
        CLPOLY_ASSERT(fac.factors.size() == 2);
    }

    // ============================================================
    // QQ 输入
    // ============================================================

    CLPOLY_TEST("factorize: QQ bivariate");
    {
        variable xq("x"), yq("y");
        polynomial_QQ f = pow(xq,2) - pow(yq,2);
        auto fac = factorize(f);
        // 验证
        polynomial_QQ prod;
        basic_monomial<grlex> m0;
        prod.push_back({m0, fac.content});
        for (auto& [fi, ei] : fac.factors)
        {
            polynomial_QQ fi_pow = fi;
            for (uint64_t e = 1; e < ei; ++e)
            {
                fi_pow = fi_pow * fi;
                fi_pow.normalization();
            }
            prod = prod * fi_pow;
            prod.normalization();
        }
        CLPOLY_ASSERT_EQ(prod, f);
        CLPOLY_ASSERT(fac.factors.size() == 2);
    }

    // ============================================================
    // D5 回归测试: 曾因 LC 分配 + 求值点确定性 bug 失败的用例
    // ============================================================

    CLPOLY_TEST("factorize: D5 bivar regression 1");
    {
        // f1 = 4x³ + y³ + 3xy - 3y², f2 = 5x³ + 5xy - 2x - 2y
        auto f1 = make_lex(4*pow(x,3) + pow(y,3) + 3*x*y - 3*pow(y,2));
        auto f2 = make_lex(5*pow(x,3) + 5*x*y - 2*x - 2*y);
        auto f = f1 * f2; f.normalization();
        auto fac = factorize(f);
        verify_factorization(f, fac, "D5 bivar 1");
        CLPOLY_ASSERT(fac.factors.size() >= 2);
    }

    CLPOLY_TEST("factorize: D5 bivar regression 2");
    {
        // f1 = 2x²y + 2xy² - 5x - 4y, f2 = -4x³ + 5y³ - 4x² + 3y
        auto f1 = make_lex(2*pow(x,2)*y + 2*x*pow(y,2) - 5*x - 4*y);
        auto f2 = make_lex(-4*pow(x,3) + 5*pow(y,3) - 4*pow(x,2) + 3*y);
        auto f = f1 * f2; f.normalization();
        auto fac = factorize(f);
        verify_factorization(f, fac, "D5 bivar 2");
        CLPOLY_ASSERT(fac.factors.size() >= 2);
    }

    CLPOLY_TEST("factorize: D5 trivar regression 1");
    {
        // f1 = 2xy - 2xz - 3y² - 3x, f2 = -x² + y - z - 2
        auto f1 = make_lex(2*x*y - 2*x*z - 3*pow(y,2) - 3*x);
        auto f2 = make_lex(-pow(x,2) + y - z - 2*polynomial_ZZ(1));
        auto f = f1 * f2; f.normalization();
        auto fac = factorize(f);
        verify_factorization(f, fac, "D5 trivar 1");
        CLPOLY_ASSERT(fac.factors.size() >= 2);
    }

    CLPOLY_TEST("factorize: D5 trivar regression 2");
    {
        // f1 = x² - xz - 3y² - x, f2 = xz - 3z² - x - z
        auto f1 = make_lex(pow(x,2) - x*z - 3*pow(y,2) - x);
        auto f2 = make_lex(x*z - 3*pow(z,2) - x - z);
        auto f = f1 * f2; f.normalization();
        auto fac = factorize(f);
        verify_factorization(f, fac, "D5 trivar 2");
        CLPOLY_ASSERT(fac.factors.size() >= 2);
    }

    CLPOLY_TEST("factorize: D5 trivar regression 3");
    {
        // f1 = -2xy - 2y² + 3yz - 2y, f2 = -2xz + 2yz + 3z² - 2y
        auto f1 = make_lex(-2*x*y - 2*pow(y,2) + 3*y*z - 2*y);
        auto f2 = make_lex(-2*x*z + 2*y*z + 3*pow(z,2) - 2*y);
        auto f = f1 * f2; f.normalization();
        auto fac = factorize(f);
        verify_factorization(f, fac, "D5 trivar 3");
        CLPOLY_ASSERT(fac.factors.size() >= 2);
    }

    // ============================================================
    // 单项式 content 提取回归测试
    // ============================================================

    CLPOLY_TEST("factorize: monomial content x divides all terms");
    {
        // g = x · (-2xy+2yz-2y+3z) · (2x-2y-3z+1), 每项含 x
        auto f1 = polynomial_ZZ(-ZZ(2)*x*y + ZZ(2)*y*z - ZZ(2)*y + ZZ(3)*z);
        auto f2 = polynomial_ZZ(ZZ(2)*x*x - ZZ(2)*x*y - ZZ(3)*x*z + x);
        auto f = make_lex(f1 * f2);
        auto fac = factorize(f);
        verify_factorization(f, fac, "mono content x");
        CLPOLY_ASSERT(fac.factors.size() >= 3);
    }

    CLPOLY_TEST("factorize: monomial content xy divides all terms");
    {
        // g = xy · (x+y+1) · (x-y+2)
        auto f = make_lex(polynomial_ZZ(x*y) *
                          polynomial_ZZ(x+y+ZZ(1)) *
                          polynomial_ZZ(x-y+ZZ(2)));
        auto fac = factorize(f);
        verify_factorization(f, fac, "mono content xy");
        CLPOLY_ASSERT(fac.factors.size() >= 4);
    }

    CLPOLY_TEST("factorize: no monomial content (has constant term)");
    {
        // g = (x+y+1)(x-y+2), 有常数项，不提取
        auto f = make_lex(polynomial_ZZ(x+y+ZZ(1)) *
                          polynomial_ZZ(x-y+ZZ(2)));
        auto fac = factorize(f);
        verify_factorization(f, fac, "no mono content");
        CLPOLY_ASSERT(fac.factors.size() >= 2);
    }

    // ============================================================
    // gamma 互素性回归测试 (D6)
    // ============================================================

    CLPOLY_TEST("factorize: gamma coprimality fail1 (L=6z)");
    {
        // gamma=6=2·3, lj=z → 需要 |z(α)| 与 6 互素
        auto f = make_lex(
            polynomial_ZZ(ZZ(3)*x*z + ZZ(2)*z*z - z + ZZ(2)) *
            polynomial_ZZ(ZZ(2)*x*x - ZZ(2)*x*y - ZZ(3)*y + ZZ(1)));
        auto fac = factorize(f);
        verify_factorization(f, fac, "gamma coprime fail1");
        CLPOLY_ASSERT(fac.factors.size() >= 2);
    }

    CLPOLY_TEST("factorize: gamma coprimality fail2 (L=6y)");
    {
        // gamma=6=2·3, lj=y → 需要 |y(α)| 与 6 互素
        auto f = make_lex(
            polynomial_ZZ(ZZ(3)*x*y - y*y - ZZ(3)*y*z + ZZ(2)) *
            polynomial_ZZ(ZZ(2)*x*x + ZZ(2)*z*z + x + ZZ(3)));
        auto fac = factorize(f);
        verify_factorization(f, fac, "gamma coprime fail2");
        CLPOLY_ASSERT(fac.factors.size() >= 2);
    }

    CLPOLY_TEST("factorize: gamma shares prime with lc factor (L=-4y²+4y)");
    {
        // gamma=-4, lc_factors=[(y,1),(y-1,1)]
        // For any integer y, one of {y, y-1} is even → gcd(4, ·) > 1
        // Power extraction can't handle this; GCD matching can
        auto f = make_lex(
            polynomial_ZZ(ZZ(2)*x*y - ZZ(2)*x + ZZ(3)) *
            polynomial_ZZ(x*y + ZZ(3)*y - ZZ(3)) *
            polynomial_ZZ(ZZ(2)*x*x - ZZ(2)*x*y + y*y));
        auto fac = factorize(f);
        verify_factorization(f, fac, "gamma shares prime");
        CLPOLY_ASSERT(fac.factors.size() >= 3);
    }

    // ============================================================
    // 随机二变量因式分解
    // ============================================================

    CLPOLY_TEST("factorize: random bivariate 2 factors");
    for (int trial = 0; trial < 10; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1_gr = random_polynomial<ZZ>({x, y}, 3, 4, {-5, 5});
        auto f2_gr = random_polynomial<ZZ>({x, y}, 3, 4, {-5, 5});
        if (f1_gr.empty() || f2_gr.empty()) continue;
        auto f_gr = f1_gr * f2_gr;
        if (f_gr.empty()) continue;
        auto f = make_lex(f_gr);
        auto fac = factorize(f);
        verify_factorization(f, fac, "random bivar");
    }

    CLPOLY_TEST("factorize: random bivariate 3 factors");
    for (int trial = 0; trial < 5; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1_gr = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
        auto f2_gr = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
        auto f3_gr = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
        if (f1_gr.empty() || f2_gr.empty() || f3_gr.empty()) continue;
        auto f_gr = f1_gr * f2_gr * f3_gr;
        if (f_gr.empty()) continue;
        auto f = make_lex(f_gr);
        auto fac = factorize(f);
        verify_factorization(f, fac, "random bivar 3fac");
    }

    // ============================================================
    // 随机三变量因式分解
    // ============================================================

    CLPOLY_TEST("factorize: random trivariate 2 factors");
    for (int trial = 0; trial < 5; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1_gr = random_polynomial<ZZ>({x, y, z}, 2, 4, {-3, 3});
        auto f2_gr = random_polynomial<ZZ>({x, y, z}, 2, 4, {-3, 3});
        if (f1_gr.empty() || f2_gr.empty()) continue;
        auto f_gr = f1_gr * f2_gr;
        if (f_gr.empty()) continue;
        auto f = make_lex(f_gr);
        auto fac = factorize(f);
        verify_factorization(f, fac, "random trivar");
    }

    // ============================================================
    // grlex dispatch: 随机
    // ============================================================

    CLPOLY_TEST("factorize: random grlex bivariate");
    for (int trial = 0; trial < 5; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_polynomial<ZZ>({x, y}, 3, 4, {-5, 5});
        auto f2 = random_polynomial<ZZ>({x, y}, 3, 4, {-5, 5});
        if (f1.empty() || f2.empty()) continue;
        polynomial_ZZ f = f1 * f2;
        if (f.empty()) continue;
        auto fac = factorize(f);
        // Verify product
        polynomial_ZZ prod;
        basic_monomial<grlex> m0;
        prod.push_back({m0, fac.content});
        for (auto& [fi, ei] : fac.factors) {
            polynomial_ZZ fi_pow = fi;
            for (uint64_t e = 1; e < ei; ++e) {
                fi_pow = fi_pow * fi;
                fi_pow.normalization();
            }
            prod = prod * fi_pow;
            prod.normalization();
        }
        CLPOLY_ASSERT_EQ(prod, f);
    }

    // ============================================================
    // 带 content 的随机因式分解
    // ============================================================

    CLPOLY_TEST("factorize: random bivariate with content");
    for (int trial = 0; trial < 5; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1_gr = random_polynomial<ZZ>({x, y}, 3, 4, {-5, 5});
        auto f2_gr = random_polynomial<ZZ>({x, y}, 3, 4, {-5, 5});
        if (f1_gr.empty() || f2_gr.empty()) continue;
        ZZ c(trial + 2);
        auto f_gr = c * f1_gr * f2_gr;
        if (f_gr.empty()) continue;
        auto f = make_lex(f_gr);
        auto fac = factorize(f);
        verify_factorization(f, fac, "random bivar content");
    }

    // ============================================================
    // 带重数的随机因式分解
    // ============================================================

    CLPOLY_TEST("factorize: random bivariate with multiplicity");
    for (int trial = 0; trial < 5; ++trial)
    {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1_gr = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
        auto f2_gr = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
        if (f1_gr.empty() || f2_gr.empty()) continue;
        auto f_gr = pow(f1_gr, 2) * f2_gr;
        if (f_gr.empty()) continue;
        auto f = make_lex(f_gr);
        auto fac = factorize(f);
        verify_factorization(f, fac, "random bivar mult");
    }

    return clpoly_test::test_summary();
}
