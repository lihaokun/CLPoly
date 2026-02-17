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

    return clpoly_test::test_summary();
}
