/**
 * @file test_multivar_helpers.cc
 * @brief Phase 1 测试: pp(), polynomial_GCD(ZZ), __taylor_coeff()
 */
#include "clpoly_test.hh"
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <iostream>

using namespace clpoly;
using PolyZZ = polynomial_<ZZ, lex>;

// Helper: 构造 lex 多项式
PolyZZ make_lex(const polynomial_ZZ& p)
{
    PolyZZ result;
    poly_convert(p, result);
    return result;
}

int main()
{
    variable x("x"), y("y"), z("z");

    // ============================================================
    CLPOLY_TEST("pp: bivariate");
    {
        // f = 2x^2*y + 4x*y + 2y = 2y * (x^2 + 2x + 1)
        auto f = make_lex(2*pow(x,2)*y + 4*x*y + 2*y);
        auto c = cont(f);
        auto p = pp(f);
        CLPOLY_ASSERT_EQ(c * p, f);
        auto expected_cont = make_lex(2*y);
        auto expected_pp = make_lex(pow(x,2) + 2*x + polynomial_ZZ(1));
        CLPOLY_ASSERT_EQ(c, expected_cont);
        CLPOLY_ASSERT_EQ(p, expected_pp);
    }

    CLPOLY_TEST("pp: trivariate");
    {
        // f = (y+z)(x^2 + x), cont(f, x) = y+z, pp(f) = x^2+x
        auto f = make_lex(pow(x,2)*y + pow(x,2)*z + x*y + x*z);
        auto c = cont(f);
        auto p = pp(f);
        CLPOLY_ASSERT_EQ(c * p, f);
    }

    CLPOLY_TEST("pp: already primitive");
    {
        auto f = make_lex(pow(x,2) + y);
        auto c = cont(f);
        auto p = pp(f);
        CLPOLY_ASSERT_EQ(c * p, f);
        CLPOLY_ASSERT(is_number(c));
    }

    CLPOLY_TEST("pp: constant");
    {
        auto f = make_lex(polynomial_ZZ(6));
        auto p = pp(f);
        CLPOLY_ASSERT(is_number(p));
        CLPOLY_ASSERT_EQ(p.front().second, ZZ(1));
    }

    // ============================================================
    CLPOLY_TEST("polynomial_GCD(ZZ): basic coprime");
    {
        // a = x^2 + x + 1, b = x - 1 (coprime)
        upolynomial_<ZZ> a({{umonomial(2), ZZ(1)}, {umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});
        upolynomial_<ZZ> b({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(-1)}});
        upolynomial_<ZZ> s, t;
        auto c = polynomial_GCD(a, b, s, t);
        // s*a + t*b = c
        auto lhs = s * a + t * b;
        lhs.normalization();
        CLPOLY_ASSERT_EQ(lhs, c);
        CLPOLY_ASSERT(!c.empty());
    }

    CLPOLY_TEST("polynomial_GCD(ZZ): coprime");
    {
        // a = x^2 + 1, b = x + 1
        upolynomial_<ZZ> a({{umonomial(2), ZZ(1)}, {umonomial(0), ZZ(1)}});
        upolynomial_<ZZ> b({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});
        upolynomial_<ZZ> s, t;
        auto c = polynomial_GCD(a, b, s, t);
        auto lhs = s * a + t * b;
        lhs.normalization();
        CLPOLY_ASSERT_EQ(lhs, c);
    }

    CLPOLY_TEST("polynomial_GCD(ZZ): non-monic coprime");
    {
        // a = 2x + 1, b = 3x + 1
        upolynomial_<ZZ> a({{umonomial(1), ZZ(2)}, {umonomial(0), ZZ(1)}});
        upolynomial_<ZZ> b({{umonomial(1), ZZ(3)}, {umonomial(0), ZZ(1)}});
        upolynomial_<ZZ> s, t;
        auto c = polynomial_GCD(a, b, s, t);
        auto lhs = s * a + t * b;
        lhs.normalization();
        CLPOLY_ASSERT_EQ(lhs, c);
    }

    CLPOLY_TEST("polynomial_GCD(ZZ): non-coprime");
    {
        // a = x^2 - 1 = (x-1)(x+1), b = x - 1, gcd = x - 1
        upolynomial_<ZZ> a({{umonomial(2), ZZ(1)}, {umonomial(0), ZZ(-1)}});
        upolynomial_<ZZ> b({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(-1)}});
        upolynomial_<ZZ> s, t;
        auto c = polynomial_GCD(a, b, s, t);
        auto lhs = s * a + t * b;
        lhs.normalization();
        CLPOLY_ASSERT_EQ(lhs, c);
        // gcd should be x - 1 (lc > 0)
        CLPOLY_ASSERT_EQ(c.front().second, ZZ(1));
        CLPOLY_ASSERT_EQ(get_deg(c), 1);
    }

    CLPOLY_TEST("polynomial_GCD(Zp): coprime");
    {
        // a = x^2 + 1, b = x + 1 in Z/7Z[x]
        uint32_t p = 7;
        upolynomial_<Zp> a({{umonomial(2), Zp(1,p)}, {umonomial(0), Zp(1,p)}});
        upolynomial_<Zp> b({{umonomial(1), Zp(1,p)}, {umonomial(0), Zp(1,p)}});
        upolynomial_<Zp> s, t;
        auto c = polynomial_GCD(a, b, s, t);
        auto lhs = s * a + t * b;
        lhs.normalization();
        CLPOLY_ASSERT_EQ(lhs, c);
        // coprime → gcd = 1
        CLPOLY_ASSERT_EQ(get_deg(c), 0);
    }

    CLPOLY_TEST("polynomial_GCD(Zp): non-coprime");
    {
        // a = x^2 - 1, b = x - 1 in Z/7Z, gcd = x - 1
        uint32_t p = 7;
        upolynomial_<Zp> a({{umonomial(2), Zp(1,p)}, {umonomial(0), Zp(-1,p)}});
        upolynomial_<Zp> b({{umonomial(1), Zp(1,p)}, {umonomial(0), Zp(-1,p)}});
        upolynomial_<Zp> s, t;
        auto c = polynomial_GCD(a, b, s, t);
        auto lhs = s * a + t * b;
        lhs.normalization();
        CLPOLY_ASSERT_EQ(lhs, c);
        // gcd = x + 6 (monic, -1 ≡ 6 mod 7)
        CLPOLY_ASSERT_EQ(get_deg(c), 1);
        CLPOLY_ASSERT_EQ(c.front().second.number(), 1);  // monic
    }

    CLPOLY_TEST("polynomial_GCD(ZZ): Bezout chain for 3 factors");
    {
        // v1 = 2x+1, v2 = 2x-1, v3 = 2x+3
        upolynomial_<ZZ> v1({{umonomial(1), ZZ(2)}, {umonomial(0), ZZ(1)}});
        upolynomial_<ZZ> v2({{umonomial(1), ZZ(2)}, {umonomial(0), ZZ(-1)}});
        upolynomial_<ZZ> v3({{umonomial(1), ZZ(2)}, {umonomial(0), ZZ(3)}});

        std::vector<upolynomial_<ZZ>> sv(3);
        ZZ denom(1);
        upolynomial_<ZZ> g_acc = v1;
        sv[0] = upolynomial_<ZZ>({{umonomial(0), ZZ(1)}});

        {
            upolynomial_<ZZ> alpha, beta;
            auto c_poly = polynomial_GCD(g_acc, v2, alpha, beta);
            // c_poly 应为常数 (coprime)
            assert(c_poly.size() == 1 && c_poly.front().first.deg() == 0);
            ZZ c_k = c_poly.front().second;
            sv[0] = sv[0] * beta;
            sv[0].normalization();
            auto tmp = upolynomial_<ZZ>({{umonomial(0), denom}});
            sv[1] = tmp * alpha;
            sv[1].normalization();
            denom *= c_k;
            g_acc = g_acc * v2;
            g_acc.normalization();
        }
        {
            upolynomial_<ZZ> alpha, beta;
            auto c_poly = polynomial_GCD(g_acc, v3, alpha, beta);
            ZZ c_k = c_poly.front().second;
            sv[0] = sv[0] * beta;
            sv[0].normalization();
            sv[1] = sv[1] * beta;
            sv[1].normalization();
            auto tmp = upolynomial_<ZZ>({{umonomial(0), denom}});
            sv[2] = tmp * alpha;
            sv[2].normalization();
            denom *= c_k;
        }

        auto vh0 = v2 * v3; vh0.normalization();
        auto vh1 = v1 * v3; vh1.normalization();
        auto vh2 = v1 * v2; vh2.normalization();

        auto sum = sv[0] * vh0 + sv[1] * vh1 + sv[2] * vh2;
        sum.normalization();
        upolynomial_<ZZ> expected({{umonomial(0), denom}});
        CLPOLY_ASSERT_EQ(sum, expected);
    }

    // ============================================================
    CLPOLY_TEST("polynomial_GCD(multivar): bivariate coprime");
    {
        // f = x^2 + y, g = x + y (coprime in Z[x,y])
        auto f = make_lex(pow(x,2) + y);
        auto g = make_lex(x + y);
        PolyZZ s, t;
        auto c = polynomial_GCD(f, g, s, t);
        // s*f + t*g = c
        auto lhs = s * f + t * g;
        lhs.normalization();
        CLPOLY_ASSERT_EQ(lhs, c);
    }

    CLPOLY_TEST("polynomial_GCD(multivar): common factor");
    {
        // f = (x+y)(x-y) = x^2 - y^2, g = (x+y)(x+1) = x^2 + x + xy + y
        // gcd should be proportional to x+y
        auto f = make_lex(pow(x,2) - pow(y,2));
        auto g = make_lex(pow(x,2) + x + x*y + y);
        PolyZZ s, t;
        auto c = polynomial_GCD(f, g, s, t);
        auto lhs = s * f + t * g;
        lhs.normalization();
        CLPOLY_ASSERT_EQ(lhs, c);
        // c should have degree 1 in x
        CLPOLY_ASSERT_EQ(degree(c, x), 1);
    }

    CLPOLY_TEST("gcd wrapper: multivar extended");
    {
        auto f = make_lex(pow(x,2) + y);
        auto g = make_lex(x + y);
        PolyZZ s, t;
        auto c1 = gcd(f, g, s, t);
        auto lhs = s * f + t * g;
        lhs.normalization();
        CLPOLY_ASSERT_EQ(lhs, c1);
    }

    CLPOLY_TEST("gcd wrapper: upolynomial extended");
    {
        upolynomial_<ZZ> a({{umonomial(2), ZZ(1)}, {umonomial(0), ZZ(1)}});
        upolynomial_<ZZ> b({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});
        upolynomial_<ZZ> s, t;
        auto c = gcd(a, b, s, t);
        auto lhs = s * a + t * b;
        lhs.normalization();
        CLPOLY_ASSERT_EQ(lhs, c);
    }

    // ============================================================
    CLPOLY_TEST("__taylor_coeff: univariate");
    {
        // f = x^3, alpha = 1
        // f = (x-1)^3 + 3(x-1)^2 + 3(x-1) + 1
        auto f = make_lex(pow(x, 3));
        ZZ alpha(1);

        auto c0 = __taylor_coeff(f, x, alpha, 0);
        auto c1 = __taylor_coeff(f, x, alpha, 1);
        auto c2 = __taylor_coeff(f, x, alpha, 2);
        auto c3 = __taylor_coeff(f, x, alpha, 3);

        auto one = make_lex(polynomial_ZZ(1));
        auto three = make_lex(polynomial_ZZ(3));
        CLPOLY_ASSERT_EQ(c0, one);
        CLPOLY_ASSERT_EQ(c1, three);
        CLPOLY_ASSERT_EQ(c2, three);
        CLPOLY_ASSERT_EQ(c3, one);
    }

    CLPOLY_TEST("__taylor_coeff: bivariate");
    {
        // f = x^2*y + x*y^2
        // Expand in (y - 0) = y:
        // j=0: assign(f, y, 0) = 0
        // j=1: f/y = x^2 + x*y, assign(that, y, 0) = x^2
        // j=2: (x^2 + x*y)/y = x, assign(that, y, 0) = x
        auto f = make_lex(pow(x,2)*y + x*pow(y,2));
        ZZ alpha(0);

        auto c0 = __taylor_coeff(f, y, alpha, 0);
        auto c1 = __taylor_coeff(f, y, alpha, 1);
        auto c2 = __taylor_coeff(f, y, alpha, 2);

        PolyZZ zero;
        CLPOLY_ASSERT_EQ(c0, zero);
        CLPOLY_ASSERT_EQ(c1, make_lex(pow(x,2)));
        CLPOLY_ASSERT_EQ(c2, make_lex(x));
    }

    CLPOLY_TEST("__taylor_coeff: reconstruction");
    {
        // f = x^2 + 2*x*y + y^2
        // Σ __taylor_coeff(f, y, 1, j) * (y-1)^j == f
        auto f = make_lex(pow(x,2) + 2*x*y + pow(y,2));
        ZZ alpha(1);

        auto y_lex = make_lex(y);
        auto one_lex = make_lex(polynomial_ZZ(1));
        auto y_minus_a = y_lex - one_lex;

        PolyZZ sum;
        PolyZZ power = one_lex;  // (y-1)^0 = 1

        for (int j = 0; j <= 2; ++j)
        {
            auto cj = __taylor_coeff(f, y, alpha, j);
            sum = sum + cj * power;
            sum.normalization();
            power = power * y_minus_a;
            power.normalization();
        }
        CLPOLY_ASSERT_EQ(sum, f);
    }

    CLPOLY_TEST("__taylor_coeff: non-zero alpha bivariate");
    {
        // f = (x + y - 2)(x - y + 3) = x^2 + x + 5y - y^2 - 6
        auto f = make_lex(pow(x,2) + x + 5*y - pow(y,2) - polynomial_ZZ(6));
        ZZ alpha(2);
        int deg_y = 2;

        auto y_lex = make_lex(y);
        auto two_lex = make_lex(polynomial_ZZ(2));
        auto one_lex = make_lex(polynomial_ZZ(1));
        auto y_minus_a = y_lex - two_lex;

        PolyZZ sum;
        PolyZZ power = one_lex;

        for (int j = 0; j <= deg_y; ++j)
        {
            auto cj = __taylor_coeff(f, y, alpha, j);
            sum = sum + cj * power;
            sum.normalization();
            power = power * y_minus_a;
            power.normalization();
        }
        CLPOLY_ASSERT_EQ(sum, f);
    }

    return clpoly_test::test_summary();
}
