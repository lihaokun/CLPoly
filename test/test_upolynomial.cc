#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

int main() {
    // ======== Construction ========
    CLPOLY_TEST("upoly_construction");
    {
        // 3x^2 + 2x + 1
        upolynomial_ZZ p({{umonomial(2), ZZ(3)}, {umonomial(1), ZZ(2)}, {umonomial(0), ZZ(1)}});
        CLPOLY_ASSERT_EQ((int)p.size(), 3);
        CLPOLY_ASSERT_EQ(degree(p), (int64_t)2);
    }

    // ======== Arithmetic ========
    CLPOLY_TEST("upoly_add");
    {
        upolynomial_ZZ p({{umonomial(2), ZZ(3)}, {umonomial(0), ZZ(1)}});   // 3x^2 + 1
        upolynomial_ZZ q({{umonomial(2), ZZ(-1)}, {umonomial(1), ZZ(2)}});  // -x^2 + 2x
        auto sum = p + q;
        // 2x^2 + 2x + 1
        CLPOLY_ASSERT_EQ(degree(sum), (int64_t)2);
        CLPOLY_ASSERT_EQ((int)sum.size(), 3);
    }

    CLPOLY_TEST("upoly_mul");
    {
        upolynomial_ZZ p({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});   // x + 1
        upolynomial_ZZ q({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(-1)}});  // x - 1
        auto product = p * q;
        // x^2 - 1
        upolynomial_ZZ expected({{umonomial(2), ZZ(1)}, {umonomial(0), ZZ(-1)}});
        CLPOLY_ASSERT_EQ(product, expected);
    }

    // ======== Evaluation ========
    CLPOLY_TEST("upoly_eval");
    {
        // 2x^3 - 3x + 1
        upolynomial_ZZ p({{umonomial(3), ZZ(2)}, {umonomial(1), ZZ(-3)}, {umonomial(0), ZZ(1)}});
        // Evaluate at x=2: 2*8 - 3*2 + 1 = 16 - 6 + 1 = 11
        ZZ result = assign(p, ZZ(2));
        CLPOLY_ASSERT_EQ(result, ZZ(11));
    }

    CLPOLY_TEST("upoly_eval_zero");
    {
        upolynomial_ZZ p({{umonomial(3), ZZ(2)}, {umonomial(1), ZZ(-3)}, {umonomial(0), ZZ(1)}});
        ZZ result = assign(p, ZZ(0));
        CLPOLY_ASSERT_EQ(result, ZZ(1));
    }

    // ======== Derivative ========
    CLPOLY_TEST("upoly_derivative");
    {
        // 3x^3 + 2x^2 - x + 5
        upolynomial_ZZ p({{umonomial(3), ZZ(3)}, {umonomial(2), ZZ(2)}, {umonomial(1), ZZ(-1)}, {umonomial(0), ZZ(5)}});
        auto dp = derivative(p);
        // 9x^2 + 4x - 1
        upolynomial_ZZ expected({{umonomial(2), ZZ(9)}, {umonomial(1), ZZ(4)}, {umonomial(0), ZZ(-1)}});
        CLPOLY_ASSERT_EQ(dp, expected);
    }

    // ======== Conversion: multivariate -> univariate ========
    CLPOLY_TEST("upoly_conversion");
    {
        variable x("x");
        polynomial_ZZ f = 2*pow(x,3) - 3*pow(x,2) + x - 5;
        upolynomial_ZZ uf;
        poly_convert(f, uf);
        CLPOLY_ASSERT_EQ(degree(uf), (int64_t)3);
        // Evaluate at x=1: 2 - 3 + 1 - 5 = -5
        ZZ val = assign(uf, ZZ(1));
        ZZ expected = ZZ(-5);
        CLPOLY_ASSERT_EQ(val, expected);
    }

    // ======== Power ========
    CLPOLY_TEST("upoly_pow");
    {
        upolynomial_ZZ p({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});  // x + 1
        auto p2 = pow(p, 2);
        auto p_times_p = p * p;
        CLPOLY_ASSERT_EQ(p2, p_times_p);
    }

    // ======== Content ========
    CLPOLY_TEST("upoly_content");
    {
        // 6x^2 + 4x + 2 -> content = 2
        upolynomial_ZZ p({{umonomial(2), ZZ(6)}, {umonomial(1), ZZ(4)}, {umonomial(0), ZZ(2)}});
        ZZ c = cont(p);
        CLPOLY_ASSERT_EQ(c, ZZ(2));
    }

    return clpoly_test::test_summary();
}
