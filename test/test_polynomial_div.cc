#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include "testdata_expected.hh"

using namespace clpoly;

int main() {
    variable x("x"), y("y"), z("z");

    // ======== Simple prem/pquo ========
    CLPOLY_TEST("prem_simple");
    {
        auto f = testdata::prem_f1();
        auto g = testdata::prem_g1();
        auto r = prem(f, g, x);
        CLPOLY_ASSERT_EQ(r, testdata::prem_result1());
    }

    CLPOLY_TEST("pquo_simple");
    {
        auto f = testdata::prem_f1();
        auto g = testdata::prem_g1();
        polynomial_ZZ r;
        auto q = pquo(r, f, g, x);
        CLPOLY_ASSERT_EQ(q, testdata::pquo_result1());
        CLPOLY_ASSERT_EQ(r, testdata::prem_result1());
    }

    // ======== Multivariate prem/pquo ========
    CLPOLY_TEST("prem_multivariate");
    {
        auto f = testdata::prem_f2();
        auto g = testdata::prem_g2();
        auto r = prem(f, g, x);
        CLPOLY_ASSERT_EQ(r, testdata::prem_result2());
    }

    CLPOLY_TEST("pquo_multivariate");
    {
        auto f = testdata::prem_f2();
        auto g = testdata::prem_g2();
        polynomial_ZZ r;
        auto q = pquo(r, f, g, x);
        CLPOLY_ASSERT_EQ(q, testdata::pquo_result2());
        CLPOLY_ASSERT_EQ(r, testdata::prem_result2());
    }

    // ======== prem/pquo identity: lc(g,x)^d * f == q*g + r ========
    CLPOLY_TEST("prem_pquo_identity_simple");
    {
        auto f = testdata::prem_f1();
        auto g = testdata::prem_g1();
        polynomial_ZZ r;
        auto q = pquo(r, f, g, x);
        auto lc = leadcoeff(g, x);
        int64_t d = degree(f, x) - degree(g, x) + 1;
        auto lhs = pow(lc, d) * f;
        auto rhs = q * g + r;
        CLPOLY_ASSERT_EQ(lhs, rhs);
    }

    CLPOLY_TEST("prem_pquo_identity_multivariate");
    {
        auto f = testdata::prem_f2();
        auto g = testdata::prem_g2();
        polynomial_ZZ r;
        auto q = pquo(r, f, g, x);
        auto lc = leadcoeff(g, x);
        int64_t d = degree(f, x) - degree(g, x) + 1;
        auto lhs = pow(lc, d) * f;
        auto rhs = q * g + r;
        CLPOLY_ASSERT_EQ(lhs, rhs);
    }

    // ======== Identity with hand-constructed polynomials ========
    CLPOLY_TEST("prem_pquo_identity_3var");
    {
        polynomial_ZZ f = pow(x,3)*y + 2*pow(x,2)*z - x*pow(y,2) + 3;
        polynomial_ZZ g = pow(x,2) + y*z;
        polynomial_ZZ r;
        auto q = pquo(r, f, g, x);
        auto lc = leadcoeff(g, x);
        int64_t d = degree(f, x) - degree(g, x) + 1;
        auto lhs = pow(lc, d) * f;
        auto rhs = q * g + r;
        CLPOLY_ASSERT_EQ(lhs, rhs);
    }

    // ======== prem when deg(f) < deg(g) ========
    CLPOLY_TEST("prem_deg_less");
    {
        polynomial_ZZ f = x + 1;
        polynomial_ZZ g = pow(x,3) + x;
        auto r = prem(f, g, x);
        CLPOLY_ASSERT_EQ(r, f);
    }

    // ======== prem by constant ========
    CLPOLY_TEST("prem_by_constant");
    {
        polynomial_ZZ f = pow(x,2) + 1;
        polynomial_ZZ g({{monomial(), ZZ(3)}});
        polynomial_ZZ zero;
        auto r = prem(f, g, x);
        CLPOLY_ASSERT_EQ(r, zero);
    }

    // ======== Exact division ========
    CLPOLY_TEST("exact_division");
    {
        polynomial_ZZ f = pow(x,2) - 1;
        polynomial_ZZ g = x + 1;
        auto result = f / g;
        polynomial_ZZ expected = x - 1;
        CLPOLY_ASSERT_EQ(result, expected);
    }

    CLPOLY_TEST("exact_division_multivariate");
    {
        polynomial_ZZ a = x*y + 1;
        polynomial_ZZ b = pow(x,2) - y;
        polynomial_ZZ product = a * b;
        CLPOLY_ASSERT_EQ(product / a, b);
        CLPOLY_ASSERT_EQ(product / b, a);
    }

    // ======== Monomial division ========
    CLPOLY_TEST("monomial_division");
    {
        polynomial_ZZ f = 6*pow(x,3)*pow(y,2) - 3*pow(x,2)*y;
        polynomial_ZZ d({{pow(x,1)*y, ZZ(3)}});
        auto result = f / d;
        polynomial_ZZ expected = 2*pow(x,2)*y - x;
        CLPOLY_ASSERT_EQ(result, expected);
    }

    return clpoly_test::test_summary();
}
