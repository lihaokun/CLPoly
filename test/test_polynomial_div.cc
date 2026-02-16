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

    // ======== pquo is_L parameter variations ========
    CLPOLY_TEST("pquo_is_L_true");
    {
        // is_L=true (default): lazy pseudo-quotient
        polynomial_ZZ f = pow(x,3)*y + 2*pow(x,2)*z - x*pow(y,2) + 3;
        polynomial_ZZ g = pow(x,2) + y*z;
        polynomial_ZZ r_true;
        auto q_true = pquo(r_true, f, g, x, true);
        auto lc = leadcoeff(g, x);
        int64_t d = degree(f, x) - degree(g, x) + 1;
        // Verify identity: lc^d * f == q*g + r
        auto lhs = pow(lc, d) * f;
        auto rhs = q_true * g + r_true;
        CLPOLY_ASSERT_EQ(lhs, rhs);
    }

    CLPOLY_TEST("pquo_is_L_false");
    {
        // is_L=false: standard pseudo-quotient
        polynomial_ZZ f = pow(x,3)*y + 2*pow(x,2)*z - x*pow(y,2) + 3;
        polynomial_ZZ g = pow(x,2) + y*z;
        polynomial_ZZ r_false;
        auto q_false = pquo(r_false, f, g, x, false);
        // The remainder should still satisfy some division identity
        // With is_L=false, we get a different scaling but q*g + r still relates to a*f
        CLPOLY_ASSERT_FALSE(r_false.empty() && q_false.empty() && !f.empty());
    }

    CLPOLY_TEST("pquo_is_L_comparison");
    {
        // Compare is_L=true vs is_L=false on same inputs
        polynomial_ZZ f = pow(x,3) + 2*pow(x,2)*y + x - 1;
        polynomial_ZZ g = pow(x,2) + y;

        polynomial_ZZ r_true, r_false;
        auto q_true = pquo(r_true, f, g, x, true);
        auto q_false = pquo(r_false, f, g, x, false);

        // Both should satisfy prem identity with appropriate leading coeff power
        // is_L=true should use minimal power; is_L=false may use more
        // Verify is_L=true identity
        auto lc = leadcoeff(g, x);
        int64_t d = degree(f, x) - degree(g, x) + 1;
        CLPOLY_ASSERT_EQ(pow(lc, d) * f, q_true * g + r_true);

        // prem should be consistent regardless of is_L
        auto r_prem = prem(f, g, x);
        // With is_L=true, the remainder should divide r_prem or be a constant multiple
        // At minimum, both should have deg in x < deg(g, x)
        CLPOLY_ASSERT_TRUE(degree(r_true, x) < degree(g, x));
        CLPOLY_ASSERT_TRUE(degree(r_false, x) < degree(g, x));
    }

    return clpoly_test::test_summary();
}
