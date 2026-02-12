#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include "testdata_expected.hh"

using namespace clpoly;

// Helper: check if a divides b (b/a has zero remainder)
bool divides(const polynomial_ZZ& a, const polynomial_ZZ& b) {
    if (a.empty()) return b.empty();
    if (b.empty()) return true;
    auto vars = get_variables(a);
    if (vars.empty()) return true;
    // Try polynomial division
    auto q = b / a;
    return q * a == b;
}

int main() {
    variable x("x"), y("y"), z("z");

    // ======== Simple GCD: known common factor ========
    CLPOLY_TEST("gcd_known_factor");
    {
        auto f = testdata::gcd_f1();
        auto g = testdata::gcd_g1();
        auto expected = testdata::gcd_result1();
        auto result = gcd(f, g);
        // GCD is unique up to sign/associate, so check divisibility
        CLPOLY_ASSERT(divides(result, f));
        CLPOLY_ASSERT(divides(result, g));
        // The degree should match
        CLPOLY_ASSERT_EQ(degree(result), degree(expected));
    }

    // ======== Coprime polynomials ========
    CLPOLY_TEST("gcd_coprime");
    {
        auto f = testdata::gcd_f2_coprime();
        auto g = testdata::gcd_g2_coprime();
        auto result = gcd(f, g);
        CLPOLY_ASSERT(is_number(result));
    }

    // ======== GCD with self ========
    CLPOLY_TEST("gcd_self");
    {
        auto f = testdata::arith_f1_simple();
        auto result = gcd(f, f);
        CLPOLY_ASSERT(divides(result, f));
        CLPOLY_ASSERT(divides(f, result));
    }

    // ======== GCD with zero ========
    CLPOLY_TEST("gcd_with_zero");
    {
        auto f = testdata::arith_f1_simple();
        polynomial_ZZ zero;
        auto result = gcd(f, zero);
        CLPOLY_ASSERT(divides(result, f));
        CLPOLY_ASSERT(divides(f, result));
    }

    // ======== Multivariate GCD ========
    CLPOLY_TEST("gcd_multivariate");
    {
        auto f = testdata::gcd_f3_mv();
        auto g = testdata::gcd_g3_mv();
        auto expected = testdata::gcd_result3_mv();
        auto result = gcd(f, g);
        CLPOLY_ASSERT(divides(result, f));
        CLPOLY_ASSERT(divides(result, g));
        CLPOLY_ASSERT_EQ(degree(result), degree(expected));
    }

    // ======== GCD algebraic identity: gcd(f*h, g*h) is divisible by h ========
    CLPOLY_TEST("gcd_identity");
    {
        polynomial_ZZ f = pow(x,2) + 1;
        polynomial_ZZ g = pow(x,3) - x + 2;
        polynomial_ZZ h = x*y + z - 1;
        auto fh = f * h;
        auto gh = g * h;
        auto result = gcd(fh, gh);
        CLPOLY_ASSERT(divides(h, result));
    }

    // ======== Complex multivariate GCD ========
    CLPOLY_TEST("gcd_complex");
    {
        auto f = testdata::complex_gcd_f();
        auto g = testdata::complex_gcd_g();
        auto common = testdata::complex_gcd_common();
        auto result = gcd(f, g);
        CLPOLY_ASSERT(divides(result, f));
        CLPOLY_ASSERT(divides(result, g));
        CLPOLY_ASSERT(divides(common, result));
    }

    // ======== Squarefree decomposition ========
    CLPOLY_TEST("squarefree_simple");
    {
        // (x+1)^2 * (x-2): not squarefree
        auto f = testdata::sqf_f1();
        CLPOLY_ASSERT_FALSE(is_squarefree(f));
    }

    CLPOLY_TEST("squarefree_higher_mult");
    {
        // (x+1)^2 * (x-1)^3: not squarefree
        auto f = testdata::sqf_f2();
        CLPOLY_ASSERT_FALSE(is_squarefree(f));
    }

    CLPOLY_TEST("squarefree_already");
    {
        // x^3 + x + 1: already squarefree
        auto f = testdata::sqf_f3_already();
        CLPOLY_ASSERT_TRUE(is_squarefree(f));
    }

    CLPOLY_TEST("squarefree_multivariate");
    {
        // (x+y)^2 * (x-y+1): not squarefree
        auto f = testdata::sqf_f4_mv();
        CLPOLY_ASSERT_FALSE(is_squarefree(f));
    }

    // ======== Squarefree factorization: product should equal original ========
    CLPOLY_TEST("squarefree_factorize_identity");
    {
        auto f = testdata::sqf_f2();
        auto factors = squarefreefactorize(f);
        CLPOLY_ASSERT(factors.size() > 0);
        // Reconstruct: product of fi^mi should equal f (up to constant)
        polynomial_ZZ product({{monomial(), ZZ(1)}});
        for (auto& p : factors) {
            product = product * pow(p.first, p.second);
        }
        // Check divisibility both ways (up to constant multiple)
        CLPOLY_ASSERT(divides(product, f));
        CLPOLY_ASSERT(divides(f, product));
    }

    CLPOLY_TEST("squarefree_factorize_each_factor_is_squarefree");
    {
        auto f = testdata::sqf_f2();
        auto factors = squarefreefactorize(f);
        for (auto& p : factors) {
            if (!is_number(p.first)) {
                CLPOLY_ASSERT_TRUE(is_squarefree(p.first));
            }
        }
    }

    return clpoly_test::test_summary();
}
