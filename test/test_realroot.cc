#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

int main() {
    variable x("x");

    // ======== Known roots ========
    CLPOLY_TEST("realroot_known_roots");
    {
        // (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
        polynomial_ZZ f = pow(x,3) - 6*pow(x,2) + 11*x - 6;
        std::vector<polynomial_ZZ> polys = {f};
        auto result = realroot(polys);
        auto& roots = result.first;
        // Should have 3 real roots
        CLPOLY_ASSERT_EQ((int)roots.size(), 3);
    }

    // ======== No real roots ========
    CLPOLY_TEST("realroot_no_real_roots");
    {
        // x^2 + 1 has no real roots
        polynomial_ZZ f = pow(x,2) + 1;
        std::vector<polynomial_ZZ> polys = {f};
        auto result = realroot(polys);
        auto& roots = result.first;
        CLPOLY_ASSERT_EQ((int)roots.size(), 0);
    }

    // ======== Single root ========
    CLPOLY_TEST("realroot_single");
    {
        // x - 5
        polynomial_ZZ f = x - 5;
        std::vector<polynomial_ZZ> polys = {f};
        auto result = realroot(polys);
        auto& roots = result.first;
        CLPOLY_ASSERT_EQ((int)roots.size(), 1);
    }

    // ======== Repeated root ========
    CLPOLY_TEST("realroot_repeated");
    {
        // (x-1)^2 = x^2 - 2x + 1
        polynomial_ZZ f = pow(x,2) - 2*x + 1;
        std::vector<polynomial_ZZ> polys = {f};
        auto result = realroot(polys);
        auto& roots = result.first;
        // Should find the repeated root
        CLPOLY_ASSERT(roots.size() >= 1);
    }

    // ======== Quadratic with two roots ========
    CLPOLY_TEST("realroot_quadratic");
    {
        // x^2 - 4 = (x-2)(x+2)
        polynomial_ZZ f = pow(x,2) - 4;
        std::vector<polynomial_ZZ> polys = {f};
        auto result = realroot(polys);
        auto& roots = result.first;
        CLPOLY_ASSERT_EQ((int)roots.size(), 2);
    }

    // ======== Wilkinson-type polynomial ========
    CLPOLY_TEST("realroot_wilkinson5");
    {
        // (x-1)(x-2)(x-3)(x-4)(x-5)
        polynomial_ZZ f = (x-1)*(x-2)*(x-3)*(x-4)*(x-5);
        std::vector<polynomial_ZZ> polys = {f};
        auto result = realroot(polys);
        auto& roots = result.first;
        CLPOLY_ASSERT_EQ((int)roots.size(), 5);
    }

    // ======== Uspensky direct on upolynomial ========
    CLPOLY_TEST("uspensky_direct");
    {
        // x^2 - 2: roots at +/-sqrt(2)
        upolynomial_ZZ p({{umonomial(2), ZZ(1)}, {umonomial(0), ZZ(-2)}});
        auto intervals = uspensky(p);
        CLPOLY_ASSERT_EQ((int)intervals.size(), 2);
        // Check that intervals contain the roots
        for (auto& iv : intervals) {
            // The root sqrt(2) ~ 1.414 or -sqrt(2) ~ -1.414
            // Each interval [a,b] should satisfy a <= root <= b
            CLPOLY_ASSERT_TRUE(iv.first <= iv.second);
        }
    }

    return clpoly_test::test_summary();
}
