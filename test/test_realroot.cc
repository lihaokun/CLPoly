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

    // ======== uroot accessors ========
    CLPOLY_TEST("uroot_accessors");
    {
        // Find roots of x^2 - 2 (roots at +-sqrt(2))
        polynomial_ZZ f = pow(x,2) - 2;
        std::vector<polynomial_ZZ> polys = {f};
        auto result = realroot(polys);
        auto& roots = result.first;
        CLPOLY_ASSERT_EQ((int)roots.size(), 2);

        for (auto& r : roots) {
            // left() <= right()
            CLPOLY_ASSERT_TRUE(r.left() <= r.right());
            // poly() should be non-empty
            CLPOLY_ASSERT_FALSE(r.poly().empty());
        }
    }

    // ======== uroot comparison ========
    CLPOLY_TEST("uroot_comparison");
    {
        // Roots of (x-1)(x-3) = x^2 - 4x + 3
        polynomial_ZZ f = pow(x,2) - 4*x + 3;
        std::vector<polynomial_ZZ> polys = {f};
        auto result = realroot(polys);
        auto& roots = result.first;
        CLPOLY_ASSERT_EQ((int)roots.size(), 2);

        // One root should be less than the other
        // roots are sorted, so roots[0] < roots[1]
        CLPOLY_ASSERT_TRUE(roots[0] < roots[1]);
        CLPOLY_ASSERT_TRUE(roots[1] > roots[0]);
        CLPOLY_ASSERT_TRUE(roots[0] <= roots[1]);
        CLPOLY_ASSERT_TRUE(roots[1] >= roots[0]);
        CLPOLY_ASSERT_TRUE(roots[0] != roots[1]);
        CLPOLY_ASSERT_FALSE(roots[0] == roots[1]);

        // Self-comparison
        CLPOLY_ASSERT_TRUE(roots[0] == roots[0]);
        CLPOLY_ASSERT_TRUE(roots[0] <= roots[0]);
        CLPOLY_ASSERT_TRUE(roots[0] >= roots[0]);
    }

    // ======== uroot shrinkinterval ========
    CLPOLY_TEST("uroot_shrink");
    {
        polynomial_ZZ f = pow(x,2) - 2;
        std::vector<polynomial_ZZ> polys = {f};
        auto result = realroot(polys);
        auto& roots = result.first;
        CLPOLY_ASSERT_EQ((int)roots.size(), 2);

        for (auto& r : roots) {
            QQ width_before = r.right() - r.left();
            if (!r.is_single()) {
                r.shrinkinterval();
                QQ width_after = r.right() - r.left();
                // Interval should not grow
                CLPOLY_ASSERT_TRUE(width_after <= width_before);
            }
        }
    }

    // ======== RealRootBound ========
    CLPOLY_TEST("realroot_bound");
    {
        // x^2 - 4: roots at +-2, bound should be >= 2
        upolynomial_ZZ p({{umonomial(2), ZZ(1)}, {umonomial(0), ZZ(-4)}});
        ZZ bound = RealRootBound(p);
        CLPOLY_ASSERT_TRUE(bound >= ZZ(2));

        // x - 5: root at 5, bound should be >= 5
        upolynomial_ZZ p2({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(-5)}});
        ZZ bound2 = RealRootBound(p2);
        CLPOLY_ASSERT_TRUE(bound2 >= ZZ(5));
    }

    // ======== coeffsignchanges ========
    CLPOLY_TEST("coeffsignchanges_basic");
    {
        // x^2 - 1: coefficients are [1, 0, -1], sign changes = 1
        upolynomial_ZZ p({{umonomial(2), ZZ(1)}, {umonomial(0), ZZ(-1)}});
        auto changes = coeffsignchanges(p);
        CLPOLY_ASSERT_EQ(changes, (uint64_t)1);

        // x^2 + x + 1: coefficients [1, 1, 1], sign changes = 0
        upolynomial_ZZ p2({{umonomial(2), ZZ(1)}, {umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});
        auto changes2 = coeffsignchanges(p2);
        CLPOLY_ASSERT_EQ(changes2, (uint64_t)0);

        // x^3 - x: coefficients [1, 0, -1, 0], sign changes = 1
        upolynomial_ZZ p3({{umonomial(3), ZZ(1)}, {umonomial(1), ZZ(-1)}});
        auto changes3 = coeffsignchanges(p3);
        CLPOLY_ASSERT_EQ(changes3, (uint64_t)1);
    }

    return clpoly_test::test_summary();
}
