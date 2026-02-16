/**
 * @file test_crosscheck_ntl.cc
 * @brief Cross-library correctness tests: CLPoly upolynomial_ZZ vs NTL ZZX.
 *
 * Random univariate polynomials are generated in CLPoly, converted to NTL,
 * operations are performed in both libraries, and results are compared.
 */
#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include "crosscheck_ntl.hh"

using namespace clpoly;

int main() {

    // ======== Round-trip: CLPoly -> NTL -> CLPoly ========
    CLPOLY_TEST("crosscheck_ntl_roundtrip");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(7, 5, {-100, 100});
        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(f);
        auto f_back = crosscheck::ntl_to_clpoly_upoly(ntl_f);
        CLPOLY_ASSERT_EQ(f, f_back);
    }

    // ======== Addition: f + g ========
    CLPOLY_TEST("crosscheck_ntl_add");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(6, 4, {-100, 100});
        auto g = random_upolynomial<ZZ>(6, 4, {-100, 100});

        auto clpoly_result = f + g;

        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
        NTL::ZZX nr = nf + ng;
        auto ntl_result = crosscheck::ntl_to_clpoly_upoly(nr);

        CLPOLY_ASSERT_EQ(clpoly_result, ntl_result);
    }

    // ======== Subtraction: f - g ========
    CLPOLY_TEST("crosscheck_ntl_sub");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(6, 4, {-100, 100});
        auto g = random_upolynomial<ZZ>(6, 4, {-100, 100});

        auto clpoly_result = f - g;

        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
        NTL::ZZX nr = nf - ng;
        auto ntl_result = crosscheck::ntl_to_clpoly_upoly(nr);

        CLPOLY_ASSERT_EQ(clpoly_result, ntl_result);
    }

    // ======== Multiplication: f * g ========
    CLPOLY_TEST("crosscheck_ntl_mul");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(5, 4, {-100, 100});
        auto g = random_upolynomial<ZZ>(5, 4, {-100, 100});

        auto clpoly_result = f * g;

        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
        NTL::ZZX nr = nf * ng;
        auto ntl_result = crosscheck::ntl_to_clpoly_upoly(nr);

        CLPOLY_ASSERT_EQ(clpoly_result, ntl_result);
    }

    // ======== GCD (up to sign) ========
    CLPOLY_TEST("crosscheck_ntl_gcd");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        // Construct f = h * a, g = h * b so gcd is at least h
        auto h = random_upolynomial<ZZ>(3, 3, {-20, 20});
        auto a = random_upolynomial<ZZ>(3, 3, {-20, 20});
        auto b = random_upolynomial<ZZ>(3, 3, {-20, 20});
        auto f = h * a;
        auto g = h * b;

        if (f.empty() || g.empty()) continue;

        // CLPoly gcd
        auto clpoly_gcd = polynomial_GCD(f, g);

        // NTL gcd
        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
        NTL::ZZX ngcd = NTL::GCD(nf, ng);
        auto ntl_gcd = crosscheck::ntl_to_clpoly_upoly(ngcd);

        // GCD is unique up to sign over ZZ:
        // 1. Degrees should match
        CLPOLY_ASSERT_EQ(degree(clpoly_gcd), degree(ntl_gcd));
        // 2. Leading coefficients should have same absolute value
        if (!clpoly_gcd.empty() && !ntl_gcd.empty()) {
            CLPOLY_ASSERT(clpoly_gcd == ntl_gcd || clpoly_gcd == -ntl_gcd);
        }
    }

    // ======== Derivative ========
    CLPOLY_TEST("crosscheck_ntl_derivative");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(7, 5, {-100, 100});

        // CLPoly derivative
        auto clpoly_deriv = derivative(f);

        // NTL derivative
        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
        NTL::ZZX nd = crosscheck::ntl_diff(nf);
        auto ntl_deriv = crosscheck::ntl_to_clpoly_upoly(nd);

        CLPOLY_ASSERT_EQ(clpoly_deriv, ntl_deriv);
    }

    // ======== Point evaluation ========
    CLPOLY_TEST("crosscheck_ntl_eval");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(6, 5, {-50, 50});
        ZZ point(trial * 3 - 5);  // evaluation points: -5, -2, 1, 4, 7

        // CLPoly evaluation
        ZZ clpoly_val = assign(f, point);

        // NTL evaluation
        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
        NTL::ZZ npoint = crosscheck::clpoly_zz_to_ntl(point);
        NTL::ZZ nval = crosscheck::ntl_eval(nf, npoint);
        ZZ ntl_val = crosscheck::ntl_zz_to_clpoly(nval);

        CLPOLY_ASSERT_EQ(clpoly_val, ntl_val);
    }

    // ======== High power: pow(f, 20) with binomials ========
    CLPOLY_TEST("crosscheck_ntl_pow20");
    for (int trial = 0; trial < 3; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        // Binomial: 2 terms, deg 1 -> pow 20 gives 21 terms
        auto f = random_upolynomial<ZZ>(1, 2, {-5, 5});
        if (f.empty()) continue;

        auto clpoly_result = pow(f, 20);

        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
        // NTL has no power() for ZZX; compute via repeated squaring
        NTL::ZZX nr;
        NTL::SetCoeff(nr, 0, NTL::ZZ(1));  // nr = 1
        for (int i = 0; i < 20; ++i) nr = nr * nf;
        auto ntl_result = crosscheck::ntl_to_clpoly_upoly(nr);

        CLPOLY_ASSERT_EQ(clpoly_result, ntl_result);
    }

    // ======== Large coefficients: manually constructed ========
    CLPOLY_TEST("crosscheck_ntl_large_coeff");
    {
        ZZ big1("999999999999999999999");      // ~10^21
        ZZ big2("-123456789012345678901234");
        ZZ big3("777777777777777777777777777");

        // f = big1*t^5 + big2*t^2 + big3
        upolynomial_ZZ f({
            {umonomial(5), big1},
            {umonomial(2), big2},
            {umonomial(0), big3}
        });
        // g = big3*t^3 + big1*t + big2
        upolynomial_ZZ g({
            {umonomial(3), big3},
            {umonomial(1), big1},
            {umonomial(0), big2}
        });

        // Test add
        {
            auto clpoly_result = f + g;
            NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
            NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
            NTL::ZZX nr = nf + ng;
            auto ntl_result = crosscheck::ntl_to_clpoly_upoly(nr);
            CLPOLY_ASSERT_EQ(clpoly_result, ntl_result);
        }
        // Test mul
        {
            auto clpoly_result = f * g;
            NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
            NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
            NTL::ZZX nr = nf * ng;
            auto ntl_result = crosscheck::ntl_to_clpoly_upoly(nr);
            CLPOLY_ASSERT_EQ(clpoly_result, ntl_result);
        }
        // Test roundtrip with product
        {
            auto product = f * g;
            NTL::ZZX np = crosscheck::clpoly_upoly_to_ntl(product);
            auto back = crosscheck::ntl_to_clpoly_upoly(np);
            CLPOLY_ASSERT_EQ(product, back);
        }
    }

    // ======== Large coefficients: random with scaling ========
    CLPOLY_TEST("crosscheck_ntl_large_coeff_random");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(6, 4, {-50, 50});
        auto g = random_upolynomial<ZZ>(6, 4, {-50, 50});
        ZZ scale = pow(ZZ(10), 18);
        auto f_big = f * scale;
        auto g_big = g * scale;

        auto clpoly_result = f_big * g_big;

        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f_big);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g_big);
        NTL::ZZX nr = nf * ng;
        auto ntl_result = crosscheck::ntl_to_clpoly_upoly(nr);

        CLPOLY_ASSERT_EQ(clpoly_result, ntl_result);
    }

    // ======== High degree: degree 80, sparse ========
    CLPOLY_TEST("crosscheck_ntl_high_degree_mul");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(80, 10, {-50, 50});
        auto g = random_upolynomial<ZZ>(80, 10, {-50, 50});

        auto clpoly_result = f * g;

        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
        NTL::ZZX nr = nf * ng;
        auto ntl_result = crosscheck::ntl_to_clpoly_upoly(nr);

        CLPOLY_ASSERT_EQ(clpoly_result, ntl_result);
    }

    // ======== High degree: derivative ========
    CLPOLY_TEST("crosscheck_ntl_high_degree_derivative");
    for (int trial = 0; trial < 3; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(100, 15, {-100, 100});

        auto clpoly_deriv = derivative(f);

        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
        NTL::ZZX nd = crosscheck::ntl_diff(nf);
        auto ntl_deriv = crosscheck::ntl_to_clpoly_upoly(nd);

        CLPOLY_ASSERT_EQ(clpoly_deriv, ntl_deriv);
    }

    // ======== High degree + large coeff: evaluation ========
    CLPOLY_TEST("crosscheck_ntl_high_degree_eval");
    for (int trial = 0; trial < 3; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(50, 8, {-50, 50});
        ZZ scale = pow(ZZ(10), 15);
        auto f_big = f * scale;
        ZZ point(trial * 2 - 2);  // -2, 0, 2

        ZZ clpoly_val = assign(f_big, point);

        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f_big);
        NTL::ZZ npoint = crosscheck::clpoly_zz_to_ntl(point);
        NTL::ZZ nval = crosscheck::ntl_eval(nf, npoint);
        ZZ ntl_val = crosscheck::ntl_zz_to_clpoly(nval);

        CLPOLY_ASSERT_EQ(clpoly_val, ntl_val);
    }

    return clpoly_test::test_summary();
}
