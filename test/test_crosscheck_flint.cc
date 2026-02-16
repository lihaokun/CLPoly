/**
 * @file test_crosscheck_flint.cc
 * @brief Cross-library correctness tests: CLPoly polynomial_ZZ vs FLINT fmpz_mpoly.
 *
 * Random multivariate polynomials are generated in CLPoly, converted to FLINT,
 * operations are performed in both libraries, and results are compared.
 */
#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include "crosscheck_flint.hh"

using namespace clpoly;

// Helper: check polynomial divisibility (a divides b)
static bool divides(const polynomial_ZZ& a, const polynomial_ZZ& b) {
    if (a.empty()) return b.empty();
    if (b.empty()) return true;
    auto q = b / a;
    return q * a == b;
}

int main() {
    variable x("x"), y("y"), z("z");
    std::vector<variable> vars2 = {x, y};
    std::vector<variable> vars3 = {x, y, z};

    // ======== Round-trip: CLPoly -> FLINT -> CLPoly ========
    CLPOLY_TEST("crosscheck_flint_roundtrip");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto& vars = (trial < 5) ? vars2 : vars3;
        auto f = random_polynomial<ZZ>(vars, 4, 5, {-50, 50});
        auto all_vars = crosscheck::collect_vars(f);
        auto fp = crosscheck::clpoly_to_flint(f, all_vars);
        auto f_back = crosscheck::flint_to_clpoly(fp);
        CLPOLY_ASSERT_EQ(f, f_back);
    }

    // ======== Addition: f + g ========
    CLPOLY_TEST("crosscheck_flint_add");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto& vars = (trial < 5) ? vars2 : vars3;
        auto f = random_polynomial<ZZ>(vars, 3, 5, {-50, 50});
        auto g = random_polynomial<ZZ>(vars, 3, 5, {-50, 50});

        // CLPoly result
        auto clpoly_result = f + g;

        // FLINT result
        auto all_vars = crosscheck::collect_vars(f, g);
        auto ff = crosscheck::clpoly_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_to_flint(g, all_vars);
        crosscheck::FlintPoly fr(all_vars);
        fmpz_mpoly_add(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== Subtraction: f - g ========
    CLPOLY_TEST("crosscheck_flint_sub");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto& vars = (trial < 5) ? vars2 : vars3;
        auto f = random_polynomial<ZZ>(vars, 3, 5, {-50, 50});
        auto g = random_polynomial<ZZ>(vars, 3, 5, {-50, 50});

        auto clpoly_result = f - g;

        auto all_vars = crosscheck::collect_vars(f, g);
        auto ff = crosscheck::clpoly_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_to_flint(g, all_vars);
        crosscheck::FlintPoly fr(all_vars);
        fmpz_mpoly_sub(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== Multiplication: f * g ========
    CLPOLY_TEST("crosscheck_flint_mul");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto& vars = (trial < 5) ? vars2 : vars3;
        auto f = random_polynomial<ZZ>(vars, 3, 4, {-50, 50});
        auto g = random_polynomial<ZZ>(vars, 3, 4, {-50, 50});

        auto clpoly_result = f * g;

        auto all_vars = crosscheck::collect_vars(f, g);
        auto ff = crosscheck::clpoly_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_to_flint(g, all_vars);
        crosscheck::FlintPoly fr(all_vars);
        fmpz_mpoly_mul(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== Power: pow(f, k) ========
    CLPOLY_TEST("crosscheck_flint_pow");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto& vars = (trial < 3) ? vars2 : vars3;
        auto f = random_polynomial<ZZ>(vars, 2, 3, {-10, 10});
        unsigned long k = (trial % 2 == 0) ? 2 : 3;

        auto clpoly_result = pow(f, k);

        auto all_vars = crosscheck::collect_vars(f);
        auto ff = crosscheck::clpoly_to_flint(f, all_vars);
        crosscheck::FlintPoly fr(all_vars);
        fmpz_mpoly_pow_ui(fr.poly, ff.poly, k, fr.ctx);
        auto flint_result = crosscheck::flint_to_clpoly(fr);

        // pow(f, 2) and pow(f, 3) crosscheck
        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== GCD (up to unit) ========
    CLPOLY_TEST("crosscheck_flint_gcd");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        // Construct f = h * a, g = h * b so gcd is at least h
        auto h = random_polynomial<ZZ>(vars2, 2, 3, {-10, 10});
        auto a = random_polynomial<ZZ>(vars2, 2, 3, {-10, 10});
        auto b = random_polynomial<ZZ>(vars2, 2, 3, {-10, 10});
        auto f = h * a;
        auto g = h * b;

        if (f.empty() || g.empty()) continue;

        // CLPoly gcd
        auto clpoly_gcd = gcd(f, g);

        // FLINT gcd
        auto all_vars = crosscheck::collect_vars(f, g);
        auto ff = crosscheck::clpoly_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_to_flint(g, all_vars);
        crosscheck::FlintPoly fgcd(all_vars);
        int ok = fmpz_mpoly_gcd(fgcd.poly, ff.poly, fg.poly, fgcd.ctx);
        CLPOLY_ASSERT(ok);

        auto flint_gcd_poly = crosscheck::flint_to_clpoly(fgcd);

        // GCD is unique up to unit (±1 over ZZ):
        // 1. Degrees should match
        CLPOLY_ASSERT_EQ(degree(clpoly_gcd), degree(flint_gcd_poly));
        // 2. Each should divide the other (up to sign)
        if (!clpoly_gcd.empty() && !flint_gcd_poly.empty()) {
            CLPOLY_ASSERT(divides(clpoly_gcd, flint_gcd_poly)
                       || divides(flint_gcd_poly, clpoly_gcd));
        }
    }

    // ======== High power: pow(f, 20) with binomials ========
    CLPOLY_TEST("crosscheck_flint_pow20");
    for (int trial = 0; trial < 3; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        // Small binomial: 2 terms, deg 1 -> pow 20 gives 21 terms, manageable
        auto f = random_polynomial<ZZ>(vars2, 1, 2, {-5, 5});
        if (f.empty()) continue;

        auto clpoly_result = pow(f, 20);

        auto all_vars = crosscheck::collect_vars(f);
        auto ff = crosscheck::clpoly_to_flint(f, all_vars);
        crosscheck::FlintPoly fr(all_vars);
        fmpz_mpoly_pow_ui(fr.poly, ff.poly, 20, fr.ctx);
        auto flint_result = crosscheck::flint_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== Large coefficients (overflow int64) ========
    CLPOLY_TEST("crosscheck_flint_large_coeff");
    {
        // Manually construct polynomials with huge coefficients
        ZZ big1("999999999999999999999");   // ~10^21
        ZZ big2("-123456789012345678901234");
        ZZ big3("777777777777777777777777777");

        // f = big1 * x^3*y + big2 * x*y^2 + big3
        polynomial_ZZ f({
            {pow(x,3)*y, big1},
            {x*pow(y,2), big2},
            {monomial(), big3}
        });
        // g = big3 * x^2 + big1 * y + big2
        polynomial_ZZ g({
            {pow(x,2), big3},
            {monomial(y), big1},
            {monomial(), big2}
        });

        auto all_vars = crosscheck::collect_vars(f, g);

        // Test add
        {
            auto clpoly_result = f + g;
            auto ff = crosscheck::clpoly_to_flint(f, all_vars);
            auto fg = crosscheck::clpoly_to_flint(g, all_vars);
            crosscheck::FlintPoly fr(all_vars);
            fmpz_mpoly_add(fr.poly, ff.poly, fg.poly, fr.ctx);
            auto flint_result = crosscheck::flint_to_clpoly(fr);
            CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
        }
        // Test mul
        {
            auto clpoly_result = f * g;
            auto ff = crosscheck::clpoly_to_flint(f, all_vars);
            auto fg = crosscheck::clpoly_to_flint(g, all_vars);
            crosscheck::FlintPoly fr(all_vars);
            fmpz_mpoly_mul(fr.poly, ff.poly, fg.poly, fr.ctx);
            auto flint_result = crosscheck::flint_to_clpoly(fr);
            CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
        }
        // Test roundtrip with product (very large coefficients after mul)
        {
            auto product = f * g;
            auto fp = crosscheck::clpoly_to_flint(product, all_vars);
            auto back = crosscheck::flint_to_clpoly(fp);
            CLPOLY_ASSERT_EQ(product, back);
        }
    }

    // ======== Large coefficients: random with scaling ========
    CLPOLY_TEST("crosscheck_flint_large_coeff_random");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial<ZZ>(vars2, 3, 4, {-50, 50});
        auto g = random_polynomial<ZZ>(vars2, 3, 4, {-50, 50});
        // Scale by 10^18 to push into big integer territory
        ZZ scale = pow(ZZ(10), 18);
        auto f_big = f * scale;
        auto g_big = g * scale;

        auto clpoly_result = f_big * g_big;

        auto all_vars = crosscheck::collect_vars(f_big, g_big);
        auto ff = crosscheck::clpoly_to_flint(f_big, all_vars);
        auto fg = crosscheck::clpoly_to_flint(g_big, all_vars);
        crosscheck::FlintPoly fr(all_vars);
        fmpz_mpoly_mul(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== High degree: deg 20, sparse ========
    CLPOLY_TEST("crosscheck_flint_high_degree");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial<ZZ>(vars2, 20, 8, {-30, 30});
        auto g = random_polynomial<ZZ>(vars2, 20, 8, {-30, 30});

        // mul
        auto clpoly_result = f * g;

        auto all_vars = crosscheck::collect_vars(f, g);
        auto ff = crosscheck::clpoly_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_to_flint(g, all_vars);
        crosscheck::FlintPoly fr(all_vars);
        fmpz_mpoly_mul(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== QQ: round-trip ========
    CLPOLY_TEST("crosscheck_flint_qq_roundtrip");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial_QQ(vars2, 3, 4, {-50, 50}, 12);
        auto all_vars = crosscheck::collect_vars_qq(f);
        auto fp = crosscheck::clpoly_qq_to_flint(f, all_vars);
        auto f_back = crosscheck::flint_qq_to_clpoly(fp);
        CLPOLY_ASSERT_EQ(f, f_back);
    }

    // ======== QQ: addition ========
    CLPOLY_TEST("crosscheck_flint_qq_add");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial_QQ(vars2, 3, 4, {-50, 50}, 12);
        auto g = random_polynomial_QQ(vars2, 3, 4, {-50, 50}, 12);

        auto clpoly_result = f + g;

        auto all_vars = crosscheck::collect_vars_qq(f, g);
        auto ff = crosscheck::clpoly_qq_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_qq_to_flint(g, all_vars);
        crosscheck::FlintQPoly fr(all_vars);
        fmpq_mpoly_add(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_qq_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== QQ: subtraction ========
    CLPOLY_TEST("crosscheck_flint_qq_sub");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial_QQ(vars2, 3, 4, {-50, 50}, 12);
        auto g = random_polynomial_QQ(vars2, 3, 4, {-50, 50}, 12);

        auto clpoly_result = f - g;

        auto all_vars = crosscheck::collect_vars_qq(f, g);
        auto ff = crosscheck::clpoly_qq_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_qq_to_flint(g, all_vars);
        crosscheck::FlintQPoly fr(all_vars);
        fmpq_mpoly_sub(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_qq_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== QQ: multiplication ========
    CLPOLY_TEST("crosscheck_flint_qq_mul");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial_QQ(vars2, 3, 4, {-30, 30}, 8);
        auto g = random_polynomial_QQ(vars2, 3, 4, {-30, 30}, 8);

        auto clpoly_result = f * g;

        auto all_vars = crosscheck::collect_vars_qq(f, g);
        auto ff = crosscheck::clpoly_qq_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_qq_to_flint(g, all_vars);
        crosscheck::FlintQPoly fr(all_vars);
        fmpq_mpoly_mul(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_qq_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    return clpoly_test::test_summary();
}
