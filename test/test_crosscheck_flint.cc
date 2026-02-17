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
    variable x("x"), y("y"), z("z"), w("w");
    std::vector<variable> vars2 = {x, y};
    std::vector<variable> vars3 = {x, y, z};
    std::vector<variable> vars4 = {x, y, z, w};

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

    // ======== Exact Division: f*g / g == f ========
    CLPOLY_TEST("crosscheck_flint_exact_div");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto& vars = (trial < 5) ? vars2 : vars3;
        auto f = random_polynomial<ZZ>(vars, 3, 4, {-50, 50});
        auto g = random_polynomial<ZZ>(vars, 2, 3, {-50, 50});
        if (f.empty() || g.empty()) continue;
        auto product = f * g;

        auto clpoly_quot = product / g;

        auto all_vars = crosscheck::collect_vars(product, g);
        auto fp = crosscheck::clpoly_to_flint(product, all_vars);
        auto fg = crosscheck::clpoly_to_flint(g, all_vars);
        crosscheck::FlintPoly fq(all_vars);
        int ok = fmpz_mpoly_divides(fq.poly, fp.poly, fg.poly, fq.ctx);
        CLPOLY_ASSERT(ok);
        auto flint_quot = crosscheck::flint_to_clpoly(fq);

        CLPOLY_ASSERT_EQ(clpoly_quot, f);
        CLPOLY_ASSERT_EQ(flint_quot, f);
    }

    // ======== Exact Division: large coefficients ========
    CLPOLY_TEST("crosscheck_flint_exact_div_large_coeff");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial<ZZ>(vars2, 3, 4, {-50, 50});
        auto g = random_polynomial<ZZ>(vars2, 2, 3, {-50, 50});
        if (f.empty() || g.empty()) continue;
        ZZ scale = pow(ZZ(10), 18);
        auto f_big = f * scale;
        auto product = f_big * g;

        auto clpoly_quot = product / g;

        auto all_vars = crosscheck::collect_vars(product, g);
        auto fp = crosscheck::clpoly_to_flint(product, all_vars);
        auto fg = crosscheck::clpoly_to_flint(g, all_vars);
        crosscheck::FlintPoly fq(all_vars);
        int ok = fmpz_mpoly_divides(fq.poly, fp.poly, fg.poly, fq.ctx);
        CLPOLY_ASSERT(ok);
        auto flint_quot = crosscheck::flint_to_clpoly(fq);

        CLPOLY_ASSERT_EQ(clpoly_quot, f_big);
        CLPOLY_ASSERT_EQ(flint_quot, f_big);
    }

    // ======== 4-variable multiplication ========
    CLPOLY_TEST("crosscheck_flint_4var_mul");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial<ZZ>(vars4, 3, 6, {-30, 30});
        auto g = random_polynomial<ZZ>(vars4, 3, 6, {-30, 30});

        auto clpoly_result = f * g;

        auto all_vars = crosscheck::collect_vars(f, g);
        auto ff = crosscheck::clpoly_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_to_flint(g, all_vars);
        crosscheck::FlintPoly fr(all_vars);
        fmpz_mpoly_mul(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== 4-variable large coeff ========
    CLPOLY_TEST("crosscheck_flint_4var_large_coeff");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial<ZZ>(vars4, 2, 5, {-30, 30});
        auto g = random_polynomial<ZZ>(vars4, 2, 5, {-30, 30});
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

    // ======== Dense polynomial multiplication ========
    CLPOLY_TEST("crosscheck_flint_dense_mul");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        // Dense: 2 vars, deg=3 => C(5,2)=10 monomials
        auto f = random_polynomial<ZZ>(vars2, 3, 10, {-30, 30});
        auto g = random_polynomial<ZZ>(vars2, 3, 10, {-30, 30});

        auto clpoly_result = f * g;

        auto all_vars = crosscheck::collect_vars(f, g);
        auto ff = crosscheck::clpoly_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_to_flint(g, all_vars);
        crosscheck::FlintPoly fr(all_vars);
        fmpz_mpoly_mul(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== Large coefficients: GCD ========
    CLPOLY_TEST("crosscheck_flint_large_coeff_gcd");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto h = random_polynomial<ZZ>(vars2, 2, 3, {-10, 10});
        auto a = random_polynomial<ZZ>(vars2, 2, 3, {-10, 10});
        auto b = random_polynomial<ZZ>(vars2, 2, 3, {-10, 10});
        if (h.empty() || a.empty() || b.empty()) continue;
        ZZ scale = pow(ZZ(10), 15);
        auto f = (h * a) * scale;
        auto g = (h * b) * scale;

        auto clpoly_gcd = gcd(f, g);

        auto all_vars = crosscheck::collect_vars(f, g);
        auto ff = crosscheck::clpoly_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_to_flint(g, all_vars);
        crosscheck::FlintPoly fgcd(all_vars);
        int ok = fmpz_mpoly_gcd(fgcd.poly, ff.poly, fg.poly, fgcd.ctx);
        CLPOLY_ASSERT(ok);
        auto flint_gcd_poly = crosscheck::flint_to_clpoly(fgcd);

        CLPOLY_ASSERT_EQ(degree(clpoly_gcd), degree(flint_gcd_poly));
        if (!clpoly_gcd.empty() && !flint_gcd_poly.empty()) {
            CLPOLY_ASSERT(divides(clpoly_gcd, flint_gcd_poly)
                       || divides(flint_gcd_poly, clpoly_gcd));
        }
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

    // ======== QQ: power ========
    CLPOLY_TEST("crosscheck_flint_qq_pow");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial_QQ(vars2, 2, 3, {-10, 10}, 8);
        unsigned long k = (trial % 2 == 0) ? 2 : 3;

        auto clpoly_result = pow(f, k);

        auto all_vars = crosscheck::collect_vars_qq(f);
        auto ff = crosscheck::clpoly_qq_to_flint(f, all_vars);
        crosscheck::FlintQPoly fr(all_vars);
        fmpq_mpoly_pow_ui(fr.poly, ff.poly, k, fr.ctx);
        auto flint_result = crosscheck::flint_qq_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== QQ: exact division ========
    CLPOLY_TEST("crosscheck_flint_qq_exact_div");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial_QQ(vars2, 3, 4, {-30, 30}, 8);
        auto g = random_polynomial_QQ(vars2, 2, 3, {-20, 20}, 6);
        if (f.empty() || g.empty()) continue;
        auto product = f * g;

        auto clpoly_quot = product / g;

        auto all_vars = crosscheck::collect_vars_qq(product, g);
        auto fp = crosscheck::clpoly_qq_to_flint(product, all_vars);
        auto fg = crosscheck::clpoly_qq_to_flint(g, all_vars);
        crosscheck::FlintQPoly fq(all_vars);
        fmpq_mpoly_div(fq.poly, fp.poly, fg.poly, fq.ctx);
        auto flint_quot = crosscheck::flint_qq_to_clpoly(fq);

        CLPOLY_ASSERT_EQ(clpoly_quot, f);
        CLPOLY_ASSERT_EQ(flint_quot, f);
    }

    // ======== QQ: 4-variable multiplication ========
    CLPOLY_TEST("crosscheck_flint_qq_4var_mul");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial_QQ(vars4, 2, 5, {-20, 20}, 8);
        auto g = random_polynomial_QQ(vars4, 2, 5, {-20, 20}, 8);

        auto clpoly_result = f * g;

        auto all_vars = crosscheck::collect_vars_qq(f, g);
        auto ff = crosscheck::clpoly_qq_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_qq_to_flint(g, all_vars);
        crosscheck::FlintQPoly fr(all_vars);
        fmpq_mpoly_mul(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_qq_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== QQ: large denominators ========
    CLPOLY_TEST("crosscheck_flint_qq_large_denom");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_polynomial_QQ(vars2, 3, 4, {-50, 50}, 100);
        auto g = random_polynomial_QQ(vars2, 3, 4, {-50, 50}, 100);

        auto clpoly_result = f * g;

        auto all_vars = crosscheck::collect_vars_qq(f, g);
        auto ff = crosscheck::clpoly_qq_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_qq_to_flint(g, all_vars);
        crosscheck::FlintQPoly fr(all_vars);
        fmpq_mpoly_mul(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_result = crosscheck::flint_qq_to_clpoly(fr);

        CLPOLY_ASSERT_EQ(clpoly_result, flint_result);
    }

    // ======== QQ: large numerators + denominators combined ========
    CLPOLY_TEST("crosscheck_flint_qq_large_num_denom");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        // Large numerators via scaling, large denominators via den_max
        auto f = random_polynomial_QQ(vars2, 2, 4, {-50, 50}, 50);
        auto g = random_polynomial_QQ(vars2, 2, 4, {-50, 50}, 50);

        // Also test addition (exercises common-denominator logic)
        auto clpoly_sum = f + g;
        auto all_vars = crosscheck::collect_vars_qq(f, g);
        auto ff = crosscheck::clpoly_qq_to_flint(f, all_vars);
        auto fg = crosscheck::clpoly_qq_to_flint(g, all_vars);
        crosscheck::FlintQPoly fr(all_vars);
        fmpq_mpoly_add(fr.poly, ff.poly, fg.poly, fr.ctx);
        auto flint_sum = crosscheck::flint_qq_to_clpoly(fr);
        CLPOLY_ASSERT_EQ(clpoly_sum, flint_sum);

        // And subtraction
        auto clpoly_diff = f - g;
        crosscheck::FlintQPoly fr2(all_vars);
        fmpq_mpoly_sub(fr2.poly, ff.poly, fg.poly, fr2.ctx);
        auto flint_diff = crosscheck::flint_qq_to_clpoly(fr2);
        CLPOLY_ASSERT_EQ(clpoly_diff, flint_diff);
    }

    // ======== Factorization: CLPoly vs FLINT (univariate) ========

    // Helper: 将 FLINT 因子列表转成排序后的 (upolynomial_ZZ, uint64_t) 列表（lc > 0）
    auto normalize_flint_factors = [](const crosscheck::FlintFactorResult& fac) {
        std::vector<std::pair<upolynomial_ZZ, uint64_t>> result;
        for (auto& [fi, ei] : fac.factors) {
            auto f = fi;
            if (!f.empty() && f.front().second < 0) f = -f;
            result.push_back({std::move(f), (uint64_t)ei});
        }
        std::sort(result.begin(), result.end());
        return result;
    };
    auto normalize_cl_factors = [](const factorization<upolynomial_<ZZ>>& fac) {
        std::vector<std::pair<upolynomial_ZZ, uint64_t>> result;
        for (auto& [fi, ei] : fac.factors) {
            auto f = fi;
            if (!f.empty() && f.front().second < 0) f = -f;
            result.push_back({std::move(f), ei});
        }
        std::sort(result.begin(), result.end());
        return result;
    };

    // 固定多项式: x^6 - 1
    CLPOLY_TEST("crosscheck_flint_factor_x6m1");
    {
        upolynomial_ZZ uf({
            {umonomial(6), ZZ(1)},
            {umonomial(0), ZZ(-1)}
        });
        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);

        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 含重因子: (x+1)^2 * (x-2)^3
    CLPOLY_TEST("crosscheck_flint_factor_with_mult");
    {
        upolynomial_ZZ xp1({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});
        upolynomial_ZZ xm2({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(-2)}});
        upolynomial_ZZ uf = pow(xp1, 2) * pow(xm2, 3);

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);

        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 随机多项式: f = f1 * f2
    CLPOLY_TEST("crosscheck_flint_factor_random");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto uf1 = random_upolynomial<ZZ>(3, 3, {-10, 10});
        auto uf2 = random_upolynomial<ZZ>(3, 3, {-10, 10});
        if (uf1.empty() || uf2.empty()) continue;
        auto uf = uf1 * uf2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);

        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 不可约: x^4 + 1
    CLPOLY_TEST("crosscheck_flint_factor_irreducible");
    {
        upolynomial_ZZ uf({{umonomial(4), ZZ(1)}, {umonomial(0), ZZ(1)}});
        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);

        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 非首一: 6x^2 + 5x + 1 = (2x+1)(3x+1)
    CLPOLY_TEST("crosscheck_flint_factor_nonmonic");
    {
        upolynomial_ZZ uf({
            {umonomial(2), ZZ(6)},
            {umonomial(1), ZZ(5)},
            {umonomial(0), ZZ(1)}
        });
        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);

        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 随机高次: deg 5+4=9, 10 trials — 精确比较因子列表
    CLPOLY_TEST("crosscheck_flint_factor_random_high_deg");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(5, 4, {-20, 20});
        auto f2 = random_upolynomial<ZZ>(4, 3, {-20, 20});
        if (f1.empty() || f2.empty()) continue;
        auto uf = f1 * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 大系数因式分解: 因子 ×10^12
    CLPOLY_TEST("crosscheck_flint_factor_large_coeff");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(3, 3, {-50, 50});
        auto f2 = random_upolynomial<ZZ>(3, 3, {-50, 50});
        if (f1.empty() || f2.empty()) continue;
        ZZ scale = pow(ZZ(10), 12);
        auto uf = (f1 * scale) * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 3 个随机因子
    CLPOLY_TEST("crosscheck_flint_factor_3_factors");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        auto f2 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        auto f3 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        if (f1.empty() || f2.empty() || f3.empty()) continue;
        auto uf = f1 * f2 * f3;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 稠密因子
    CLPOLY_TEST("crosscheck_flint_factor_dense");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(4, 5, {-20, 20});
        auto f2 = random_upolynomial<ZZ>(3, 4, {-20, 20});
        if (f1.empty() || f2.empty()) continue;
        auto uf = f1 * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 非首一 + 大系数: 手动构造
    CLPOLY_TEST("crosscheck_flint_factor_nonmonic_large");
    {
        ZZ big("100000000000000000");  // 10^17
        upolynomial_ZZ f1({{umonomial(1), big}, {umonomial(0), ZZ(1)}});
        upolynomial_ZZ f2({{umonomial(2), ZZ(1)}, {umonomial(0), -big}});
        auto uf = f1 * f2;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 4 个因子
    CLPOLY_TEST("crosscheck_flint_factor_4_factors");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        auto f2 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        auto f3 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        auto f4 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        if (f1.empty() || f2.empty() || f3.empty() || f4.empty()) continue;
        auto uf = f1 * f2 * f3 * f4;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 高次: deg 15+ (7+8=15)
    CLPOLY_TEST("crosscheck_flint_factor_very_high_deg");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(7, 5, {-10, 10});
        auto f2 = random_upolynomial<ZZ>(8, 5, {-10, 10});
        if (f1.empty() || f2.empty()) continue;
        auto uf = f1 * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 含重因子: 随机 f1^2 * f2
    CLPOLY_TEST("crosscheck_flint_factor_random_mult");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        auto f2 = random_upolynomial<ZZ>(3, 3, {-10, 10});
        if (f1.empty() || f2.empty()) continue;
        auto uf = pow(f1, 2) * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 显式含因子 x: x * f1 * f2
    CLPOLY_TEST("crosscheck_flint_factor_with_x");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(3, 3, {-20, 20});
        auto f2 = random_upolynomial<ZZ>(3, 3, {-20, 20});
        if (f1.empty() || f2.empty()) continue;
        upolynomial_ZZ x_mono({{umonomial(1), ZZ(1)}});
        auto uf = x_mono * f1 * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 自然大系数: 10^20 级
    CLPOLY_TEST("crosscheck_flint_factor_huge_coeff");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        ZZ base = pow(ZZ(10), 20);
        ZZ a1 = base + ZZ(trial * 7 + 3);
        ZZ a0 = base - ZZ(trial * 13 + 1);
        // (a1*x + a0) * (x^2 + trial+1)
        upolynomial_ZZ f1({{umonomial(1), a1}, {umonomial(0), a0}});
        upolynomial_ZZ f2({{umonomial(2), ZZ(1)}, {umonomial(0), ZZ(trial + 1)}});
        auto uf = f1 * f2;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    return clpoly_test::test_summary();
}
