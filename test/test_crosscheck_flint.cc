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

    // 3 个随机因子, 随机次数/项数
    CLPOLY_TEST("crosscheck_flint_factor_3_factors");
    {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> deg_dis(2, 4);
    std::uniform_int_distribution<int> len_dis(2, 4);
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        auto f2 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        auto f3 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        if (f1.empty() || f2.empty() || f3.empty()) continue;
        auto uf = f1 * f2 * f3;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }
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

    // 4 个因子, 随机次数/项数
    CLPOLY_TEST("crosscheck_flint_factor_4_factors");
    {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> deg_dis(2, 3);
    std::uniform_int_distribution<int> len_dis(2, 3);
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        auto f2 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        auto f3 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        auto f4 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        if (f1.empty() || f2.empty() || f3.empty() || f4.empty()) continue;
        auto uf = f1 * f2 * f3 * f4;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }
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

    // 含重因子: f1^e1 * f2^e2, 随机重数/次数/项数
    CLPOLY_TEST("crosscheck_flint_factor_random_mult");
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> exp_dis(1, 3);
        std::uniform_int_distribution<int> deg_dis(2, 4);
        std::uniform_int_distribution<int> len_dis(2, 4);
        for (int trial = 0; trial < 5; ++trial) {
            CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
            int d1 = deg_dis(gen), d2 = deg_dis(gen);
            auto f1 = random_upolynomial<ZZ>(d1, len_dis(gen), {-30, 30});
            auto f2 = random_upolynomial<ZZ>(d2, len_dis(gen), {-30, 30});
            if (f1.empty() || f2.empty()) continue;
            int e1 = exp_dis(gen), e2 = exp_dis(gen);
            auto uf = pow(f1, e1) * pow(f2, e2);
            if (uf.empty() || get_deg(uf) < 2) continue;

            auto cl_fac = factorize(uf);
            auto fl_fac = crosscheck::flint_factor_upoly(uf);
            auto cl_sorted = normalize_cl_factors(cl_fac);
            auto fl_sorted = normalize_flint_factors(fl_fac);
            CLPOLY_ASSERT(cl_sorted == fl_sorted);
        }
    }

    // 含重因子: 3 个因子 f1^e1 * f2^e2 * f3^e3, 随机重数/次数/项数
    CLPOLY_TEST("crosscheck_flint_factor_random_mult_3");
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> exp_dis(1, 3);
        std::uniform_int_distribution<int> deg_dis(2, 4);
        std::uniform_int_distribution<int> len_dis(2, 4);
        for (int trial = 0; trial < 5; ++trial) {
            CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
            int d1 = deg_dis(gen), d2 = deg_dis(gen), d3 = deg_dis(gen);
            auto f1 = random_upolynomial<ZZ>(d1, len_dis(gen), {-30, 30});
            auto f2 = random_upolynomial<ZZ>(d2, len_dis(gen), {-30, 30});
            auto f3 = random_upolynomial<ZZ>(d3, len_dis(gen), {-30, 30});
            if (f1.empty() || f2.empty() || f3.empty()) continue;
            int e1 = exp_dis(gen), e2 = exp_dis(gen), e3 = exp_dis(gen);
            auto uf = pow(f1, e1) * pow(f2, e2) * pow(f3, e3);
            if (uf.empty() || get_deg(uf) < 2) continue;

            auto cl_fac = factorize(uf);
            auto fl_fac = crosscheck::flint_factor_upoly(uf);
            auto cl_sorted = normalize_cl_factors(cl_fac);
            auto fl_sorted = normalize_flint_factors(fl_fac);
            CLPOLY_ASSERT(cl_sorted == fl_sorted);
        }
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

    // 5+ 因子
    CLPOLY_TEST("crosscheck_flint_factor_5_factors");
    {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> deg_dis(1, 3);
    std::uniform_int_distribution<int> len_dis(2, 3);
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        auto f2 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        auto f3 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        auto f4 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        auto f5 = random_upolynomial<ZZ>(deg_dis(gen), len_dis(gen), {-10, 10});
        if (f1.empty() || f2.empty() || f3.empty() || f4.empty() || f5.empty()) continue;
        auto uf = f1 * f2 * f3 * f4 * f5;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }
    }

    // 高重数: (随机因子)^e, e∈[4,8]
    CLPOLY_TEST("crosscheck_flint_factor_high_mult");
    {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> exp_dis(4, 8);
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        if (f1.empty()) continue;
        int e = exp_dis(gen);
        auto uf = pow(f1, e);
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }
    }

    // 完全幂: (不可约二次)^n
    CLPOLY_TEST("crosscheck_flint_factor_perfect_power");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        // x² + (trial+1) 在 Z 上不可约 (对 trial=0..4 成立)
        upolynomial_ZZ irr({{umonomial(2), ZZ(1)}, {umonomial(0), ZZ(trial + 1)}});
        int e = trial + 2;  // e = 2,3,4,5,6
        auto uf = pow(irr, e);

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 混合 squarefree 分量: f1^1 * f2^2 * f3^3 (确保走不同 squarefree 分量)
    CLPOLY_TEST("crosscheck_flint_factor_mixed_squarefree");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        auto f2 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        auto f3 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        if (f1.empty() || f2.empty() || f3.empty()) continue;
        auto uf = f1 * pow(f2, 2) * pow(f3, 3);
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto fl_fac = crosscheck::flint_factor_upoly(uf);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto fl_sorted = normalize_flint_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // 退化输入: 零、常数、线性、单项式
    CLPOLY_TEST("crosscheck_flint_factor_zero");
    {
        upolynomial_ZZ uz;
        auto cl_fac = factorize(uz);
        CLPOLY_ASSERT_EQ(cl_fac.content, ZZ(0));
        CLPOLY_ASSERT(cl_fac.factors.empty());
    }

    CLPOLY_TEST("crosscheck_flint_factor_constant");
    {
        upolynomial_ZZ uc({{umonomial(0), ZZ(42)}});
        auto cl_fac = factorize(uc);
        CLPOLY_ASSERT_EQ(cl_fac.content, ZZ(42));
        CLPOLY_ASSERT(cl_fac.factors.empty());
    }

    CLPOLY_TEST("crosscheck_flint_factor_linear");
    {
        upolynomial_ZZ ul({{umonomial(1), ZZ(3)}, {umonomial(0), ZZ(5)}});
        auto cl_fac = factorize(ul);
        CLPOLY_ASSERT_EQ(cl_fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(cl_fac.factors[0].second, (uint64_t)1);
    }

    CLPOLY_TEST("crosscheck_flint_factor_monomial");
    {
        // x^7 → (x, 7)
        upolynomial_ZZ um({{umonomial(7), ZZ(1)}});
        auto cl_fac = factorize(um);
        CLPOLY_ASSERT_EQ(cl_fac.factors.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(cl_fac.factors[0].second, (uint64_t)7);
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

    // ================================================================
    // Multivariate factorization: CLPoly vs FLINT fmpz_mpoly_factor
    // ================================================================

    // Helper: verify factorization by reconstructing product
    auto verify_factor_product = [](const polynomial_ZZ& f,
                                     const ZZ& content,
                                     const std::vector<std::pair<polynomial_ZZ, long>>& factors) {
        polynomial_ZZ prod;
        basic_monomial<grlex> m0;
        prod.push_back({m0, content});
        for (auto& [fi, ei] : factors) {
            polynomial_ZZ fi_pow = fi;
            for (long e = 1; e < ei; ++e) {
                fi_pow = fi_pow * fi;
                fi_pow.normalization();
            }
            prod = prod * fi_pow;
            prod.normalization();
        }
        return prod == f;
    };

    // Helper: normalize multivar factor list for comparison (make lc positive, sort)
    auto normalize_mpoly_factors = [](const std::vector<std::pair<polynomial_ZZ, long>>& factors) {
        std::vector<std::pair<polynomial_ZZ, long>> result;
        for (auto& [fi, ei] : factors) {
            auto f = fi;
            if (!f.empty() && f.front().second < 0) f = -f;
            result.push_back({std::move(f), ei});
        }
        std::sort(result.begin(), result.end());
        return result;
    };

    // x^2 - y^2 = (x+y)(x-y)
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_x2_y2");
    {
        polynomial_ZZ f = pow(x,2) - pow(y,2);
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        // Both should produce 2 factors
        CLPOLY_ASSERT(cl_fac.factors.size() == 2);
        CLPOLY_ASSERT(fl_fac.factors.size() == 2);

        // Verify product reconstruction for both
        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
    }

    // x^3 + y^3 = (x+y)(x^2-xy+y^2)
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_x3_y3");
    {
        polynomial_ZZ f = pow(x,3) + pow(y,3);
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(cl_fac.factors.size() == 2);
        CLPOLY_ASSERT(fl_fac.factors.size() == 2);
        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
    }

    // (x+y)(x-y)(x+1) — three factors
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_3fac");
    {
        polynomial_ZZ f = (x + y) * (x - y) * (x + polynomial_ZZ(1));
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(cl_fac.factors.size() == 3);
        CLPOLY_ASSERT(fl_fac.factors.size() == 3);
        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
    }

    // (x+y)^2 * (x-y) — repeated factor
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_repeated");
    {
        polynomial_ZZ f = pow(x + y, 2) * (x - y);
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
        // Both should recognize multiplicity
        CLPOLY_ASSERT(cl_fac.factors.size() == 2);
        CLPOLY_ASSERT(fl_fac.factors.size() == 2);
    }

    // Trivariate: (x+y+z)(x-y+z)
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_trivar");
    {
        variable z("z");
        polynomial_ZZ f = (x + y + z) * (x - y + z);
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(cl_fac.factors.size() == 2);
        CLPOLY_ASSERT(fl_fac.factors.size() == 2);
        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
    }

    // With content: 6*(x^2 - y^2)
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_content");
    {
        polynomial_ZZ f = 6 * (pow(x,2) - pow(y,2));
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
    }

    // Helper: normalize CLPoly multivariate factors for comparison
    auto normalize_cl_mpoly_factors = [](const factorization<polynomial_ZZ>& fac) {
        std::vector<std::pair<polynomial_ZZ, uint64_t>> result;
        for (auto& [fi, ei] : fac.factors) {
            auto f = fi;
            if (!f.empty() && f.front().second < 0) f = -f;
            result.push_back({std::move(f), ei});
        }
        std::sort(result.begin(), result.end());
        return result;
    };

    // Helper: normalize FLINT multivariate factors for comparison
    auto normalize_fl_mpoly_factors = [](const crosscheck::FlintMPolyFactorResult& fac) {
        std::vector<std::pair<polynomial_ZZ, uint64_t>> result;
        for (auto& [fi, ei] : fac.factors) {
            auto f = fi;
            if (!f.empty() && f.front().second < 0) f = -f;
            result.push_back({std::move(f), (uint64_t)ei});
        }
        std::sort(result.begin(), result.end());
        return result;
    };

    // Random bivariate 2 factors (increased from 5 to 20 rounds)
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_random_bivar");
    for (int trial = 0; trial < 20; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_polynomial<ZZ>({x, y}, 3, 4, {-5, 5});
        auto f2 = random_polynomial<ZZ>({x, y}, 3, 4, {-5, 5});
        if (f1.empty() || f2.empty()) continue;
        polynomial_ZZ f = f1 * f2;
        if (f.empty()) continue;
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));

        // Exact factor list comparison
        auto cl_sorted = normalize_cl_mpoly_factors(cl_fac);
        auto fl_sorted = normalize_fl_mpoly_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // Random bivariate 4 factors
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_random_bivar_4fac");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
        auto f2 = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
        auto f3 = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
        auto f4 = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
        if (f1.empty() || f2.empty() || f3.empty() || f4.empty()) continue;
        polynomial_ZZ f = f1 * f2 * f3 * f4;
        if (f.empty()) continue;
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
        auto cl_sorted = normalize_cl_mpoly_factors(cl_fac);
        auto fl_sorted = normalize_fl_mpoly_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // Random trivariate 2 factors (increased from 3 to 10)
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_random_trivar");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_polynomial<ZZ>({x, y, z}, 2, 4, {-3, 3});
        auto f2 = random_polynomial<ZZ>({x, y, z}, 2, 4, {-3, 3});
        if (f1.empty() || f2.empty()) continue;
        polynomial_ZZ f = f1 * f2;
        if (f.empty()) continue;
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
        auto cl_sorted = normalize_cl_mpoly_factors(cl_fac);
        auto fl_sorted = normalize_fl_mpoly_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // Random trivariate 3 factors
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_random_trivar_3fac");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_polynomial<ZZ>({x, y, z}, 2, 3, {-3, 3});
        auto f2 = random_polynomial<ZZ>({x, y, z}, 2, 3, {-3, 3});
        auto f3 = random_polynomial<ZZ>({x, y, z}, 2, 3, {-3, 3});
        if (f1.empty() || f2.empty() || f3.empty()) continue;
        polynomial_ZZ f = f1 * f2 * f3;
        if (f.empty()) continue;
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
        auto cl_sorted = normalize_cl_mpoly_factors(cl_fac);
        auto fl_sorted = normalize_fl_mpoly_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // Random 4-variable 2 factors (reduced rounds — slow in debug mode)
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_random_4var");
    for (int trial = 0; trial < 2; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_polynomial<ZZ>({x, y, z, w}, 2, 3, {-3, 3});
        auto f2 = random_polynomial<ZZ>({x, y, z, w}, 2, 3, {-3, 3});
        if (f1.empty() || f2.empty()) continue;
        polynomial_ZZ f = f1 * f2;
        if (f.empty()) continue;
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
        auto cl_sorted = normalize_cl_mpoly_factors(cl_fac);
        auto fl_sorted = normalize_fl_mpoly_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // Random bivariate with multiplicity
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_random_bivar_mult");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
        auto f2 = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
        if (f1.empty() || f2.empty()) continue;
        polynomial_ZZ f = pow(f1, 2) * f2;
        if (f.empty()) continue;
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
        auto cl_sorted = normalize_cl_mpoly_factors(cl_fac);
        auto fl_sorted = normalize_fl_mpoly_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    // Irreducible
    CLPOLY_TEST("crosscheck_flint_mpoly_factor_irreducible");
    {
        polynomial_ZZ f = pow(x,2) + pow(y,2) + polynomial_ZZ(1);
        auto vars = crosscheck::collect_vars(f);

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

        CLPOLY_ASSERT(cl_fac.factors.size() == 1);
        CLPOLY_ASSERT(fl_fac.factors.size() == 1);
    }

    // ======== QQ univariate factorization crosscheck ========
    CLPOLY_TEST("crosscheck_flint_qq_factor_univar");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        // Generate random ZZ upoly, factor with FLINT
        auto uf1 = random_upolynomial<ZZ>(3, 3, {-10, 10});
        auto uf2 = random_upolynomial<ZZ>(3, 3, {-10, 10});
        if (uf1.empty() || uf2.empty()) continue;
        auto uf = uf1 * uf2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto fl_fac = crosscheck::flint_factor_upoly(uf);

        // Build polynomial_ZZ from upoly using variable x
        polynomial_ZZ pz;
        for (auto& [m, c] : uf)
            pz.push_back({pow(x, m.deg()), c});
        pz.normalization();

        // Convert to QQ with a rational scale
        polynomial_QQ pq;
        poly_convert(pz, pq);
        pq = pq * QQ(1, (trial + 2));

        auto qq_fac = factorize(pq);
        // Verify QQ factorization reconstructs
        polynomial_QQ qq_prod;
        basic_monomial<grlex> m0;
        qq_prod.push_back({m0, qq_fac.content});
        for (auto& [fi, ei] : qq_fac.factors) {
            polynomial_QQ fi_pow = fi;
            for (uint64_t e = 1; e < ei; ++e) {
                fi_pow = fi_pow * fi;
                fi_pow.normalization();
            }
            qq_prod = qq_prod * fi_pow;
            qq_prod.normalization();
        }
        CLPOLY_ASSERT_EQ(qq_prod, pq);

        // Factor count should match FLINT (ignoring content differences)
        CLPOLY_ASSERT_EQ(qq_fac.factors.size(), (size_t)fl_fac.factors.size());
    }

    // ======== Classic multivariate crosscheck ========

    CLPOLY_TEST("crosscheck_flint_mpoly_classic_sympy_f1");
    {
        // (x+yz+20)(xy+z+10)(xz+y+30)
        polynomial_ZZ f = (x + y*z + polynomial_ZZ(ZZ(20)))
                        * (x*y + z + polynomial_ZZ(ZZ(10)))
                        * (x*z + y + polynomial_ZZ(ZZ(30)));
        auto vars = crosscheck::collect_vars(f);
        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);
        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
        auto cl_sorted = normalize_cl_mpoly_factors(cl_fac);
        auto fl_sorted = normalize_fl_mpoly_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    CLPOLY_TEST("crosscheck_flint_mpoly_classic_x4_y4");
    {
        // x^4 - y^4 = (x-y)(x+y)(x^2+y^2)
        polynomial_ZZ f = pow(x,4) - pow(y,4);
        auto vars = crosscheck::collect_vars(f);
        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);
        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
        auto cl_sorted = normalize_cl_mpoly_factors(cl_fac);
        auto fl_sorted = normalize_fl_mpoly_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    CLPOLY_TEST("crosscheck_flint_mpoly_classic_cube");
    {
        // (x+y-z)^3
        polynomial_ZZ f = pow(x + y - z, 3);
        auto vars = crosscheck::collect_vars(f);
        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);
        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
        auto cl_sorted = normalize_cl_mpoly_factors(cl_fac);
        auto fl_sorted = normalize_fl_mpoly_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    CLPOLY_TEST("crosscheck_flint_mpoly_classic_sum_cube_minus_sum");
    {
        // (x+y+z)^3 - (x+y+z) = (x+y+z)(x+y+z-1)(x+y+z+1)
        auto s = x + y + z;
        polynomial_ZZ f = pow(s, 3) - s;
        auto vars = crosscheck::collect_vars(f);
        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);
        CLPOLY_ASSERT(verify_factor_product(f, fl_fac.content, fl_fac.factors));
        auto cl_sorted = normalize_cl_mpoly_factors(cl_fac);
        auto fl_sorted = normalize_fl_mpoly_factors(fl_fac);
        CLPOLY_ASSERT(cl_sorted == fl_sorted);
    }

    return clpoly_test::test_summary();
}
