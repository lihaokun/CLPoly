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

    // ======== Exact Division: f*g / g == f ========
    CLPOLY_TEST("crosscheck_ntl_exact_div");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(5, 4, {-50, 50});
        auto g = random_upolynomial<ZZ>(3, 3, {-50, 50});
        if (f.empty() || g.empty()) continue;
        auto product = f * g;

        auto clpoly_quot = product / g;

        NTL::ZZX np = crosscheck::clpoly_upoly_to_ntl(product);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
        NTL::ZZX nq;
        NTL::divide(nq, np, ng);
        auto ntl_quot = crosscheck::ntl_to_clpoly_upoly(nq);

        CLPOLY_ASSERT_EQ(clpoly_quot, f);
        CLPOLY_ASSERT_EQ(ntl_quot, f);
    }

    // ======== Exact Division: large coefficients ========
    CLPOLY_TEST("crosscheck_ntl_exact_div_large_coeff");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(4, 4, {-50, 50});
        auto g = random_upolynomial<ZZ>(3, 3, {-50, 50});
        if (f.empty() || g.empty()) continue;
        ZZ scale = pow(ZZ(10), 18);
        auto f_big = f * scale;
        auto product = f_big * g;

        auto clpoly_quot = product / g;

        NTL::ZZX np = crosscheck::clpoly_upoly_to_ntl(product);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
        NTL::ZZX nq;
        NTL::divide(nq, np, ng);
        auto ntl_quot = crosscheck::ntl_to_clpoly_upoly(nq);

        CLPOLY_ASSERT_EQ(clpoly_quot, f_big);
        CLPOLY_ASSERT_EQ(ntl_quot, f_big);
    }

    // ======== Dense polynomial multiplication ========
    CLPOLY_TEST("crosscheck_ntl_dense_mul");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        // Dense: len close to deg+1 (all coefficients present)
        auto f = random_upolynomial<ZZ>(8, 9, {-50, 50});
        auto g = random_upolynomial<ZZ>(8, 9, {-50, 50});

        auto clpoly_result = f * g;

        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
        NTL::ZZX nr = nf * ng;
        auto ntl_result = crosscheck::ntl_to_clpoly_upoly(nr);

        CLPOLY_ASSERT_EQ(clpoly_result, ntl_result);
    }

    // ======== Dense polynomial with large coefficients ========
    CLPOLY_TEST("crosscheck_ntl_dense_large_coeff");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(6, 7, {-50, 50});
        auto g = random_upolynomial<ZZ>(6, 7, {-50, 50});
        ZZ scale = pow(ZZ(10), 20);
        auto f_big = f * scale;
        auto g_big = g * scale;

        auto clpoly_result = f_big * g_big;

        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f_big);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g_big);
        NTL::ZZX nr = nf * ng;
        auto ntl_result = crosscheck::ntl_to_clpoly_upoly(nr);

        CLPOLY_ASSERT_EQ(clpoly_result, ntl_result);
    }

    // ======== Large coefficients: GCD ========
    CLPOLY_TEST("crosscheck_ntl_large_coeff_gcd");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto h = random_upolynomial<ZZ>(2, 3, {-20, 20});
        auto a = random_upolynomial<ZZ>(2, 3, {-20, 20});
        auto b = random_upolynomial<ZZ>(2, 3, {-20, 20});
        if (h.empty() || a.empty() || b.empty()) continue;
        ZZ scale = pow(ZZ(10), 15);
        auto f = (h * a) * scale;
        auto g = (h * b) * scale;

        auto clpoly_gcd = polynomial_GCD(f, g);

        NTL::ZZX nf = crosscheck::clpoly_upoly_to_ntl(f);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
        NTL::ZZX ngcd = NTL::GCD(nf, ng);
        auto ntl_gcd = crosscheck::ntl_to_clpoly_upoly(ngcd);

        CLPOLY_ASSERT_EQ(degree(clpoly_gcd), degree(ntl_gcd));
        if (!clpoly_gcd.empty() && !ntl_gcd.empty()) {
            CLPOLY_ASSERT(clpoly_gcd == ntl_gcd || clpoly_gcd == -ntl_gcd);
        }
    }

    // ======== Large coefficients: exact division + GCD combined ========
    CLPOLY_TEST("crosscheck_ntl_large_coeff_div_gcd");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f = random_upolynomial<ZZ>(4, 4, {-30, 30});
        auto g = random_upolynomial<ZZ>(3, 3, {-30, 30});
        if (f.empty() || g.empty()) continue;
        ZZ scale = pow(ZZ(10), 20);
        auto f_big = f * scale;
        auto product = f_big * g;

        // Division
        auto clpoly_quot = product / g;
        NTL::ZZX np = crosscheck::clpoly_upoly_to_ntl(product);
        NTL::ZZX ng = crosscheck::clpoly_upoly_to_ntl(g);
        NTL::ZZX nq;
        NTL::divide(nq, np, ng);
        auto ntl_quot = crosscheck::ntl_to_clpoly_upoly(nq);
        CLPOLY_ASSERT_EQ(clpoly_quot, ntl_quot);

        // GCD of product and g should divide g
        auto clpoly_gcd = polynomial_GCD(product, g);
        NTL::ZZX ngcd = NTL::GCD(np, ng);
        auto ntl_gcd = crosscheck::ntl_to_clpoly_upoly(ngcd);
        CLPOLY_ASSERT_EQ(degree(clpoly_gcd), degree(ntl_gcd));
    }

    // ======== Factorization: CLPoly vs NTL ========

    // Helper: NTL 因子转成排序后的 (upolynomial_ZZ, uint64_t) 列表（lc > 0）
    auto normalize_ntl_factors = [](const crosscheck::NTLFactorResult& fac) {
        std::vector<std::pair<upolynomial_ZZ, uint64_t>> result;
        for (auto& [fi, ei] : fac.factors) {
            auto f = crosscheck::ntl_to_clpoly_upoly(fi);
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

    // 固定多项式: x^6 - 1 = (x-1)(x+1)(x^2-x+1)(x^2+x+1)
    CLPOLY_TEST("crosscheck_ntl_factor_x6m1");
    {
        upolynomial_ZZ uf({
            {umonomial(6), ZZ(1)},
            {umonomial(0), ZZ(-1)}
        });
        auto cl_fac = factorize(uf);

        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);

        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 含重因子: (x+1)^2 * (x-2)^3
    CLPOLY_TEST("crosscheck_ntl_factor_with_mult");
    {
        upolynomial_ZZ xp1({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});
        upolynomial_ZZ xm2({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(-2)}});
        upolynomial_ZZ uf = pow(xp1, 2) * pow(xm2, 3);

        auto cl_fac = factorize(uf);

        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);

        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 随机多项式: 构造 f = f1 * f2, 排序后逐因子比较
    CLPOLY_TEST("crosscheck_ntl_factor_random");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(3, 3, {-10, 10});
        auto f2 = random_upolynomial<ZZ>(3, 3, {-10, 10});
        if (f1.empty() || f2.empty()) continue;
        auto uf = f1 * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);

        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);

        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 不可约: x^4 + 1
    CLPOLY_TEST("crosscheck_ntl_factor_irreducible");
    {
        upolynomial_ZZ uf({{umonomial(4), ZZ(1)}, {umonomial(0), ZZ(1)}});
        auto cl_fac = factorize(uf);

        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);

        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 随机多项式: 高次因子 (deg 5+4=9)
    CLPOLY_TEST("crosscheck_ntl_factor_random_high_deg");
    for (int trial = 0; trial < 10; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(5, 4, {-20, 20});
        auto f2 = random_upolynomial<ZZ>(4, 3, {-20, 20});
        if (f1.empty() || f2.empty()) continue;
        auto uf = f1 * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 大系数因式分解: 因子 ×10^12 后乘积
    CLPOLY_TEST("crosscheck_ntl_factor_large_coeff");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(3, 3, {-50, 50});
        auto f2 = random_upolynomial<ZZ>(3, 3, {-50, 50});
        if (f1.empty() || f2.empty()) continue;
        ZZ scale = pow(ZZ(10), 12);
        auto uf = (f1 * scale) * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 3 个随机因子
    CLPOLY_TEST("crosscheck_ntl_factor_3_factors");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        auto f2 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        auto f3 = random_upolynomial<ZZ>(2, 2, {-10, 10});
        if (f1.empty() || f2.empty() || f3.empty()) continue;
        auto uf = f1 * f2 * f3;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 稠密因子: deg=4 len=5 (几乎全项)
    CLPOLY_TEST("crosscheck_ntl_factor_dense");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(4, 5, {-20, 20});
        auto f2 = random_upolynomial<ZZ>(3, 4, {-20, 20});
        if (f1.empty() || f2.empty()) continue;
        auto uf = f1 * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 非首一 + 大系数: 手动构造
    CLPOLY_TEST("crosscheck_ntl_factor_nonmonic_large");
    {
        ZZ big("100000000000000000");  // 10^17
        // (big*x + 1)(x^2 - big)
        upolynomial_ZZ f1({{umonomial(1), big}, {umonomial(0), ZZ(1)}});
        upolynomial_ZZ f2({{umonomial(2), ZZ(1)}, {umonomial(0), -big}});
        auto uf = f1 * f2;

        auto cl_fac = factorize(uf);
        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 4 个因子
    CLPOLY_TEST("crosscheck_ntl_factor_4_factors");
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
        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 高次: deg 15+ (7+8=15)
    CLPOLY_TEST("crosscheck_ntl_factor_very_high_deg");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(7, 5, {-10, 10});
        auto f2 = random_upolynomial<ZZ>(8, 5, {-10, 10});
        if (f1.empty() || f2.empty()) continue;
        auto uf = f1 * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 含重因子: f1^e1 * f2^e2, 随机重数
    CLPOLY_TEST("crosscheck_ntl_factor_random_mult");
    {
        int exps[][2] = {{2,1}, {3,1}, {2,2}, {1,3}, {2,3}};
        for (int trial = 0; trial < 5; ++trial) {
            CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
            auto f1 = random_upolynomial<ZZ>(2, 2, {-10, 10});
            auto f2 = random_upolynomial<ZZ>(2, 2, {-10, 10});
            if (f1.empty() || f2.empty()) continue;
            int e1 = exps[trial][0], e2 = exps[trial][1];
            auto uf = pow(f1, e1) * pow(f2, e2);
            if (uf.empty() || get_deg(uf) < 2) continue;

            auto cl_fac = factorize(uf);
            auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
            auto ntl_fac = crosscheck::ntl_factor(ntl_f);
            auto cl_sorted = normalize_cl_factors(cl_fac);
            auto ntl_sorted = normalize_ntl_factors(ntl_fac);
            CLPOLY_ASSERT(cl_sorted == ntl_sorted);
        }
    }

    // 含重因子: 3 个因子 f1^e1 * f2^e2 * f3^e3
    CLPOLY_TEST("crosscheck_ntl_factor_random_mult_3");
    {
        int exps[][3] = {{2,1,1}, {1,2,1}, {1,1,2}, {2,2,1}, {2,1,3}};
        for (int trial = 0; trial < 5; ++trial) {
            CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
            auto f1 = random_upolynomial<ZZ>(2, 2, {-10, 10});
            auto f2 = random_upolynomial<ZZ>(2, 2, {-10, 10});
            auto f3 = random_upolynomial<ZZ>(2, 2, {-10, 10});
            if (f1.empty() || f2.empty() || f3.empty()) continue;
            int e1 = exps[trial][0], e2 = exps[trial][1], e3 = exps[trial][2];
            auto uf = pow(f1, e1) * pow(f2, e2) * pow(f3, e3);
            if (uf.empty() || get_deg(uf) < 2) continue;

            auto cl_fac = factorize(uf);
            auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
            auto ntl_fac = crosscheck::ntl_factor(ntl_f);
            auto cl_sorted = normalize_cl_factors(cl_fac);
            auto ntl_sorted = normalize_ntl_factors(ntl_fac);
            CLPOLY_ASSERT(cl_sorted == ntl_sorted);
        }
    }

    // 显式含因子 x: x * f1 * f2
    CLPOLY_TEST("crosscheck_ntl_factor_with_x");
    for (int trial = 0; trial < 5; ++trial) {
        CLPOLY_TEST_SECTION("trial_" + std::to_string(trial));
        auto f1 = random_upolynomial<ZZ>(3, 3, {-20, 20});
        auto f2 = random_upolynomial<ZZ>(3, 3, {-20, 20});
        if (f1.empty() || f2.empty()) continue;
        upolynomial_ZZ x_mono({{umonomial(1), ZZ(1)}});
        auto uf = x_mono * f1 * f2;
        if (uf.empty() || get_deg(uf) < 2) continue;

        auto cl_fac = factorize(uf);
        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    // 自然大系数: 10^20 级
    CLPOLY_TEST("crosscheck_ntl_factor_huge_coeff");
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
        auto ntl_f = crosscheck::clpoly_upoly_to_ntl(uf);
        auto ntl_fac = crosscheck::ntl_factor(ntl_f);
        auto cl_sorted = normalize_cl_factors(cl_fac);
        auto ntl_sorted = normalize_ntl_factors(ntl_fac);
        CLPOLY_ASSERT(cl_sorted == ntl_sorted);
    }

    return clpoly_test::test_summary();
}
