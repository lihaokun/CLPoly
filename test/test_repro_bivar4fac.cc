/**
 * Mass random test for multivar factorization (bivariate 4 factors).
 *
 * Usage: test_repro_bivar4fac [N_TRIALS] [-v]
 *   N_TRIALS: number of random trials (default 10)
 *   -v: verbose, print all polynomials
 */
#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include <iostream>
#include <random>
#include <cstring>
#include <cstdlib>

using namespace clpoly;

// Generate a random bivariate polynomial of degree d with coefficients in [-C, C]
polynomial_ZZ random_bivar(polynomial_ZZ px, polynomial_ZZ py, int d, int C, std::mt19937& rng) {
    std::uniform_int_distribution<int> coeff_dist(-C, C);
    polynomial_ZZ result;
    std::vector<polynomial_ZZ> xpow(d+1), ypow(d+1);
    xpow[0] = polynomial_ZZ({{{},(ZZ)1}}); ypow[0] = xpow[0];
    for (int k = 1; k <= d; k++) {
        xpow[k] = xpow[k-1] * px; xpow[k].normalization();
        ypow[k] = ypow[k-1] * py; ypow[k].normalization();
    }
    for (int i = 0; i <= d; i++) {
        for (int j = 0; j <= d - i; j++) {
            int c = coeff_dist(rng);
            if (c == 0) continue;
            auto term = polynomial_ZZ({{{},(ZZ)c}}) * xpow[i] * ypow[j];
            term.normalization();
            result = result + term;
            result.normalization();
        }
    }
    return result;
}

int main(int argc, char* argv[]) {
    int N_TRIALS = 10;
    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0)
            verbose = true;
        else
            N_TRIALS = atoi(argv[i]);
    }

    variable x("x"), y("y");
    int pass = 0, fail = 0;

    for (int trial = 0; trial < N_TRIALS; trial++) {
        std::mt19937 rng(trial * 12345 + 9876);

        polynomial_ZZ px(x), py(y);
        polynomial_ZZ factors[4];
        for (int i = 0; i < 4; i++)
            factors[i] = random_bivar(px, py, 1 + (trial % 2), 5, rng);

        bool skip = false;
        for (int i = 0; i < 4; i++)
            if (factors[i].empty() || is_number(factors[i]) || get_deg(factors[i]) == 0)
                skip = true;
        if (skip) {
            if (verbose) std::cerr << "trial " << trial << " skip" << std::endl;
            continue;
        }

        polynomial_ZZ f = factors[0] * factors[1] * factors[2] * factors[3];

        if (verbose) {
            std::cerr << "trial " << trial << std::endl;
            for (int i = 0; i < 4; i++)
                std::cerr << "  input[" << i << "] = " << factors[i] << std::endl;
        }

        auto fac = factorize(f);

        polynomial_ZZ product({{{},(ZZ)fac.content}});
        for (auto& [fi, ei] : fac.factors)
            for (uint64_t e = 0; e < ei; e++) {
                product = product * fi;
                product.normalization();
            }

        bool ok = (product == f || product == -f);

        if (!ok) {
            std::cerr << "FAIL trial " << trial << std::endl;
            for (int i = 0; i < 4; i++)
                std::cerr << "  input[" << i << "] = " << factors[i] << std::endl;
            fail++;
        } else {
            if (verbose) std::cerr << "  PASS" << std::endl;
            pass++;
        }
    }

    CLPOLY_TEST("mass_random_bivar4fac");
    std::cout << "Trials: " << (pass + fail) << "  Pass: " << pass << "  Fail: " << fail << std::endl;
    CLPOLY_ASSERT_EQ(fail, 0);

    // Deterministic reproducer
    CLPOLY_TEST("repro_bivar4fac_trial80");
    {
        polynomial_ZZ one({{{},(ZZ)1}});
        polynomial_ZZ two({{{},(ZZ)2}});
        polynomial_ZZ three({{{},(ZZ)3}});
        polynomial_ZZ f1 = -two*x*y + y*y - one;
        polynomial_ZZ f2 = -three*x*y - two*y + one;
        polynomial_ZZ f3 = x*x + x*y + two*x;
        polynomial_ZZ f4 = three*x*y - three*y - two;
        polynomial_ZZ f = f1 * f2 * f3 * f4;
        auto fac = factorize(f);
        size_t nontrivial = 0;
        for (auto& [fi, ei] : fac.factors)
            if (get_deg(fi) > 0) ++nontrivial;
        CLPOLY_ASSERT_EQ(nontrivial, (size_t)5);
    }

    return clpoly_test::test_summary();
}
