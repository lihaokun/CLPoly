/**
 * Mass random test for multivar factorization (trivariate 3 factors).
 *
 * Usage: test_repro_trivar3fac [N_TRIALS] [-v]
 *   N_TRIALS: number of random trials (default 10)
 *   -v: verbose, print all polynomials
 */
#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include <iostream>
#include <cstring>

using namespace clpoly;

int main(int argc, char* argv[]) {
    int N_TRIALS = 10;
    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-v") == 0)
            verbose = true;
        else
            N_TRIALS = atoi(argv[i]);
    }

    variable x("x"), y("y"), z("z");
    int pass = 0, fail = 0;

    for (int trial = 0; trial < N_TRIALS; trial++) {
        auto f1 = random_polynomial<ZZ>({x, y, z}, 2, 3, {-3, 3});
        auto f2 = random_polynomial<ZZ>({x, y, z}, 2, 3, {-3, 3});
        auto f3 = random_polynomial<ZZ>({x, y, z}, 2, 3, {-3, 3});
        if (f1.empty() || f2.empty() || f3.empty()) {
            if (verbose) std::cerr << "trial " << trial << " skip" << std::endl;
            continue;
        }
        polynomial_ZZ f = f1 * f2 * f3;
        if (f.empty()) {
            if (verbose) std::cerr << "trial " << trial << " skip" << std::endl;
            continue;
        }

        std::cerr << "trial " << trial << " ..." << std::flush;

        if (verbose) {
            std::cerr << std::endl;
            std::cerr << "  f1 = " << f1 << std::endl;
            std::cerr << "  f2 = " << f2 << std::endl;
            std::cerr << "  f3 = " << f3 << std::endl;
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
            std::cerr << " FAIL" << std::endl;
            std::cerr << "  f1 = " << f1 << std::endl;
            std::cerr << "  f2 = " << f2 << std::endl;
            std::cerr << "  f3 = " << f3 << std::endl;
            fail++;
            if (fail >= 3) break;
        } else {
            std::cerr << " ok" << std::endl;
            pass++;
        }
    }

    CLPOLY_TEST("mass_random_trivar3fac");
    std::cout << "Trials: " << (pass + fail) << "  Pass: " << pass << "  Fail: " << fail << std::endl;
    CLPOLY_ASSERT_EQ(fail, 0);

    return clpoly_test::test_summary();
}
