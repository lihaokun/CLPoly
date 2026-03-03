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

    // Deterministic reproducer: Wang LCC non-divisor hang (trial 387)
    // L = lc(f,x) = -6y²(3y+1), γ=-6. CLPoly 旧 coprime 检查
    // gcd(Eⱼ, cs·γ)==1 对此多项式数学上不可满足，导致无限循环。
    // 修复：non-divisor 累积素因子剥离（GCL §8.7 condition 3）
    CLPOLY_TEST("repro_trivar3fac_nondivisor_hang");
    {
        polynomial_ZZ one({{{},(ZZ)1}});
        polynomial_ZZ two({{{},(ZZ)2}});
        polynomial_ZZ three({{{},(ZZ)3}});
        polynomial_ZZ f1 = -three*x*y - x + two;
        polynomial_ZZ f2 = three*x*y + three*y*z - y;
        polynomial_ZZ f3 = two*x*y - three*y*z - three;
        polynomial_ZZ f = f1 * f2 * f3;
        auto fac = factorize(f);
        polynomial_ZZ product({{{},(ZZ)fac.content}});
        for (auto& [fi, ei] : fac.factors)
            for (uint64_t e = 0; e < ei; e++) {
                product = product * fi;
                product.normalization();
            }
        CLPOLY_ASSERT(product == f || product == -f);
    }

    return clpoly_test::test_summary();
}
