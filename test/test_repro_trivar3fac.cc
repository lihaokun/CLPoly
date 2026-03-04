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
#include <cstdlib>

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

    // M2: 纯变量 LC 因子 + 小素数 γ（valuation 提取测试）
    // lc(f,x) 含 y^k 因子，需要拆分到多个提升因子
    CLPOLY_TEST("repro_trivar3fac_pure_var_lc");
    {
        polynomial_ZZ one({{{},(ZZ)1}});
        polynomial_ZZ two({{{},(ZZ)2}});
        // f1 = y*x + z, f2 = y*x - 1, f3 = 2*x + y + z
        // lc(f,x) = 2*y^2，γ=2，LC 因子含 y^2 需拆分为 y*y
        polynomial_ZZ f1 = y*x + z;
        polynomial_ZZ f2 = y*x - one;
        polynomial_ZZ f3 = two*x + y + z;
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

    // M2: LC 含多个不可约因子需分配到不同提升因子
    CLPOLY_TEST("repro_trivar3fac_multi_lc_factors");
    {
        polynomial_ZZ one({{{},(ZZ)1}});
        polynomial_ZZ two({{{},(ZZ)2}});
        polynomial_ZZ three({{{},(ZZ)3}});
        // f1 = y*x + 1, f2 = (y+z)*x - 2, f3 = 3*x + z - 1
        // lc(f,x) = 3*y*(y+z)，不可约因子 {y, y+z} 需分配到不同因子
        polynomial_ZZ f1 = y*x + one;
        polynomial_ZZ f2 = (y+z)*x - two;
        polynomial_ZZ f3 = three*x + z - one;
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

    // M2: 常数 LC（q==1 提前退出路径）
    // 所有因子的 lc(fi, x) 均为常数 → LC 分配平凡
    CLPOLY_TEST("repro_trivar3fac_const_lc");
    {
        polynomial_ZZ one({{{},(ZZ)1}});
        polynomial_ZZ two({{{},(ZZ)2}});
        // f1 = x + y + z, f2 = x - y + 1, f3 = x + 2z - 1
        // lc(f,x) = 1，LC 分配平凡
        polynomial_ZZ f1 = x + y + z;
        polynomial_ZZ f2 = x - y + one;
        polynomial_ZZ f3 = x + two*z - one;
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
