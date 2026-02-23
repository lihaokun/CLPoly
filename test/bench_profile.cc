// bench_profile.cc — 手动计时 profiling
// 编译：make _build/profile/bench_profile  (uses -O2 -pg -DCLPOLY_PROFILE)
// 运行：_build/profile/bench_profile

#define CLPOLY_PROFILE
#include "clpoly/clpoly.hh"
#include "clpoly/polynomial_factorize.hh"
#include <iostream>
#include <random>
#include <chrono>

using namespace clpoly;

static upolynomial_ZZ make_random_poly(std::mt19937& rng, int n_factors)
{
    std::uniform_int_distribution<int> coeff(-9, 9);
    upolynomial_ZZ f({{0, ZZ(1)}});
    for (int i = 0; i < n_factors; ++i) {
        int a = 0; while (a == 0) a = coeff(rng);
        int b = 0; while (b == 0) b = coeff(rng);
        f = f * upolynomial_ZZ({{1, ZZ(a)}, {0, ZZ(b)}});
    }
    return f;
}

static upolynomial_ZZ make_wilkinson(int n)
{
    upolynomial_ZZ f({{1, ZZ(1)}, {0, ZZ(-1)}});
    for (int k = 2; k <= n; ++k)
        f = f * upolynomial_ZZ({{1, ZZ(1)}, {0, ZZ(-k)}});
    return f;
}

static void run_case(const char* label, const upolynomial_ZZ& poly, int reps)
{
    reset_factorize_profile();
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < reps; ++i) {
        volatile auto r = factorize(poly); (void)r;
    }
    auto t1 = std::chrono::steady_clock::now();
    double wall_ms = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / 1000.0;
    std::printf("──────────────────────────────────────────────\n");
    std::printf("Case: %s  (x%d, wall=%.1fms, %.2fms/call)\n",
        label, reps, wall_ms, wall_ms / reps);
    print_factorize_profile();
}

int main()
{
    std::mt19937 rng(42);
    auto poly3 = make_random_poly(rng, 3);   // ~deg15 (3 fac)
    auto poly5 = make_random_poly(rng, 5);   // ~deg29 (5 fac)
    auto w20   = make_wilkinson(20);          // W(20)

    // uni-70: small-coeff 70-factor poly
    std::mt19937 rng2(123);
    std::uniform_int_distribution<int> sm(-3, 3);
    upolynomial_ZZ poly70({{0, ZZ(1)}});
    for (int i = 0; i < 70; ++i) {
        int a = 0; while (a == 0) a = sm(rng2);
        int b = 0; while (b == 0) b = sm(rng2);
        poly70 = poly70 * upolynomial_ZZ({{1, ZZ(a)}, {0, ZZ(b)}});
    }

    run_case("~deg15 (3 fac)", poly3, 2000);
    run_case("~deg29 (5 fac)", poly5,  500);
    run_case("Wilkinson W(20)", w20,   500);
    run_case("uni-70 factors",  poly70, 100);

    return 0;
}
