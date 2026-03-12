#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <iostream>
#include <iomanip>

using namespace clpoly;

int main() {
    // 用例 3: large-coeff 6 factors (deg10)
    upolynomial_ZZ f1({{2, ZZ(1)}, {1, ZZ(500)}, {0, ZZ(1)}});
    upolynomial_ZZ f2({{2, ZZ(1)}, {1, ZZ(-300)}, {0, ZZ(7)}});
    upolynomial_ZZ f3({{1, ZZ(1)}, {0, ZZ(200)}});
    upolynomial_ZZ f4({{1, ZZ(1)}, {0, ZZ(-400)}});
    upolynomial_ZZ f5({{1, ZZ(2)}, {0, ZZ(999)}});
    upolynomial_ZZ f6({{2, ZZ(1)}, {1, ZZ(123)}, {0, ZZ(-456)}});
    auto f = f1 * f2 * f3 * f4 * f5 * f6;

    std::cout << "f: deg=" << get_deg(f) << " terms=" << f.size()
              << " lc=" << f.front().second << std::endl;
    std::cout << "true factors = 6" << std::endl;

    ZZ lc_f = f.front().second;
    int64_t deg_f = get_deg(f);
    std::mt19937 rng(42);

    // 扫描小素数
    std::cout << "\n--- Small primes (p=2..500) ---" << std::endl;
    std::cout << std::setw(8) << "p" << std::setw(8) << "r" << std::setw(10) << "status" << std::endl;
    int best_r = 999; uint64_t best_p = 0;
    for (uint64_t p = 2; p < 500; p = next_prime_64(p)) {
        ZZ lc_mod;
        ZZ::fdiv_r(lc_mod, lc_f, ZZ(p));
        if (!lc_mod) { continue; }

        auto fp = polynomial_mod(f, p);
        if (fp.empty() || get_deg(fp) != deg_f) continue;

        auto fp_deriv = derivative(fp);
        if (fp_deriv.empty()) continue;
        auto g = polynomial_GCD(fp, fp_deriv);
        if (get_deg(g) > 0) { std::cout << std::setw(8) << p << "  not squarefree" << std::endl; continue; }

        __upoly_make_monic(fp);
        auto ddf = __ddf_Zp(fp);
        int r = 0;
        for (auto& [gk, dk] : ddf) {
            std::vector<upolynomial_<Zp>> edf_out;
            __edf_Zp(edf_out, gk, dk, rng);
            r += edf_out.size();
        }
        std::cout << std::setw(8) << p << std::setw(8) << r;
        if (r == 6) std::cout << "  *** optimal";
        std::cout << std::endl;
        if (r < best_r) { best_r = r; best_p = p; }
    }
    std::cout << "Best small: p=" << best_p << " r=" << best_r << std::endl;

    // 扫描大素数
    std::cout << "\n--- Large primes (top 20 from 2^64-59) ---" << std::endl;
    uint64_t lp = UINT64_MAX - 58;
    for (int i = 0; i < 20; ++i, lp = prev_prime_64(lp)) {
        ZZ lc_mod;
        ZZ::fdiv_r(lc_mod, lc_f, ZZ(lp));
        if (!lc_mod) { std::cout << std::setw(22) << lp << "  lc=0 skip" << std::endl; continue; }

        auto fp = polynomial_mod(f, lp);
        if (fp.empty() || get_deg(fp) != deg_f) continue;

        auto fp_deriv = derivative(fp);
        if (fp_deriv.empty()) continue;
        auto g = polynomial_GCD(fp, fp_deriv);
        if (get_deg(g) > 0) { std::cout << std::setw(22) << lp << "  not squarefree" << std::endl; continue; }

        __upoly_make_monic(fp);
        auto ddf = __ddf_Zp(fp);
        int r = 0;
        for (auto& [gk, dk] : ddf) {
            std::vector<upolynomial_<Zp>> edf_out;
            __edf_Zp(edf_out, gk, dk, rng);
            r += edf_out.size();
        }
        std::cout << std::setw(22) << lp << std::setw(8) << r;
        if (r == 6) std::cout << "  *** optimal";
        std::cout << std::endl;
    }

    return 0;
}
