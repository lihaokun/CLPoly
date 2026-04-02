/**
 * @file bench_factorize_limits.cc
 * @brief 探测 factorize() 在不同维度上的性能边界
 *
 * 维度：系数范围 × 变量数 × 因子数 × 次数
 * 目标：找到 crosscheck 可承受的上限（单次 <1s）
 */
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>

using namespace clpoly;

template<typename Fn>
double time_ms(Fn fn) {
    auto t0 = std::chrono::high_resolution_clock::now();
    fn();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

void header(const std::string& title) {
    std::cout << "\n==== " << title << " ====" << std::endl;
}

int main()
{
    variable x("x"), y("y"), z("z"), w("w");
    std::cout << "==== factorize() performance boundary exploration ====" << std::endl;

    // ---- 1. 单变量：系数范围 ----
    header("Univariate: coefficient range (deg~12, 4 factors)");
    {
        int coeff_ranges[] = {5, 50, 500, 5000, 50000, 500000};
        for (int C : coeff_ranges) {
            // 4 个 deg3 因子
            auto f1 = random_upolynomial<ZZ>(3, 3, {-C, C});
            auto f2 = random_upolynomial<ZZ>(3, 3, {-C, C});
            auto f3 = random_upolynomial<ZZ>(3, 3, {-C, C});
            auto f4 = random_upolynomial<ZZ>(3, 3, {-C, C});
            if (f1.empty() || f2.empty() || f3.empty() || f4.empty()) continue;
            auto f = f1 * f2 * f3 * f4;

            double ms = time_ms([&]{ volatile auto r = factorize(f); });
            std::cout << "  coeff[-" << std::setw(6) << C << "," << std::setw(6) << C << "]"
                      << "  deg=" << std::setw(3) << get_deg(f)
                      << "  terms=" << std::setw(3) << f.size()
                      << "  " << std::fixed << std::setprecision(1) << std::setw(8) << ms << "ms"
                      << std::endl;
        }
    }

    // ---- 2. 单变量：因子数 ----
    header("Univariate: factor count (linear factors, coeff[-10,10])");
    {
        int nfacs[] = {5, 10, 15, 20, 30, 40, 50};
        for (int n : nfacs) {
            upolynomial_ZZ f({{1, ZZ(1)}, {0, ZZ(-1)}});
            for (int k = 2; k <= n; ++k) {
                upolynomial_ZZ lin({{1, ZZ(1)}, {0, ZZ(-k)}});
                f = f * lin;
            }
            double ms = time_ms([&]{ volatile auto r = factorize(f); });
            std::cout << "  " << std::setw(2) << n << " linear factors"
                      << "  deg=" << std::setw(3) << get_deg(f)
                      << "  " << std::fixed << std::setprecision(1) << std::setw(8) << ms << "ms"
                      << std::endl;
            if (ms > 5000) break;
        }
    }

    // ---- 3. 单变量：次数（少因子）----
    header("Univariate: degree (2 factors, coeff[-10,10])");
    {
        int degs[] = {5, 10, 15, 20, 30, 40};
        for (int d : degs) {
            auto f1 = random_upolynomial<ZZ>(d, d, {-10, 10});
            auto f2 = random_upolynomial<ZZ>(d, d, {-10, 10});
            if (f1.empty() || f2.empty()) continue;
            auto f = f1 * f2;

            double ms = time_ms([&]{ volatile auto r = factorize(f); });
            std::cout << "  2 x deg" << std::setw(2) << d
                      << "  total_deg=" << std::setw(3) << get_deg(f)
                      << "  terms=" << std::setw(3) << f.size()
                      << "  " << std::fixed << std::setprecision(1) << std::setw(8) << ms << "ms"
                      << std::endl;
            if (ms > 5000) break;
        }
    }

    // ---- 4. 双变量：系数范围 ----
    header("Bivariate: coefficient range (2 factors, deg~3)");
    {
        int coeff_ranges[] = {5, 50, 500, 5000, 50000};
        for (int C : coeff_ranges) {
            auto f1 = random_polynomial<ZZ>({x, y}, 3, 4, {-C, C});
            auto f2 = random_polynomial<ZZ>({x, y}, 3, 4, {-C, C});
            if (f1.empty() || f2.empty()) continue;
            auto f = f1 * f2;

            double ms = time_ms([&]{ volatile auto r = factorize(f); });
            std::cout << "  coeff[-" << std::setw(5) << C << "," << std::setw(5) << C << "]"
                      << "  terms=" << std::setw(3) << f.size()
                      << "  " << std::fixed << std::setprecision(1) << std::setw(8) << ms << "ms"
                      << std::endl;
            if (ms > 5000) break;
        }
    }

    // ---- 5. 双变量：因子数 ----
    header("Bivariate: factor count (deg~2 each, coeff[-3,3])");
    {
        int nfacs[] = {2, 3, 4, 5, 6, 8, 10};
        for (int n : nfacs) {
            polynomial_ZZ f = polynomial_ZZ(1);
            bool ok = true;
            for (int i = 0; i < n; ++i) {
                auto fi = random_polynomial<ZZ>({x, y}, 2, 3, {-3, 3});
                if (fi.empty()) { ok = false; break; }
                f = f * fi;
                f.normalization();
            }
            if (!ok) continue;

            double ms = time_ms([&]{ volatile auto r = factorize(f); });
            std::cout << "  " << std::setw(2) << n << " factors"
                      << "  terms=" << std::setw(4) << f.size()
                      << "  " << std::fixed << std::setprecision(1) << std::setw(8) << ms << "ms"
                      << std::endl;
            if (ms > 5000) break;
        }
    }

    // ---- 6. 三变量：因子数 ----
    header("Trivariate: factor count (deg~2 each, coeff[-3,3])");
    {
        int nfacs[] = {2, 3, 4, 5, 6};
        for (int n : nfacs) {
            polynomial_ZZ f = polynomial_ZZ(1);
            bool ok = true;
            for (int i = 0; i < n; ++i) {
                auto fi = random_polynomial<ZZ>({x, y, z}, 2, 3, {-3, 3});
                if (fi.empty()) { ok = false; break; }
                f = f * fi;
                f.normalization();
            }
            if (!ok) continue;

            double ms = time_ms([&]{ volatile auto r = factorize(f); });
            std::cout << "  " << std::setw(2) << n << " factors"
                      << "  terms=" << std::setw(4) << f.size()
                      << "  " << std::fixed << std::setprecision(1) << std::setw(8) << ms << "ms"
                      << std::endl;
            if (ms > 5000) break;
        }
    }

    // ---- 7. 四变量 ----
    header("4-variable: factor count (deg~2 each, coeff[-3,3])");
    {
        int nfacs[] = {2, 3, 4};
        for (int n : nfacs) {
            polynomial_ZZ f = polynomial_ZZ(1);
            bool ok = true;
            for (int i = 0; i < n; ++i) {
                auto fi = random_polynomial<ZZ>({x, y, z, w}, 2, 3, {-3, 3});
                if (fi.empty()) { ok = false; break; }
                f = f * fi;
                f.normalization();
            }
            if (!ok) continue;

            double ms = time_ms([&]{ volatile auto r = factorize(f); });
            std::cout << "  " << std::setw(2) << n << " factors"
                      << "  terms=" << std::setw(4) << f.size()
                      << "  " << std::fixed << std::setprecision(1) << std::setw(8) << ms << "ms"
                      << std::endl;
            if (ms > 5000) break;
        }
    }

    // ---- 8. 双变量大系数 + 多因子 ----
    header("Bivariate: large coeff + multiple factors");
    {
        struct Config { int nfac; int C; int deg; };
        Config cfgs[] = {
            {2, 100, 3}, {3, 100, 2}, {4, 100, 2},
            {2, 1000, 3}, {3, 1000, 2},
            {2, 10000, 3}, {3, 10000, 2},
        };
        for (auto& cfg : cfgs) {
            polynomial_ZZ f = polynomial_ZZ(1);
            bool ok = true;
            for (int i = 0; i < cfg.nfac; ++i) {
                auto fi = random_polynomial<ZZ>({x, y}, cfg.deg, cfg.deg+1, {-cfg.C, cfg.C});
                if (fi.empty()) { ok = false; break; }
                f = f * fi;
                f.normalization();
            }
            if (!ok) continue;

            double ms = time_ms([&]{ volatile auto r = factorize(f); });
            std::cout << "  " << cfg.nfac << "fac deg" << cfg.deg
                      << " C=" << std::setw(5) << cfg.C
                      << "  terms=" << std::setw(4) << f.size()
                      << "  " << std::fixed << std::setprecision(1) << std::setw(8) << ms << "ms"
                      << std::endl;
        }
    }

    return 0;
}
