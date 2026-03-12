/**
 * @file bench_bivar_deg.cc
 * @brief 探测双变量因式分解在更高次数下的性能
 */
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <chrono>
#include <iostream>
#include <iomanip>

using namespace clpoly;

template<typename Fn>
double time_ms(Fn fn) {
    auto t0 = std::chrono::high_resolution_clock::now();
    fn();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

int main()
{
    variable x("x"), y("y"), z("z");
    std::cout << "==== Bivariate: degree exploration ====" << std::endl;

    // 2 factors, varying degree
    std::cout << "\n--- 2 factors, coeff[-5,5] ---" << std::endl;
    int degs[] = {3, 4, 5, 6, 7, 8, 10};
    for (int d : degs) {
        auto f1 = random_polynomial<ZZ>({x, y}, d, d+2, {-5, 5});
        auto f2 = random_polynomial<ZZ>({x, y}, d, d+2, {-5, 5});
        if (f1.empty() || f2.empty()) continue;
        auto f = f1 * f2;
        f.normalization();

        double ms = time_ms([&]{ volatile auto r = factorize(f); });
        std::cout << "  deg" << d << "*deg" << d
                  << "  total_terms=" << std::setw(4) << f.size()
                  << "  " << std::fixed << std::setprecision(1) << std::setw(9) << ms << "ms"
                  << std::endl;
        if (ms > 5000) break;
    }

    // 3 factors, varying degree
    std::cout << "\n--- 3 factors, coeff[-5,5] ---" << std::endl;
    int degs3[] = {2, 3, 4, 5, 6};
    for (int d : degs3) {
        auto f1 = random_polynomial<ZZ>({x, y}, d, d+2, {-5, 5});
        auto f2 = random_polynomial<ZZ>({x, y}, d, d+2, {-5, 5});
        auto f3 = random_polynomial<ZZ>({x, y}, d, d+2, {-5, 5});
        if (f1.empty() || f2.empty() || f3.empty()) continue;
        auto f = f1 * f2 * f3;
        f.normalization();

        double ms = time_ms([&]{ volatile auto r = factorize(f); });
        std::cout << "  deg" << d << "*3"
                  << "  total_terms=" << std::setw(4) << f.size()
                  << "  " << std::fixed << std::setprecision(1) << std::setw(9) << ms << "ms"
                  << std::endl;
        if (ms > 5000) break;
    }

    // 2 factors, large coeff, varying degree
    std::cout << "\n--- 2 factors, coeff[-500,500] ---" << std::endl;
    for (int d : degs) {
        auto f1 = random_polynomial<ZZ>({x, y}, d, d+2, {-500, 500});
        auto f2 = random_polynomial<ZZ>({x, y}, d, d+2, {-500, 500});
        if (f1.empty() || f2.empty()) continue;
        auto f = f1 * f2;
        f.normalization();

        double ms = time_ms([&]{ volatile auto r = factorize(f); });
        std::cout << "  deg" << d << "*deg" << d
                  << "  total_terms=" << std::setw(4) << f.size()
                  << "  " << std::fixed << std::setprecision(1) << std::setw(9) << ms << "ms"
                  << std::endl;
        if (ms > 5000) break;
    }

    // Trivariate: varying degree
    std::cout << "\n==== Trivariate: degree exploration ====" << std::endl;
    std::cout << "\n--- 2 factors, coeff[-5,5] ---" << std::endl;
    int tdegs[] = {2, 3, 4, 5};
    for (int d : tdegs) {
        auto f1 = random_polynomial<ZZ>({x, y, z}, d, d+2, {-5, 5});
        auto f2 = random_polynomial<ZZ>({x, y, z}, d, d+2, {-5, 5});
        if (f1.empty() || f2.empty()) continue;
        auto f = f1 * f2;
        f.normalization();

        double ms = time_ms([&]{ volatile auto r = factorize(f); });
        std::cout << "  deg" << d << "*deg" << d
                  << "  total_terms=" << std::setw(4) << f.size()
                  << "  " << std::fixed << std::setprecision(1) << std::setw(9) << ms << "ms"
                  << std::endl;
        if (ms > 5000) break;
    }

    return 0;
}
