/**
 * @file bench_prime_mode.cc
 * @brief 端到端 benchmark: SMALL vs LARGE 素数模式切换边界实验
 *
 * 通过 __g_use_large_prime 全局开关控制素数选取模式，
 * 直接调用完整 factorize() 测量端到端性能。
 *
 * Build (release):
 *   g++ -O3 -DNDEBUG -I./ test/bench_prime_mode.cc -o _build/release/bin/bench_prime_mode lib/libclpoly.a -lgmpxx -lgmp
 */
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <random>

using namespace clpoly;

namespace {

struct Timer {
    std::chrono::high_resolution_clock::time_point t0;
    Timer() : t0(std::chrono::high_resolution_clock::now()) {}
    double ms() const {
        return std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - t0).count();
    }
};

struct Result {
    double ms;
    size_t nfactors;
    bool correct;
};

bool verify(const upolynomial_<ZZ>& f, const factorization<upolynomial_<ZZ>>& fac)
{
    upolynomial_<ZZ> prod;
    prod.push_back(std::make_pair(umonomial(0), fac.content));
    for (auto& [fi, ei] : fac.factors) {
        auto fp = pow(fi, (int)ei);
        pair_vec_multiplies(prod.data(), prod.data(), fp.data(), prod.comp());
        prod.normalization();
    }
    return prod == f;
}

Result bench(const upolynomial_<ZZ>& f, bool large, int reps)
{
    __g_use_large_prime = large;
    Result r;
    double total = 0;
    for (int i = 0; i < reps; ++i) {
        Timer t;
        auto fac = factorize(f);
        total += t.ms();
        if (i == 0) {
            r.nfactors = fac.factors.size();
            r.correct = verify(f, fac);
        }
    }
    r.ms = total / reps;
    __g_use_large_prime = false;
    return r;
}

// 构造辅助
upolynomial_<ZZ> linear(int64_t a) {
    upolynomial_<ZZ> f;
    f.push_back(std::make_pair(umonomial(1), ZZ(1)));
    f.push_back(std::make_pair(umonomial(0), ZZ(a)));
    return f;
}

upolynomial_<ZZ> affine(int64_t lc, int64_t c0) {
    upolynomial_<ZZ> f;
    f.push_back(std::make_pair(umonomial(1), ZZ(lc)));
    f.push_back(std::make_pair(umonomial(0), ZZ(c0)));
    return f;
}

upolynomial_<ZZ> mul_all(std::initializer_list<upolynomial_<ZZ>> fs)
{
    upolynomial_<ZZ> prod;
    prod.push_back(std::make_pair(umonomial(0), ZZ(1)));
    for (auto& fi : fs) {
        pair_vec_multiplies(prod.data(), prod.data(), fi.data(), prod.comp());
        prod.normalization();
    }
    return prod;
}

void header()
{
    std::cout << std::left << std::setw(48) << "Case"
              << std::right
              << std::setw(10) << "SMALL"
              << std::setw(7) << "nfac"
              << std::setw(10) << "LARGE"
              << std::setw(7) << "nfac"
              << std::setw(8) << "S/L"
              << std::setw(8) << "best"
              << std::endl;
    std::cout << std::string(98, '-') << std::endl;
}

void row(const std::string& name, Result s, Result l)
{
    double ratio = (l.ms > 0.001) ? s.ms / l.ms : 999.0;
    const char* w = (ratio < 0.9) ? "SMALL" : (ratio > 1.1) ? "LARGE" : "tie";
    std::cout << std::left << std::setw(48) << name
              << std::right << std::fixed << std::setprecision(3)
              << std::setw(10) << s.ms
              << std::setw(7) << s.nfactors
              << std::setw(10) << l.ms
              << std::setw(7) << l.nfactors
              << std::setw(8) << std::setprecision(2) << ratio
              << std::setw(8) << w;
    if (!s.correct) std::cout << " [S-FAIL]";
    if (!l.correct) std::cout << " [L-FAIL]";
    std::cout << std::endl;
}

} // anon

int main()
{
    std::mt19937 rng(12345);
    std::cout << "==== Prime Mode Benchmark: SMALL vs LARGE ====\n";

    // §1. 系数规模 (deg5, 4 distinct factors)
    // f = (x+α)(x-α)(x+2)(3x+7)  [all distinct roots when α≥3]
    std::cout << "\n-- §1: coeff magnitude (deg5=4fac, (x+a)(x-a)(x+2)(3x+7)) --\n";
    header();
    for (int64_t a : {3, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 1000000}) {
        auto f = mul_all({linear(a), linear(-a), linear(2), affine(3, 7)});
        int reps = (a <= 100) ? 30 : 10;
        row("alpha=" + std::to_string(a), bench(f, false, reps), bench(f, true, reps));
    }

    // §2. 系数规模 (deg~10, 5 random quadratic factors)
    std::cout << "\n-- §2: coeff magnitude (deg~10, 5 random quadratic factors) --\n";
    header();
    for (int C : {5, 10, 50, 100, 500, 1000, 5000, 10000}) {
        double sm = 0, lg = 0; size_t sn = 0, ln = 0; bool sok = true, lok = true;
        int runs = 5;
        for (int t = 0; t < runs; ++t) {
            upolynomial_<ZZ> prod;
            prod.push_back(std::make_pair(umonomial(0), ZZ(1)));
            for (int i = 0; i < 5; ++i) {
                auto fi = random_upolynomial<ZZ>(2, 2, {-C, C});
                if (fi.empty() || get_deg(fi) < 1) { fi = linear(rng() % 100 + 1); }
                pair_vec_multiplies(prod.data(), prod.data(), fi.data(), prod.comp());
                prod.normalization();
            }
            auto rs = bench(prod, false, 1);
            auto rl = bench(prod, true, 1);
            sm += rs.ms; lg += rl.ms;
            sn = rs.nfactors; ln = rl.nfactors;
            if (!rs.correct) sok = false;
            if (!rl.correct) lok = false;
        }
        row("deg~10 5fac C=" + std::to_string(C),
            Result{sm/runs, sn, sok}, Result{lg/runs, ln, lok});
    }

    // §3. 次数 (大系数线性因子)
    std::cout << "\n-- §3: degree (linear factors, roots 100,150,200,...) --\n";
    header();
    for (int n : {4, 6, 8, 10, 15, 20, 25, 30}) {
        upolynomial_<ZZ> f;
        f.push_back(std::make_pair(umonomial(0), ZZ(1)));
        for (int i = 0; i < n; ++i) {
            auto fi = linear(100 + i * 50);
            pair_vec_multiplies(f.data(), f.data(), fi.data(), f.comp());
            f.normalization();
        }
        int reps = (n <= 10) ? 5 : 1;
        row("deg" + std::to_string(n) + " " + std::to_string(n) + "fac large",
            bench(f, false, reps), bench(f, true, reps));
    }

    // §4. 次数 (小系数二次因子)
    std::cout << "\n-- §4: degree (quadratic factors, C=5) --\n";
    header();
    for (int n : {2, 3, 4, 5, 6, 8, 10}) {
        double sm = 0, lg = 0; size_t sn = 0, ln = 0; bool sok = true, lok = true;
        int runs = 5;
        for (int t = 0; t < runs; ++t) {
            upolynomial_<ZZ> prod;
            prod.push_back(std::make_pair(umonomial(0), ZZ(1)));
            for (int i = 0; i < n; ++i) {
                auto fi = random_upolynomial<ZZ>(2, 2, {-5, 5});
                if (fi.empty() || get_deg(fi) < 1) { fi = linear(rng() % 10 + 1); }
                pair_vec_multiplies(prod.data(), prod.data(), fi.data(), prod.comp());
                prod.normalization();
            }
            auto rs = bench(prod, false, 1);
            auto rl = bench(prod, true, 1);
            sm += rs.ms; lg += rl.ms;
            sn = rs.nfactors; ln = rl.nfactors;
            if (!rs.correct) sok = false;
            if (!rl.correct) lok = false;
        }
        row("deg" + std::to_string(n*2) + " " + std::to_string(n) + "fac C=5",
            Result{sm/runs, sn, sok}, Result{lg/runs, ln, lok});
    }

    // §5. 经典用例
    std::cout << "\n-- §5: classic --\n";
    header();
    for (int n : {10, 15, 20, 25}) {
        upolynomial_<ZZ> f;
        f.push_back(std::make_pair(umonomial(0), ZZ(1)));
        for (int i = 1; i <= n; ++i) {
            auto fi = linear(-i);
            pair_vec_multiplies(f.data(), f.data(), fi.data(), f.comp());
            f.normalization();
        }
        row("Wilkinson W(" + std::to_string(n) + ")",
            bench(f, false, 3), bench(f, true, 3));
    }
    // Swinnerton-Dyer S3
    {
        upolynomial_<ZZ> f;
        f.push_back(std::make_pair(umonomial(8), ZZ(1)));
        f.push_back(std::make_pair(umonomial(6), ZZ(-40)));
        f.push_back(std::make_pair(umonomial(4), ZZ(352)));
        f.push_back(std::make_pair(umonomial(2), ZZ(-960)));
        f.push_back(std::make_pair(umonomial(0), ZZ(576)));
        row("Swinnerton-Dyer S3", bench(f, false, 5), bench(f, true, 5));
    }

    // §6. B2 风格大系数用例
    std::cout << "\n-- §6: B2-style large coeff --\n";
    header();
    {
        // B2 原始
        upolynomial_<ZZ> f;
        f.push_back(std::make_pair(umonomial(4), ZZ(96000)));
        f.push_back(std::make_pair(umonomial(3), ZZ("38464001")));
        f.push_back(std::make_pair(umonomial(2), ZZ("-61401631999")));
        f.push_back(std::make_pair(umonomial(1), ZZ("-24616960640000")));
        f.push_back(std::make_pair(umonomial(0), ZZ("-24555520640000")));
        row("B2 original (deg4)", bench(f, false, 20), bench(f, true, 20));
    }
    {
        // 6 因子 deg8
        auto f = mul_all({
            linear(800), linear(-800), linear(1200), linear(-1200),
            affine(5, 7), affine(11, 17)
        });
        row("6fac deg6 alpha~1200", bench(f, false, 10), bench(f, true, 10));
    }
    {
        // 8 因子 deg8
        auto f = mul_all({
            linear(500), linear(-500), linear(1000), linear(-1000),
            linear(2000), linear(-2000), affine(7, 11), affine(13, 19)
        });
        row("8fac deg8 alpha~2000", bench(f, false, 5), bench(f, true, 5));
    }
    {
        // 大 lc 4因子
        auto f = mul_all({
            affine(96000, 38368001), linear(1), linear(-800), linear(800)
        });
        row("4fac deg4 lc=96000", bench(f, false, 20), bench(f, true, 20));
    }

    // §7. 极端大系数 — 推高系数膨胀以寻找大素数优势边界
    std::cout << "\n-- §7: extreme large coeff (pushing coefficient inflation) --\n";
    header();
    // 7a: deg10 = 4 quadratic + 2 linear, 大根
    for (int64_t a : {100, 1000, 10000, 100000, 1000000}) {
        upolynomial_<ZZ> f1, f2, f3, f4, f5, f6;
        // x^2 + a*x + 1
        f1.push_back(std::make_pair(umonomial(2), ZZ(1)));
        f1.push_back(std::make_pair(umonomial(1), ZZ(a)));
        f1.push_back(std::make_pair(umonomial(0), ZZ(1)));
        // x^2 - a*x + 3
        f2.push_back(std::make_pair(umonomial(2), ZZ(1)));
        f2.push_back(std::make_pair(umonomial(1), ZZ(-a)));
        f2.push_back(std::make_pair(umonomial(0), ZZ(3)));
        // x^2 + (a/2)*x + 7
        f3.push_back(std::make_pair(umonomial(2), ZZ(1)));
        f3.push_back(std::make_pair(umonomial(1), ZZ(a/2)));
        f3.push_back(std::make_pair(umonomial(0), ZZ(7)));
        // x^2 - 1
        f4.push_back(std::make_pair(umonomial(2), ZZ(1)));
        f4.push_back(std::make_pair(umonomial(0), ZZ(-1)));
        // x + a
        f5.push_back(std::make_pair(umonomial(1), ZZ(1)));
        f5.push_back(std::make_pair(umonomial(0), ZZ(a)));
        // x - a
        f6.push_back(std::make_pair(umonomial(1), ZZ(1)));
        f6.push_back(std::make_pair(umonomial(0), ZZ(-a)));
        auto f = mul_all({f1, f2, f3, f4, f5, f6});
        int reps = (a <= 1000) ? 10 : 3;
        row("deg10 6fac a=" + std::to_string(a), bench(f, false, reps), bench(f, true, reps));
    }
    // 7b: deg15 = 5 quadratic + 5 linear, 大根
    for (int64_t a : {1000, 10000, 100000, 1000000}) {
        auto f = mul_all({
            linear(a), linear(-a), linear(a/2), linear(-a/3), linear(a/5 + 1),
            affine(1, a*a), affine(1, -a*a + 3), affine(1, a*a/2 + 7),
            affine(3, a + 1), affine(7, a - 1)
        });
        int reps = (a <= 1000) ? 5 : 2;
        row("deg10 10fac a=" + std::to_string(a), bench(f, false, reps), bench(f, true, reps));
    }
    // 7c: 非首一大 lc, deg8
    for (int64_t lc : {(int64_t)100, (int64_t)10000, (int64_t)1000000, (int64_t)100000000}) {
        auto f = mul_all({
            affine(lc, 1), affine(lc + 1, -1),
            linear(1000), linear(-1000),
            affine(3, 7), affine(11, 17),
            linear(500), linear(-500)
        });
        row("deg8 8fac lc=" + std::to_string(lc), bench(f, false, 5), bench(f, true, 5));
    }
    // 7d: 系数 10^18 级别 (接近 int64 极限)
    {
        int64_t big = 1000000000LL; // 10^9
        auto f = mul_all({
            linear(big), linear(-big), linear(big/2), linear(-big/3),
            affine(3, big + 1), affine(7, big - 1)
        });
        row("deg6 6fac a=1e9", bench(f, false, 3), bench(f, true, 3));
    }
    {
        int64_t big = 1000000000LL;
        auto f = mul_all({
            linear(big), linear(-big), linear(big/2), linear(-big/3),
            affine(3, big + 1), affine(7, big - 1),
            linear(big/5 + 1), linear(-big/7)
        });
        row("deg8 8fac a=1e9", bench(f, false, 3), bench(f, true, 3));
    }
    {
        int64_t big = 1000000000LL;
        auto f = mul_all({
            linear(big), linear(-big), linear(big/2), linear(-big/3),
            affine(3, big + 1), affine(7, big - 1),
            linear(big/5 + 1), linear(-big/7),
            affine(11, big/11), affine(13, -big/13)
        });
        row("deg10 10fac a=1e9", bench(f, false, 2), bench(f, true, 2));
    }

    std::cout << "\nDone." << std::endl;
    return 0;
}
