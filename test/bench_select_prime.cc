/**
 * @file bench_select_prime.cc
 * @brief 实验：对比小素数 vs 大素数在 DDF+EDF 中的开销
 *
 * 对同一多项式，分别测量：
 *   A) 3 个小素数（从 p=2 顺序枚举）的 DDF+EDF 总时间
 *   B) 3 个随机 63-bit 素数的 DDF+EDF 总时间
 *   C) 单个大素数 p = 2^64 - 59 的 DDF+EDF 时间
 *
 * 同时记录每种方案选出的 r（因子数），以评估 r 最小化效果。
 *
 * Build: make test/bench_select_prime && _build/release/bin/bench_select_prime
 */
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <string>
#include <cassert>

using namespace clpoly;

// --- 复用 polynomial_factorize_zp.hh 中的内部函数 ---
// 这些函数在 clpoly namespace 中，但是 inline 在 .hh 中

struct ddf_edf_result {
    uint64_t prime;
    size_t r;           // 因子数
    double ms;          // 耗时
};

// 对单个素数 p 做一次完整的 mod + sqfree_check + DDF + EDF
// 返回 {p, r, ms}；若素数不合法返回 r=0
ddf_edf_result try_one_prime(const upolynomial_<ZZ>& f, uint64_t p)
{
    int64_t deg_f = get_deg(f);
    ZZ lc_f = f.front().second;
    std::mt19937 rng(42);

    auto t0 = std::chrono::high_resolution_clock::now();

    // lc(f) mod p == 0?
    ZZ lc_mod;
    ZZ::fdiv_r(lc_mod, lc_f, ZZ(p));
    if (!lc_mod) return {p, 0, 0.0};

    // f mod p
    auto fp = polynomial_mod(f, p);
    if (fp.empty() || get_deg(fp) != deg_f)
        return {p, 0, 0.0};

    // squarefree check: gcd(fp, fp') == 1
    auto fp_deriv = derivative(fp);
    if (fp_deriv.empty()) return {p, 0, 0.0};
    auto g = polynomial_GCD(fp, fp_deriv);
    if (get_deg(g) > 0) return {p, 0, 0.0};

    // DDF + EDF
    __upoly_make_monic(fp);
    auto ddf = __ddf_Zp(fp);

    std::vector<upolynomial_<Zp>> irr_factors;
    for (auto& [gk, dk] : ddf) {
        std::vector<upolynomial_<Zp>> edf_out;
        __edf_Zp(edf_out, gk, dk, rng);
        for (auto& hi : edf_out)
            irr_factors.push_back(std::move(hi));
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    return {p, irr_factors.size(), ms};
}

// 方案 A: 3 个小素数（顺序枚举）
struct scheme_result {
    double total_ms;
    size_t best_r;
    std::vector<ddf_edf_result> trials;
};

scheme_result scheme_small_primes(const upolynomial_<ZZ>& f, size_t max_tries = 3)
{
    scheme_result res{0.0, SIZE_MAX, {}};
    size_t tried = 0;
    for (uint64_t p = 2; tried < max_tries; p = next_prime_64(p)) {
        auto r = try_one_prime(f, p);
        if (r.r == 0) continue;  // 不合法，跳过
        ++tried;
        res.total_ms += r.ms;
        res.trials.push_back(r);
        if (r.r < res.best_r) res.best_r = r.r;
    }
    return res;
}

// 方案 B: 3 个随机 63-bit 素数
scheme_result scheme_random_large_primes(const upolynomial_<ZZ>& f, size_t max_tries = 3)
{
    scheme_result res{0.0, SIZE_MAX, {}};
    // 随机 63-bit 素数: 从 [2^62, 2^63] 范围随机取，再找下一个素数
    std::mt19937_64 rng(12345);
    std::uniform_int_distribution<uint64_t> dist(1ULL << 62, (1ULL << 63) - 1);

    size_t tried = 0;
    size_t attempts = 0;
    while (tried < max_tries && attempts < 100) {
        uint64_t candidate = dist(rng);
        uint64_t p = next_prime_64(candidate);
        ++attempts;
        auto r = try_one_prime(f, p);
        if (r.r == 0) continue;
        ++tried;
        res.total_ms += r.ms;
        res.trials.push_back(r);
        if (r.r < res.best_r) res.best_r = r.r;
    }
    return res;
}

// 方案 C: 单个大素数 2^64 - 59
scheme_result scheme_single_large(const upolynomial_<ZZ>& f)
{
    scheme_result res{0.0, SIZE_MAX, {}};
    uint64_t p = UINT64_MAX - 58;  // 2^64 - 59
    auto r = try_one_prime(f, p);
    if (r.r > 0) {
        res.total_ms = r.ms;
        res.best_r = r.r;
        res.trials.push_back(r);
    }
    return res;
}

void print_scheme(const std::string& name, const scheme_result& res)
{
    std::cout << "  " << std::left << std::setw(30) << name
              << std::right << std::setw(10) << std::fixed << std::setprecision(3) << res.total_ms << "ms"
              << "  r_best=" << res.best_r;
    std::cout << "  [";
    for (size_t i = 0; i < res.trials.size(); ++i) {
        if (i > 0) std::cout << ", ";
        std::cout << "p=";
        if (res.trials[i].prime > (1ULL << 32))
            std::cout << "2^" << (int)(log2((double)res.trials[i].prime) + 0.5);
        else
            std::cout << res.trials[i].prime;
        std::cout << ":r=" << res.trials[i].r
                  << ":" << std::fixed << std::setprecision(3) << res.trials[i].ms << "ms";
    }
    std::cout << "]" << std::endl;
}

void run_experiment(const std::string& label, const upolynomial_<ZZ>& f, int repeats = 5)
{
    std::cout << "\n--- " << label << " (deg=" << get_deg(f)
              << ", terms=" << f.size() << ") ---" << std::endl;

    // 多次运行取平均
    scheme_result avg_a{}, avg_b{}, avg_c{};
    scheme_result last_a, last_b, last_c;

    for (int rep = 0; rep < repeats; ++rep) {
        last_a = scheme_small_primes(f);
        last_b = scheme_random_large_primes(f);
        last_c = scheme_single_large(f);
        avg_a.total_ms += last_a.total_ms;
        avg_b.total_ms += last_b.total_ms;
        avg_c.total_ms += last_c.total_ms;
    }
    avg_a.total_ms /= repeats; avg_a.best_r = last_a.best_r; avg_a.trials = last_a.trials;
    avg_b.total_ms /= repeats; avg_b.best_r = last_b.best_r; avg_b.trials = last_b.trials;
    avg_c.total_ms /= repeats; avg_c.best_r = last_c.best_r; avg_c.trials = last_c.trials;

    print_scheme("A: 3 small primes (2,3,5..)", avg_a);
    print_scheme("B: 3 random 63-bit primes", avg_b);
    print_scheme("C: single p=2^64-59", avg_c);

    // 对比
    if (avg_a.total_ms > 0 && avg_c.total_ms > 0) {
        std::cout << "  >> C/A time ratio: " << std::fixed << std::setprecision(2)
                  << avg_c.total_ms / avg_a.total_ms << "x" << std::endl;
    }
    if (avg_a.total_ms > 0 && avg_b.total_ms > 0) {
        std::cout << "  >> B/A time ratio: " << std::fixed << std::setprecision(2)
                  << avg_b.total_ms / avg_a.total_ms << "x" << std::endl;
    }
}

int main()
{
    variable x("x");
    std::cout << "==== __select_prime: small vs large prime DDF+EDF experiment ====" << std::endl;

    // 用例 1: ~deg15 (3 factors)
    {
        upolynomial_ZZ uf1({{5, ZZ(1)}, {3, ZZ(-2)}, {1, ZZ(3)}, {0, ZZ(1)}});
        upolynomial_ZZ uf2({{4, ZZ(1)}, {2, ZZ(1)}, {0, ZZ(-1)}});
        upolynomial_ZZ uf3({{6, ZZ(2)}, {3, ZZ(-1)}, {1, ZZ(4)}, {0, ZZ(-3)}});
        auto f = uf1 * uf2 * uf3;
        run_experiment("~deg15 (3 factors)", f);
    }

    // 用例 2: ~deg29 (5 factors)
    {
        upolynomial_ZZ uf1({{5, ZZ(1)}, {3, ZZ(-2)}, {1, ZZ(3)}, {0, ZZ(1)}});
        upolynomial_ZZ uf2({{4, ZZ(1)}, {2, ZZ(1)}, {0, ZZ(-1)}});
        upolynomial_ZZ uf3({{6, ZZ(2)}, {3, ZZ(-1)}, {1, ZZ(4)}, {0, ZZ(-3)}});
        upolynomial_ZZ uf4({{8, ZZ(1)}, {4, ZZ(-3)}, {2, ZZ(2)}, {0, ZZ(5)}});
        upolynomial_ZZ uf5({{6, ZZ(3)}, {3, ZZ(2)}, {1, ZZ(-1)}, {0, ZZ(7)}});
        auto f = uf1 * uf2 * uf3 * uf4 * uf5;
        run_experiment("~deg29 (5 factors)", f);
    }

    // 用例 3: Wilkinson W(15)
    {
        upolynomial_ZZ wilk({{1, ZZ(1)}, {0, ZZ(-1)}});
        for (int k = 2; k <= 15; ++k) {
            upolynomial_ZZ lin({{1, ZZ(1)}, {0, ZZ(-k)}});
            wilk = wilk * lin;
        }
        run_experiment("Wilkinson W(15)", wilk);
    }

    // 用例 4: Wilkinson W(25)
    {
        upolynomial_ZZ wilk({{1, ZZ(1)}, {0, ZZ(-1)}});
        for (int k = 2; k <= 25; ++k) {
            upolynomial_ZZ lin({{1, ZZ(1)}, {0, ZZ(-k)}});
            wilk = wilk * lin;
        }
        run_experiment("Wilkinson W(25)", wilk);
    }

    // 用例 5: B2 大系数多项式
    {
        upolynomial_<ZZ> f;
        f.push_back({umonomial(4), ZZ(96000)});
        f.push_back({umonomial(3), ZZ("38464001")});
        f.push_back({umonomial(2), ZZ("-61401631999")});
        f.push_back({umonomial(1), ZZ("-24616960640000")});
        f.push_back({umonomial(0), ZZ("-24555520640000")});
        run_experiment("B2 large-coeff (deg4, 4 factors)", f, 20);
    }

    // 用例 6: x^24-1 (cyclotomic, many factors)
    {
        upolynomial_ZZ f({{24, ZZ(1)}, {0, ZZ(-1)}});
        run_experiment("x^24-1 (cyclotomic)", f);
    }

    // 用例 7: Wilkinson W(40)
    {
        upolynomial_ZZ f({{1, ZZ(1)}, {0, ZZ(-1)}});
        for (int k = 2; k <= 40; ++k) {
            upolynomial_ZZ lin({{1, ZZ(1)}, {0, ZZ(-k)}});
            f = f * lin;
        }
        run_experiment("Wilkinson W(40)", f, 1);
    }

    return 0;
}
