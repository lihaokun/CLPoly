/**
 * @file bench_phase2_cost.cc
 * @brief 实验：Phase 2 Hensel 提升开销 vs 大素数 DDF 开销
 *
 * 对 B2 大系数多项式和更大的合成用例，测量：
 *   1. 小素数下 Phase 1 + Phase 2 的 Hensel+重组 总时间
 *   2. 大素数下 DDF+EDF + Phase 1 Hensel+重组 总时间（无 Phase 2）
 *   3. 完整 factorize() 端到端时间
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

struct timing_result {
    double select_ms;     // DDF+EDF
    double hensel_p1_ms;  // Phase 1 Hensel lift
    double recom_p1_ms;   // Phase 1 recombine
    double hensel_p2_ms;  // Phase 2 Hensel lift (0 if no Phase 2)
    double recom_p2_ms;   // Phase 2 recombine (0 if no Phase 2)
    double total_ms;
    uint64_t prime;
    int r;
    int a_h;
    int a_mig;
    bool phase2_triggered;
    size_t result_factors;
};

// 对指定素数做完整 Hensel+重组流程，测量各阶段时间
timing_result run_with_prime(const upolynomial_<ZZ>& f, uint64_t p)
{
    timing_result res{};
    res.prime = p;
    int64_t deg_f = get_deg(f);
    ZZ lc_f = f.front().second;
    std::mt19937 rng(42);

    // --- select_prime 阶段：mod + sqfree + DDF + EDF ---
    auto t0 = std::chrono::high_resolution_clock::now();

    ZZ lc_mod;
    ZZ::fdiv_r(lc_mod, lc_f, ZZ(p));
    if (!lc_mod) { res.select_ms = -1; return res; }

    auto fp = polynomial_mod(f, p);
    if (fp.empty() || get_deg(fp) != deg_f) { res.select_ms = -1; return res; }

    auto fp_deriv = derivative(fp);
    if (fp_deriv.empty()) { res.select_ms = -1; return res; }
    auto g = polynomial_GCD(fp, fp_deriv);
    if (get_deg(g) > 0) { res.select_ms = -1; return res; }

    __upoly_make_monic(fp);
    auto ddf = __ddf_Zp(fp);

    std::vector<upolynomial_<Zp>> factors;
    for (auto& [gk, dk] : ddf) {
        std::vector<upolynomial_<Zp>> edf_out;
        __edf_Zp(edf_out, gk, dk, rng);
        for (auto& hi : edf_out)
            factors.push_back(std::move(hi));
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    res.select_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    res.r = (int)factors.size();

    if (res.r <= 1) {
        res.result_factors = 1;
        res.total_ms = res.select_ms;
        return res;
    }

    // --- 精度计算 ---
    auto [a_h, a_mig] = __heuristic_starting_precision(f, res.r, p);
    res.a_h = a_h;
    res.a_mig = a_mig;

    // --- Phase 1: Hensel lift to a_h ---
    auto t2 = std::chrono::high_resolution_clock::now();
    auto [lifted_h, m_h] = __hensel_lift(f, factors, p, a_h);
    auto t3 = std::chrono::high_resolution_clock::now();
    res.hensel_p1_ms = std::chrono::duration<double, std::milli>(t3 - t2).count();

    // --- Phase 1: recombine ---
    auto result = __vanhoeij_recombine(f, lifted_h, m_h);
    auto t4 = std::chrono::high_resolution_clock::now();
    res.recom_p1_ms = std::chrono::duration<double, std::milli>(t4 - t3).count();

    // --- Phase 2? ---
    if ((int)result.size() < res.r && a_h < a_mig) {
        res.phase2_triggered = true;

        auto t5 = std::chrono::high_resolution_clock::now();
        auto [lifted_mig, m_mig] = __hensel_lift(f, factors, p);
        auto t6 = std::chrono::high_resolution_clock::now();
        res.hensel_p2_ms = std::chrono::duration<double, std::milli>(t6 - t5).count();

        auto result2 = __vanhoeij_recombine(f, lifted_mig, m_mig);
        auto t7 = std::chrono::high_resolution_clock::now();
        res.recom_p2_ms = std::chrono::duration<double, std::milli>(t7 - t6).count();

        res.result_factors = result2.size();
    } else {
        res.phase2_triggered = false;
        res.result_factors = result.size();
    }

    auto tend = std::chrono::high_resolution_clock::now();
    res.total_ms = std::chrono::duration<double, std::milli>(tend - t0).count();
    return res;
}

void print_result(const std::string& label, const timing_result& r)
{
    std::cout << "  " << std::left << std::setw(22) << label;
    if (r.select_ms < 0) {
        std::cout << "  (prime invalid)" << std::endl;
        return;
    }
    std::cout << std::right << std::fixed << std::setprecision(3)
              << "  p=";
    if (r.prime > (1ULL << 32))
        std::cout << "2^" << (int)(log2((double)r.prime) + 0.5);
    else
        std::cout << r.prime;

    std::cout << "  r=" << r.r
              << "  a_h=" << r.a_h << "  a_mig=" << r.a_mig
              << "  | DDF=" << std::setw(7) << r.select_ms << "ms"
              << "  H1=" << std::setw(7) << r.hensel_p1_ms << "ms"
              << "  R1=" << std::setw(7) << r.recom_p1_ms << "ms";
    if (r.phase2_triggered) {
        std::cout << "  H2=" << std::setw(7) << r.hensel_p2_ms << "ms"
                  << "  R2=" << std::setw(7) << r.recom_p2_ms << "ms";
    }
    std::cout << "  | total=" << std::setw(7) << r.total_ms << "ms"
              << "  factors=" << r.result_factors
              << (r.phase2_triggered ? " [P2]" : "")
              << std::endl;
}

void run_experiment(const std::string& name, const upolynomial_<ZZ>& f, int repeats = 10)
{
    std::cout << "\n--- " << name << " (deg=" << get_deg(f)
              << ", terms=" << f.size() << ") ---" << std::endl;

    // 选几个代表性素数
    uint64_t primes[] = {
        7,                    // 典型小素数（__select_prime 常选中）
        101,                  // 中等素数
        65521,                // 16-bit 最大素数附近
        UINT64_MAX - 58       // 2^64 - 59
    };

    for (uint64_t p : primes) {
        // 多次运行取平均
        timing_result avg{};
        timing_result last{};
        int valid_runs = 0;

        for (int rep = 0; rep < repeats; ++rep) {
            last = run_with_prime(f, p);
            if (last.select_ms < 0) break;
            avg.select_ms += last.select_ms;
            avg.hensel_p1_ms += last.hensel_p1_ms;
            avg.recom_p1_ms += last.recom_p1_ms;
            avg.hensel_p2_ms += last.hensel_p2_ms;
            avg.recom_p2_ms += last.recom_p2_ms;
            avg.total_ms += last.total_ms;
            ++valid_runs;
        }
        if (valid_runs > 0) {
            avg.select_ms /= valid_runs;
            avg.hensel_p1_ms /= valid_runs;
            avg.recom_p1_ms /= valid_runs;
            avg.hensel_p2_ms /= valid_runs;
            avg.recom_p2_ms /= valid_runs;
            avg.total_ms /= valid_runs;
            avg.prime = last.prime;
            avg.r = last.r;
            avg.a_h = last.a_h;
            avg.a_mig = last.a_mig;
            avg.phase2_triggered = last.phase2_triggered;
            avg.result_factors = last.result_factors;
            print_result(std::to_string(p), avg);
        } else {
            std::cout << "  p=" << p << "  (invalid)" << std::endl;
        }
    }
}

int main()
{
    std::cout << "==== Phase 2 cost experiment: small vs large prime Hensel ====" << std::endl;

    // 用例 1: B2 大系数 (deg4, 4 factors)
    {
        upolynomial_<ZZ> f;
        f.push_back({umonomial(4), ZZ(96000)});
        f.push_back({umonomial(3), ZZ("38464001")});
        f.push_back({umonomial(2), ZZ("-61401631999")});
        f.push_back({umonomial(1), ZZ("-24616960640000")});
        f.push_back({umonomial(0), ZZ("-24555520640000")});
        run_experiment("B2 large-coeff (deg4)", f, 50);
    }

    // 用例 2: 更大的大系数多项式 — 模拟方向 5 的 alpha=500 代入
    // (x+500)(x-500)(x+1)(3x+7) 展开后系数到 ~10^6 量级
    {
        upolynomial_ZZ f1({{1, ZZ(1)}, {0, ZZ(500)}});
        upolynomial_ZZ f2({{1, ZZ(1)}, {0, ZZ(-500)}});
        upolynomial_ZZ f3({{1, ZZ(1)}, {0, ZZ(1)}});
        upolynomial_ZZ f4({{1, ZZ(3)}, {0, ZZ(7)}});
        auto f = f1 * f2 * f3 * f4;
        run_experiment("simulated alpha=500 (deg4)", f, 50);
    }

    // 用例 3: 更高次 — 6 个因子，大系数
    // (x^2+500x+1)(x^2-300x+7)(x+200)(x-400)(2x+999)(x^2+123x-456)
    {
        upolynomial_ZZ f1({{2, ZZ(1)}, {1, ZZ(500)}, {0, ZZ(1)}});
        upolynomial_ZZ f2({{2, ZZ(1)}, {1, ZZ(-300)}, {0, ZZ(7)}});
        upolynomial_ZZ f3({{1, ZZ(1)}, {0, ZZ(200)}});
        upolynomial_ZZ f4({{1, ZZ(1)}, {0, ZZ(-400)}});
        upolynomial_ZZ f5({{1, ZZ(2)}, {0, ZZ(999)}});
        upolynomial_ZZ f6({{2, ZZ(1)}, {1, ZZ(123)}, {0, ZZ(-456)}});
        auto f = f1 * f2 * f3 * f4 * f5 * f6;
        run_experiment("large-coeff 6 factors (deg10)", f, 20);
    }

    // 用例 4: ~deg15 典型小系数（对比基线）
    {
        upolynomial_ZZ uf1({{5, ZZ(1)}, {3, ZZ(-2)}, {1, ZZ(3)}, {0, ZZ(1)}});
        upolynomial_ZZ uf2({{4, ZZ(1)}, {2, ZZ(1)}, {0, ZZ(-1)}});
        upolynomial_ZZ uf3({{6, ZZ(2)}, {3, ZZ(-1)}, {1, ZZ(4)}, {0, ZZ(-3)}});
        auto f = uf1 * uf2 * uf3;
        run_experiment("small-coeff (deg15, baseline)", f, 20);
    }

    // 用例 5: Wilkinson W(15)
    {
        upolynomial_ZZ wilk({{1, ZZ(1)}, {0, ZZ(-1)}});
        for (int k = 2; k <= 15; ++k) {
            upolynomial_ZZ lin({{1, ZZ(1)}, {0, ZZ(-k)}});
            wilk = wilk * lin;
        }
        run_experiment("Wilkinson W(15)", wilk, 10);
    }

    return 0;
}
