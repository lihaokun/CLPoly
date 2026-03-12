/**
 * @file test_factorize_trace.cc
 * @brief 深入追踪 Case 1 不完全分解 bug 的根因
 *
 * Case 1: f = 10707840*x^12 - 4358016*x^11 - ... + 4287360*x
 * CLPoly: 4 factors (deg6 未分裂)
 * FLINT:  5 factors (两个 deg3)
 */
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <clpoly/polynomial_factorize_zp.hh>
#include <iostream>
#include <iomanip>

using namespace clpoly;

// 复制自 polynomial_factorize_univar.hh 的内部函数声明
// 为了调试，我们需要直接调用这些函数

int main()
{
    // ===== 构造 Case 1 多项式 =====
    upolynomial_<ZZ> f;
    f.push_back({umonomial(12), ZZ("10707840")});
    f.push_back({umonomial(11), ZZ("-4358016")});
    f.push_back({umonomial(10), ZZ("-17120232")});
    f.push_back({umonomial(9), ZZ("-1940088")});
    f.push_back({umonomial(8), ZZ("3144483")});
    f.push_back({umonomial(7), ZZ("-20309937")});
    f.push_back({umonomial(6), ZZ("-15166461")});
    f.push_back({umonomial(5), ZZ("11473425")});
    f.push_back({umonomial(4), ZZ("-2040858")});
    f.push_back({umonomial(3), ZZ("-7493748")});
    f.push_back({umonomial(2), ZZ("6025272")});
    f.push_back({umonomial(1), ZZ("4287360")});

    std::cout << "f = " << f << std::endl;
    std::cout << "deg(f) = " << get_deg(f) << ", terms = " << f.size() << std::endl;
    std::cout << "lc(f) = " << f.front().second << std::endl;

    // ===== Step 1: 提取内容和本原化 =====
    ZZ ct = cont(f);
    std::cout << "\n=== Step 1: content & primitive ===\n";
    std::cout << "content(f) = " << ct << std::endl;

    upolynomial_<ZZ> f_prim = f;
    if (f_prim.front().second < ZZ(0)) ct = -ct;
    for (auto& term : f_prim)
        term.second /= ct;
    std::cout << "f_prim = " << f_prim << std::endl;
    std::cout << "lc(f_prim) = " << f_prim.front().second << std::endl;

    // ===== Step 2: 无平方分解 =====
    std::cout << "\n=== Step 2: squarefree factorization ===\n";
    variable __x("x");
    polynomial_<ZZ,lex> poly_prim;
    poly_convert(f_prim, poly_prim, __x);
    auto sqf = squarefreefactorize(poly_prim);
    std::cout << "Number of squarefree factors: " << sqf.size() << std::endl;
    for (auto& [sqf_fac, mult] : sqf) {
        upolynomial_<ZZ> usqf;
        poly_convert(sqf_fac, usqf);
        std::cout << "  mult=" << mult << " deg=" << get_deg(usqf) << " : " << usqf << std::endl;
    }

    // ===== Step 3: 素数选择（手动模拟 __select_prime）=====
    std::cout << "\n=== Step 3: prime selection (manual trace) ===\n";

    // 对本原无平方部分做素数选择
    // 如果 sqf 只有一个因子，就用那个
    upolynomial_<ZZ> f_sqfree;
    for (auto& [sqf_fac, mult] : sqf) {
        if (is_number(sqf_fac)) continue;
        poly_convert(sqf_fac, f_sqfree);
    }
    if (f_sqfree.empty()) f_sqfree = f_prim;

    std::cout << "f_sqfree = " << f_sqfree << std::endl;
    std::cout << "deg(f_sqfree) = " << get_deg(f_sqfree) << std::endl;
    std::cout << "lc(f_sqfree) = " << f_sqfree.front().second << std::endl;

    ZZ lc_f = f_sqfree.front().second;
    int64_t deg_f = get_deg(f_sqfree);

    // 逐素数追踪
    size_t max_tries = 3;
    std::mt19937 rng(42);

    struct prime_info {
        uint64_t p;
        size_t nfactors;
        std::vector<int> factor_degs;
    };
    std::vector<prime_info> tried_primes;

    std::cout << "\nScanning primes:\n";
    for (uint64_t p = 2, tried = 0; tried < max_tries + 10 && p < 200; p = next_prime_64(p))
    {
        // lc mod p
        ZZ lc_mod;
        ZZ::fdiv_r(lc_mod, lc_f, ZZ(p));
        if (!lc_mod) {
            std::cout << "  p=" << p << ": skip (lc ≡ 0 mod p)\n";
            continue;
        }

        // f mod p
        auto fp = polynomial_mod(f_sqfree, p);
        if (fp.empty() || get_deg(fp) != deg_f) {
            std::cout << "  p=" << p << ": skip (deg drops)\n";
            continue;
        }

        // squarefree check
        auto fp_deriv = derivative(fp);
        if (fp_deriv.empty()) {
            std::cout << "  p=" << p << ": skip (f'=0)\n";
            continue;
        }
        auto g = polynomial_GCD(fp, fp_deriv);
        if (get_deg(g) > 0) {
            std::cout << "  p=" << p << ": skip (not squarefree, gcd deg=" << get_deg(g) << ")\n";
            continue;
        }

        ++tried;

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

        prime_info pi;
        pi.p = p;
        pi.nfactors = irr_factors.size();
        for (auto& fi : irr_factors)
            pi.factor_degs.push_back(get_deg(fi));
        std::sort(pi.factor_degs.begin(), pi.factor_degs.end());

        std::cout << "  p=" << p << " (tried #" << tried << "): r=" << pi.nfactors << ", degs=[";
        for (size_t i = 0; i < pi.factor_degs.size(); ++i) {
            if (i > 0) std::cout << ",";
            std::cout << pi.factor_degs[i];
        }
        std::cout << "]";

        tried_primes.push_back(pi);

        // Mark the one that __select_prime would choose
        if (tried <= max_tries)
            std::cout << " <-- tried by __select_prime";
        std::cout << std::endl;
    }

    // ===== Step 4: 确认 __select_prime 的实际选择 =====
    // 找 tried_primes 前3个中 nfactors 最小的
    if (!tried_primes.empty()) {
        size_t best_idx = 0;
        for (size_t i = 1; i < std::min((size_t)3, tried_primes.size()); ++i) {
            if (tried_primes[i].nfactors < tried_primes[best_idx].nfactors)
                best_idx = i;
        }
        auto& best = tried_primes[best_idx];
        std::cout << "\n__select_prime would choose: p=" << best.p
                  << ", r=" << best.nfactors << ", degs=[";
        for (size_t i = 0; i < best.factor_degs.size(); ++i) {
            if (i > 0) std::cout << ",";
            std::cout << best.factor_degs[i];
        }
        std::cout << "]\n";

        // ===== Step 5: Hensel 精度计算 =====
        std::cout << "\n=== Step 5: Hensel precision ===\n";
        int r = best.nfactors;
        uint64_t p = best.p;
        double logp = std::log((double)p);
        int min_b = (int)ZZ(p).sizeinbase(2);
        int N = (int)f_sqfree.size() - 1;
        double a_h_d = std::ceil(
            (2.5 * r + min_b) * std::log(2.0) / logp
            + std::log((double)(N + 1)) / (2.0 * logp));
        int a_h_raw = std::max(1, (int)a_h_d);

        // Mignotte bound
        auto binomial = [](int64_t n, int64_t k) -> ZZ {
            if (k < 0 || k > n) return ZZ(0);
            if (k == 0 || k == n) return ZZ(1);
            if (k > n - k) k = n - k;
            ZZ result(1);
            for (int64_t i = 0; i < k; ++i) { result *= ZZ(n - i); result /= ZZ(i + 1); }
            return result;
        };
        auto upoly_norm_l2_sq = [](const upolynomial_<ZZ>& ff) -> ZZ {
            ZZ s(0); for (auto& t : ff) s += t.second * t.second; return s;
        };
        auto isqrt_ceil = [](const ZZ& n) -> ZZ {
            if (n <= ZZ(0)) return ZZ(0);
            size_t bits = n.sizeinbase(2);
            ZZ x(1); x <<= (bits + 1) / 2;
            while (true) { ZZ x1 = (x + n / x) / ZZ(2); if (x1 >= x) break; x = x1; }
            if (x * x < n) x += 1;
            return x;
        };

        ZZ B_mig = binomial(deg_f, deg_f / 2) * isqrt_ceil(upoly_norm_l2_sq(f_sqfree));
        ZZ lc_abs = lc_f; if (lc_abs < ZZ(0)) lc_abs = -lc_abs;
        ZZ target = ZZ(2) * lc_abs * B_mig;
        int a_mig = 0;
        ZZ pa(1);
        while (pa <= target) { pa *= ZZ(p); ++a_mig; }
        int a_h = std::min(a_mig, a_h_raw);

        std::cout << "r = " << r << ", p = " << p << std::endl;
        std::cout << "B_mig = " << B_mig << std::endl;
        std::cout << "target = 2*|lc|*B_mig = " << target << std::endl;
        std::cout << "a_h_raw = " << a_h_raw << ", a_mig = " << a_mig << std::endl;
        std::cout << "a_h (used for Phase 1) = " << a_h << std::endl;
        std::cout << "Phase 2 trigger condition: result.size() < " << r
                  << " AND " << a_h << " < " << a_mig
                  << " => " << (a_h < a_mig ? "CAN trigger" : "CANNOT trigger (a_h == a_mig)")
                  << std::endl;
    }

    // ===== Step 6: 验证真因子在模 p 下的分解 =====
    std::cout << "\n=== Step 6: true factors mod p ===\n";
    // 真因子（FLINT 的结果）
    upolynomial_<ZZ> h0, h1, h2, h3, h4;
    // x
    h0.push_back({umonomial(1), ZZ(1)});
    // 65x^2-61x-80
    h1.push_back({umonomial(2), ZZ(65)});
    h1.push_back({umonomial(1), ZZ(-61)});
    h1.push_back({umonomial(0), ZZ(-80)});
    // 24x^3+9x+14
    h2.push_back({umonomial(3), ZZ(24)});
    h2.push_back({umonomial(1), ZZ(9)});
    h2.push_back({umonomial(0), ZZ(14)});
    // 44x^3-13x^2+29
    h3.push_back({umonomial(3), ZZ(44)});
    h3.push_back({umonomial(2), ZZ(-13)});
    h3.push_back({umonomial(0), ZZ(29)});
    // 52x^3+43x^2-44
    h4.push_back({umonomial(3), ZZ(52)});
    h4.push_back({umonomial(2), ZZ(43)});
    h4.push_back({umonomial(0), ZZ(-44)});

    // 对每个 tried prime，显示真因子的模分解
    for (auto& pi : tried_primes) {
        if (pi.p > 50) break;  // 只看小素数
        uint64_t p = pi.p;
        std::cout << "\np=" << p << ":\n";
        auto show_factor_mod_p = [&](const std::string& name, const upolynomial_<ZZ>& h) {
            auto hp = polynomial_mod(h, p);
            __upoly_make_monic(hp);
            auto ddf = __ddf_Zp(hp);
            std::cout << "  " << name << " mod " << p << ": ";
            for (auto& [gk, dk] : ddf) {
                std::cout << dk << "^" << (get_deg(gk)/dk) << " ";
            }
            std::cout << "(h mod p = " << hp << ")" << std::endl;
        };
        show_factor_mod_p("h0=x", h0);
        show_factor_mod_p("h1=65x^2-61x-80", h1);
        show_factor_mod_p("h2=24x^3+9x+14", h2);
        show_factor_mod_p("h3=44x^3-13x^2+29", h3);
        show_factor_mod_p("h4=52x^3+43x^2-44", h4);

        // 计算总模因子数
        int total_mod_factors = 0;
        auto count_mf = [&](const upolynomial_<ZZ>& h) {
            auto hp = polynomial_mod(h, p);
            __upoly_make_monic(hp);
            auto ddf = __ddf_Zp(hp);
            for (auto& [gk, dk] : ddf) total_mod_factors += get_deg(gk)/dk;
        };
        count_mf(h0); count_mf(h1); count_mf(h2); count_mf(h3); count_mf(h4);
        std::cout << "  Total mod factors (from true factors): " << total_mod_factors
                  << " (reported r=" << pi.nfactors << ")" << std::endl;
    }

    // ===== Step 7: 实际调用 factorize =====
    std::cout << "\n=== Step 7: actual factorize result ===\n";
    auto result = factorize(f);
    std::cout << "content = " << result.content << std::endl;
    std::cout << "factors (" << result.factors.size() << "):" << std::endl;
    for (auto& [fi, ei] : result.factors)
        std::cout << "  [e=" << ei << "] deg=" << get_deg(fi) << " : " << fi << std::endl;

    // 验证乘积
    upolynomial_<ZZ> prod;
    prod.push_back(std::make_pair(umonomial(0), result.content));
    for (auto& [fi, ei] : result.factors) {
        auto fi_pow = pow(fi, (int)ei);
        pair_vec_multiplies(prod.data(), prod.data(), fi_pow.data(), prod.comp());
        prod.normalization();
    }
    std::cout << "Product == f? " << (prod == f ? "YES" : "NO") << std::endl;

    // ===== 验证 deg6 是否可分 =====
    for (auto& [fi, ei] : result.factors) {
        if (get_deg(fi) == 6) {
            std::cout << "\n=== Checking deg6 factor ===\n";
            std::cout << "deg6 = " << fi << std::endl;

            // 尝试 factorize 这个 deg6
            std::cout << "Re-factorizing deg6 factor...\n";
            auto re_result = factorize(fi);
            std::cout << "Re-factorize content = " << re_result.content << std::endl;
            std::cout << "Re-factorize factors (" << re_result.factors.size() << "):" << std::endl;
            for (auto& [gi, ge] : re_result.factors)
                std::cout << "  [e=" << ge << "] deg=" << get_deg(gi) << " : " << gi << std::endl;

            // 尝试除以已知的因子
            std::cout << "\nTrial divide deg6 by h3 (44x^3-13x^2+29):" << std::endl;
            upolynomial_<ZZ> q, rem;
            pair_vec_div(q.data(), rem.data(), fi.data(), h3.data(), fi.comp());
            std::cout << "  q = " << q << std::endl;
            std::cout << "  r = " << rem << std::endl;
            std::cout << "  Divides? " << (rem.empty() ? "YES" : "NO") << std::endl;
        }
    }

    return 0;
}
