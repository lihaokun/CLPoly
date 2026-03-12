/**
 * @file test_factorize_stress.cc
 * @brief 因式分解大规模随机压力测试
 *
 * 覆盖维度：系数范围 × 次数 × 因子数 × 变量数
 * 验证方式：乘积还原 + FLINT 交叉因子数对比
 *
 * Usage: test_factorize_stress [trial_scale] [-v]
 *   trial_scale: 试验数量倍率 (default 1)
 *   -v: verbose，输出每个失败用例的详情
 *
 * Build: make test/test_factorize_stress
 */
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include "clpoly_test.hh"
#include "crosscheck_flint.hh"
#include <cstdlib>
#include <iomanip>

using namespace clpoly;

// ---- 验证辅助 ----

// 单变量：乘积还原
static bool verify_upoly_product(const upolynomial_<ZZ>& f,
                                  const factorization<upolynomial_<ZZ>>& fac)
{
    upolynomial_<ZZ> prod;
    prod.push_back(std::make_pair(umonomial(0), fac.content));
    for (auto& [fi, ei] : fac.factors) {
        auto fi_pow = pow(fi, (int)ei);
        pair_vec_multiplies(prod.data(), prod.data(), fi_pow.data(), prod.comp());
        prod.normalization();
    }
    return prod == f;
}

// 多变量：乘积还原
static bool verify_mpoly_product(const polynomial_ZZ& f,
                                  const factorization<polynomial_ZZ>& fac)
{
    basic_monomial<grlex> m0;
    polynomial_ZZ prod;
    prod.push_back(std::make_pair(m0, fac.content));
    for (auto& [fi, ei] : fac.factors) {
        auto fi_pow = fi;
        for (uint64_t e = 1; e < ei; ++e) {
            fi_pow = fi_pow * fi;
            fi_pow.normalization();
        }
        prod = prod * fi_pow;
        prod.normalization();
    }
    return prod == f;
}

// 归一化因子列表：首一化（单变量 lc > 0），按因子排序
static std::vector<std::pair<upolynomial_ZZ, uint64_t>>
normalize_upoly_factors(const std::vector<std::pair<upolynomial_ZZ, uint64_t>>& factors)
{
    std::vector<std::pair<upolynomial_ZZ, uint64_t>> result;
    for (auto& [fi, ei] : factors) {
        auto f = fi;
        if (!f.empty() && f.front().second < 0) f = -f;
        result.push_back({std::move(f), ei});
    }
    std::sort(result.begin(), result.end());
    return result;
}

// FLINT 交叉对比：归一化因子集合一致
static bool crosscheck_upoly_factors(const upolynomial_<ZZ>& f,
                                      const factorization<upolynomial_<ZZ>>& cl_fac)
{
    auto fl_fac = crosscheck::flint_factor_upoly(f);

    // 归一化 CLPoly 因子
    std::vector<std::pair<upolynomial_ZZ, uint64_t>> cl_norm;
    for (auto& [fi, ei] : cl_fac.factors) {
        auto g = fi;
        if (!g.empty() && g.front().second < 0) g = -g;
        cl_norm.push_back({std::move(g), ei});
    }
    std::sort(cl_norm.begin(), cl_norm.end());

    // 归一化 FLINT 因子
    std::vector<std::pair<upolynomial_ZZ, uint64_t>> fl_norm;
    for (auto& [fi, ei] : fl_fac.factors) {
        auto g = fi;
        if (!g.empty() && g.front().second < 0) g = -g;
        fl_norm.push_back({std::move(g), (uint64_t)ei});
    }
    std::sort(fl_norm.begin(), fl_norm.end());

    return cl_norm == fl_norm;
}

static bool crosscheck_mpoly_factors(const polynomial_ZZ& f,
                                      const factorization<polynomial_ZZ>& cl_fac)
{
    auto vars = crosscheck::collect_vars(f);
    auto fl_fac = crosscheck::flint_factor_mpoly(f, vars);

    // 归一化 CLPoly
    std::vector<std::pair<polynomial_ZZ, uint64_t>> cl_norm;
    for (auto& [fi, ei] : cl_fac.factors) {
        auto g = fi;
        if (!g.empty() && g.front().second < 0) g = -g;
        cl_norm.push_back({std::move(g), ei});
    }
    std::sort(cl_norm.begin(), cl_norm.end());

    // 归一化 FLINT
    std::vector<std::pair<polynomial_ZZ, uint64_t>> fl_norm;
    for (auto& [fi, ei] : fl_fac.factors) {
        auto g = fi;
        if (!g.empty() && g.front().second < 0) g = -g;
        fl_norm.push_back({std::move(g), (uint64_t)ei});
    }
    std::sort(fl_norm.begin(), fl_norm.end());

    return cl_norm == fl_norm;
}

// ---- 统计 ----

struct Stats {
    int total = 0;
    int passed = 0;
    int product_fail = 0;
    int cross_fail = 0;
    int skipped = 0;
};

static void print_stats(const std::string& name, const Stats& s) {
    std::cout << "  " << std::left << std::setw(50) << name
              << std::right
              << "  total=" << std::setw(4) << s.total
              << "  pass=" << std::setw(4) << s.passed
              << "  prod_fail=" << s.product_fail
              << "  cross_fail=" << s.cross_fail
              << "  skip=" << s.skipped
              << std::endl;
    if (s.product_fail > 0 || s.cross_fail > 0) {
        std::cout << "  *** FAILURES DETECTED ***" << std::endl;
    }
}

int main(int argc, char* argv[])
{
    int trial_scale = 1;
    bool verbose = false;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-v") verbose = true;
        else trial_scale = std::max(1, atoi(argv[i]));
    }

    variable x("x"), y("y"), z("z"), w("w");
    int total_failures = 0;

    std::cout << "==== Factorize stress test (scale=" << trial_scale << ") ====" << std::endl;

    // ================================================================
    // §1. 单变量：大系数
    // ================================================================
    std::cout << "\n-- Univariate: large coefficients --" << std::endl;
    {
        int coeff_ranges[] = {100, 1000, 10000, 100000};
        int nfacs_list[] = {2, 3, 4};
        for (int C : coeff_ranges) {
            for (int nfac : nfacs_list) {
                Stats s;
                std::string name = "uni C=" + std::to_string(C) + " " + std::to_string(nfac) + "fac deg3";
                for (int trial = 0; trial < 20 * trial_scale; ++trial) {
                    upolynomial_<ZZ> f;
                    f.push_back({umonomial(0), ZZ(1)});
                    bool ok = true;
                    for (int i = 0; i < nfac; ++i) {
                        auto fi = random_upolynomial<ZZ>(3, 3, {-C, C});
                        if (fi.empty()) { ok = false; break; }
                        pair_vec_multiplies(f.data(), f.data(), fi.data(), f.comp());
                        f.normalization();
                    }
                    if (!ok || f.empty() || get_deg(f) < 2) { s.skipped++; s.total++; continue; }

                    s.total++;
                    auto cl_fac = factorize(f);

                    if (!verify_upoly_product(f, cl_fac)) {
                        s.product_fail++;
                        if (verbose) std::cerr << "PRODUCT FAIL: " << name << " trial " << trial << " f=" << f << std::endl;
                        continue;
                    }
                    if (!crosscheck_upoly_factors(f, cl_fac)) {
                        s.cross_fail++;
                        if (verbose) std::cerr << "NFACTOR FAIL: " << name << " trial " << trial << " f=" << f
                                               << " cl=" << cl_fac.factors.size() << std::endl;
                        continue;
                    }
                    s.passed++;
                }
                print_stats(name, s);
                total_failures += s.product_fail + s.cross_fail;
            }
        }
    }

    // ================================================================
    // §2. 单变量：高次
    // ================================================================
    std::cout << "\n-- Univariate: high degree --" << std::endl;
    {
        struct Cfg { int deg; int nfac; int C; };
        Cfg cfgs[] = {
            {5, 2, 10}, {5, 3, 10}, {5, 4, 10},
            {8, 2, 10}, {8, 3, 10},
            {10, 2, 10},
            {5, 2, 500}, {8, 2, 500},
        };
        for (auto& cfg : cfgs) {
            Stats s;
            std::string name = "uni deg" + std::to_string(cfg.deg) + " " +
                               std::to_string(cfg.nfac) + "fac C=" + std::to_string(cfg.C);
            for (int trial = 0; trial < 15 * trial_scale; ++trial) {
                upolynomial_<ZZ> f;
                f.push_back({umonomial(0), ZZ(1)});
                bool ok = true;
                for (int i = 0; i < cfg.nfac; ++i) {
                    auto fi = random_upolynomial<ZZ>(cfg.deg, cfg.deg, {-cfg.C, cfg.C});
                    if (fi.empty()) { ok = false; break; }
                    pair_vec_multiplies(f.data(), f.data(), fi.data(), f.comp());
                    f.normalization();
                }
                if (!ok || f.empty() || get_deg(f) < 2) { s.skipped++; s.total++; continue; }

                s.total++;
                auto cl_fac = factorize(f);

                if (!verify_upoly_product(f, cl_fac)) {
                    s.product_fail++;
                    if (verbose) std::cerr << "PRODUCT FAIL: " << name << " trial " << trial << std::endl;
                    continue;
                }
                if (!crosscheck_upoly_factors(f, cl_fac)) {
                    s.cross_fail++;
                    if (verbose) std::cerr << "NFACTOR FAIL: " << name << " trial " << trial
                                           << " cl=" << cl_fac.factors.size() << std::endl;
                    continue;
                }
                s.passed++;
            }
            print_stats(name, s);
            total_failures += s.product_fail + s.cross_fail;
        }
    }

    // ================================================================
    // §3. 单变量：重数
    // ================================================================
    std::cout << "\n-- Univariate: with multiplicity --" << std::endl;
    {
        Stats s;
        std::string name = "uni multiplicity (various)";
        for (int trial = 0; trial < 30 * trial_scale; ++trial) {
            auto f1 = random_upolynomial<ZZ>(3, 3, {-20, 20});
            auto f2 = random_upolynomial<ZZ>(3, 3, {-20, 20});
            if (f1.empty() || f2.empty()) { s.skipped++; s.total++; continue; }
            int e1 = 1 + (trial % 3);  // 1, 2, 3
            int e2 = 1 + ((trial / 3) % 2);  // 1, 2
            auto f = pow(f1, e1) * pow(f2, e2);
            if (f.empty() || get_deg(f) < 2) { s.skipped++; s.total++; continue; }

            s.total++;
            auto cl_fac = factorize(f);
            if (!verify_upoly_product(f, cl_fac)) {
                s.product_fail++;
                if (verbose) std::cerr << "PRODUCT FAIL: " << name << " trial " << trial << std::endl;
                continue;
            }
            if (!crosscheck_upoly_factors(f, cl_fac)) {
                s.cross_fail++;
                if (verbose) std::cerr << "NFACTOR FAIL: " << name << " trial " << trial
                                       << " cl=" << cl_fac.factors.size() << std::endl;
                continue;
            }
            s.passed++;
        }
        print_stats(name, s);
        total_failures += s.product_fail + s.cross_fail;
    }

    // ================================================================
    // §4. 双变量：大系数
    // ================================================================
    std::cout << "\n-- Bivariate: large coefficients --" << std::endl;
    {
        struct Cfg { int nfac; int deg; int C; };
        Cfg cfgs[] = {
            {2, 3, 100}, {2, 3, 500}, {2, 3, 1000},
            {3, 2, 100}, {3, 2, 500}, {3, 2, 1000},
            {4, 2, 100}, {4, 2, 500},
            {2, 5, 50},  {2, 5, 200},
        };
        for (auto& cfg : cfgs) {
            Stats s;
            std::string name = "bivar " + std::to_string(cfg.nfac) + "fac deg" +
                               std::to_string(cfg.deg) + " C=" + std::to_string(cfg.C);
            for (int trial = 0; trial < 15 * trial_scale; ++trial) {
                polynomial_ZZ f = polynomial_ZZ(1);
                bool ok = true;
                for (int i = 0; i < cfg.nfac; ++i) {
                    auto fi = random_polynomial<ZZ>({x, y}, cfg.deg, cfg.deg+1, {-cfg.C, cfg.C});
                    if (fi.empty()) { ok = false; break; }
                    f = f * fi;
                    f.normalization();
                }
                if (!ok || f.empty()) { s.skipped++; s.total++; continue; }

                s.total++;
                auto cl_fac = factorize(f);

                if (!verify_mpoly_product(f, cl_fac)) {
                    s.product_fail++;
                    if (verbose) std::cerr << "PRODUCT FAIL: " << name << " trial " << trial << std::endl;
                    continue;
                }
                if (!crosscheck_mpoly_factors(f, cl_fac)) {
                    s.cross_fail++;
                    if (verbose) std::cerr << "NFACTOR FAIL: " << name << " trial " << trial
                                           << " cl=" << cl_fac.factors.size() << std::endl;
                    continue;
                }
                s.passed++;
            }
            print_stats(name, s);
            total_failures += s.product_fail + s.cross_fail;
        }
    }

    // ================================================================
    // §5. 双变量：高次
    // ================================================================
    std::cout << "\n-- Bivariate: high degree --" << std::endl;
    {
        struct Cfg { int nfac; int deg; int C; };
        Cfg cfgs[] = {
            {2, 5, 5}, {2, 6, 5}, {2, 7, 5}, {2, 8, 5},
            {3, 4, 5}, {3, 5, 5},
            {2, 5, 100}, {2, 6, 100},
        };
        for (auto& cfg : cfgs) {
            Stats s;
            std::string name = "bivar " + std::to_string(cfg.nfac) + "fac deg" +
                               std::to_string(cfg.deg) + " C=" + std::to_string(cfg.C);
            for (int trial = 0; trial < 10 * trial_scale; ++trial) {
                polynomial_ZZ f = polynomial_ZZ(1);
                bool ok = true;
                for (int i = 0; i < cfg.nfac; ++i) {
                    auto fi = random_polynomial<ZZ>({x, y}, cfg.deg, cfg.deg+2, {-cfg.C, cfg.C});
                    if (fi.empty()) { ok = false; break; }
                    f = f * fi;
                    f.normalization();
                }
                if (!ok || f.empty()) { s.skipped++; s.total++; continue; }

                s.total++;
                auto cl_fac = factorize(f);

                if (!verify_mpoly_product(f, cl_fac)) {
                    s.product_fail++;
                    if (verbose) std::cerr << "PRODUCT FAIL: " << name << " trial " << trial << std::endl;
                    continue;
                }
                if (!crosscheck_mpoly_factors(f, cl_fac)) {
                    s.cross_fail++;
                    if (verbose) std::cerr << "NFACTOR FAIL: " << name << " trial " << trial
                                           << " cl=" << cl_fac.factors.size() << std::endl;
                    continue;
                }
                s.passed++;
            }
            print_stats(name, s);
            total_failures += s.product_fail + s.cross_fail;
        }
    }

    // ================================================================
    // §6. 三变量
    // ================================================================
    std::cout << "\n-- Trivariate --" << std::endl;
    {
        struct Cfg { int nfac; int deg; int C; };
        Cfg cfgs[] = {
            {2, 2, 5}, {2, 3, 5},
            {3, 2, 5}, {3, 2, 50},
            {4, 2, 3},
            {2, 2, 100}, {2, 3, 50},
        };
        for (auto& cfg : cfgs) {
            Stats s;
            std::string name = "trivar " + std::to_string(cfg.nfac) + "fac deg" +
                               std::to_string(cfg.deg) + " C=" + std::to_string(cfg.C);
            for (int trial = 0; trial < 10 * trial_scale; ++trial) {
                polynomial_ZZ f = polynomial_ZZ(1);
                bool ok = true;
                for (int i = 0; i < cfg.nfac; ++i) {
                    auto fi = random_polynomial<ZZ>({x, y, z}, cfg.deg, cfg.deg+2, {-cfg.C, cfg.C});
                    if (fi.empty()) { ok = false; break; }
                    f = f * fi;
                    f.normalization();
                }
                if (!ok || f.empty()) { s.skipped++; s.total++; continue; }

                s.total++;
                auto cl_fac = factorize(f);

                if (!verify_mpoly_product(f, cl_fac)) {
                    s.product_fail++;
                    if (verbose) std::cerr << "PRODUCT FAIL: " << name << " trial " << trial << std::endl;
                    continue;
                }
                if (!crosscheck_mpoly_factors(f, cl_fac)) {
                    s.cross_fail++;
                    if (verbose) std::cerr << "NFACTOR FAIL: " << name << " trial " << trial
                                           << " cl=" << cl_fac.factors.size() << std::endl;
                    continue;
                }
                s.passed++;
            }
            print_stats(name, s);
            total_failures += s.product_fail + s.cross_fail;
        }
    }

    // ================================================================
    // §7. 四变量
    // ================================================================
    std::cout << "\n-- 4-variable --" << std::endl;
    {
        struct Cfg { int nfac; int deg; int C; };
        Cfg cfgs[] = {
            {2, 2, 5}, {2, 2, 50},
            {3, 2, 3},
        };
        for (auto& cfg : cfgs) {
            Stats s;
            std::string name = "4var " + std::to_string(cfg.nfac) + "fac deg" +
                               std::to_string(cfg.deg) + " C=" + std::to_string(cfg.C);
            for (int trial = 0; trial < 5 * trial_scale; ++trial) {
                polynomial_ZZ f = polynomial_ZZ(1);
                bool ok = true;
                for (int i = 0; i < cfg.nfac; ++i) {
                    auto fi = random_polynomial<ZZ>({x, y, z, w}, cfg.deg, cfg.deg+2, {-cfg.C, cfg.C});
                    if (fi.empty()) { ok = false; break; }
                    f = f * fi;
                    f.normalization();
                }
                if (!ok || f.empty()) { s.skipped++; s.total++; continue; }

                s.total++;
                auto cl_fac = factorize(f);

                if (!verify_mpoly_product(f, cl_fac)) {
                    s.product_fail++;
                    if (verbose) std::cerr << "PRODUCT FAIL: " << name << " trial " << trial << std::endl;
                    continue;
                }
                if (!crosscheck_mpoly_factors(f, cl_fac)) {
                    s.cross_fail++;
                    if (verbose) std::cerr << "NFACTOR FAIL: " << name << " trial " << trial
                                           << " cl=" << cl_fac.factors.size() << std::endl;
                    continue;
                }
                s.passed++;
            }
            print_stats(name, s);
            total_failures += s.product_fail + s.cross_fail;
        }
    }

    // ================================================================
    // §8. 双变量：重数
    // ================================================================
    std::cout << "\n-- Bivariate: with multiplicity --" << std::endl;
    {
        Stats s;
        std::string name = "bivar multiplicity (various)";
        for (int trial = 0; trial < 20 * trial_scale; ++trial) {
            auto f1 = random_polynomial<ZZ>({x, y}, 2, 3, {-10, 10});
            auto f2 = random_polynomial<ZZ>({x, y}, 2, 3, {-10, 10});
            if (f1.empty() || f2.empty()) { s.skipped++; s.total++; continue; }
            int e1 = 1 + (trial % 3);
            int e2 = 1 + ((trial / 3) % 2);
            auto f = pow(f1, e1) * pow(f2, e2);
            f.normalization();
            if (f.empty()) { s.skipped++; s.total++; continue; }

            s.total++;
            auto cl_fac = factorize(f);
            if (!verify_mpoly_product(f, cl_fac)) {
                s.product_fail++;
                if (verbose) std::cerr << "PRODUCT FAIL: " << name << " trial " << trial << std::endl;
                continue;
            }
            if (!crosscheck_mpoly_factors(f, cl_fac)) {
                s.cross_fail++;
                if (verbose) std::cerr << "NFACTOR FAIL: " << name << " trial " << trial
                                       << " cl=" << cl_fac.factors.size() << std::endl;
                continue;
            }
            s.passed++;
        }
        print_stats(name, s);
        total_failures += s.product_fail + s.cross_fail;
    }

    // ================================================================
    // Summary
    // ================================================================
    std::cout << "\n==============================" << std::endl;
    if (total_failures == 0) {
        std::cout << "ALL PASSED" << std::endl;
    } else {
        std::cout << "TOTAL FAILURES: " << total_failures << std::endl;
    }
    return total_failures > 0 ? 1 : 0;
}
