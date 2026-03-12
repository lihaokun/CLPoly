#include <clpoly/clpoly.hh>
#include <chrono>
#include <iostream>
#include <random>
#include <iomanip>

using namespace clpoly;

// 旧实现：稀疏 Euclid GCD（从 polynomial_gcd.hh 复制的原始代码）
static int64_t sparse_polynomial_GCD(upolynomial_<Zp>& Pout,
                              const upolynomial_<Zp>& G,
                              const upolynomial_<Zp>& F,
                              const Zp& Lc_gcd,
                              int64_t deg)
{
    assert(!F.empty() && !G.empty());
    upolynomial_<Zp> Pout_, Pout_1, Pout_2;
    Pout = G;
    pair_vec_div(Pout_2.data(), Pout_1.data(), F.data(), Pout.data(), Pout.comp());
    while (!Pout_1.empty()) {
        std::swap(Pout.data(), Pout_.data());
        std::swap(Pout.data(), Pout_1.data());
        pair_vec_div(Pout_2.data(), Pout_1.data(), Pout_.data(), Pout.data(), Pout.comp());
    }
    Zp lc_inv = Pout.front().second.inv() * Lc_gcd;
    for (auto& i : Pout)
        i.second *= lc_inv;
    int64_t deg_ = Pout.front().first.deg();
    return (deg_ <= deg) ? deg_ : -1;
}

// 新实现：稠密 Euclid GCD（当前 __polynomial_GCD）
static int64_t dense_polynomial_GCD(upolynomial_<Zp>& Pout,
                             const upolynomial_<Zp>& G,
                             const upolynomial_<Zp>& F,
                             const Zp& Lc_gcd,
                             int64_t deg)
{
    assert(!F.empty() && !G.empty());
    uint64_t p = F.front().second.prime();
    dense_upoly_zp f_d(F, p), g_d(G, p);
    dense_upoly_zp result;
    dense_upoly_zp::gcd(result, f_d, g_d);
    uint64_t scale = result.nmod_mul(
        result.nmod_inv(result.lead()), Lc_gcd.number());
    result.scalar_mul(scale);
    int64_t d = result.deg();
    if (d > deg) return -1;
    Pout = result.to_upoly();
    return d;
}

// 生成随机 Zp 多项式（稠密，度 = deg）
static upolynomial_<Zp> random_upoly_zp(int64_t deg, uint64_t prime, std::mt19937_64& rng)
{
    upolynomial_<Zp> p;
    std::uniform_int_distribution<uint64_t> dist(1, prime - 1);
    for (int64_t i = deg; i >= 0; --i) {
        uint64_t c = dist(rng);
        p.push_back({umonomial(i), Zp(c, prime)});
    }
    return p;
}

// 生成有公因子的测试用例: F = H*A, G = H*B
struct TestCase {
    upolynomial_<Zp> F, G;
    int64_t common_deg;
    int64_t total_deg_f, total_deg_g;
};

static TestCase make_test(int64_t deg_a, int64_t deg_b, int64_t deg_common,
                          uint64_t prime, std::mt19937_64& rng)
{
    auto H = random_upoly_zp(deg_common, prime, rng);
    auto A = random_upoly_zp(deg_a, prime, rng);
    auto B = random_upoly_zp(deg_b, prime, rng);

    // F = H * A, G = H * B (用稠密乘法)
    dense_upoly_zp h_d(H, prime), a_d(A, prime), b_d(B, prime);
    dense_upoly_zp f_d, g_d;
    dense_upoly_zp::mul(f_d, h_d, a_d);
    dense_upoly_zp::mul(g_d, h_d, b_d);

    TestCase tc;
    tc.F = f_d.to_upoly();
    tc.G = g_d.to_upoly();
    tc.common_deg = deg_common;
    tc.total_deg_f = deg_a + deg_common;
    tc.total_deg_g = deg_b + deg_common;
    return tc;
}

struct BenchResult {
    double sparse_us;
    double dense_us;
};

static BenchResult bench_one(const TestCase& tc, uint64_t prime, int repeats)
{
    Zp lc_gcd(1, prime);
    int64_t deg = tc.common_deg;

    // 预热
    {
        upolynomial_<Zp> out;
        dense_polynomial_GCD(out, tc.F, tc.G, lc_gcd, deg + 10);
    }

    // 稀疏
    double sparse_total = 0;
    for (int r = 0; r < repeats; ++r) {
        upolynomial_<Zp> out;
        auto t0 = std::chrono::high_resolution_clock::now();
        sparse_polynomial_GCD(out, tc.F, tc.G, lc_gcd, deg + 10);
        auto t1 = std::chrono::high_resolution_clock::now();
        sparse_total += std::chrono::duration<double, std::micro>(t1 - t0).count();
    }

    // 稠密
    double dense_total = 0;
    for (int r = 0; r < repeats; ++r) {
        upolynomial_<Zp> out;
        auto t0 = std::chrono::high_resolution_clock::now();
        dense_polynomial_GCD(out, tc.F, tc.G, lc_gcd, deg + 10);
        auto t1 = std::chrono::high_resolution_clock::now();
        dense_total += std::chrono::duration<double, std::micro>(t1 - t0).count();
    }

    return {sparse_total / repeats, dense_total / repeats};
}

int main()
{
    uint64_t prime = UINT64_C(18446744073709551557); // 2^64 - 59
    std::mt19937_64 rng(42);

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "=== dense_upoly_zp vs sparse upolynomial_<Zp>: Zp Euclid GCD ===\n";
    std::cout << "Prime: 2^64 - 59\n\n";

    std::cout << std::setw(40) << std::left << "用例"
              << std::setw(14) << std::right << "稀疏(us)"
              << std::setw(14) << "稠密(us)"
              << std::setw(10) << "加速"
              << "\n";
    std::cout << std::string(78, '-') << "\n";

    struct Config {
        const char* name;
        int64_t deg_a, deg_b, deg_common;
        int repeats;
    };

    Config configs[] = {
        {"deg100+common50",       50,  50,  50,  500},
        {"deg200+common100",     100, 100, 100,  200},
        {"deg500+common250",     250, 250, 250,   50},
        {"deg1000+common500",    500, 500, 500,   20},
        {"deg2000+common1000",  1000,1000,1000,    5},
        {"deg100 coprime",       100, 100,   0, 1000},
        {"deg500 coprime",       500, 500,   0,  100},
        {"deg1000 coprime",     1000,1000,   0,   20},
        {"deg2000 coprime",     2000,2000,   0,    5},
    };

    for (auto& cfg : configs) {
        auto tc = make_test(cfg.deg_a, cfg.deg_b, cfg.deg_common, prime, rng);
        auto res = bench_one(tc, prime, cfg.repeats);
        double speedup = res.sparse_us / res.dense_us;

        std::cout << std::setw(40) << std::left << cfg.name
                  << std::setw(12) << std::right << res.sparse_us
                  << std::setw(14) << res.dense_us
                  << std::setw(8) << speedup << "x"
                  << "\n";
    }

    std::cout << "\n";
    return 0;
}
