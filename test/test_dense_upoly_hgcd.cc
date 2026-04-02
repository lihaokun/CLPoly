#include <clpoly/dense_upoly_zp.hh>
#include "clpoly_test.hh"
#include <random>

using namespace clpoly;

// Standalone Euclid GCD for reference
void euclid_gcd(dense_upoly_zp& G, const dense_upoly_zp& A, const dense_upoly_zp& B) {
    dense_upoly_zp a = A, b = B;
    if (a.deg() < b.deg()) std::swap(a, b);
    dense_upoly_zp q, r;
    while (!b.empty()) {
        dense_upoly_zp::divrem(q, r, a, b);
        a = std::move(b);
        b = std::move(r);
    }
    G = std::move(a);
}

// Generate random dense_upoly_zp of given degree
dense_upoly_zp make_random_poly(std::mt19937_64& rng, size_t deg, uint64_t p) {
    upolynomial_<Zp> sparse;
    for (int64_t i = deg; i >= 0; --i) {
        uint64_t c = rng() % p;
        if (c == 0 && i == (int64_t)deg) c = 1;
        if (c != 0)
            sparse.push_back({umonomial(i), Zp(c, p)});
    }
    return dense_upoly_zp(sparse, p);
}

int main() {
    uint64_t p = 18446744073709551557ULL;  // 2^64 - 59

    // ======== Basic HGCD correctness: coprime inputs ========
    CLPOLY_TEST("hgcd_gcd_coprime");
    {
        std::mt19937_64 rng(100);
        auto a = make_random_poly(rng, 200, p);
        auto b = make_random_poly(rng, 180, p);
        dense_upoly_zp g_euclid, g_hgcd;
        euclid_gcd(g_euclid, a, b);
        dense_upoly_zp::gcd(g_hgcd, a, b);
        CLPOLY_ASSERT_EQ(g_euclid.deg(), g_hgcd.deg());
    }

    // ======== HGCD with common factor (the original failing case) ========
    CLPOLY_TEST("hgcd_gcd_common_factor_deg170");
    {
        std::mt19937_64 rng(212);
        auto common = make_random_poly(rng, 170, p);
        auto f1 = make_random_poly(rng, 170, p);
        auto f2 = make_random_poly(rng, 170, p);
        dense_upoly_zp a, b;
        dense_upoly_zp::mul(a, common, f1);
        dense_upoly_zp::mul(b, common, f2);

        dense_upoly_zp g_euclid, g_hgcd;
        euclid_gcd(g_euclid, a, b);
        dense_upoly_zp::gcd(g_hgcd, a, b);

        CLPOLY_ASSERT_EQ(g_euclid.deg(), g_hgcd.deg());
        CLPOLY_ASSERT(g_hgcd.deg() >= 170);
    }

    // ======== Sweep over many seeds and common degrees ========
    CLPOLY_TEST("hgcd_gcd_sweep_500_seeds");
    {
        int failures = 0;
        int total = 0;
        for (int seed = 0; seed < 500; ++seed) {
            for (int common_deg : {10, 50, 100, 170}) {
                std::mt19937_64 rng(seed * 1000 + common_deg);
                auto common = make_random_poly(rng, common_deg, p);
                auto f1 = make_random_poly(rng, common_deg, p);
                auto f2 = make_random_poly(rng, common_deg, p);

                dense_upoly_zp a, b;
                dense_upoly_zp::mul(a, common, f1);
                dense_upoly_zp::mul(b, common, f2);

                dense_upoly_zp g_euclid, g_hgcd;
                euclid_gcd(g_euclid, a, b);
                dense_upoly_zp::gcd(g_hgcd, a, b);

                ++total;
                if (g_euclid.deg() != g_hgcd.deg()) {
                    ++failures;
                }
            }
        }
        CLPOLY_ASSERT_EQ(failures, 0);
    }

    // ======== GCD below HGCD cutoff (pure Euclid path) ========
    CLPOLY_TEST("hgcd_gcd_below_cutoff");
    {
        for (int deg = 5; deg < 100; deg += 15) {
            std::mt19937_64 rng(deg);
            auto common = make_random_poly(rng, deg, p);
            auto f1 = make_random_poly(rng, deg, p);
            auto f2 = make_random_poly(rng, deg, p);
            dense_upoly_zp a, b;
            dense_upoly_zp::mul(a, common, f1);
            dense_upoly_zp::mul(b, common, f2);

            dense_upoly_zp g_euclid, g_hgcd;
            euclid_gcd(g_euclid, a, b);
            dense_upoly_zp::gcd(g_hgcd, a, b);
            CLPOLY_ASSERT_EQ(g_euclid.deg(), g_hgcd.deg());
        }
    }

    // ======== GCD of coprime polynomials = constant ========
    CLPOLY_TEST("hgcd_gcd_coprime_result_constant");
    {
        std::mt19937_64 rng(42);
        auto a = make_random_poly(rng, 300, p);
        auto b = make_random_poly(rng, 280, p);
        dense_upoly_zp g;
        dense_upoly_zp::gcd(g, a, b);
        // Most random pairs are coprime
        CLPOLY_ASSERT_EQ(g.deg(), (int64_t)0);
    }

    // ======== GCD where one divides the other ========
    CLPOLY_TEST("hgcd_gcd_one_divides_other");
    {
        std::mt19937_64 rng(77);
        auto common = make_random_poly(rng, 200, p);
        auto f1 = make_random_poly(rng, 100, p);
        dense_upoly_zp a;
        dense_upoly_zp::mul(a, common, f1);

        dense_upoly_zp g;
        dense_upoly_zp::gcd(g, a, common);
        CLPOLY_ASSERT_EQ(g.deg(), common.deg());
    }

    // ======== GCD with equal inputs ========
    CLPOLY_TEST("hgcd_gcd_equal_inputs");
    {
        std::mt19937_64 rng(123);
        auto a = make_random_poly(rng, 250, p);
        dense_upoly_zp g;
        dense_upoly_zp::gcd(g, a, a);
        CLPOLY_ASSERT_EQ(g.deg(), a.deg());
    }

    // ======== GCD with high-degree inputs ========
    CLPOLY_TEST("hgcd_gcd_high_degree");
    {
        std::mt19937_64 rng(999);
        auto common = make_random_poly(rng, 300, p);
        auto f1 = make_random_poly(rng, 200, p);
        auto f2 = make_random_poly(rng, 200, p);
        dense_upoly_zp a, b;
        dense_upoly_zp::mul(a, common, f1);
        dense_upoly_zp::mul(b, common, f2);

        dense_upoly_zp g_euclid, g_hgcd;
        euclid_gcd(g_euclid, a, b);
        dense_upoly_zp::gcd(g_hgcd, a, b);
        CLPOLY_ASSERT_EQ(g_euclid.deg(), g_hgcd.deg());
    }

    // ======== GCD with small prime ========
    CLPOLY_TEST("hgcd_gcd_small_prime");
    {
        uint64_t sp = 1000000007ULL;
        std::mt19937_64 rng(555);
        auto make_sp = [&](size_t deg) -> dense_upoly_zp {
            upolynomial_<Zp> sparse;
            for (int64_t i = deg; i >= 0; --i) {
                uint64_t c = rng() % sp;
                if (c == 0 && i == (int64_t)deg) c = 1;
                if (c != 0)
                    sparse.push_back({umonomial(i), Zp(c, sp)});
            }
            return dense_upoly_zp(sparse, sp);
        };
        auto common = make_sp(150);
        auto f1 = make_sp(150);
        auto f2 = make_sp(150);
        dense_upoly_zp a, b;
        dense_upoly_zp::mul(a, common, f1);
        dense_upoly_zp::mul(b, common, f2);

        dense_upoly_zp g_euclid, g_hgcd;
        euclid_gcd(g_euclid, a, b);
        dense_upoly_zp::gcd(g_hgcd, a, b);
        CLPOLY_ASSERT_EQ(g_euclid.deg(), g_hgcd.deg());
    }

    // ======== GCD at HGCD cutoff boundary ========
    CLPOLY_TEST("hgcd_gcd_cutoff_boundary");
    {
        // Test degrees near GCD_CUTOFF (340) and HGCD_CUTOFF (100)
        for (int common_deg : {80, 100, 120, 170, 200}) {
            std::mt19937_64 rng(common_deg * 7);
            auto common = make_random_poly(rng, common_deg, p);
            auto f1 = make_random_poly(rng, common_deg, p);
            auto f2 = make_random_poly(rng, common_deg, p);
            dense_upoly_zp a, b;
            dense_upoly_zp::mul(a, common, f1);
            dense_upoly_zp::mul(b, common, f2);

            dense_upoly_zp g_euclid, g_hgcd;
            euclid_gcd(g_euclid, a, b);
            dense_upoly_zp::gcd(g_hgcd, a, b);
            CLPOLY_ASSERT_EQ(g_euclid.deg(), g_hgcd.deg());
        }
    }

    return clpoly_test::test_summary();
}
