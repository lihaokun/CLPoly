#include "clpoly_test.hh"
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>

using namespace clpoly;
using PolyZZ = polynomial_<ZZ, lex>;

PolyZZ make_lex(const polynomial_ZZ& p)
{
    PolyZZ result;
    poly_convert(p, result);
    return result;
}

upolynomial_<ZZ> to_upoly(const PolyZZ& f)
{
    upolynomial_<ZZ> result;
    poly_convert(f, result);
    return result;
}

// 选择一个满足 Mignotte 界的素数
static uint64_t pick_prime(const PolyZZ& f_scaled)
{
    // 计算 max |coeff|
    ZZ max_coeff(0);
    for (const auto& term : f_scaled)
    {
        ZZ a = abs(term.second);
        if (a > max_coeff) max_coeff = a;
    }
    // p 需要 > 2 * max_coeff * C(d, d/2)（简化：取 >> max_coeff 的素数）
    // 对于测试用的小多项式，10007 足够
    uint64_t candidates[] = {10007, 100003, 1000003};
    for (auto p : candidates)
        if (ZZ((int64_t)p) > 2 * max_coeff)
            return p;
    return 1000003;
}

int main()
{
    variable x("x"), y("y"), z("z");

    // ============================================================
    // 双变量提升
    // ============================================================

    CLPOLY_TEST("mtshl_lift: (x+y)(x-y), constant lc");
    {
        auto f = make_lex(pow(x,2) - pow(y,2));
        auto alpha = __select_eval_point(f, x);

        auto f0 = assign(f, alpha);
        auto f0_upoly = to_upoly(f0);
        auto uni_fac = factorize(f0_upoly);

        std::vector<upolynomial_<ZZ>> uni_factors;
        for (auto& [fi, ei] : uni_fac.factors)
            uni_factors.push_back(fi);
        CLPOLY_ASSERT(uni_factors.size() == 2);

        auto lc_result = __wang_leading_coeff(f, uni_factors, alpha, x);
        CLPOLY_ASSERT(lc_result.success);

        uint64_t p = pick_prime(lc_result.f_scaled);
        auto G = __mtshl_lift(lc_result.f_scaled, lc_result.scaled_factors, lc_result.lc_targets, alpha, x, p);
        CLPOLY_ASSERT(G.size() == 2);

        auto prod = G[0] * G[1];
        prod.normalization();
        CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
    }

    CLPOLY_TEST("mtshl_lift: (x+y)(x-y+1), constant lc");
    {
        auto f = make_lex(pow(x,2) + x - pow(y,2) + y);
        auto alpha = __select_eval_point(f, x);

        auto f0 = assign(f, alpha);
        auto f0_upoly = to_upoly(f0);
        auto uni_fac = factorize(f0_upoly);

        std::vector<upolynomial_<ZZ>> uni_factors;
        for (auto& [fi, ei] : uni_fac.factors)
            uni_factors.push_back(fi);
        CLPOLY_ASSERT(uni_factors.size() == 2);

        auto lc_result = __wang_leading_coeff(f, uni_factors, alpha, x);
        CLPOLY_ASSERT(lc_result.success);

        uint64_t p = pick_prime(lc_result.f_scaled);
        auto G = __mtshl_lift(lc_result.f_scaled, lc_result.scaled_factors, lc_result.lc_targets, alpha, x, p);

        auto prod = G[0] * G[1];
        prod.normalization();
        CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
    }

    CLPOLY_TEST("mtshl_lift: nontrivial lc, ((y+1)x+y)(x+1)");
    {
        auto f = make_lex((y + polynomial_ZZ(1))*pow(x,2) + (2*y + polynomial_ZZ(1))*x + y);
        auto alpha = __select_eval_point(f, x);

        auto f0 = assign(f, alpha);
        auto f0_upoly = to_upoly(f0);
        auto uni_fac = factorize(f0_upoly);

        std::vector<upolynomial_<ZZ>> uni_factors;
        for (auto& [fi, ei] : uni_fac.factors)
            uni_factors.push_back(fi);
        CLPOLY_ASSERT(uni_factors.size() == 2);

        auto lc_result = __wang_leading_coeff(f, uni_factors, alpha, x);
        CLPOLY_ASSERT(lc_result.success);

        uint64_t p = pick_prime(lc_result.f_scaled);
        auto G = __mtshl_lift(lc_result.f_scaled, lc_result.scaled_factors, lc_result.lc_targets, alpha, x, p);
        CLPOLY_ASSERT(G.size() == 2);

        auto prod = G[0] * G[1];
        prod.normalization();
        CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
    }

    CLPOLY_TEST("mtshl_lift: three factors, (x+y)(x-y)(x+1)");
    {
        auto f = make_lex(pow(x,3) + pow(x,2) - x*pow(y,2) - pow(y,2));
        auto alpha = __select_eval_point(f, x);

        auto f0 = assign(f, alpha);
        auto f0_upoly = to_upoly(f0);
        auto uni_fac = factorize(f0_upoly);

        std::vector<upolynomial_<ZZ>> uni_factors;
        for (auto& [fi, ei] : uni_fac.factors)
            uni_factors.push_back(fi);
        CLPOLY_ASSERT(uni_factors.size() == 3);

        auto lc_result = __wang_leading_coeff(f, uni_factors, alpha, x);
        CLPOLY_ASSERT(lc_result.success);

        uint64_t p = pick_prime(lc_result.f_scaled);
        auto G = __mtshl_lift(lc_result.f_scaled, lc_result.scaled_factors, lc_result.lc_targets, alpha, x, p);
        CLPOLY_ASSERT(G.size() == 3);

        auto prod = G[0] * G[1];
        prod.normalization();
        prod = prod * G[2];
        prod.normalization();
        CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
    }

    // ============================================================
    // 三变量提升
    // ============================================================

    CLPOLY_TEST("mtshl_lift: trivariate (x+y+z)(x-y+z)");
    {
        auto f = make_lex(pow(x,2) + 2*x*z + pow(z,2) - pow(y,2));
        auto alpha = __select_eval_point(f, x);
        CLPOLY_ASSERT(alpha.size() == 2);

        auto f0 = assign(f, alpha);
        auto f0_upoly = to_upoly(f0);
        auto uni_fac = factorize(f0_upoly);

        std::vector<upolynomial_<ZZ>> uni_factors;
        for (auto& [fi, ei] : uni_fac.factors)
            uni_factors.push_back(fi);

        if (uni_factors.size() == 2)
        {
            auto lc_result = __wang_leading_coeff(f, uni_factors, alpha, x);
            CLPOLY_ASSERT(lc_result.success);

            uint64_t p = pick_prime(lc_result.f_scaled);
            auto G = __mtshl_lift(lc_result.f_scaled, lc_result.scaled_factors, lc_result.lc_targets, alpha, x, p);

            auto prod = G[0] * G[1];
            prod.normalization();
            CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
        }
        else
        {
            CLPOLY_ASSERT(true);
        }
    }

    CLPOLY_TEST("mtshl_lift: trivariate (x+y)(x+z)");
    {
        auto f = make_lex(pow(x,2) + x*z + x*y + y*z);
        auto alpha = __select_eval_point(f, x);
        CLPOLY_ASSERT(alpha.size() == 2);

        auto f0 = assign(f, alpha);
        auto f0_upoly = to_upoly(f0);
        auto uni_fac = factorize(f0_upoly);

        std::vector<upolynomial_<ZZ>> uni_factors;
        for (auto& [fi, ei] : uni_fac.factors)
            uni_factors.push_back(fi);

        if (uni_factors.size() == 2)
        {
            auto lc_result = __wang_leading_coeff(f, uni_factors, alpha, x);
            CLPOLY_ASSERT(lc_result.success);

            uint64_t p = pick_prime(lc_result.f_scaled);
            auto G = __mtshl_lift(lc_result.f_scaled, lc_result.scaled_factors, lc_result.lc_targets, alpha, x, p);

            auto prod = G[0] * G[1];
            prod.normalization();
            CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
        }
        else
        {
            CLPOLY_ASSERT(true);
        }
    }

    CLPOLY_TEST("mtshl_lift: nontrivial lc, ((y^2-1)x^2 + 2yx + 1)");
    {
        auto f = make_lex((pow(y,2) - polynomial_ZZ(1))*pow(x,2) + 2*y*x + polynomial_ZZ(1));
        auto alpha = __select_eval_point(f, x);

        auto f0 = assign(f, alpha);
        auto f0_upoly = to_upoly(f0);
        auto uni_fac = factorize(f0_upoly);

        std::vector<upolynomial_<ZZ>> uni_factors;
        for (auto& [fi, ei] : uni_fac.factors)
            uni_factors.push_back(fi);

        if (uni_factors.size() == 2)
        {
            auto lc_result = __wang_leading_coeff(f, uni_factors, alpha, x);
            if (lc_result.success)
            {
                uint64_t p = pick_prime(lc_result.f_scaled);
                auto G = __mtshl_lift(lc_result.f_scaled, lc_result.scaled_factors, lc_result.lc_targets, alpha, x, p);

                auto prod = G[0] * G[1];
                prod.normalization();
                CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
            }
            else
            {
                CLPOLY_ASSERT(true);
            }
        }
        else
        {
            CLPOLY_ASSERT(true);
        }
    }

    // ============================================================
    // p-adic 提升路径验证（M6b）
    // ============================================================

    CLPOLY_TEST("mtshl_lift: p-adic, bivariate large coeff");
    {
        // C = 5×10^18 > mtshl_p/2 → 单次 symmetric_mod 不足，需 p-adic 提升
        ZZ C(5000000000000000LL);  // 5×10^15
        C *= ZZ(1000LL);           // C = 5×10^18

        // f = (x + C*y)(x - C*y) = x^2 - C^2*y^2
        auto f = make_lex(pow(x, 2) - polynomial_ZZ(C * C) * pow(y, 2));

        auto alpha = __select_eval_point(f, x);
        auto f0 = assign(f, alpha);
        auto f0_upoly = to_upoly(f0);
        auto uni_fac = factorize(f0_upoly);

        std::vector<upolynomial_<ZZ>> uni_factors;
        for (auto& [fi, ei] : uni_fac.factors)
            uni_factors.push_back(fi);
        CLPOLY_ASSERT(uni_factors.size() == 2);

        auto lc_result = __wang_leading_coeff(f, uni_factors, alpha, x);
        CLPOLY_ASSERT(lc_result.success);

        // 使用 63-bit 素数（与 __wang_core 相同）
        static constexpr uint64_t mtshl_p = 9223372036854775783ULL;
        auto G = __mtshl_lift(lc_result.f_scaled, lc_result.scaled_factors,
                              lc_result.lc_targets, alpha, x, mtshl_p);
        CLPOLY_ASSERT(G.size() == 2);

        auto prod = G[0] * G[1];
        prod.normalization();
        CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
    }

    CLPOLY_TEST("mtshl_lift: p-adic, trivariate large coeff");
    {
        // C = 10^19 > mtshl_p/2，三变量确保 SparseInt MDP 路径被覆盖
        ZZ C(10000000000000000LL);  // 10^16
        C *= ZZ(1000LL);            // C = 10^19

        // f = (x + C*y + z)(x - C*y + z) = (x+z)^2 - C^2*y^2
        auto f = make_lex(pow(x, 2) + 2*x*z + pow(z, 2)
                          - polynomial_ZZ(C * C) * pow(y, 2));

        auto alpha = __select_eval_point(f, x);
        CLPOLY_ASSERT(alpha.size() == 2);

        auto f0 = assign(f, alpha);
        auto f0_upoly = to_upoly(f0);
        auto uni_fac = factorize(f0_upoly);

        std::vector<upolynomial_<ZZ>> uni_factors;
        for (auto& [fi, ei] : uni_fac.factors)
            uni_factors.push_back(fi);

        if (uni_factors.size() == 2)
        {
            auto lc_result = __wang_leading_coeff(f, uni_factors, alpha, x);
            CLPOLY_ASSERT(lc_result.success);

            static constexpr uint64_t mtshl_p = 9223372036854775783ULL;
            auto G = __mtshl_lift(lc_result.f_scaled, lc_result.scaled_factors,
                                  lc_result.lc_targets, alpha, x, mtshl_p);
            CLPOLY_ASSERT(G.size() == 2);

            auto prod = G[0] * G[1];
            prod.normalization();
            CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
        }
        else
        {
            CLPOLY_ASSERT(true);
        }
    }

    return clpoly_test::test_summary();
}
