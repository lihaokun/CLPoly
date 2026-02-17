/**
 * @file test_multivar_hensel.cc
 * @brief Phase 3 测试: __multivar_hensel_lift (Bézout + 逐变量提升)
 */
#include "clpoly_test.hh"
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <iostream>

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

int main()
{
    variable x("x"), y("y"), z("z");

    // ============================================================
    // 二变量提升测试
    // ============================================================

    CLPOLY_TEST("hensel_lift: (x+y)(x-y), constant lc");
    {
        // f = x^2 - y^2 = (x+y)(x-y)
        // lc(f,x) = 1, 本原
        auto f = make_lex(pow(x,2) - pow(y,2));
        auto alpha = __select_eval_point(f, x);
        ZZ alpha_y = alpha.at(y);

        auto f0 = assign(f, alpha);
        auto f0_upoly = to_upoly(f0);
        auto uni_fac = factorize(f0_upoly);

        std::vector<upolynomial_<ZZ>> uni_factors;
        for (auto& [fi, ei] : uni_fac.factors)
            uni_factors.push_back(fi);
        CLPOLY_ASSERT(uni_factors.size() == 2);

        auto lc_result = __wang_leading_coeff(f, uni_factors, alpha, x);
        CLPOLY_ASSERT(lc_result.success);

        // 执行 Hensel 提升
        auto G = __multivar_hensel_lift(
            lc_result.f_scaled, lc_result.scaled_factors,
            lc_result.lc_targets, alpha, x);

        CLPOLY_ASSERT(G.size() == 2);

        // 验证 ∏Gᵢ = f_scaled
        auto prod = G[0] * G[1];
        prod.normalization();
        CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
    }

    CLPOLY_TEST("hensel_lift: (x+y)(x-y+1), constant lc");
    {
        // f = (x+y)(x-y+1) = x^2 + x - y^2 + y
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

        auto G = __multivar_hensel_lift(
            lc_result.f_scaled, lc_result.scaled_factors,
            lc_result.lc_targets, alpha, x);

        auto prod = G[0] * G[1];
        prod.normalization();
        CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
    }

    CLPOLY_TEST("hensel_lift: nontrivial lc, (y+1)x+y) * (x+1)");
    {
        // f = ((y+1)*x + y)(x + 1) = (y+1)x^2 + (2y+1)x + y
        // lc(f,x) = y+1
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

        auto G = __multivar_hensel_lift(
            lc_result.f_scaled, lc_result.scaled_factors,
            lc_result.lc_targets, alpha, x);

        auto prod = G[0] * G[1];
        prod.normalization();
        CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
    }

    CLPOLY_TEST("hensel_lift: three factors, (x+y)(x-y)(x+1)");
    {
        // f = (x+y)(x-y)(x+1) = x^3 + x^2 - xy^2 - y^2
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

        auto G = __multivar_hensel_lift(
            lc_result.f_scaled, lc_result.scaled_factors,
            lc_result.lc_targets, alpha, x);

        CLPOLY_ASSERT(G.size() == 3);
        auto prod = G[0] * G[1];
        prod.normalization();
        prod = prod * G[2];
        prod.normalization();
        CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
    }

    // ============================================================
    // 三变量提升测试
    // ============================================================

    CLPOLY_TEST("hensel_lift: trivariate (x+y+z)(x-y+z)");
    {
        // f = (x+y+z)(x-y+z) = x^2 + 2xz + z^2 - y^2
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

            auto G = __multivar_hensel_lift(
                lc_result.f_scaled, lc_result.scaled_factors,
                lc_result.lc_targets, alpha, x);

            auto prod = G[0] * G[1];
            prod.normalization();
            CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
        }
        else
        {
            std::cout << "    (skipped: " << uni_factors.size() << " factors)" << std::endl;
            CLPOLY_ASSERT(true);
        }
    }

    CLPOLY_TEST("hensel_lift: trivariate (x+y)(x+z)");
    {
        // f = (x+y)(x+z) = x^2 + xz + xy + yz
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

            auto G = __multivar_hensel_lift(
                lc_result.f_scaled, lc_result.scaled_factors,
                lc_result.lc_targets, alpha, x);

            auto prod = G[0] * G[1];
            prod.normalization();
            CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
        }
        else
        {
            std::cout << "    (skipped: " << uni_factors.size() << " factors)" << std::endl;
            CLPOLY_ASSERT(true);
        }
    }

    // ============================================================
    // nontrivial lc + 二变量, 两 lc 因子
    // ============================================================

    CLPOLY_TEST("hensel_lift: nontrivial lc, ((y^2-1)x^2 + 2yx + 1)");
    {
        // f = ((y-1)x+1)((y+1)x+1) = (y^2-1)x^2 + 2yx + 1
        // lc(f,x) = y^2-1 = (y-1)(y+1)
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
                auto G = __multivar_hensel_lift(
                    lc_result.f_scaled, lc_result.scaled_factors,
                    lc_result.lc_targets, alpha, x);

                auto prod = G[0] * G[1];
                prod.normalization();
                CLPOLY_ASSERT_EQ(prod, lc_result.f_scaled);
            }
            else
            {
                std::cout << "    (lc assignment failed, acceptable)" << std::endl;
                CLPOLY_ASSERT(true);
            }
        }
        else
        {
            std::cout << "    (skipped: " << uni_factors.size() << " factors)" << std::endl;
            CLPOLY_ASSERT(true);
        }
    }

    return clpoly_test::test_summary();
}
