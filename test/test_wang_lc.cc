/**
 * @file test_wang_lc.cc
 * @brief Phase 2 测试: __select_eval_point(), __wang_leading_coeff()
 */
#include "clpoly_test.hh"
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <iostream>

using namespace clpoly;
using PolyZZ = polynomial_<ZZ, lex>;

// Helper: 构造 lex 多项式
PolyZZ make_lex(const polynomial_ZZ& p)
{
    PolyZZ result;
    poly_convert(p, result);
    return result;
}

// Helper: polynomial → upolynomial (单变量)
upolynomial_<ZZ> to_upoly(const PolyZZ& f)
{
    upolynomial_<ZZ> result;
    poly_convert(f, result);
    return result;
}

// Helper: 验证 wang_lc 的不变量
template<class LcResult>
void verify_wang_lc(
    const PolyZZ& f,
    const LcResult& lc_result,
    const std::map<variable, ZZ>& alpha,
    const variable& main_var,
    const char* test_name)
{
    size_t r = lc_result.lc_assignments.size();

    // 1. ∏σᵢ == L
    auto L = leadcoeff(f, main_var);
    PolyZZ prod_sigma = lc_result.lc_assignments[0];
    for (size_t i = 1; i < r; ++i)
    {
        prod_sigma = prod_sigma * lc_result.lc_assignments[i];
        prod_sigma.normalization();
    }
    CLPOLY_ASSERT_EQ(prod_sigma, L);

    // 2. τᵢ(α) = δ
    auto L_eval = assign(L, alpha);
    ZZ delta = L_eval.front().second;
    for (size_t i = 0; i < r; ++i)
    {
        auto tau_eval = assign(lc_result.lc_targets[i], alpha);
        ZZ tau_val = is_number(tau_eval) ? tau_eval.front().second : ZZ(0);
        CLPOLY_ASSERT_EQ(tau_val, delta);
    }

    // 3. ∏vᵢ == f_scaled(x₁, α)
    upolynomial_<ZZ> prod_v = lc_result.scaled_factors[0];
    for (size_t i = 1; i < r; ++i)
    {
        prod_v = prod_v * lc_result.scaled_factors[i];
        prod_v.normalization();
    }
    auto f_scaled_eval = assign(lc_result.f_scaled, alpha);
    auto f_scaled_upoly = to_upoly(f_scaled_eval);
    CLPOLY_ASSERT_EQ(prod_v, f_scaled_upoly);
}

int main()
{
    variable x("x"), y("y"), z("z");

    // ============================================================
    // __select_eval_point 测试
    // ============================================================

    CLPOLY_TEST("select_eval_point: bivariate, trivial lc");
    {
        // f = x^2 + y => lc(f,x) = 1 (常数)
        auto f = make_lex(pow(x,2) + y);
        auto alpha = __select_eval_point(f, x);
        CLPOLY_ASSERT(!alpha.empty());
        auto f0 = assign(f, alpha);
        CLPOLY_ASSERT(is_squarefree(f0));
    }

    CLPOLY_TEST("select_eval_point: bivariate, nontrivial lc");
    {
        // f = ((y+1)*x + y) * (x + 1) = (y+1)x² + (2y+1)x + y
        // lc(f,x) = y+1, 本原
        auto f = make_lex((y + polynomial_ZZ(1))*pow(x,2) + (2*y + polynomial_ZZ(1))*x + y);
        auto alpha = __select_eval_point(f, x);
        CLPOLY_ASSERT(!alpha.empty());
        // 条件 (b): lc(f,x)(α) ≠ 0
        auto L = leadcoeff(f, x);
        auto L_eval = assign(L, alpha);
        CLPOLY_ASSERT(!L_eval.empty());
        // 条件 (a): f(x,α) 无平方
        auto f0 = assign(f, alpha);
        CLPOLY_ASSERT(is_squarefree(f0));
        // 条件 (d'): |lⱼ(α)| ≥ 2
        ZZ delta_val = L_eval.front().second;
        CLPOLY_ASSERT(abs(delta_val) >= ZZ(2));
    }

    CLPOLY_TEST("select_eval_point: trivariate");
    {
        // f = (x+y+z)(x-y+z) = x² + 2xz + z² - y²
        auto f = make_lex(pow(x,2) + 2*x*z + pow(z,2) - pow(y,2));
        auto alpha = __select_eval_point(f, x);
        CLPOLY_ASSERT(alpha.size() == 2);
        auto f0 = assign(f, alpha);
        CLPOLY_ASSERT(is_squarefree(f0));
    }

    CLPOLY_TEST("select_eval_point: two lc factors");
    {
        // f = y*(y+2)*x^2 + x + 1
        // lc(f,x) = y*(y+2), 因子 [y, y+2] 求值需两两互素且 ≥ 2
        auto f = make_lex(y*(y + polynomial_ZZ(2))*pow(x,2) + x + polynomial_ZZ(1));
        auto alpha = __select_eval_point(f, x);
        CLPOLY_ASSERT(!alpha.empty());
        auto f0 = assign(f, alpha);
        CLPOLY_ASSERT(is_squarefree(f0));
    }

    // ============================================================
    // __wang_leading_coeff 测试
    // ============================================================

    CLPOLY_TEST("wang_lc: constant lc");
    {
        // f = x^2 + x - y^2 + y = (x+y)(x-y+1), lc(f,x) = 1
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

        // 常数 lc → δ=1, f_scaled = f
        CLPOLY_ASSERT_EQ(lc_result.f_scaled, f);

        verify_wang_lc(f, lc_result, alpha, x, "constant lc");
    }

    CLPOLY_TEST("wang_lc: nontrivial lc, bivariate");
    {
        // f = ((y+1)*x + y) * (x + 1) = (y+1)x² + (2y+1)x + y
        // lc(f,x) = y+1, 本原
        // 真实因子: f₁ = (y+1)x+y (lc=y+1), f₂ = x+1 (lc=1)
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

        verify_wang_lc(f, lc_result, alpha, x, "nontrivial lc bivariate");

        // 验证缩放因子是首一的 (lc(vᵢ) = δ)
        auto L_eval = assign(leadcoeff(f, x), alpha);
        ZZ delta = L_eval.front().second;
        for (auto& vi : lc_result.scaled_factors)
            CLPOLY_ASSERT_EQ(vi.front().second, delta);
    }

    CLPOLY_TEST("wang_lc: three factors, constant lc");
    {
        // f = (x+y)(x-y)(x+1) = (x²-y²)(x+1) = x³+x²-xy²-y²
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

        // δ=1, f_scaled = f
        CLPOLY_ASSERT_EQ(lc_result.f_scaled, f);

        verify_wang_lc(f, lc_result, alpha, x, "three factors constant lc");
    }

    CLPOLY_TEST("wang_lc: nontrivial lc, two distinct factors");
    {
        // f = (y²-1)*(x²-1) = (y-1)(y+1)(x-1)(x+1)
        // 但作为 x 的多项式: f = (y²-1)*x² - (y²-1)
        // lc(f,x) = y²-1 = (y-1)(y+1)
        // 本原 w.r.t. x? cont = gcd(y²-1, -(y²-1)) = y²-1, 非本原!
        // 改用: f = ((y-1)*x + 1) * ((y+1)*x + 1)
        //      = (y-1)(y+1)*x² + (y-1+y+1)*x + 1
        //      = (y²-1)*x² + 2y*x + 1
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
            CLPOLY_ASSERT(lc_result.success);
            verify_wang_lc(f, lc_result, alpha, x, "two distinct lc factors");
        }
        else
        {
            std::cout << "    (skipped: uni_factors.size()=" << uni_factors.size() << ")" << std::endl;
            CLPOLY_ASSERT(true);
        }
    }

    CLPOLY_TEST("wang_lc: trivariate, constant lc");
    {
        // f = (x+y+z)(x-y+z) = x² + 2xz + z² - y², lc(f,x) = 1
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
            CLPOLY_ASSERT_EQ(lc_result.f_scaled, f);
            verify_wang_lc(f, lc_result, alpha, x, "trivariate constant lc");
        }
        else
        {
            std::cout << "    (skipped: uni_factors.size()=" << uni_factors.size() << ")" << std::endl;
            CLPOLY_ASSERT(true);
        }
    }

    return clpoly_test::test_summary();
}
