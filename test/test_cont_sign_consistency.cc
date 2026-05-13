/**
 * @file test_cont_sign_consistency.cc
 * @brief Issue #17 回归测试：cont() 始终非负（对齐 Maple/SymPy/FLINT 约定）
 *
 * 修正方案：docs/fixes/2026-05-12-cont-single-term-sign.md
 *
 * Bug：cont() 对单项多项式返回带符号值（c_0 直接返回），与多项情形
 * （gcd 折正）行为不一致。多变量 cont 在 single-group 输入下同样问题。
 *
 * 修复：
 * - polynomial_gcd.hh:1067 (univariate ZZ): init `ZZ c = 0`，让 gcd(0, c_i)
 *   走 GMP `mpz_gcd` 取绝对值，自动归正。
 * - upolynomial.hh:220 (template Tc): 同样模式。
 * - polynomial_gcd.hh:581+ (multivariate): return 前加 sign normalize。
 */
#include <clpoly/clpoly.hh>
#include <iostream>

using namespace clpoly;

static int passed = 0, total = 0;

#define CHECK(cond, label) do {                                     \
    ++total;                                                        \
    if (cond) { ++passed; std::cout << "  [PASS] " << label << std::endl; } \
    else { std::cout << "  [FAIL] " << label << std::endl; }        \
} while(0)

int main()
{
    variable x("x"), y("y");

    // ============================================================
    // §1. Univariate cont — 主回归（issue #17 复现）
    // ============================================================
    std::cout << "=== Univariate cont sign consistency ===" << std::endl;
    {
        auto p = polynomial_<ZZ>(x);
        ZZ c = cont(p);
        std::cout << "    cont(x) = " << c << " (expect 1)" << std::endl;
        CHECK(c == ZZ(1), "cont_pos_single");
    }
    {
        // bug-revealing case：单项 -x，修前 cont=-1，修后 cont=1
        auto p = polynomial_<ZZ>(-x);
        ZZ c = cont(p);
        std::cout << "    cont(-x) = " << c << " (expect 1 [fix])" << std::endl;
        CHECK(c == ZZ(1), "cont_neg_single_x");
    }
    {
        // bug-revealing case：单项 -2x，修前 -2，修后 2
        auto p = polynomial_<ZZ>(-2*x);
        ZZ c = cont(p);
        std::cout << "    cont(-2*x) = " << c << " (expect 2 [fix])" << std::endl;
        CHECK(c == ZZ(2), "cont_neg_single_2x");
    }
    {
        // 多项情形已正确，回归 sanity
        auto p = polynomial_<ZZ>(-x+1);
        ZZ c = cont(p);
        std::cout << "    cont(-x+1) = " << c << " (expect 1)" << std::endl;
        CHECK(c == ZZ(1), "cont_neg_multi_minus_x_plus_1");
    }
    {
        auto p = polynomial_<ZZ>(-2*x+2);
        ZZ c = cont(p);
        std::cout << "    cont(-2*x+2) = " << c << " (expect 2)" << std::endl;
        CHECK(c == ZZ(2), "cont_neg_multi_minus_2x_plus_2");
    }
    {
        // 常数多项式（本质是单项 mono=[]）：修前 cont=-6，修后 6
        auto p = polynomial_<ZZ>(ZZ(-6));
        ZZ c = cont(p);
        std::cout << "    cont(-6 const) = " << c << " (expect 6 [fix])" << std::endl;
        CHECK(c == ZZ(6), "cont_constant_negative");
    }
    {
        // 零多项式（empty）
        polynomial_<ZZ> p;
        ZZ c = cont(p);
        std::cout << "    cont(0 empty) = " << c << " (expect 1 [convention])" << std::endl;
        CHECK(c == ZZ(1), "cont_zero_poly");
    }

    // ============================================================
    // §2. Multivariate cont — single-group 修复（用 polynomial_<ZZ, lex>）
    // ============================================================
    std::cout << "\n=== Multivariate cont (lex order) sign consistency ===" << std::endl;
    {
        // 单项多变量 -x：mvpoly view, single-group bug
        polynomial_<ZZ, lex> p;
        poly_convert(polynomial_<ZZ>(-x), p);
        auto c = cont(p);
        std::cout << "    cont(-x mv lex) = " << c << " (expect first coef > 0 [fix])" << std::endl;
        CHECK(!c.empty() && c.front().second > ZZ(0), "mv_cont_neg_single");
    }
    {
        // -2*x*y: single x-deg 1 group
        polynomial_<ZZ, lex> p;
        poly_convert(polynomial_<ZZ>(-2*x*y), p);
        auto c = cont(p);
        std::cout << "    cont(-2xy mv lex) = " << c << " (expect first coef > 0 [fix])" << std::endl;
        CHECK(!c.empty() && c.front().second > ZZ(0), "mv_cont_neg_term_xy");
    }
    {
        // multi-group sanity (polynomial_GCD 已归正路径)
        polynomial_<ZZ, lex> p;
        poly_convert(polynomial_<ZZ>(-x+1), p);
        auto c = cont(p);
        std::cout << "    cont(-x+1 mv lex) = " << c << " (expect first coef > 0)" << std::endl;
        CHECK(!c.empty() && c.front().second > ZZ(0), "mv_cont_multi_group_sanity");
    }

    // ============================================================
    // §3. 真乘法重建：(F / cont(F)) * cont(F) = F（守恒不变量）
    //     多变量路径：cont(F) 返回 polynomial（首项 lex 之外的多项式 GCD）
    // ============================================================
    std::cout << "\n=== Invariant: (F / cont(F)) * cont(F) == F ===" << std::endl;
    {
        // 用 lex 序：cont(lex poly) 返回 polynomial_<ZZ,lex>（去首变元 GCD）
        auto test_invariant = [](const polynomial_<ZZ,lex>& F, const std::string& name) {
            if (F.empty()) return;
            auto c = cont(F);           // polynomial_<ZZ,lex>
            auto quotient = F / c;      // polynomial / polynomial
            auto reconstructed = quotient * c;
            CHECK(reconstructed == F, "F = (F/cont)*cont : " + name);
            std::cout << "    " << name << ": F=" << F
                      << ", cont=" << c << ", F/cont=" << quotient << std::endl;
        };
        auto make_lex = [](const polynomial_<ZZ>& g) {
            polynomial_<ZZ,lex> p;
            poly_convert(g, p);
            return p;
        };
        test_invariant(make_lex(polynomial_<ZZ>(x)),          "x");
        test_invariant(make_lex(polynomial_<ZZ>(-x)),         "-x");
        test_invariant(make_lex(polynomial_<ZZ>(-2*x)),       "-2x");
        test_invariant(make_lex(polynomial_<ZZ>(-x+1)),       "-x+1");
        test_invariant(make_lex(polynomial_<ZZ>(-2*x+2)),     "-2x+2");
        test_invariant(make_lex(polynomial_<ZZ>(-2*x*y+3*x)), "-2xy+3x");
    }

    // ============================================================
    // Summary
    // ============================================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "Summary: " << passed << "/" << total << " passed" << std::endl;
    return passed == total ? 0 : 1;
}
