#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

// 辅助：构建 Zp 单变量多项式 from {(deg, coeff)} pairs
static upolynomial_<Zp> make_upzp(
    std::initializer_list<std::pair<int64_t, int>> terms, uint32_t p)
{
    upolynomial_<Zp> poly;
    for (auto& [d, c] : terms)
    {
        Zp val(c, p);
        if (val.number() != 0)
            poly.push_back({umonomial(d), val});
    }
    return poly;
}

// 辅助：验证 Σ sigma[i] * ∏_{l≠i} F[l] = c（单变量）
static bool verify_univar_mdp(
    const std::vector<upolynomial_<Zp>>& F,
    const upolynomial_<Zp>& c,
    const std::vector<upolynomial_<Zp>>& sigma)
{
    int r = (int)F.size();
    upolynomial_<Zp> lhs;
    for (int i = 0; i < r; i++)
    {
        // bi = ∏_{l≠i} F[l]
        upolynomial_<Zp> bi;
        bool first = true;
        for (int l = 0; l < r; l++)
        {
            if (l == i) continue;
            if (first) { bi = F[l]; first = false; }
            else { bi = bi * F[l]; bi.normalization(); }
        }
        auto term = sigma[i] * bi;
        term.normalization();
        if (lhs.empty())
            lhs = term;
        else
        {
            lhs = lhs + term;
            lhs.normalization();
        }
    }
    return lhs == c;
}

// 辅助：验证 Σ result[i] * ∏_{l≠i} F[l] = c（多变量 Zp）
template<class Comp>
static bool verify_multivar_mdp(
    const std::vector<polynomial_<Zp, Comp>>& F,
    const polynomial_<Zp, Comp>& c,
    const std::vector<polynomial_<Zp, Comp>>& result)
{
    int r = (int)F.size();
    auto comp_ptr = c.comp_ptr();
    polynomial_<Zp, Comp> lhs(comp_ptr);
    for (int i = 0; i < r; i++)
    {
        polynomial_<Zp, Comp> bi(comp_ptr);
        bool first = true;
        for (int l = 0; l < r; l++)
        {
            if (l == i) continue;
            if (first) { bi = F[l]; first = false; }
            else { bi = bi * F[l]; bi.normalization(); }
        }
        auto term = result[i] * bi;
        term.normalization();
        lhs = lhs + term;
        lhs.normalization();
    }
    return lhs == c;
}

// 辅助：构建 Zp 多变量多项式（lex 序）
// terms: {(e1, e2, coeff), ...}
template<class Comp>
static polynomial_<Zp, Comp> make_bivar_zp(
    const Comp* comp_ptr,
    const variable& x1, const variable& x2,
    std::initializer_list<std::tuple<int,int,int>> terms, uint32_t p)
{
    polynomial_<Zp, Comp> poly(comp_ptr);
    for (auto& [e1, e2, c] : terms)
    {
        basic_monomial<Comp> m(comp_ptr);
        if (e1 > 0) m.push_back({x1, e1});
        if (e2 > 0) m.push_back({x2, e2});
        m.normalization();
        Zp val(c, p);
        if (val.number() != 0) poly.push_back({m, val});
    }
    poly.normalization();
    return poly;
}

int main()
{
    // ================================================================
    // §M2-a  __mtshl_zp_univar_mdp（单变量 Bézout 链）
    // ================================================================

    CLPOLY_TEST("univar_mdp_r2");
    {
        // p=101, F = [x+2, x+3], c = 1（deg(c)=0 < deg(∏F)=2）
        uint32_t p = 101;
        std::vector<upolynomial_<Zp>> F = {
            make_upzp({{1, 1}, {0, 2}}, p),  // x + 2
            make_upzp({{1, 1}, {0, 3}}, p),  // x + 3
        };
        upolynomial_<Zp> c = make_upzp({{0, 1}}, p);
        std::vector<upolynomial_<Zp>> sigma;
        bool ok = __mtshl_zp_univar_mdp(F, c, sigma);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)sigma.size(), 2);
        CLPOLY_ASSERT_TRUE(verify_univar_mdp(F, c, sigma));
    }

    CLPOLY_TEST("univar_mdp_r2_linear_c");
    {
        // p=101, F = [x+1, x+100], c = 5x + 3（deg(c)=1 < deg(∏F)=2）
        uint32_t p = 101;
        std::vector<upolynomial_<Zp>> F = {
            make_upzp({{1, 1}, {0, 1}}, p),    // x + 1
            make_upzp({{1, 1}, {0, 100}}, p),  // x - 1
        };
        upolynomial_<Zp> c = make_upzp({{1, 5}, {0, 3}}, p);
        std::vector<upolynomial_<Zp>> sigma;
        bool ok = __mtshl_zp_univar_mdp(F, c, sigma);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)sigma.size(), 2);
        CLPOLY_ASSERT_TRUE(verify_univar_mdp(F, c, sigma));
        // deg(sigma[i]) < deg(F[i]) = 1 → sigma[i] are constants
        CLPOLY_ASSERT_TRUE(get_deg(sigma[0]) <= 0);
        CLPOLY_ASSERT_TRUE(get_deg(sigma[1]) <= 0);
    }

    CLPOLY_TEST("univar_mdp_r3");
    {
        // p=101, F = [x+1, x+2, x+3], c = x^2 + 1（deg=2 < deg(∏F)=3）
        uint32_t p = 101;
        std::vector<upolynomial_<Zp>> F = {
            make_upzp({{1, 1}, {0, 1}}, p),
            make_upzp({{1, 1}, {0, 2}}, p),
            make_upzp({{1, 1}, {0, 3}}, p),
        };
        upolynomial_<Zp> c = make_upzp({{2, 1}, {0, 1}}, p);
        std::vector<upolynomial_<Zp>> sigma;
        bool ok = __mtshl_zp_univar_mdp(F, c, sigma);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)sigma.size(), 3);
        CLPOLY_ASSERT_TRUE(verify_univar_mdp(F, c, sigma));
    }

    CLPOLY_TEST("univar_mdp_r4");
    {
        // p=97, F = [x+1, x+2, x+3, x+4], c = 2x^3 + x + 5（deg=3 < deg(∏F)=4）
        uint32_t p = 97;
        std::vector<upolynomial_<Zp>> F = {
            make_upzp({{1, 1}, {0, 1}}, p),
            make_upzp({{1, 1}, {0, 2}}, p),
            make_upzp({{1, 1}, {0, 3}}, p),
            make_upzp({{1, 1}, {0, 4}}, p),
        };
        upolynomial_<Zp> c = make_upzp({{3, 2}, {1, 1}, {0, 5}}, p);
        std::vector<upolynomial_<Zp>> sigma;
        bool ok = __mtshl_zp_univar_mdp(F, c, sigma);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)sigma.size(), 4);
        CLPOLY_ASSERT_TRUE(verify_univar_mdp(F, c, sigma));
    }

    CLPOLY_TEST("univar_mdp_zero_c");
    {
        uint32_t p = 101;
        std::vector<upolynomial_<Zp>> F = {
            make_upzp({{1, 1}, {0, 1}}, p),
            make_upzp({{1, 1}, {0, 2}}, p),
        };
        upolynomial_<Zp> c;
        std::vector<upolynomial_<Zp>> sigma;
        bool ok = __mtshl_zp_univar_mdp(F, c, sigma);
        CLPOLY_ASSERT_TRUE(ok);
        for (auto& si : sigma)
            CLPOLY_ASSERT_TRUE(si.empty());
    }

    // ================================================================
    // §M2-b  __mtshl_multi_bdp（双变量 Taylor 提升）
    // 测试策略：先构造已知 σi，再计算 c = Σ σi * bi，
    //           然后用 multi_bdp 恢复 σi，验证 Σ result*bi = c
    // ================================================================

    CLPOLY_TEST("multi_bdp_r2_bivar");
    {
        // p=101, F1 = x1+x2, F2 = x1-x2+1
        // 已知解: σ1 = x2, σ2 = 3（deg < deg(F[i],x1)=1, 均为 x2 多项式）
        // c = σ1*b1 + σ2*b2 = x2*(x1-x2+1) + 3*(x1+x2)
        //   = x1*x2-x2^2+x2 + 3*x1+3*x2 = x1*x2+3*x1-x2^2+4*x2
        uint32_t p = 101;
        variable x1("x1"), x2("x2");
        lex comp_lex;

        std::vector<polynomial_<Zp, lex>> F = {
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,1,1}}, p),              // x1 + x2
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,1,100}, {0,0,1}}, p),   // x1 - x2 + 1
        };
        // c = x1*x2 + 3*x1 - x2^2 + 4*x2
        polynomial_<Zp, lex> c = make_bivar_zp(&comp_lex, x1, x2,
            {{1,1,1}, {1,0,3}, {0,2,100}, {0,1,4}}, p);  // -x2^2 = 100*x2^2 (mod 101)
        Zp alpha2(2, p);

        std::vector<polynomial_<Zp, lex>> result;
        bool ok = __mtshl_multi_bdp(F, c, x1, x2, alpha2, result);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)result.size(), 2);
        CLPOLY_ASSERT_TRUE(verify_multivar_mdp(F, c, result));
    }

    CLPOLY_TEST("multi_bdp_r2_constant_c");
    {
        // F1 = x1+1, F2 = x1+2（无 x2 项，纯单变量情形）
        // c = 5x1 + 3（deg < deg(∏F)=2）
        // multi_bdp 退化为纯单变量 MDP（Taylor 提升无事可做）
        uint32_t p = 101;
        variable x1("x1"), x2("x2");
        lex comp_lex;

        std::vector<polynomial_<Zp, lex>> F = {
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,1}}, p),   // x1 + 1
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,2}}, p),   // x1 + 2
        };
        polynomial_<Zp, lex> c = make_bivar_zp(&comp_lex, x1, x2,
            {{1,0,5}, {0,0,3}}, p);
        Zp alpha2(0, p);

        std::vector<polynomial_<Zp, lex>> result;
        bool ok = __mtshl_multi_bdp(F, c, x1, x2, alpha2, result);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)result.size(), 2);
        CLPOLY_ASSERT_TRUE(verify_multivar_mdp(F, c, result));
    }

    CLPOLY_TEST("multi_bdp_r3_bivar");
    {
        // r=3, p=97
        // F1 = x1+1, F2 = x1+x2, F3 = x1+2*x2+3
        // 已知解: σ1 = 2, σ2 = x2+1, σ3 = 3
        // b1 = F2*F3, b2 = F1*F3, b3 = F1*F2
        // 构造 c = Σ σi*bi
        uint32_t p = 97;
        variable x1("x1"), x2("x2");
        lex comp_lex;

        std::vector<polynomial_<Zp, lex>> F = {
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,1}}, p),             // x1 + 1
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,1,1}}, p),             // x1 + x2
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,1,2}, {0,0,3}}, p),    // x1 + 2*x2 + 3
        };

        // 已知解
        polynomial_<Zp, lex> sig0 = make_bivar_zp(&comp_lex, x1, x2, {{0,0,2}}, p);           // 2
        polynomial_<Zp, lex> sig1 = make_bivar_zp(&comp_lex, x1, x2, {{0,1,1}, {0,0,1}}, p);  // x2 + 1
        polynomial_<Zp, lex> sig2 = make_bivar_zp(&comp_lex, x1, x2, {{0,0,3}}, p);           // 3

        // 计算 bi
        auto b0 = F[1] * F[2]; b0.normalization();
        auto b1 = F[0] * F[2]; b1.normalization();
        auto b2 = F[0] * F[1]; b2.normalization();

        // c = σ0*b0 + σ1*b1 + σ2*b2
        auto t0 = sig0 * b0; t0.normalization();
        auto t1 = sig1 * b1; t1.normalization();
        auto t2 = sig2 * b2; t2.normalization();
        auto c = t0 + t1; c.normalization();
        c = c + t2; c.normalization();

        Zp alpha2(3, p);

        std::vector<polynomial_<Zp, lex>> result;
        bool ok = __mtshl_multi_bdp(F, c, x1, x2, alpha2, result);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)result.size(), 3);
        CLPOLY_ASSERT_TRUE(verify_multivar_mdp(F, c, result));
    }

    CLPOLY_TEST("multi_bdp_zero_c");
    {
        uint32_t p = 101;
        variable x1("x1"), x2("x2");
        lex comp_lex;

        std::vector<polynomial_<Zp, lex>> F = {
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,1}}, p),
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,2}}, p),
        };
        polynomial_<Zp, lex> c(&comp_lex);
        Zp alpha2(0, p);

        std::vector<polynomial_<Zp, lex>> result;
        bool ok = __mtshl_multi_bdp(F, c, x1, x2, alpha2, result);
        CLPOLY_ASSERT_TRUE(ok);
        for (auto& ri : result)
            CLPOLY_ASSERT_TRUE(ri.empty());
    }

    return clpoly_test::test_summary();
}
