#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

// 辅助：构建 Zp 多变量多项式（lex 序）
// terms: {(e1, e2, coeff), ...} for bivariate
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

// 辅助：构建三变量 Zp 多项式（lex 序）
// terms: {(e1, e2, e3, coeff)}
template<class Comp>
static polynomial_<Zp, Comp> make_trivar_zp(
    const Comp* comp_ptr,
    const variable& x1, const variable& x2, const variable& x3,
    std::initializer_list<std::tuple<int,int,int,int>> terms, uint32_t p)
{
    polynomial_<Zp, Comp> poly(comp_ptr);
    for (auto& [e1, e2, e3, c] : terms)
    {
        basic_monomial<Comp> m(comp_ptr);
        if (e1 > 0) m.push_back({x1, e1});
        if (e2 > 0) m.push_back({x2, e2});
        if (e3 > 0) m.push_back({x3, e3});
        m.normalization();
        Zp val(c, p);
        if (val.number() != 0) poly.push_back({m, val});
    }
    poly.normalization();
    return poly;
}

// 辅助：验证 Σ result[i] * ∏_{l≠i} F[l] = c
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

int main()
{
    // ================================================================
    // §M3  __mtshl_sparse_int（稀疏插值 MDP）
    // 测试策略：构造已知 σi + forms，计算 c = Σ σi * bi，
    //           用 sparse_int 恢复 σi，验证 Σ result*bi = c
    //           + 与 multi_bdp 交叉验证
    // ================================================================

    CLPOLY_TEST("sparse_int_r2_bivar");
    {
        // p=101, F1=x1+x2, F2=x1-x2+1
        // 已知解: σ1=x2, σ2=3
        // forms[0] = {x2}, forms[1] = {1}（常数项）
        uint32_t p = 101;
        variable x1("x1"), x2("x2");
        lex comp_lex;

        std::vector<polynomial_<Zp, lex>> F = {
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,1,1}}, p),              // x1 + x2
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,1,100}, {0,0,1}}, p),   // x1 - x2 + 1
        };
        polynomial_<Zp, lex> sig0 = make_bivar_zp(&comp_lex, x1, x2, {{0,1,1}}, p);    // x2
        polynomial_<Zp, lex> sig1 = make_bivar_zp(&comp_lex, x1, x2, {{0,0,3}}, p);    // 3

        // b0 = F[1], b1 = F[0]
        auto t0 = sig0 * F[1]; t0.normalization();
        auto t1 = sig1 * F[0]; t1.normalization();
        auto c = t0 + t1; c.normalization();

        // forms: σ0 的支撑 = {x2}，σ1 的支撑 = {1}
        basic_monomial<lex> m_x2(&comp_lex);
        m_x2.push_back({x2, 1}); m_x2.normalization();
        basic_monomial<lex> m_const(&comp_lex);
        m_const.normalization();

        std::vector<std::vector<basic_monomial<lex>>> forms = {
            {m_x2},     // forms[0]
            {m_const},  // forms[1]
        };

        std::vector<polynomial_<Zp, lex>> result;
        bool ok = __mtshl_sparse_int(F, c, forms, x1, {x2}, p, result);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)result.size(), 2);
        CLPOLY_ASSERT_TRUE(verify_multivar_mdp(F, c, result));
    }

    CLPOLY_TEST("sparse_int_r2_bivar_multi_term");
    {
        // σ0 = 2*x2 + 5, σ1 = x2 + 1（多项解）
        uint32_t p = 101;
        variable x1("x1"), x2("x2");
        lex comp_lex;

        std::vector<polynomial_<Zp, lex>> F = {
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,1}}, p),   // x1 + 1
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,2}}, p),   // x1 + 2
        };
        polynomial_<Zp, lex> sig0 = make_bivar_zp(&comp_lex, x1, x2,
            {{0,1,2}, {0,0,5}}, p);  // 2*x2 + 5
        polynomial_<Zp, lex> sig1 = make_bivar_zp(&comp_lex, x1, x2,
            {{0,1,1}, {0,0,1}}, p);  // x2 + 1

        auto t0 = sig0 * F[1]; t0.normalization();
        auto t1 = sig1 * F[0]; t1.normalization();
        auto c = t0 + t1; c.normalization();

        // forms: σ0 支撑 = {x2, 1}，σ1 支撑 = {x2, 1}
        basic_monomial<lex> m_x2(&comp_lex);
        m_x2.push_back({x2, 1}); m_x2.normalization();
        basic_monomial<lex> m_const(&comp_lex);
        m_const.normalization();

        std::vector<std::vector<basic_monomial<lex>>> forms = {
            {m_x2, m_const},
            {m_x2, m_const},
        };

        std::vector<polynomial_<Zp, lex>> result;
        bool ok = __mtshl_sparse_int(F, c, forms, x1, {x2}, p, result);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)result.size(), 2);
        CLPOLY_ASSERT_TRUE(verify_multivar_mdp(F, c, result));
    }

    CLPOLY_TEST("sparse_int_r3_bivar");
    {
        // r=3, p=97
        // F1=x1+1, F2=x1+x2, F3=x1+2*x2+3
        // σ0=2, σ1=x2+1, σ2=3
        uint32_t p = 97;
        variable x1("x1"), x2("x2");
        lex comp_lex;

        std::vector<polynomial_<Zp, lex>> F = {
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,1}}, p),
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,1,1}}, p),
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,1,2}, {0,0,3}}, p),
        };

        polynomial_<Zp, lex> sig0 = make_bivar_zp(&comp_lex, x1, x2, {{0,0,2}}, p);
        polynomial_<Zp, lex> sig1 = make_bivar_zp(&comp_lex, x1, x2, {{0,1,1}, {0,0,1}}, p);
        polynomial_<Zp, lex> sig2 = make_bivar_zp(&comp_lex, x1, x2, {{0,0,3}}, p);

        auto b0 = F[1] * F[2]; b0.normalization();
        auto b1 = F[0] * F[2]; b1.normalization();
        auto b2 = F[0] * F[1]; b2.normalization();
        auto c = sig0 * b0; c.normalization();
        auto t1 = sig1 * b1; t1.normalization();
        auto t2 = sig2 * b2; t2.normalization();
        c = c + t1; c.normalization();
        c = c + t2; c.normalization();

        // forms
        basic_monomial<lex> m_x2(&comp_lex);
        m_x2.push_back({x2, 1}); m_x2.normalization();
        basic_monomial<lex> m_const(&comp_lex);
        m_const.normalization();

        std::vector<std::vector<basic_monomial<lex>>> forms = {
            {m_const},
            {m_x2, m_const},
            {m_const},
        };

        std::vector<polynomial_<Zp, lex>> result;
        bool ok = __mtshl_sparse_int(F, c, forms, x1, {x2}, p, result);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)result.size(), 3);
        CLPOLY_ASSERT_TRUE(verify_multivar_mdp(F, c, result));
    }

    CLPOLY_TEST("sparse_int_zero_c");
    {
        uint32_t p = 101;
        variable x1("x1"), x2("x2");
        lex comp_lex;

        std::vector<polynomial_<Zp, lex>> F = {
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,1}}, p),
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,2}}, p),
        };
        polynomial_<Zp, lex> c(&comp_lex);

        basic_monomial<lex> m_const(&comp_lex);
        m_const.normalization();
        std::vector<std::vector<basic_monomial<lex>>> forms = {
            {m_const}, {m_const},
        };

        std::vector<polynomial_<Zp, lex>> result;
        bool ok = __mtshl_sparse_int(F, c, forms, x1, {x2}, p, result);
        CLPOLY_ASSERT_TRUE(ok);
        for (auto& ri : result)
            CLPOLY_ASSERT_TRUE(ri.empty());
    }

    CLPOLY_TEST("sparse_int_empty_forms");
    {
        // forms 全空 → 直接返回全零
        uint32_t p = 101;
        variable x1("x1"), x2("x2");
        lex comp_lex;

        std::vector<polynomial_<Zp, lex>> F = {
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,1}}, p),
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,0,2}}, p),
        };
        polynomial_<Zp, lex> c(&comp_lex);  // c=0 也是

        std::vector<std::vector<basic_monomial<lex>>> forms = {{}, {}};

        std::vector<polynomial_<Zp, lex>> result;
        bool ok = __mtshl_sparse_int(F, c, forms, x1, {x2}, p, result);
        CLPOLY_ASSERT_TRUE(ok);
        for (auto& ri : result)
            CLPOLY_ASSERT_TRUE(ri.empty());
    }

    CLPOLY_TEST("sparse_int_trivar_r2");
    {
        // 三变量: F[i] ∈ Zp[x1,x2,x3]
        // F1 = x1+x2+x3, F2 = x1-x2+1
        // σ0 = x2+2*x3, σ1 = x3+5
        uint32_t p = 101;
        variable x1("x1"), x2("x2"), x3("x3");
        lex comp_lex;

        auto F1 = make_trivar_zp(&comp_lex, x1, x2, x3,
            {{1,0,0,1}, {0,1,0,1}, {0,0,1,1}}, p);        // x1+x2+x3
        auto F2 = make_trivar_zp(&comp_lex, x1, x2, x3,
            {{1,0,0,1}, {0,1,0,100}, {0,0,0,1}}, p);      // x1-x2+1

        auto sig0 = make_trivar_zp(&comp_lex, x1, x2, x3,
            {{0,1,0,1}, {0,0,1,2}}, p);                    // x2+2*x3
        auto sig1 = make_trivar_zp(&comp_lex, x1, x2, x3,
            {{0,0,1,1}, {0,0,0,5}}, p);                    // x3+5

        std::vector<polynomial_<Zp, lex>> F = {F1, F2};
        auto t0 = sig0 * F2; t0.normalization();
        auto t1 = sig1 * F1; t1.normalization();
        auto c = t0 + t1; c.normalization();

        // forms
        basic_monomial<lex> m_x2(&comp_lex);
        m_x2.push_back({x2, 1}); m_x2.normalization();
        basic_monomial<lex> m_x3(&comp_lex);
        m_x3.push_back({x3, 1}); m_x3.normalization();
        basic_monomial<lex> m_const(&comp_lex);
        m_const.normalization();

        std::vector<std::vector<basic_monomial<lex>>> forms = {
            {m_x2, m_x3},       // σ0 支撑: x2, x3
            {m_x3, m_const},    // σ1 支撑: x3, 1
        };

        std::vector<polynomial_<Zp, lex>> result;
        bool ok = __mtshl_sparse_int(F, c, forms, x1, {x2, x3}, p, result);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)result.size(), 2);
        CLPOLY_ASSERT_TRUE(verify_multivar_mdp(F, c, result));
    }

    // 交叉验证：sparse_int 与 multi_bdp 在双变量情形结果一致
    CLPOLY_TEST("sparse_int_vs_multi_bdp_crosscheck");
    {
        uint32_t p = 101;
        variable x1("x1"), x2("x2");
        lex comp_lex;

        std::vector<polynomial_<Zp, lex>> F = {
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,1,1}}, p),              // x1+x2
            make_bivar_zp(&comp_lex, x1, x2, {{1,0,1}, {0,1,100}, {0,0,1}}, p),   // x1-x2+1
        };
        polynomial_<Zp, lex> sig0 = make_bivar_zp(&comp_lex, x1, x2, {{0,1,1}}, p);
        polynomial_<Zp, lex> sig1 = make_bivar_zp(&comp_lex, x1, x2, {{0,0,3}}, p);
        auto t0 = sig0 * F[1]; t0.normalization();
        auto t1 = sig1 * F[0]; t1.normalization();
        auto c = t0 + t1; c.normalization();

        // sparse_int
        basic_monomial<lex> m_x2(&comp_lex);
        m_x2.push_back({x2, 1}); m_x2.normalization();
        basic_monomial<lex> m_const(&comp_lex);
        m_const.normalization();
        std::vector<std::vector<basic_monomial<lex>>> forms = {
            {m_x2}, {m_const},
        };
        std::vector<polynomial_<Zp, lex>> result_sparse;
        bool ok1 = __mtshl_sparse_int(F, c, forms, x1, {x2}, p, result_sparse);
        CLPOLY_ASSERT_TRUE(ok1);

        // multi_bdp
        Zp alpha2(2, p);
        std::vector<polynomial_<Zp, lex>> result_dense;
        bool ok2 = __mtshl_multi_bdp(F, c, x1, x2, alpha2, result_dense);
        CLPOLY_ASSERT_TRUE(ok2);

        // 两者应给出等价解（Σ result*bi = c）
        CLPOLY_ASSERT_TRUE(verify_multivar_mdp(F, c, result_sparse));
        CLPOLY_ASSERT_TRUE(verify_multivar_mdp(F, c, result_dense));

        // 由于 MDP 解唯一（deg(σi,x1) < deg(F[i],x1)），结果应完全相同
        CLPOLY_ASSERT_EQ((int)result_sparse.size(), (int)result_dense.size());
        for (int i = 0; i < (int)result_sparse.size(); i++)
            CLPOLY_ASSERT_TRUE(result_sparse[i] == result_dense[i]);
    }

    return clpoly_test::test_summary();
}
