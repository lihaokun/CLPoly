#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

// ================================================================
// 辅助：从 UPZp 中按次数提取系数
// ================================================================
static uint64_t get_coeff(const upolynomial_<Zp>& p, int64_t deg)
{
    for (const auto& term : p)
        if (term.first.deg() == deg)
            return term.second.number();
    return 0;
}

int main()
{
    // ================================================================
    // §M1-A  __si_vandermonde_solve
    // ================================================================

    // ---- s=1 --------------------------------------------------------
    CLPOLY_TEST("vandermonde_s1");
    {
        // v[1] = c1 * θ1^1；设 p=101, θ1=5, c1=3 → v[1]=15
        uint64_t p = 101;
        std::vector<Zp> values = { Zp(15, p) };
        std::vector<Zp> thetas = { Zp(5,  p) };
        std::vector<Zp> coeffs;
        bool ok = __si_vandermonde_solve(values, thetas, coeffs);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)coeffs.size(), 1);
        // c1 = v[1]/θ1 = 15/5 = 3
        CLPOLY_ASSERT_EQ((int)coeffs[0].number(), 3);
    }

    // ---- s=2 --------------------------------------------------------
    CLPOLY_TEST("vandermonde_s2");
    {
        // p=101, θ=(5,7), c=(3,4)
        // v[1] = 3*5 + 4*7 = 15+28 = 43
        // v[2] = 3*25 + 4*49 = 75+196 = 271 ≡ 69 (mod 101)
        uint64_t p = 101;
        std::vector<Zp> values = { Zp(43, p), Zp(69, p) };
        std::vector<Zp> thetas = { Zp(5,  p), Zp(7,  p) };
        std::vector<Zp> coeffs;
        bool ok = __si_vandermonde_solve(values, thetas, coeffs);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)coeffs.size(), 2);
        CLPOLY_ASSERT_EQ((int)coeffs[0].number(), 3);
        CLPOLY_ASSERT_EQ((int)coeffs[1].number(), 4);
    }

    // ---- s=3 --------------------------------------------------------
    CLPOLY_TEST("vandermonde_s3");
    {
        // p=101, θ=(2,3,5), c=(1,2,3)
        // v[l] = 1*2^l + 2*3^l + 3*5^l
        // v[1] = 2 + 6 + 15 = 23
        // v[2] = 4 + 18 + 75 = 97
        // v[3] = 8 + 54 + 375 = 437 ≡ 437-4*101=33 (mod 101)
        uint64_t p = 101;
        std::vector<Zp> values = { Zp(23, p), Zp(97, p), Zp(33, p) };
        std::vector<Zp> thetas = { Zp(2,  p), Zp(3,  p), Zp(5,  p) };
        std::vector<Zp> coeffs;
        bool ok = __si_vandermonde_solve(values, thetas, coeffs);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)coeffs.size(), 3);
        CLPOLY_ASSERT_EQ((int)coeffs[0].number(), 1);
        CLPOLY_ASSERT_EQ((int)coeffs[1].number(), 2);
        CLPOLY_ASSERT_EQ((int)coeffs[2].number(), 3);
    }

    // ---- 碰撞 → 奇异 -----------------------------------------------
    CLPOLY_TEST("vandermonde_collision");
    {
        // θ = (5, 5)：矩阵行 [1,1] 和 [5,5] 线性相关 → 奇异
        uint64_t p = 101;
        std::vector<Zp> values = { Zp(10, p), Zp(50, p) };
        std::vector<Zp> thetas = { Zp(5,  p), Zp(5,  p) };
        std::vector<Zp> coeffs;
        bool ok = __si_vandermonde_solve(values, thetas, coeffs);
        CLPOLY_ASSERT_FALSE(ok);
    }

    // ---- 逆向验证：随机系数 round-trip --------------------------------
    CLPOLY_TEST("vandermonde_roundtrip");
    {
        // p=97, s=4, c=(2,5,7,3), θ=(3,4,6,8)
        // 先正向计算 v[l]，再用 vandermonde_solve 恢复 c
        uint64_t p = 97;
        std::vector<Zp> thetas = { Zp(3,p), Zp(4,p), Zp(6,p), Zp(8,p) };
        std::vector<Zp> c_orig = { Zp(2,p), Zp(5,p), Zp(7,p), Zp(3,p) };
        int s = 4;

        // 正向: v[l] = Σ c_t * θ_t^l, l=1..s
        std::vector<Zp> values(s, Zp(0,p));
        for (int l = 0; l < s; l++)
            for (int t = 0; t < s; t++)
                values[l] = values[l] + c_orig[t] * pow(thetas[t], (int64_t)(l + 1));

        std::vector<Zp> coeffs;
        bool ok = __si_vandermonde_solve(values, thetas, coeffs);
        CLPOLY_ASSERT_TRUE(ok);
        CLPOLY_ASSERT_EQ((int)coeffs.size(), s);
        for (int t = 0; t < s; t++)
            CLPOLY_ASSERT_EQ(coeffs[t].number(), c_orig[t].number());
    }

    // ================================================================
    // §M1-B  __si_theta_array_eval
    // ================================================================

    // ---- 单变量退化情形（aux_vars 为空）-----------------------------
    CLPOLY_TEST("theta_array_eval_univar");
    {
        // f = 3*x^2 + 2*x + 1 ∈ Zp[x]（无 aux 变量）
        // θ-array 退化：θ_m=1 (乘积为空), running_m=1^l=1, images[l-1]=f(x)
        uint64_t p = 101;
        variable x("x");
        univariate_priority_order comp(x);
        polynomial_<Zp, univariate_priority_order> f(&comp);
        {
            basic_monomial<univariate_priority_order> m2(&comp);
            m2.push_back({x, 2}); m2.normalization();
            f.push_back({m2, Zp(3, p)});
            basic_monomial<univariate_priority_order> m1(&comp);
            m1.push_back({x, 1}); m1.normalization();
            f.push_back({m1, Zp(2, p)});
            basic_monomial<univariate_priority_order> m0(&comp);
            m0.normalization();
            f.push_back({m0, Zp(1, p)});
            f.normalization();
        }
        std::vector<variable> aux_vars;
        std::vector<Zp>       sparse_betas;
        std::vector<upolynomial_<Zp>> images;
        __si_theta_array_eval(f, x, aux_vars, sparse_betas, 3, images);

        CLPOLY_ASSERT_EQ((int)images.size(), 3);
        // 每个 image 应为 3*x^2 + 2*x + 1（running_m^l=1^l=1）
        for (int l = 0; l < 3; l++)
        {
            CLPOLY_ASSERT_EQ((int)get_coeff(images[l], 2), 3);
            CLPOLY_ASSERT_EQ((int)get_coeff(images[l], 1), 2);
            CLPOLY_ASSERT_EQ((int)get_coeff(images[l], 0), 1);
        }
    }

    // ---- 双变量（x1 主变量，x2 辅助变量）---------------------------
    CLPOLY_TEST("theta_array_eval_bivar");
    {
        // f = 3*x1^2*x2 + 2*x1 + 1*x2^2 ∈ Zp[x1, x2]
        // aux_vars=[x2], β2=3, s=3, p=101
        // θ_m for each term:
        //   3*x1^2*x2: θ = β2^1 = 3
        //   2*x1:      θ = β2^0 = 1
        //   1*x2^2:    θ = β2^2 = 9
        // l=1: 3*3*x1^2 + 2*1*x1 + 1*9 = 9*x1^2 + 2*x1 + 9
        // l=2: 3*9*x1^2 + 2*1*x1 + 1*81 = 27*x1^2 + 2*x1 + 81
        // l=3: 3*27*x1^2 + 2*1*x1 + 1*729 = 81*x1^2 + 2*x1 + 22
        //       (729 mod 101 = 729 - 7*101 = 22)
        uint64_t p = 101;
        variable x1("x1"), x2("x2");
        lex comp_lex;  // lex = lex_<less>
        polynomial_<Zp, lex> f(&comp_lex);
        {
            basic_monomial<lex> m(&comp_lex);
            m.push_back({x1, 2}); m.push_back({x2, 1}); m.normalization();
            f.push_back({m, Zp(3, p)});
        }
        {
            basic_monomial<lex> m(&comp_lex);
            m.push_back({x1, 1}); m.normalization();
            f.push_back({m, Zp(2, p)});
        }
        {
            basic_monomial<lex> m(&comp_lex);
            m.push_back({x2, 2}); m.normalization();
            f.push_back({m, Zp(1, p)});
        }
        f.normalization();

        std::vector<variable> aux_vars    = { x2 };
        std::vector<Zp>       sparse_betas = { Zp(3, p) };
        std::vector<upolynomial_<Zp>> images;

        __si_theta_array_eval(f, x1, aux_vars, sparse_betas, 3, images);

        CLPOLY_ASSERT_EQ((int)images.size(), 3);

        // l=1: 9*x1^2 + 2*x1 + 9
        CLPOLY_ASSERT_EQ((int)get_coeff(images[0], 2), 9);
        CLPOLY_ASSERT_EQ((int)get_coeff(images[0], 1), 2);
        CLPOLY_ASSERT_EQ((int)get_coeff(images[0], 0), 9);

        // l=2: 27*x1^2 + 2*x1 + 81
        CLPOLY_ASSERT_EQ((int)get_coeff(images[1], 2), 27);
        CLPOLY_ASSERT_EQ((int)get_coeff(images[1], 1), 2);
        CLPOLY_ASSERT_EQ((int)get_coeff(images[1], 0), 81);

        // l=3: 81*x1^2 + 2*x1 + 22
        CLPOLY_ASSERT_EQ((int)get_coeff(images[2], 2), 81);
        CLPOLY_ASSERT_EQ((int)get_coeff(images[2], 1), 2);
        CLPOLY_ASSERT_EQ((int)get_coeff(images[2], 0), 22);
    }

    // ---- 三变量（x1 主变量，x2,x3 辅助变量）-----------------------
    CLPOLY_TEST("theta_array_eval_trivar");
    {
        // f = x1*x2*x3 + x2^2 + x3 ∈ Zp[x1,x2,x3]
        // aux_vars=[x2,x3], β2=2, β3=5, s=2, p=97
        // θ_m:
        //   x1*x2*x3: θ = β2^1 * β3^1 = 2*5 = 10
        //   x2^2:     θ = β2^2 * β3^0 = 4
        //   x3:       θ = β2^0 * β3^1 = 5
        // l=1:
        //   1*(10^1)*x1 + 1*(4^1)*x1^0 + 1*(5^1)*x1^0
        //   = 10*x1 + (4+5) = 10*x1 + 9
        // l=2:
        //   1*(10^2)*x1 + 1*(4^2)*x1^0 + 1*(5^2)*x1^0
        //   = 100*x1 + (16+25)
        //   100 mod 97 = 3, 41 mod 97 = 41
        //   = 3*x1 + 41
        uint64_t p = 97;
        variable x1("x1"), x2("x2"), x3("x3");
        lex comp_lex;
        polynomial_<Zp, lex> f(&comp_lex);
        {
            basic_monomial<lex> m(&comp_lex);
            m.push_back({x1, 1}); m.push_back({x2, 1}); m.push_back({x3, 1}); m.normalization();
            f.push_back({m, Zp(1, p)});
        }
        {
            basic_monomial<lex> m(&comp_lex);
            m.push_back({x2, 2}); m.normalization();
            f.push_back({m, Zp(1, p)});
        }
        {
            basic_monomial<lex> m(&comp_lex);
            m.push_back({x3, 1}); m.normalization();
            f.push_back({m, Zp(1, p)});
        }
        f.normalization();

        std::vector<variable> aux_vars    = { x2, x3 };
        std::vector<Zp>       sparse_betas = { Zp(2, p), Zp(5, p) };
        std::vector<upolynomial_<Zp>> images;

        __si_theta_array_eval(f, x1, aux_vars, sparse_betas, 2, images);

        CLPOLY_ASSERT_EQ((int)images.size(), 2);

        // l=1: 10*x1 + 9
        CLPOLY_ASSERT_EQ((int)get_coeff(images[0], 1), 10);
        CLPOLY_ASSERT_EQ((int)get_coeff(images[0], 0), 9);

        // l=2: 3*x1 + 41  (100 mod 97 = 3)
        CLPOLY_ASSERT_EQ((int)get_coeff(images[1], 1), 3);
        CLPOLY_ASSERT_EQ((int)get_coeff(images[1], 0), 41);
    }

    return clpoly_test::test_summary();
}
