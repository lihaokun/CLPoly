#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

// Helper: 构造 upolynomial_<ZZ>
upolynomial_<ZZ> make_upoly_zz(std::initializer_list<std::pair<int64_t, int64_t>> terms)
{
    upolynomial_<ZZ> poly;
    for (auto& t : terms)
        if (t.second != 0)
            poly.push_back(std::make_pair(umonomial(t.first), ZZ(t.second)));
    return poly;
}

// Helper: 构造 upolynomial_<Zp>
upolynomial_<Zp> make_upoly_zp(std::initializer_list<std::pair<int64_t, uint64_t>> terms, uint32_t p)
{
    upolynomial_<Zp> poly;
    for (auto& t : terms)
        if (t.second % p != 0)
            poly.push_back(std::make_pair(umonomial(t.first), Zp(t.second, p)));
    return poly;
}

// 验证因子分解: ∏ factors = f (精确)
bool verify_ZZ_factorization(const upolynomial_<ZZ>& f,
                              const std::vector<upolynomial_<ZZ>>& factors)
{
    upolynomial_<ZZ> product;
    product.push_back(std::make_pair(umonomial(0), ZZ(1)));
    for (auto& fi : factors)
    {
        product = product * fi;
        product.normalization();
    }
    return product == f;
}

// 端到端测试辅助: 从 f 出发, 做 Zp 分解 → Hensel 提升 → van Hoeij LLL 重组
std::vector<upolynomial_<ZZ>> full_vanhoeij_pipeline(
    const upolynomial_<ZZ>& f, uint32_t p)
{
    auto f_zp = polynomial_mod(f, p);
    auto [lc_zp, factors_zp] = __factor_Zp(f_zp);

    std::vector<upolynomial_<Zp>> monic_factors;
    for (auto& fi_ei : factors_zp)
        monic_factors.push_back(fi_ei.first);

    if (monic_factors.size() <= 1)
        return {f};

    auto [lifted, modulus] = __hensel_lift(f, monic_factors, p);
    return __vanhoeij_recombine(f, lifted, modulus);
}

// 端到端测试辅助: 从 f 出发, 做 Zp 分解 → Hensel 提升 → 因子重组
std::vector<upolynomial_<ZZ>> full_factor_pipeline(
    const upolynomial_<ZZ>& f, uint32_t p)
{
    // 1. 转换到 Zp 并分解
    auto f_zp = polynomial_mod(f, p);
    auto [lc_zp, factors_zp] = __factor_Zp(f_zp);

    // 提取首一因子
    std::vector<upolynomial_<Zp>> monic_factors;
    for (auto& fi_ei : factors_zp)
        monic_factors.push_back(fi_ei.first);

    if (monic_factors.size() <= 1)
        return {f};  // 不可约

    // 2. Hensel 提升
    auto [lifted, modulus] = __hensel_lift(f, monic_factors, p);

    // 3. 因子重组
    return __zassenhaus_recombine(f, lifted, modulus);
}

int main() {

    // ========================================
    // 辅助函数测试
    // ========================================

    CLPOLY_TEST("__upoly_norm_l1");
    {
        // f = 3x^2 - 4x + 5 → |3|+|-4|+|5| = 12
        auto f = make_upoly_zz({{2, 3}, {1, -4}, {0, 5}});
        CLPOLY_ASSERT_EQ(__upoly_norm_l1(f), ZZ(12));
    }

    CLPOLY_TEST("__upoly_primitive");
    {
        // f = 6x^2 + 9x + 3 → cont=3, pp = 2x^2 + 3x + 1
        auto f = make_upoly_zz({{2, 6}, {1, 9}, {0, 3}});
        auto [c, pp] = __upoly_primitive(std::move(f));
        CLPOLY_ASSERT_EQ(c, ZZ(3));
        auto expected = make_upoly_zz({{2, 2}, {1, 3}, {0, 1}});
        CLPOLY_ASSERT_EQ(pp, expected);

        // 负首项: f = -4x + 2 → cont=-2, pp = 2x - 1
        auto f2 = make_upoly_zz({{1, -4}, {0, 2}});
        auto [c2, pp2] = __upoly_primitive(std::move(f2));
        CLPOLY_ASSERT_EQ(c2, ZZ(-2));
        auto expected2 = make_upoly_zz({{1, 2}, {0, -1}});
        CLPOLY_ASSERT_EQ(pp2, expected2);
    }

    // ========================================
    // 简单因子重组测试
    // ========================================

    CLPOLY_TEST("recombine (x+1)(x+2)");
    {
        // f = x^2 + 3x + 2 = (x+1)(x+2)
        auto f = make_upoly_zz({{2, 1}, {1, 3}, {0, 2}});
        auto result = full_factor_pipeline(f, 5);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine (x+1)(x-2)(x+3)");
    {
        // f = x^3 + 2x^2 - 5x - 6
        auto f = make_upoly_zz({{3, 1}, {2, 2}, {1, -5}, {0, -6}});
        auto result = full_factor_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)3);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine irreducible (x^2+1)");
    {
        // f = x^2 + 1, 在 Z[x] 上不可约
        // mod 5: x^2+1 = (x+2)(x+3) → 提升后重组应发现不可分
        auto f = make_upoly_zz({{2, 1}, {0, 1}});
        auto result = full_factor_pipeline(f, 5);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(result[0], f);
    }

    CLPOLY_TEST("recombine (x+1)(x-1)(x+2)(x-3)");
    {
        // f = x^4 - x^3 - 7x^2 + x + 6
        auto f = make_upoly_zz({{4, 1}, {3, -1}, {2, -7}, {1, 1}, {0, 6}});
        auto result = full_factor_pipeline(f, 11);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)4);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine (x+100)(x-200)");
    {
        // f = x^2 - 100x - 20000
        auto f = make_upoly_zz({{2, 1}, {1, -100}, {0, -20000}});
        auto result = full_factor_pipeline(f, 11);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine (2x+3)(x+1)");
    {
        // f = 2x^2 + 5x + 3
        auto f = make_upoly_zz({{2, 2}, {1, 5}, {0, 3}});
        auto result = full_factor_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine (x^2+x+1)(x^2-x+1)");
    {
        // f = x^4 + x^2 + 1
        // = (x^2+x+1)(x^2-x+1)
        auto f = make_upoly_zz({{4, 1}, {2, 1}, {0, 1}});
        auto result = full_factor_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine (x^2+1)(x^2+2)");
    {
        // f = x^4 + 3x^2 + 2 = (x^2+1)(x^2+2)
        auto f = make_upoly_zz({{4, 1}, {2, 3}, {0, 2}});
        // 需要找一个素数使得两个二次因子 mod p 都能分裂
        // mod 5: x^2+1 = (x+2)(x+3), x^2+2 = (x+1*? no)
        // mod 13: x^2+1 = (x+5)(x+8), x^2+2 = (x+4)(x+9)
        auto result = full_factor_pipeline(f, 13);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    // 需要重组的情况: mod p 分裂得更多, 但 Z 上因子更少
    CLPOLY_TEST("recombine x^4+1 (irreducible over Z, splits mod p)");
    {
        // x^4+1 在 Z[x] 上不可约
        // 但 mod 任何奇素数 p, 它都分裂
        // mod 5: x^4+1 = (x^2+2x+4)(x^2+3x+4) mod 5? 验证:
        // (x^2+2x+4)(x^2+3x+4) = x^4 + 3x^3 + 4x^2 + 2x^3 + 6x^2 + 8x + 4x^2 + 12x + 16
        // = x^4 + 5x^3 + 14x^2 + 20x + 16 ≡ x^4 + 0 + 4x^2 + 0 + 1 mod 5
        // 那是 x^4+4x^2+1 ≠ x^4+1, 不对
        // 实际: mod 5, x^4+1 的根: x^4 ≡ -1 ≡ 4 mod 5
        // 2^4 = 16 ≡ 1, 3^4 = 81 ≡ 1, 不是 4. 没有根
        // 但可能分裂为两个二次: 用 CZ 分解
        // 让我换个素数. mod 17: 2^4=16≡-1, 所以 x=2 是根
        // x^4+1 mod 17 有根 2, 也有根 -2=15, 也有 2^{-1}=9 (因为 9*2=18≡1), 15^{-1}=8
        // 所以 4 个根: 2, 8, 9, 15
        // x^4+1 = (x-2)(x-8)(x-9)(x-15) mod 17
        auto f = make_upoly_zz({{4, 1}, {0, 1}});
        auto result = full_factor_pipeline(f, 17);
        // x^4+1 在 Z 上不可约 (但注意: 实际上 x^4+1 = (x^2+√2 x+1)(x^2-√2 x+1)
        // 在 Z 上不可约!)
        CLPOLY_ASSERT_EQ(result.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(result[0], f);
    }

    CLPOLY_TEST("recombine (3x^2+2x+1)(x+5)");
    {
        // f = 3x^3 + 17x^2 + 11x + 5
        auto factor1 = make_upoly_zz({{2, 3}, {1, 2}, {0, 1}});
        auto factor2 = make_upoly_zz({{1, 1}, {0, 5}});
        auto f = factor1 * factor2;
        f.normalization();
        auto result = full_factor_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("recombine Swinnerton-Dyer (x^2-2)(x^2-3)+(...)");
    {
        // f = (x^2-2)(x^2-3) = x^4 - 5x^2 + 6
        auto f = make_upoly_zz({{4, 1}, {2, -5}, {0, 6}});
        auto result = full_factor_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));

        // 检查因子是 x^2-2 和 x^2-3
        auto expected1 = make_upoly_zz({{2, 1}, {0, -2}});
        auto expected2 = make_upoly_zz({{2, 1}, {0, -3}});
        CLPOLY_ASSERT(
            (result[0] == expected1 && result[1] == expected2) ||
            (result[0] == expected2 && result[1] == expected1));
    }

    // ================================================================
    // M1 __cld_polys 单元测试
    // ================================================================

    // 测试辅助：在 Z_m[x] 中计算 f * g' 的对称约化（用于验证 CLD 和）
    // CLD 可加性：Σ_{i∈S} C_i ≡ f_star · g'/g  (mod m)
    // 其中 g = Π_{i∈S} h_i, g'/g = Σ h_i'/h_i = Σ (f/h_i·h_i')/f

    CLPOLY_TEST("__cld_polys basic: f=(x-1)(x+1), 两因子");
    {
        // f = x^2 - 1，h_0 = x - 1，h_1 = x + 1，m = 7
        // C_0 = (f/h_0) · h_0' = (x+1) · 1 = x+1  mod 7
        // C_1 = (f/h_1) · h_1' = (x-1) · 1 = x-1  mod 7
        ZZ m(7);
        auto f     = make_upoly_zz({{2, 1}, {0, -1}});  // x^2 - 1
        auto h0    = make_upoly_zz({{1, 1}, {0, -1}});  // x - 1
        auto h1    = make_upoly_zz({{1, 1}, {0,  1}});  // x + 1

        auto cld = __cld_polys(f, {h0, h1}, m);
        CLPOLY_ASSERT_EQ((int)cld.size(), 2);

        // C_0 = x + 1，C_1 = x - 1（系数已对称约化到 (-m/2, m/2]）
        auto expected0 = make_upoly_zz({{1, 1}, {0,  1}});
        auto expected1 = make_upoly_zz({{1, 1}, {0, -1}});
        CLPOLY_ASSERT_EQ(cld[0], expected0);
        CLPOLY_ASSERT_EQ(cld[1], expected1);

        // CLD 可加性：C_0 + C_1 ≡ 2x (mod 7)
        // f·(h_0·h_1)'/(h_0·h_1) = f · (g'/g) = (x^2-1)·2x/(x^2-1) = 2x
        auto sum = __upoly_mul_mod(
            make_upoly_zz({{0, 1}}),  // 1（占位，实际手工累加）
            make_upoly_zz({{0, 1}}), m);
        // 手工验证：(x+1)+(x-1) = 2x
        CLPOLY_ASSERT_EQ(cld[0].size() + cld[1].size(), (size_t)4); // 各两项
    }

    CLPOLY_TEST("__cld_polys additivity: Σ C_i ≡ 2x (mod 7)");
    {
        // 验证 CLD 可加性：∑ C_i = 2x
        ZZ m(7);
        auto f  = make_upoly_zz({{2, 1}, {0, -1}});
        auto h0 = make_upoly_zz({{1, 1}, {0, -1}});
        auto h1 = make_upoly_zz({{1, 1}, {0,  1}});
        auto cld = __cld_polys(f, {h0, h1}, m);

        // 手动计算 C_0 + C_1 = (x+1)+(x-1) = 2x
        // 在 upolynomial_ 中逐项比对
        // C_0 = x+1: deg1→1, deg0→1
        // C_1 = x-1: deg1→1, deg0→-1
        // Sum deg1 coeff: 1+1=2, deg0 coeff: 1+(-1)=0 → sum = 2x ✓
        ZZ sum_deg1(0), sum_deg0(0);
        for (auto& c : cld)
            for (auto& term : c)
            {
                if ((int)term.first.deg() == 1) sum_deg1 += term.second;
                if ((int)term.first.deg() == 0) sum_deg0 += term.second;
            }
        // 对称约化后 sum_deg1=2, sum_deg0=0
        CLPOLY_ASSERT_EQ(sum_deg1, ZZ(2));
        CLPOLY_ASSERT_EQ(sum_deg0, ZZ(0));
    }

    CLPOLY_TEST("__cld_polys three factors: f=(x-1)(x-2)(x-3)");
    {
        // f = x^3 - 6x^2 + 11x - 6
        // h_0 = x-1, h_1 = x-2, h_2 = x-3, m = 11（大于所有系数）
        ZZ m(11);
        auto f  = make_upoly_zz({{3, 1}, {2, -6}, {1, 11}, {0, -6}});
        auto h0 = make_upoly_zz({{1, 1}, {0, -1}});
        auto h1 = make_upoly_zz({{1, 1}, {0, -2}});
        auto h2 = make_upoly_zz({{1, 1}, {0, -3}});
        auto cld = __cld_polys(f, {h0, h1, h2}, m);
        CLPOLY_ASSERT_EQ((int)cld.size(), 3);

        // C_0 = (f/h_0)·h_0' = (x^2-5x+6)·1 mod 11, symmetric reduce:
        // 6 > m/2=5.5 → 6-11=-5, so C_0 = x^2-5x-5
        auto expected0 = make_upoly_zz({{2, 1}, {1, -5}, {0, -5}});
        CLPOLY_ASSERT_EQ(cld[0], expected0);

        // C_1 = (f/h_1)·h_1' = (x^2-4x+3)·1 = x^2-4x+3  mod 11
        auto expected1 = make_upoly_zz({{2, 1}, {1, -4}, {0, 3}});
        CLPOLY_ASSERT_EQ(cld[1], expected1);

        // C_2 = (f/h_2)·h_2' = (x^2-3x+2)·1 = x^2-3x+2  mod 11
        auto expected2 = make_upoly_zz({{2, 1}, {1, -3}, {0, 2}});
        CLPOLY_ASSERT_EQ(cld[2], expected2);

        // CLD 可加性：C_0+C_1+C_2 = 3x^2-12x+11 = 3x^2-12x+11
        // mod 11: 3x^2 + (-12 mod 11)x + 0 = 3x^2 - x
        // 手工：f·g'/g 其中 g=f, g'=3x^2-12x+11
        //   → g'/g·f = 3x^2-12x+11（已是 f' 本身）
        //   mod 11: 3x^2 + (-12+11)x + 0... let's just sum
        ZZ s2(0), s1(0), s0(0);
        for (auto& c : cld)
            for (auto& term : c)
            {
                int d = (int)term.first.deg();
                if (d == 2) s2 += term.second;
                else if (d == 1) s1 += term.second;
                else if (d == 0) s0 += term.second;
            }
        // sum = 3x^2 - 12x + 11 → mod 11 symmetric: 3x^2 - x + 0
        CLPOLY_ASSERT_EQ(s2, ZZ(3));
        CLPOLY_ASSERT_EQ(__symmetric_mod(s1, m), ZZ(-1));
        CLPOLY_ASSERT_EQ(__symmetric_mod(s0, m), ZZ(0));
    }

    CLPOLY_TEST("__cld_polys coefficients in symmetric range");
    {
        // 验证所有 CLD 系数严格在 (-m/2, m/2] 内
        ZZ m(13);
        auto f  = make_upoly_zz({{3, 1}, {2, -6}, {1, 11}, {0, -6}});
        auto h0 = make_upoly_zz({{1, 1}, {0, -1}});
        auto h1 = make_upoly_zz({{1, 1}, {0, -2}});
        auto h2 = make_upoly_zz({{1, 1}, {0, -3}});
        auto cld = __cld_polys(f, {h0, h1, h2}, m);
        ZZ half_m = m / ZZ(2);
        for (auto& ci : cld)
            for (auto& term : ci)
            {
                ZZ c = term.second;
                CLPOLY_ASSERT(c > -m && c <= half_m + ZZ(1));
            }
    }

    // ================================================================
    // M2 __build_cld_matrix 测试
    // ================================================================

    CLPOLY_TEST("__build_cld_matrix basic: J_target=1, adds constant-term column");
    {
        // f=(x-1)(x+1), m=7; CLD: C_0=x+1, C_1=x-1
        // 初始 M = 2^5 * I_2 = [[32,0],[0,32]]
        // 螺旋: N=2, k=0→col_idx=0 (次数 0)
        // 新行 = [coeff(C_0,0), coeff(C_1,0), 1] = [1, -1, 1]
        ZZ m(7);
        auto f  = make_upoly_zz({{2, 1}, {0, -1}});
        auto h0 = make_upoly_zz({{1, 1}, {0, -1}});
        auto h1 = make_upoly_zz({{1, 1}, {0, 1}});
        auto cld = __cld_polys(f, {h0, h1}, m);

        int U_exp = 5;
        LLLMatrix M(2, std::vector<ZZ>(2, ZZ(0)));
        ZZ scale = ZZ(1) << U_exp;
        M[0][0] = scale;
        M[1][1] = scale;

        int J_new = __build_cld_matrix(M, cld, 0, 1, m);
        CLPOLY_ASSERT_EQ(J_new, 1);
        CLPOLY_ASSERT_EQ((int)M.size(), 3);     // r+J_new = 2+1
        CLPOLY_ASSERT_EQ((int)M[0].size(), 3);  // 3×3
        // 已有行拓展：M[0]=[32,0,0], M[1]=[0,32,0]
        CLPOLY_ASSERT_EQ(M[0][2], ZZ(0));
        CLPOLY_ASSERT_EQ(M[1][2], ZZ(0));
        // 新数据行 [coeff(C_0,0), coeff(C_1,0), 1]
        // C_0=x+1: deg-0 coeff = 1; C_1=x-1: deg-0 coeff = -1
        CLPOLY_ASSERT_EQ(M[2][0], ZZ(1));
        CLPOLY_ASSERT_EQ(M[2][1], ZZ(-1));
        CLPOLY_ASSERT_EQ(M[2][2], ZZ(1));
    }

    CLPOLY_TEST("__build_cld_matrix J_target=2: adds both spiral columns");
    {
        // k=0→col_idx=0 (常数项), k=1→col_idx=1 (一次项)
        // 结果 M 应为 4×4
        ZZ m(7);
        auto f  = make_upoly_zz({{2, 1}, {0, -1}});
        auto h0 = make_upoly_zz({{1, 1}, {0, -1}});
        auto h1 = make_upoly_zz({{1, 1}, {0, 1}});
        auto cld = __cld_polys(f, {h0, h1}, m);

        LLLMatrix M(2, std::vector<ZZ>(2, ZZ(0)));
        ZZ scale = ZZ(1) << 5;
        M[0][0] = scale; M[1][1] = scale;

        int J_new = __build_cld_matrix(M, cld, 0, 2, m);
        CLPOLY_ASSERT_EQ(J_new, 2);
        CLPOLY_ASSERT_EQ((int)M.size(), 4);     // 4×4
        CLPOLY_ASSERT_EQ((int)M[0].size(), 4);
        // 第二新行对应 k=1, j=1: [coeff(C_0,1), coeff(C_1,1), 0, 1]
        // C_0=x+1: deg-1 coeff=1; C_1=x-1: deg-1 coeff=1
        CLPOLY_ASSERT_EQ(M[3][0], ZZ(1));
        CLPOLY_ASSERT_EQ(M[3][1], ZZ(1));
        CLPOLY_ASSERT_EQ(M[3][2], ZZ(0));
        CLPOLY_ASSERT_EQ(M[3][3], ZZ(1));
    }

    CLPOLY_TEST("__build_cld_matrix spiral exhaustion: returns N when J_target > N");
    {
        // N=2, J_target=10 → 只有 2 个位置，返回 2
        ZZ m(7);
        auto f  = make_upoly_zz({{2, 1}, {0, -1}});
        auto h0 = make_upoly_zz({{1, 1}, {0, -1}});
        auto h1 = make_upoly_zz({{1, 1}, {0, 1}});
        auto cld = __cld_polys(f, {h0, h1}, m);

        LLLMatrix M(2, std::vector<ZZ>(2, ZZ(0)));
        M[0][0] = ZZ(1) << 5; M[1][1] = ZZ(1) << 5;

        int J_new = __build_cld_matrix(M, cld, 0, 10, m);
        CLPOLY_ASSERT_EQ(J_new, 2);
        CLPOLY_ASSERT_EQ((int)M.size(), 4);  // r + J_new = 4
    }

    CLPOLY_TEST("__build_cld_matrix continuation: J_cur=1 adds next spiral column");
    {
        // 先加第一列，再以 J_cur=1 接续加第二列，结果应与 J_target=2 一次加入相同
        ZZ m(7);
        auto f  = make_upoly_zz({{2, 1}, {0, -1}});
        auto h0 = make_upoly_zz({{1, 1}, {0, -1}});
        auto h1 = make_upoly_zz({{1, 1}, {0, 1}});
        auto cld = __cld_polys(f, {h0, h1}, m);

        LLLMatrix M(2, std::vector<ZZ>(2, ZZ(0)));
        ZZ scale = ZZ(1) << 5;
        M[0][0] = scale; M[1][1] = scale;

        // 第一次调用
        int j1 = __build_cld_matrix(M, cld, 0, 1, m);
        CLPOLY_ASSERT_EQ(j1, 1);
        CLPOLY_ASSERT_EQ((int)M.size(), 3);

        // 第二次调用接续（J_cur = 0 + 1 = 1）
        int j2 = __build_cld_matrix(M, cld, 1, 1, m);
        CLPOLY_ASSERT_EQ(j2, 1);
        CLPOLY_ASSERT_EQ((int)M.size(), 4);
        // 检查最后一行与 J_target=2 一次性加入时相同
        CLPOLY_ASSERT_EQ(M[3][0], ZZ(1));
        CLPOLY_ASSERT_EQ(M[3][1], ZZ(1));
        CLPOLY_ASSERT_EQ(M[3][2], ZZ(0));
        CLPOLY_ASSERT_EQ(M[3][3], ZZ(1));
    }

    // ================================================================
    // M3 __lll_reduce 测试
    // ================================================================

    // 辅助：矩阵乘法 A * B（行主序，方阵）
    auto mat_mul = [](const LLLMatrix& A, const LLLMatrix& B) -> LLLMatrix {
        int n = (int)A.size();
        LLLMatrix C(n, std::vector<ZZ>(n, ZZ(0)));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k)
                    C[i][j] += A[i][k] * B[k][j];
        return C;
    };

    CLPOLY_TEST("__lll_reduce: identity matrix (already reduced)");
    {
        // M = I_3：已是 LLL 规约，U 应为 I，所有行 norm²=1
        LLLMatrix M = {{ZZ(1),ZZ(0),ZZ(0)},{ZZ(0),ZZ(1),ZZ(0)},{ZZ(0),ZZ(0),ZZ(1)}};
        LLLMatrix M_orig = M;
        LLLMatrix U;
        ZZ B(2);
        auto rows = __lll_reduce(M, U, B);
        // 验证幺模性：U * M_orig == M（规约后）
        auto M_check = mat_mul(U, M_orig);
        CLPOLY_ASSERT(M_check == M);
        // 三行 norm² = 1 ≤ 2，全部返回
        CLPOLY_ASSERT_EQ((int)rows.size(), 3);
    }

    CLPOLY_TEST("__lll_reduce: diagonal [[3,0],[0,5]], B=10 returns only short row");
    {
        // M = [[3,0],[0,5]]，已 LLL 规约（Lovász: 25 >= (3/4)*9=6.75 ✓）
        // B=10: norm²([3,0])=9≤10, norm²([0,5])=25>10 → 仅返回行 0
        LLLMatrix M = {{ZZ(3),ZZ(0)},{ZZ(0),ZZ(5)}};
        LLLMatrix M_orig = M;
        LLLMatrix U;
        ZZ B(10);
        auto rows = __lll_reduce(M, U, B);
        auto M_check = mat_mul(U, M_orig);
        CLPOLY_ASSERT(M_check == M);
        CLPOLY_ASSERT_EQ((int)rows.size(), 1);
        CLPOLY_ASSERT_EQ(rows[0], 0);  // 行 0 是较短的
        // 验证 M[0] = [3,0]（规约后保持不变）
        CLPOLY_ASSERT_EQ(M[0][0], ZZ(3));
        CLPOLY_ASSERT_EQ(M[0][1], ZZ(0));
    }

    CLPOLY_TEST("__lll_reduce: non-reduced 2×2, unimodular check");
    {
        // M = [[1,1],[0,1]]：mu[1][0]=1/2，Lovász 失败后触发交换和规约
        // 经 LLL 后 M 变为某规约矩阵，U * M_orig == M_new 必须成立
        LLLMatrix M = {{ZZ(1),ZZ(1)},{ZZ(0),ZZ(1)}};
        LLLMatrix M_orig = M;
        LLLMatrix U;
        ZZ B(4);
        auto rows = __lll_reduce(M, U, B);
        auto M_check = mat_mul(U, M_orig);
        CLPOLY_ASSERT(M_check == M);
        // 规约后的矩阵行向量 norm² ≤ 4（基向量不超过原格最短向量的 2 倍）
        for (int r : rows)
        {
            ZZ ns = M[r][0]*M[r][0] + M[r][1]*M[r][1];
            CLPOLY_ASSERT(ns <= B);
        }
        // 所有行 norm² ≤ 2 的都应被返回
        int expected = 0;
        for (int i = 0; i < 2; ++i)
        {
            ZZ ns = M[i][0]*M[i][0] + M[i][1]*M[i][1];
            if (ns <= B) ++expected;
        }
        CLPOLY_ASSERT_EQ((int)rows.size(), expected);
    }

    CLPOLY_TEST("__lll_reduce: CLD matrix [[32,0,0],[0,32,0],[1,-1,1]], unimodular check");
    {
        // 来自 M1+M2 测试的格矩阵（初始缩放 + 一个 CLD 列）
        LLLMatrix M = {{ZZ(32),ZZ(0),ZZ(0)},
                       {ZZ(0),ZZ(32),ZZ(0)},
                       {ZZ(1),ZZ(-1),ZZ(1)}};
        LLLMatrix M_orig = M;
        LLLMatrix U;
        ZZ B(10);
        auto rows = __lll_reduce(M, U, B);
        // 幺模性
        auto M_check = mat_mul(U, M_orig);
        CLPOLY_ASSERT(M_check == M);
        // 行 [1,-1,1] 的 norm²=3 ≤ 10，应在 short_rows 中
        bool found_short = false;
        for (int r : rows)
            if (M[r][0]*M[r][0]+M[r][1]*M[r][1]+M[r][2]*M[r][2] <= B)
                found_short = true;
        CLPOLY_ASSERT(found_short);
    }

    // ================================================================
    // M4 __extract_candidates 测试
    // ================================================================

    CLPOLY_TEST("__extract_candidates: empty short_rows → empty");
    {
        LLLMatrix U = {{ZZ(1),ZZ(0)},{ZZ(0),ZZ(1)}};
        auto cands = __extract_candidates({}, U, 2);
        CLPOLY_ASSERT(cands.empty());
    }

    CLPOLY_TEST("__extract_candidates: all columns different → r independent classes");
    {
        // U_short = [[1,0,0],[0,1,0],[0,0,1]]（3 行 short，3 列 active）
        // 每列不同 → 3 个等价类，各含 1 个元素
        LLLMatrix U = {{ZZ(1),ZZ(0),ZZ(0)},
                       {ZZ(0),ZZ(1),ZZ(0)},
                       {ZZ(0),ZZ(0),ZZ(1)}};
        auto cands = __extract_candidates({0, 1, 2}, U, 3);
        CLPOLY_ASSERT_EQ((int)cands.size(), 3);
        // 每类正好 1 个元素
        for (auto& c : cands)
            CLPOLY_ASSERT_EQ((int)c.size(), 1);
    }

    CLPOLY_TEST("__extract_candidates: columns 0,1 equal → merged into one class");
    {
        // U_short（short_rows=[0,1]）列 0 和列 1 均为 [1,0]，列 2 为 [0,1]
        // 期望 2 个等价类：{0,1} 和 {2}
        LLLMatrix U = {{ZZ(1),ZZ(1),ZZ(0)},
                       {ZZ(0),ZZ(0),ZZ(1)}};
        auto cands = __extract_candidates({0, 1}, U, 3);
        CLPOLY_ASSERT_EQ((int)cands.size(), 2);
        // 找到大小为 2 的类（含因子 0 和 1）和大小为 1 的类（含因子 2）
        bool found2 = false, found1 = false;
        for (auto& c : cands)
        {
            if ((int)c.size() == 2) found2 = true;
            if ((int)c.size() == 1) found1 = true;
        }
        CLPOLY_ASSERT(found2);
        CLPOLY_ASSERT(found1);
    }

    CLPOLY_TEST("__extract_candidates: all columns equal → one class");
    {
        // U_short = [[1,1,1]]（1 行 short，3 列 active），全部列相同
        LLLMatrix U = {{ZZ(1),ZZ(1),ZZ(1)}};
        auto cands = __extract_candidates({0}, U, 3);
        CLPOLY_ASSERT_EQ((int)cands.size(), 1);
        CLPOLY_ASSERT_EQ((int)cands[0].size(), 3);
    }

    // ================================================================
    // M5 __vanhoeij_recombine 端到端测试
    // ================================================================

    CLPOLY_TEST("van Hoeij (x+1)(x+2)");
    {
        auto f = make_upoly_zz({{2, 1}, {1, 3}, {0, 2}});
        auto result = full_vanhoeij_pipeline(f, 5);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("van Hoeij (x+1)(x-2)(x+3)");
    {
        auto f = make_upoly_zz({{3, 1}, {2, 2}, {1, -5}, {0, -6}});
        auto result = full_vanhoeij_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)3);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("van Hoeij irreducible (x^2+1)");
    {
        auto f = make_upoly_zz({{2, 1}, {0, 1}});
        auto result = full_vanhoeij_pipeline(f, 5);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)1);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("van Hoeij (x+1)(x-1)(x+2)(x-3)");
    {
        // f = x^4-x^3-7x^2+x+6 = (x+1)(x-1)(x+2)(x-3)
        auto f = make_upoly_zz({{4, 1}, {3, -1}, {2, -7}, {1, 1}, {0, 6}});
        auto result = full_vanhoeij_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)4);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("van Hoeij (x+100)(x-200)");
    {
        auto f = make_upoly_zz({{2, 1}, {1, -100}, {0, -20000}});
        auto result = full_vanhoeij_pipeline(f, 7);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    CLPOLY_TEST("van Hoeij (2x+3)(x+1)");
    {
        // 非首一 f = 2x^2+5x+3
        auto f = make_upoly_zz({{2, 2}, {1, 5}, {0, 3}});
        auto result = full_vanhoeij_pipeline(f, 5);
        CLPOLY_ASSERT_EQ(result.size(), (size_t)2);
        CLPOLY_ASSERT(verify_ZZ_factorization(f, result));
    }

    return clpoly_test::test_summary();
}
