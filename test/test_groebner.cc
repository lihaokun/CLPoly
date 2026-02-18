#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

int main() {
    variable x("x"), y("y"), z("z");

    // ======== M1: monomial LCM ========

    CLPOLY_TEST("lcm_basic");
    {
        monomial m1({{x, 3}, {y, 2}});
        monomial m2({{x, 1}, {y, 4}, {z, 2}});
        auto l = lcm(m1, m2);
        CLPOLY_ASSERT_EQ(l.deg(x), (int64_t)3);  // max(3,1)
        CLPOLY_ASSERT_EQ(l.deg(y), (int64_t)4);  // max(2,4)
        CLPOLY_ASSERT_EQ(l.deg(z), (int64_t)2);  // max(0,2)
        CLPOLY_ASSERT_EQ(l.deg(), (int64_t)9);    // 3+4+2
    }

    CLPOLY_TEST("lcm_divisibility");
    {
        // lcm(m1,m2) must be divisible by both m1 and m2
        monomial m1({{x, 3}, {y, 2}});
        monomial m2({{x, 1}, {y, 4}, {z, 2}});
        auto l = lcm(m1, m2);
        monomial q;
        CLPOLY_ASSERT_TRUE(is_divexact(q, l, m1));
        CLPOLY_ASSERT_TRUE(is_divexact(q, l, m2));
    }

    CLPOLY_TEST("lcm_consistent_with_gcd");
    {
        // lcm(m1,m2) == m1*m2/gcd(m1,m2)
        monomial m1({{x, 3}, {y, 2}});
        monomial m2({{x, 1}, {y, 4}, {z, 2}});
        auto l = lcm(m1, m2);
        auto g = gcd(m1, m2);
        auto expected = m1 * m2 / g;
        CLPOLY_ASSERT_EQ(l, expected);
    }

    CLPOLY_TEST("lcm_identity");
    {
        monomial m1({{x, 2}, {y, 3}});
        // lcm(m, 1) = m
        monomial one;
        CLPOLY_ASSERT_EQ(lcm(m1, one), m1);
        CLPOLY_ASSERT_EQ(lcm(one, m1), m1);
        // lcm(1, 1) = 1
        CLPOLY_ASSERT_EQ(lcm(one, one), one);
    }

    CLPOLY_TEST("lcm_self");
    {
        monomial m1({{x, 2}, {y, 3}});
        CLPOLY_ASSERT_EQ(lcm(m1, m1), m1);
    }

    CLPOLY_TEST("lcm_disjoint_variables");
    {
        // gcd = 1, lcm = m1*m2
        monomial m1({{x, 2}});
        monomial m2({{y, 3}});
        auto l = lcm(m1, m2);
        CLPOLY_ASSERT_EQ(l, m1 * m2);
        CLPOLY_ASSERT_EQ(l.deg(), (int64_t)5);
    }

    CLPOLY_TEST("lcm_coprime_detection");
    {
        // gcd(m1,m2) = 1 iff lcm.deg() == m1.deg() + m2.deg()
        monomial m1({{x, 2}});
        monomial m2({{y, 3}});
        auto l = lcm(m1, m2);
        CLPOLY_ASSERT_EQ(l.deg(), m1.deg() + m2.deg());  // coprime

        monomial m3({{x, 1}, {y, 1}});
        auto l2 = lcm(m1, m3);
        CLPOLY_ASSERT_TRUE(l2.deg() < m1.deg() + m3.deg());  // not coprime
    }

    CLPOLY_TEST("lcm_three_variables");
    {
        monomial m1({{x, 5}, {y, 1}, {z, 3}});
        monomial m2({{x, 2}, {y, 4}, {z, 1}});
        auto l = lcm(m1, m2);
        CLPOLY_ASSERT_EQ(l.deg(x), (int64_t)5);
        CLPOLY_ASSERT_EQ(l.deg(y), (int64_t)4);
        CLPOLY_ASSERT_EQ(l.deg(z), (int64_t)3);
        CLPOLY_ASSERT_EQ(l.deg(), (int64_t)12);
    }

    // ======== M2: S-polynomial ========

    CLPOLY_TEST("spoly_leading_term_cancels");
    {
        // S(f,g) should have its leading term cancelled
        // f = x^2 + y, g = x*y + x
        polynomial_QQ f = pow(x,2) + y;
        polynomial_QQ g = x*y + x;
        auto s = s_polynomial(f, g);
        // L = lcm(x^2, x*y) = x^2*y
        // S = (y)*f - (x)*g = x^2*y + y^2 - x^2*y - x^2 = y^2 - x^2
        polynomial_QQ expected = pow(y,2) - pow(x,2);
        CLPOLY_ASSERT_EQ(s, expected);
    }

    CLPOLY_TEST("spoly_same_polynomial");
    {
        // S(f,f) = 0
        polynomial_QQ f = pow(x,2) + y + 1;
        auto s = s_polynomial(f, f);
        CLPOLY_ASSERT_TRUE(s.empty());
    }

    CLPOLY_TEST("spoly_coprime_leading_terms");
    {
        // If LM(f) and LM(g) are coprime, S-poly is "large" (no simplification)
        // f = x^2, g = y^2
        polynomial_QQ f = pow(x,2);
        polynomial_QQ g = pow(y,2);
        auto s = s_polynomial(f, g);
        // L = x^2*y^2
        // S = y^2*f - x^2*g = x^2*y^2 - x^2*y^2 = 0
        CLPOLY_ASSERT_TRUE(s.empty());
    }

    CLPOLY_TEST("spoly_with_coefficients");
    {
        // f = 2x^2 + 3y, g = 5x*y + 7
        // S = (1/2)*y*(2x^2+3y) - (1/5)*x*(5xy+7)
        //   = x^2*y + 3y^2/2 - x^2*y - 7x/5
        //   = 3y^2/2 - 7x/5
        polynomial_QQ f = QQ(2)*pow(x,2) + QQ(3)*y;
        polynomial_QQ g = QQ(5)*x*y + QQ(7);
        auto s = s_polynomial(f, g);
        polynomial_QQ expected = QQ(3,2)*pow(y,2) - QQ(7,5)*x;
        CLPOLY_ASSERT_EQ(s, expected);
    }

    CLPOLY_TEST("spoly_linear");
    {
        // f = x + y, g = x - y
        // L = x, m_f = 1, m_g = 1
        // S = f - g = 2y
        polynomial_QQ f = x + y;
        polynomial_QQ g = x - y;
        auto s = s_polynomial(f, g);
        polynomial_QQ expected = QQ(2)*y;
        CLPOLY_ASSERT_EQ(s, expected);
    }

    CLPOLY_TEST("spoly_three_variables");
    {
        // f = x*y*z + x, g = x*y + z
        // LM(f)=xyz, LM(g)=xy, L=lcm=xyz
        // m_f=1, m_g=z
        // S = f - z*g = xyz + x - xyz - z^2 = x - z^2
        polynomial_QQ f = x*y*z + x;
        polynomial_QQ g = x*y + z;
        auto s = s_polynomial(f, g);
        polynomial_QQ expected = x - pow(z,2);
        CLPOLY_ASSERT_EQ(s, expected);
    }

    // ======== M3: Normal Form ========

    CLPOLY_TEST("nf_reduces_to_zero");
    {
        // NF(f, {f}) = 0
        polynomial_QQ f = pow(x,2) + y;
        std::vector<polynomial_QQ> G = {f};
        auto nf = normal_form(f, G);
        CLPOLY_ASSERT_TRUE(nf.empty());
    }

    CLPOLY_TEST("nf_constant_unchanged");
    {
        // NF(c, G) = c if c is a constant not divisible by any LM(G)
        polynomial_QQ c = polynomial_QQ(QQ(5));
        polynomial_QQ g = pow(x,2) + y;
        std::vector<polynomial_QQ> G = {g};
        auto nf = normal_form(c, G);
        CLPOLY_ASSERT_EQ(nf, c);
    }

    CLPOLY_TEST("nf_single_divisor");
    {
        // f = x^3 + x^2 + x + 1, G = {x^2 + 1}
        // x^3 + x^2 + x + 1 = x*(x^2+1) + (x^2+1) - x*(1) ... let's verify
        // x^3 / x^2 -> quot_m = x, c = 1
        // h = x^3+x^2+x+1 - x*(x^2+1) = x^2+x+1 - x = x^2+1... no
        // h = x^3+x^2+x+1 - x*(x^2+1) = x^3+x^2+x+1 - x^3 - x = x^2+1
        // x^2 / x^2 -> quot_m = 1, c = 1
        // h = x^2+1 - (x^2+1) = 0
        polynomial_QQ f = pow(x,3) + pow(x,2) + x + 1;
        polynomial_QQ g = pow(x,2) + 1;
        std::vector<polynomial_QQ> G = {g};
        auto nf = normal_form(f, G);
        CLPOLY_ASSERT_TRUE(nf.empty());
    }

    CLPOLY_TEST("nf_remainder");
    {
        // f = x^2 + x + 1, G = {x + 1}
        // x^2 / x -> quot_m = x, h = x^2+x+1 - x*(x+1) = x+1-x = 1... no
        // h = x^2+x+1 - x*(x+1) = x^2+x+1 - x^2 - x = 1
        polynomial_QQ f = pow(x,2) + x + 1;
        polynomial_QQ g = x + 1;
        std::vector<polynomial_QQ> G = {g};
        auto nf = normal_form(f, G);
        polynomial_QQ expected = polynomial_QQ(QQ(1));
        CLPOLY_ASSERT_EQ(nf, expected);
    }

    CLPOLY_TEST("nf_multivariate");
    {
        // f = x^2*y + x*y^2 + y^2
        // G = {x*y - 1, y^2 - 1}
        // grlex order: x^2*y > x*y^2 > y^2
        // step1: x^2*y / x*y -> quot_m=x, c=1, h = f - x*(xy-1) = f - x^2*y + x = x*y^2 + y^2 + x
        // step2: x*y^2 / x*y -> quot_m=y, c=1, h = x*y^2+y^2+x - y*(xy-1) = y^2+x+y
        // step3: y^2 / y^2 -> quot_m=1, c=1, h = y^2+x+y - (y^2-1) = x+y+1
        // step4: x not divisible by xy or y^2, move to r. h = y+1
        // step5: y not divisible by xy or y^2, move to r. h = 1
        // step6: 1 not divisible, move to r. h = 0
        // NF = x + y + 1
        polynomial_QQ f = pow(x,2)*y + x*pow(y,2) + pow(y,2);
        polynomial_QQ g1 = x*y - 1;
        polynomial_QQ g2 = pow(y,2) - 1;
        std::vector<polynomial_QQ> G = {g1, g2};
        auto nf = normal_form(f, G);
        polynomial_QQ expected = x + y + 1;
        CLPOLY_ASSERT_EQ(nf, expected);
    }

    CLPOLY_TEST("nf_empty_basis");
    {
        // NF(f, {}) = f
        polynomial_QQ f = pow(x,2) + y;
        std::vector<polynomial_QQ> G;
        auto nf = normal_form(f, G);
        CLPOLY_ASSERT_EQ(nf, f);
    }

    CLPOLY_TEST("nf_zero_polynomial");
    {
        // NF(0, G) = 0
        polynomial_QQ zero;
        polynomial_QQ g = pow(x,2) + 1;
        std::vector<polynomial_QQ> G = {g};
        auto nf = normal_form(zero, G);
        CLPOLY_ASSERT_TRUE(nf.empty());
    }

    // ======== M4+M5: Buchberger ========

    // 辅助: 检查 GB 基本性质
    // 1. 所有生成元对 GB 的 NF 为 0
    // 2. 每个 GB 元素首一
    // 3. 无冗余（无 LM(gi) | LM(gj)）
    auto check_gb = [](const std::vector<polynomial_QQ>& gb,
                       const std::vector<polynomial_QQ>& gens,
                       const std::string& label) {
        // 每个生成元 NF = 0
        for (size_t i = 0; i < gens.size(); ++i)
        {
            auto nf = normal_form(gens[i], gb);
            if (!nf.empty()) {
                std::cerr << "  " << label << ": gen[" << i << "] NF != 0" << std::endl;
                return false;
            }
        }
        // 首一
        for (size_t i = 0; i < gb.size(); ++i)
        {
            if (gb[i].empty() || gb[i].front().second != QQ(1)) {
                std::cerr << "  " << label << ": gb[" << i << "] not monic" << std::endl;
                return false;
            }
        }
        // 无冗余
        for (size_t i = 0; i < gb.size(); ++i)
            for (size_t j = 0; j < gb.size(); ++j) {
                if (i == j) continue;
                monomial tmp;
                if (is_divexact(tmp, gb[i].front().first, gb[j].front().first)) {
                    std::cerr << "  " << label << ": redundant LM" << std::endl;
                    return false;
                }
            }
        return true;
    };

    CLPOLY_TEST("gb_single_generator");
    {
        // GB({x^2+1}) = {x^2+1}
        polynomial_QQ f = pow(x,2) + 1;
        std::vector<polynomial_QQ> gens = {f};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(gb[0], f);
    }

    CLPOLY_TEST("gb_trivial_one");
    {
        // GB({1}) = {1}
        polynomial_QQ one = polynomial_QQ(QQ(1));
        std::vector<polynomial_QQ> gens = {one};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)1);
    }

    CLPOLY_TEST("gb_linear");
    {
        // GB({x+y-1, x-y}) grlex = {y - 1/2, x - 1/2}
        // Mathematica: {2y-1, 2x-1} -> monic: {y-1/2, x-1/2}
        polynomial_QQ f1 = x + y - 1;
        polynomial_QQ f2 = x - y;
        std::vector<polynomial_QQ> gens = {f1, f2};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)2);
        CLPOLY_ASSERT_TRUE(check_gb(gb, gens, "gb_linear"));
    }

    CLPOLY_TEST("gb_circle_hyperbola");
    {
        // I = <x^2+y^2-1, xy-1>
        // Mathematica grlex: {xy-1, x^2+y^2-1, y^3-y+x} (3 elements)
        polynomial_QQ f1 = pow(x,2) + pow(y,2) - 1;
        polynomial_QQ f2 = x*y - 1;
        std::vector<polynomial_QQ> gens = {f1, f2};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)3);
        CLPOLY_ASSERT_TRUE(check_gb(gb, gens, "gb_circle_hyperbola"));
    }

    CLPOLY_TEST("gb_x2_y_x3_x");
    {
        // I = <x^2-y, x^3-x>
        // Mathematica grlex: {y^2-y, xy-x, x^2-y} (3 elements)
        polynomial_QQ f1 = pow(x,2) - y;
        polynomial_QQ f2 = pow(x,3) - x;
        std::vector<polynomial_QQ> gens = {f1, f2};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)3);
        CLPOLY_ASSERT_TRUE(check_gb(gb, gens, "gb_x2_y_x3_x"));
    }

    CLPOLY_TEST("gb_symmetric_3var");
    {
        // I = <x+y+z, xy+yz+zx, xyz-1>
        // Mathematica grlex: {x+y+z, y^2+yz+z^2, z^3-1} (3 elements)
        polynomial_QQ f1 = x + y + z;
        polynomial_QQ f2 = x*y + y*z + z*x;
        polynomial_QQ f3 = x*y*z - 1;
        std::vector<polynomial_QQ> gens = {f1, f2, f3};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)3);
        CLPOLY_ASSERT_TRUE(check_gb(gb, gens, "gb_symmetric_3var"));
    }

    CLPOLY_TEST("gb_duplicate_generators");
    {
        // GB({f, f}) = GB({f})
        polynomial_QQ f = pow(x,2) + y;
        std::vector<polynomial_QQ> gens = {f, f};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(gb[0], f);
    }

    CLPOLY_TEST("gb_with_zero");
    {
        // 零多项式应被过滤
        polynomial_QQ f = pow(x,2) + 1;
        polynomial_QQ zero;
        std::vector<polynomial_QQ> gens = {zero, f, zero};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)1);
        CLPOLY_ASSERT_EQ(gb[0], f);
    }

    CLPOLY_TEST("gb_ideal_membership");
    {
        // f1, f2 ∈ <GB>, verify via NF
        polynomial_QQ f1 = pow(x,2) + pow(y,2) - 1;
        polynomial_QQ f2 = x*y - 1;
        std::vector<polynomial_QQ> gens = {f1, f2};
        auto gb = groebner_basis(gens);
        // x^3 - x should be in the ideal (can be verified)
        // More importantly: any linear combination should reduce to 0
        auto nf1 = normal_form(f1, gb);
        auto nf2 = normal_form(f2, gb);
        CLPOLY_ASSERT_TRUE(nf1.empty());
        CLPOLY_ASSERT_TRUE(nf2.empty());
        // Something NOT in the ideal should have nonzero NF
        polynomial_QQ g = x + y;
        auto nfg = normal_form(g, gb);
        CLPOLY_ASSERT_FALSE(nfg.empty());
    }

    CLPOLY_TEST("gb_exact_match_linear");
    {
        // GB({x+y-1, x-y}) = {x - 1/2, y - 1/2}
        polynomial_QQ f1 = x + y - 1;
        polynomial_QQ f2 = x - y;
        std::vector<polynomial_QQ> gens = {f1, f2};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)2);
        // 约化 GB 唯一，验证内容
        polynomial_QQ ex1 = polynomial_QQ(x) - QQ(1,2);
        polynomial_QQ ex2 = polynomial_QQ(y) - QQ(1,2);
        // gb 按某种顺序排列，检查集合相等
        bool match = (gb[0] == ex1 && gb[1] == ex2) ||
                     (gb[0] == ex2 && gb[1] == ex1);
        CLPOLY_ASSERT_TRUE(match);
    }

    // ======== M6: ZZ 入口 ========

    CLPOLY_TEST("gb_zz_basic");
    {
        // ZZ 输入: {2x + 3y, 4x - y}
        polynomial_ZZ f1 = 2*x + 3*y;
        polynomial_ZZ f2 = 4*x - y;
        std::vector<polynomial_ZZ> gens = {f1, f2};
        auto gb = groebner_basis(gens);
        // 结果应为 ZZ, 本原, lc > 0
        for (auto& g : gb)
        {
            CLPOLY_ASSERT_FALSE(g.empty());
            CLPOLY_ASSERT_TRUE(g.front().second > 0);
        }
        // 验证生成元在理想中（转 QQ 检查 NF）
        std::vector<polynomial_QQ> gb_qq;
        for (auto& g : gb) {
            polynomial_QQ gq;
            poly_convert(g, gq);
            gb_qq.push_back(gq);
        }
        for (auto& f : gens) {
            polynomial_QQ fq;
            poly_convert(f, fq);
            CLPOLY_ASSERT_TRUE(normal_form(fq, gb_qq).empty());
        }
    }

    CLPOLY_TEST("gb_zz_circle_hyperbola");
    {
        // ZZ 版本: {x^2+y^2-1, xy-1}
        polynomial_ZZ f1 = pow(x,2) + pow(y,2) - 1;
        polynomial_ZZ f2 = x*y - 1;
        std::vector<polynomial_ZZ> gens = {f1, f2};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)3);
        for (auto& g : gb) {
            CLPOLY_ASSERT_TRUE(g.front().second > 0);
            // 验证本原: 系数 GCD = 1
            ZZ c = abs(g.front().second);
            for (auto& term : g)
                c = gcd(c, abs(term.second));
            CLPOLY_ASSERT_EQ(c, ZZ(1));
        }
    }

    // ======== M6: Ideal 类 ========

    CLPOLY_TEST("ideal_groebner_basis");
    {
        polynomial_QQ f1 = pow(x,2) + pow(y,2) - 1;
        polynomial_QQ f2 = x*y - 1;
        Ideal<QQ> I = {f1, f2};
        auto& gb = I.groebner_basis();
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)3);
    }

    CLPOLY_TEST("ideal_contains");
    {
        polynomial_QQ f1 = pow(x,2) + pow(y,2) - 1;
        polynomial_QQ f2 = x*y - 1;
        Ideal<QQ> I = {f1, f2};
        CLPOLY_ASSERT_TRUE(I.contains(f1));
        CLPOLY_ASSERT_TRUE(I.contains(f2));
        CLPOLY_ASSERT_FALSE(I.contains(x + y));
    }

    CLPOLY_TEST("ideal_reduce");
    {
        polynomial_QQ f1 = pow(x,2) + pow(y,2) - 1;
        polynomial_QQ f2 = x*y - 1;
        Ideal<QQ> I = {f1, f2};
        auto r = I.reduce(f1);
        CLPOLY_ASSERT_TRUE(r.empty());
        auto r2 = I.reduce(x + y);
        CLPOLY_ASSERT_FALSE(r2.empty());
    }

    CLPOLY_TEST("ideal_zz");
    {
        polynomial_ZZ f1 = pow(x,2) - y;
        polynomial_ZZ f2 = pow(x,3) - x;
        Ideal<ZZ> I = {f1, f2};
        auto& gb = I.groebner_basis();
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)3);
    }

    // ======== 交叉测试: 与 Mathematica GroebnerBasis[] 精确对比 ========

    // 辅助: 将 GB vector 转为 set 比较（约化 GB 唯一，但 vector 顺序不定）
    auto gb_set_eq = [](const std::vector<polynomial_QQ>& gb,
                        const std::vector<polynomial_QQ>& expected) -> bool {
        if (gb.size() != expected.size()) return false;
        for (auto& e : expected) {
            bool found = false;
            for (auto& g : gb)
                if (g == e) { found = true; break; }
            if (!found) return false;
        }
        return true;
    };

    CLPOLY_TEST("cross_linear");
    {
        // Mathematica: GroebnerBasis[{x+y-1,x-y},{x,y},MonomialOrder->DegreeLexicographic]
        // = {-1+2y, -1+2x} → monic: {y-1/2, x-1/2}
        polynomial_QQ f1 = x + y - 1;
        polynomial_QQ f2 = x - y;
        auto gb = groebner_basis(std::vector<polynomial_QQ>{f1, f2});
        std::vector<polynomial_QQ> expected = {
            polynomial_QQ(x) - QQ(1,2),
            polynomial_QQ(y) - QQ(1,2)
        };
        CLPOLY_ASSERT_TRUE(gb_set_eq(gb, expected));
    }

    CLPOLY_TEST("cross_circle_hyperbola");
    {
        // Mathematica: {xy-1, x^2+y^2-1, x-y+y^3} (all monic)
        polynomial_QQ f1 = pow(x,2) + pow(y,2) - 1;
        polynomial_QQ f2 = x*y - 1;
        auto gb = groebner_basis(std::vector<polynomial_QQ>{f1, f2});
        std::vector<polynomial_QQ> expected = {
            x*y - 1,
            pow(x,2) + pow(y,2) - 1,
            pow(y,3) + x - y       // y^3 - y + x, grlex: y^3 > x > y
        };
        CLPOLY_ASSERT_TRUE(gb_set_eq(gb, expected));
    }

    CLPOLY_TEST("cross_x2y_x3x");
    {
        // Mathematica: {y^2-y, xy-x, x^2-y} (all monic)
        polynomial_QQ f1 = pow(x,2) - y;
        polynomial_QQ f2 = pow(x,3) - x;
        auto gb = groebner_basis(std::vector<polynomial_QQ>{f1, f2});
        std::vector<polynomial_QQ> expected = {
            pow(y,2) - y,
            x*y - x,
            pow(x,2) - y
        };
        CLPOLY_ASSERT_TRUE(gb_set_eq(gb, expected));
    }

    CLPOLY_TEST("cross_symmetric_3var");
    {
        // Mathematica: {x+y+z, -y^2-yz-z^2, 1-z^3}
        //   → monic: {x+y+z, y^2+yz+z^2, z^3-1}
        polynomial_QQ f1 = x + y + z;
        polynomial_QQ f2 = x*y + y*z + z*x;
        polynomial_QQ f3 = x*y*z - 1;
        auto gb = groebner_basis(std::vector<polynomial_QQ>{f1, f2, f3});
        std::vector<polynomial_QQ> expected = {
            x + y + z,
            pow(y,2) + y*z + pow(z,2),
            pow(z,3) - 1
        };
        CLPOLY_ASSERT_TRUE(gb_set_eq(gb, expected));
    }

    CLPOLY_TEST("cross_nonlinear2");
    {
        // Mathematica: {-1+y, 1+x} → monic: {y-1, x+1}
        polynomial_QQ f1 = QQ(2)*x + QQ(3)*y - 1;
        polynomial_QQ f2 = x - y + 2;
        auto gb = groebner_basis(std::vector<polynomial_QQ>{f1, f2});
        std::vector<polynomial_QQ> expected = {
            polynomial_QQ(y) - 1,
            polynomial_QQ(x) + 1
        };
        CLPOLY_ASSERT_TRUE(gb_set_eq(gb, expected));
    }

    CLPOLY_TEST("cross_katsura2");
    {
        // Mathematica: 4 elements
        // gb[1] = x + 2y + 2z - 1           (monic, LC=1)
        // gb[2] = -y - 4z + 10yz + 12z^2    (LT=10yz, monic: yz + 6/5*z^2 - 1/10*y - 2/5*z)
        // gb[3] = -y + 5y^2 + z - 3z^2      (LT=5y^2, monic: y^2 - 3/5*z^2 - 1/5*y + 1/5*z)
        // gb[4] = 7y + 3z - 79z^2 + 210z^3  (LT=210z^3, monic: z^3 - 79/210*z^2 + 1/30*y + 1/70*z)
        polynomial_QQ f1 = x + QQ(2)*y + QQ(2)*z - 1;
        polynomial_QQ f2 = pow(x,2) + QQ(2)*pow(y,2) + QQ(2)*pow(z,2) - x;
        polynomial_QQ f3 = QQ(2)*x*y + QQ(2)*y*z - y;
        auto gb = groebner_basis(std::vector<polynomial_QQ>{f1, f2, f3});
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)4);

        std::vector<polynomial_QQ> expected = {
            x + QQ(2)*y + QQ(2)*z - 1,
            y*z + QQ(6,5)*pow(z,2) - QQ(1,10)*y - QQ(2,5)*z,
            pow(y,2) - QQ(3,5)*pow(z,2) - QQ(1,5)*y + QQ(1,5)*z,
            pow(z,3) - QQ(79,210)*pow(z,2) + QQ(1,30)*y + QQ(1,70)*z
        };
        CLPOLY_ASSERT_TRUE(gb_set_eq(gb, expected));
    }

    // ======== 随机性质测试 ========

    // 验证 GB 性质的辅助函数
    auto verify_gb_properties = [](const std::vector<polynomial_QQ>& gb,
                                   const std::vector<polynomial_QQ>& gens,
                                   const std::string& label) -> bool {
        if (gb.empty()) {
            // 空 GB 意味着理想是零理想，所有生成元应为零
            for (auto& f : gens)
                if (!f.empty()) return false;
            return true;
        }

        // 性质 1: GB 中每个元素首一
        for (size_t i = 0; i < gb.size(); i++)
            if (gb[i].front().second != QQ(1)) return false;

        // 性质 2: 无冗余 (没有 LM(gi) | LM(gj) for i≠j)
        for (size_t i = 0; i < gb.size(); i++)
            for (size_t j = 0; j < gb.size(); j++) {
                if (i == j) continue;
                monomial q;
                if (is_divexact(q, gb[i].front().first, gb[j].front().first))
                    return false;
            }

        // 性质 3: 每个生成元 mod GB = 0
        for (auto& f : gens)
            if (!normal_form(f, gb).empty()) return false;

        // 性质 4: Buchberger 准则 - 所有 S-多项式 mod GB = 0
        for (size_t i = 0; i < gb.size(); i++)
            for (size_t j = i+1; j < gb.size(); j++) {
                auto sp = s_polynomial(gb[i], gb[j]);
                if (!normal_form(sp, gb).empty()) return false;
            }

        return true;
    };

    CLPOLY_TEST("random_gb_2var_qq");
    for (int round = 0; round < 20; round++) {
        CLPOLY_TEST_SECTION("round_" + std::to_string(round));
        auto f1 = random_polynomial_QQ({x, y}, 3, 3, {-10, 10}, 5);
        auto f2 = random_polynomial_QQ({x, y}, 3, 3, {-10, 10}, 5);
        if (f1.empty() && f2.empty()) continue;
        std::vector<polynomial_QQ> gens = {f1, f2};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_TRUE(verify_gb_properties(gb, gens, "2var_qq_" + std::to_string(round)));
    }

    CLPOLY_TEST("random_gb_3var_qq");
    for (int round = 0; round < 10; round++) {
        CLPOLY_TEST_SECTION("round_" + std::to_string(round));
        auto f1 = random_polynomial_QQ({x, y, z}, 2, 3, {-5, 5}, 3);
        auto f2 = random_polynomial_QQ({x, y, z}, 2, 3, {-5, 5}, 3);
        auto f3 = random_polynomial_QQ({x, y, z}, 2, 2, {-5, 5}, 3);
        if (f1.empty() && f2.empty() && f3.empty()) continue;
        std::vector<polynomial_QQ> gens = {f1, f2, f3};
        auto gb = groebner_basis(gens);
        CLPOLY_ASSERT_TRUE(verify_gb_properties(gb, gens, "3var_qq_" + std::to_string(round)));
    }

    CLPOLY_TEST("random_gb_2var_zz");
    for (int round = 0; round < 15; round++) {
        CLPOLY_TEST_SECTION("round_" + std::to_string(round));
        auto f1 = random_polynomial<ZZ>({x, y}, 3, 3, {-10, 10});
        auto f2 = random_polynomial<ZZ>({x, y}, 3, 3, {-10, 10});
        if (f1.empty() && f2.empty()) continue;
        std::vector<polynomial_ZZ> gens = {f1, f2};
        auto gb = groebner_basis(gens);
        // ZZ 结果: 本原, lc > 0
        for (auto& g : gb) {
            CLPOLY_ASSERT_TRUE(g.front().second > 0);
            ZZ c = abs(g.front().second);
            for (auto& term : g) c = gcd(c, abs(term.second));
            CLPOLY_ASSERT_EQ(c, ZZ(1));
        }
        // 转 QQ 验证 GB 性质
        std::vector<polynomial_QQ> gb_qq, gens_qq;
        for (auto& g : gb) { polynomial_QQ gq; poly_convert(g, gq); gb_qq.push_back(gq); }
        for (auto& f : gens) { polynomial_QQ fq; poly_convert(f, fq); gens_qq.push_back(fq); }
        // 生成元 mod GB = 0
        for (auto& fq : gens_qq)
            CLPOLY_ASSERT_TRUE(normal_form(fq, gb_qq).empty());
    }

    CLPOLY_TEST("random_gb_3var_zz");
    for (int round = 0; round < 5; round++) {
        CLPOLY_TEST_SECTION("round_" + std::to_string(round));
        auto f1 = random_polynomial<ZZ>({x, y, z}, 2, 3, {-5, 5});
        auto f2 = random_polynomial<ZZ>({x, y, z}, 2, 3, {-5, 5});
        if (f1.empty() && f2.empty()) continue;
        std::vector<polynomial_ZZ> gens = {f1, f2};
        auto gb = groebner_basis(gens);
        for (auto& g : gb) {
            CLPOLY_ASSERT_TRUE(g.front().second > 0);
        }
        std::vector<polynomial_QQ> gb_qq, gens_qq;
        for (auto& g : gb) { polynomial_QQ gq; poly_convert(g, gq); gb_qq.push_back(gq); }
        for (auto& f : gens) { polynomial_QQ fq; poly_convert(f, fq); gens_qq.push_back(fq); }
        for (auto& fq : gens_qq)
            CLPOLY_ASSERT_TRUE(normal_form(fq, gb_qq).empty());
    }

    // ======== lex 序交叉测试: 与 Mathematica 精确对比 ========

    // lex 序辅助: set 比较
    auto gb_set_eq_lex = [](const std::vector<polynomial_<QQ, lex>>& gb,
                            const std::vector<polynomial_<QQ, lex>>& expected) -> bool {
        if (gb.size() != expected.size()) return false;
        for (auto& e : expected) {
            bool found = false;
            for (auto& g : gb)
                if (g == e) { found = true; break; }
            if (!found) return false;
        }
        return true;
    };

    // lex 下构造多项式的辅助
    auto lex_var = [](variable v) -> polynomial_<QQ, lex> {
        polynomial_<QQ, lex> p;
        basic_monomial<lex> m;
        m.push_back({v, 1});
        p.data().push_back({m, QQ(1)});
        return p;
    };
    auto lex_pow = [](variable v, int64_t n) -> polynomial_<QQ, lex> {
        polynomial_<QQ, lex> p;
        basic_monomial<lex> m;
        m.push_back({v, n});
        p.data().push_back({m, QQ(1)});
        return p;
    };
    auto lex_const = [](QQ c) -> polynomial_<QQ, lex> {
        polynomial_<QQ, lex> p;
        basic_monomial<lex> m;  // empty = degree 0
        p.data().push_back({m, c});
        return p;
    };

    polynomial_<QQ, lex> lx = lex_var(x), ly = lex_var(y), lz = lex_var(z);

    CLPOLY_TEST("lex_circle_hyperbola");
    {
        // Mathematica lex: {1-y^2+y^4, x-y+y^3}
        //   → monic: {y^4-y^2+1, y^3-y+x}
        // Note: lex order x > y, so y^3 < x in lex but we need to check actual CLPoly ordering
        polynomial_<QQ, lex> f1 = lx*lx + ly*ly - lex_const(QQ(1));
        polynomial_<QQ, lex> f2 = lx*ly - lex_const(QQ(1));
        auto gb = groebner_basis(std::vector<polynomial_<QQ, lex>>{f1, f2});
        // lex GB should have 2 elements (elimination: last one is univariate in y)
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)2);

        // Expected: {y^4-y^2+1, x-y+y^3}
        polynomial_<QQ, lex> e1 = lex_pow(y,4) - lex_pow(y,2) + lex_const(QQ(1));
        polynomial_<QQ, lex> e2 = lx + lex_pow(y,3) - ly;
        CLPOLY_ASSERT_TRUE(gb_set_eq_lex(gb, {e1, e2}));
    }

    CLPOLY_TEST("lex_symmetric_3var");
    {
        // Mathematica lex: {-1+z^3, y^2+yz+z^2, x+y+z}
        //   → monic: {z^3-1, y^2+yz+z^2, x+y+z}
        polynomial_<QQ, lex> f1 = lx + ly + lz;
        polynomial_<QQ, lex> f2 = lx*ly + ly*lz + lz*lx;
        polynomial_<QQ, lex> f3 = lx*ly*lz - lex_const(QQ(1));
        auto gb = groebner_basis(std::vector<polynomial_<QQ, lex>>{f1, f2, f3});
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)3);

        polynomial_<QQ, lex> e1 = lex_pow(z,3) - lex_const(QQ(1));
        polynomial_<QQ, lex> e2 = lex_pow(y,2) + ly*lz + lex_pow(z,2);
        polynomial_<QQ, lex> e3 = lx + ly + lz;
        CLPOLY_ASSERT_TRUE(gb_set_eq_lex(gb, {e1, e2, e3}));
    }

    CLPOLY_TEST("lex_x2y_x3x");
    {
        // Mathematica lex: {-y+y^2, -x+xy, x^2-y}
        //   → monic: {y^2-y, xy-x, x^2-y}
        polynomial_<QQ, lex> f1 = lx*lx - ly;
        polynomial_<QQ, lex> f2 = lx*lx*lx - lx;
        auto gb = groebner_basis(std::vector<polynomial_<QQ, lex>>{f1, f2});
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)3);

        polynomial_<QQ, lex> e1 = lex_pow(y,2) - ly;
        polynomial_<QQ, lex> e2 = lx*ly - lx;
        polynomial_<QQ, lex> e3 = lx*lx - ly;
        CLPOLY_ASSERT_TRUE(gb_set_eq_lex(gb, {e1, e2, e3}));
    }

    CLPOLY_TEST("lex_katsura2");
    {
        // Mathematica lex: 3 elements
        // z + z^2 - 40z^3 + 84z^4  → monic: 84z^4 - 40z^3 + z^2 + z
        //   divide by 84: z^4 - 40/84*z^3 + 1/84*z^2 + 1/84*z
        //   = z^4 - 10/21*z^3 + 1/84*z^2 + 1/84*z
        // 7y + 3z - 79z^2 + 210z^3  → monic in leading term...
        //   lex order: y > z, so LT = 7y → monic: y + 3/7*z - 79/7*z^2 + 30*z^3
        // -7 + 7x + 8z + 158z^2 - 420z^3  → LT = 7x, monic: x + 8/7*z + 158/7*z^2 - 60*z^3 - 1
        polynomial_<QQ, lex> f1 = lx + QQ(2)*ly + QQ(2)*lz - lex_const(QQ(1));
        polynomial_<QQ, lex> f2 = lx*lx + QQ(2)*ly*ly + QQ(2)*lz*lz - lx;
        polynomial_<QQ, lex> f3 = QQ(2)*lx*ly + QQ(2)*ly*lz - ly;
        auto gb = groebner_basis(std::vector<polynomial_<QQ, lex>>{f1, f2, f3});
        CLPOLY_ASSERT_EQ(gb.size(), (size_t)3);

        // Verify elimination property: last element should be univariate in z
        // (In lex x > y > z, the GB triangulates)
        // Also verify generators reduce to 0
        auto nf1 = normal_form(f1, gb);
        auto nf2 = normal_form(f2, gb);
        auto nf3 = normal_form(f3, gb);
        CLPOLY_ASSERT_TRUE(nf1.empty());
        CLPOLY_ASSERT_TRUE(nf2.empty());
        CLPOLY_ASSERT_TRUE(nf3.empty());

        // Exact match with Mathematica (monic forms)
        polynomial_<QQ, lex> e1 = lx + QQ(8,7)*lz + QQ(158,7)*lex_pow(z,2)
                                  - QQ(60)*lex_pow(z,3) - lex_const(QQ(1));
        polynomial_<QQ, lex> e2 = ly + QQ(3,7)*lz - QQ(79,7)*lex_pow(z,2)
                                  + QQ(30)*lex_pow(z,3);
        polynomial_<QQ, lex> e3 = lex_pow(z,4) - QQ(10,21)*lex_pow(z,3)
                                  + QQ(1,84)*lex_pow(z,2) + QQ(1,84)*lz;
        CLPOLY_ASSERT_TRUE(gb_set_eq_lex(gb, {e1, e2, e3}));
    }

    return clpoly_test::test_summary();
}
