#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

int main() {
    variable x("x"), y("y"), z("z");

    // ======== degree ========
    CLPOLY_TEST("degree_total");
    {
        polynomial_ZZ f = 3*pow(x,3)*pow(y,2) + x*y - 5;
        CLPOLY_ASSERT_EQ(degree(f), (int64_t)5);  // x^3*y^2 has total degree 5
    }

    CLPOLY_TEST("degree_in_variable");
    {
        polynomial_ZZ f = 3*pow(x,3)*pow(y,2) + x*y - 5;
        CLPOLY_ASSERT_EQ(degree(f, x), (int64_t)3);
        CLPOLY_ASSERT_EQ(degree(f, y), (int64_t)2);
        CLPOLY_ASSERT_EQ(degree(f, z), (int64_t)0);  // z not in f
    }

    CLPOLY_TEST("degree_zero_poly");
    {
        polynomial_ZZ zero;
        CLPOLY_ASSERT_EQ(degree(zero), (int64_t)0);
        CLPOLY_ASSERT_EQ(degree(zero, x), (int64_t)0);
    }

    CLPOLY_TEST("degree_constant");
    {
        polynomial_ZZ c({{monomial(), ZZ(7)}});
        CLPOLY_ASSERT_EQ(degree(c), (int64_t)0);
        CLPOLY_ASSERT_EQ(degree(c, x), (int64_t)0);
    }

    // ======== leadcoeff ========
    CLPOLY_TEST("leadcoeff_univariate");
    {
        polynomial_ZZ f = 5*pow(x,3) - 2*pow(x,2) + x + 7;
        auto lc = leadcoeff(f, x);
        // leadcoeff(5x^3 - 2x^2 + x + 7, x) = 5
        CLPOLY_ASSERT_TRUE(is_number(lc));
        CLPOLY_ASSERT_EQ(lc.front().second, ZZ(5));
    }

    CLPOLY_TEST("leadcoeff_multivariate");
    {
        polynomial_ZZ f = pow(x,2)*y + 3*pow(x,2) + x*pow(y,2);
        auto lc = leadcoeff(f, x);
        // leadcoeff(x^2*y + 3x^2 + x*y^2, x) = y + 3
        polynomial_ZZ expected = y + 3;
        CLPOLY_ASSERT_EQ(lc, expected);
    }

    // ======== Regression: leadcoeff 通用版本 monomial 规范化 ========
    CLPOLY_TEST("leadcoeff_normalization");
    {
        // Bug: leadcoeff(O,F,v) 通用版直接 m[...].second=0
        //      留下零指数项 (v,0) 且 __deg 未更新
        //      导致结果 monomial 不规范，与正常构造的多项式不等

        // 3 变量，提取 x 的 leading coefficient
        polynomial_ZZ f = 3*pow(x,2)*pow(y,2) + pow(x,2)*z + 5*pow(x,2) + 2*x*y + 7;
        auto lc = leadcoeff(f, x);
        // lc(f, x) = 3*y^2 + z + 5 （x^2 的系数）
        polynomial_ZZ expected = 3*pow(y,2) + z + 5;
        CLPOLY_ASSERT_EQ(lc, expected);

        // 结果不应包含变量 x
        auto vars = get_variables(lc);
        for (auto& vp : vars) {
            CLPOLY_ASSERT_NE(vp.first, x);
        }

        // 结果中各 monomial 的 deg() 应与手工构造的一致
        // 若有零指数残留或 __deg 未更新，排序和 deg 都会错
        auto it_lc = lc.begin();
        auto it_exp = expected.begin();
        for (; it_lc != lc.end() && it_exp != expected.end(); ++it_lc, ++it_exp) {
            CLPOLY_ASSERT_EQ(it_lc->first.deg(), it_exp->first.deg());
            CLPOLY_ASSERT_EQ(it_lc->first, it_exp->first);
        }
    }

    CLPOLY_TEST("leadcoeff_result_usable_in_arithmetic");
    {
        // leadcoeff 的结果应能正常参与后续运算
        polynomial_ZZ f = pow(x,3)*pow(y,2)*z + 2*pow(x,3)*y + pow(x,2) - 1;
        auto lc = leadcoeff(f, x);  // y^2*z + 2*y
        polynomial_ZZ g = y + 1;

        // 如果 lc 内部 monomial 不规范，乘法结果会出错
        auto product = lc * g;
        polynomial_ZZ expected_lc = pow(y,2)*z + 2*y;
        auto expected_product = expected_lc * g;
        CLPOLY_ASSERT_EQ(product, expected_product);
    }

    // ======== is_number ========
    CLPOLY_TEST("is_number");
    {
        polynomial_ZZ zero;
        polynomial_ZZ c({{monomial(), ZZ(42)}});
        polynomial_ZZ f = pow(x,2) + 1;

        CLPOLY_ASSERT_TRUE(is_number(zero));
        CLPOLY_ASSERT_TRUE(is_number(c));
        CLPOLY_ASSERT_FALSE(is_number(f));
    }

    // ======== get_variables ========
    CLPOLY_TEST("get_variables");
    {
        polynomial_ZZ f = pow(x,2)*y + z + 1;
        auto vars = get_variables(f);
        CLPOLY_ASSERT_EQ((int)vars.size(), 3);
    }

    CLPOLY_TEST("get_variables_constant");
    {
        polynomial_ZZ c({{monomial(), ZZ(5)}});
        auto vars = get_variables(c);
        CLPOLY_ASSERT_EQ((int)vars.size(), 0);
    }

    // ======== assign ========
    CLPOLY_TEST("assign_substitute");
    {
        polynomial_ZZ f = pow(x,2) + 2*x + 1;
        auto result = assign(f, x, ZZ(3));
        // x=3: 9 + 6 + 1 = 16
        CLPOLY_ASSERT_TRUE(is_number(result));
        CLPOLY_ASSERT_EQ(result.front().second, ZZ(16));
    }

    CLPOLY_TEST("assign_partial");
    {
        polynomial_ZZ f = pow(x,2)*y + x*z + 1;
        auto result = assign(f, y, ZZ(2));
        // 2*x^2 + x*z + 1
        polynomial_ZZ expected = 2*pow(x,2) + x*z + 1;
        CLPOLY_ASSERT_EQ(result, expected);
    }

    CLPOLY_TEST("assign_variable_not_present");
    {
        polynomial_ZZ f = pow(x,2) + 1;
        auto result = assign(f, y, ZZ(5));
        CLPOLY_ASSERT_EQ(result, f);
    }

    // ======== derivative ========
    CLPOLY_TEST("derivative_univariate");
    {
        polynomial_ZZ f = 3*pow(x,4) - 2*pow(x,2) + x - 5;
        auto df = derivative(f, x);
        polynomial_ZZ expected = 12*pow(x,3) - 4*x + 1;
        CLPOLY_ASSERT_EQ(df, expected);
    }

    CLPOLY_TEST("derivative_multivariate");
    {
        polynomial_ZZ f = pow(x,3)*pow(y,2) + x*y + pow(y,3);
        auto dfx = derivative(f, x);
        polynomial_ZZ expected_x = 3*pow(x,2)*pow(y,2) + y;
        CLPOLY_ASSERT_EQ(dfx, expected_x);

        auto dfy = derivative(f, y);
        polynomial_ZZ expected_y = 2*pow(x,3)*y + x + 3*pow(y,2);
        CLPOLY_ASSERT_EQ(dfy, expected_y);
    }

    CLPOLY_TEST("derivative_constant_is_zero");
    {
        polynomial_ZZ c({{monomial(), ZZ(42)}});
        auto dc = derivative(c, x);
        polynomial_ZZ zero;
        CLPOLY_ASSERT_EQ(dc, zero);
    }

    CLPOLY_TEST("derivative_zero_is_zero");
    {
        polynomial_ZZ zero;
        auto dz = derivative(zero, x);
        CLPOLY_ASSERT_EQ(dz, zero);
    }

    // ======== coeff ========
    CLPOLY_TEST("coeff_extraction");
    {
        polynomial_ZZ f = 3*pow(x,2) + 2*x + 1;
        auto coeffs = coeff(f, x);
        // coeffs[0] is coefficient of x^0, coeffs[1] of x^1, coeffs[2] of x^2
        CLPOLY_ASSERT_EQ((int)coeffs.size(), 3);
    }

    // ======== Equality and comparison ========
    CLPOLY_TEST("equality");
    {
        polynomial_ZZ f = pow(x,2) + 2*x + 1;
        polynomial_ZZ g = pow(x,2) + 2*x + 1;
        polynomial_ZZ h = pow(x,2) + 2*x;
        CLPOLY_ASSERT_EQ(f, g);
        CLPOLY_ASSERT_NE(f, h);
    }

    CLPOLY_TEST("copy_construct");
    {
        polynomial_ZZ f = pow(x,3) - x + 1;
        polynomial_ZZ g(f);
        CLPOLY_ASSERT_EQ(f, g);
    }

    // ======== poly_convert: ordering ========
    CLPOLY_TEST("poly_convert_ordering");
    {
        polynomial_ZZ f = pow(x,2)*y + x*pow(z,2) + y*z + 1;
        lex_<custom_var_order> mo(custom_var_order({z, y, x}));
        polynomial_<ZZ, lex_<custom_var_order>> f_lex(&mo);
        poly_convert(f, f_lex);
        CLPOLY_ASSERT_EQ((int)f_lex.size(), (int)f.size());
        // Convert back and check equality
        polynomial_ZZ f_back;
        poly_convert(f_lex, f_back);
        CLPOLY_ASSERT_EQ(f_back, f);
    }

    // ======== polynomial_QQ ========
    CLPOLY_TEST("polynomial_QQ_division");
    {
        polynomial_QQ p1(1);
        polynomial_QQ p2(2);
        auto result = p1 / p2;
        CLPOLY_ASSERT_TRUE(is_number(result));
        CLPOLY_ASSERT_FALSE(result.empty());
    }

    // ======== Hash ========
    CLPOLY_TEST("polynomial_hash");
    {
        polynomial_ZZ f = pow(x,2) + 2*x + 1;
        polynomial_ZZ g = pow(x,2) + 2*x + 1;
        auto hf = std::hash<polynomial_ZZ>()(f);
        auto hg = std::hash<polynomial_ZZ>()(g);
        CLPOLY_ASSERT_EQ(hf, hg);
    }

    // ======== get_variables for vector of polynomials ========
    CLPOLY_TEST("get_variables_vector");
    {
        polynomial_ZZ f1 = pow(x,2) + y;
        polynomial_ZZ f2 = pow(z,3) + 1;
        polynomial_ZZ f3 = x*y;
        std::vector<polynomial_ZZ> polys = {f1, f2, f3};
        auto vars = get_variables(polys);
        CLPOLY_ASSERT_EQ((int)vars.size(), 3);  // x, y, z

        // Single-polynomial vector
        std::vector<polynomial_ZZ> single = {f1};
        auto vars1 = get_variables(single);
        CLPOLY_ASSERT_EQ((int)vars1.size(), 2);  // x, y

        // Empty vector
        std::vector<polynomial_ZZ> empty_vec;
        auto vars_empty = get_variables(empty_vec);
        CLPOLY_ASSERT_EQ((int)vars_empty.size(), 0);
    }

    // ======== assign with map (multi-variable substitution) ========
    CLPOLY_TEST("assign_multi");
    {
        polynomial_ZZ f = pow(x,2)*y + x*z + y*z + 1;
        std::map<variable, ZZ> subs;
        subs[x] = ZZ(2);
        subs[y] = ZZ(3);
        auto result = assign(f, subs);
        // x=2, y=3: 4*3 + 2*z + 3*z + 1 = 12 + 5*z + 1 = 5*z + 13
        polynomial_ZZ expected = 5*z + 13;
        CLPOLY_ASSERT_EQ(result, expected);

        // Substitute all variables
        subs[z] = ZZ(-1);
        auto full_result = assign(f, subs);
        // x=2, y=3, z=-1: 12 - 2 - 3 + 1 = 8
        CLPOLY_ASSERT_TRUE(is_number(full_result));
        CLPOLY_ASSERT_EQ(full_result.front().second, ZZ(8));

        // Empty map: no substitution
        std::map<variable, ZZ> empty_map;
        auto unchanged = assign(f, empty_map);
        CLPOLY_ASSERT_EQ(unchanged, f);
    }

    // ======== coeff: univariate coefficient extraction ========
    CLPOLY_TEST("coeff_univariate");
    {
        // f = 5x^3 - 2x^2 + 3x + 7
        polynomial_ZZ f = 5*pow(x,3) - 2*pow(x,2) + 3*x + 7;
        auto coeffs = coeff(f, x);
        // coeffs should have 4 entries: coeff of x^0 through x^3
        CLPOLY_ASSERT_EQ((int)coeffs.size(), 4);

        // Verify each coefficient
        CLPOLY_ASSERT_EQ(coeffs[0].front().second, ZZ(7));   // x^0
        CLPOLY_ASSERT_EQ(coeffs[1].front().second, ZZ(3));   // x^1
        CLPOLY_ASSERT_EQ(coeffs[2].front().second, ZZ(-2));  // x^2
        CLPOLY_ASSERT_EQ(coeffs[3].front().second, ZZ(5));   // x^3
    }

    return clpoly_test::test_summary();
}
