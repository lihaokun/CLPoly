/**
 * @file test_coeff.cc
 * @brief 测试各种序下的系数计算正确性（统一使用多项式 5*z*x*y + 2*x*y + 3*z*y + 4*y + 8）
 */
#include "clpoly_test.hh"
#include <clpoly/clpoly.hh>
#include <iostream>

int main() {
    using namespace clpoly;

    // 定义变量
    variable x("x"), y("y"), z("z");

    // ---------- 1. 分次字典序 (grlex) 测试 ----------
    // polynomial_ZZ 默认使用分次字典序
    polynomial_ZZ poly_grlex = 5*z*x*y + 2*x*y + 3*z*y + 4*y + 8;

    CLPOLY_TEST("coeff in grlex order with respect to z");
    {
        auto coeffs = coeff(poly_grlex, z);
        CLPOLY_ASSERT_EQ(coeffs.size(), 2);

        polynomial_ZZ expected0 = 2*x*y + 4*y + 8;
        CLPOLY_ASSERT_EQ(coeffs[0], expected0);

        polynomial_ZZ expected1 = 5*x*y + 3*y;
        CLPOLY_ASSERT_EQ(coeffs[1], expected1);
    }

    CLPOLY_TEST("coeff in grlex order with respect to y");
    {
        auto coeffs = coeff(poly_grlex, y);
        CLPOLY_ASSERT_EQ(coeffs.size(), 2);

        polynomial_ZZ expected0 = 8;
        CLPOLY_ASSERT_EQ(coeffs[0], expected0);

        polynomial_ZZ expected1 = 5*z*x + 2*x + 3*z + 4;
        CLPOLY_ASSERT_EQ(coeffs[1], expected1);
    }

    CLPOLY_TEST("coeff in grlex order with respect to x");
    {
        auto coeffs = coeff(poly_grlex, x);
        CLPOLY_ASSERT_EQ(coeffs.size(), 2);

        polynomial_ZZ expected0 = 3*z*y + 4*y + 8;
        CLPOLY_ASSERT_EQ(coeffs[0], expected0);

        polynomial_ZZ expected1 = 5*z*y + 2*y;
        CLPOLY_ASSERT_EQ(coeffs[1], expected1);
    }

    // ---------- 2. 字典序 (lex) 测试 ----------
    // 自定义变量序: x < y < z
    lex_<custom_var_order> mo_lex;
    mo_lex = lex_<custom_var_order>(custom_var_order({x, y, z}));
    polynomial_<ZZ, lex_<custom_var_order>> poly_lex(&mo_lex);
    poly_lex = 5*z*x*y + 2*x*y + 3*z*y + 4*y + 8;

    CLPOLY_TEST("coeff in lex order with respect to z");
    {
        auto coeffs = coeff(poly_lex, z);
        CLPOLY_ASSERT_EQ(coeffs.size(), 2);

        polynomial_<ZZ, lex_<custom_var_order>> expected0(&mo_lex);
        expected0 = 2*x*y + 4*y + 8;
        CLPOLY_ASSERT_EQ(coeffs[0], expected0);

        polynomial_<ZZ, lex_<custom_var_order>> expected1(&mo_lex);
        expected1 = 5*x*y + 3*y;
        CLPOLY_ASSERT_EQ(coeffs[1], expected1);
    }

    CLPOLY_TEST("coeff in lex order with respect to y");
    {
        auto coeffs = coeff(poly_lex, y);
        CLPOLY_ASSERT_EQ(coeffs.size(), 2);

        polynomial_<ZZ, lex_<custom_var_order>> expected0(&mo_lex);
        expected0 = 8;
        CLPOLY_ASSERT_EQ(coeffs[0], expected0);

        polynomial_<ZZ, lex_<custom_var_order>> expected1(&mo_lex);
        expected1 = 5*z*x + 2*x + 3*z + 4;
        CLPOLY_ASSERT_EQ(coeffs[1], expected1);
    }

    CLPOLY_TEST("coeff in lex order with respect to x");
    {
        auto coeffs = coeff(poly_lex, x);
        CLPOLY_ASSERT_EQ(coeffs.size(), 2);

        polynomial_<ZZ, lex_<custom_var_order>> expected0(&mo_lex);
        expected0 = 3*z*y + 4*y + 8;
        CLPOLY_ASSERT_EQ(coeffs[0], expected0);

        polynomial_<ZZ, lex_<custom_var_order>> expected1(&mo_lex);
        expected1 = 5*z*y + 2*y;
        CLPOLY_ASSERT_EQ(coeffs[1], expected1);
    }

    // 字典序下测试无参 coeff(poly) 默认使用首变量 x
    CLPOLY_TEST("coeff in lex order without variable (default main variable x)");
    {
        auto coeffs_default = coeff(poly_lex);          // 无参版本，应默认按首变量 x 展开
        auto coeffs_x = coeff(poly_lex, x);             // 显式指定 x
        CLPOLY_ASSERT_EQ(coeffs_default.size(), coeffs_x.size());
        for (size_t i = 0; i < coeffs_default.size(); ++i) {
            CLPOLY_ASSERT_EQ(coeffs_default[i], coeffs_x[i]);
        }
    }

    // ---------- 3. 单变量优先序 (univariate priority) 测试 ----------
    // 以 x 为优先变量（假设优先变量是 x）
    univariate_priority_order uni_order_x(x);
    polynomial_<ZZ, lex_<univariate_order>> poly_uni_x(&uni_order_x);
    poly_uni_x = 5*z*x*y + 2*x*y + 3*z*y + 4*y + 8;

    CLPOLY_TEST("coeff in univariate priority (x first) with respect to z");
    {
        auto coeffs = coeff(poly_uni_x, z);
        CLPOLY_ASSERT_EQ(coeffs.size(), 2);

        polynomial_<ZZ, lex_<univariate_order>> expected0(&uni_order_x);
        expected0 = 2*x*y + 4*y + 8;
        CLPOLY_ASSERT_EQ(coeffs[0], expected0);

        polynomial_<ZZ, lex_<univariate_order>> expected1(&uni_order_x);
        expected1 = 5*x*y + 3*y;
        CLPOLY_ASSERT_EQ(coeffs[1], expected1);
    }

    CLPOLY_TEST("coeff in univariate priority (x first) with respect to x");
    {
        auto coeffs = coeff(poly_uni_x, x);
        CLPOLY_ASSERT_EQ(coeffs.size(), 2);

        polynomial_<ZZ, lex_<univariate_order>> expected0(&uni_order_x);
        expected0 = 3*z*y + 4*y + 8;
        CLPOLY_ASSERT_EQ(coeffs[0], expected0);

        polynomial_<ZZ, lex_<univariate_order>> expected1(&uni_order_x);
        expected1 = 5*z*y + 2*y;
        CLPOLY_ASSERT_EQ(coeffs[1], expected1);
    }

    CLPOLY_TEST("coeff in univariate priority (x first) with respect to y");
    {
        auto coeffs = coeff(poly_uni_x, y);
        CLPOLY_ASSERT_EQ(coeffs.size(), 2);

        polynomial_<ZZ, lex_<univariate_order>> expected0(&uni_order_x);
        expected0 = 8;
        CLPOLY_ASSERT_EQ(coeffs[0], expected0);

        polynomial_<ZZ, lex_<univariate_order>> expected1(&uni_order_x);
        expected1 = 5*z*x + 2*x + 3*z + 4;
        CLPOLY_ASSERT_EQ(coeffs[1], expected1);
    }

    // 单变量优先序下测试无参 coeff(poly) 默认使用首变量 x
    CLPOLY_TEST("coeff in univariate priority without variable (default main variable x)");
    {
        auto coeffs_default = coeff(poly_uni_x);        // 无参版本，应默认按首变量 x 展开
        auto coeffs_x = coeff(poly_uni_x, x);           // 显式指定 x
        CLPOLY_ASSERT_EQ(coeffs_default.size(), coeffs_x.size());
        for (size_t i = 0; i < coeffs_default.size(); ++i) {
            CLPOLY_ASSERT_EQ(coeffs_default[i], coeffs_x[i]);
        }
    }

        // ---------- 4. 常数多项式测试 ----------
    polynomial_ZZ const_poly_grlex = 10;   // 分次字典序

    // 字典序常数多项式
    polynomial_<ZZ, lex_<custom_var_order>> const_poly_lex(&mo_lex);
    const_poly_lex = 10;

    // 单变量优先序常数多项式
    polynomial_<ZZ, lex_<univariate_order>> const_poly_uni(&uni_order_x);
    const_poly_uni = 10;

    CLPOLY_TEST("coeff of constant polynomial with respect to any variable");
    {
        // 测试分次字典序（polynomial_ZZ）
        auto coeffs_z = coeff(const_poly_grlex, z);
        CLPOLY_ASSERT_EQ(coeffs_z.size(), 1);
        polynomial_ZZ expected = 10;
        CLPOLY_ASSERT_EQ(coeffs_z[0], expected);

        auto coeffs_y = coeff(const_poly_grlex, y);
        CLPOLY_ASSERT_EQ(coeffs_y.size(), 1);
        CLPOLY_ASSERT_EQ(coeffs_y[0], expected);

            
        auto coeffs_x = coeff(const_poly_grlex, x);
        CLPOLY_ASSERT_EQ(coeffs_x.size(), 1);
        CLPOLY_ASSERT_EQ(coeffs_x[0], expected);

        // 测试字典序无参 coeff (不指定变量)
        auto coeffs_lex = coeff(const_poly_lex);
        CLPOLY_ASSERT_EQ(coeffs_lex.size(), 1);
        CLPOLY_ASSERT_EQ(coeffs_lex[0], const_poly_lex);

        // 测试单变量优先序 coeff (不指定变量)
        auto coeffs_uni = coeff(const_poly_uni);
        CLPOLY_ASSERT_EQ(coeffs_uni.size(), 1);
        CLPOLY_ASSERT_EQ(coeffs_uni[0], const_poly_uni);
    }

    return clpoly_test::test_summary();
}