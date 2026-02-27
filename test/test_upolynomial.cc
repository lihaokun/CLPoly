#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

int main() {
    // ======== Construction ========
    CLPOLY_TEST("upoly_construction");
    {
        // 3x^2 + 2x + 1
        upolynomial_ZZ p({{umonomial(2), ZZ(3)}, {umonomial(1), ZZ(2)}, {umonomial(0), ZZ(1)}});
        CLPOLY_ASSERT_EQ((int)p.size(), 3);
        CLPOLY_ASSERT_EQ(degree(p), (int64_t)2);
    }

    // ======== Arithmetic ========
    CLPOLY_TEST("upoly_add");
    {
        upolynomial_ZZ p({{umonomial(2), ZZ(3)}, {umonomial(0), ZZ(1)}});   // 3x^2 + 1
        upolynomial_ZZ q({{umonomial(2), ZZ(-1)}, {umonomial(1), ZZ(2)}});  // -x^2 + 2x
        auto sum = p + q;
        // 2x^2 + 2x + 1
        CLPOLY_ASSERT_EQ(degree(sum), (int64_t)2);
        CLPOLY_ASSERT_EQ((int)sum.size(), 3);
    }

    CLPOLY_TEST("upoly_mul");
    {
        upolynomial_ZZ p({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});   // x + 1
        upolynomial_ZZ q({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(-1)}});  // x - 1
        auto product = p * q;
        // x^2 - 1
        upolynomial_ZZ expected({{umonomial(2), ZZ(1)}, {umonomial(0), ZZ(-1)}});
        CLPOLY_ASSERT_EQ(product, expected);
    }

    // ======== Evaluation ========
    CLPOLY_TEST("upoly_eval");
    {
        // 2x^3 - 3x + 1
        upolynomial_ZZ p({{umonomial(3), ZZ(2)}, {umonomial(1), ZZ(-3)}, {umonomial(0), ZZ(1)}});
        // Evaluate at x=2: 2*8 - 3*2 + 1 = 16 - 6 + 1 = 11
        ZZ result = assign(p, ZZ(2));
        CLPOLY_ASSERT_EQ(result, ZZ(11));
    }

    CLPOLY_TEST("upoly_eval_zero");
    {
        upolynomial_ZZ p({{umonomial(3), ZZ(2)}, {umonomial(1), ZZ(-3)}, {umonomial(0), ZZ(1)}});
        ZZ result = assign(p, ZZ(0));
        CLPOLY_ASSERT_EQ(result, ZZ(1));
    }

    // ======== Derivative ========
    CLPOLY_TEST("upoly_derivative");
    {
        // 3x^3 + 2x^2 - x + 5
        upolynomial_ZZ p({{umonomial(3), ZZ(3)}, {umonomial(2), ZZ(2)}, {umonomial(1), ZZ(-1)}, {umonomial(0), ZZ(5)}});
        auto dp = derivative(p);
        // 9x^2 + 4x - 1
        upolynomial_ZZ expected({{umonomial(2), ZZ(9)}, {umonomial(1), ZZ(4)}, {umonomial(0), ZZ(-1)}});
        CLPOLY_ASSERT_EQ(dp, expected);
    }

    // ======== Conversion: multivariate -> univariate ========
    CLPOLY_TEST("upoly_conversion");
    {
        variable x("x");
        polynomial_ZZ f = 2*pow(x,3) - 3*pow(x,2) + x - 5;
        upolynomial_ZZ uf;
        poly_convert(f, uf);
        CLPOLY_ASSERT_EQ(degree(uf), (int64_t)3);
        // Evaluate at x=1: 2 - 3 + 1 - 5 = -5
        ZZ val = assign(uf, ZZ(1));
        ZZ expected = ZZ(-5);
        CLPOLY_ASSERT_EQ(val, expected);
    }

    // ======== Power ========
    CLPOLY_TEST("upoly_pow");
    {
        upolynomial_ZZ p({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});  // x + 1
        auto p2 = pow(p, 2);
        auto p_times_p = p * p;
        CLPOLY_ASSERT_EQ(p2, p_times_p);
    }

    // ======== Content ========
    CLPOLY_TEST("upoly_content");
    {
        // 6x^2 + 4x + 2 -> content = 2
        upolynomial_ZZ p({{umonomial(2), ZZ(6)}, {umonomial(1), ZZ(4)}, {umonomial(0), ZZ(2)}});
        ZZ c = cont(p);
        CLPOLY_ASSERT_EQ(c, ZZ(2));
    }

    // ======== Subtraction ========
    CLPOLY_TEST("upoly_subtraction");
    {
        upolynomial_ZZ p({{umonomial(2), ZZ(3)}, {umonomial(0), ZZ(1)}});   // 3x^2 + 1
        upolynomial_ZZ q({{umonomial(2), ZZ(1)}, {umonomial(1), ZZ(2)}});   // x^2 + 2x
        auto diff = p - q;
        // 2x^2 - 2x + 1
        upolynomial_ZZ expected({{umonomial(2), ZZ(2)}, {umonomial(1), ZZ(-2)}, {umonomial(0), ZZ(1)}});
        CLPOLY_ASSERT_EQ(diff, expected);

        // p - p == 0
        auto zero_result = p - p;
        upolynomial_ZZ zero;
        CLPOLY_ASSERT_EQ(zero_result, zero);

        // p - 0 == p
        CLPOLY_ASSERT_EQ(p - zero, p);
    }

    // ======== Division ========
    CLPOLY_TEST("upoly_division");
    {
        // (x^2 - 1) / (x + 1) = x - 1
        upolynomial_ZZ p({{umonomial(2), ZZ(1)}, {umonomial(0), ZZ(-1)}});   // x^2 - 1
        upolynomial_ZZ q({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});    // x + 1
        auto quotient = p / q;
        upolynomial_ZZ expected({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(-1)}});  // x - 1
        CLPOLY_ASSERT_EQ(quotient, expected);

        // Verify: quotient * divisor == dividend
        CLPOLY_ASSERT_EQ(quotient * q, p);

        // Division of (x^3 + 1) by (x + 1) = x^2 - x + 1
        upolynomial_ZZ p2({{umonomial(3), ZZ(1)}, {umonomial(0), ZZ(1)}});
        auto q2 = p2 / q;
        CLPOLY_ASSERT_EQ(q2 * q, p2);
    }

    // ======== is_number ========
    CLPOLY_TEST("upoly_is_number");
    {
        // Constant polynomial
        upolynomial_ZZ c({{umonomial(0), ZZ(42)}});
        CLPOLY_ASSERT_TRUE(is_number(c));

        // Zero polynomial
        upolynomial_ZZ zero;
        CLPOLY_ASSERT_TRUE(is_number(zero));

        // Non-constant polynomial
        upolynomial_ZZ p({{umonomial(1), ZZ(1)}, {umonomial(0), ZZ(1)}});
        CLPOLY_ASSERT_FALSE(is_number(p));

        // Higher degree
        upolynomial_ZZ p2({{umonomial(3), ZZ(2)}, {umonomial(0), ZZ(1)}});
        CLPOLY_ASSERT_FALSE(is_number(p2));
    }

    // ======== Reverse conversion: univariate -> multivariate ========
    CLPOLY_TEST("upoly_reverse_convert");
    {
        variable x("x");
        // Create a multivariate polynomial with single variable
        polynomial_ZZ f = 3*pow(x,3) - 2*pow(x,2) + x - 5;
        // Convert to univariate
        upolynomial_ZZ uf;
        poly_convert(f, uf);
        CLPOLY_ASSERT_EQ(degree(uf), (int64_t)3);

        // Evaluate both at same point and compare
        ZZ val_u = assign(uf, ZZ(4));
        auto val_m = assign(f, x, ZZ(4));
        // 3*64 - 2*16 + 4 - 5 = 192 - 32 + 4 - 5 = 159
        CLPOLY_ASSERT_EQ(val_u, ZZ(159));
        CLPOLY_ASSERT_TRUE(is_number(val_m));
        CLPOLY_ASSERT_EQ(val_m.front().second, ZZ(159));

        // Convert between ZZ and QQ univariate
        upolynomial_<QQ> uq;
        poly_convert(uf, uq);
        CLPOLY_ASSERT_EQ(degree(uq), degree(uf));
        // Convert back
        upolynomial_ZZ uf2;
        poly_convert(uq, uf2);
        CLPOLY_ASSERT_EQ(uf2, uf);
    }

    // ======== poly_convert: Zp -> ZZ ========
    CLPOLY_TEST("poly_convert_zp_to_zz_basic");
    {
        uint64_t p = 97;
        // 3x^2 + 5x + 7（系数均小于 p，直接恢复）
        upolynomial_<Zp> zp_poly;
        zp_poly.push_back({umonomial(2), Zp(3, p)});
        zp_poly.push_back({umonomial(1), Zp(5, p)});
        zp_poly.push_back({umonomial(0), Zp(7, p)});

        upolynomial_ZZ zz_poly;
        poly_convert(zp_poly, zz_poly);

        CLPOLY_ASSERT_EQ(degree(zz_poly), (int64_t)2);
        CLPOLY_ASSERT_EQ((int)zz_poly.size(), 3);
        CLPOLY_ASSERT_EQ(zz_poly[0].second, ZZ(3));
        CLPOLY_ASSERT_EQ(zz_poly[1].second, ZZ(5));
        CLPOLY_ASSERT_EQ(zz_poly[2].second, ZZ(7));
    }

    CLPOLY_TEST("poly_convert_zp_to_zz_zero");
    {
        uint64_t p = 101;
        upolynomial_<Zp> zp_poly;  // 零多项式
        upolynomial_ZZ zz_poly;
        poly_convert(zp_poly, zz_poly);
        CLPOLY_ASSERT_TRUE(zz_poly.empty());
    }

    CLPOLY_TEST("poly_convert_zp_to_zz_constant");
    {
        uint64_t p = 101;
        upolynomial_<Zp> zp_poly;
        zp_poly.push_back({umonomial(0), Zp(42, p)});
        upolynomial_ZZ zz_poly;
        poly_convert(zp_poly, zz_poly);
        CLPOLY_ASSERT_EQ((int)zz_poly.size(), 1);
        CLPOLY_ASSERT_EQ(zz_poly[0].second, ZZ(42));
    }

    CLPOLY_TEST("poly_convert_zp_to_zz_roundtrip");
    {
        // ZZ -> Zp (polynomial_mod) -> ZZ (poly_convert)，系数在 [0, p-1] 内时为恒等
        uint64_t p = 97;
        upolynomial_ZZ orig({{umonomial(3), ZZ(10)}, {umonomial(1), ZZ(3)}, {umonomial(0), ZZ(1)}});
        auto zp = polynomial_mod(orig, p);
        upolynomial_ZZ recovered;
        poly_convert(zp, recovered);
        CLPOLY_ASSERT_EQ(orig, recovered);
    }

    CLPOLY_TEST("poly_convert_zp_to_zz_zero_coeff_dropped");
    {
        // 系数恰好为 p 的倍数，polynomial_mod 应丢弃该项
        uint64_t p = 7;
        // 7x^2 + 3x + 1 → mod 7 → 3x + 1
        upolynomial_ZZ orig({{umonomial(2), ZZ(7)}, {umonomial(1), ZZ(3)}, {umonomial(0), ZZ(1)}});
        auto zp = polynomial_mod(orig, p);
        upolynomial_ZZ recovered;
        poly_convert(zp, recovered);
        upolynomial_ZZ expected({{umonomial(1), ZZ(3)}, {umonomial(0), ZZ(1)}});
        CLPOLY_ASSERT_EQ(recovered, expected);
    }

    // ======== Scalar multiplication ========
    CLPOLY_TEST("upoly_scalar_mul");
    {
        upolynomial_ZZ p({{umonomial(2), ZZ(3)}, {umonomial(1), ZZ(-2)}, {umonomial(0), ZZ(1)}});
        // p * 5
        auto r = p * ZZ(5);
        upolynomial_ZZ expected({{umonomial(2), ZZ(15)}, {umonomial(1), ZZ(-10)}, {umonomial(0), ZZ(5)}});
        CLPOLY_ASSERT_EQ(r, expected);
        // 5 * p (commutative)
        CLPOLY_ASSERT_EQ(ZZ(5) * p, expected);
        // p * 0 == 0
        CLPOLY_ASSERT_EQ(p * ZZ(0), upolynomial_ZZ());
        // p * 1 == p
        CLPOLY_ASSERT_EQ(p * ZZ(1), p);
        // p * (-1) == -p
        CLPOLY_ASSERT_EQ(p * ZZ(-1), -p);
    }

    return clpoly_test::test_summary();
}
