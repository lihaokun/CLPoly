#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include "testdata_expected.hh"

using namespace clpoly;

int main() {
    variable x("x"), y("y"), z("z");

    // ======== Simple: single variable ========
    CLPOLY_TEST("arith_simple_add");
    {
        auto f = testdata::arith_f1_simple();
        auto g = testdata::arith_g1_simple();
        CLPOLY_ASSERT_EQ(f + g, testdata::arith_add_simple());
    }

    CLPOLY_TEST("arith_simple_sub");
    {
        auto f = testdata::arith_f1_simple();
        auto g = testdata::arith_g1_simple();
        CLPOLY_ASSERT_EQ(f - g, testdata::arith_sub_simple());
    }

    CLPOLY_TEST("arith_simple_mul");
    {
        auto f = testdata::arith_f1_simple();
        auto g = testdata::arith_g1_simple();
        CLPOLY_ASSERT_EQ(f * g, testdata::arith_mul_simple());
    }

    CLPOLY_TEST("arith_simple_pow2");
    {
        auto f = testdata::arith_f1_simple();
        CLPOLY_ASSERT_EQ(pow(f, 2), testdata::arith_pow2_simple());
        CLPOLY_ASSERT_EQ(pow(f, 2), f * f);
    }

    CLPOLY_TEST("arith_simple_pow3");
    {
        auto f = testdata::arith_f1_simple();
        CLPOLY_ASSERT_EQ(pow(f, 3), testdata::arith_pow3_simple());
        CLPOLY_ASSERT_EQ(pow(f, 3), f * f * f);
    }

    // ======== Medium: two variables ========
    CLPOLY_TEST("arith_medium_mul");
    {
        auto f = testdata::arith_f1_medium();
        auto g = testdata::arith_g1_medium();
        CLPOLY_ASSERT_EQ(f * g, testdata::arith_mul_medium());
    }

    // ======== Complex: three variables ========
    CLPOLY_TEST("arith_complex_mul");
    {
        auto f = testdata::arith_f1_complex();
        auto g = testdata::arith_g1_complex();
        CLPOLY_ASSERT_EQ(f * g, testdata::arith_mul_complex());
    }

    CLPOLY_TEST("arith_complex_random_mul");
    {
        auto f = testdata::complex_f1();
        auto g = testdata::complex_g1();
        CLPOLY_ASSERT_EQ(f * g, testdata::complex_mul1());
    }

    // ======== Algebraic identities ========
    CLPOLY_TEST("arith_commutativity");
    {
        auto f = testdata::arith_f1_medium();
        auto g = testdata::arith_g1_medium();
        CLPOLY_ASSERT_EQ(f + g, g + f);
        CLPOLY_ASSERT_EQ(f * g, g * f);
    }

    CLPOLY_TEST("arith_associativity");
    {
        polynomial_ZZ f = 3*pow(x,2) + x - 1;
        polynomial_ZZ g = pow(x,2)*y + 2;
        polynomial_ZZ h = y*pow(x,3) - x;
        CLPOLY_ASSERT_EQ((f + g) + h, f + (g + h));
        CLPOLY_ASSERT_EQ((f * g) * h, f * (g * h));
    }

    CLPOLY_TEST("arith_distributivity");
    {
        auto f = testdata::arith_f1_complex();
        auto g = testdata::arith_g1_complex();
        polynomial_ZZ h = pow(x,2)*y - z + 3;
        CLPOLY_ASSERT_EQ(f * (g + h), f * g + f * h);
    }

    // ======== Zero and identity ========
    CLPOLY_TEST("arith_zero_identity");
    {
        auto f = testdata::arith_f1_simple();
        polynomial_ZZ zero;
        CLPOLY_ASSERT_EQ(f + zero, f);
        CLPOLY_ASSERT_EQ(f - f, zero);
        CLPOLY_ASSERT_EQ(f * zero, zero);
    }

    CLPOLY_TEST("arith_one_identity");
    {
        auto f = testdata::arith_f1_medium();
        polynomial_ZZ one = polynomial_ZZ({{monomial(), ZZ(1)}});
        CLPOLY_ASSERT_EQ(f * one, f);
        CLPOLY_ASSERT_EQ(pow(f, 0), one);
        CLPOLY_ASSERT_EQ(pow(f, 1), f);
    }

    CLPOLY_TEST("arith_negation");
    {
        auto f = testdata::arith_f1_simple();
        polynomial_ZZ zero;
        CLPOLY_ASSERT_EQ(f + (-f), zero);
        CLPOLY_ASSERT_EQ(-(-f), f);
    }

    // ======== Scalar operations ========
    CLPOLY_TEST("arith_scalar_mul");
    {
        polynomial_ZZ f = pow(x,3) - 2*x + 1;
        polynomial_ZZ expected = 3*pow(x,3) - 6*x + 3;
        CLPOLY_ASSERT_EQ(f * ZZ(3), expected);
    }

    CLPOLY_TEST("arith_scalar_add");
    {
        polynomial_ZZ f = pow(x,2) + x;
        polynomial_ZZ expected = pow(x,2) + x + 5;
        CLPOLY_ASSERT_EQ(f + 5, expected);
    }

    // ======== Constant polynomial ========
    CLPOLY_TEST("arith_constant");
    {
        polynomial_ZZ c({{monomial(), ZZ(7)}});
        CLPOLY_ASSERT_EQ(pow(c, 3), polynomial_ZZ({{monomial(), ZZ(343)}}));
    }

    // ======== pow identity: pow(f,2) == f*f, pow(f,3) == f*f*f ========
    CLPOLY_TEST("arith_pow_identity");
    {
        auto f = testdata::arith_f1_complex();
        CLPOLY_ASSERT_EQ(pow(f, 2), f * f);
        CLPOLY_ASSERT_EQ(pow(f, 3), f * f * f);
    }

    // ======== Scalar reverse operations: scalar op polynomial ========
    CLPOLY_TEST("polynomial_scalar_reverse");
    {
        polynomial_ZZ f = pow(x,2) + x + 1;

        // scalar - polynomial
        auto diff = 10 - f;
        polynomial_ZZ expected_diff = -pow(x,2) - x + 9;
        CLPOLY_ASSERT_EQ(diff, expected_diff);

        // scalar * polynomial
        auto prod = ZZ(3) * f;
        polynomial_ZZ expected_prod = 3*pow(x,2) + 3*x + 3;
        CLPOLY_ASSERT_EQ(prod, expected_prod);

        // scalar + polynomial
        auto sum = 5 + f;
        polynomial_ZZ expected_sum = pow(x,2) + x + 6;
        CLPOLY_ASSERT_EQ(sum, expected_sum);

        // Commutativity: scalar*f == f*scalar
        CLPOLY_ASSERT_EQ(ZZ(3) * f, f * ZZ(3));
        CLPOLY_ASSERT_EQ(5 + f, f + 5);
    }

    // ======== Unary operators ========
    CLPOLY_TEST("polynomial_unary");
    {
        polynomial_ZZ f = pow(x,3) - 2*x + 1;

        // Unary +
        polynomial_ZZ pos_f = +f;
        CLPOLY_ASSERT_EQ(pos_f, f);

        // Unary -
        polynomial_ZZ neg_f = -f;
        polynomial_ZZ expected_neg = -pow(x,3) + 2*x - 1;
        CLPOLY_ASSERT_EQ(neg_f, expected_neg);

        // +(-f) == -f
        CLPOLY_ASSERT_EQ(+neg_f, neg_f);

        // Unary on zero
        polynomial_ZZ zero;
        CLPOLY_ASSERT_EQ(+zero, zero);
        CLPOLY_ASSERT_EQ(-zero, zero);

        // Unary on constant
        polynomial_ZZ c({{monomial(), ZZ(5)}});
        polynomial_ZZ neg_c({{monomial(), ZZ(-5)}});
        CLPOLY_ASSERT_EQ(-c, neg_c);
    }

    // ======== pow(variable, n) to construct polynomial ========
    CLPOLY_TEST("polynomial_variable_pow");
    {
        // pow(x, n) creates a monomial, usable in polynomial expressions
        auto p1 = pow(x, 5);
        CLPOLY_ASSERT_EQ(p1.deg(), (int64_t)5);
        CLPOLY_ASSERT_EQ(p1.deg(x), (int64_t)5);

        // pow(x, 0) = empty monomial (constant 1 in monomial context)
        auto p0 = pow(x, 0);
        CLPOLY_ASSERT_TRUE(p0.empty());
        CLPOLY_ASSERT_EQ(p0.deg(), (int64_t)0);

        // pow(x, 1)
        auto p1x = pow(x, 1);
        CLPOLY_ASSERT_EQ(p1x.deg(), (int64_t)1);

        // Use in polynomial arithmetic
        polynomial_ZZ f = pow(x, 4) + pow(x, 2) + 1;
        polynomial_ZZ g = pow(x, 2) + 1;
        auto product = f - g;
        polynomial_ZZ expected = pow(x, 4);
        CLPOLY_ASSERT_EQ(product, expected);

        // Multi-variable pow
        polynomial_ZZ h = pow(x, 2) * pow(y, 3) + pow(z, 4);
        CLPOLY_ASSERT_EQ(degree(h, x), (int64_t)2);
        CLPOLY_ASSERT_EQ(degree(h, y), (int64_t)3);
        CLPOLY_ASSERT_EQ(degree(h, z), (int64_t)4);
    }

    return clpoly_test::test_summary();
}
