#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include "testdata_expected.hh"

using namespace clpoly;

int main() {
    variable x("x"), y("y"), z("z");

    // ======== Resultant: simple univariate ========
    CLPOLY_TEST("resultant_simple");
    {
        // x^2+x+1 and x^3-1 share roots (x^2+x+1 divides x^3-1)
        auto f = testdata::res_f1();
        auto g = testdata::res_g1();
        auto res = resultant(f, g, x);
        auto expected_val = testdata::res_result1_val();
        // Result should be a constant (number)
        CLPOLY_ASSERT_TRUE(is_number(res));
        if (!res.empty())
            CLPOLY_ASSERT_EQ(res.front().second, expected_val);
        else
            CLPOLY_ASSERT_EQ(ZZ(0), expected_val);
    }

    // ======== Resultant: multivariate ========
    CLPOLY_TEST("resultant_multivariate");
    {
        auto f = testdata::res_f2();
        auto g = testdata::res_g2();
        auto res = resultant(f, g, x);
        auto expected = testdata::res_result2();
        CLPOLY_ASSERT_EQ(res, expected);
    }

    // ======== Resultant: common factor => resultant is 0 ========
    CLPOLY_TEST("resultant_common_factor");
    {
        auto f = testdata::res_f3_common();
        auto g = testdata::res_g3_common();
        auto res = resultant(f, g, x);
        polynomial_ZZ zero;
        CLPOLY_ASSERT_EQ(res, zero);
    }

    // ======== Resultant: if gcd(f,g) != 1 then resultant == 0 ========
    CLPOLY_TEST("resultant_gcd_nontriv_implies_zero");
    {
        polynomial_ZZ common = pow(x,2) + y;
        polynomial_ZZ f = common * (x + 1);
        polynomial_ZZ g = common * (x - y + 2);
        auto res = resultant(f, g, x);
        polynomial_ZZ zero;
        CLPOLY_ASSERT_EQ(res, zero);
    }

    // ======== Resultant: coprime polynomials => resultant != 0 ========
    CLPOLY_TEST("resultant_coprime_nonzero");
    {
        polynomial_ZZ f = pow(x,2) + 1;
        polynomial_ZZ g = x + 2;
        auto res = resultant(f, g, x);
        CLPOLY_ASSERT_TRUE(is_number(res));
        CLPOLY_ASSERT_FALSE(res.empty());
        // res(x^2+1, x+2) = 2^2+1 = 5
        CLPOLY_ASSERT_EQ(res.front().second, ZZ(5));
    }

    // ======== Resultant: complex 3-variable ========
    CLPOLY_TEST("resultant_complex");
    {
        auto f = testdata::complex_f1();
        auto g = testdata::complex_g1();
        auto res = resultant(f, g, x);
        auto expected = testdata::complex_res_result();
        CLPOLY_ASSERT_EQ(res, expected);
    }

    // ======== Discriminant: repeated root => discriminant == 0 ========
    CLPOLY_TEST("discriminant_repeated_root");
    {
        auto f = testdata::disc_f1_repeated();
        auto disc = discriminant(f, x);
        auto expected_val = testdata::disc_result1_val();
        CLPOLY_ASSERT_TRUE(is_number(disc));
        if (!disc.empty())
            CLPOLY_ASSERT_EQ(disc.front().second, expected_val);
        else
            CLPOLY_ASSERT_EQ(ZZ(0), expected_val);
    }

    // ======== Discriminant: no repeated roots ========
    CLPOLY_TEST("discriminant_no_repeat");
    {
        auto f = testdata::disc_f2();
        auto disc = discriminant(f, x);
        auto expected_val = testdata::disc_result2_val();
        CLPOLY_ASSERT_TRUE(is_number(disc));
        CLPOLY_ASSERT_EQ(disc.front().second, expected_val);
    }

    // ======== Discriminant: multivariate ========
    CLPOLY_TEST("discriminant_multivariate");
    {
        auto f = testdata::disc_f3_mv();
        auto disc = discriminant(f, x);
        auto expected = testdata::disc_result3_mv();
        CLPOLY_ASSERT_EQ(disc, expected);
    }

    // ======== Discriminant identity: f has repeated root <=> disc(f) == 0 ========
    CLPOLY_TEST("discriminant_identity");
    {
        // f = (x-1)^2*(x+2) has repeated root
        polynomial_ZZ f = pow(x,3) - 3*x + 2;
        auto disc = discriminant(f, x);
        polynomial_ZZ zero;
        CLPOLY_ASSERT_EQ(disc, zero);

        // g = x^2 - 2 has no repeated root
        polynomial_ZZ g = pow(x,2) - 2;
        auto disc_g = discriminant(g, x);
        CLPOLY_ASSERT_NE(disc_g, zero);
    }

    // ======== Subresultant ========
    CLPOLY_TEST("subresultant_basic");
    {
        // Use coprime polynomials so resultant is nonzero
        polynomial_ZZ f = pow(x,3) + x + 1;
        polynomial_ZZ g = pow(x,2) + x + 1;
        auto sres = subresultant(f, g, x);
        // sres should be a non-empty vector
        CLPOLY_ASSERT(sres.size() > 0);
        // S_0 = resultant for coprime polynomials
        auto res = resultant(f, g, x);
        // The first entry (index 0) should be S_0 = resultant
        CLPOLY_ASSERT_EQ(sres[0], res);
    }

    return clpoly_test::test_summary();
}
