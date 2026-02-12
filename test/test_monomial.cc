#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"

using namespace clpoly;

int main() {
    variable x("x"), y("y"), z("z");

    // ======== Construction ========
    CLPOLY_TEST("monomial_default");
    {
        monomial m;
        CLPOLY_ASSERT_TRUE(m.empty());
        CLPOLY_ASSERT_EQ(m.deg(), (int64_t)0);
    }

    CLPOLY_TEST("monomial_from_variable");
    {
        monomial m = pow(x, 3);
        CLPOLY_ASSERT_EQ(m.deg(), (int64_t)3);
        CLPOLY_ASSERT_FALSE(m.empty());
    }

    CLPOLY_TEST("monomial_multivariate");
    {
        monomial m = pow(x, 2) * pow(y, 3);
        CLPOLY_ASSERT_EQ(m.deg(), (int64_t)5);
        CLPOLY_ASSERT_EQ(m.deg(x), (int64_t)2);
        CLPOLY_ASSERT_EQ(m.deg(y), (int64_t)3);
        CLPOLY_ASSERT_EQ(m.deg(z), (int64_t)0);
    }

    // ======== Multiplication ========
    CLPOLY_TEST("monomial_multiply");
    {
        monomial m1 = pow(x, 2) * y;
        monomial m2 = x * pow(y, 3);
        monomial product = m1 * m2;
        CLPOLY_ASSERT_EQ(product.deg(x), (int64_t)3);
        CLPOLY_ASSERT_EQ(product.deg(y), (int64_t)4);
    }

    // ======== Power ========
    CLPOLY_TEST("monomial_pow");
    {
        monomial m = x * y;
        monomial m2 = pow(m, 3);
        CLPOLY_ASSERT_EQ(m2.deg(x), (int64_t)3);
        CLPOLY_ASSERT_EQ(m2.deg(y), (int64_t)3);
        CLPOLY_ASSERT_EQ(m2.deg(), (int64_t)6);
    }

    // ======== GCD of monomials ========
    CLPOLY_TEST("monomial_gcd");
    {
        monomial m1({{x, 3}, {y, 2}});
        monomial m2({{x, 1}, {y, 4}, {z, 2}});
        auto g = gcd(m1, m2);
        CLPOLY_ASSERT_EQ(g.deg(x), (int64_t)1);
        CLPOLY_ASSERT_EQ(g.deg(y), (int64_t)2);
        CLPOLY_ASSERT_EQ(g.deg(z), (int64_t)0);
    }

    // ======== Ordering: grlex (default) ========
    CLPOLY_TEST("monomial_ordering_grlex");
    {
        monomial m1 = pow(x, 3);           // deg 3
        monomial m2 = pow(x, 2) * y;       // deg 3
        monomial m3 = pow(x, 2);           // deg 2
        // In grlex, higher total degree comes first
        grlex comp;
        CLPOLY_ASSERT_TRUE(comp(m1, m3));  // deg 3 > deg 2
        CLPOLY_ASSERT_TRUE(comp(m2, m3));  // deg 3 > deg 2
    }

    // ======== Empty monomial (constant) ========
    CLPOLY_TEST("monomial_empty");
    {
        monomial m;
        CLPOLY_ASSERT_TRUE(m.empty());
        CLPOLY_ASSERT_EQ(m.deg(), (int64_t)0);
        // Multiplying by empty should give the other
        monomial m2 = pow(x, 2);
        monomial product = m * m2;
        CLPOLY_ASSERT_EQ(product.deg(), m2.deg());
    }

    // ======== Division ========
    CLPOLY_TEST("monomial_division");
    {
        monomial m1 = pow(x, 3) * pow(y, 2);
        monomial m2 = pow(x, 1) * y;
        monomial q = m1 / m2;
        CLPOLY_ASSERT_EQ(q.deg(x), (int64_t)2);
        CLPOLY_ASSERT_EQ(q.deg(y), (int64_t)1);
    }

    // ======== Regression: pow(multivar_monomial) must update __deg ========
    CLPOLY_TEST("monomial_pow_multivar_deg");
    {
        // Bug: pow(basic_monomial,n) 修改指数后未调用 re_deg()
        monomial m3v = pow(x, 2) * y * pow(z, 3);  // deg = 6
        monomial m3v_p2 = pow(m3v, 2);              // should be deg 12
        CLPOLY_ASSERT_EQ(m3v_p2.deg(x), (int64_t)4);
        CLPOLY_ASSERT_EQ(m3v_p2.deg(y), (int64_t)2);
        CLPOLY_ASSERT_EQ(m3v_p2.deg(z), (int64_t)6);
        CLPOLY_ASSERT_EQ(m3v_p2.deg(), (int64_t)12);

        // pow(m, 0) should give empty monomial with deg 0
        monomial m0 = pow(m3v, 0);
        CLPOLY_ASSERT_TRUE(m0.empty());
        CLPOLY_ASSERT_EQ(m0.deg(), (int64_t)0);

        // grlex ordering must be consistent after pow
        grlex comp;
        monomial small = pow(x, 2);  // deg 2
        CLPOLY_ASSERT_TRUE(comp(m3v_p2, small));  // deg 12 > deg 2
    }

    // ======== Push_back and normalization ========
    CLPOLY_TEST("monomial_push_back_normalization");
    {
        monomial m;
        m.push_back({x, 2});
        m.push_back({y, 3});
        CLPOLY_ASSERT_TRUE(m.is_normal());
        CLPOLY_ASSERT_EQ(m.deg(), (int64_t)5);

        // push_back out of order -> not normal
        monomial m2;
        m2.push_back({y, 3});
        m2.push_back({x, 2});
        CLPOLY_ASSERT_FALSE(m2.is_normal());
        m2.normalization();
        CLPOLY_ASSERT_TRUE(m2.is_normal());
        CLPOLY_ASSERT_EQ(m, m2);
    }

    // ======== str() ========
    CLPOLY_TEST("monomial_str");
    {
        monomial m = pow(x, 2) * y;
        std::string s = m.str();
        CLPOLY_ASSERT_FALSE(s.empty());
    }

    // ======== Move semantics ========
    CLPOLY_TEST("monomial_move");
    {
        monomial m1 = pow(x, 3) * pow(y, 2);
        int64_t orig_deg = m1.deg();
        monomial m2(std::move(m1));
        CLPOLY_ASSERT_EQ(m2.deg(), orig_deg);
        CLPOLY_ASSERT_TRUE(m1.empty());

        monomial m3;
        m3 = std::move(m2);
        CLPOLY_ASSERT_EQ(m3.deg(), orig_deg);
        CLPOLY_ASSERT_TRUE(m2.empty());
    }

    // ======== Construction from vector ========
    CLPOLY_TEST("monomial_from_vector");
    {
        std::vector<std::pair<variable, int64_t>> v = {{x, 1}, {y, 2}};
        monomial m(std::move(v));
        CLPOLY_ASSERT_EQ(m.deg(x), (int64_t)1);
        CLPOLY_ASSERT_EQ(m.deg(y), (int64_t)2);
        CLPOLY_ASSERT_EQ(m.deg(), (int64_t)3);
    }

    // ======== Equality ========
    CLPOLY_TEST("monomial_equality");
    {
        monomial m1 = pow(x, 2) * y;
        monomial m2 = pow(x, 2) * y;
        monomial m3 = pow(x, 2) * pow(y, 2);
        CLPOLY_ASSERT_TRUE(m1 == m2);
        CLPOLY_ASSERT_TRUE(m1 != m3);
    }

    return clpoly_test::test_summary();
}
