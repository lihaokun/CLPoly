#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include <map>
#include <unordered_map>

using namespace clpoly;

int main() {
    // ======== Construction ========
    CLPOLY_TEST("variable_default");
    {
        variable v;
        CLPOLY_ASSERT_EQ(v.serial(), (std::size_t)0);
    }

    CLPOLY_TEST("variable_named");
    {
        variable v("testvar_a");
        CLPOLY_ASSERT_EQ(v.name(), std::string("testvar_a"));
        CLPOLY_ASSERT(v.serial() > 0);
    }

    // ======== Equality / Comparison ========
    CLPOLY_TEST("variable_equality");
    {
        variable v1("testvar_eq1");
        variable v2("testvar_eq1");
        variable v3("testvar_eq2");
        CLPOLY_ASSERT_TRUE(v1 == v2);
        CLPOLY_ASSERT_TRUE(v1 != v3);
    }

    CLPOLY_TEST("variable_ordering");
    {
        variable v1("testvar_ord1");
        variable v2("testvar_ord2");
        // ordering is by serial number
        CLPOLY_ASSERT_TRUE((v1 < v2) || (v2 < v1));
        CLPOLY_ASSERT_FALSE(v1 < v1);
    }

    // ======== Copy ========
    CLPOLY_TEST("variable_copy");
    {
        variable v1("testvar_cp");
        variable v2 = v1;
        CLPOLY_ASSERT_EQ(v1, v2);
        CLPOLY_ASSERT_EQ(v1.serial(), v2.serial());
    }

    // ======== Get by serial ========
    CLPOLY_TEST("variable_get_by_serial");
    {
        variable v1("testvar_ser");
        auto serial = v1.serial();
        variable v2 = variable::get_variable(serial);
        CLPOLY_ASSERT_EQ(v1, v2);
    }

    // ======== Get by name (idempotent) ========
    CLPOLY_TEST("variable_get_by_name");
    {
        variable v1("testvar_nm");
        variable v2 = variable::get_variable("testvar_nm");
        CLPOLY_ASSERT_EQ(v1, v2);
    }

    // ======== Used in std::map ========
    CLPOLY_TEST("variable_in_map");
    {
        variable v1("testvar_m1");
        variable v2("testvar_m2");
        std::map<variable, int> m;
        m[v1] = 10;
        m[v2] = 20;
        CLPOLY_ASSERT_EQ(m[v1], 10);
        CLPOLY_ASSERT_EQ(m[v2], 20);
    }

    // ======== Used in std::unordered_map (hashing) ========
    CLPOLY_TEST("variable_in_unordered_map");
    {
        variable v1("testvar_um1");
        variable v2("testvar_um2");
        std::unordered_map<variable, int> m;
        m[v1] = 100;
        m[v2] = 200;
        CLPOLY_ASSERT_EQ(m[v1], 100);
        CLPOLY_ASSERT_EQ(m[v2], 200);
    }

    // ======== Deletion and reuse ========
    CLPOLY_TEST("variable_delete");
    {
        variable v = variable::new_variable("testvar_del");
        auto serial = v.serial();
        variable::del_variable("testvar_del");
        // After deletion, creating a new variable may reuse the serial
        variable v2 = variable::new_variable("testvar_del2");
        // v2 should have some valid serial
        CLPOLY_ASSERT(v2.serial() > 0);
        // Clean up
        variable::del_variable("testvar_del2");
    }

    return clpoly_test::test_summary();
}
