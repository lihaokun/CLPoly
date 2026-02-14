#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include <sstream>

using namespace clpoly;

int main() {
    // ================================================================
    //  ZZ (mpz_class) tests
    // ================================================================

    CLPOLY_TEST("ZZ_construction");
    {
        ZZ a(0);
        CLPOLY_ASSERT_EQ(a, ZZ(0));

        ZZ b(123);
        CLPOLY_ASSERT_EQ(b, ZZ(123));

        ZZ c("-999");
        CLPOLY_ASSERT_EQ(c, ZZ(-999));

        ZZ big("123456789012345678901234567890");
        CLPOLY_ASSERT_EQ(big, ZZ("123456789012345678901234567890"));
    }

    CLPOLY_TEST("ZZ_arithmetic");
    {
        ZZ a(10), b(3);
        CLPOLY_ASSERT_EQ(ZZ(a + b), ZZ(13));
        CLPOLY_ASSERT_EQ(ZZ(a - b), ZZ(7));
        CLPOLY_ASSERT_EQ(ZZ(a * b), ZZ(30));
        CLPOLY_ASSERT_EQ(ZZ(a / b), ZZ(3));
    }

    CLPOLY_TEST("ZZ_pow");
    {
        CLPOLY_ASSERT_EQ(pow(ZZ(2), 10), ZZ(1024));
        CLPOLY_ASSERT_EQ(pow(ZZ(3), 0), ZZ(1));
        CLPOLY_ASSERT_EQ(pow(ZZ(-2), 3), ZZ(-8));
    }

    CLPOLY_TEST("ZZ_addmul");
    {
        ZZ op(5), op1(3), op2(4);
        addmul(op, op1, op2);
        CLPOLY_ASSERT_EQ(op, ZZ(5 + 3 * 4));
    }

    CLPOLY_TEST("ZZ_submul");
    {
        ZZ op(20), op1(3), op2(4);
        submul(op, op1, op2);
        CLPOLY_ASSERT_EQ(op, ZZ(20 - 3 * 4));
    }

    CLPOLY_TEST("ZZ_div_floor");
    {
        ZZ op, op1(7), op2(3);
        __div(op, op1, op2);
        CLPOLY_ASSERT_EQ(op, ZZ(2));

        ZZ op3, op4(-7), op5(3);
        __div(op3, op4, op5);
        CLPOLY_ASSERT_EQ(op3, ZZ(-3)); // floor division
    }

    CLPOLY_TEST("ZZ_div_with_remainder");
    {
        ZZ q, r, a(17), b(5);
        __div(q, r, a, b);
        CLPOLY_ASSERT_EQ(a, ZZ(q * b + r));
        CLPOLY_ASSERT_EQ(q, ZZ(3));
        CLPOLY_ASSERT_EQ(r, ZZ(2));

        ZZ q2, r2, a2(-17), b2(5);
        __div(q2, r2, a2, b2);
        CLPOLY_ASSERT_EQ(a2, ZZ(q2 * b2 + r2));
    }

    CLPOLY_TEST("ZZ_set_zero");
    {
        ZZ a(42);
        set_zero(a);
        CLPOLY_ASSERT_EQ(a, ZZ(0));
    }

    CLPOLY_TEST("ZZ_zore_check");
    {
        CLPOLY_ASSERT_TRUE(zore_check<ZZ>()(ZZ(0)));
        CLPOLY_ASSERT_FALSE(zore_check<ZZ>()(ZZ(1)));
        CLPOLY_ASSERT_FALSE(zore_check<ZZ>()(ZZ(-1)));
    }

    CLPOLY_TEST("ZZ_hash_value");
    {
        ZZ a(12345), b(12345), c(54321);
        CLPOLY_ASSERT_EQ(hash_value(a), hash_value(b));
        // Different values should (very likely) have different hashes
        CLPOLY_ASSERT_NE(hash_value(a), hash_value(c));
    }

    CLPOLY_TEST("ZZ_sizeinbase");
    {
        ZZ a(255);
        CLPOLY_ASSERT_TRUE(sizeinbase(a, 2) >= 8);
        ZZ b(1000);
        CLPOLY_ASSERT_EQ(sizeinbase(b, 10), (size_t)4);
    }

    // ================================================================
    //  QQ (mpq_class) tests
    // ================================================================

    CLPOLY_TEST("QQ_construction");
    {
        QQ a(0);
        CLPOLY_ASSERT_EQ(a, QQ(0));

        QQ b(3, 4);
        CLPOLY_ASSERT_EQ(b, QQ(3, 4));

        // Auto-canonicalization: 6/8 -> 3/4
        QQ c(6, 8);
        c.canonicalize();
        CLPOLY_ASSERT_EQ(c, QQ(3, 4));
    }

    CLPOLY_TEST("QQ_arithmetic");
    {
        QQ a(1, 2), b(1, 3);
        QQ sum = a + b;
        sum.canonicalize();
        CLPOLY_ASSERT_EQ(sum, QQ(5, 6));

        QQ diff = a - b;
        diff.canonicalize();
        CLPOLY_ASSERT_EQ(diff, QQ(1, 6));

        QQ prod = a * b;
        prod.canonicalize();
        CLPOLY_ASSERT_EQ(prod, QQ(1, 6));

        QQ quot = a / b;
        quot.canonicalize();
        CLPOLY_ASSERT_EQ(quot, QQ(3, 2));
    }

    CLPOLY_TEST("QQ_pow");
    {
        QQ r = pow(QQ(2, 3), (uint64_t)3);
        r.canonicalize();
        CLPOLY_ASSERT_EQ(r, QQ(8, 27));
    }

    CLPOLY_TEST("QQ_div_from_ZZ");
    {
        QQ op;
        __div(op, ZZ(3), ZZ(7));
        op.canonicalize();
        CLPOLY_ASSERT_EQ(op, QQ(3, 7));
    }

    CLPOLY_TEST("QQ_zore_check");
    {
        CLPOLY_ASSERT_TRUE(zore_check<QQ>()(QQ(0)));
        CLPOLY_ASSERT_FALSE(zore_check<QQ>()(QQ(1, 2)));
    }

    CLPOLY_TEST("QQ_hash_value");
    {
        QQ a(3, 4), b(3, 4), c(4, 5);
        CLPOLY_ASSERT_EQ(hash_value(a), hash_value(b));
        CLPOLY_ASSERT_NE(hash_value(a), hash_value(c));
    }

    // ================================================================
    //  Zp (modular arithmetic) tests
    // ================================================================

    CLPOLY_TEST("Zp_construction");
    {
        Zp a(7);
        CLPOLY_ASSERT_EQ(a.number(), (uint64_t)0);
        CLPOLY_ASSERT_EQ(a.prime(), (uint32_t)7);

        Zp b((uint64_t)10, 7);
        CLPOLY_ASSERT_EQ(b.number(), (uint64_t)3);

        Zp c((int64_t)10, 7);
        CLPOLY_ASSERT_EQ(c.number(), (uint64_t)3);

        Zp d(10, 7);
        CLPOLY_ASSERT_EQ(d.number(), (uint64_t)3);

        Zp e(ZZ(10), 7);
        CLPOLY_ASSERT_EQ(e.number(), (uint64_t)3);
    }

    CLPOLY_TEST("Zp_negative_construction");
    {
        Zp a(-1, 7);
        CLPOLY_ASSERT_EQ(a.number(), (uint64_t)6);

        Zp b((int64_t)-3, 7);
        CLPOLY_ASSERT_EQ(b.number(), (uint64_t)4);
    }

    CLPOLY_TEST("Zp_arithmetic");
    {
        Zp a(3, 7), b(5, 7);
        Zp sum = a + b;
        CLPOLY_ASSERT_EQ(sum.number(), (uint64_t)1); // (3+5)%7=1

        Zp diff = a - b;
        CLPOLY_ASSERT_EQ(diff.number(), (uint64_t)5); // (3-5+7)%7=5

        Zp prod = a * b;
        CLPOLY_ASSERT_EQ(prod.number(), (uint64_t)1); // (3*5)%7=15%7=1

        Zp quot = a / b;
        // a/b = a * inv(b); verify by checking quot * b == a
        Zp check = quot * b;
        CLPOLY_ASSERT_EQ(check.number(), a.number());
    }

    CLPOLY_TEST("Zp_inv");
    {
        uint32_t p = 7;
        for (uint64_t i = 1; i < p; ++i) {
            Zp a(i, p);
            Zp ai = a.inv();
            Zp prod = a * ai;
            CLPOLY_ASSERT_EQ(prod.number(), (uint64_t)1);
        }
    }

    CLPOLY_TEST("Zp_pow_fermat");
    {
        uint32_t p = 13;
        for (uint64_t i = 1; i < p; ++i) {
            Zp a(i, p);
            Zp result = pow(a, (int64_t)(p - 1));
            CLPOLY_ASSERT_EQ(result.number(), (uint64_t)1);
        }
    }

    CLPOLY_TEST("Zp_addmul");
    {
        Zp op(2, 7), op1(3, 7), op2(4, 7);
        addmul(op, op1, op2);
        CLPOLY_ASSERT_EQ(op.number(), (uint64_t)((2 + 3 * 4) % 7));
    }

    CLPOLY_TEST("Zp_submul");
    {
        Zp op(6, 7), op1(3, 7), op2(4, 7);
        submul(op, op1, op2);
        // 6 + 7 - (3*4)%7 = 6 + 7 - 5 = 8, 8%7 = 1
        CLPOLY_ASSERT_EQ(op.number(), (uint64_t)((6 + 7 - (3 * 4) % 7) % 7));
    }

    CLPOLY_TEST("Zp_assignment");
    {
        Zp a(7);
        a = (int64_t)10;
        CLPOLY_ASSERT_EQ(a.number(), (uint64_t)3);

        a = ZZ(15);
        CLPOLY_ASSERT_EQ(a.number(), (uint64_t)1); // 15%7=1
    }

    CLPOLY_TEST("Zp_normalization");
    {
        Zp a(3, 7);
        a.number() = 20; // bypass modular reduction
        a.normalization();
        CLPOLY_ASSERT_EQ(a.number(), (uint64_t)(20 % 7));
    }

    CLPOLY_TEST("Zp_hash_value");
    {
        Zp a(3, 7), b(3, 7), c(4, 7);
        CLPOLY_ASSERT_EQ(hash_value(a), hash_value(b));
        CLPOLY_ASSERT_NE(hash_value(a), hash_value(c));
    }

    CLPOLY_TEST("Zp_output");
    {
        Zp a(5, 7);
        std::ostringstream ss;
        ss << a;
        CLPOLY_ASSERT_EQ(ss.str(), std::string("5"));
    }

    return clpoly_test::test_summary();
}
