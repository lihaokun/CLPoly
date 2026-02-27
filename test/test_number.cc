#include <clpoly/clpoly.hh>
#include "clpoly_test.hh"
#include <sstream>
#include <climits>
#include <stdexcept>

using namespace clpoly;

int main() {
    // ================================================================
    //  ZZ construction
    // ================================================================

    CLPOLY_TEST("ZZ_construction");
    {
        ZZ a;
        CLPOLY_ASSERT_EQ(a, ZZ(0));

        ZZ b(0);
        CLPOLY_ASSERT_EQ(b, ZZ(0));

        ZZ c(123);
        CLPOLY_ASSERT_EQ(c, ZZ(123));

        ZZ d(-999);
        CLPOLY_ASSERT_EQ(d, ZZ(-999));

        ZZ e("-999");
        CLPOLY_ASSERT_EQ(e, ZZ(-999));

        ZZ big("123456789012345678901234567890");
        CLPOLY_ASSERT_EQ(big, ZZ("123456789012345678901234567890"));

        // std::string constructor
        ZZ f(std::string("42"));
        CLPOLY_ASSERT_EQ(f, ZZ(42));
    }

    // Regression test for issue #3
    CLPOLY_TEST("ZZ_int64_construction");
    {
        int64_t v1 = 1234567890123LL;
        ZZ a(v1);
        CLPOLY_ASSERT_EQ(a.get_si(), v1);

        long long v2 = -9876543210LL;
        ZZ b(v2);
        CLPOLY_ASSERT_EQ(b.get_si(), v2);

        ZZ c = a + b;
        ZZ expected(v1 + v2);
        CLPOLY_ASSERT_EQ(c, expected);

        ZZ d;
        d = 42LL;
        CLPOLY_ASSERT_EQ(d, ZZ(42));
    }

    CLPOLY_TEST("ZZ_unsigned_construction");
    {
        ZZ a((unsigned int)42);
        CLPOLY_ASSERT_EQ(a, ZZ(42));

        ZZ b((unsigned long)123456);
        CLPOLY_ASSERT_EQ(b, ZZ(123456));

        ZZ c((unsigned long long)99999);
        CLPOLY_ASSERT_EQ(c, ZZ(99999));

        // Large unsigned long long > INT64_MAX → promote
        unsigned long long big_val = static_cast<unsigned long long>(INT64_MAX) + 1ULL;
        ZZ d(big_val);
        CLPOLY_ASSERT_TRUE(d > ZZ(0));
        CLPOLY_ASSERT_TRUE(d > ZZ(INT64_MAX));
    }

    CLPOLY_TEST("ZZ_copy_move");
    {
        ZZ a(42);
        ZZ b(a);  // copy
        CLPOLY_ASSERT_EQ(a, b);

        ZZ big("999999999999999999999999999999");
        ZZ big2(big);  // copy big
        CLPOLY_ASSERT_EQ(big, big2);

        ZZ c(std::move(b));  // move
        CLPOLY_ASSERT_EQ(c, ZZ(42));

        ZZ d;
        d = a;  // copy assign
        CLPOLY_ASSERT_EQ(d, ZZ(42));

        ZZ e;
        e = std::move(d);  // move assign
        CLPOLY_ASSERT_EQ(e, ZZ(42));
    }

    // ================================================================
    //  ZZ state query
    // ================================================================

    CLPOLY_TEST("ZZ_is_zero_is_one");
    {
        CLPOLY_ASSERT_TRUE(ZZ(0).is_zero());
        CLPOLY_ASSERT_FALSE(ZZ(1).is_zero());
        CLPOLY_ASSERT_FALSE(ZZ(-1).is_zero());

        CLPOLY_ASSERT_TRUE(ZZ(1).is_one());
        CLPOLY_ASSERT_FALSE(ZZ(0).is_one());
        CLPOLY_ASSERT_FALSE(ZZ(2).is_one());
        CLPOLY_ASSERT_FALSE(ZZ(-1).is_one());
    }

    CLPOLY_TEST("ZZ_sgn");
    {
        CLPOLY_ASSERT_EQ(sgn(ZZ(42)), 1);
        CLPOLY_ASSERT_EQ(sgn(ZZ(0)), 0);
        CLPOLY_ASSERT_EQ(sgn(ZZ(-42)), -1);

        // big integer sgn
        CLPOLY_ASSERT_EQ(sgn(ZZ("999999999999999999999999999999")), 1);
        CLPOLY_ASSERT_EQ(sgn(ZZ("-999999999999999999999999999999")), -1);
    }

    CLPOLY_TEST("ZZ_bool");
    {
        CLPOLY_ASSERT_FALSE(static_cast<bool>(ZZ(0)));
        CLPOLY_ASSERT_TRUE(static_cast<bool>(ZZ(1)));
        CLPOLY_ASSERT_TRUE(static_cast<bool>(ZZ(-1)));
        CLPOLY_ASSERT_TRUE(static_cast<bool>(ZZ("99999999999999999999999")));
    }

    // ================================================================
    //  ZZ comparison
    // ================================================================

    CLPOLY_TEST("ZZ_comparison");
    {
        // small vs small
        CLPOLY_ASSERT_TRUE(ZZ(1) < ZZ(2));
        CLPOLY_ASSERT_FALSE(ZZ(2) < ZZ(1));
        CLPOLY_ASSERT_FALSE(ZZ(2) < ZZ(2));

        CLPOLY_ASSERT_TRUE(ZZ(2) > ZZ(1));
        CLPOLY_ASSERT_TRUE(ZZ(1) <= ZZ(2));
        CLPOLY_ASSERT_TRUE(ZZ(2) <= ZZ(2));
        CLPOLY_ASSERT_TRUE(ZZ(2) >= ZZ(1));
        CLPOLY_ASSERT_TRUE(ZZ(2) >= ZZ(2));

        CLPOLY_ASSERT_TRUE(ZZ(-5) < ZZ(5));
        CLPOLY_ASSERT_TRUE(ZZ(-5) < ZZ(-3));

        // small vs big
        ZZ big("999999999999999999999999999999");
        CLPOLY_ASSERT_TRUE(ZZ(1) < big);
        CLPOLY_ASSERT_TRUE(big > ZZ(1));
        CLPOLY_ASSERT_TRUE(ZZ(-1) < big);

        // big vs big
        ZZ big2("999999999999999999999999999998");
        CLPOLY_ASSERT_TRUE(big2 < big);
        CLPOLY_ASSERT_TRUE(big > big2);
        CLPOLY_ASSERT_TRUE(big >= big);
        CLPOLY_ASSERT_EQ(big, big);
    }

    // ================================================================
    //  ZZ arithmetic
    // ================================================================

    CLPOLY_TEST("ZZ_arithmetic");
    {
        ZZ a(10), b(3);
        CLPOLY_ASSERT_EQ(a + b, ZZ(13));
        CLPOLY_ASSERT_EQ(a - b, ZZ(7));
        CLPOLY_ASSERT_EQ(a * b, ZZ(30));
        CLPOLY_ASSERT_EQ(a / b, ZZ(3));  // truncated division
    }

    CLPOLY_TEST("ZZ_modulo");
    {
        CLPOLY_ASSERT_EQ(ZZ(17) % ZZ(5), ZZ(2));
        CLPOLY_ASSERT_EQ(ZZ(-17) % ZZ(5), ZZ(-2));  // truncated
        CLPOLY_ASSERT_EQ(ZZ(20) % ZZ(5), ZZ(0));
    }

    CLPOLY_TEST("ZZ_negate");
    {
        CLPOLY_ASSERT_EQ(-ZZ(42), ZZ(-42));
        CLPOLY_ASSERT_EQ(-ZZ(-42), ZZ(42));
        CLPOLY_ASSERT_EQ(-ZZ(0), ZZ(0));

        // INT64_MIN: -INT64_MIN overflows int64_t, must promote
        ZZ neg_min = -ZZ(INT64_MIN);
        ZZ pos_min("9223372036854775808");  // 2^63
        CLPOLY_ASSERT_EQ(neg_min, pos_min);
    }

    CLPOLY_TEST("ZZ_compound_assignment");
    {
        ZZ a(10);
        a += ZZ(5);
        CLPOLY_ASSERT_EQ(a, ZZ(15));

        a -= ZZ(3);
        CLPOLY_ASSERT_EQ(a, ZZ(12));

        a *= ZZ(4);
        CLPOLY_ASSERT_EQ(a, ZZ(48));

        a /= ZZ(6);
        CLPOLY_ASSERT_EQ(a, ZZ(8));

        a %= ZZ(3);
        CLPOLY_ASSERT_EQ(a, ZZ(2));
    }

    // ================================================================
    //  ZZ overflow / promote-demote (core feature)
    // ================================================================

    CLPOLY_TEST("ZZ_overflow_add");
    {
        ZZ a(INT64_MAX);
        ZZ b(1);
        ZZ c = a + b;
        // Should promote to big integer
        CLPOLY_ASSERT_TRUE(c > ZZ(INT64_MAX));
        CLPOLY_ASSERT_EQ(c - b, a);

        // Negative overflow
        ZZ d(INT64_MIN);
        ZZ e = d - b;
        CLPOLY_ASSERT_TRUE(e < ZZ(INT64_MIN));
        CLPOLY_ASSERT_EQ(e + b, d);
    }

    CLPOLY_TEST("ZZ_overflow_mul");
    {
        ZZ a(INT64_MAX);
        ZZ b(2);
        ZZ c = a * b;
        CLPOLY_ASSERT_TRUE(c > ZZ(INT64_MAX));
        CLPOLY_ASSERT_EQ(c / b, a);
    }

    CLPOLY_TEST("ZZ_overflow_div_min");
    {
        // INT64_MIN / -1 = -INT64_MIN, overflows
        ZZ a(INT64_MIN);
        ZZ b(-1);
        ZZ c = a / b;
        CLPOLY_ASSERT_EQ(c, ZZ("9223372036854775808"));
    }

    CLPOLY_TEST("ZZ_demote");
    {
        // big → small: compute something big then bring it back to small range
        ZZ big = ZZ(INT64_MAX) + ZZ(100);
        CLPOLY_ASSERT_TRUE(big > ZZ(INT64_MAX));
        ZZ small = big - ZZ(INT64_MAX);
        CLPOLY_ASSERT_EQ(small, ZZ(100));
    }

    CLPOLY_TEST("ZZ_mixed_small_big");
    {
        ZZ big("100000000000000000000");  // 10^20
        ZZ small(42);

        ZZ sum = big + small;
        CLPOLY_ASSERT_EQ(sum, ZZ("100000000000000000042"));

        ZZ diff = big - small;
        CLPOLY_ASSERT_EQ(diff, ZZ("99999999999999999958"));

        ZZ prod = big * small;
        CLPOLY_ASSERT_EQ(prod, ZZ("4200000000000000000000"));

        ZZ quot = big / small;
        ZZ rem = big % small;
        CLPOLY_ASSERT_EQ(quot * small + rem, big);
    }

    CLPOLY_TEST("ZZ_big_big_arithmetic");
    {
        ZZ a("123456789012345678901234567890");
        ZZ b("987654321098765432109876543210");
        ZZ sum = a + b;
        CLPOLY_ASSERT_EQ(sum, ZZ("1111111110111111111011111111100"));

        ZZ prod = a * b;
        ZZ check = prod / a;
        CLPOLY_ASSERT_EQ(check, b);
    }

    // ================================================================
    //  ZZ special functions
    // ================================================================

    CLPOLY_TEST("ZZ_pow");
    {
        CLPOLY_ASSERT_EQ(pow(ZZ(2), 10), ZZ(1024));
        CLPOLY_ASSERT_EQ(pow(ZZ(3), 0), ZZ(1));
        CLPOLY_ASSERT_EQ(pow(ZZ(-2), 3), ZZ(-8));
        CLPOLY_ASSERT_EQ(pow(ZZ(0), 5), ZZ(0));
        CLPOLY_ASSERT_EQ(pow(ZZ(1), 100), ZZ(1));

        // Large power → big integer
        ZZ p = pow(ZZ(2), 64);
        CLPOLY_ASSERT_EQ(p, ZZ("18446744073709551616"));
    }

    CLPOLY_TEST("ZZ_gcd");
    {
        CLPOLY_ASSERT_EQ(gcd(ZZ(12), ZZ(8)), ZZ(4));
        CLPOLY_ASSERT_EQ(gcd(ZZ(17), ZZ(13)), ZZ(1));
        CLPOLY_ASSERT_EQ(gcd(ZZ(0), ZZ(5)), ZZ(5));
        CLPOLY_ASSERT_EQ(gcd(ZZ(5), ZZ(0)), ZZ(5));
        CLPOLY_ASSERT_EQ(gcd(ZZ(-12), ZZ(8)), ZZ(4));
        CLPOLY_ASSERT_EQ(gcd(ZZ(12), ZZ(-8)), ZZ(4));

        // big integer gcd
        ZZ a("123456789012345678901234567890");
        ZZ b("987654321098765432109876543210");
        ZZ g = gcd(a, b);
        CLPOLY_ASSERT_EQ(a % g, ZZ(0));
        CLPOLY_ASSERT_EQ(b % g, ZZ(0));
    }

    CLPOLY_TEST("ZZ_lcm");
    {
        CLPOLY_ASSERT_EQ(lcm(ZZ(4), ZZ(6)), ZZ(12));
        CLPOLY_ASSERT_EQ(lcm(ZZ(3), ZZ(7)), ZZ(21));
        CLPOLY_ASSERT_EQ(lcm(ZZ(0), ZZ(5)), ZZ(0));
        CLPOLY_ASSERT_EQ(lcm(ZZ(5), ZZ(0)), ZZ(0));

        // lcm(a,b) * gcd(a,b) == |a*b|
        ZZ a(12), b(8);
        CLPOLY_ASSERT_EQ(lcm(a, b) * gcd(a, b), abs(a * b));
    }

    CLPOLY_TEST("ZZ_abs");
    {
        CLPOLY_ASSERT_EQ(abs(ZZ(42)), ZZ(42));
        CLPOLY_ASSERT_EQ(abs(ZZ(-42)), ZZ(42));
        CLPOLY_ASSERT_EQ(abs(ZZ(0)), ZZ(0));

        // INT64_MIN
        ZZ a = abs(ZZ(INT64_MIN));
        CLPOLY_ASSERT_EQ(a, ZZ("9223372036854775808"));
    }

    CLPOLY_TEST("ZZ_addmul");
    {
        ZZ op(5), op1(3), op2(4);
        addmul(op, op1, op2);
        CLPOLY_ASSERT_EQ(op, ZZ(5 + 3 * 4));

        // addmul with big numbers
        ZZ big("100000000000000000000");
        ZZ x(999999999), y(999999999);
        addmul(big, x, y);
        CLPOLY_ASSERT_EQ(big, ZZ("100000000000000000000") + ZZ(999999999) * ZZ(999999999));
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

        // Positive case: floor == truncated
        ZZ op6, op7(10), op8(3);
        __div(op6, op7, op8);
        CLPOLY_ASSERT_EQ(op6, ZZ(3));
    }

    CLPOLY_TEST("ZZ_div_with_remainder");
    {
        ZZ q, r, a(17), b(5);
        __div(q, r, a, b);
        CLPOLY_ASSERT_EQ(a, q * b + r);
        CLPOLY_ASSERT_EQ(q, ZZ(3));
        CLPOLY_ASSERT_EQ(r, ZZ(2));

        ZZ q2, r2, a2(-17), b2(5);
        __div(q2, r2, a2, b2);
        CLPOLY_ASSERT_EQ(a2, q2 * b2 + r2);
        // floor division: remainder has same sign as divisor
        CLPOLY_ASSERT_TRUE(r2 >= ZZ(0));
    }

    // ================================================================
    //  ZZ swap
    // ================================================================

    CLPOLY_TEST("ZZ_swap");
    {
        ZZ a(42), b(99);
        a.swap(b);
        CLPOLY_ASSERT_EQ(a, ZZ(99));
        CLPOLY_ASSERT_EQ(b, ZZ(42));

        // swap with big integer
        ZZ c(7), d("999999999999999999999999999999");
        swap(c, d);  // free function
        CLPOLY_ASSERT_EQ(c, ZZ("999999999999999999999999999999"));
        CLPOLY_ASSERT_EQ(d, ZZ(7));
    }

    // ================================================================
    //  ZZ IO
    // ================================================================

    CLPOLY_TEST("ZZ_output");
    {
        std::ostringstream ss;
        ss << ZZ(42);
        CLPOLY_ASSERT_EQ(ss.str(), std::string("42"));

        std::ostringstream ss2;
        ss2 << ZZ(-123);
        CLPOLY_ASSERT_EQ(ss2.str(), std::string("-123"));

        std::ostringstream ss3;
        ss3 << ZZ(0);
        CLPOLY_ASSERT_EQ(ss3.str(), std::string("0"));

        // big integer output
        std::ostringstream ss4;
        ZZ big("123456789012345678901234567890");
        ss4 << big;
        CLPOLY_ASSERT_EQ(ss4.str(), std::string("123456789012345678901234567890"));
    }

    CLPOLY_TEST("ZZ_input");
    {
        std::istringstream is("12345");
        ZZ a;
        is >> a;
        CLPOLY_ASSERT_EQ(a, ZZ(12345));

        std::istringstream is2("-999");
        ZZ b;
        is2 >> b;
        CLPOLY_ASSERT_EQ(b, ZZ(-999));

        std::istringstream is3("123456789012345678901234567890");
        ZZ c;
        is3 >> c;
        CLPOLY_ASSERT_EQ(c, ZZ("123456789012345678901234567890"));
    }

    // ================================================================
    //  ZZ misc (existing tests kept)
    // ================================================================

    CLPOLY_TEST("ZZ_set_zero");
    {
        ZZ a(42);
        set_zero(a);
        CLPOLY_ASSERT_EQ(a, ZZ(0));

        // set_zero on big integer
        ZZ big("999999999999999999999999");
        set_zero(big);
        CLPOLY_ASSERT_EQ(big, ZZ(0));
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
        CLPOLY_ASSERT_NE(hash_value(a), hash_value(c));

        // big integer hash consistency
        ZZ big1("123456789012345678901234567890");
        ZZ big2("123456789012345678901234567890");
        CLPOLY_ASSERT_EQ(hash_value(big1), hash_value(big2));
    }

    CLPOLY_TEST("ZZ_sizeinbase");
    {
        ZZ a(255);
        CLPOLY_ASSERT_TRUE(sizeinbase(a, 2) >= 8);
        ZZ b(1000);
        CLPOLY_ASSERT_EQ(sizeinbase(b, 10), (size_t)4);
        CLPOLY_ASSERT_EQ(sizeinbase(ZZ(0), 10), (size_t)1);
    }

    CLPOLY_TEST("ZZ_fdiv_ui");
    {
        CLPOLY_ASSERT_EQ(ZZ(10).fdiv_ui(3), (uint64_t)1);
        CLPOLY_ASSERT_EQ(ZZ(-10).fdiv_ui(3), (uint64_t)2);  // floor mod
        CLPOLY_ASSERT_EQ(ZZ(0).fdiv_ui(7), (uint64_t)0);
        CLPOLY_ASSERT_EQ(ZZ(21).fdiv_ui(7), (uint64_t)0);
    }

    CLPOLY_TEST("ZZ_shift");
    {
        // --- operator<<= / operator>>= (compound, small integers) ---
        ZZ a(1);
        a <<= 4;
        CLPOLY_ASSERT_EQ(a, ZZ(16));
        a >>= 2;
        CLPOLY_ASSERT_EQ(a, ZZ(4));

        // --- binary operator<< / operator>> (non-modifying) ---
        // lvalue path (const ZZ&)
        ZZ x(3);
        CLPOLY_ASSERT_EQ(x << 3,  ZZ(24));
        CLPOLY_ASSERT_EQ(x >> 1,  ZZ(1));
        CLPOLY_ASSERT_EQ(x,       ZZ(3));       // x 不变

        // rvalue path (ZZ&&) — 典型用例：临时对象直接位移，无拷贝
        CLPOLY_ASSERT_EQ(ZZ(1) << 10, ZZ(1024));
        CLPOLY_ASSERT_EQ(ZZ(128) >> 3, ZZ(16));

        // --- 大整数路径（超出 int64 小整数范围）---
        // 1 << 64 → 需要 mpz 路径
        ZZ big = ZZ(1) << 64;
        CLPOLY_ASSERT_EQ(big, ZZ("18446744073709551616"));
        // 再右移回来
        CLPOLY_ASSERT_EQ(big >> 64, ZZ(1));

        // --- 负数右移（算术移位 / floor toward -∞）---
        CLPOLY_ASSERT_EQ(ZZ(-8) >> 1, ZZ(-4));
        CLPOLY_ASSERT_EQ(ZZ(-1) >> 1, ZZ(-1));  // floor(-0.5) = -1

        // --- 边界：移位量 >= 63（小整数路径特殊处理）---
        CLPOLY_ASSERT_EQ(ZZ(1)  >> 63, ZZ(0));
        CLPOLY_ASSERT_EQ(ZZ(-1) >> 63, ZZ(-1));
        CLPOLY_ASSERT_EQ(ZZ(1)  >> 100, ZZ(0));  // 大移位量走 mpz 路径

        // --- 零移位 ---
        CLPOLY_ASSERT_EQ(ZZ(42) << 0, ZZ(42));
        CLPOLY_ASSERT_EQ(ZZ(42) >> 0, ZZ(42));

        // --- 模拟 M5 中的 van Hoeij bound 计算：ZZ(r+1) * (ZZ(1) << (2*U_exp)) ---
        int r = 10, U_exp = 4;
        ZZ B = ZZ(r + 1) * (ZZ(1) << (2 * U_exp));
        CLPOLY_ASSERT_EQ(B, ZZ(11) * ZZ(256));   // 11 * 2^8
    }

    // ================================================================
    //  QQ construction
    // ================================================================

    CLPOLY_TEST("QQ_construction");
    {
        QQ a;
        CLPOLY_ASSERT_EQ(a, QQ(0));

        QQ b(0);
        CLPOLY_ASSERT_EQ(b, QQ(0));

        QQ c(42);
        CLPOLY_ASSERT_EQ(c, QQ(42));

        QQ d(42L);
        CLPOLY_ASSERT_EQ(d, QQ(42));

        QQ e(42LL);
        CLPOLY_ASSERT_EQ(e, QQ(42));

        QQ f(ZZ(42));
        CLPOLY_ASSERT_EQ(f, QQ(42));

        QQ g(3, 4);
        CLPOLY_ASSERT_EQ(g, QQ(3, 4));

        // Auto-canonicalization: 6/8 -> 3/4
        QQ h(6, 8);
        CLPOLY_ASSERT_EQ(h, QQ(3, 4));

        // ZZ numerator and denominator
        QQ i(ZZ(10), ZZ(4));
        CLPOLY_ASSERT_EQ(i, QQ(5, 2));

        // Negative denominator → canonicalize to positive
        QQ j(3, -4);
        CLPOLY_ASSERT_EQ(j, QQ(-3, 4));
    }

    CLPOLY_TEST("QQ_zero_denominator");
    {
        bool threw = false;
        try {
            QQ bad(1, 0);
        } catch (const std::domain_error&) {
            threw = true;
        }
        CLPOLY_ASSERT_TRUE(threw);
    }

    CLPOLY_TEST("QQ_copy_move");
    {
        QQ a(3, 4);
        QQ b(a);  // copy
        CLPOLY_ASSERT_EQ(a, b);

        QQ c(std::move(b));  // move
        CLPOLY_ASSERT_EQ(c, QQ(3, 4));

        QQ d;
        d = a;  // copy assign
        CLPOLY_ASSERT_EQ(d, QQ(3, 4));

        QQ e;
        e = std::move(d);  // move assign
        CLPOLY_ASSERT_EQ(e, QQ(3, 4));

        QQ f;
        f = 42LL;
        CLPOLY_ASSERT_EQ(f, QQ(42));
    }

    // ================================================================
    //  QQ accessors
    // ================================================================

    CLPOLY_TEST("QQ_accessors");
    {
        QQ a(3, 4);
        CLPOLY_ASSERT_EQ(a.num(), ZZ(3));
        CLPOLY_ASSERT_EQ(a.den(), ZZ(4));
        CLPOLY_ASSERT_EQ(a.get_num(), ZZ(3));
        CLPOLY_ASSERT_EQ(a.get_den(), ZZ(4));

        QQ b(-5, 7);
        CLPOLY_ASSERT_EQ(b.num(), ZZ(-5));
        CLPOLY_ASSERT_EQ(b.den(), ZZ(7));

        QQ c(0);
        CLPOLY_ASSERT_EQ(c.num(), ZZ(0));
        CLPOLY_ASSERT_EQ(c.den(), ZZ(1));
    }

    // ================================================================
    //  QQ state query
    // ================================================================

    CLPOLY_TEST("QQ_is_zero");
    {
        CLPOLY_ASSERT_TRUE(QQ(0).is_zero());
        CLPOLY_ASSERT_FALSE(QQ(1, 2).is_zero());
        CLPOLY_ASSERT_FALSE(QQ(-1).is_zero());
    }

    CLPOLY_TEST("QQ_bool");
    {
        CLPOLY_ASSERT_FALSE(static_cast<bool>(QQ(0)));
        CLPOLY_ASSERT_TRUE(static_cast<bool>(QQ(1, 2)));
        CLPOLY_ASSERT_TRUE(static_cast<bool>(QQ(-3, 4)));
    }

    CLPOLY_TEST("QQ_sgn");
    {
        CLPOLY_ASSERT_EQ(sgn(QQ(3, 4)), 1);
        CLPOLY_ASSERT_EQ(sgn(QQ(0)), 0);
        CLPOLY_ASSERT_EQ(sgn(QQ(-3, 4)), -1);
    }

    // ================================================================
    //  QQ comparison
    // ================================================================

    CLPOLY_TEST("QQ_comparison");
    {
        CLPOLY_ASSERT_TRUE(QQ(1, 3) < QQ(1, 2));
        CLPOLY_ASSERT_FALSE(QQ(1, 2) < QQ(1, 3));
        CLPOLY_ASSERT_FALSE(QQ(1, 2) < QQ(1, 2));

        CLPOLY_ASSERT_TRUE(QQ(1, 2) > QQ(1, 3));
        CLPOLY_ASSERT_TRUE(QQ(1, 3) <= QQ(1, 2));
        CLPOLY_ASSERT_TRUE(QQ(1, 2) <= QQ(1, 2));
        CLPOLY_ASSERT_TRUE(QQ(1, 2) >= QQ(1, 3));
        CLPOLY_ASSERT_TRUE(QQ(1, 2) >= QQ(1, 2));

        CLPOLY_ASSERT_TRUE(QQ(-1, 2) < QQ(1, 2));
        CLPOLY_ASSERT_TRUE(QQ(-3, 4) < QQ(-1, 4));

        CLPOLY_ASSERT_EQ(QQ(2, 4), QQ(1, 2));
        CLPOLY_ASSERT_NE(QQ(1, 2), QQ(1, 3));
    }

    // ================================================================
    //  QQ arithmetic
    // ================================================================

    CLPOLY_TEST("QQ_arithmetic");
    {
        QQ a(1, 2), b(1, 3);
        CLPOLY_ASSERT_EQ(a + b, QQ(5, 6));
        CLPOLY_ASSERT_EQ(a - b, QQ(1, 6));
        CLPOLY_ASSERT_EQ(a * b, QQ(1, 6));
        CLPOLY_ASSERT_EQ(a / b, QQ(3, 2));
    }

    CLPOLY_TEST("QQ_negate");
    {
        CLPOLY_ASSERT_EQ(-QQ(3, 4), QQ(-3, 4));
        CLPOLY_ASSERT_EQ(-QQ(-3, 4), QQ(3, 4));
        CLPOLY_ASSERT_EQ(-QQ(0), QQ(0));
    }

    CLPOLY_TEST("QQ_compound_assignment");
    {
        QQ a(1, 2);
        a += QQ(1, 3);
        CLPOLY_ASSERT_EQ(a, QQ(5, 6));

        a -= QQ(1, 6);
        CLPOLY_ASSERT_EQ(a, QQ(2, 3));

        a *= QQ(3, 4);
        CLPOLY_ASSERT_EQ(a, QQ(1, 2));

        a /= QQ(1, 4);
        CLPOLY_ASSERT_EQ(a, QQ(2));
    }

    CLPOLY_TEST("QQ_div_by_zero");
    {
        bool threw = false;
        try {
            QQ a(1, 2);
            QQ b(0);
            QQ c = a / b;
            (void)c;
        } catch (const std::domain_error&) {
            threw = true;
        }
        CLPOLY_ASSERT_TRUE(threw);
    }

    // ================================================================
    //  QQ special functions
    // ================================================================

    CLPOLY_TEST("QQ_pow");
    {
        CLPOLY_ASSERT_EQ(pow(QQ(2, 3), (uint64_t)3), QQ(8, 27));
        CLPOLY_ASSERT_EQ(pow(QQ(2, 3), (uint64_t)0), QQ(1));
        CLPOLY_ASSERT_EQ(pow(QQ(-1, 2), (uint64_t)2), QQ(1, 4));
    }

    CLPOLY_TEST("QQ_abs");
    {
        CLPOLY_ASSERT_EQ(abs(QQ(3, 4)), QQ(3, 4));
        CLPOLY_ASSERT_EQ(abs(QQ(-3, 4)), QQ(3, 4));
        CLPOLY_ASSERT_EQ(abs(QQ(0)), QQ(0));
    }

    CLPOLY_TEST("QQ_div_from_ZZ");
    {
        QQ op;
        __div(op, ZZ(3), ZZ(7));
        CLPOLY_ASSERT_EQ(op, QQ(3, 7));

        QQ op2;
        QQ r2;
        __div(op2, r2, ZZ(10), ZZ(3));
        CLPOLY_ASSERT_EQ(op2, QQ(10, 3));
        CLPOLY_ASSERT_EQ(r2, QQ(0));
    }

    // ================================================================
    //  QQ swap
    // ================================================================

    CLPOLY_TEST("QQ_swap");
    {
        QQ a(1, 2), b(3, 4);
        a.swap(b);
        CLPOLY_ASSERT_EQ(a, QQ(3, 4));
        CLPOLY_ASSERT_EQ(b, QQ(1, 2));

        swap(a, b);  // free function
        CLPOLY_ASSERT_EQ(a, QQ(1, 2));
        CLPOLY_ASSERT_EQ(b, QQ(3, 4));
    }

    // ================================================================
    //  QQ IO
    // ================================================================

    CLPOLY_TEST("QQ_output");
    {
        std::ostringstream ss;
        ss << QQ(3, 4);
        CLPOLY_ASSERT_EQ(ss.str(), std::string("3/4"));

        std::ostringstream ss2;
        ss2 << QQ(5);
        CLPOLY_ASSERT_EQ(ss2.str(), std::string("5"));

        std::ostringstream ss3;
        ss3 << QQ(0);
        CLPOLY_ASSERT_EQ(ss3.str(), std::string("0"));

        std::ostringstream ss4;
        ss4 << QQ(-7, 3);
        CLPOLY_ASSERT_EQ(ss4.str(), std::string("-7/3"));
    }

    CLPOLY_TEST("QQ_input");
    {
        std::istringstream is("3/4");
        QQ a;
        is >> a;
        CLPOLY_ASSERT_EQ(a, QQ(3, 4));

        std::istringstream is2("42");
        QQ b;
        is2 >> b;
        CLPOLY_ASSERT_EQ(b, QQ(42));
    }

    // ================================================================
    //  QQ misc (existing tests kept)
    // ================================================================

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

        // Canonical form: 6/8 and 3/4 should have same hash
        QQ d(6, 8);
        CLPOLY_ASSERT_EQ(hash_value(a), hash_value(d));
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
        uint64_t p = 7;
        for (uint64_t i = 1; i < p; ++i) {
            Zp a(i, p);
            Zp ai = a.inv();
            Zp prod = a * ai;
            CLPOLY_ASSERT_EQ(prod.number(), (uint64_t)1);
        }
    }

    CLPOLY_TEST("Zp_pow_fermat");
    {
        uint64_t p = 13;
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

    // ================================================================
    //  Zp: set_zero (回归: pair_vec_multiplies 默认构造 Zp 后调 set_zero 崩溃)
    // ================================================================

    CLPOLY_TEST("Zp_set_zero");
    {
        // 正常 Zp
        Zp a(5, 7);
        set_zero(a);
        CLPOLY_ASSERT_EQ(a.number(), (uint64_t)0);
        CLPOLY_ASSERT_EQ(a.prime(), (uint32_t)7);

        // _p==0 的默认构造 Zp (回归: 之前 set_zero 触发 assert/_p==0 除零 UB)
        Zp b;
        CLPOLY_ASSERT_EQ(b.prime(), (uint32_t)0);
        set_zero(b);
        CLPOLY_ASSERT_EQ(b.number(), (uint64_t)0);
    }

    CLPOLY_TEST("Zp_zore_check");
    {
        zore_check<Zp> chk;
        CLPOLY_ASSERT_TRUE(chk(Zp(0, 7)));
        CLPOLY_ASSERT_TRUE(!chk(Zp(3, 7)));
    }

    // ================================================================
    //  Zp: operator= _p==0 时赋零 (回归: 消除除零 UB)
    // ================================================================

    CLPOLY_TEST("Zp_assign_zero_uninitialized");
    {
        Zp a;  // _p==0
        a = (int64_t)0;
        CLPOLY_ASSERT_EQ(a.number(), (uint64_t)0);
        CLPOLY_ASSERT_EQ(a.prime(), (uint32_t)0);

        Zp b;
        b = ZZ(0);
        CLPOLY_ASSERT_EQ(b.number(), (uint64_t)0);
    }

    // ================================================================
    //  Zp: addmul/submul 从 _p==0 累加器开始 (pair_vec_multiplies 场景)
    // ================================================================

    CLPOLY_TEST("Zp_addmul_from_default");
    {
        // 模拟 pair_vec_multiplies 的使用模式:
        // T2 k; set_zero(k); addmul(k, a, b);
        Zp k;
        set_zero(k);
        Zp a(3, 7), b(5, 7);
        addmul(k, a, b);
        CLPOLY_ASSERT_EQ(k.prime(), (uint32_t)7);
        CLPOLY_ASSERT_EQ(k.number(), (uint64_t)(15 % 7));  // 3*5=15, 15%7=1

        // 再累加一次
        Zp c(4, 7), d(2, 7);
        addmul(k, c, d);
        CLPOLY_ASSERT_EQ(k.number(), (uint64_t)((1 + 8) % 7));  // 1+4*2=9, 9%7=2
    }

    CLPOLY_TEST("Zp_submul_from_default");
    {
        Zp k;
        set_zero(k);
        Zp a(3, 7), b(5, 7);
        addmul(k, a, b);  // k=1
        Zp c(2, 7), d(3, 7);
        submul(k, c, d);
        // 1 + 7 - (2*3)%7 = 1 + 7 - 6 = 2
        CLPOLY_ASSERT_EQ(k.number(), (uint64_t)2);
    }

    // ================================================================
    //  Zp: 复合运算、边界值
    // ================================================================

    CLPOLY_TEST("Zp_compound_assignment");
    {
        Zp a(3, 7);
        a += Zp(5, 7);
        CLPOLY_ASSERT_EQ(a.number(), (uint64_t)1);  // (3+5)%7=1

        a *= Zp(4, 7);
        CLPOLY_ASSERT_EQ(a.number(), (uint64_t)4);  // (1*4)%7=4

        a /= Zp(2, 7);
        // 4 * inv(2,7) = 4 * 4 = 16 % 7 = 2
        CLPOLY_ASSERT_EQ(a.number(), (uint64_t)2);
    }

    CLPOLY_TEST("Zp_negation");
    {
        Zp a(3, 7);
        Zp neg = -a;
        CLPOLY_ASSERT_EQ(neg.number(), (uint64_t)4);  // 7-3=4
        Zp sum = a + neg;
        CLPOLY_ASSERT_EQ(sum.number(), (uint64_t)0);

        // -0 == 0
        Zp z(0, 7);
        CLPOLY_ASSERT_EQ((-z).number(), (uint64_t)0);
    }

    CLPOLY_TEST("Zp_comparison_with_int");
    {
        Zp a(3, 7);
        CLPOLY_ASSERT_TRUE(a == 3);
        CLPOLY_ASSERT_TRUE(3 == a);
        CLPOLY_ASSERT_TRUE(a != 4);
        CLPOLY_ASSERT_TRUE(4 != a);
        CLPOLY_ASSERT_TRUE(a == 10);  // 10 % 7 == 3
        CLPOLY_ASSERT_TRUE(a >= 0);
        CLPOLY_ASSERT_TRUE(a >= 3);
        CLPOLY_ASSERT_TRUE(!(a >= 4));
    }

    CLPOLY_TEST("Zp_scalar_mul");
    {
        Zp a(3, 7);
        Zp r1 = a * (int64_t)5;
        CLPOLY_ASSERT_EQ(r1.number(), (uint64_t)1);  // (3*5)%7=1
        Zp r2 = (int64_t)5 * a;
        CLPOLY_ASSERT_EQ(r2.number(), (uint64_t)1);

        // 负标量
        Zp r3 = a * (int64_t)(-1);
        CLPOLY_ASSERT_EQ(r3.number(), (uint64_t)4);  // 3*(-1) -> 3*(7-1)=3*6=18%7=4
    }

    CLPOLY_TEST("Zp_bool_conversion");
    {
        Zp a(0, 7);
        CLPOLY_ASSERT_TRUE(!bool(a));
        Zp b(3, 7);
        CLPOLY_ASSERT_TRUE(bool(b));
    }

    CLPOLY_TEST("Zp_uint64_conversion");
    {
        Zp a(5, 7);
        uint64_t v = (uint64_t)a;
        CLPOLY_ASSERT_EQ(v, (uint64_t)5);
    }

    // ================================================================
    //  Zp 多项式乘法 (回归: pair_vec_multiplies 堆路径)
    // ================================================================

    CLPOLY_TEST("Zp_upoly_mul_heap_path");
    {
        // 两个 ≥2 项的 Zp 多项式相乘，触发堆乘法路径
        // 之前在此路径上 set_zero 导致崩溃
        uint64_t p = 7;
        upolynomial_<Zp> a({
            {umonomial(2), Zp(3, p)},
            {umonomial(1), Zp(2, p)},
            {umonomial(0), Zp(1, p)}
        });  // 3x^2 + 2x + 1
        upolynomial_<Zp> b({
            {umonomial(2), Zp(1, p)},
            {umonomial(1), Zp(4, p)},
            {umonomial(0), Zp(5, p)}
        });  // x^2 + 4x + 5
        auto c = a * b;
        // (3x^2+2x+1)(x^2+4x+5) = 3x^4+14x^3+24x^2+14x+5
        // mod 7: 3x^4 + 0x^3 + 3x^2 + 0x + 5
        CLPOLY_ASSERT_EQ(get_deg(c), (int64_t)4);
        CLPOLY_ASSERT_EQ(c.front().second.number(), (uint64_t)3);   // lc = 3
        CLPOLY_ASSERT_EQ(c.back().second.number(), (uint64_t)5);    // 常数项 = 5
        CLPOLY_ASSERT_EQ(c.size(), (size_t)3);  // 只有 x^4, x^2, x^0 项非零
        // assign 求值验证: a(2)*b(2) == c(2)
        Zp pt(2, p);
        Zp va = assign(a, pt);
        Zp vb = assign(b, pt);
        Zp vc = assign(c, pt);
        CLPOLY_ASSERT_EQ(vc.number(), (va * vb).number());
    }

    // ================================================================
    // 63-bit 素数测试
    // ================================================================

    CLPOLY_TEST("Zp 63-bit prime basic arithmetic");
    {
        uint64_t p = (uint64_t)9223372036854775783ULL;  // 最大 < 2^63 的素数
        uint64_t va = 12345678901234567ULL, vb = 98765432109876543ULL;
        Zp a(va, p);
        Zp b(vb, p);

        // 加法
        Zp c = a + b;
        CLPOLY_ASSERT_EQ(c.number(), a.number() + b.number());  // 和 < p，无需归约

        // 乘法 (与 __int128 直接计算对比)
        Zp d = a * b;
        uint64_t expected = (unsigned __int128)a.number() * b.number() % p;
        CLPOLY_ASSERT_EQ(d.number(), expected);

        // 逆元
        Zp e = a.inv();
        CLPOLY_ASSERT_EQ((a * e).number(), (uint64_t)1);

        // 负值构造
        Zp f(-1, p);
        CLPOLY_ASSERT_EQ(f.number(), p - 1);

        // ZZ 构造
        Zp g(ZZ("123456789012345678901234567890"), p);
        uint64_t expected_g = ZZ("123456789012345678901234567890").fdiv_ui(p);
        CLPOLY_ASSERT_EQ(g.number(), expected_g);
    }

    CLPOLY_TEST("Zp 63-bit prime mul stress");
    {
        uint64_t p = (uint64_t)9223372036854775783ULL;
        // 测试边界值乘法
        Zp a(p - 1, p);  // 最大值
        Zp b(p - 1, p);
        Zp c = a * b;
        uint64_t expected = (unsigned __int128)(p - 1) * (p - 1) % p;
        CLPOLY_ASSERT_EQ(c.number(), expected);

        // 0 和 1 边界
        Zp zero(0, p);
        Zp one(1, p);
        CLPOLY_ASSERT_EQ((a * zero).number(), (uint64_t)0);
        CLPOLY_ASSERT_EQ((a * one).number(), a.number());

        // 多次乘法累积
        Zp acc(1, p);
        uint64_t bval = 123456789ULL;
        Zp base(bval, p);
        for (int i = 0; i < 100; i++)
            acc *= base;
        // 验证 acc * base^{-100} == 1
        Zp base_inv = base.inv();
        Zp dec = acc;
        for (int i = 0; i < 100; i++)
            dec *= base_inv;
        CLPOLY_ASSERT_EQ(dec.number(), (uint64_t)1);
    }

    CLPOLY_TEST("Zp 63-bit prime division and pow");
    {
        uint64_t p = (uint64_t)9223372036854775783ULL;
        Zp a(42, p);
        Zp b(7, p);
        Zp c = a / b;
        CLPOLY_ASSERT_EQ((c * b).number(), a.number());

        // Fermat 小定理：a^(p-1) = 1
        // p-1 太大，用 pow(a, p-2) == a.inv() 代替
        Zp a_inv = a.inv();
        Zp a_inv2 = pow(a, (int64_t)(p - 2));
        CLPOLY_ASSERT_EQ(a_inv.number(), a_inv2.number());
    }

    CLPOLY_TEST("Zp 63-bit prime add/sub boundary");
    {
        uint64_t p = (uint64_t)9223372036854775783ULL;
        Zp a(p - 1, p);
        Zp b(1, p);

        // (p-1) + 1 = 0 mod p
        CLPOLY_ASSERT_EQ((a + b).number(), (uint64_t)0);

        // 0 - 1 = p - 1
        Zp zero(0, p);
        CLPOLY_ASSERT_EQ((zero - b).number(), p - 1);

        // 取负
        CLPOLY_ASSERT_EQ((-a).number(), (uint64_t)1);
        CLPOLY_ASSERT_EQ((-zero).number(), (uint64_t)0);
    }

    return clpoly_test::test_summary();
}
