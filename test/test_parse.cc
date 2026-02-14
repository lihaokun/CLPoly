#include <clpoly/parse/basic.hh>
#include "clpoly_test.hh"

using namespace clpoly::parse;

int main() {
    // ================================================================
    //  is_blank
    // ================================================================

    CLPOLY_TEST("is_blank");
    {
        CLPOLY_ASSERT_TRUE(is_blank(' '));
        CLPOLY_ASSERT_TRUE(is_blank('\n'));
        CLPOLY_ASSERT_TRUE(is_blank('\r'));
        CLPOLY_ASSERT_FALSE(is_blank('a'));
        CLPOLY_ASSERT_FALSE(is_blank('0'));
        CLPOLY_ASSERT_FALSE(is_blank('('));
    }

    // ================================================================
    //  is_parentheses
    // ================================================================

    CLPOLY_TEST("is_parentheses");
    {
        CLPOLY_ASSERT_TRUE(is_parentheses('('));
        CLPOLY_ASSERT_TRUE(is_parentheses(')'));
        CLPOLY_ASSERT_FALSE(is_parentheses('['));
        CLPOLY_ASSERT_FALSE(is_parentheses('a'));
        CLPOLY_ASSERT_FALSE(is_parentheses(' '));
    }

    // ================================================================
    //  is_note
    // ================================================================

    CLPOLY_TEST("is_note");
    {
        CLPOLY_ASSERT_TRUE(is_note(';'));
        CLPOLY_ASSERT_FALSE(is_note('#'));
        CLPOLY_ASSERT_FALSE(is_note('a'));
        CLPOLY_ASSERT_FALSE(is_note(' '));
    }

    // ================================================================
    //  is_quoted
    // ================================================================

    CLPOLY_TEST("is_quoted");
    {
        CLPOLY_ASSERT_TRUE(is_quoted('|'));
        CLPOLY_ASSERT_FALSE(is_quoted('"'));
        CLPOLY_ASSERT_FALSE(is_quoted('\''));
        CLPOLY_ASSERT_FALSE(is_quoted('a'));
    }

    // ================================================================
    //  is_number
    // ================================================================

    CLPOLY_TEST("is_number");
    {
        for (char c = '0'; c <= '9'; ++c) {
            CLPOLY_ASSERT_TRUE(is_number(c));
        }
        CLPOLY_ASSERT_FALSE(is_number('a'));
        CLPOLY_ASSERT_FALSE(is_number(' '));
        CLPOLY_ASSERT_FALSE(is_number('/'));
    }

    // ================================================================
    //  is_characters
    // ================================================================

    CLPOLY_TEST("is_characters");
    {
        // All 17 special characters
        const char specials[] = {
            '~', '!', '@', '$', '%', '^', '&', '*',
            '_', '-', '+', '=', '<', '>', '.', '?', '/'
        };
        for (char c : specials) {
            CLPOLY_ASSERT_TRUE(is_characters(c));
        }
        // Letters and digits should return false
        CLPOLY_ASSERT_FALSE(is_characters('a'));
        CLPOLY_ASSERT_FALSE(is_characters('Z'));
        CLPOLY_ASSERT_FALSE(is_characters('0'));
        CLPOLY_ASSERT_FALSE(is_characters('9'));
        CLPOLY_ASSERT_FALSE(is_characters(' '));
        CLPOLY_ASSERT_FALSE(is_characters('('));
        CLPOLY_ASSERT_FALSE(is_characters(';'));
    }

    return clpoly_test::test_summary();
}
