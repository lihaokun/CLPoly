/**
 * @file clpoly_test.hh
 * @brief CLPoly lightweight test framework - no external dependencies
 */
#ifndef CLPOLY_TEST_HH
#define CLPOLY_TEST_HH

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

namespace clpoly_test {

struct TestResult {
    int passed = 0;
    int failed = 0;
    std::vector<std::string> failures;
};

inline TestResult& global_result() {
    static TestResult r;
    return r;
}

inline std::string& current_test() {
    static std::string name;
    return name;
}

inline std::string& current_section() {
    static std::string name;
    return name;
}

inline std::string location_str(const char* file, int line) {
    std::ostringstream ss;
    ss << file << ":" << line;
    return ss.str();
}

inline void record_pass() {
    global_result().passed++;
}

inline void record_fail(const char* file, int line, const std::string& msg) {
    auto& r = global_result();
    r.failed++;
    std::ostringstream ss;
    ss << "  FAIL [" << current_test();
    if (!current_section().empty())
        ss << "/" << current_section();
    ss << "] " << file << ":" << line << " - " << msg;
    r.failures.push_back(ss.str());
    std::cerr << ss.str() << std::endl;
}

inline int test_summary() {
    auto& r = global_result();
    std::cout << "\n========================================" << std::endl;
    std::cout << "  Total: " << (r.passed + r.failed)
              << "  Passed: " << r.passed
              << "  Failed: " << r.failed << std::endl;
    if (r.failed > 0) {
        std::cout << "  FAILURES:" << std::endl;
        for (auto& f : r.failures)
            std::cout << f << std::endl;
    }
    std::cout << "========================================" << std::endl;
    return r.failed > 0 ? 1 : 0;
}

} // namespace clpoly_test

#define CLPOLY_TEST(name) \
    clpoly_test::current_test() = (name); \
    clpoly_test::current_section().clear(); \
    std::cout << "[ TEST ] " << (name) << std::endl;

#define CLPOLY_TEST_SECTION(name) \
    clpoly_test::current_section() = (name);

#define CLPOLY_ASSERT(cond) \
    do { \
        if (cond) { \
            clpoly_test::record_pass(); \
        } else { \
            clpoly_test::record_fail(__FILE__, __LINE__, \
                std::string("Assertion failed: ") + #cond); \
        } \
    } while(0)

#define CLPOLY_ASSERT_EQ(a, b) \
    do { \
        auto __clpoly_a = (a); \
        auto __clpoly_b = (b); \
        if (__clpoly_a == __clpoly_b) { \
            clpoly_test::record_pass(); \
        } else { \
            std::ostringstream __clpoly_ss; \
            __clpoly_ss << #a << " == " << #b \
                        << "\n    actual:   " << __clpoly_a \
                        << "\n    expected: " << __clpoly_b; \
            clpoly_test::record_fail(__FILE__, __LINE__, __clpoly_ss.str()); \
        } \
    } while(0)

#define CLPOLY_ASSERT_NE(a, b) \
    do { \
        auto __clpoly_a = (a); \
        auto __clpoly_b = (b); \
        if (__clpoly_a != __clpoly_b) { \
            clpoly_test::record_pass(); \
        } else { \
            std::ostringstream __clpoly_ss; \
            __clpoly_ss << #a << " != " << #b \
                        << "\n    both are: " << __clpoly_a; \
            clpoly_test::record_fail(__FILE__, __LINE__, __clpoly_ss.str()); \
        } \
    } while(0)

#define CLPOLY_ASSERT_TRUE(cond)  CLPOLY_ASSERT(cond)
#define CLPOLY_ASSERT_FALSE(cond) CLPOLY_ASSERT(!(cond))

#endif
