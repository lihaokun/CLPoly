// 背靠背测试 C++ 端：调用真实 CLPoly 函数，输出结果
#include <iostream>
#include <cstdint>
#include "../../clpoly/polynomial_factorize_zp.hh"

using namespace clpoly;

// 输出 Zp
void print_zp(const Zp& z) {
    std::cout << z.number() << " " << z.prime();
}

// 输出多项式
void print_poly(const upolynomial_<Zp>& f) {
    std::cout << f.size();
    for (auto& term : f) {
        std::cout << " " << term.first.deg() << " " << term.second.number() ;
    }
}

// ============================================================
// 测试 __make_zp
// ============================================================
void test_make_zp() {
    struct { int64_t val; uint64_t p; } cases[] = {
        {0, 3}, {1, 3}, {2, 3}, {7, 13}, {-1, 5}, {12, 13}, {0, 17}, {100, 17}
    };
    std::cout << "make_zp " << sizeof(cases)/sizeof(cases[0]) << std::endl;
    for (auto& c : cases) {
        auto result = __make_zp(c.val, c.p);
        std::cout << c.val << " " << c.p << " -> ";
        print_zp(result);
        std::cout << std::endl;
    }
}

// ============================================================
// 测试 __extract_pth_root
// ============================================================
void test_extract_pth_root() {
    // f = x^6 + 2x^3 + 1 over F_3 → pth_root = x^2 + 2x + 1
    uint64_t p = 3;
    upolynomial_<Zp> f;
    f.push_back({umonomial(6), Zp(1, p)});
    f.push_back({umonomial(3), Zp(2, p)});
    f.push_back({umonomial(0), Zp(1, p)});

    std::cout << "extract_pth_root 1" << std::endl;
    auto result = __extract_pth_root(f);
    std::cout << "input ";  print_poly(f);  std::cout << std::endl;
    std::cout << "output "; print_poly(result); std::cout << std::endl;
}

// ============================================================
// 测试 __upoly_subtract_x
// ============================================================
void test_subtract_x() {
    uint64_t p = 5;
    // f = x^2 + 3x + 1
    upolynomial_<Zp> f;
    f.push_back({umonomial(2), Zp(1, p)});
    f.push_back({umonomial(1), Zp(3, p)});
    f.push_back({umonomial(0), Zp(1, p)});

    std::cout << "subtract_x 1" << std::endl;
    auto result = __upoly_subtract_x(f, p);
    std::cout << "input ";  print_poly(f);  std::cout << std::endl;
    std::cout << "output "; print_poly(result); std::cout << std::endl;
}

// ============================================================
// 测试 __upoly_make_monic
// ============================================================
void test_make_monic() {
    uint64_t p = 7;
    // f = 3x^2 + 5x + 2
    upolynomial_<Zp> f;
    f.push_back({umonomial(2), Zp(3, p)});
    f.push_back({umonomial(1), Zp(5, p)});
    f.push_back({umonomial(0), Zp(2, p)});

    std::cout << "make_monic 1" << std::endl;
    std::cout << "input "; print_poly(f); std::cout << std::endl;
    __upoly_make_monic(f);
    std::cout << "output "; print_poly(f); std::cout << std::endl;
}

int main() {
    test_make_zp();
    test_extract_pth_root();
    test_subtract_x();
    test_make_monic();
    return 0;
}
