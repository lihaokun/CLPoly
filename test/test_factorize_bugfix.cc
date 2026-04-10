/**
 * @file test_factorize_bugfix.cc
 * @brief 验证 lc-baking 首一归一化修复的两个复现用例
 *
 * 修复前：CLPoly 返回 3 个因子（两个不可约三次被合并为一个六次）
 * 修复后：应返回 4 个不可约三次因子
 */
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include <iostream>
#include <cassert>

using namespace clpoly;

// 验证因式分解：乘积还原 + 无 degree >= 4 的因子
// Bug 表现：两个不可约三次被合并为一个六次因子
// 修复后：所有因子 degree <= 3
bool verify(const std::string& label, const upolynomial_<ZZ>& f,
            int max_factor_deg)
{
    std::cout << "=== " << label << " ===" << std::endl;
    std::cout << "  deg=" << get_deg(f) << " terms=" << f.size() << std::endl;

    auto result = factorize(f);

    int max_deg = 0;
    std::cout << "  content=" << result.content
              << ", factors=" << result.factors.size() << std::endl;
    for (auto& [fi, ei] : result.factors) {
        int d = get_deg(fi);
        if (d > max_deg) max_deg = d;
        std::cout << "    [e=" << ei << " deg=" << d << "] " << fi << std::endl;
    }

    // 验证乘积
    upolynomial_<ZZ> prod;
    prod.push_back(std::make_pair(umonomial(0), result.content));
    for (auto& [fi, ei] : result.factors) {
        auto fi_pow = pow(fi, (int)ei);
        pair_vec_multiplies(prod.data(), prod.data(), fi_pow.data(), prod.comp());
        prod.normalization();
    }
    bool product_ok = (prod == f);
    std::cout << "  product == f? " << (product_ok ? "YES" : "NO") << std::endl;

    bool deg_ok = (max_deg <= max_factor_deg);
    std::cout << "  max factor deg=" << max_deg << " <= " << max_factor_deg << "? "
              << (deg_ok ? "YES" : "NO (BUG: incomplete factorization)") << std::endl;

    bool ok = product_ok && deg_ok;
    std::cout << (ok ? "  PASS" : "  *** FAIL ***") << std::endl;
    return ok;
}

int main()
{
    int passed = 0, total = 0;

    // Case 1: 4 个 deg3 因子，系数范围 C=100
    // 修复前: 3 因子（1 个 deg6 + 2 个 deg3）
    // 修复后: 4 个 deg3
    {
        upolynomial_<ZZ> f;
        f.push_back({umonomial(12), ZZ("10707840")});
        f.push_back({umonomial(11), ZZ("-4358016")});
        f.push_back({umonomial(10), ZZ("-17120232")});
        f.push_back({umonomial(9), ZZ("-1940088")});
        f.push_back({umonomial(8), ZZ("3144483")});
        f.push_back({umonomial(7), ZZ("-20309937")});
        f.push_back({umonomial(6), ZZ("-15166461")});
        f.push_back({umonomial(5), ZZ("11473425")});
        f.push_back({umonomial(4), ZZ("-2040858")});
        f.push_back({umonomial(3), ZZ("-7493748")});
        f.push_back({umonomial(2), ZZ("6025272")});
        f.push_back({umonomial(1), ZZ("4287360")});
        ++total;
        if (verify("Case 1: C=100 4fac deg3", f, 3)) ++passed;
    }

    // Case 2: 4 个 deg3 因子，系数范围 C=10000
    {
        upolynomial_<ZZ> f;
        f.push_back({umonomial(12), ZZ("226923284090880")});
        f.push_back({umonomial(11), ZZ("350116202211840")});
        f.push_back({umonomial(10), ZZ("-1455995722245888")});
        f.push_back({umonomial(9), ZZ("-1949032086005432")});
        f.push_back({umonomial(8), ZZ("768244069014269")});
        f.push_back({umonomial(7), ZZ("2086218141435743")});
        f.push_back({umonomial(6), ZZ("2623576334221089")});
        f.push_back({umonomial(5), ZZ("-1360439581386935")});
        f.push_back({umonomial(4), ZZ("-731715361913790")});
        f.push_back({umonomial(3), ZZ("-922193000618706")});
        f.push_back({umonomial(2), ZZ("492400886704650")});
        ++total;
        if (verify("Case 2: C=10000 4fac deg3", f, 3)) ++passed;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Total: " << total << "  Passed: " << passed
              << "  Failed: " << (total - passed) << std::endl;
    std::cout << "========================================" << std::endl;

    return (passed == total) ? 0 : 1;
}
