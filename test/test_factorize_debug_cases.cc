/**
 * @file test_factorize_debug_cases.cc
 * @brief 手工验证 stress test 中的 cross_fail 用例
 */
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include "crosscheck_flint.hh"
#include <iostream>

using namespace clpoly;

void check(const std::string& label, const upolynomial_<ZZ>& f)
{
    std::cout << "\n=== " << label << " ===" << std::endl;
    std::cout << "f = " << f << std::endl;
    std::cout << "deg=" << get_deg(f) << " terms=" << f.size() << std::endl;

    auto cl = factorize(f);
    auto fl = crosscheck::flint_factor_upoly(f);

    std::cout << "\nCLPoly: content=" << cl.content << ", " << cl.factors.size() << " factors:" << std::endl;
    for (auto& [fi, ei] : cl.factors)
        std::cout << "  [e=" << ei << "] " << fi << std::endl;

    std::cout << "\nFLINT:  content=" << fl.content << ", " << fl.factors.size() << " factors:" << std::endl;
    for (auto& [fi, ei] : fl.factors)
        std::cout << "  [e=" << ei << "] " << fi << std::endl;

    // 验证乘积
    upolynomial_<ZZ> cl_prod;
    cl_prod.push_back(std::make_pair(umonomial(0), cl.content));
    for (auto& [fi, ei] : cl.factors) {
        auto fi_pow = pow(fi, (int)ei);
        pair_vec_multiplies(cl_prod.data(), cl_prod.data(), fi_pow.data(), cl_prod.comp());
        cl_prod.normalization();
    }
    std::cout << "\nCLPoly product == f? " << (cl_prod == f ? "YES" : "NO") << std::endl;
}

int main()
{
    // 用例 1: uni C=100 4fac deg3 trial 18
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
        check("Case 1: C=100 4fac deg3", f);
    }

    // 用例 2: uni C=10000 4fac deg3 trial 9
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
        check("Case 2: C=10000 4fac deg3", f);
    }

    return 0;
}
