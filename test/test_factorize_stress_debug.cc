/**
 * @file test_factorize_stress_debug.cc
 * @brief 复现 nfactor_fail 的详细调试
 */
#include <clpoly/clpoly.hh>
#include <clpoly/polynomial_factorize.hh>
#include "crosscheck_flint.hh"
#include <iostream>
#include <iomanip>

using namespace clpoly;

void check_upoly(const std::string& label, const upolynomial_<ZZ>& f)
{
    std::cout << "\n--- " << label << " ---" << std::endl;
    std::cout << "f = " << f << std::endl;
    std::cout << "deg = " << get_deg(f) << ", terms = " << f.size() << std::endl;

    auto cl_fac = factorize(f);
    auto fl_fac = crosscheck::flint_factor_upoly(f);

    std::cout << "CLPoly: content=" << cl_fac.content
              << ", " << cl_fac.factors.size() << " factors:" << std::endl;
    for (auto& [fi, ei] : cl_fac.factors)
        std::cout << "  (" << fi << ")^" << ei << std::endl;

    std::cout << "FLINT:  content=" << fl_fac.content
              << ", " << fl_fac.factors.size() << " factors:" << std::endl;
    for (auto& [fi, ei] : fl_fac.factors)
        std::cout << "  (" << fi << ")^" << ei << std::endl;

    // 验证乘积
    upolynomial_<ZZ> cl_prod;
    cl_prod.push_back(std::make_pair(umonomial(0), cl_fac.content));
    for (auto& [fi, ei] : cl_fac.factors) {
        auto fi_pow = pow(fi, (int)ei);
        pair_vec_multiplies(cl_prod.data(), cl_prod.data(), fi_pow.data(), cl_prod.comp());
        cl_prod.normalization();
    }
    std::cout << "CLPoly product " << (cl_prod == f ? "OK" : "MISMATCH") << std::endl;
}

int main()
{
    // 复现三类失败: uni C=100 3fac, uni C=100000 3fac, uni deg8 3fac
    // 用固定种子尝试找到失败用例

    std::cout << "==== Reproducing nfactor failures ====" << std::endl;

    // 扫描 uni C=100 3fac deg3
    std::cout << "\n== uni C=100 3fac deg3 ==" << std::endl;
    for (int trial = 0; trial < 20; ++trial) {
        auto f1 = random_upolynomial<ZZ>(3, 3, {-100, 100});
        auto f2 = random_upolynomial<ZZ>(3, 3, {-100, 100});
        auto f3 = random_upolynomial<ZZ>(3, 3, {-100, 100});
        if (f1.empty() || f2.empty() || f3.empty()) continue;
        upolynomial_<ZZ> f;
        f.push_back(std::make_pair(umonomial(0), ZZ(1)));
        pair_vec_multiplies(f.data(), f.data(), f1.data(), f.comp()); f.normalization();
        pair_vec_multiplies(f.data(), f.data(), f2.data(), f.comp()); f.normalization();
        pair_vec_multiplies(f.data(), f.data(), f3.data(), f.comp()); f.normalization();
        if (f.empty() || get_deg(f) < 2) continue;

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_upoly(f);
        if (cl_fac.factors.size() != fl_fac.factors.size()) {
            check_upoly("trial " + std::to_string(trial), f);
        }
    }

    // 扫描 uni deg8 3fac C=10
    std::cout << "\n== uni deg8 3fac C=10 ==" << std::endl;
    for (int trial = 0; trial < 15; ++trial) {
        auto f1 = random_upolynomial<ZZ>(8, 8, {-10, 10});
        auto f2 = random_upolynomial<ZZ>(8, 8, {-10, 10});
        auto f3 = random_upolynomial<ZZ>(8, 8, {-10, 10});
        if (f1.empty() || f2.empty() || f3.empty()) continue;
        upolynomial_<ZZ> f;
        f.push_back(std::make_pair(umonomial(0), ZZ(1)));
        pair_vec_multiplies(f.data(), f.data(), f1.data(), f.comp()); f.normalization();
        pair_vec_multiplies(f.data(), f.data(), f2.data(), f.comp()); f.normalization();
        pair_vec_multiplies(f.data(), f.data(), f3.data(), f.comp()); f.normalization();
        if (f.empty() || get_deg(f) < 2) continue;

        auto cl_fac = factorize(f);
        auto fl_fac = crosscheck::flint_factor_upoly(f);
        if (cl_fac.factors.size() != fl_fac.factors.size()) {
            check_upoly("trial " + std::to_string(trial), f);
        }
    }

    return 0;
}
