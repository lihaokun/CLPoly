/**
 * @file polynomial_type.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义polynomial_type
 * 
 */
#ifndef CLPOLY_POLYNOMIAL_TYPE__HH
#define CLPOLY_POLYNOMIAL_TYPE__HH
#include <clpoly/monomial.hh>
#include <clpoly/number.hh>
#include <clpoly/basic_polynomial.hh>
namespace clpoly{
    template <class Tc,class comp=grlex>
    using polynomial_=basic_polynomial<basic_monomial<comp>,Tc,comp>;
    using polynomial_ZZ=polynomial_<ZZ>;
    using polynomial_QQ=polynomial_<QQ>;
}
#endif