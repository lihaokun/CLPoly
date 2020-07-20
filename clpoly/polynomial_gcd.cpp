/*
Module Name:
    polynomial_gcd.hh
Abstract:
    定义polynomial_gcd
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_POLYNOMIAL_GCD_HH
#define CLPOLY_POLYNOMIAL_GCD_HH
#include <clpoly/polynomial.hh>
#include <boost/math/special_functions/prime.hpp>

namespace clpoly{    
    
    template <class Tc,class comp>
    inline polynomial_<Tc,comp> polynomial_GCD(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F)
    {
        assert(G.comp_ptr()==G.comp_ptr() || G.comp()==F.comp());
        polynomial_<Tc,comp> G_lc(G.comp_ptr());
        auto vars=G.variables();
        auto F_vars=F.variables();
        _variables_pair_marge(F_vars.begin(),F_vars.end(),vars,G.comp());
        
    }
}
#endif