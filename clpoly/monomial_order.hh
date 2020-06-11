/*
Module Name:
    basic_monomial.hh
Abstract:
    定义类：monomial
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_MONOMIAL_ORDER_HH
#define CLPOLY_MONOMIAL_ORDER_HH
#include <functional>
#include "basic_monomial.hh"
namespace clpoly{
    struct grlex
    {
        constexpr bool operator()(const variable & v1,const variable & v2) const {return v1<v2;}
        inline bool operator()(const basic_monomial<grlex> &m1,const basic_monomial<grlex> &m2)const {return m1>m2;}
        constexpr bool operator==(const grlex &g1)const {return true;}
        constexpr bool operator!=(const grlex &g1)const {return false;}
    };
    struct grevlex
    {
        constexpr bool operator()(const variable & v1,const variable & v2) const {return v1>v2;}
        inline bool operator()(const basic_monomial<grevlex> &m1,const basic_monomial<grevlex> &m2)const {return m1>m2;}
        constexpr bool operator==(const grevlex &g1)const {return true;}
        constexpr bool operator!=(const grevlex &g1)const {return false;}
    };
    // template <class comp>
    // constexpr const int64_t get_up_deg(const basic_monomial<comp>& m1){return m1.empty()?0:m1.back().second;}
    struct univariate_priority_order;
    constexpr const int64_t get_up_deg(const basic_monomial<univariate_priority_order>& m1);
    struct univariate_priority_order
    {
        variable v;
        univariate_priority_order():v(){}
        univariate_priority_order(variable _v):v(_v){}
        constexpr bool operator()(const variable & v1,const variable & v2) const {return(v1!=v && (v2==v || v1<v2));}
        inline bool operator()(const basic_monomial<univariate_priority_order> &m1,const basic_monomial<univariate_priority_order> &m2)const 
        {
            auto d1=get_up_deg(m1);
            auto d2=get_up_deg(m2);
            return (d1>d2 || (d1==d2 && m1>m2));
        }
        constexpr bool operator==(const univariate_priority_order &g1)const {return v==g1.v;}
        constexpr bool operator!=(const univariate_priority_order &g1)const {return v!=g1.v;}
        constexpr const variable & var() const{return v;}
    };
    constexpr const int64_t get_up_deg(const basic_monomial<univariate_priority_order>& m1)
    {return (m1.empty() || m1.back().first!=m1.comp().var())?0:m1.back().second;}
    
}
#endif
