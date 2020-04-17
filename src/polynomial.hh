/*
Module Name:
    atomic_polynomial.hh
Abstract:
    定义类：polynomial
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_POLYNOMIAL_HH
#define CLPOLY_POLYNOMIAL_HH
#include "basic.hh"
#include "variable.hh"
#include "monomial.hh"
#include "number.hh"
#include "basic_polynomial.hh"
#include <list>
#include <string>
namespace clpoly{
    template <class Tc,class comp_v=lex,class comp_m=grlex>
    using polynomial_=basic_polynomial<basic_monomial<comp_v>,Tc,comp_m>;
    using polynomial_ZZ=polynomial_<ZZ>;

    monomial operator* (const variable & v1,const variable & v2)
    {
        if (monomial::comp(v1,v2))  
            return monomial({{v1,1},{v2,1}});
        return monomial({{v2,1},{v1,1}});          
    }
    monomial operator* (const monomial & v1,const variable & v2)
    {
        return v1*monomial(v2);          
    }

    polynomial_ZZ operator* (const ZZ & v1,const monomial & v2)
    {
        return polynomial_ZZ(v2)*v1;          
    }

    polynomial_ZZ operator* (const ZZ & v1,const variable & v2)
    {
        return polynomial_ZZ(v2)*v1;          
    }

    polynomial_ZZ operator* (int64_t v1,const variable & v2)
    {
        return polynomial_ZZ(v2)*ZZ(v1);          
    }
    polynomial_ZZ operator* (const polynomial_ZZ &  v1,const variable & v2)
    {
        return v1*polynomial_ZZ(v2);          
    }

    polynomial_ZZ operator+ (const monomial & m1,const monomial & m2)
    {
        if (polynomial_ZZ::comp(m1,m2))
            return polynomial_ZZ({{m1,1},{m2,1}});
        return polynomial_ZZ({{m2,1},{m1,1}});
    }
    polynomial_ZZ operator+ (const monomial & m,const polynomial_ZZ & p)
    {
        return p+polynomial_ZZ(m);
    }
    polynomial_ZZ operator+ (const polynomial_ZZ & p,const monomial & m)
    {
        return p+polynomial_ZZ(m);
    }

    polynomial_ZZ operator+ (const monomial & m,const ZZ & p)
    {
        return polynomial_ZZ({{m,1},{{},p}});
    }

    polynomial_ZZ operator+ (const ZZ & m,const polynomial_ZZ & p)
    {
        return p+m;
    }

    template<class Tc,class comp_v,class comp_m>
    std::list<std::pair<variable,size_t>> get_variables(const polynomial_<Tc,comp_v,comp_m>& p)
    {
        std::list<std::pair<variable,size_t>> l;
        typename std::list<std::pair<variable,size_t>>::iterator l_ptr;
        typename basic_monomial<comp_v>::const_iterator m_ptr;

        for (const auto & i:p)
        {
            if (!i.first.empty())
            {
                m_ptr=i.first.begin();
                if (!l.empty())
                {
                    l_ptr=l.begin();
                    for(;m_ptr!=i.first.end() && i.first.comp(m_ptr->first,l_ptr->first);l.push_front(*(m_ptr++)));
                    while(l_ptr!=l.end() && m_ptr!=i.first.end())
                    {
                        if (i.first.comp(m_ptr->first,l_ptr->first))
                        {
                            l.insert(l_ptr,*(m_ptr++));
                        }
                        else
                        {
                            if (m_ptr->first==l_ptr->first)//equal_to 
                            {
                                if (m_ptr->second>l_ptr->second) //greater
                                    l_ptr->second=m_ptr->second;
                                ++l_ptr;++m_ptr;
                            }    
                            else
                                ++l_ptr;
                        }
                        
                    } 
                }
                for(;m_ptr!=i.first.end();l.push_back(*(m_ptr++)));  
            }
        }
        return l;
    }
    template <class Tc,class comp_v>
    int64_t get_deg(const polynomial_<Tc,comp_v,graded<basic_monomial<comp_v>>> & p)
    {
        if (p.empty())
            return 0;
        return p.front().first.deg(); 
    }
    template <class Tc,class comp_v,class comp_m>
    int64_t get_deg(const polynomial_<Tc,comp_v,comp_m> & p)
    {
        if (p.empty())
            return 0;
        int64_t deg=p.front().first.deg();
        for(auto &i:p)
            deg=std::max(i.first.deg(),deg);
        return deg;
    }
}
#endif