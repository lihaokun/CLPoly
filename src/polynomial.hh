/*
Module Name:
    atomic_polynomial.hh
Abstract:
    定义polynomial
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
    template <class Tc,class comp=grlex>
    using polynomial_=basic_polynomial<basic_monomial<comp>,Tc,comp>;
    using polynomial_ZZ=polynomial_<ZZ>;

    monomial operator* (const variable & v1,const variable & v2)
    {
        if (monomial::compare_type()(v1,v2))  
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
        return polynomial_ZZ(monomial(v2))*v1;          
    }

    polynomial_ZZ operator* (int64_t v1,const variable & v2)
    {
        return polynomial_ZZ(monomial(v2))*ZZ(v1);          
    }
    
    polynomial_ZZ operator* (const polynomial_ZZ &  v1,const variable & v2)
    {
        return v1*polynomial_ZZ(monomial(v2));          
    }

    polynomial_ZZ operator+ (const monomial & m1,const monomial & m2)
    {
        if (polynomial_ZZ::compare_type()(m1,m2))
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

    template<class Tc,class comp>
    std::list<std::pair<variable,size_t>> get_variables(const polynomial_<Tc,comp>& p)
    {
        std::list<std::pair<variable,size_t>> l;
        typename std::list<std::pair<variable,size_t>>::iterator l_ptr;
        typename basic_monomial<comp>::const_iterator m_ptr;

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

    template <class Tc>
    int64_t get_deg(const polynomial_<Tc> & p)
    {
        if (p.empty())
            return 0;
        return p.front().first.deg(); 
    }

    template <class Tc,class comp>
    int64_t get_deg(const polynomial_<Tc,comp> & p)
    {
        if (p.empty())
            return 0;
        int64_t deg=p.front().first.deg();
        for(auto &i:p)
            deg=std::max(i.first.deg(),deg);
        return deg;
    }

    template<class T1,class T2,class comp1,class comp2>
    void poly_convert(const polynomial_<T1,comp1>& p1,polynomial_<T2,comp2> & p2)
    {
        p2.clear();
        basic_monomial<comp2> m(p2.comp_ptr());
        for (auto &i:p1)
        {
            m=i.first.data();
            p2.push_back({std::move(m),i.second});
        }
        p2.normalization();
    }
    template<class T1,class T2,class comp1,class comp2>
    void poly_convert(polynomial_<T1,comp1>&& p1,polynomial_<T2,comp2> & p2)
    {
        p2.clear();
        basic_monomial<comp2> m(p2.comp_ptr());
        for (auto &i:p1)
        {
            m=std::move(i.first.data());
            p2.push_back({std::move(m),std::move(i.second)});
        }
        p2.normalization();
        p1.clear();
    }

    template <class Tc,class comp>
    inline polynomial_<Tc,comp> prem(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v)
    {
        assert(G.comp_ptr()!=F.comp_ptr());
        univariate_first_order comp_v(v);
        polynomial_<Tc,univariate_first_order>  G1(&comp_v);
        polynomial_<Tc,univariate_first_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        polynomial_<Tc,univariate_first_order>  O1(&comp_v);
        prem(G1,F1,O1);
        polynomial_<Tc,comp> O(G.comp_ptr());
        poly_convert(std::move(O1),O);
        return O;
    }
    // template <class Tc>
    // void __onestep__prem(
    //     std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & G,
    //     std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & F,
    //     std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & O
    // )
    // {
        
    // }
    // template <class Tc>
    // void prem(const polynomial_<Tc,univariate_first_order>&G,
    //           const polynomial_<Tc,univariate_first_order>&F,
    //           polynomial_<Tc,univariate_first_order>&O
    //         )
    // {
    //     O.clear();
    //     if (F.empty() || get_first_deg(F.begin()->first)==0)
    //         return void();
        
    //     const variable & var=G.comp().variable;
    //     std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_G0;
    //     //const std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> &tmp_G=G.data();
    //     std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_F0;
    //     std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_F_;
    //     std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_S;
    //     std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_M;

    //     while ()
    //     {
            
    //     }
        
    // }
    

    // template <std::size_t v,class Tc>
    // polynomial_<Tc,vfvcomp<v>,vfmcomp<v>> prem(const polynomial_<Tc,vfvcomp<v>,vfmcomp<v>> &G,const polynomial_<Tc,vfvcomp<v>,vfmcomp<v>> & F)
    // {
    //     static variable var=get_variables(v);
    //     std::vector<std::pair<basic_monomial<vfvcomp<v>>,Tc>> tmp_G0;
    //     std::vector<std::pair<basic_monomial<vfvcomp<v>>,Tc>> tmp_G=G.data();
    //     std::vector<std::pair<basic_monomial<vfvcomp<v>>,Tc>> tmp_F0;
    //     std::vector<std::pair<basic_monomial<vfvcomp<v>>,Tc>> tmp_F_;
    //     std::vector<std::pair<basic_monomial<vfvcomp<v>>,Tc>> tmp_S;
    //     std::vector<std::pair<basic_monomial<vfvcomp<v>>,Tc>> tmp_M;
    // }
}
#endif