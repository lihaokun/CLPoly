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

//#define x**y pow(x,y)
namespace clpoly{
    template <class Tc,class comp=grlex>
    using polynomial_=basic_polynomial<basic_monomial<comp>,Tc,comp>;
    using polynomial_ZZ=polynomial_<ZZ>;

    // template <class T>
    // monomial pow(const basic_monomial<T> & m,uint64_t n)
    // {
    //     return m.power(n);
    // }
    // template <class T1,class T2,class T3>
    // monomial pow(const basic_polynomial<T1,T2,T3> & m,uint64_t n)
    // {
    //     return m.power(n);
    // }
    
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
        assert(G.comp_ptr()==F.comp_ptr());
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
    template <class Tc,class compare>
    void __onestep__prem(
        const std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & G,
        int64_t f_deg,
        const std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & F,
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & O,
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & F1,
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & G1,
        compare comp 

    )
    {

        auto G_ptr=G.begin();
        auto G_end=G.end();
        int64_t g_deg=get_first_deg(G_ptr->first);
        G1.clear();
        G1.reserve(G.size());
        while(G_ptr!=G.end() && get_first_deg(G_ptr->first)==g_deg)
        {
            G1.push_back(*(G_ptr++));
            G1.back().first.begin()->second+=g_deg-f_deg;
        }
        O.clear();
        if(G_ptr==G.end())
        {
            return void();
        }
        pair_vec_multiplies(F1,G1,F,comp);
        auto F_ptr=F1.begin();
        auto F_end=F1.end();
        O.reserve((G_end-G_ptr)+F1.size());
        Tc tmp;
        while (G_ptr!=G.end() && F_ptr!=F_end)
        {
            if (comp(G_ptr->first,F_ptr->first))
            {
                O.push_back(*(G_ptr++));
            }
            else{
                if ( G_ptr->first==F_ptr->first) //equal_to
                {
                    if (!zore_check<Tc>()(tmp=G_ptr->second-F_ptr->second)) //minus
                        O.emplace_back(F_ptr->first,tmp);
                    ++G_ptr;++F_ptr;
                }
                else
                {
                    O.emplace_back(F_ptr->first,-F_ptr->second); //negate
                    ++F_ptr;
                }
            }
        }
    
        while (F_ptr!=F_end)
        {
            O.emplace_back(F_ptr->first,-F_ptr->second);//negate
            ++F_ptr;
        }
        while(G_ptr!=G.end())
        {
            O.push_back(*(G_ptr++));
        }  
    }
    template <class Tc>
    void prem(const polynomial_<Tc,univariate_first_order>&G,
              const polynomial_<Tc,univariate_first_order>&F,
              polynomial_<Tc,univariate_first_order>&O
            )
    {
        O.data().clear();
        const variable & var=G.comp().var();
        int64_t f_deg;
        int64_t g_deg;
        if (G.empty()||  F.empty() ||
            (f_deg=get_first_deg(F.begin()->first,var))==0 || 
            (g_deg=get_first_deg(G.begin()->first,var))==0)
            return void();
        if (g_deg<f_deg)
        {
            O=G;
            return void();
        }
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_G;
        //tmp_G.reserve(G.size());
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_G1;
        tmp_G1.reserve(G.size());
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_F_;
        tmp_F_.reserve(F.size());
        for (auto & i:F)
            if (get_first_deg(i.first,var)!=f_deg)
                tmp_F_.push_back(i);
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_F1;
        O.data().reserve(G.size());
        __onestep__prem(G.data(),f_deg,tmp_F_,O.data(),tmp_F1,tmp_G1,G.comp());
        while (get_first_deg(O.begin()->first,var)>=f_deg)
        {
            std::cout<<"prem_:"<<O<<std::endl;
            tmp_G=std::move(O.data());
            __onestep__prem(tmp_G,f_deg,tmp_F_,O.data(),tmp_F1,tmp_G1,G.comp());    
        }
    }
}
#endif