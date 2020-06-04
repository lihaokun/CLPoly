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
#include <random>

namespace clpoly{
    template <class Tc,class comp=grlex>
    using polynomial_=basic_polynomial<basic_monomial<comp>,Tc,comp>;
    using polynomial_ZZ=polynomial_<ZZ>;
    using polynomial_QQ=polynomial_<QQ>;
    
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
    
    template<class Tc>
    polynomial_<Tc> operator+ (const monomial & m,const polynomial_<Tc> & p)
    {
        return p+polynomial_<Tc>(m);
    }

    template<class Tc>
    polynomial_<Tc> operator+ (const polynomial_<Tc> & p,const monomial & m)
    {
        return p+polynomial_<Tc>(m);
    }

    polynomial_ZZ operator+ (const monomial & m,const ZZ & p)
    {
        return polynomial_ZZ({{m,1},{{},p}});
    }

    template<class Tc>
    polynomial_<Tc> operator+ (const ZZ & m,const polynomial_<Tc> & p)
    {
        return p+m;
    }
    template<class Tc>
    polynomial_<Tc> random_polynomial(const std::vector<variable> & v,uint64_t deg,double p,int up,int down)
    {
        std::vector<std::pair<monomial,Tc>> p1;
        std::random_device rd; 
        std::mt19937 gen(rd());
        auto size=v.size();
        std::bernoulli_distribution mp(p);
        std::uniform_int_distribution<> dis(down, up);
        std::vector<std::pair<variable,int64_t>>  m;
        for (auto &i:v)
            m.emplace_back(i,0);
        for(int d=deg;d>=0;--d)
        {
            
            if (d)
            {
                int64_t i0=0;
                int64_t sum=d;
                m.begin()->second=d;
                while(1)
                {
                    if (mp(gen))
                    {
                        p1.emplace_back(monomial(m),Tc(dis(gen)));
                    }
                    while (i0>=0)
                    {
                        if (m[i0].second && i0<size-1)
                        {
                            --m[i0].second;
                            ++i0;
                            m[i0].second=d-sum+1;
                            sum=d;
                            break;
                        }
                        else
                        {
                            sum-=m[i0].second;
                            m[i0].second=0;
                            --i0;
                        }
                    }   
                    if (i0<0)   
                        break;

                }
            }
            else
            {
                p1.push_back({{},Tc(dis(gen))});
            }
            
        }
        return polynomial_<Tc>(p1);
    }
    template<class Tc,class comp>
    std::list<std::pair<variable,int64_t>> get_variables(const polynomial_<Tc,comp>& p)
    {
        std::list<std::pair<variable,int64_t>> l;
        typename std::list<std::pair<variable,int64_t>>::iterator l_ptr;
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
        prem(O1,G1,F1);
        polynomial_<Tc,comp> O(G.comp_ptr());
        poly_convert(std::move(O1),O);
        return O;
    }
    template <class Tc,class compare>
    void __onestep__prem(
        const std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & G,
        const std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & F0,
        int64_t f_deg,
        const std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & F,
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & O,
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & F1,
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & G1,
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> & G2,
        compare comp 
    )
    {

        auto G_ptr=G.begin();
        auto G_end=G.end();
        int64_t g_deg=get_first_deg(G_ptr->first);
        G1.clear();
        G1.reserve(G.size());
        while(G_ptr!=G_end && get_first_deg(G_ptr->first)==g_deg)
        {
            G1.push_back(*(G_ptr++));
            G1.back().first.begin()->second=g_deg-f_deg;
        }
        if(G_ptr==G_end)
        {
            return void();
        }
        pair_vec_multiplies(F1,G1,F,comp);
        auto F_ptr=F1.begin();
        auto F_end=F1.end();
        G2.clear();
        G2.reserve((G_end-G_ptr));
        G1.clear();
        while(G_ptr!=G_end)
            G2.push_back(*(G_ptr++));
        pair_vec_multiplies(G1,G2,F0,comp);
        pair_vec_sub(O,G1,F1,comp);
    }
    template <class Tc>
    void prem(
            polynomial_<Tc,univariate_first_order>&O,
            const polynomial_<Tc,univariate_first_order>&G,
            const polynomial_<Tc,univariate_first_order>&F
            )
    {
        O.clear();
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
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_G1;
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_G2;
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_F0;
        tmp_F0.reserve(F.size());
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_F_;
        tmp_F_.reserve(F.size());
        for (auto & i:F)
            if (get_first_deg(i.first,var)!=f_deg)
                tmp_F_.push_back(i);
            else
            {
                tmp_F0.push_back(i);
                tmp_F0.back().first.begin()->second=0;
            }
        std::vector<std::pair<basic_monomial<univariate_first_order>,Tc>> tmp_F1;
        O.data().reserve(G.size());
        __onestep__prem(G.data(),tmp_F0,f_deg,tmp_F_,O.data(),tmp_F1,tmp_G1,tmp_G2,G.comp());
        //std::cout<<"prem_:"<<O<<std::endl;
        while (get_first_deg(O.begin()->first,var)>=f_deg)
        {
            tmp_G=std::move(O.data());
            __onestep__prem(tmp_G,tmp_F0,f_deg,tmp_F_,O.data(),tmp_F1,tmp_G1,tmp_G2,G.comp());    
            //std::cout<<"prem_:"<<O<<std::endl;
        }
    }
    
    template <class Tc,class comp>
    inline polynomial_<Tc,comp> resultant(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v)
    {
        assert(G.comp_ptr()==F.comp_ptr());
        univariate_first_order comp_v(v);
        polynomial_<Tc,univariate_first_order>  G1(&comp_v);
        polynomial_<Tc,univariate_first_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        polynomial_<Tc,univariate_first_order>  O1(&comp_v);
        resultant(O1,G1,F1);
        polynomial_<Tc,comp> O(G.comp_ptr());
        poly_convert(std::move(O1),O);
        return O;
    }
    template <class Tc>
    void resultant
        (   
            polynomial_<Tc,univariate_first_order>&O,
            const polynomial_<Tc,univariate_first_order>&F,
            const polynomial_<Tc,univariate_first_order>&G
            
        )
    {
        O.clear();
        const variable & var=G.comp().var();
        int64_t l;
        int64_t m;
        if (G.empty()||  F.empty() )
            return void();
        if ( (m=get_first_deg(F.begin()->first,var))<(l=get_first_deg(G.begin()->first,var)))
        {
            resultant(O1,G1,F1);
            return void();
        }
        

        
    }
}
#endif