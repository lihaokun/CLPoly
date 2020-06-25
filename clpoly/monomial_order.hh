/*
Module Name:
    monomial_order.hh
Abstract:
    定义类：单项序和单项压缩算法
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_MONOMIAL_ORDER_HH
#define CLPOLY_MONOMIAL_ORDER_HH
#include <functional>
#include <cassert>
#include "basic_monomial.hh"
namespace clpoly{
    struct lex
    {
        constexpr bool operator()(const variable & v1,const variable & v2) const {return v1<v2;}
        //constexpr bool operator()(uint64_t v1,uint64_t v2) const {return v1<v2;}
        inline bool operator()(const basic_monomial<lex> &m1,const basic_monomial<lex> &m2)const {return pair_vec_comp(m1.data(),m2.data(),m1.comp());}
        constexpr bool operator==(const lex &g1)const {return true;}
        constexpr bool operator!=(const lex &g1)const {return false;}
        
    };
    struct grlex
    {
        constexpr bool operator()(const variable & v1,const variable & v2) const {return v1<v2;}
        //constexpr bool operator()(uint64_t v1,uint64_t v2) const {return v1>v2;}
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
            assert(m1.comp().v==m2.comp().v);
            auto d1=get_up_deg(m1);
            auto d2=get_up_deg(m2);
            return (d1>d2 || (d1==d2 && pair_vec_comp(m1.data(),m2.data(),m1.comp())));
        }
        constexpr bool operator==(const univariate_priority_order &g1)const {return v==g1.v;}
        constexpr bool operator!=(const univariate_priority_order &g1)const {return v!=g1.v;}
        constexpr const variable & var() const{return v;}
    };
    constexpr const int64_t get_up_deg(const basic_monomial<univariate_priority_order>& m1)
    {return (m1.empty() || m1.back().first!=m1.comp().var())?0:m1.back().second;}


/*compression*/

    template<class Tc1,class Tc2,class compare,class compare2>
    constexpr bool __is_monomial_multiplies_compression(
        const std::vector<std::pair<basic_monomial<compare>,Tc1>> & v1_,
        const std::vector<std::pair<basic_monomial<compare>,Tc2>> & v2_,
        const compare2 & comp,
        std::list<variable>& vars
    )
    {
        return false;
    }
    template<class compare>
    constexpr uint64_t _monomial_compression(const basic_monomial<compare> & m,const std::list<variable>& vars)
    {
        return 0;
    }
    template<class compare>
    constexpr void _monomial_decompression(uint64_t mc,basic_monomial<compare> & m,const std::list<variable>& vars,const compare * comp)
    {}

    template<class Tc1,class Tc2>
    bool __is_monomial_multiplies_compression(
        const std::vector<std::pair<basic_monomial<grlex>,Tc1>> & v1_,
        const std::vector<std::pair<basic_monomial<grlex>,Tc2>> & v2_,
        const grlex & comp,
        std::list<variable>& vars
    )
    {
        for (auto & i:v1_)
            for (auto & j:i.first)
                if (j.second<0)
                    return false;
        for (auto & i:v2_)
            for (auto & j:i.first)
                if (j.second<0)
                    return false;
        vars.clear();
        __pair_vec_variables(v1_,vars);
        __pair_vec_variables(v2_,vars);
        if (vars.size()>1 && v1_.begin()->first.deg()+v2_.begin()->first.deg()>=(uint64_t(1)<<(64/(vars.size()+1)))) return false;
        return true;
    }
    template<>
    inline uint64_t _monomial_compression(const basic_monomial<grlex> & m,const std::list<variable>& vars)
    {
        uint64_t mc=0;
        if (m.empty() || vars.empty() )return mc;
        if (vars.size()>1)
        {
            uint l=64/(vars.size()+1);
            mc=m.deg();
            auto m_ptr=m.begin();
            for (auto &i:vars)
            {
                //if (m_ptr==m.end())

                mc<<=l;
                if (m_ptr!=m.end()  && i == m_ptr->first)
                {
                    mc+=m_ptr->second;
                    ++m_ptr;
                }
            }
        }
        else
        {
            mc=m.deg();
        }
        return mc;
    }
    template<>
    inline void _monomial_decompression(uint64_t mc,basic_monomial<grlex> & m,const std::list<variable>& vars,const grlex * comp_ptr)
    {
        m.clear();
        m.comp(comp_ptr);
        if (!mc || vars.empty()) return void();
        m.reserve(vars.size());
        if (vars.size()>1)
        {
            uint l=64/(vars.size()+1);
            uint ll=l*vars.size();
            uint64_t mod=(uint64_t(1)<<(l*(vars.size()+1)))-(uint64_t(1)<<ll);
            m.deg()=(mc & mod)>>ll;
            uint64_t part;
            for (auto &i:vars)
            {
                mc<<=l;
                part=(mc & mod)>>ll;
                if (part!=0)
                    m.data().push_back({i,part});
            }
        }
        else
        {
            m.push_back({*vars.begin(),mc});
        }
    }

    template<class Tc1,class Tc2>
    bool __is_monomial_multiplies_compression(
        const std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc1>> & v1_,
        const std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc2>> & v2_,
        const univariate_priority_order & comp,
        std::list<variable>& vars
    )
    {
        uint64_t deg=0;
        for (auto & i:v1_)
            for (auto & j:i.first)
                if (j.second<0)
                    return false;
                else if(j.second>deg)
                    deg=j.second;

                
        for (auto & i:v2_)
            for (auto & j:i.first)
                if (j.second<0)
                    return false;
                else if(j.second>deg)
                    deg=j.second;
        vars.clear();
        __pair_vec_variables(v1_,vars);
        __pair_vec_variables(v2_,vars);
        if (vars.size()>1  && deg>=(uint64_t(1)<<(64/vars.size()))) return false;
        return true;
    }
    template<>
    inline uint64_t _monomial_compression(const basic_monomial<univariate_priority_order> & m,const std::list<variable>& vars)
    {
        uint64_t mc=0;
        if (m.empty() || vars.empty() )return mc;
        uint l=64/(vars.size());
        variable v=m.comp().v;
        auto m_ptr=m.begin();
        if (vars.back()==v)
        {
            mc=get_up_deg(m);
        }
        for (auto &i:vars)
        {
            if (i!=v)
            {
                mc<<=l;
                if (m_ptr!=m.end()  && i == m_ptr->first)
                {
                    mc+=m_ptr->second;
                    ++m_ptr;
                }
            }
        }
        return mc;
    }
    template<>
    inline void _monomial_decompression(uint64_t mc,basic_monomial<univariate_priority_order> & m,const std::list<variable>& vars,const univariate_priority_order * comp_ptr)
    {
        m.clear();
        m.comp(comp_ptr);
        if (!mc || vars.empty()) return void();
        m.reserve(vars.size());
        variable v=comp_ptr->v;
        uint l=64/(vars.size());
        uint ll=l*(vars.size()-1);
        uint64_t mod=(uint64_t(1)<<(l*(vars.size())))-(uint64_t(1)<<ll);
        uint64_t deg;
        if (vars.back()==v)
        {
            deg=(mc & mod)>>ll;
            mc<<=l;
        }
        uint64_t part;
        for (auto &i:vars)
        {
            if (i!=v)
            {
                part=(mc & mod)>>ll;
                if (part!=0)
                    m.push_back({i,part});
                mc<<=l;
            }
            else
            {
                if (deg !=0)
                    m.push_back({v,deg});
            }   
        }

    }

}
#endif
