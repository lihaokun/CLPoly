/**
 * @file polynomial.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义polynomial和相关的运算
*/
#ifndef CLPOLY_POLYNOMIAL_HH
#define CLPOLY_POLYNOMIAL_HH

#include <clpoly/polynomial_.hh>
#include <list>
#include <map>
#include <string>

namespace clpoly{

    
/*======================================主要函数=========================================*/
    // template <class Tc,class comp>
    // int64_t degree(const polynomial_<Tc,comp> & p);
    // template <class Tc,class comp>
    // int64_t degree(const polynomial_<Tc,comp> & p,const variable & v);
    // template<class Tc>
    // template<class T1,class T2,class comp1,class comp2>
    // void poly_convert(const polynomial_<T1,comp1>& p_in,polynomial_<T2,comp2> & p_out);
    // template <class Tc,class comp>
    // inline polynomial_<Tc,comp> leadcoeff(const polynomial_<Tc,comp>&F);
    // template <class Tc,class comp>
    // inline polynomial_<Tc,comp> prem(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v);
    // template <class Tc,class comp>
    // inline polynomial_<Tc,comp> resultant(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v);
    // template <class Tc, class To>
    // polynomial_<Tc,To> assign(const polynomial_<Tc,To>& p,const variable & v,const Tc & c);
    // template <class compare>
    // polynomial_<Zp,compare> polynomial_mod(polynomial_<ZZ,compare> && p, uint32_t prime)
/*======================================实现===============================================*/
   
 
    

    template <class Tc,class comp>
    int64_t degree(const polynomial_<Tc,comp> & p)
    {
        return p.degree();
    }
    template <class Tc,class comp>
    int64_t degree(const polynomial_<Tc,comp> & p,const variable & v)
    {
        // auto & l=p.variables();
        // for (auto &i:l)
        //     if (l.first==v)
        //         return l.second;
        // return 0;
        if (p.empty())
            return 0;
        auto ptr=p.begin();
        int64_t deg=(ptr++)->first.deg(v),tmp;
        for (;ptr!=p.end();++ptr)
            if ((tmp=ptr->first.deg(v))>deg)
                deg=tmp;
        return deg;
    }
    
    template <class Tc>
    constexpr int64_t get_up_deg(const polynomial_<Tc,univariate_priority_order>& p)
    {
        return p.empty()?0:get_up_deg(p.begin()->first);
    }
    template <class var_order>
    constexpr int64_t get_first_deg(const basic_monomial<lex_<var_order>>& m)
    {
        return m.empty()?0:m.front().second;
    }
    template <class var_order>
    constexpr variable get_first_var(const basic_monomial<lex_<var_order>>& m)
    {
        return m.empty()?variable():m.front().first;
    }
    template <class Tc,class var_order>
    constexpr int64_t get_first_deg(const polynomial_<Tc,lex_<var_order>>& p)
    {
        return p.empty()?0:(p.front().first.empty()?0:p.front().first.front().second);
    }
    template <class Tc,class var_order>
    constexpr variable get_first_var(const polynomial_<Tc,lex_<var_order>>& p)
    {
        return p.empty()?variable():(p.front().first.empty()?variable():p.front().first.front().first);
    }

    template <class Tc,class var_order>
    std::pair<variable,int64_t> get_last_var_deg(const polynomial_<Tc,lex_<var_order>>& p)
    {
        variable v=get_first_var(p);
        int64_t d=get_first_deg(p);
        for (auto &i:p)
        {
            if (!i.first.empty())
            {
                if (p.comp(v,i.first.back().first))
                {
                    v=i.first.back().first;
                    d=i.first.back().second;
                }
                else if (v==i.first.back().first && d<i.first.back().second)
                {
                    d=i.first.back().second;
                }
            }
        }
        return {v,d};
    }
    template <class Tc,class var_order>
    int64_t degree(const polynomial_<Tc,lex_<var_order>> & p,const variable & v)
    {
        if (p.empty())
            return 0;
        if (get_first_var(p)==v)
            return get_first_deg(p);
        if (p.comp()(v,get_first_var(p)))
            return 0;
        auto ptr=p.begin();
        int64_t deg=(ptr++)->first.deg(v),tmp;
        for (;ptr!=p.end();++ptr)
            if ((tmp=ptr->first.deg(v))>deg)
                deg=tmp;
        return deg;
    }

    template<class comp>
    basic_monomial<comp> __change_monomial_var_deg(const basic_monomial<comp> &m,variable v,uint64_t d)
    {
        basic_monomial<comp>  mout(m.comp_ptr());
        auto ptr=m.begin();
        bool s=true;
        for (;ptr!=m.end();++ptr)
        {
            if (v==ptr->first)
            {
                if (d!=0)
                    mout.push_back({v,d});
                ++ptr;
                s=false;
                break;
            }
            if (m.comp(v,ptr->first))
            {
                if (d!=0)
                    mout.push_back({v,d});
                s=false;
                break;
            }
            mout.push_back(*ptr);
        }
        for(;ptr!=m.end();++ptr)
            mout.push_back(*ptr);
        if (s && d!=0)
            mout.push_back({v,d});
        return mout;
    }

    inline basic_monomial<univariate_priority_order> __change_up_monomial_var_deg(const basic_monomial<univariate_priority_order> &m,uint64_t d)
    {
        basic_monomial<univariate_priority_order>  mout(m.comp_ptr());
        auto ptr=m.begin();
        if (d==0)
        {
            if (get_up_deg(m)==0)
                return m;
            
            ++ptr;
        }
        else
        {
            mout.push_back({get_up_var(m.comp()),d});
            if (get_up_deg(m)!=0)
                ++ptr;
        }
        for(;ptr!=m.end();++ptr)
            mout.push_back(*ptr);
        return mout;
    }

    template<class Tc>
    std::vector<polynomial_<Tc,univariate_priority_order>>  coeff(const polynomial_<Tc,univariate_priority_order> &F)
    {
        int64_t m=get_up_deg(F);
        polynomial_<Tc,univariate_priority_order>  tmp(F.comp_ptr());
        std::vector<polynomial_<Tc,univariate_priority_order>>  coF(get_up_deg(F)+1,tmp);
        for (auto &i:F)
        {
            if (get_up_deg(i.first)!=m)
            {
                coF[m]=std::move(tmp);
                m=get_up_deg(i.first);
            }
            tmp.push_back({__change_up_monomial_var_deg(i.first,0),i.second});
        }
        coF[m]=std::move(tmp);
        return coF;
    }
    
    template<class Tc,class comp>
    std::vector<polynomial_<Tc,comp>>  coeff(const polynomial_<Tc,comp> &F, variable v)
    {
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order> F_(&comp_v);
        F_=F;
        auto coeff_l=coeff(F_);
        std::vector<polynomial_<Tc,comp>>  coF(coeff_l.size(),polynomial_<Tc,comp>(F.comp_ptr()));
        for (size_t i=0;i<coeff_l.size();++i)
        {
            coF[i]=coeff_l[i];
        }
        return coF;
    }

    template <class Tc>
    inline int64_t leadcoeff(polynomial_<Tc,univariate_priority_order>&O,const polynomial_<Tc,univariate_priority_order>&F)
    {
        assert(&O!=&F);
        int64_t d=get_up_deg(F);
        if (d && !F.empty())
        {
            O.clear();
            O.comp(F.comp_ptr());
            for (auto & i:F)
            {
                if (get_up_deg(i.first)!=d)
                    return d;
                O.push_back({__change_up_monomial_var_deg(i.first,0),i.second});
                // O.push_back(i);
                // O.back().first.pop_back();
            }
        }
        else
            O=F;
        return d;
    }

    template<class var_order>
    polynomial_<ZZ,lex_<var_order>> leadcoeff(const polynomial_<ZZ,lex_<var_order>> &F_)
    {
        polynomial_<ZZ,lex_<var_order>>  lc(F_.comp_ptr());
        auto v=get_first_var(F_);
        int64_t deg=get_first_deg(F_);
        basic_monomial<lex_<var_order>> m(F_.comp_ptr());
        for (auto &i:F_)
        {
            if ((!i.first.empty() &&  i.first.front().first==v && i.first.front().second == deg))
            {
                m.clear();
                m.reserve(i.first.size());
                auto ptr=i.first.begin();
                ++ptr;
                for (;ptr!=i.first.end();++ptr)
                {
                    m.push_back(*ptr);
                }
                lc.push_back({std::move(m),i.second});
            }
            else{
                break;
             }
        }
        return lc;
    }
    template <class Tc,class comp >
    inline int64_t leadcoeff(polynomial_<Tc,comp>&O,const polynomial_<Tc,comp>&F,const variable & v)
    {
        assert(&O!=&F);
        int64_t deg=0,tmp;
        for (auto &i:F)
            if ((tmp=i.first.deg(v))>deg)
                deg=tmp;
        if (deg)
        {
            O.clear();
            O.comp(F.comp_ptr());
            typename basic_monomial<comp>::const_iterator tmp_I;
            for (auto &i:F)
                if ((tmp_I=i.first.find(v))!=i.first.end() && tmp_I->second==deg)
                {
                    O.push_back({__change_monomial_var_deg(i.first,v,0),i.second});
                }
            O.normalization();
        }
        else 
            O=F;
        return deg;
    }
    template <class Tc,class comp >
    inline polynomial_<Tc,comp> leadcoeff(const polynomial_<Tc,comp>&F,const variable & v)
    {
        polynomial_<Tc,comp> O(F.comp_ptr());
        leadcoeff(O,F,v);
        return O;
    }
    


    template <class Tc,class comp>
    inline polynomial_<Tc,comp> prem(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v,bool is_L=true)
    {
        if (F.size()==0)
        {
//#ifndef NDEBUG  
            throw std::invalid_argument("Error:polynomial prem div 0.");
//#endif            
            return polynomial_<Tc,comp>();
        }
        polynomial_<Tc,comp> O(G.comp_ptr());
        if (is_number(F))
            return  O;
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order>  G1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        polynomial_<Tc,univariate_priority_order>  O1(&comp_v);
        prem(O1,G1,F1,v,is_L);
        poly_convert(std::move(O1),O);
        return O;
    }
    template <class Tc,class var_order>
    inline polynomial_<Tc,lex_<var_order>> prem(const polynomial_<Tc,lex_<var_order>> &G,const polynomial_<Tc,lex_<var_order>> & F,const variable & v,bool is_L=true)
    {
        if (F.size()==0)
        {
//#ifndef NDEBUG  
            throw std::invalid_argument("Error:polynomial prem div 0.");
//#endif            
            return polynomial_<Tc,lex_<var_order>>();
        }
        polynomial_<Tc,lex_<var_order>> O(G.comp_ptr());
        if (is_number(F))
            return O;
        if (comp_consistent(F.comp(),G.comp()) && (v==get_first_var(F) || F.comp()(v,get_first_var(F)))&& (is_number(G) || v==get_first_var(G) || F.comp()(v,get_first_var(G))))
        {
            prem(O,G,F,v,is_L);
            return O;
        }
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order>  G1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        polynomial_<Tc,univariate_priority_order>  O1(&comp_v);
        prem(O1,G1,F1,v,is_L);
        poly_convert(std::move(O1),O);
        return O;
    }
    


    template <class Tc,class var_order>
    void __onestep__prem_v1(
        const std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> & G,
        const std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> & F0,
        int64_t f_deg,
        const std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> & F,
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> & O,
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> & F1,
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> & G1,
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> & G2,
        variable v,
        const lex_<var_order> & comp 
    )
    {
        auto G_ptr=G.begin();
        auto G_end=G.end();
        int64_t g_deg=get_first_deg(G_ptr->first);
        G1.clear();
        G1.reserve(G.size());
        while(G_ptr!=G_end && get_first_deg(G_ptr->first)==g_deg && get_first_var(G_ptr->first)==v)
        {
            G1.push_back({__change_monomial_var_deg(G_ptr->first,v,g_deg-f_deg),G_ptr->second});
            ++G_ptr;
        }

        pair_vec_multiplies(F1,G1,F,comp);
        if(G_ptr==G_end)
        {
            pair_vec_negate(F1);
            std::swap(F1,O);
        }
        else{
            G2.clear();
            G2.reserve((G_end-G_ptr));
            while(G_ptr!=G_end)
                G2.push_back(*(G_ptr++));
            pair_vec_multiplies(O,G2,F0,comp);
            pair_vec_sub(G2,O,F1,comp);
            std::swap(G2,O);
        }
    }
    
    template <class Tc,class var_order>
    void prem(
            polynomial_<Tc,lex_<var_order>>&O,
            const polynomial_<Tc,lex_<var_order>>&G,
            const polynomial_<Tc,lex_<var_order>>&F,
            variable v,
            bool is_L=true
            )
    {
        const lex_<var_order> &  comp=F.comp();
        
        if (F.size()==0)
        {
//#ifndef NDEBUG  
            throw std::invalid_argument("Error:polynomial prem div 0.");
//#endif            
            return void();
        }
        O.comp(&comp);
        O.clear();
        if (is_number(F))
            return void();
        // std::cout<< F<<std::endl;
        // std::cout<< G<<std::endl;
        // std::cout<< v<<std::endl;
        
        assert((v==get_first_var(F) || comp(v,get_first_var(F)))&& (is_number(G) || v==get_first_var(G) || comp(v,get_first_var(G))));
        int64_t f_deg=get_first_deg(F);
        int64_t g_deg=get_first_deg(G);
        if (get_first_var(F)!=v)
            return void();
        if (get_first_var(G)!=v ||g_deg<f_deg)
        {
            O=G;
            return void();
        }
        int64_t d=g_deg-f_deg;
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_G;
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_G1;
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_G2;
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_F0;
        tmp_F0.reserve(F.size());
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_F_;
        tmp_F_.reserve(F.size());
        auto ptr=F.begin();
        for (;ptr!=F.end();++ptr)
            if (get_first_var(ptr->first)==v &&get_first_deg(ptr->first)==f_deg)
            {
                 tmp_F0.push_back({__change_monomial_var_deg(ptr->first,v,0),ptr->second});
            }
            else
                break;
        for (;ptr!=F.end();++ptr)
            tmp_F_.push_back(*ptr);    
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_F1;
        O.data().reserve(G.size());
        __onestep__prem_v1(G.data(),tmp_F0,f_deg,tmp_F_,O.data(),tmp_F1,tmp_G1,tmp_G2,v,comp);
        // std::cout<<"prem_v1_:"<<O<<std::endl;
        while (get_first_var(O)==v && get_first_deg(O)>=f_deg)
        {
            std::swap(tmp_G,O.data());
            __onestep__prem_v1(tmp_G,tmp_F0,f_deg,tmp_F_,O.data(),tmp_F1,tmp_G1,tmp_G2,v,comp);    
            --d;
            // std::cout<<"prem_v1_:"<<O<<std::endl;
        }
        if (is_L && d>0)
        {
            if (d==1)
            {
                pair_vec_multiplies(tmp_G,O.data(),tmp_F0,comp);
                std::swap(tmp_G,O.data());
            }
            else
            {
                pair_vec_power(tmp_G1,tmp_F0,d,comp);
                pair_vec_multiplies(tmp_G,O.data(),tmp_G1,comp);
                std::swap(tmp_G,O.data());
            }
            
        }
        
        //L.data()=std::move(tmp_F0);
        // L=pow(L,d);
    }


    template <class Tc,class var_order>
    void pquo(
            polynomial_<Tc,lex_<var_order>>&Q,
            polynomial_<Tc,lex_<var_order>>&R,
            const polynomial_<Tc,lex_<var_order>>&G,
            const polynomial_<Tc,lex_<var_order>>&F,
            variable v,
            bool is_L=true
            )
    {
        const lex_<var_order> &  comp=F.comp();
        
        if (F.size()==0)
        {
//#ifndef NDEBUG  
            throw std::invalid_argument("Error:polynomial prem div 0.");
//#endif            
            return void();
        }
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_G;
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_G1;
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_G2;
        R.comp(&comp);
        R.clear();
        Q.comp(&comp);
        Q.clear();
       
        if (is_number(F)||get_first_var(F)!=v)
        {
            int64_t d=(v==get_first_var(G))?get_first_deg(G):0;
            if (d)
            {
                if(d==1)
                {
                    pair_vec_multiplies(Q.data(),G.data(),F.data(),comp);
                }
                else
                {
                    pair_vec_power(tmp_G1,F.data(),d,comp);
                    pair_vec_multiplies(Q.data(),G.data(),tmp_G1,comp);
                }
            }
            else
                Q=G;
            return void();
        }
        // std::cout<< F<<std::endl;
        // std::cout<< G<<std::endl;
        // std::cout<< v<<std::endl;
        
        assert((v==get_first_var(F) || comp(v,get_first_var(F)))&& (is_number(G) || v==get_first_var(G) || comp(v,get_first_var(G))));
        int64_t f_deg=get_first_deg(F);
        int64_t g_deg=get_first_deg(G);
        int64_t d=g_deg-f_deg;  
        if (get_first_var(G)!=v ||g_deg<f_deg)
        {
            R=G;
            return void();
        }

        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_F0;
        tmp_F0.reserve(F.size());
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_F_;
        tmp_F_.reserve(F.size());
        auto ptr=F.begin();
        for (;ptr!=F.end();++ptr)
            if (get_first_var(ptr->first)==v &&get_first_deg(ptr->first)==f_deg)
            {
                 tmp_F0.push_back({__change_monomial_var_deg(ptr->first,v,0),ptr->second});
            }
            else
                break;
        for (;ptr!=F.end();++ptr)
            tmp_F_.push_back(*ptr);    
        std::vector<std::pair<basic_monomial<lex_<var_order>>,Tc>> tmp_F1;
        R.data().reserve(G.size());
        __onestep__prem_v1(G.data(),tmp_F0,f_deg,tmp_F_,R.data(),tmp_F1,tmp_G1,tmp_G2,v,comp);
        std::swap(tmp_G1,Q.data());
        // std::cout<<"prem_v1_:"<<R<<std::endl;
        while (get_first_var(R)==v && get_first_deg(R)>=f_deg)
        {
            std::swap(tmp_G,R.data());
            __onestep__prem_v1(tmp_G,tmp_F0,f_deg,tmp_F_,R.data(),tmp_F1,tmp_G1,tmp_G2,v,comp);  
            pair_vec_multiplies(tmp_G2,tmp_F0,Q.data(),comp);
            tmp_G2.reserve(tmp_G2.size()+tmp_G1.size());
            for (auto &i:tmp_G1)
                tmp_G2.push_back(std::move(i));
            std::swap(tmp_G2,Q.data());
            --d;
            // std::cout<<"prem_v1_:"<<R<<std::endl;
        }
        if (is_L && d>0)
        {
            if (d==1)
            {
                pair_vec_multiplies(tmp_G,R.data(),tmp_F0,comp);
                std::swap(tmp_G,R.data());
                pair_vec_multiplies(tmp_G,Q.data(),tmp_F0,comp);
                std::swap(tmp_G,Q.data());
            }
            else
            {
                pair_vec_power(tmp_G1,tmp_F0,d,comp);
                pair_vec_multiplies(tmp_G,R.data(),tmp_G1,comp);
                std::swap(tmp_G,R.data());
                pair_vec_multiplies(tmp_G,Q.data(),tmp_G1,comp);
                std::swap(tmp_G,Q.data());
            }
            
        }
        
        //L.data()=std::move(tmp_F0);
        // L=pow(L,d);
    }
    
 
    template <class Tc,class comp>
    inline polynomial_<Tc,comp> pquo(polynomial_<Tc,comp> &R,const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v,bool is_L=true)
    {
        if (F.size()==0)
        {
//#ifndef NDEBUG  
            throw std::invalid_argument("Error:polynomial prem div 0.");
//#endif            
            return polynomial_<Tc,comp>();
        }
        polynomial_<Tc,comp> Q(G.comp_ptr());
        R.clear();
        R.comp(G.comp_ptr());
        if (is_number(F))
        {
            auto d=degree(G,v);
            if (d)
                Q=G*pow(F,d);
            else
                Q=G;
            return Q;

        }
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order>  G1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        polynomial_<Tc,univariate_priority_order>  Q1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  R1(&comp_v);
        pquo(Q1,R1,G1,F1,v,is_L);
        poly_convert(std::move(Q1),Q);
        poly_convert(std::move(R1),R);
        return Q;
    }
    template <class Tc,class var_order>
    inline polynomial_<Tc,lex_<var_order>> pquo(polynomial_<Tc,lex_<var_order>> &R,const polynomial_<Tc,lex_<var_order>> &G,const polynomial_<Tc,lex_<var_order>> & F,const variable & v,bool is_L=true)
    {
        if (F.size()==0)
        {
//#ifndef NDEBUG  
            throw std::invalid_argument("Error:polynomial prem div 0.");
//#endif            
            return polynomial_<Tc,lex_<var_order>>();
        }
        polynomial_<Tc,lex_<var_order>> Q(G.comp_ptr());
        R.clear();
        R.comp(G.comp_ptr());
        if (is_number(F))
        {
            auto d=degree(G,v);
            if (d)
                Q=G*pow(F,d);
            else
                Q=G;
            return Q;

        }
        if (comp_consistent(F.comp(),G.comp()) && (v==get_first_var(F) || F.comp()(v,get_first_var(F)))&& (is_number(G) || v==get_first_var(G) || F.comp()(v,get_first_var(G))))
        {
            pquo(Q,R,G,F,v,is_L);
            return Q;
        }
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order>  G1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        polynomial_<Tc,univariate_priority_order>  Q1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  R1(&comp_v);
        pquo(Q1,R1,G1,F1,v,is_L);
        poly_convert(std::move(Q1),Q);
        poly_convert(std::move(R1),R);
        return Q;
    }
    template <class Tc,class comp>
    inline polynomial_<Tc,comp> pquo(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v,bool is_L=true)
    {
        polynomial_<Tc,comp> R(G.comp_ptr());
        return  pquo(R,G,F,v,is_L);
    }

    template <class compare>
    polynomial_<Zp,compare> polynomial_mod(const polynomial_<ZZ,compare> & p, uint32_t prime)
    {
        polynomial_<Zp,compare> new_p(p.comp_ptr());
        // std::cout<<p<<" mod "<<prime;
        Zp coeff(prime);
        for (auto & i:p)
        {
            coeff=i.second; 
            if (coeff)
                new_p.push_back({i.first,std::move(coeff)});
        }
        // std::cout<<" ="<<new_p<<std::endl;
        return new_p;
    }
    template <class compare>
    polynomial_<Zp,compare> polynomial_mod(polynomial_<ZZ,compare> && p, uint32_t prime)
    {
        polynomial_<Zp,compare> new_p(p.comp_ptr());
        Zp coeff(prime);
        for (auto & i:p)
        {
            coeff=i.second; 
            if (coeff)
            new_p.push_back({std::move(i.first),std::move(coeff)});
        }
        p.clear();
        return new_p;
    }

    template <class Tc, class To>
    polynomial_<Tc,To> assign(const polynomial_<Tc,To>& p,const variable & v,const Tc & c)
    {
        polynomial_<Tc,To> Pout(p.comp_ptr());
        basic_monomial<To> m(p.comp_ptr());
        basic_monomial<To> m1(p.comp_ptr());
        bool f=true;
        Tc z, z1;
        for (auto &i:p)
        {

            m.clear();m.reserve(i.first.size());
            z=i.second;
            for (auto& j:i.first)
                if (j.first==v)
                {
                    z*=pow(c,j.second);
                }
                else
                    m.push_back(j);
            if (!f && m!=m1)
            {
                if (z1)
                    Pout.push_back({std::move(m1),std::move(z1)});
                m1=std::move(m);
                z1=std::move(z);
            }
            else if (f)
            {
                m1=std::move(m);
                z1=std::move(z);
                f=false;
            }
            else
            {
                z1+=z;
            }

            
        }
        if (z1)
            Pout.push_back({std::move(m1),std::move(z1)});
        Pout.normalization();
        return Pout;
    }
    
    template <class Tc, class To>
    polynomial_<Tc,To> assign(const polynomial_<Tc,To>& p,const std::map<variable,Tc> & ass_list)
    {
        polynomial_<Tc,To> Pout(p.comp_ptr());
        basic_monomial<To> m(p.comp_ptr());
        basic_monomial<To> m1(p.comp_ptr());
        bool f=true;
        Tc z, z1;
        for (auto &i:p)
        {

            m.clear();m.reserve(i.first.size());
            z=i.second;
            for (auto& j:i.first)
            {
                auto ptr=ass_list.find(j.first);
                if (ptr!=ass_list.end())
                {
                    z*=pow(ptr->second,j.second);
                }
                else
                    m.push_back(j);
            }
            if (!f && m!=m1)
            {
                if (z1)
                    Pout.push_back({std::move(m1),std::move(z1)});
                m1=std::move(m);
                z1=std::move(z);
            }
            else if (f)
            {
                m1=std::move(m);
                z1=std::move(z);
                f=false;
            }
            else
            {
                z1+=z;
            }

            
        }
        if (z1)
            Pout.push_back({std::move(m1),std::move(z1)});
        Pout.normalization();
        return Pout;
    }
    template <class Tc1,class Tc2,class Tc3, class To>
    polynomial_<Tc1,To> assign(const polynomial_<Tc2,To>& p,const std::map<variable,Tc3> & ass_list)
    {
        polynomial_<Tc1,To> Pout(p.comp_ptr());
        basic_monomial<To> m(p.comp_ptr());
        basic_monomial<To> m1(p.comp_ptr());
        bool f=true;
        Tc1 z, z1;
        for (auto &i:p)
        {

            m.clear();m.reserve(i.first.size());
            z=i.second;
            for (auto& j:i.first)
            {
                auto ptr=ass_list.find(j.first);
                if (ptr!=ass_list.end())
                {
                    z*=pow(ptr->second,j.second);
                }
                else
                    m.push_back(j);
            }
            if (!f && m!=m1)
            {
                if (z1)
                    Pout.push_back({std::move(m1),std::move(z1)});
                m1=std::move(m);
                z1=std::move(z);
            }
            else if (f)
            {
                m1=std::move(m);
                z1=std::move(z);
                f=false;
            }
            else
            {
                z1+=z;
            }

            
        }
        if (z1)
            Pout.push_back({std::move(m1),std::move(z1)});
        Pout.normalization();
        return Pout;
    }
    template <class Tc,class Tm>
    polynomial_<Tc,Tm>  derivative(const polynomial_<Tc,Tm> & p,variable x)
    {
        polynomial_<Tc,Tm> Pout(p.comp_ptr());
        basic_monomial<Tm> m(p.comp_ptr());
        int64_t b;
        for (auto &i:p)
        {
            m.clear();
            b=0;
            for (auto &j:i.first)
            {
                if (j.first==x)
                {
                    b=j.second;
                    if (j.second-1)
                        m.push_back({x,j.second-1});
                }
                else
                {
                    m.push_back(j);
                }
            }
            if (b)
            {
                Pout.push_back({std::move(m),b*i.second});
            }
            
        }
        return Pout;
    }
    template <class Tc,class var_order>
    polynomial_<Tc,lex_<var_order>>  derivative(const polynomial_<Tc,lex_<var_order>> & p)
    {
        polynomial_<Tc,lex_<var_order>> Pout(p.comp_ptr());
        basic_monomial<lex_<var_order>> m(p.comp_ptr());
        auto x=get_first_var(p);
        for (auto &i:p)
        {
            
            if (!i.first.empty() && x==i.first.front().first)
            {
                auto ptr=i.first.begin();
                m.clear();
                if (ptr->second-1)
                    m.push_back({x,ptr->second-1});
                for (++ptr;ptr!=i.first.end();++ptr) m.push_back(*ptr);
                Pout.push_back({std::move(m),i.first.front().second*i.second});
            }
            
        }
        return Pout;
    }

    template <class Tc,class Tm>
    constexpr bool is_number(const polynomial_<Tc,Tm> & F)
    {
        return (F.empty() || F.size()==1 && F.front().first.empty());
    }

}
#endif