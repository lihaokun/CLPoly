/*
Module Name:
    polynomial.hh
Abstract:
    定义polynomial和相关的运算
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_POLYNOMIAL_HH
#define CLPOLY_POLYNOMIAL_HH

#include <clpoly/polynomial_.hh>
#include <list>
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
    // polynomial_<Tc,To> association(const polynomial_<Tc,To>& p,const variable & v,const Tc & c);
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
        for (;ptr!=ptr.end();++ptr)
            if ((tmp=ptr->first.deg(v))>deg)
                deg=tmp;
        return deg;
    }
    
    template <class Tc>
    constexpr int64_t get_up_deg(const polynomial_<Tc,univariate_priority_order>& p)
    {
        return p.empty()?0:get_up_deg(p.begin()->first);
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
    constexpr std::pair<variable,int64_t> get_last_var_deg(const polynomial_<Tc,lex_<var_order>>& p)
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
    template <class Tc>
    int64_t degree(const polynomial_<Tc,univariate_priority_order> & p,const variable & v)
    {
        if (get_up_var(p.comp())==v)
            return get_up_deg(p);
        if (p.empty())
            return 0;
        auto ptr=p.begin();
        int64_t deg=(ptr++)->first.deg(v),tmp;
        for (;ptr!=ptr.end();++ptr)
            if ((tmp=ptr->first.deg(v))>deg)
                deg=tmp;
        return deg;
    }


    basic_monomial<univariate_priority_order> __change_up_monomial_var_deg(const basic_monomial<univariate_priority_order> &m,uint64_t d)
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
            basic_monomial<comp> m(F.comp_ptr());
            typename basic_monomial<comp>::const_iterator tmp_I;
            for (auto &i:F)
                if ((tmp_I=i.first.find(v))!=i.first.end() && tmp_I->second==deg)
                {
                    m=i.first;
                    m[tmp_I-i.first.begin()].second=0;
                    O.push_back({std::move(m),i.second});
                }
            O.normalization();
        }
        else 
            O=F;
        return deg;
    }
    
    template <class Tc>
    inline polynomial_<Tc,univariate_priority_order> leadcoeff(const polynomial_<Tc,univariate_priority_order>&F)
    {
        polynomial_<Tc,univariate_priority_order> O(F.comp_ptr());
        leadcoeff(O,F);
        return O;
    }
    
    template <class Tc,class comp>
    inline polynomial_<Tc,comp> leadcoeff(const polynomial_<Tc,comp>&F)
    {
        polynomial_<Tc,comp> O(F.comp_ptr());
        leadcoeff(O,F);
        return O;
    }

    // #define prem_v1 prem
    
    template <class Tc>
    inline polynomial_<Tc,univariate_priority_order> prem(const polynomial_<Tc,univariate_priority_order> &G,const polynomial_<Tc,univariate_priority_order> & F)
    {
        polynomial_<Tc,univariate_priority_order>  O1(&G.comp_ptr());
        prem(O1,G,F);
        return O1;
    }
    template <class Tc,class comp>
    inline polynomial_<Tc,comp> prem(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v)
    {
        if (F.size()==0)
        {
//#ifndef NDEBUG  
            throw std::invalid_argument("Error:polynomial prem div 0.");
//#endif            
            return polynomial_<Tc,comp>();
        }
        polynomial_<Tc,comp> O(G.comp_ptr());
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order>  G1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        polynomial_<Tc,univariate_priority_order>  O1(&comp_v);
        prem(O1,G1,F1);
        poly_convert(std::move(O1),O);
        return O;
    }

    


    template <class Tc,class compare>
    void __onestep__prem_v1(
        const std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & G,
        const std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & F0,
        int64_t f_deg,
        const std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & F,
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & O,
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & F1,
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & G1,
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & G2,
        compare comp 
    )
    {

        auto G_ptr=G.begin();
        auto G_end=G.end();
        int64_t g_deg=get_up_deg(G_ptr->first);
        G1.clear();
        G1.reserve(G.size());
        while(G_ptr!=G_end && get_up_deg(G_ptr->first)==g_deg)
        {
            G1.push_back({__change_up_monomial_var_deg(G_ptr->first,g_deg-f_deg),G_ptr->second});
            ++G_ptr;
        }
        // if (g_deg-f_deg)
        // {
        //     while(G_ptr!=G_end && get_up_deg(G_ptr->first)==g_deg)
        //     {
        //         G1.push_back(*(G_ptr++));
        //         G1.back().first.back().second=g_deg-f_deg;
        //         G1.back().first.deg()-=f_deg;
        //     }
        // }
        // else
        // {
        //     while(G_ptr!=G_end && get_up_deg(G_ptr->first)==g_deg)
        //     {
        //         G1.push_back(*(G_ptr++));
        //         G1.back().first.pop_back();
        //     }
        // }

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
            pair_vec_multiplies(G1,G2,F0,comp);
            pair_vec_sub(O,G1,F1,comp);
        }
    }
    
    template <class Tc>
    void prem(
            polynomial_<Tc,univariate_priority_order>&O,
            //polynomial_<Tc,univariate_priority_order>&L,
            const polynomial_<Tc,univariate_priority_order>&G,
            const polynomial_<Tc,univariate_priority_order>&F,
            bool is_L=true
            )
    {
        if (F.size()==0)
        {
//#ifndef NDEBUG  
            throw std::invalid_argument("Error:polynomial prem div 0.");
//#endif            
            return void();
        }
        O.clear();
        int64_t f_deg=get_up_deg(F);
        int64_t g_deg=get_up_deg(G);
        if (f_deg<=0)
            return void();
        if (g_deg<f_deg)
        {
            O=G;
            return void();
        }
        int64_t d=g_deg-f_deg;
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_G;
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_G1;
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_G2;
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_F0;
        tmp_F0.reserve(F.size());
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_F_;
        tmp_F_.reserve(F.size());
        for (auto & i:F)
            if (get_up_deg(i.first)!=f_deg)
                tmp_F_.push_back(i);
            else
            {
                 tmp_F0.push_back({__change_up_monomial_var_deg(i.first,0),i.second});
                // tmp_F0.push_back(i);
                // tmp_F0.back().first.pop_back();
            }
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_F1;
        O.data().reserve(G.size());
        __onestep__prem_v1(G.data(),tmp_F0,f_deg,tmp_F_,O.data(),tmp_F1,tmp_G1,tmp_G2,G.comp());
        // std::cout<<"prem_v1_:"<<O<<std::endl;
        while (get_up_deg(O)>=f_deg)
        {
            std::swap(tmp_G,O.data());
            __onestep__prem_v1(tmp_G,tmp_F0,f_deg,tmp_F_,O.data(),tmp_F1,tmp_G1,tmp_G2,G.comp());    
            --d;
            // std::cout<<"prem_v1_:"<<O<<std::endl;
        }
        if (is_L && d>0)
        {
            if (d==1)
            {
                pair_vec_multiplies(tmp_G,O.data(),tmp_F0,G.comp());
                std::swap(tmp_G,O.data());
            }
            else
            {
                pair_vec_power(tmp_G1,tmp_F0,d,G.comp());
                pair_vec_multiplies(tmp_G,O.data(),tmp_G1,G.comp());
                std::swap(tmp_G,O.data());
            }
            
        }
        
        //L.data()=std::move(tmp_F0);
        // L=pow(L,d);
    }
    
    template <class Tc,class comp>
    inline polynomial_<Tc,comp> resultant(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v)
    {
        assert(comp_consistent(G.comp(),F.comp()));
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order>  G1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        polynomial_<Tc,univariate_priority_order>  O1(&comp_v);
        resultant(O1,G1,F1);
        polynomial_<Tc,comp> O(G.comp_ptr());
        poly_convert(std::move(O1),O);
        return O;
    }
    
    
    template <class Tc>
    void resultant
        (   
            polynomial_<Tc,univariate_priority_order>&O,
            const polynomial_<Tc,univariate_priority_order>&F,
            const polynomial_<Tc,univariate_priority_order>&G
            
        )
    {
        O.clear();
        const univariate_priority_order &comp=G.comp();
        int64_t l=get_up_deg(G);
        int64_t m=get_up_deg(F);
        if (m<0 || l<0 || G.empty()||F.empty())
            return void();
        if (m<l)
        {
            resultant(O,G,F);
            if ((l&1) &&(m&1))
                pair_vec_negate(O.data());
            return void();
        }
        if (m==0)
        {
            O=pow(F,l);
            return void();
        }
        if (l==0)
        {
            O=pow(G,m);
            return  void();
        }

        polynomial_<Tc,univariate_priority_order> S_j_1(&comp);
        polynomial_<Tc,univariate_priority_order> S_j(&comp);
        polynomial_<Tc,univariate_priority_order> S_r_1(&comp);
        polynomial_<Tc,univariate_priority_order> S_r(&comp);
        
        polynomial_<Tc,univariate_priority_order> R_(&comp);
        polynomial_<Tc,univariate_priority_order> tmp1(&comp);
        polynomial_<Tc,univariate_priority_order> tmp2(&comp);
        int64_t j,r,u;
        // std::cout<<"F="<<F<<std::endl;
        // std::cout<<"G="<<G<<std::endl;
        
        if (l<m)
        {
            j=m-1;
            if(j==l)
            {
                prem(S_j,F,G);
                --j;
                //  std::cout<<j<<":"<<S_j<<std::endl;
                S_j_1=G;
            }
            else{
                leadcoeff(R_,G);
                if (j-l<2)
                    pair_vec_multiplies(S_r.data(),R_.data(),G.data(),comp);
                else 
                {
                    pair_vec_power(tmp1.data(),R_.data(),j-l,comp);
                    pair_vec_multiplies(S_r.data(),tmp1.data(),G.data(),comp);
                }
                //  std::cout<<l<<":"<<S_r<<std::endl;
                prem(S_r_1,F,G);
                if (j-l & 1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<l-1<<":"<<S_r_1<<std::endl;
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
                j=l-1;
            }

        }else
        {
            j=m;
            prem(S_j,F,G);
            pair_vec_negate(S_j.data());
            //  std::cout<<j-1<<":"<<S_j<<std::endl;
            if (!(--j))
            {
                O.data()=std::move(S_j.data());
                return void();
            }
            r=get_up_deg(S_j);
            if (r<0)
                return void();
            if (r<j)
            {
                
                if (r!=1)
                {
                    leadcoeff(R_,S_j);
                    if (j-r<2)
                        pair_vec_multiplies(S_r.data(),R_.data(),S_j.data(),comp);
                    else 
                    {
                        pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                        pair_vec_multiplies(S_r.data(),tmp1.data(),S_j.data(),comp);
                    }
                }
                
                //  std::cout<<r<<":"<<S_r<<std::endl;
                if (!r)
                {
                    O.data()=std::move(S_r.data());
                    return void();
                }
                prem(tmp1,G,S_j);
                //std::cout<<"G"<<":"<<G<<std::endl;
                //std::cout<<"S_"<<j<<":"<<S_j<<std::endl;
                
                //  std::cout<<r-1<<":"<<tmp1<<std::endl;
                leadcoeff(R_,G);
                pair_vec_div(S_r_1.data(),tmp1.data(),R_.data(),comp);
                if ((j-r) & 1==1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<r-1<<":"<<S_r_1<<std::endl;
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
                j=r-1;
            }
            else
            {
                prem(tmp1,G,S_j);
                leadcoeff(R_,G);
                pair_vec_div(S_j_1.data(),tmp1.data(),R_.data(),comp);
                //  std::cout<<j-1<<":"<<S_j_1<<std::endl;
                --j;
                swap(S_j.data(),S_j_1.data());
            }
        }
        while (j>0)
        {
            r=get_up_deg(S_j);
            if (r<0)
                return void();
            if (r<j)
            {
                
                //std::cout<<"R_j"<<":"<<tmp1<<std::endl;
                leadcoeff(R_,S_j_1);
                //std::cout<<"R_j+1"<<":"<<R_<<std::endl;
                if (r!=1)
                {
                    leadcoeff(tmp1,S_j);
                    if (j-r<2)
                    {
                        pair_vec_multiplies(tmp2.data(),tmp1.data(),S_j.data(),comp);
                        pair_vec_div(S_r.data(),tmp2.data(),R_.data(),comp);
                    }
                    else 
                    {
                        pair_vec_power(S_r.data(),tmp1.data(),j-r,comp);
                        pair_vec_multiplies(tmp2.data(),S_r.data(),S_j.data(),comp);
                        pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                        pair_vec_div(S_r.data(),tmp2.data(),tmp1.data(),comp);
                    }
                }
                else
                {
                    pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                }
                
                //  std::cout<<r<<":"<<S_r<<std::endl;
                if (!r)
                {
                    O.data()=std::move(S_r.data());
                    return void();
                }
                pair_vec_multiplies(tmp2.data(),R_.data(),R_.data(),comp);
                if (j-r<2)
                {
                    pair_vec_multiplies(tmp1.data(),R_.data(),tmp2.data(),comp);    
                    prem(tmp2,S_j_1,S_j);
                    pair_vec_div(S_r_1.data(),tmp2.data(),tmp1.data(),comp);    
                }        
                else
                {
                    pair_vec_multiplies(R_.data(),tmp2.data(),tmp1.data(),comp);
                    prem(tmp2,S_j_1,S_j);
                    pair_vec_div(S_r_1.data(),tmp2.data(),R_.data(),comp);
                }
                if ((j-r) & 1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<r-1<<":"<<S_r_1<<std::endl;
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
                j=r-1;
            }
            else
            {
                prem(tmp1,S_j_1,S_j);
                // std::cout<<"tmp1:"<<tmp1<<std::endl;
                leadcoeff(R_,S_j_1);
                // std::cout<<"R_:"<<R_<<std::endl;
                pair_vec_multiplies(tmp2.data(),R_.data(),R_.data(),comp);
                // std::cout<<"tmp2:"<<tmp2<<std::endl;
                pair_vec_div(S_j_1.data(),tmp1.data(),tmp2.data(),comp);
                --j;
                //  std::cout<<j<<":"<<S_j_1<<std::endl;
                swap(S_j.data(),S_j_1.data());
            }
            
        }
        O.data()=std::move(S_j.data());
        
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
    polynomial_<Tc,To> association(const polynomial_<Tc,To>& p,const variable & v,const Tc & c)
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