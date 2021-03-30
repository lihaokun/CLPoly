/**
 * @file polynomial_gcd.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义polynomial_gcd
*/
#ifndef CLPOLY_POLYNOMIAL_GCD_HH
#define CLPOLY_POLYNOMIAL_GCD_HH
#include <clpoly/polynomial.hh>
#include <clpoly/upolynomial.hh>
#include <boost/math/special_functions/prime.hpp>
#include <cmath>
#include <cstdint>
#include <cassert>
#include <random>

namespace clpoly{ 

    template <class comp>
    basic_monomial<comp> gcd(const basic_monomial<comp> &m1,const basic_monomial<comp> &m2)
    {
        assert(comp_consistent(m1.comp(),m2.comp()));
        basic_monomial<comp> m(m1.comp_ptr());
        if (m1.empty()|| m2.empty())
        {
            return m;
        }
        m.reserve(m1.size());
        auto m1_ptr=m1.begin();
        auto m2_ptr=m2.begin();
        while (m1_ptr!=m1.end() && m2_ptr!=m2.end())
        {
            if (m1.comp(m1_ptr->first,m2_ptr->first))
            {
                ++m1_ptr;
            }
            else
            {
                if (m1_ptr->first==m2_ptr->first)
                {
                    m.push_back({m1_ptr->first,std::min(m1_ptr->second,m2_ptr->second)});
                    ++m1_ptr;
                }
                ++m2_ptr;
            }
        }
        return m;
    }
    template<class var_order>
    polynomial_<ZZ,lex_<var_order>> cont(const polynomial_<ZZ,lex_<var_order>> &F_);
    template<class var_order>
    polynomial_<ZZ,lex_<var_order>> leadcoeff(const polynomial_<ZZ,lex_<var_order>> &F_);
    template<class var_order>
    int64_t  __polynomial_GCD(       polynomial_<Zp,lex_<var_order>> & Pout,
                            const polynomial_<Zp,lex_<var_order>> & F,
                            const polynomial_<Zp,lex_<var_order>> & G,
                            const polynomial_<Zp,lex_<var_order>> & Lc_gcd,
                            int64_t deg);
    template<class var_order>
    polynomial_<ZZ,lex_<var_order>> polynomial_GCD(polynomial_<ZZ,lex_<var_order>>  F, polynomial_<ZZ,lex_<var_order>>  G);
    template <class compare>
    inline polynomial_<ZZ,compare> polynomial_GCD(const polynomial_<ZZ,compare>  &F, const polynomial_<ZZ,compare> & G)
    {
        polynomial_<ZZ,lex> F_,G_;
        poly_convert(F,F_);
        poly_convert(G,G_);
        polynomial_<ZZ,compare>  Pout(F.comp_ptr());
        poly_convert(polynomial_GCD(std::move(F_),std::move(G_)),Pout);
        return Pout;
    }
    template <class comp>
    inline polynomial_<ZZ,comp> gcd(const polynomial_<ZZ,comp> &F,const polynomial_<ZZ,comp> & G)
    {
        return polynomial_GCD(F,G);
    }
    template<class var_order>
    bool is_squarefree (const polynomial_<ZZ,lex_<var_order>> &  F)
    {
        if (F.empty())
            return false;
        if (is_number(F))
            return true;
        auto f_cont=cont(F);
        if (is_squarefree(f_cont))
        {
            auto F_=F/f_cont;
            auto F_1=polynomial_GCD(F_,derivative(F_));
            if (is_number(F_1))
                return true;
        }
        return false;
    }
    template<class comp>
    bool is_squarefree (const polynomial_<ZZ,comp> &  F)
    {
        polynomial_<ZZ,lex> F_;
        poly_convert(F,F_);
        return is_squarefree(F_);
    }
    template<class var_order>
    std::vector<std::pair<polynomial_<ZZ,lex_<var_order>>,uint64_t>> squarefree (const polynomial_<ZZ,lex_<var_order>> &  F)
    {
        // std::cout<<"F:"<<F<<std::endl;
        if (F.empty())
            return {};
        if (is_number(F))
            return {{F,1}};
        auto f_cont=cont(F);
        auto F_=F/f_cont;
        auto lst=squarefree(f_cont);
        polynomial_<ZZ,lex_<var_order>> F_1(F.comp_ptr()),F_2(F.comp_ptr()),F_3(F.comp_ptr());
        for(uint64_t i=1;;++i)
        {
            // std::cout<<"F_:"<<F_<<std::endl;
            // std::cout<<"D(F_):"<<derivative(F_)<<std::endl;
            
            F_1=polynomial_GCD(F_,derivative(F_));
            // std::cout<<"F_1:"<<F_1<<std::endl;
            if (is_number(F_1))
            {
                if (i>1 && F_!=F_3 )
                {
                    lst.push_back({F_3/F_,i-1});
                }
                lst.push_back({F_,i});
                break;
            }
            F_2=F_/F_1;
            // std::cout<<"F_2:"<<F_2<<std::endl;
            if (F_2!=F_3)
            {
                if (i>1)
                {
                    lst.push_back({F_3/F_2,i-1});
                }
                F_3=std::move(F_2);
                
            }
            F_=std::move(F_1);
        }
        return lst;
    }
    template <class comp>
    std::vector<std::pair<polynomial_<ZZ,comp>,uint64_t>> squarefree (const polynomial_<ZZ,comp> &  F)
    {
        polynomial_<ZZ,lex> F_;
        poly_convert(F,F_);
        auto lst=squarefree(F_);
        std::vector<std::pair<polynomial_<ZZ,comp>,uint64_t>> lst_(lst.size(),{polynomial_<ZZ,comp>(F.comp_ptr()),0});
        for (uint64_t i=0;i<lst.size();++i)
        {
            poly_convert(lst[i].first,lst_[i].first);
            lst_[i].second=lst[i].second;
        }
        return lst_;
    }
    template<class var_order>
    std::vector<polynomial_<ZZ,lex_<var_order>>> squarefreebasis (std::vector<polynomial_<ZZ,lex_<var_order>>>  F)
    {
        std::vector<polynomial_<ZZ,lex_<var_order>>> lst,lst_;
        for (auto &i:F)
        {
            lst_.clear();   
            auto l_=squarefree(i);
            if (l_.empty())
                continue;
            auto l_ptr=l_.begin();
            for (++l_ptr;l_ptr!=l_.end();++l_ptr)
            {
                i=std::move(l_ptr->first);
                for (auto &j:lst)
                {
                    auto F1=polynomial_GCD(j,i);
                    if (!is_number(F1))
                    {
                        i=i/F1;
                        if (F1!=j)
                        {
                            j=j/F1;
                            lst_.push_back(std::move(F1));
                        }
                        if (is_number(i))
                            break;
                    }
                }
                if (!is_number(i))
                    lst_.push_back(std::move(i));
            }
            lst.reserve(lst.size()+lst_.size());
            for (auto &j:lst_)
            {
                 lst.push_back(std::move(j));
            }
        }
        return lst;
    }
    
    template <class comp>
    std::vector<polynomial_<ZZ,comp>> squarefreebasis (const std::vector<polynomial_<ZZ,comp>> &  F)
    {
        std::vector<polynomial_<ZZ,comp>> L;
        if (F.empty())
            return L;
        std::vector<polynomial_<ZZ,lex>> lst;
        polynomial_<ZZ,lex> p;
        lst.resize(F.size());
        for (auto &i:F)
        {
            poly_convert(i,p);
            lst.push_back(std::move(p));
        }
        lst=squarefreebasis(lst);
        L.reserve(lst.size());
        polynomial_<ZZ,comp> p1(F.front().comp_ptr());
        for (auto &i:lst)
        {
            poly_convert(i,p1);
            L.push_back(std::move(p1));
        }
        return L;
    }
    





    template<class var_order>
    polynomial_<ZZ,lex_<var_order>> polynomial_GCD(polynomial_<ZZ,lex_<var_order>>  F, polynomial_<ZZ,lex_<var_order>>  G)
    {
        assert(comp_consistent(F.comp(),G.comp()));
        if (F.empty())
            return G;
        if (G.empty())
            return F;
        if (F.size()==1 || G.size()==1)
        {
            auto ptr=F.begin();
            auto m=ptr->first;
            auto c=ptr->second;
            for(++ptr;ptr!=F.end();++ptr)
            {
                m=gcd(m,ptr->first);
                c=gcd(c,ptr->second);
            }
            for(ptr=G.begin();ptr!=G.end();++ptr)
            {
                m=gcd(m,ptr->first);
                c=gcd(c,ptr->second);
            }
            polynomial_<ZZ,lex_<var_order>> Pout(F.comp_ptr());
            Pout={{m,c}};
            return Pout;
            
        }
        
        variable ffv,gfv;
        while ((ffv=get_first_var(F))!=(gfv=get_first_var(G)) || !ffv.serial())
        {
            if (!ffv.serial() || !gfv.serial())
            {
                ZZ cont=0;
                for(auto &i:F)
                    cont=gcd(cont,i.second);
                for(auto &i:G)
                    cont=gcd(cont,i.second);
                polynomial_<ZZ,lex_<var_order>> Pout(F.comp_ptr());
                Pout={{basic_monomial<lex_<var_order>>(F.comp_ptr()),cont}};
                return Pout;    
            }
            if (F.comp(ffv,gfv))
            {
                F=cont(F);
            }
            else
            {
                G=cont(G);
            }
        }
        polynomial_<ZZ,lex_<var_order>> F_cont=cont(F);
        polynomial_<ZZ,lex_<var_order>> G_cont=cont(G);
        polynomial_<ZZ,lex_<var_order>> cont_gcd=polynomial_GCD(cont(F),cont(G));

        F=F/F_cont;G=G/G_cont;
        polynomial_<ZZ,lex_<var_order>> lc_gcd=polynomial_GCD(leadcoeff(F),leadcoeff(G));
       
        
        // std::cout<<"F:"<<F<<std::endl;
        // std::cout<<"G:"<<G<<std::endl;
        // std::cout<<"cont_gcd:"<<cont_gcd<<std::endl;
        // std::cout<<"lc_gcd:"<<lc_gcd<<std::endl;
        
        std::uint32_t tmp_x=std::max(degree(F),degree(G));
        if (tmp_x<2) tmp_x=2;
        std::uint32_t p_index=tmp_x/std::log(tmp_x);
        std::uint32_t prime=boost::math::prime(p_index);
        while (prime <tmp_x)
        {
            prime=boost::math::prime(++p_index);
        }

        polynomial_<ZZ,lex_<var_order>> Pout_(F.comp_ptr()),tmp_Pout_(F.comp_ptr()),R(F.comp_ptr());
        polynomial_<Zp,lex_<var_order>> Pout_mod(F.comp_ptr()),f_p(F.comp_ptr()),g_p(F.comp_ptr()),lc_gcd_p(F.comp_ptr());
        polynomial_<ZZ,lex_<var_order>> tmp(F.comp_ptr());
        ZZ Pout_prime;
        std::int64_t Pout_d=INT64_MAX;
        std::int64_t tmp_Pout_d=INT64_MAX;
        
        while (1)
        {
            
            while (F.begin()->second % prime ==0 || G.begin()->second % prime ==0)
            {
                prime=boost::math::prime(++p_index);
            }
            f_p=polynomial_mod(F,prime);
            g_p=polynomial_mod(G,prime);
            lc_gcd_p=polynomial_mod(lc_gcd,prime);

            // std::cout<<"p:"<<prime<<std::endl;
            // std::cout<<"f_p:"<<f_p<<std::endl;
            // std::cout<<"g_p:"<<g_p<<std::endl;
            // std::cout<<"lc_gcd_p:"<<lc_gcd_p<<std::endl;

            tmp_Pout_d=__polynomial_GCD(Pout_mod,f_p,g_p,lc_gcd_p,Pout_d);
            if (tmp_Pout_d==-1)
            {
                prime=boost::math::prime(++p_index);
                continue;
            }
            
            // std::cout<<"poly_mod:"<<Pout_mod<<std::endl;
            
            if (tmp_Pout_d < Pout_d)
            {
                Pout_d=tmp_Pout_d;
                Pout_prime=prime;
                poly_convert(Pout_mod,Pout_);
                for (auto &i:Pout_)
                {
                    i.second%=Pout_prime;
                    if (i.second>Pout_prime/2)
                    {
                        i.second-=Pout_prime;
                    }
                }
                
                // std::cout<<"Pout_:"<<Pout_<<std::endl;
                
            }
            else
            {
                Zp tmp_inv(Pout_prime,prime);
                tmp_inv=tmp_inv.inv();
                tmp_Pout_.clear();
                tmp_Pout_.reserve(Pout_.size());
                auto Pout_ptr=Pout_.begin();
                auto Pout_end=Pout_.end();
                auto Pm_ptr=Pout_mod.begin();
                auto Pm_end=Pout_mod.end();
                while (Pout_ptr!=Pout_end && Pm_ptr!=Pm_end)
                {
                    if (F.comp(Pout_ptr->first,Pm_ptr->first))
                    {
                        tmp_Pout_.push_back({Pout_ptr->first,Pout_ptr->second-Pout_ptr->second*tmp_inv.number()*Pout_prime});
                        ++Pout_ptr;
                    }
                    else
                    {
                        if (Pout_ptr->first==Pm_ptr->first)
                        {
                            tmp_Pout_.push_back({std::move(Pm_ptr->first),Pout_ptr->second+
                            (Pm_ptr->second.number()-Pout_ptr->second)*tmp_inv.number()*Pout_prime});
                            ++Pout_ptr;
                        }
                        else
                        {
                            tmp_Pout_.push_back({std::move(Pm_ptr->first),Pm_ptr->second.number()*tmp_inv.number()*Pout_prime});
                            
                        }
                        ++Pm_ptr;
                    }
                    
                }
                while (Pout_ptr!=Pout_end)
                {
                    tmp_Pout_.push_back({std::move(Pout_ptr->first),Pout_ptr->second-Pout_ptr->second*tmp_inv.number()*Pout_prime});
                    ++Pout_ptr;
                }
                while (Pm_ptr!=Pm_end)
                {
                    tmp_Pout_.push_back({std::move(Pm_ptr->first),Pm_ptr->second.number()*tmp_inv.number()*Pout_prime});
                    ++Pm_ptr;
                }
                Pout_prime*=prime; 
                for (auto &i:tmp_Pout_)
                {
                    i.second%=Pout_prime;
                    if (i.second>Pout_prime/2)
                    {
                        i.second-=Pout_prime;
                    }
                    else if (i.second<-Pout_prime/2)
                    {
                        i.second+=Pout_prime;
                    }
                }
                
                // std::cout<<"Pout_:"<<tmp_Pout_<<std::endl;
                
                if (tmp_Pout_==Pout_)
                {
                    if (get_first_var(Pout_)==get_first_var(F))
                    {
                        auto cont_=cont(Pout_);
                        // std::cout<<"cont_:"<<cont_<<std::endl;
                        pair_vec_div(tmp.data(),tmp_Pout_.data(),cont_.data(),F.comp());
                        std::swap(tmp,tmp_Pout_);
                        
                        // std::cout<<tmp_Pout_<<std::endl;
                        
                        pair_vec_div(tmp.data(),R.data(),F.data(),tmp_Pout_.data(),F.comp());
                        if (R.empty())
                        {
                            pair_vec_div(tmp.data(),R.data(),G.data(),tmp_Pout_.data(),F.comp());
                            if (R.empty())
                            {
                                pair_vec_multiplies(tmp.data(),tmp_Pout_.data(),cont_gcd.data(),F.comp());
                                return tmp;
                            }
                        }
                    }
                    else
                    {
                        return cont_gcd;
                    }
                }
                swap(tmp_Pout_.data(),Pout_.data());
                       
            }
            prime=boost::math::prime(++p_index);
        }  

        
    }

    

   
    template<class var_order>
    polynomial_<ZZ,lex_<var_order>> cont(const polynomial_<ZZ,lex_<var_order>> &F_)
    {
        // std::cout<<F_<<std::endl;
        polynomial_<ZZ,lex_<var_order>>  cont(F_.comp_ptr()), tmp(F_.comp_ptr());
        auto v=get_first_var(F_);
        int64_t deg=get_first_deg(F_);
        int64_t tmp_deg=deg;
        basic_monomial<lex_<var_order>> m(F_.comp_ptr());
        for (auto &i:F_)
        {
            if ((!i.first.empty() &&  i.first.front().first==v && i.first.front().second == tmp_deg) ||  
                ((i.first.empty() || i.first.front().first!=v)  && tmp_deg==0))
            {
                m.clear();
                m.reserve(i.first.size());
                auto ptr=i.first.begin();
                if (tmp_deg)
                    ++ptr;
                for (;ptr!=i.first.end();++ptr)
                {
                    m.push_back(*ptr);
                }
                tmp.push_back({std::move(m),i.second});
            }
            else{
                if (tmp_deg==deg)
                {
                    cont=std::move(tmp);
                    // std::cout<<"c:"<<cont<<std::endl;
                }
                else
                {
                    cont=polynomial_GCD(cont,tmp);
                    // std::cout<<"c:"<<tmp<<std::endl;
                    // std::cout<<"c:"<<cont<<std::endl;
                }
                tmp_deg=(!i.first.empty() && i.first.front().first==v )?i.first.front().second : 0;
                tmp.clear();
                m.clear();
                m.reserve(i.first.size());
                auto ptr=i.first.begin();
                if (tmp_deg)
                    ++ptr;
                for (;ptr!=i.first.end();++ptr)
                {
                    m.push_back(*ptr);
                }
                tmp.push_back({std::move(m),i.second});
            }
        }
        if (tmp_deg==deg)
        {
            cont=std::move(tmp);
            // std::cout<<"c:"<<cont<<std::endl;
        }
        else
        {
            cont=polynomial_GCD(cont,tmp);
            // std::cout<<"c:"<<tmp<<std::endl;
            // std::cout<<"c:"<<cont<<std::endl;
        }
        return cont;
    } 
    // polynomial_<ZZ,univariate_priority_order> cont(const polynomial_<ZZ,univariate_priority_order> &F_)
    // {
    //     auto & v_order=F_.comp();
    //     polynomial_<ZZ,univariate_priority_order>  cont(&v_order),
    //                                                cont_(&v_order),
    //                                                tmp(&v_order);
    //     int64_t deg=get_up_deg(F_);
    //     int64_t tmp_deg=deg;
    //     for (auto &i:F_)
    //     {
    //         if (get_up_deg(i.first)==tmp_deg)
    //         {
    //             tmp.push_back(i);
    //             if (tmp_deg)
    //                 tmp.back().first.pop_back();
    //         }
    //         else{
    //             if (tmp_deg==deg)
    //             {
    //                 cont=std::move(tmp);
                   
    //             }
    //             else
    //             {
    //                 cont=polynomial_GCD(cont,tmp);
                    
    //             }
    //             tmp_deg=get_up_deg(i.first);
    //             tmp.clear();
    //             tmp.push_back(i);
    //             if (tmp_deg)
    //                 tmp.back().first.pop_back();
    //         }
    //     }
    //     if (tmp_deg==deg)
    //     {
    //         cont=std::move(tmp);
    //     }
    //     else
    //     {
    //         cont=polynomial_GCD(cont,tmp);
    //     }
    //     return cont;
    // } 
    
    template<class var_order>
    int64_t  __polynomial_GCD(       polynomial_<Zp,lex_<var_order>> & Pout,
                            const polynomial_<Zp,lex_<var_order>> & F,
                            const polynomial_<Zp,lex_<var_order>> & G,
                            const polynomial_<Zp,lex_<var_order>> & Lc_gcd,
                            int64_t deg)
    {
        int64_t deg_;
        variable fvf,gvf;
        assert(!F.empty() && !G.empty() && !Lc_gcd.empty());
        fvf=get_first_var(F);
        gvf=get_first_var(G);
        assert(fvf==gvf);

        if (!fvf.serial())
        {
            assert(Lc_gcd.size()==1 && Lc_gcd.begin()->first.empty());
            Pout=Lc_gcd;
            return 0;
        }
        auto flvd=get_last_var_deg(F);
        auto glvd=get_last_var_deg(G);
        polynomial_<Zp,lex_<var_order>> Pout_;
        polynomial_<Zp,lex_<var_order>> Pout_1;
        polynomial_<Zp,lex_<var_order>> Pout_2;
        if (fvf==flvd.first && gvf==glvd.first)
        {
            Pout=G;
            pair_vec_div(Pout_2.data(),Pout_1.data(),F.data(),Pout.data(),Pout.comp());
            //std::cout<<Pout_1<<std::endl;
            while(!Pout_1.empty())
            {
                swap(Pout.data(),Pout_.data());
                swap(Pout.data(),Pout_1.data());
                pair_vec_div(Pout_2.data(),Pout_1.data(),Pout_.data(),Pout.data(),Pout.comp());
                //std::cout<<Pout_1<<std::endl;
            }
            assert(Lc_gcd.size()==1 && Lc_gcd.begin()->first.empty());
            Zp lc_inv=Pout.front().second.inv()*Lc_gcd.begin()->second;
            for (auto &i:Pout)
                i.second*=lc_inv;
            deg_=get_first_deg(Pout);
            if (deg_<=deg)
                return deg_;
            else
                return -1;
        }
        else
        {
            int64_t f_d=get_first_deg(F);
            int64_t g_d=get_first_deg(G);
            uint32_t prime=F.begin()->second.prime();
            Zp p_(prime);
            auto & comp=F.comp();
            polynomial_<Zp,lex_<var_order>> F_v(F.comp_ptr());
            polynomial_<Zp,lex_<var_order>> G_v(F.comp_ptr());
            polynomial_<Zp,lex_<var_order>> lc_v(F.comp_ptr());
            int64_t num_s=0;
            int64_t v_d;
            variable v;
            if (flvd.first==glvd.first)
            {
                v=flvd.first;
                v_d=std::max(flvd.second,glvd.second)+1;
            }
            else
            {
                v=comp(flvd.first,glvd.first)?glvd.first:flvd.first;
                v_d=1;
            }
            std::vector<std::pair<basic_monomial<lex_<var_order>>,std::vector<Zp>>> _Pout,tmp_Pout;
            std::vector<Zp> points;
            std::vector<Zp> tmp_p;
            points.reserve(v_d);
            std::random_device rd; 
            std::mt19937 gen(rd());  
            std::vector<int> v_bool(prime,1);  
            uint32_t p_tmp;
            for (int32_t i=0;i<prime;++i)
            {
                // p_.number()=i;
                std::uniform_int_distribution<uint64_t> dis(1, prime-i);
                p_tmp=dis(gen);
                uint64_t j_tmp=0;
                for (;p_tmp>0;p_tmp-=v_bool[j_tmp],++j_tmp);
                v_bool[j_tmp-1]=0;p_.number()=j_tmp-1;

                F_v=association(F,v,p_);
                G_v=association(G,v,p_);
                lc_v=association(Lc_gcd,v,p_);
                //std::cout<<v<<"->"<<p_<<std::endl;
                //std::cout<<"F_v:"<<F_v<<std::endl;
                //std::cout<<"G_v:"<<G_v<<std::endl;
                //std::cout<<"lc_v:"<<Lc_gcd<<std::endl;

                if (get_first_deg(F_v)==f_d && get_first_var(F_v)==fvf && get_first_deg(G_v)==g_d && get_first_var(G_v)==gvf )
                {
                    deg_=__polynomial_GCD(Pout_,F_v,G_v,lc_v,deg);
                    //std::cout<<"v:"<<v<<" deg:"<<deg_<<" GDD:"<<Pout_<<std::endl;
                    if (deg_!=-1 && deg_<=deg)
                    {
                        if (deg_<deg || !num_s)
                        {
                            deg=deg_;
                            num_s=1;
                            _Pout.clear();
                            for (auto &i:Pout_)
                            {
                                _Pout.push_back({std::move(i.first),{i.second}});
                            }
                            points={p_};
                        }
                        else
                        {
                            ++num_s;
                            tmp_Pout.clear();
                            tmp_Pout.reserve(Pout_.size());
                            points.push_back(p_);
                            auto _P_ptr=_Pout.begin();
                            auto _P_end=_Pout.end();
                            auto P_ptr=Pout_.begin();
                            auto P_end=Pout_.end();
                            
                            while (_P_ptr!=_P_end && P_ptr!=P_end)
                            {
                                if (comp(P_ptr->first,_P_ptr->first))
                                {
                                    tmp_p.clear();
                                    tmp_p.resize(num_s,Zp(0,prime));
                                    tmp_p.back()=P_ptr->second;
                                    tmp_Pout.push_back({std::move(P_ptr->first),std::move(tmp_p)});
                                    ++P_ptr;
                                }
                                else{
                                    if (P_ptr->first==_P_ptr->first)
                                    {
                                        tmp_Pout.push_back(std::move(*_P_ptr));
                                        tmp_Pout.back().second.push_back(P_ptr->second);
                                        ++P_ptr;
                                    }else
                                    {
                                        tmp_Pout.push_back(std::move(*_P_ptr));
                                        tmp_Pout.back().second.push_back(Zp(0,prime));
                                    }
                                    ++_P_ptr;
                                }
                            }
                            while (P_ptr!=P_end)
                            {
                                tmp_p.clear();
                                tmp_p.resize(num_s,Zp(0,prime));
                                tmp_p.back()=P_ptr->second;
                                tmp_Pout.push_back({std::move(P_ptr->first),std::move(tmp_p)});
                                ++P_ptr;
                            }
                            while (_P_ptr!=_P_end)
                            {
                                tmp_Pout.push_back(std::move(*_P_ptr));
                                tmp_Pout.back().second.push_back(Zp(0,prime));
                                ++_P_ptr;
                            }  
                            swap(tmp_Pout,_Pout);
                            
                        }
                        if (num_s==v_d)
                        {
                            //std::cout<<std::endl;
                            std::vector<Zp> lag;
                            lag.resize(v_d*v_d);
                            Zp tmp_inv;
                            for (int64_t i=0;i<v_d;++i)
                            {
                                lag[i*v_d]=Zp(1,prime);
                                for (int64_t j=1;j<v_d;++j)
                                    lag[i*v_d+j]=Zp(0,prime);
                                for (int64_t j=0;j<v_d;++j)
                                    if (i!=j)
                                    {
                                        tmp_inv=(points[i]-points[j]).inv();
                                        for (int64_t k=v_d-1;k>0;--k)
                                        {
                                            lag[i*v_d+k]=(lag[i*v_d+k-1]-lag[i*v_d+k]*points[j])*tmp_inv;
                                        }
                                        lag[i*v_d]=(-lag[i*v_d]*points[j])*tmp_inv;

                                    }
                                         
                            }
                            Pout.clear();
                            Pout.reserve(_Pout.size());
                            basic_monomial<lex_<var_order>> m(F.comp_ptr());
                            for (auto &i:_Pout)
                            {
                                m=i.first;
                                m.push_back({v,0});
                                
                                for (int64_t j=v_d-1;j>=0;--j)
                                {
                                    p_.number()=0;
                                   
                                    for (int64_t k=0;k<v_d;++k)
                                    {
                                        p_+=i.second[k]*lag[k*v_d+j];
                                    }
                                    if (p_)
                                    {
                                        if (j)
                                        {
                                            m.back().second=j;
                                            m.deg()+=j;
                                            Pout.push_back({m,p_});
                                            m.deg()-=j;
                                        }
                                        else
                                        {
                                            Pout.push_back({std::move(i.first),p_});
                                        }
                                        
                                    }  

                                }

                            }
                            return deg;
                        }
                        
                        
                    }
                    

                }
                if (prime<=i+v_d-num_s)
                {
                    return -1;
                }
            }
            
            return -1;

        }

    }
    ZZ cont(const upolynomial_<ZZ> &G)
    {
        if (G.empty())
            return 1;
        auto ptr=G.begin();
        ZZ c=(ptr++)->second;
        for (;ptr!=G.end();++ptr)
        {
            c=gcd(c,ptr->second);
        }
        return c;

    }

    umonomial gcd(const umonomial & G,const umonomial& F)
    {
        return umonomial(std::min(G.deg(),F.deg()));
    }
     int64_t  __polynomial_GCD(upolynomial_<Zp> &Pout,
                            const upolynomial_<Zp> &G,
                            const upolynomial_<Zp> &F,
                            const Zp & Lc_gcd,
                            int64_t deg)
    {
        int64_t deg_;
        assert(!F.empty() && !G.empty() );
        upolynomial_<Zp> Pout_;
        upolynomial_<Zp> Pout_1;
        upolynomial_<Zp> Pout_2;
        Pout=G;
        pair_vec_div(Pout_2.data(),Pout_1.data(),F.data(),Pout.data(),Pout.comp());
        //std::cout<<Pout_1<<std::endl;
        while(!Pout_1.empty())
        {
            swap(Pout.data(),Pout_.data());
            swap(Pout.data(),Pout_1.data());
            pair_vec_div(Pout_2.data(),Pout_1.data(),Pout_.data(),Pout.data(),Pout.comp());
            //std::cout<<Pout_1<<std::endl;
        }
        Zp lc_inv=Pout.front().second.inv()*Lc_gcd;
        for (auto &i:Pout)
            i.second*=lc_inv;
        deg_=Pout.front().first.deg();
        if (deg_<=deg)
            return deg_;
        else
            return -1;
    }
    upolynomial_<ZZ>  polynomial_GCD(upolynomial_<ZZ> G,upolynomial_<ZZ> F)
    {
         if (F.empty())
            return G;
        if (G.empty())
            return F;
        if (F.size()==1 || G.size()==1)
        {
            auto ptr=F.begin();
            auto m=ptr->first;
            auto c=ptr->second;
            for(++ptr;ptr!=F.end();++ptr)
            {
                m=gcd(m,ptr->first);
                c=gcd(c,ptr->second);
            }
            for(ptr=G.begin();ptr!=G.end();++ptr)
            {
                m=gcd(m,ptr->first);
                c=gcd(c,ptr->second);
            }
            upolynomial_<ZZ> Pout;
            Pout={{m,c}};
            return Pout;
            
        }

        ZZ f_cont=cont(F);
        ZZ g_cont=cont(G);
        ZZ cont_gcd=gcd(f_cont,g_cont);
        for (auto &i:F)
            i.second/=f_cont;
        for (auto &i:G)
            i.second/=g_cont;
            
        std::uint32_t p_index=0;
        std::uint32_t prime=boost::math::prime(p_index);
        upolynomial_<ZZ> Pout_,tmp_Pout_,R;
        upolynomial_<Zp> Pout_mod,f_p,g_p;
        upolynomial_<ZZ> tmp;
        ZZ lc_gcd=gcd(F.begin()->second,G.begin()->second );
        Zp lc_gcd_p;
        ZZ Pout_prime;
        std::int64_t Pout_d=INT64_MAX;
        std::int64_t tmp_Pout_d=INT64_MAX;
        
        while (1)
        {
            
            while (F.begin()->second % prime ==0 || G.begin()->second % prime ==0)
            {
                prime=boost::math::prime(++p_index);
            }
            f_p=polynomial_mod(F,prime);
            g_p=polynomial_mod(G,prime);
            lc_gcd_p=Zp(lc_gcd,prime);

            // std::cout<<"p:"<<prime<<std::endl;
            // std::cout<<"f_p:"<<f_p<<std::endl;
            // std::cout<<"g_p:"<<g_p<<std::endl;
            // std::cout<<"lc_gcd_p:"<<lc_gcd_p<<std::endl;

            tmp_Pout_d=__polynomial_GCD(Pout_mod,f_p,g_p,lc_gcd_p,Pout_d);
            if (tmp_Pout_d==-1)
            {
                prime=boost::math::prime(++p_index);
                continue;
            }
            
            // std::cout<<"poly_mod:"<<Pout_mod<<std::endl;
            
            if (tmp_Pout_d < Pout_d)
            {
                Pout_d=tmp_Pout_d;
                Pout_prime=prime;
                // poly_convert(Pout_mod,Pout_);
                Pout_.clear();
                Pout_.reserve(Pout_mod.size());
                for (auto & i:Pout_mod)
                    Pout_.push_back({i.first,i.second.number()});
                

                for (auto &i:Pout_)
                {
                    i.second%=Pout_prime;
                    if (i.second>Pout_prime/2)
                    {
                        i.second-=Pout_prime;
                    }
                }
                
                // std::cout<<"Pout_:"<<Pout_<<std::endl;
                
            }
            else
            {
                Zp tmp_inv(Pout_prime,prime);
                tmp_inv=tmp_inv.inv();
                tmp_Pout_.clear();
                tmp_Pout_.reserve(Pout_.size());
                auto Pout_ptr=Pout_.begin();
                auto Pout_end=Pout_.end();
                auto Pm_ptr=Pout_mod.begin();
                auto Pm_end=Pout_mod.end();
                while (Pout_ptr!=Pout_end && Pm_ptr!=Pm_end)
                {
                    if (F.comp(Pout_ptr->first,Pm_ptr->first))
                    {
                        tmp_Pout_.push_back({Pout_ptr->first,Pout_ptr->second-Pout_ptr->second*tmp_inv.number()*Pout_prime});
                        ++Pout_ptr;
                    }
                    else
                    {
                        if (Pout_ptr->first==Pm_ptr->first)
                        {
                            tmp_Pout_.push_back({std::move(Pm_ptr->first),Pout_ptr->second+
                            (Pm_ptr->second.number()-Pout_ptr->second)*tmp_inv.number()*Pout_prime});
                            ++Pout_ptr;
                        }
                        else
                        {
                            tmp_Pout_.push_back({std::move(Pm_ptr->first),Pm_ptr->second.number()*tmp_inv.number()*Pout_prime});
                            
                        }
                        ++Pm_ptr;
                    }
                    
                }
                while (Pout_ptr!=Pout_end)
                {
                    tmp_Pout_.push_back({std::move(Pout_ptr->first),Pout_ptr->second-Pout_ptr->second*tmp_inv.number()*Pout_prime});
                    ++Pout_ptr;
                }
                while (Pm_ptr!=Pm_end)
                {
                    tmp_Pout_.push_back({std::move(Pm_ptr->first),Pm_ptr->second.number()*tmp_inv.number()*Pout_prime});
                    ++Pm_ptr;
                }
                Pout_prime*=prime; 
                for (auto &i:tmp_Pout_)
                {
                    i.second%=Pout_prime;
                    if (i.second>Pout_prime/2)
                    {
                        i.second-=Pout_prime;
                    }
                    else if (i.second<-Pout_prime/2)
                    {
                        i.second+=Pout_prime;
                    }
                }
                
                // std::cout<<"Pout_:"<<tmp_Pout_<<std::endl;
                
                if (tmp_Pout_==Pout_)
                {
                    if (Pout_d>0)
                    {
                        auto cont_=cont(Pout_);
                        // std::cout<<"cont_:"<<cont_<<std::endl;
                        // pair_vec_div(tmp.data(),tmp_Pout_.data(),cont_.data(),F.comp());
                        tmp=tmp_Pout_/upolynomial_<ZZ>({{0,cont_}});
                        std::swap(tmp,tmp_Pout_);
                        pair_vec_div(tmp.data(),R.data(),F.data(),tmp_Pout_.data(),F.comp());
                        if (R.empty())
                        {
                            pair_vec_div(tmp.data(),R.data(),G.data(),tmp_Pout_.data(),F.comp());
                            if (R.empty())
                            {
                                // pair_vec_multiplies(tmp.data(),tmp_Pout_.data(),cont_gcd.data(),F.comp());
                                tmp=tmp_Pout_*upolynomial_<ZZ>({{0,cont_gcd}});
                                return tmp;
                            }
                        }
                    }
                    else
                    {
                        return {{0,cont_gcd}};
                    }
                }
                swap(tmp_Pout_.data(),Pout_.data());
                       
            }
            prime=boost::math::prime(++p_index);
        }  

    }
   
}
#endif