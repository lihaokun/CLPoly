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
    
    polynomial_<ZZ,lex> cont(const polynomial_<ZZ,lex> &F_);
    polynomial_<ZZ,lex> leadcoeff(const polynomial_<ZZ,lex> &F_);
    int64_t  __polynomial_GCD(       polynomial_<Zp,lex> & Pout,
                            const polynomial_<Zp,lex> & F,
                            const polynomial_<Zp,lex> & G,
                            const polynomial_<Zp,lex> & Lc_gcd,
                            int64_t deg);
    
    polynomial_<ZZ,lex> polynomial_GCD(polynomial_<ZZ,lex>  F, polynomial_<ZZ,lex>  G);
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

    std::vector<std::pair<polynomial_<ZZ,lex>,uint64_t>> squarefree (const polynomial_<ZZ,lex> &  F)
    {
        // std::cout<<"F:"<<F<<std::endl;
        if (is_number(F))
            return {{F,1}};
        auto f_cont=cont(F);
        auto F_=F/f_cont;
        auto lst=squarefree(f_cont);
        polynomial_<ZZ,lex> F_1,F_2,F_3;
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
    
    // void __polynomial_GCD_get_var(
    //             variable& v,
    //             uint64_t& vd,
    //             const polynomial_<Zp,univariate_priority_order>& F)
    // {
    //     // const auto & v1=F.comp().var();
    //     v=variable();vd=0;
    //     // std::pair<variable,uint64_t> v1;
    //     for(auto & i:F)
    //     {
    //         if (get_up_deg(i.first))
    //         {
    //             if (i.first.size()>1)
    //             {
    //                 const auto & v1=i.first[i.first.size()-2];
    //                 if (F.comp(v,v1.first))
    //                 {
    //                     v=v1.first;
    //                     vd=v1.second;
    //                 }
    //                 else if (v=v1.first && vd<v1.second)
    //             }
    //         }
    //         else
    //         {
                
    //         }
            
    //     }
    //     for(auto & i:G)
    //     {
    //         if (!i.first.empty() && F.comp(vd,i.first.back().first))
    //         {
    //             vd=i.first.back().first;
    //         }
    //     }
    // } 
    // int64_t _polynomial_GCD(polynomial_<Zp,univariate_priority_order>& Pout,
    //                         const polynomial_<Zp,univariate_priority_order>& F,
    //                         const polynomial_<Zp,univariate_priority_order>& G,
    //                         const polynomial_<Zp,univariate_priority_order>& lc,
    //                         typename std::list<std::pair<variable,int64_t>>::reverse_iterator vars_,
    //                         typename std::list<std::pair<variable,int64_t>>::reverse_iterator vars_end, int64_t deg)
    // {
    //     int64_t deg_;
    //     variable v=vars_->first;
    //     int64_t v_d=vars_->second+1;
        
    //     const univariate_priority_order & comp =F.comp(); 
    //     polynomial_<Zp,univariate_priority_order> Pout_(&comp);
    //     polynomial_<Zp,univariate_priority_order> Pout_1(&comp);
    //     polynomial_<Zp,univariate_priority_order> Pout_2(&comp);
        
    //     if (vars_==vars_end)
    //     {
    //         Pout=G;
    //         pair_vec_div(Pout_2.data(),Pout_1.data(),F.data(),Pout.data(),comp);
    //         std::cout<<Pout_1<<std::endl;
    //         while(!Pout_1.empty())
    //         {
    //             swap(Pout.data(),Pout_.data());
    //             swap(Pout.data(),Pout_1.data());
    //             pair_vec_div(Pout_2.data(),Pout_1.data(),Pout_.data(),Pout.data(),comp);
    //             std::cout<<Pout_1<<std::endl;
    //         }
    //         assert(lc.size()==1 && lc.begin()->first.empty());
    //         Zp lc_inv=(Pout.begin()->second).inv()*lc.begin()->second;
    //         for (auto &i:Pout)
    //             i.second*=lc_inv;
    //         deg_=get_up_deg(Pout);
    //         if (deg_>=0 && deg_<=deg)
    //             return deg_;
    //         else
    //             return -1;
    //     }
    //     else
    //     {
    //         int64_t f_d=get_up_deg(F);
    //         int64_t g_d=get_up_deg(G);
    //         uint32_t prime=F.begin()->second.prime();
    //         Zp p_(prime);
    //         polynomial_<Zp,univariate_priority_order> F_v(&comp);
    //         polynomial_<Zp,univariate_priority_order> G_v(&comp);
    //         polynomial_<Zp,univariate_priority_order> lc_v(&comp);
    //         int64_t num_s=0;
    //         ++vars_;
    //         std::vector<std::pair<basic_monomial<univariate_priority_order>,std::vector<Zp>>> _Pout,tmp_Pout;
    //         std::vector<Zp> points;
    //         std::vector<Zp> tmp_p;
    //         points.reserve(v_d);
    //         for (int64_t i=0;i<prime;++i)
    //         {
    //             p_.number()=i;
    //             F_v=association(F,v,p_);
    //             G_v=association(G,v,p_);
    //             lc_v=association(lc,v,p_);
    //             std::cout<<v<<"->"<<p_<<std::endl;
    //             std::cout<<"F_v:"<<F_v<<std::endl;
    //             std::cout<<"G_v:"<<G_v<<std::endl;
                
    //             if (get_up_deg(F_v)==f_d && get_up_deg(G_v)==g_d)
    //             {
    //                 deg_=_polynomial_GCD(Pout_,F_v,G_v,lc_v,vars_,vars_end,deg);
    //                 std::cout<<"v:"<<v<<" deg:"<<deg_<<" GDD:"<<Pout_<<std::endl;
    //                 if (deg_>=0 && deg_<=deg)
    //                 {
    //                     if (deg_<deg || !num_s)
    //                     {
    //                         deg=deg_;
    //                         num_s=1;
    //                         _Pout.clear();
    //                         for (auto &i:Pout_)
    //                         {
    //                             _Pout.push_back({std::move(i.first),{i.second}});
    //                         }
    //                         points={p_};
    //                     }
    //                     else
    //                     {
    //                         ++num_s;
    //                         tmp_Pout.clear();
    //                         tmp_Pout.reserve(Pout_.size());
    //                         points.push_back(p_);
    //                         auto _P_ptr=_Pout.begin();
    //                         auto _P_end=_Pout.end();
    //                         auto P_ptr=Pout_.begin();
    //                         auto P_end=Pout_.end();
                            
    //                         while (_P_ptr!=_P_end && P_ptr!=P_end)
    //                         {
    //                             if (comp(P_ptr->first,_P_ptr->first))
    //                             {
    //                                 tmp_p.clear();
    //                                 tmp_p.resize(num_s,Zp(0,prime));
    //                                 tmp_p.back()=P_ptr->second;
    //                                 tmp_Pout.push_back({std::move(P_ptr->first),std::move(tmp_p)});
    //                                 ++P_ptr;
    //                             }
    //                             else{
    //                                 if (P_ptr->first==_P_ptr->first)
    //                                 {
    //                                     tmp_Pout.push_back(std::move(*_P_ptr));
    //                                     tmp_Pout.back().second.push_back(P_ptr->second);
    //                                     ++P_ptr;
    //                                 }else
    //                                 {
    //                                     tmp_Pout.push_back(std::move(*_P_ptr));
    //                                     tmp_Pout.back().second.push_back(Zp(0,prime));
    //                                 }
    //                                 ++_P_ptr;
    //                             }
    //                         }
    //                         while (P_ptr!=P_end)
    //                         {
    //                             tmp_p.clear();
    //                             tmp_p.resize(num_s,Zp(0,prime));
    //                             tmp_p.back()=P_ptr->second;
    //                             tmp_Pout.push_back({std::move(P_ptr->first),std::move(tmp_p)});
    //                             ++P_ptr;
    //                         }
    //                         while (_P_ptr!=_P_end)
    //                         {
    //                             tmp_Pout.push_back(std::move(*_P_ptr));
    //                             tmp_Pout.back().second.push_back(Zp(0,prime));
    //                             ++_P_ptr;
    //                         }  
    //                         swap(tmp_Pout,_Pout);
                            
    //                     }
    //                     if (num_s==v_d)
    //                     {
    //                         for(auto &i:_Pout)
    //                         {
    //                             std::cout<<"{";
    //                             for(auto &j:i.second)
    //                             {
    //                                 std::cout<<j<<",";
    //                             }
    //                             std::cout<<"}";                                
    //                             std::cout<<i.first<<"+";
    //                         }
    //                         std::cout<<std::endl;
    //                         std::vector<Zp> lag;
    //                         lag.resize(v_d*v_d);
    //                         Zp tmp_inv;
    //                         for (int64_t i=0;i<v_d;++i)
    //                         {
    //                             lag[i*v_d]=Zp(1,prime);
    //                             for (int64_t j=1;j<v_d;++j)
    //                                 lag[i*v_d+j]=Zp(0,prime);
    //                             for (int64_t j=0;j<v_d;++j)
    //                                 if (i!=j)
    //                                 {
    //                                     tmp_inv=(points[i]-points[j]).inv();
    //                                     for (int64_t k=v_d-1;k>0;--k)
    //                                     {
    //                                         lag[i*v_d+k]=(lag[i*v_d+k-1]-lag[i*v_d+k]*points[j])*tmp_inv;
    //                                     }
    //                                     lag[i*v_d]=(-lag[i*v_d]*points[j])*tmp_inv;

    //                                 }
                                         
    //                         }
    //                         Pout.clear();
    //                         Pout.reserve(_Pout.size());
    //                         basic_monomial<univariate_priority_order> m(&comp);
    //                         for (auto &i:_Pout)
    //                         {
    //                             m.clear();
    //                             m.reserve(i.first.size());
    //                             basic_monomial<univariate_priority_order>::iterator m_ptr;
    //                             for (auto &j:i.first)
    //                             {
    //                                 if (j.first==comp.var())
    //                                     m.push_back({v,0});
    //                                 m.push_back(j);
    //                             }
    //                             if (get_up_deg(m))
    //                             {
    //                                 m_ptr=m.end();
    //                                 --m_ptr;--m_ptr;
    //                             } 
    //                             else
    //                             {
    //                                 m.push_back({v,0});
    //                                 m_ptr=m.end();
    //                                 --m_ptr;
    //                             }
                                
    //                             for (int64_t j=v_d-1;j>=0;--j)
    //                             {
    //                                 p_.number()=0;
                                   
    //                                 for (int64_t k=0;k<v_d;++k)
    //                                 {
    //                                     p_+=i.second[k]*lag[k*v_d+j];
    //                                 }
    //                                 if (p_)
    //                                 {
    //                                     if (j)
    //                                     {
    //                                         m_ptr->second=j;
    //                                         m.deg()+=j;
    //                                         Pout.push_back({m,p_});
    //                                         m.deg()-=j;
    //                                     }
    //                                     else
    //                                     {
    //                                         Pout.push_back({std::move(i.first),p_});
    //                                     }
                                        
    //                                 }  

    //                             }

    //                         }
    //                         return deg;
    //                     }
                        
                        
    //                 }
                    

    //             }
    //             if (prime<=i+v_d-num_s)
    //             {
    //                 return -1;
    //             }
    //         }
    //         return -1;
    //     }
        

    // }







    polynomial_<ZZ,lex> polynomial_GCD(polynomial_<ZZ,lex>  F, polynomial_<ZZ,lex>  G)
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
            
            return {{m,c}};
            
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
                polynomial_<ZZ,lex> Pout;
                Pout.push_back({basic_monomial<lex>(),cont});
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
        polynomial_<ZZ,lex> F_cont=cont(F);
        polynomial_<ZZ,lex> G_cont=cont(G);
        polynomial_<ZZ,lex> cont_gcd=polynomial_GCD(cont(F),cont(G));

        F=F/F_cont;G=G/G_cont;
        polynomial_<ZZ,lex> lc_gcd=polynomial_GCD(leadcoeff(F),leadcoeff(G));
       
        
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

        polynomial_<ZZ,lex> Pout_,tmp_Pout_,R;
        polynomial_<Zp,lex> Pout_mod,f_p,g_p,lc_gcd_p;
        polynomial_<ZZ,lex> tmp;
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


    // template <class comp>
    // inline polynomial_<ZZ,comp> polynomial_GCD(const polynomial_<ZZ,comp> &F,const polynomial_<ZZ,comp> & G)
    // {
    //     assert(comp_consistent(G.comp(),F.comp()));
    //     if (F.empty())
    //         return G;
    //     if (G.empty())
    //         return F;
    //     if (G.size()==1 && G.begin()->first.empty() || F.size()==1 && F.begin()->first.empty())
    //     {
    //         ZZ cont=0;
    //         for(auto &i:F)
    //             cont=gcd(cont,i.second);
    //         for(auto &i:G)
    //             cont=gcd(cont,i.second);
    //         polynomial_<ZZ,comp> Pout(F.comp_ptr());
    //         Pout.push_back({basic_monomial<comp>(F.comp_ptr()),cont});
    //         return Pout;    
    //     }


    //     variable v_;
    //     int64_t v_d=INT64_MIN;

    //     for(auto &i:F)
    //         for (auto &j:i.first)
    //             if(j.second>v_d)
    //             {
    //                 v_=j.first;v_d=j.second;
    //             }
        
    //     variable v__;
    //     int64_t v_d_=INT64_MIN;
    //     for(auto &i:G)
    //         for (auto &j:i.first)
    //             if(j.second>v_d_)
    //             {
    //                 v__=j.first;v_d_=j.second;
    //             }
    //     if (v_d_<v_d)
    //     {
    //         v_=v__;v_d=v_d_;
    //     }
    //     univariate_priority_order v_order(v_);
    //     polynomial_<ZZ,univariate_priority_order> F_(&v_order);
    //     polynomial_<ZZ,univariate_priority_order> G_(&v_order);
    //     poly_convert(F,F_);
    //     poly_convert(G,G_);
    //     polynomial_<ZZ,univariate_priority_order>  cont(&v_order),cont_(&v_order),lc_gcd(&v_order),tmp(&v_order);
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
                    
    //                 lc_gcd=tmp;
    //                 std::cout<<"lc_gcd="<<lc_gcd<<std::endl;
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
    //         lc_gcd=tmp;
    //         cont=std::move(tmp);
    //     }
    //     else
    //     {
    //         cont=polynomial_GCD(cont,tmp);
    //     }
    //     tmp.clear();
    //     deg=get_up_deg(G_);tmp_deg=deg;
        
    //     for (auto &i:G_)
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
    //                 std::cout<<"lc_gcd="<<lc_gcd<<std::endl;
    //                 lc_gcd=polynomial_GCD(lc_gcd,tmp);
    //                 std::cout<<"lc_gcd="<<lc_gcd<<std::endl;
    //             }
    //             cont=polynomial_GCD(cont,tmp);
    //             tmp_deg=get_up_deg(i.first);
    //             tmp.clear();
    //             tmp.push_back(i);
    //             if (tmp_deg)
    //                 tmp.back().first.pop_back();
    //         }
    //     }
    //     if (tmp_deg==deg)
    //     {
    //         lc_gcd=polynomial_GCD(lc_gcd,tmp);
            
    //     }
    //     cont=polynomial_GCD(cont,tmp);
    //     if (!get_up_deg(F_) || !get_up_deg(G_))
    //     {
    //         polynomial_<ZZ,comp> Pout(F.comp_ptr());
    //         poly_convert(cont,Pout);
    //         return Pout;
    //     }
    //     std::cout<<"F_="<<F_<<std::endl;
    //     std::cout<<"G_="<<G_<<std::endl;
    //     std::cout<<"cont="<<cont<<std::endl;
    //     std::cout<<"lc_gcd="<<lc_gcd<<std::endl;
        
    //     auto vars=G_.variables();
    //     auto F_vars=F_.variables();
    //     auto vars_ptr=vars.begin();
    //     auto F_vars_ptr=F_vars.begin();
    //     while (vars_ptr!=vars.end() && F_vars_ptr!=F_vars.end())
    //     {
    //         if (v_order(vars_ptr->first,F_vars_ptr->first))
    //         {
    //             ++vars_ptr;
    //         }
    //         else
    //         {
    //             if (vars_ptr->first==F_vars_ptr->first )
    //             {
    //                 if ( vars_ptr->second>F_vars_ptr->second)
    //                     vars_ptr->second=F_vars_ptr->second;
    //                 ++vars_ptr;++F_vars_ptr;
    //             }
    //             else
    //             {
    //                 vars.insert(vars_ptr,std::move(*F_vars_ptr));
    //                 ++F_vars_ptr;
    //             }
    //         }
            
    //     }
    //     while ( F_vars_ptr!=F_vars.end())
    //     {
    //         vars.push_back(std::move(*F_vars_ptr));
    //         ++F_vars_ptr;
    //     }

    //     std::uint32_t tmp_x=vars.size()*v_d;
    //     if (tmp_x<2) tmp_x=2;
    //     vars.pop_back();
    //     std::uint32_t p_index=tmp_x/std::log(tmp_x);
    //     std::uint32_t prime=boost::math::prime(p_index);
    //     while (prime <v_d)
    //     {
    //         prime=boost::math::prime(++p_index);
    //     }

    //     polynomial_<ZZ,univariate_priority_order> Pout_(&v_order);
    //     polynomial_<ZZ,univariate_priority_order> tmp_Pout_(&v_order);
    //     polynomial_<ZZ,univariate_priority_order> R(&v_order);
        
    //     ZZ Pout_prime;
    //     std::int64_t Pout_d=INT64_MAX;
    //     std::int64_t tmp_Pout_d=INT64_MAX;
    //     polynomial_<Zp,univariate_priority_order> Pout_mod(&v_order);
    //     polynomial_<Zp,univariate_priority_order> f_p(&v_order),g_p(&v_order);


        
    //     while (1)
    //     {
            
    //         while (F_.begin()->second % prime ==0 || G_.begin()->second % prime ==0)
    //         {
    //             prime=boost::math::prime(++p_index);
    //         }
    //         f_p=polynomial_mod(F_,prime);
    //         g_p=polynomial_mod(G_,prime);
    //         // std::cout<<"p:"<<prime<<std::endl;
    //         // std::cout<<"f_p:"<<f_p<<std::endl;
    //         // std::cout<<"g_p:"<<g_p<<std::endl;
    //         tmp_Pout_d=_polynomial_GCD(Pout_mod,f_p,g_p,polynomial_mod(lc_gcd,prime),vars.rbegin(),vars.rend(),Pout_d);
    //         if (tmp_Pout_d<0)
    //         {
    //             prime=boost::math::prime(++p_index);
    //             continue;
    //         }
    //         if (tmp_Pout_d < Pout_d)
    //         {
    //             Pout_d=tmp_Pout_d;
    //             Pout_prime=prime;
    //             poly_convert(Pout_mod,Pout_);
    //             for (auto &i:Pout_)
    //             {
    //                 i.second%=Pout_prime;
    //                 if (i.second>Pout_prime/2)
    //                 {
    //                     i.second-=Pout_prime;
    //                 }
    //             }
    //             std::cout<<Pout_<<std::endl;
                
    //         }
    //         else
    //         {
    //             Zp tmp_inv(Pout_prime,prime);
    //             tmp_inv=tmp_inv.inv();
    //             tmp_Pout_.clear();
    //             tmp_Pout_.reserve(Pout_.size());
    //             auto Pout_ptr=Pout_.begin();
    //             auto Pout_end=Pout_.end();
    //             auto Pm_ptr=Pout_mod.begin();
    //             auto Pm_end=Pout_mod.end();
    //             while (Pout_ptr!=Pout_end && Pm_ptr!=Pm_end)
    //             {
    //                 if (v_order(Pout_ptr->first,Pm_ptr->first))
    //                 {
    //                     tmp_Pout_.push_back({Pout_ptr->first,Pout_ptr->second-Pout_ptr->second*tmp_inv.number()*Pout_prime});
    //                     ++Pout_ptr;
    //                 }
    //                 else
    //                 {
    //                     if (Pout_ptr->first==Pm_ptr->first)
    //                     {
    //                         tmp_Pout_.push_back({std::move(Pm_ptr->first),Pout_ptr->second+
    //                         (Pm_ptr->second.number()-Pout_ptr->second)*tmp_inv.number()*Pout_prime});
    //                         ++Pout_ptr;
    //                     }
    //                     else
    //                     {
    //                         tmp_Pout_.push_back({std::move(Pm_ptr->first),Pm_ptr->second.number()*tmp_inv.number()*Pout_prime});
                            
    //                     }
    //                     ++Pm_ptr;
    //                 }
                    
    //             }
    //             while (Pout_ptr!=Pout_end)
    //             {
    //                 tmp_Pout_.push_back({std::move(Pout_ptr->first),Pout_ptr->second-Pout_ptr->second*tmp_inv.number()*Pout_prime});
    //                 ++Pout_ptr;
    //             }
    //             while (Pm_ptr!=Pm_end)
    //             {
    //                 tmp_Pout_.push_back({std::move(Pm_ptr->first),Pm_ptr->second.number()*tmp_inv.number()*Pout_prime});
    //                 ++Pm_ptr;
    //             }
    //             Pout_prime*=prime; 
    //             for (auto &i:tmp_Pout_)
    //             {
    //                 i.second%=Pout_prime;
    //                 if (i.second>Pout_prime/2)
    //                 {
    //                     i.second-=Pout_prime;
    //                 }
    //                 else if (i.second<-Pout_prime/2)
    //                 {
    //                     i.second+=Pout_prime;
    //                 }
    //             }
    //             std::cout<<tmp_Pout_<<std::endl;
    //             if (tmp_Pout_==Pout_)
    //             {
    //                 cont_.clear();
    //                 tmp.clear();
    //                 tmp_deg=get_up_deg(tmp_Pout_);
    //                 for (auto &i:tmp_Pout_)
    //                 {
    //                     if (get_up_deg(i.first)==tmp_deg)
    //                     {
    //                         tmp.push_back(i);
    //                         if (tmp_deg)
    //                             tmp.back().first.pop_back();
    //                     }
    //                     else
    //                     {
    //                         cont_=polynomial_GCD(cont_,tmp);
    //                         tmp.clear();
    //                         tmp_deg=get_up_deg(i.first);
    //                         tmp.push_back(i);
    //                         if (tmp_deg)
    //                             tmp.back().first.pop_back();
    //                     }
    //                 }
    //                 cont_=polynomial_GCD(cont_,tmp);
    //                 std::cout<<cont_<< std::endl;
    //                 pair_vec_div(tmp.data(),tmp_Pout_.data(),cont_.data(),v_order);
    //                 std::swap(tmp,tmp_Pout_);
    //                 std::cout<<tmp_Pout_<<std::endl;
    //                 pair_vec_div(tmp.data(),R.data(),F_.data(),tmp_Pout_.data(),v_order);
    //                 if (R.empty())
    //                 {
    //                     pair_vec_div(tmp.data(),R.data(),G_.data(),tmp_Pout_.data(),v_order);
    //                     if (R.empty())
    //                     {
    //                         polynomial_<ZZ,comp> Pout(F.comp_ptr());
    //                         pair_vec_multiplies(tmp.data(),tmp_Pout_.data(),cont.data(),v_order);
    //                         poly_convert(tmp,Pout);
    //                         return Pout;
    //                     }
    //                 }
    //             }
    //             swap(tmp_Pout_.data(),Pout_.data());
                       
    //         }
    //         prime=boost::math::prime(++p_index);
    //     }  

        
    // }
    polynomial_<ZZ,lex> leadcoeff(const polynomial_<ZZ,lex> &F_)
    {
        polynomial_<ZZ,lex>  lc;
        auto v=get_first_var(F_);
        int64_t deg=get_first_deg(F_);
        basic_monomial<lex> m;
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
    polynomial_<ZZ,lex> cont(const polynomial_<ZZ,lex> &F_)
    {
        // std::cout<<F_<<std::endl;
        polynomial_<ZZ,lex>  cont, tmp;
        auto v=get_first_var(F_);
        int64_t deg=get_first_deg(F_);
        int64_t tmp_deg=deg;
        basic_monomial<lex> m;
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
    polynomial_<ZZ,univariate_priority_order> cont(const polynomial_<ZZ,univariate_priority_order> &F_)
    {
        auto & v_order=F_.comp();
        polynomial_<ZZ,univariate_priority_order>  cont(&v_order),
                                                   cont_(&v_order),
                                                   tmp(&v_order);
        int64_t deg=get_up_deg(F_);
        int64_t tmp_deg=deg;
        for (auto &i:F_)
        {
            if (get_up_deg(i.first)==tmp_deg)
            {
                tmp.push_back(i);
                if (tmp_deg)
                    tmp.back().first.pop_back();
            }
            else{
                if (tmp_deg==deg)
                {
                    cont=std::move(tmp);
                   
                }
                else
                {
                    cont=polynomial_GCD(cont,tmp);
                    
                }
                tmp_deg=get_up_deg(i.first);
                tmp.clear();
                tmp.push_back(i);
                if (tmp_deg)
                    tmp.back().first.pop_back();
            }
        }
        if (tmp_deg==deg)
        {
            cont=std::move(tmp);
        }
        else
        {
            cont=polynomial_GCD(cont,tmp);
        }
        return cont;
    } 
    

    int64_t  __polynomial_GCD(       polynomial_<Zp,lex> & Pout,
                            const polynomial_<Zp,lex> & F,
                            const polynomial_<Zp,lex> & G,
                            const polynomial_<Zp,lex> & Lc_gcd,
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
        polynomial_<Zp,lex> Pout_;
        polynomial_<Zp,lex> Pout_1;
        polynomial_<Zp,lex> Pout_2;
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
            Zp lc_inv=Pout.front().second.inv()*Lc_gcd.begin()->second;;
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
            polynomial_<Zp,lex> F_v;
            polynomial_<Zp,lex> G_v;
            polynomial_<Zp,lex> lc_v;
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
            std::vector<std::pair<basic_monomial<lex>,std::vector<Zp>>> _Pout,tmp_Pout;
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
                            basic_monomial<lex> m;
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
}
#endif