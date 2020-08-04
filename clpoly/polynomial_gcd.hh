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

namespace clpoly{    

    template <class compare>
    polynomial_<Zp,compare> polynomial_mod(const polynomial_<ZZ,compare> & p, uint32_t prime)
    {
        polynomial_<Zp,compare> new_p(p.comp_ptr());
        Zp coeff(prime);
        for (auto & i:p)
        {
            coeff=i.second; 
            if (coeff)
                new_p.push_back({i.first,std::move(coeff)});
        }
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

    int64_t _polynomial_GCD(polynomial_<Zp,univariate_priority_order>& Pout,
                            const polynomial_<Zp,univariate_priority_order>& F,
                            const polynomial_<Zp,univariate_priority_order>& G,
                            typename std::list<std::pair<variable,int64_t>>::iterator vars_,
                            const variable& v_, int64_t deg)
    {
        int64_t deg_;
        variable v=vars_->first;
        int64_t v_d=vars_->second+1;
        polynomial_<Zp,univariate_priority_order> Pout_;
        polynomial_<Zp,univariate_priority_order> Pout_1;
        polynomial_<Zp,univariate_priority_order> Pout_2;
        const univariate_priority_order & comp =F.comp(); 
        if (v==v_)
        {
            //deg_=_upolynomial_GCD(Pout,F,G);
            Pout=G;
            pair_vec_div(Pout_2.data(),Pout_1.data(),F.data(),Pout.data(),comp);
            //std::cout<<Pout_1<<std::endl;
            while(!Pout_1.empty())
            {
                swap(Pout.data(),Pout_.data());
                swap(Pout.data(),Pout_1.data());
                pair_vec_div(Pout_2.data(),Pout_1.data(),Pout_.data(),Pout.data(),comp);
                //std::cout<<Pout_1<<std::endl;
            }
            deg_=get_up_deg(Pout);
            if (deg_!=0 && deg_<=deg)
                return deg_;
            else
                return 0;
        }
        else
        {
            int64_t f_d=get_up_deg(F);
            int64_t g_d=get_up_deg(G);
            uint32_t prime=F.begin()->second.prime();
            Zp p_(prime);
            polynomial_<Zp,univariate_priority_order> F_v;
            polynomial_<Zp,univariate_priority_order> G_v;
            
            int64_t num_s=0;
            ++vars_;
            std::vector<std::pair<basic_monomial<univariate_priority_order>,std::vector<Zp>>> _Pout,tmp_Pout;
            std::vector<Zp> points;
            std::vector<Zp> tmp_p;
            points.reserve(v_d);
            for (int64_t i=0;i<prime;++i)
            {
                p_.number()=i;
                F_v=association(F,v,p_);
                G_v=association(G,v,p_);
                if (get_up_deg(F_v)==f_d && get_up_deg(G_v)==g_d)
                {
                    deg_=_polynomial_GCD(Pout_,F_v,G_v,vars_,v_,deg);
                    if (deg_!=0 && deg_<=deg)
                    {
                        if (deg_<deg)
                        {
                            deg=deg_;
                            num_s=1;
                            _Pout.clear();
                            for (auto &i:Pout_)
                            {
                                _Pout.push_back({std::move(i.first),{i.second}});
                            }
                            //points.clear();
                            points={p_};
                        }
                        else
                        {
                            ++num_s;
                            tmp_Pout.clear();
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
                                    tmp_p.resize(num_s);
                                    tmp_p.back()=P_ptr->second;
                                    tmp_Pout.push_back({std::move(P_ptr->first),std::move(tmp_p)});
                                    ++P_ptr;
                                }
                                else{
                                    tmp_Pout.push_back(std::move(*_P_ptr));
                                    if (P_ptr->first==_P_ptr->first)
                                    {
                                        tmp_Pout.back().second.push_back(P_ptr->second);
                                        ++P_ptr;
                                    }else
                                    {
                                        tmp_Pout.back().second.push_back(Zp(0,prime));
                                    }
                                    ++_P_ptr;
                                }
                            }
                            while (P_ptr!=P_end)
                            {
                                tmp_p.clear();
                                tmp_p.resize(num_s);
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
                            std::vector<Zp> lag;
                            lag.resize(v_d*v_d);
                            Zp tmp_inv;
                            for (int64_t i=0;i<v_d;++i)
                            {
                                lag[i*v_d]=Zp(1,prime);
                                
                                for (int64_t j=0;j<v_d;++j)
                                    if (i!=j)
                                    {
                                        tmp_inv=(points[i]-points[j]).inv();
                                        for (int64_t k=v_d-1;k>0;++k)
                                        {
                                            lag[i*v_d+k]=(lag[i*v_d+k-1]-lag[i*v_d+k]*points[j])*tmp_inv;
                                        }
                                        lag[i*v_d]=(-lag[i*v_d]*points[j])*tmp_inv;

                                    }
                                         
                            }
                            Pout.clear();
                            Pout.resize(tmp_Pout.size());
                            basic_monomial<univariate_priority_order> m(&comp);
                            for (auto &i:tmp_Pout)
                            {
                                for (int64_t j=v_d-1;j>=0;++j)
                                {
                                    p_.number()=0;
                                    m=i.first;
                                    m.data().insert(m.data().begin(),{v,0});
                                    for (int64_t k=0;k<v_d;++k)
                                    {
                                        p_+=i.second[k]*lag[k*v_d+j];
                                    }
                                    if (p_)
                                    {
                                        if (j)
                                        {
                                            m.begin()->second=j;
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
                if (prime<=i+num_s)
                {
                    return 0;
                }
            }
            return 0;
        }
        

    }

    template <class Tc,class comp>
    inline polynomial_<Tc,comp> polynomial_GCD(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F)
    {
        assert(G.comp_ptr()==G.comp_ptr() || G.comp()==F.comp());
        // polynomial_<Tc,comp> G_lc(F.comp_ptr());
        // polynomial_<Tc,comp> F_lc(F.comp_ptr());
        if (F.empty() || G.empty())
            return {};
        
        Tc tmp;

        //auto vars=G.variables();
        if (G.size()==1 && G.begin()->first.empty())
        {
            tmp=G.begin()->second;
            for (auto & i:F)
                tmp=gcd(i.second,tmp);
            return {{{},tmp}};
        }
        //auto F_vars=F.variables();
        if (F.size()==1 && F.begin()->first.empty())
        {
            tmp=F.begin()->second;
            for (auto & i:G)
                tmp=gcd(i.second,tmp);
            return {{{},tmp}};
        }
        
        // _variables_pair_marge(F_vars.begin(),F_vars.end(),vars,G.comp());

        

        variable v_;
        int64_t v_d=INT64_MIN;

        for(auto &i:F)
            for (auto &j:i.first)
                if(j.second>v_d)
                {
                    v_=j.first;v_d=j.second;
                }
        
        for(auto &i:G)
            for (auto &j:i.first)
                if(j.second>v_d)
                {
                    v_=j.first;v_d=j.second;
                }
        
        univariate_priority_order v_order(v_);
        polynomial_<Tc,univariate_priority_order> F_(&v_order);
        polynomial_<Tc,univariate_priority_order> G_(&v_order);
        poly_convert(F,F_);
        poly_convert(G,G_);


        // int64_t f_d=get_up_deg(F_);
        // int64_t g_d=get_up_deg(G_);
        // polynomial_<Tc,univariate_priority_order> F_lc(&v_order);
        // polynomial_<Tc,univariate_priority_order> G_lc(&v_order);
        // for(auto &i:F_)
        //     if (get_up_deg(i.first)==f_d)
        //     {
        //         F_lc.push_back(i);
        //         F_lc.first.pop_back();
        //     }
        //     else break;
        // for(auto &i:G_)
        //     if (get_up_deg(i.first)==g_d)
        //     {
        //         G_lc.push_back(i);
        //         G_lc.first.pop_back();
        //     }
        //     else break;
        

        auto vars=G.variables();
        auto F_vars=F.variables();
        _variables_pair_marge(F_vars.begin(),F_vars.end(),vars,G.comp());

        std::uint32_t tmp_x=vars.size()*v_d;
        std::int32_t p_index=2;//tmp_x/std::log(tmp_x);
        std::uint32_t prime=boost::math::prime(p_index);
        while (prime <v_d)
        {
            prime=boost::math::prime(++p_index);
        }

        polynomial_<Tc,univariate_priority_order> Pout_(&v_order);
        polynomial_<Tc,univariate_priority_order> tmp_Pout_(&v_order);
        Tc Pout_prime;
        std::int64_t Pout_d=INT64_MAX;
        std::int64_t tmp_Pout_d=INT64_MAX;
        polynomial_<Zp,univariate_priority_order> Pout_mod(&v_order);
        polynomial_<Zp,univariate_priority_order> f_p(&v_order),g_p(&v_order);
        Zp tmp_inv;

        while (1)
        {
            
            while (F_.begin()->second % prime ==0 || G_.begin()->second % prime ==0)
            {
                prime=boost::math::prime(++p_index);
            }
            f_p=polynomial_mod(F_,prime);
            g_p=polynomial_mod(G_,prime);
            std::cout<<f_p<<std::endl;
            std::cout<<g_p<<std::endl;
            while (!(tmp_Pout_d=_polynomial_GCD(Pout_mod,f_p,g_p,vars.begin(),v_,Pout_d)))
            {
                prime=boost::math::prime(++p_index);
                f_p=polynomial_mod(F_,prime);
                g_p=polynomial_mod(G_,prime);
                std::cout<<f_p<<std::endl;
                std::cout<<g_p<<std::endl;
            }
            std::cout<<Pout_mod<<std::endl;
            if (tmp_Pout_d < Pout_d)
            {
                Pout_d=tmp_Pout_d;
                Pout_prime=prime;
                poly_convert(Pout_mod,Pout_);
                std::cout<<Pout_<<std::endl;
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
                    if (v_order(Pout_ptr->first,Pm_ptr->first))
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
                }
                std::cout<<tmp_Pout_<<std::endl;
                if (tmp_Pout_==Pout_)
                {
                    for (auto &i:tmp_Pout_)
                    {
                        if (i.second>Pout_prime/2)
                        {
                            i.second-=Pout_prime;
                            i.second%=Pout_prime;
                        }
                    }
                    std::cout<<tmp_Pout_<<std::endl;
                    polynomial_<Tc,univariate_priority_order> R(&v_order);
                    pair_vec_div(tmp_Pout_.data(),R.data(),F_.data(),Pout_.data(),v_order);
                    if (R.empty())
                    {
                        pair_vec_div(tmp_Pout_.data(),R.data(),G_.data(),Pout_.data(),v_order);
                        if (R.empty())
                        {
                            polynomial_<Tc,comp> Pout(F.comp_ptr());
                            poly_convert(Pout_,Pout);
                            return Pout;
                        }
                    }
                }
                swap(tmp_Pout_.data(),Pout_.data());
                       
            }
            prime=boost::math::prime(++p_index);
        }  

        
    }

}
#endif