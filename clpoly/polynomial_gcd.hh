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


    int64_t _polynomial_GCD(polynomial_<Zp,univariate_priority_order>& Pout,
                            const polynomial_<Zp,univariate_priority_order>& F,
                            const polynomial_<Zp,univariate_priority_order>& G,
                            const polynomial_<Zp,univariate_priority_order>& lc,
                            typename std::list<std::pair<variable,int64_t>>::reverse_iterator vars_,
                            typename std::list<std::pair<variable,int64_t>>::reverse_iterator vars_end, int64_t deg)
    {
        int64_t deg_;
        variable v=vars_->first;
        int64_t v_d=vars_->second+1;
        const univariate_priority_order & comp =F.comp(); 
        polynomial_<Zp,univariate_priority_order> Pout_(&comp);
        polynomial_<Zp,univariate_priority_order> Pout_1(&comp);
        polynomial_<Zp,univariate_priority_order> Pout_2(&comp);
        
        if (vars_==vars_end)
        {
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
            assert(lc.size()==1 && lc.begin()->first.empty());
            Zp lc_inv=(Pout.begin()->second).inv()*lc.begin()->second;
            for (auto &i:Pout)
                i.second*=lc_inv;
            deg_=get_up_deg(Pout);
            if (deg_>=0 && deg_<=deg)
                return deg_;
            else
                return -1;
        }
        else
        {
            int64_t f_d=get_up_deg(F);
            int64_t g_d=get_up_deg(G);
            uint32_t prime=F.begin()->second.prime();
            Zp p_(prime);
            polynomial_<Zp,univariate_priority_order> F_v(&comp);
            polynomial_<Zp,univariate_priority_order> G_v(&comp);
            polynomial_<Zp,univariate_priority_order> lc_v(&comp);
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
                lc_v=association(lc,v,p_);
                //std::cout<<v<<"->"<<p_<<std::endl;
                //std::cout<<"F_v:"<<F_v<<std::endl;
                //std::cout<<"G_v:"<<G_v<<std::endl;
                
                if (get_up_deg(F_v)==f_d && get_up_deg(G_v)==g_d)
                {
                    deg_=_polynomial_GCD(Pout_,F_v,G_v,lc_v,vars_,vars_end,deg);
                    //std::cout<<"v:"<<v<<" deg:"<<deg_<<" GDD:"<<Pout_<<std::endl;
                    if (deg_>=0 && deg_<=deg)
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
                            for(auto &i:_Pout)
                            {
                                //std::cout<<"{";
                                for(auto &j:i.second)
                                {
                                    //std::cout<<j<<",";
                                }
                                //std::cout<<"}";                                
                                //std::cout<<i.first<<"+";
                            }
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
                            basic_monomial<univariate_priority_order> m(&comp);
                            for (auto &i:_Pout)
                            {
                                m.clear();
                                m.reserve(i.first.size());
                                basic_monomial<univariate_priority_order>::iterator m_ptr;
                                for (auto &j:i.first)
                                {
                                    if (j.first==comp.var())
                                        m.push_back({v,0});
                                    m.push_back(j);
                                }
                                if (get_up_deg(m))
                                {
                                    m_ptr=m.end();
                                    --m_ptr;--m_ptr;
                                } 
                                else
                                {
                                    m.push_back({v,0});
                                    m_ptr=m.end();
                                    --m_ptr;
                                }
                                
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
                                            m_ptr->second=j;
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

    template <class Tc,class comp>
    inline polynomial_<Tc,comp> polynomial_GCD(const polynomial_<Tc,comp> &F,const polynomial_<Tc,comp> & G)
    {
        assert(G.comp_ptr()==G.comp_ptr() || G.comp()==F.comp());
        if (F.empty())
            return G;
        if (G.empty())
            return F;
        if (G.size()==1 && G.begin()->first.empty() || F.size()==1 && F.begin()->first.empty())
        {
            Tc cont=0;
            for(auto &i:F)
                cont=gcd(cont,i.second);
            for(auto &i:G)
                cont=gcd(cont,i.second);
            polynomial_<Tc,comp> Pout(F.comp_ptr());
            Pout.push_back({basic_monomial<comp>(F.comp_ptr()),cont});
            return Pout;    
        }


        variable v_;
        int64_t v_d=INT64_MIN;

        for(auto &i:F)
            for (auto &j:i.first)
                if(j.second>v_d)
                {
                    v_=j.first;v_d=j.second;
                }
        
        variable v__;
        int64_t v_d_=INT64_MIN;
        for(auto &i:G)
            for (auto &j:i.first)
                if(j.second>v_d_)
                {
                    v__=j.first;v_d_=j.second;
                }
        if (v_d_<v_d)
        {
            v_=v__;v_d=v_d_;
        }
        univariate_priority_order v_order(v_);
        polynomial_<Tc,univariate_priority_order> F_(&v_order);
        polynomial_<Tc,univariate_priority_order> G_(&v_order);
        poly_convert(F,F_);
        poly_convert(G,G_);
        polynomial_<Tc,univariate_priority_order>  cont(&v_order),cont_(&v_order),lc_gcd(&v_order),tmp(&v_order);
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
                    
                    lc_gcd=tmp;
                    //std::cout<<"lc_gcd="<<lc_gcd<<std::endl;
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
            lc_gcd=tmp;
            cont=std::move(tmp);
        }
        else
        {
            cont=polynomial_GCD(cont,tmp);
        }
        tmp.clear();
        deg=get_up_deg(G_);tmp_deg=deg;
        
        for (auto &i:G_)
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
                    //std::cout<<"lc_gcd="<<lc_gcd<<std::endl;
                    lc_gcd=polynomial_GCD(lc_gcd,tmp);
                    //std::cout<<"lc_gcd="<<lc_gcd<<std::endl;
                }
                cont=polynomial_GCD(cont,tmp);
                tmp_deg=get_up_deg(i.first);
                tmp.clear();
                tmp.push_back(i);
                if (tmp_deg)
                    tmp.back().first.pop_back();
            }
        }
        if (tmp_deg==deg)
        {
            lc_gcd=polynomial_GCD(lc_gcd,tmp);
            
        }
        cont=polynomial_GCD(cont,tmp);
        if (!get_up_deg(F_) || !get_up_deg(G_))
        {
            polynomial_<Tc,comp> Pout(F.comp_ptr());
            poly_convert(cont,Pout);
            return Pout;
        }
        //std::cout<<"F_="<<F_<<std::endl;
        //std::cout<<"G_="<<G_<<std::endl;
        //std::cout<<"cont="<<cont<<std::endl;
        //std::cout<<"lc_gcd="<<lc_gcd<<std::endl;
        
        auto vars=G_.variables();
        auto F_vars=F_.variables();
        auto vars_ptr=vars.begin();
        auto F_vars_ptr=F_vars.begin();
        while (vars_ptr!=vars.end() && F_vars_ptr!=F_vars.end())
        {
            if (v_order(vars_ptr->first,F_vars_ptr->first))
            {
                ++vars_ptr;
            }
            else
            {
                if (vars_ptr->first==F_vars_ptr->first )
                {
                    if ( vars_ptr->second>F_vars_ptr->second)
                        vars_ptr->second=F_vars_ptr->second;
                    ++vars_ptr;++F_vars_ptr;
                }
                else
                {
                    vars.insert(vars_ptr,std::move(*F_vars_ptr));
                    ++F_vars_ptr;
                }
            }
            
        }
        while ( F_vars_ptr!=F_vars.end())
        {
            vars.push_back(std::move(*F_vars_ptr));
            ++F_vars_ptr;
        }

        std::uint32_t tmp_x=vars.size()*v_d;
        if (tmp_x<2) tmp_x=2;
        vars.pop_back();
        std::uint32_t p_index=tmp_x/std::log(tmp_x);
        std::uint32_t prime=boost::math::prime(p_index);
        while (prime <v_d)
        {
            prime=boost::math::prime(++p_index);
        }

        polynomial_<Tc,univariate_priority_order> Pout_(&v_order);
        polynomial_<Tc,univariate_priority_order> tmp_Pout_(&v_order);
        polynomial_<Tc,univariate_priority_order> R(&v_order);
        
        Tc Pout_prime;
        std::int64_t Pout_d=INT64_MAX;
        std::int64_t tmp_Pout_d=INT64_MAX;
        polynomial_<Zp,univariate_priority_order> Pout_mod(&v_order);
        polynomial_<Zp,univariate_priority_order> f_p(&v_order),g_p(&v_order);


        
        while (1)
        {
            
            while (F_.begin()->second % prime ==0 || G_.begin()->second % prime ==0)
            {
                prime=boost::math::prime(++p_index);
            }
            f_p=polynomial_mod(F_,prime);
            g_p=polynomial_mod(G_,prime);
            // std::cout<<"p:"<<prime<<std::endl;
            // std::cout<<"f_p:"<<f_p<<std::endl;
            // std::cout<<"g_p:"<<g_p<<std::endl;
            tmp_Pout_d=_polynomial_GCD(Pout_mod,f_p,g_p,polynomial_mod(lc_gcd,prime),vars.rbegin(),vars.rend(),Pout_d);
            if (tmp_Pout_d<0)
            {
                prime=boost::math::prime(++p_index);
                continue;
            }
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
                //std::cout<<Pout_<<std::endl;
                
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
                    if (i.second>Pout_prime/2)
                    {
                        i.second-=Pout_prime;
                    }
                    else if (i.second<-Pout_prime/2)
                    {
                        i.second+=Pout_prime;
                    }
                }
                //std::cout<<tmp_Pout_<<std::endl;
                if (tmp_Pout_==Pout_)
                {
                    cont_.clear();
                    tmp.clear();
                    tmp_deg=get_up_deg(tmp_Pout_);
                    for (auto &i:tmp_Pout_)
                    {
                        if (get_up_deg(i.first)==tmp_deg)
                        {
                            tmp.push_back(i);
                            if (tmp_deg)
                                tmp.back().first.pop_back();
                        }
                        else
                        {
                            cont_=polynomial_GCD(cont_,tmp);
                            tmp.clear();
                            tmp_deg=get_up_deg(i.first);
                            tmp.push_back(i);
                            if (tmp_deg)
                                tmp.back().first.pop_back();
                        }
                    }
                    cont_=polynomial_GCD(cont_,tmp);
                    //std::cout<<cont_<< std::endl;
                    pair_vec_div(tmp.data(),tmp_Pout_.data(),cont_.data(),v_order);
                    std::swap(tmp,tmp_Pout_);
                    //std::cout<<tmp_Pout_<<std::endl;
                    pair_vec_div(tmp.data(),R.data(),F_.data(),tmp_Pout_.data(),v_order);
                    if (R.empty())
                    {
                        pair_vec_div(tmp.data(),R.data(),G_.data(),tmp_Pout_.data(),v_order);
                        if (R.empty())
                        {
                            polynomial_<Tc,comp> Pout(F.comp_ptr());
                            pair_vec_multiplies(tmp.data(),tmp_Pout_.data(),cont.data(),v_order);
                            poly_convert(tmp,Pout);
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