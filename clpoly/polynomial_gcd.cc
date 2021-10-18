/**
 * @file polynomial_gcd.cc
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义polynomial_gcd
*/
#include<clpoly/polynomial_gcd.hh>
namespace clpoly{ 
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
                std::swap(tmp_Pout_.data(),Pout_.data());
                       
            }
            prime=boost::math::prime(++p_index);
        }  

    }
}