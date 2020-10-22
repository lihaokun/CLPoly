/*
Module Name:
    random.hh
Abstract:
Author:
    haokun li
Notes:
*/

#ifndef CLPOLY_RANDOM_HH
#define CLPOLY_RANDOM_HH

#include <clpoly/polynomial.hh>
#include <list>
#include <string>
#include <random>


namespace clpoly{   
    template<class Tc>
    polynomial_<Tc> random_polynomial(const std::vector<variable> & v,uint64_t deg,double p,int up,int down)
    {
        if (down > up) std::swap(up,down);
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

    template <class T>
    std::vector<T> RandomSample(const  std::vector<T> & l,uint64_t n)
    {
        
        std::random_device rd; 
        std::mt19937 gen(rd());
        auto size=l.size();
        std::vector<T>  lout;
        lout.reserve(std::max(n,size));
        std::vector<bool> bl;
        bl.reserve(size);
        uint64_t select;
        for (uint64_t i=size;i>0;--i)
            bl.push_back(true);
        for (uint64_t i=std::max(n,size);i>0;--i)
        {
            if (i==1)
            {
                select = 1;
            }
            else{
                std::uniform_int_distribution<uint64_t> dis(1, i);
                select=dis(gen);
            }
            for (uint64_t j=0;j<size;++j)
            {
                if (bl[j])
                {
                    --select;
                    if (!select)
                    {
                        bl[j]=false;
                        lout.push_back(l[j]);
                        break;
                    }
                }
            }
        }
        return lout;
    }
}
#endif