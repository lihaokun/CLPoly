/**
 * @file random.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 
 * 
 */

#ifndef CLPOLY_RANDOM_HH
#define CLPOLY_RANDOM_HH

#include <clpoly/polynomial.hh>
#include <clpoly/associatedgraph.hh>
#include <list>
#include <string>
#include <random>
#include <algorithm>


namespace clpoly{
    inline std::vector<size_t> random_select(size_t m,size_t n,bool no_repeat=true) // 0,..,m-1 中选 n 个
    {
        std::vector<size_t> v;
        if (no_repeat && n>=m)
        {
            n=m;
            for(size_t i=0;i!=n;++i)
                v.push_back(i);
            return v;
        
        }
        std::random_device rd; 
        std::mt19937 gen(rd());
        size_t tmp=no_repeat?m-n:m-1;
        std::uniform_int_distribution<> dis(0, tmp);
        for (size_t i=0;i!=n;++i)
        {
            v.push_back(dis(gen));
        }
        std::sort(v.begin(),v.end());
        if (no_repeat)
        {
            for (size_t i=0;i!=n;++i)
            {
                v[i]+=i;
            }
        }
        return v;
    }
    
    template<class Tc>
    polynomial_<Tc> random_polynomial(const std::vector<variable> & v,uint64_t deg,uint64_t len,std::pair<int,int>  coeff,bool is_add_num=false)
    {
        if (coeff.first > coeff.second) std::swap(coeff.first,coeff.second);
        std::vector<std::pair<monomial,Tc>> p1;
        std::random_device rd; 
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(coeff.first, coeff.second);
        std::uniform_int_distribution<> deg_r(0, deg);
        std::vector<std::pair<variable,int64_t>>  m;
        for (auto &i:v)
            m.emplace_back(i,0);
        {
            auto tmp=random_select(deg+1,v.size()-1,false);
            tmp.push_back(deg);
            size_t deg_=0;
            for (size_t i=0;i!=v.size();++i)
            {
                m[i].second=tmp[i]-deg_;
                deg_=tmp[i];
            }
            p1.emplace_back(monomial(m),Tc(dis(gen)));
        }
        for (size_t i=0;i<len-1;++i)
        {
            auto tmp=random_select(deg+1,v.size(),false);
            size_t deg_=0;
            for (size_t i=0;i!=v.size();++i)
            {
                m[i].second=tmp[i]-deg_;
                deg_=tmp[i];
            }
            p1.emplace_back(monomial(m),Tc(dis(gen)));
        }
        if (is_add_num)
            p1.push_back({{},Tc(dis(gen))});
        return polynomial_<Tc>(p1);
    }

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
                if (mp(gen))
                {
                    p1.push_back({{},Tc(dis(gen))});
                }
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

    template <class node>
    graph<node> random_graph(const std::vector<node> & nodes,double p)
    {
        graph<node> G;
        for (auto & i:nodes)
            G.add_node(i);
        size_t n=G.size();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::bernoulli_distribution d(p);
        for (size_t i=0;i<n;++i)
        {
            for (size_t j=i+1;j<n;++j)
                if (d(gen))
                    G.add_edge_index(i,j);
        }
        return G;
    }
    /**
     * random_polynomials<Tc>
     * 生成一组多项式
     * 变量相关图的稀疏性是 p1
     * 项稀疏性是 p2^std::max(n-1,1) n为这个多项式变量数
    **/
    template<class Tc>
    std::vector<polynomial_<Tc>> random_polynomials(const std::vector<variable> & vars,uint64_t deg,double p1,double p2,int up,int down)
    {
        auto G=random_graph<variable>(vars,p1);
        // std::cout<<G<<std::endl;
        std::vector<polynomial_<Tc>> l;
        std::vector<variable> v_l;
        
        for (auto &v:G.nodes())
        {
            v_l.clear();
            v_l.push_back(v);
            auto i=G.index(v);
            double p=1;
            for (auto &j:G.adjacency_list()[i])
            {
                if (j!=i)
                {
                    bool b=true;
                    auto node_j=G.node(j);
                    for (auto &k:v_l)
                    {
                        if (!G.is_edge(k,node_j))
                        {
                            b=false;
                            break;
                        }
                    }
                    if (b)
                    {
                        v_l.push_back(node_j);
                        p*=p2;
                    }
                }
            }
            if(v_l.size()>1)
                l.push_back(random_polynomial<clpoly::ZZ>(v_l,deg,p,up,down));
        }
        return l;
    }
}
#endif