/**
 * @file realroot.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 实根隔离
*/

#ifndef CLPOLY_REALROOT_HH
#define CLPOLY_REALROOT_HH
#include <clpoly/polynomial.hh>
#include <clpoly/upolynomial.hh>
#include <vector>
namespace clpoly{
    const upolynomial_<ZZ> __upolynomial_x_plus_1={{1,1},{0,1}};
    uint64_t coeffsignchanges(const upolynomial_<ZZ>& G)
    {
        if (G.size()<=1)
            return 0;
        auto ptr1=G.begin();
        auto ptr2=G.begin();
        ++ptr2;
        uint64_t v=0;
        for(;ptr2!=G.end();++ptr1,++ptr2)
        {
            if (sgn(ptr1->second)*sgn(ptr2->second)<0)
                ++v;
        }
        return v;
    }
    upolynomial_<ZZ> _upolynomial_1toinf(const upolynomial_<ZZ>& G)
    {
        if (G.empty())
            return G;
        auto m=G.front().first.deg();
        upolynomial_<ZZ> g,g_;
        for (auto &i:G)
        {
            g_=pow(__upolynomial_x_plus_1,m-i.first.deg());
            for (auto &j:g_)
            {
                j.second*=i.second;
            }
            g=g+g_;
        }
        return g;

    }
    upolynomial_<ZZ> _upolynomial_Bto1(upolynomial_<ZZ> G,const ZZ &B)
    {
        if (G.empty())
            return G;
        auto m=G.front().first.deg();
        for (auto &i:G)
        {
            i.second*=pow(B,i.first.deg());
        }
        return G;

    }
    upolynomial_<ZZ> _upolynomial_01to012(const upolynomial_<ZZ>& G)
    {
        if (G.empty())
            return G;
        auto m=G.front().first.deg();
        upolynomial_<ZZ> g;
        g=G;
        ZZ Z=2;
        for (auto &i:g)
        {
            i.second*=pow(Z,m-i.first.deg());
        }
        return g;
    }
    
    upolynomial_<ZZ> _upolynomial_01to121(const upolynomial_<ZZ>& G)
    {
        if (G.empty())
            return G;
        auto m=G.front().first.deg();
        upolynomial_<ZZ> g,g_;
        ZZ Z=2;
        for (auto &i:G)
        {
            g_=pow(__upolynomial_x_plus_1,i.first.deg());
            ZZ z1=pow(Z,m-i.first.deg());
            for (auto &j:g_)
            {
                j.second*=i.second*z1;
            }
            g=g+g_;
        }
        return g;

    }
    
    // void subcontraction_root_interval(upolynomial_<ZZ> G,QQ& B,QQ & E,const  QQ&  width=0)
    // {
    //     assert(E>B && width>=0);
    //     if (!width || (E-B<=width))
    //         return void();
    //     upolynomial_<ZZ> G_,G1;
    //     while (E-B>width)
    //     {
    //         QQ mid=(B+E)/2;
    //         if (!association<QQ,ZZ,QQ>(G,mid))
    //         {
    //             B=E=mid;
    //             return void();
    //         }
    //         G1=_upolynomial_01to012(G);
    //         G_=_upolynomial_1toinf(G1);
    //         if (coeffsignchanges(G_))
    //         {
    //             E=mid;
    //             G=G1;
    //             continue;
    //         }
    //         G1=_upolynomial_01to121(G);
    //         B=mid;
    //     }
        
        
    // } 
    
    void subuspensky(const upolynomial_<ZZ>& G, std::vector<std::pair<QQ,QQ>>& l,const QQ & B=0,const QQ & E=1)
    {
        // std::cout<<"G:"<<G<<std::endl;
        upolynomial_<ZZ> G_=_upolynomial_1toinf(G);
        uint64_t v=coeffsignchanges(G_);
        // std::cout<<B<<" "<<E<<" "<<v<<std::endl;
        if (v==0)   return void();
        if (v==1 )
        {
            l.push_back({B,E});
            return void();
        }
        QQ mid=(B+E)/2;
        if (!association<QQ,ZZ,QQ>(G,mid))
        {
            l.push_back({mid,mid});
        }
        subuspensky(_upolynomial_01to012(G),l,B,mid);
        subuspensky( _upolynomial_01to121(G),l,mid,E);
    } 

    ZZ RealRootBound(const upolynomial_<ZZ>& G)
    {
        if (G.empty())
            return 0;
        auto ptr=G.begin();
        ZZ B=abs(ptr->second);
        ZZ c;
        for (++ptr;ptr!=G.end();++ptr)
        {
            if((c=abs(ptr->second))>B)
            {
                B=c;
            }
        
        }   
        
        B=B/abs(G.front().second)+3;
        // ++B;
        ZZ C=pow(ZZ(2),sizeinbase(B,2)-1);
        if (B==C)
            return B;
        return  2*C;
    }
    std::vector<std::pair<QQ,QQ>> uspensky(const upolynomial_<ZZ>& G)
    {
        std::vector<std::pair<QQ,QQ>> l,l1,l2;
        if (G.empty()|| G.front().first.empty())
            return l;
        ZZ B=RealRootBound(G);
        std::cout<<B<<std::endl;
        subuspensky(_upolynomial_Bto1(G,-B),l1,0,B);
        subuspensky(_upolynomial_Bto1(G,B),l2,0,B);
        l.reserve(l1.size()+l2.size());
        for (auto i=l1.rbegin();i!=l1.rend();++i)
            l.push_back({-i->second,-i->first});
        for (auto &i:l2)
            l.push_back(i);
        
        return l;
    }


    

}
#endif