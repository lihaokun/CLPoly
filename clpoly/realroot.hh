/*
Module Name:
    realroot.hh
Abstract:
    实根隔离
Author:
    haokun li
Notes:
*/

#ifndef CLPOLY_REALROOT_HH
#define CLPOLY_REALROOT_HH
#include <clpoly/polynomial.hh>
#include <clpoly/upolynomial.hh>
#include <vector>
namespace clpoly{
    upolynomial_<ZZ> __upolynomial_x_plus_1={{1,1},{0,1}};
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
    
    void subcontraction_root_interval(upolynomial_<ZZ> G,QQ& B,QQ & E,const  QQ&  width=0)
    {
        assert(E>B && width>=0);
        if (!width || (E-B<=width))
            return void();
        upolynomial_<ZZ> G_,G1;
        while (E-B>width)
        {
            QQ mid=(B+E)/2;
            if (!association<QQ,ZZ,QQ>(G,mid))
            {
                B=E=mid;
                return void();
            }
            G1=_upolynomial_01to012(G);
            G_=_upolynomial_1toinf(G1);
            if (coeffsignchanges(G_))
            {
                E=mid;
                G=G1;
                continue;
            }
            G1=_upolynomial_01to121(G);
            B=mid;
        }
        
        
    } 
    
    void subuspensky(const upolynomial_<ZZ>& G, std::vector<std::pair<QQ,QQ>>& l,const QQ & B=0,const QQ & E=1,const  QQ&  width=0)
    {
        upolynomial_<ZZ> G_=_upolynomial_1toinf(G);
        uint64_t v=coeffsignchanges(G_);
        if (v==0)   return void();
        if (v==1 )
        {
            if (!width || (E-B<=width))
            {
                l.push_back({B,E});
                return void();
            }
            else{
                auto B1=B;
                auto E1=E;
                subcontraction_root_interval(G,B1, E1,width);
                l.push_back({B1,E1});
                return void();
            }
            
        }
        QQ mid=(B+E)/2;
        if (!association<QQ,ZZ,QQ>(G,mid))
        {
            l.push_back({mid,mid});
        }
        subuspensky(_upolynomial_01to012(G),l,B,mid,width);
        subuspensky( _upolynomial_01to121(G),l,mid,E,width);
    } 
    
 

}
#endif