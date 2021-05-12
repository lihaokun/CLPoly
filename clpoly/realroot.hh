/**
 * @file realroot.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 实根隔离
*/

#ifndef CLPOLY_REALROOT_HH
#define CLPOLY_REALROOT_HH
#include <clpoly/polynomial.hh>
#include <clpoly/upolynomial.hh>
#include <clpoly/polynomial_gcd.hh>
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

    class uroot
    {
        public:
            std::vector<upolynomial_ZZ>* upolys;
            std::map<upolynomial_ZZ,size_t>* upolymap;
            size_t poly_index;
            QQ left;
            QQ right;
            uroot(size_t _poly_index,QQ _l,QQ _r,std::vector<upolynomial_ZZ>* _upolys,std::map<upolynomial_ZZ,size_t>* _upolymap)
            :poly_index(_poly_index),left(_l),right(_r),upolys(_upolys),upolymap(_upolymap){}
            const upolynomial_ZZ & poly()
            {
                return (*upolys)[poly_index];
            }
            int64_t add_poly(clpoly::upolynomial_ZZ  p)
            {
                assert(!p.empty());
                int64_t s=1;
                if (p.front().second<0)
                {
                    p=-p;
                    s=-1;
                }
                auto it=this->upolymap->find(p);
                size_t index;
                if (it==this->upolymap->end())
                {
                    index=upolys->size();
                    (*upolymap)[p]=upolys->size();
                    upolys->push_back(std::move(p));
                }
                else
                    index=it->second;
                return index*s;
            }
            
    };
    upolynomial_<ZZ> _upolynomial_Rtoab(const upolynomial_<ZZ>& G,const QQ &a,const QQ& b)
    {
        if (G.empty())
            return G;
        ZZ den=lcm(a.get_den(),b.get_den());
        ZZ _gcd=gcd(a.get_den(),b.get_den());
        
        ZZ a1=a.get_num()*b.get_den()/_gcd;
        ZZ b1=b.get_num()*a.get_den()/_gcd;
        auto m=G.front().first.deg();
        upolynomial_<ZZ>  pnum={{1,a1},{0,b1}};
        upolynomial_<ZZ>  pden={{1,den},{0,den}};
        upolynomial_<ZZ> g,g_,g_1;
        for (auto &i:G)
        {
            g_=pow(pden,m-i.first.deg())*pow(pnum,i.first.deg());
            for (auto &j:g_)
            {
                j.second*=i.second;
            }
            g=g+g_;
        }
        return g;

    }
    void _uroot_check(uroot* r,const QQ & mid)
    {
        if (mid>=r->left || mid <=r->right)
            return void();
        QQ rig_ass=association<QQ,ZZ,QQ>(r->poly(),r->right);
        QQ mid_ass=association<QQ,ZZ,QQ>(r->poly(),mid);
        QQ left_ass=association<QQ,ZZ,QQ>(r->poly(),r->left);
        if (mid_ass==0)
        {
            r->left=r->right=mid;
            return void();
        }
        if (rig_ass==0)
        {
            if (left_ass==0) // 不应该进入
            {
                auto p=_upolynomial_Rtoab(r->poly(),r->left,mid);
                uint64_t v=coeffsignchanges(p);
                if (v==0)
                {
                    r->left=mid;return void();
                }
                r->right=mid;return void();
            }
            else
            {
                if (sgn(mid_ass)==sgn(left_ass))
                {
                    r->right=mid;return void();
                }
                r->left=mid;return void();
            }
        }
        if (sgn(mid_ass)==sgn(rig_ass))
        {
            r->left=mid;return void();  
        }
        r->right=mid;return void();
    }
    bool _uroot_check(const upolynomial_<ZZ>& G,const QQ &a,const QQ& b)
    {
        QQ rig_ass=association<QQ,ZZ,QQ>(G,b);
        if (a==b && rig_ass==0)
            return true;
        QQ left_ass=association<QQ,ZZ,QQ>(G,a);
        if (sgn(rig_ass)*sgn(left_ass)<0)
            return true;
        if (sgn(rig_ass)*sgn(left_ass)>0)
            return false;
        auto p=_upolynomial_Rtoab(G,a,b);
        uint64_t v=coeffsignchanges(p);
        if (v==0)
        {
            return false;
        }
        return true;    
    }
    int comp(uroot * r1,uroot* r2) // 1:r1<r2;0:r1=r2;-1:r1>r2
    {
        assert(r1->upolymap==r2->upolymap && r1->upolys==r2->upolys);
        int status=1;
        if (r1->left>r2->left)
        {
            std::swap(r1,r2);status=-1;
        }
        if (r1->right<r2->left)
            return status;
        if (r1->left!=r2->left || r1->right!=r2->right)
        {
            _uroot_check(r1,r2->left);
            _uroot_check(r2,r1->right);
        }
        if (r1->left!=r2->left || r1->right!=r2->right)
            return status;
        if (r1->poly_index==r2->poly_index)
            return 0;
        auto p=polynomial_GCD(r1->poly(),r2->poly());
        if (_uroot_check(p,r1->left,r1->right))
        {
            auto index=r1->add_poly(p);
            r1->poly_index=abs(index);
            r2->poly_index=abs(index);
            return 0;
        }
        else{
            auto index=r1->add_poly(r1->poly()/p);
            r1->poly_index=abs(index);
            index=r2->add_poly(r2->poly()/p);
            r2->poly_index=abs(index);
        }
        while (r1->left!=r2->left || r1->right!=r2->right)
        {
            _uroot_check(r1,(r1->left+r1->right)/2);
            _uroot_check(r2,(r2->left+r2->right)/2);
        }
        if (r1->left<=r2->left)
            return status;
        return -status;
    }
    

}
#endif