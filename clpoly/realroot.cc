/**
 * @file realroot.cc
 * @author 李昊坤 (ker@pm.me)
 * @brief 实根隔离
*/
#include <clpoly/realroot.hh>
namespace clpoly{
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
        // auto m=G.front().first.deg();
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
    void subuspensky(const upolynomial_<ZZ>& G, std::vector<std::pair<QQ,QQ>>& l,const QQ & B,const QQ & E)
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
        if (!assign<QQ,ZZ,QQ>(G,mid))
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
        // std::cout<<B<<std::endl;
        subuspensky(_upolynomial_Bto1(G,-B),l1,0,B);
        subuspensky(_upolynomial_Bto1(G,B),l2,0,B);
        l.reserve(l1.size()+l2.size());
        for (auto i=l1.rbegin();i!=l1.rend();++i)
            l.push_back({-i->second,-i->first});
        for (auto &i:l2)
            l.push_back(i);
        
        return l;
    }

    std::vector<upolynomial_<ZZ>> _upolynomial_1toinf(const std::vector<upolynomial_<ZZ>> &Gs)
    {
        std::vector<upolynomial_<ZZ>> Gout;
        Gout.reserve(Gs.size());
        for (auto & G:Gs)
        {
            if (G.empty())
            {
                Gout.push_back(G);
                continue;
            }
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
            Gout.push_back(g);
        }
        return Gout;
    }
    std::vector<upolynomial_<ZZ>>  _upolynomial_Bto1(std::vector<upolynomial_<ZZ>> Gs,const ZZ &B)
    {
        for (auto &G:Gs)
        {
            if (G.empty())
                continue;
            // auto m=G.front().first.deg();
            for (auto &i:G)
            {
                i.second*=pow(B,i.first.deg());
            }
        }
        return Gs;

    }
    void subuspensky(const std::vector<upolynomial_<ZZ>>& G,const std::vector<size_t>& I, std::vector<std::pair<QQ,QQ>>& l,std::vector<size_t>& index,const QQ & B,const QQ & E)
    {
        // std::cout<<B<<" "<<E<<std::endl;
        auto G_=_upolynomial_1toinf(G);
        
        std::vector<uint64_t> v;v.reserve(G_.size());
        uint64_t v0=0;
        uint64_t v1=0;
        size_t i_=0;
        for (size_t i=0;i<G_.size();++i)
        {
            v.push_back(coeffsignchanges(G_[i]));
            switch (v.back())
            {
            case 0:
                ++v0;
                break;
            case 1:
                i_=i;
                ++v1;
                break;

            }
        }
        // uint64_t v=coeffsignchanges(G_);
        // std::cout<<B<<" "<<E<<" "<<v<<std::endl;
        if (v0==G.size())   return void();
        if (v1==1 && v0==G.size()-1)
        {
            l.push_back({B,E});
            index.push_back(I[i_]);
            return void();
        }
        
        QQ mid=(B+E)/2;
        std::vector<upolynomial_<ZZ>> G_1;
        std::vector<size_t> I_1;
        G_1.reserve(G.size());
        I_1.reserve(G.size());
        for (size_t i=0;i<G_.size();++i)
        {
            if (!assign<QQ,ZZ,QQ>(G[i],mid))
            {
                l.push_back({mid,mid});
                index.push_back(I[i]);
                if (v[i]==1)
                {
                    v[i]=0;
                }
                break;
            }
        }
        for (size_t i=0;i<G.size();++i)
        {
            if (v[i]>0)
            {
                G_1.push_back(_upolynomial_01to012(G[i]));
                I_1.push_back(I[i]);
            }
        }
        subuspensky(G_1,I_1,l,index,B,mid);
        G_1.clear();
        for (size_t i=0;i<G.size();++i)
        {
            if (v[i]>0)
            {
                G_1.push_back(_upolynomial_01to121(G[i]));
            }
        }
        subuspensky(G_1,I_1,l,index,mid,E);
    } 

    std::pair<std::vector<std::pair<QQ,QQ>>,std::vector<size_t>> uspensky(std::vector<upolynomial_<ZZ>> G)//输入无平方基
    {
        std::vector<std::pair<QQ,QQ>> l,l1,l2;
        std::vector<size_t> I;
        I.reserve(G.size());
        std::vector<size_t> index1,index2;
        {
            size_t tmp=0;
            for (size_t i=0;i<G.size();++i)
            {
                if (!is_number(G[i]))
                {
                    I.push_back(i);
                    if (i!=tmp)
                        G[tmp]=std::move(G[i]);
                    ++tmp;

                }
            }
            if (tmp==0)
                return {{},{}};
            G.resize(tmp);
        }

        ZZ B=0;
        for (auto &i:G)
        {
            auto tmp=RealRootBound(i);
            if (tmp>B)
                B=std::move(tmp);
        }
        // std::cout<<B<<std::endl;
        
        subuspensky(_upolynomial_Bto1(G,-B),I,l1,index1,0,B);
        subuspensky(_upolynomial_Bto1(G,B),I,l2,index2,0,B);
        l.reserve(l1.size()+l2.size());
        std::vector<size_t> index(index1.rbegin(),index1.rend());
        
        for (auto i=l1.rbegin();i!=l1.rend();++i)
        {
            l.push_back({-i->second,-i->first});
        }
        for (size_t i=0;i<G.size();++i)
        {
            if (!assign<QQ,ZZ,QQ>(G[i],0))
            {
                l.push_back({0,0});
                index.push_back(I[i]);
                break;
            }
        }
        index.insert(index.end(),index2.begin(),index2.end());
        l.insert(l.end(),l2.begin(),l2.end());
        
        return {l,index};
    }
    
    
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
        QQ rig_ass=assign<QQ,ZZ,QQ>(r->poly,r->right);
        QQ mid_ass=assign<QQ,ZZ,QQ>(r->poly,mid);
        QQ left_ass=assign<QQ,ZZ,QQ>(r->poly,r->left);
        if (mid_ass==0)
        {
            r->left=r->right=mid;
            return void();
        }
        if (rig_ass==0)
        {
            if (left_ass==0) // 不应该进入
            {
                auto p=_upolynomial_Rtoab(r->poly,r->left,mid);
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
        QQ rig_ass=assign<QQ,ZZ,QQ>(G,b);
        if (a==b && rig_ass==0)
            return true;
        QQ left_ass=assign<QQ,ZZ,QQ>(G,a);
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
    int  uroot::comp(uroot * r1,uroot* r2) // 1:r1<r2;0:r1=r2;-1:r1>r2
    {
        if (r1->is_inf=-1)
            return int(r2->is_inf!=-1);
        if (r1->is_inf=1)
            return -int(r2->is_inf!=1);
        if (r2->is_inf)
            return r2->is_inf;    
        // if (r1->upolymap!=r2->upolymap && r1->upolys!=r2->upolys)
        //     return 3;
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
        if (r1->poly==r2->poly)
            return 0;
        auto p=polynomial_GCD(r1->poly,r2->poly);
        if (_uroot_check(p,r1->left,r1->right))
        {
            r1->poly=p;
            r2->poly=p;
            return 0;
        }
        else{
            r1->poly=r1->poly/p;
            r2->poly=r2->poly/p;
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