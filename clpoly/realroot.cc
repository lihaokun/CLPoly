/**
 * @file realroot.cc
 * @author 李昊坤 (ker@pm.me)
 * @brief 实根隔离
*/
#include <clpoly/realroot.hh>
namespace clpoly{
    const upolynomial_<ZZ> __upolynomial_x_plus_1={{1,1},{0,1}};
    const QQ __QQ_1_2=QQ(1,2);
    
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
    bool uspensky_shrink(const upolynomial_<ZZ>& G,QQ & B,QQ & E,bool lb=false,bool rb=false);
    void subuspensky(const upolynomial_<ZZ>& G, std::vector<std::pair<QQ,QQ>>& l,const QQ & B,const QQ & E)
    {
        // std::cout<<"G:"<<G<<std::endl;
        upolynomial_<ZZ> G_=_upolynomial_1toinf(G);
        uint64_t v=coeffsignchanges(G_);
        // std::cout<<B<<" "<<E<<" "<<v<<std::endl;
        if (v==0)   return void();
        if (v==1 )
        {
            QQ B1=B;
            QQ E1=E;
            uspensky_shrink(G,B1,E1);
            l.push_back({std::move(B1),std::move(E1)});
            return void();
        }
        QQ mid=(B+E)/2;
        subuspensky(_upolynomial_01to012(G),l,B,mid);
        if (!assign<QQ,ZZ,QQ>(G,__QQ_1_2))
        {
            l.push_back({mid,mid});
        }
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
    bool uspensky_shrink(const upolynomial_<ZZ>& G,QQ & B,QQ & E,bool lb,bool rb)
    {
        auto G_=_upolynomial_1toinf(G);
        auto v=coeffsignchanges(G_);
        if (v==0)
            return false;
        if (E-B<1 && (lb && rb))
            return true;
        QQ m=(E+B)/2;
        if (!assign<QQ,ZZ,QQ>(G,__QQ_1_2))
        {
            B=m;E=m;return true;
        }
        if (uspensky_shrink(_upolynomial_01to012(G),B,m,lb,true))
        {
            E=m;
            return true;   
        }
        
        uspensky_shrink(_upolynomial_01to121(G),m,E,true,rb);
        B=m;
        return true;

    }
    upolynomial_<ZZ> _upolynomial_Rtoab(const upolynomial_<ZZ>& G,const QQ &a,const QQ& b)
    {
        if (G.empty() || is_number(G))
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
    
    void uspensky_shrink_(const upolynomial_<ZZ>&p,  QQ & B, QQ & E)
    {
        int BS=sgn(assign<QQ,ZZ,QQ>(p,B));
        int ES=sgn(assign<QQ,ZZ,QQ>(p,E));
        bool l,r;
        QQ ml=(B+E)/2;
        QQ mr=ml;
        int msl=sgn(assign<QQ,ZZ,QQ>(p,ml));
        int msr=msl;
        if (msl==0)
        {
            B=ml;E=ml;
            return;
        }
        // std::cout<<"uspensky_shrink_ "<<B<<" "<<E<<std::endl;
        while (1)
        {
            // std::cout<<"ml,mr:"<<ml<<" "<<mr<<std::endl;
        
            if (BS!=msl)
            {
                QQ ml1=(ml+B)/2;
                // std::cout<<"ml1:"<<ml1<<std::endl;
                int msl1=sgn(assign<QQ,ZZ,QQ>(p,ml1));
                if (msl1==0)
                {
                    B=ml1;E=ml1;
                    // std::cout<<"new B E "<<B<<" "<<E<<std::endl;
                    
       
                    return;
                }
                if (msl1!=msl)
                {
                    B=ml1;E=ml;
                    msr=msl;msl=msl1;
                    break;
                }
                ml=ml1;msl=msl1;
            }
            if (ES!=msr)
            {
                QQ mr1=(mr+E)/2;
                int msr1=sgn(assign<QQ,ZZ,QQ>(p,mr1));
                if (msr1==0)
                {
                    B=mr1;E=mr1;
                    // std::cout<<"new B E "<<B<<" "<<E<<std::endl;
       
                    return;
                }
                if (msr1!=msr)
                {
                    B=mr;E=mr1;
                    msl=msr;msr=msr1;
                    break;
                }
                mr=mr1;msr=msr1;
            }
        }
        assert(B<=E);
        // std::cout<<"new B E "<<B<<" "<<E<<std::endl;
       
        while (E-B>1)
        {
            QQ m=(B+E)/2;
            int ms=sgn(assign<QQ,ZZ,QQ>(p,m));
            if (ms!=msr)
            {
                B=m;
                msl=ms;
            }
            else{
                E=m;
                msr=ms;    
            }
        }
    }

    void subuspensky(const std::vector<upolynomial_<ZZ>>& G,const QQ & B,const QQ & E,std::vector<std::pair<QQ,QQ>>& l,std::vector<size_t>& index)
    {
        std::list<std::pair<std::pair<QQ,QQ>,std::vector<size_t>>> queue;
        {
            std::vector<size_t> tmp;
            for (size_t i=0;i!=G.size();++i)
            {
                tmp.push_back(i);
            }
            queue.push_back({{B,E},tmp});
        }
        // std::vector<uint64_t> v;v.resize(G.size());
        while (!queue.empty())
        {
            // std::cout<<"new"<<std::endl;
            
            QQ a=std::move(queue.front().first.first);
            QQ b=std::move(queue.front().first.second);
            std::vector<size_t> p_s=std::move(queue.front().second);
            queue.pop_front();
            
            if (a==b)
            {
                l.push_back({std::move(a),std::move(b)});
                index.push_back(p_s[0]);
                continue;
            }

            // std::cout<<a<<" "<<b<<" "<<p_s.size()<<std::endl;
            uint64_t v0=0;
            uint64_t v1=0;
            size_t i_=0;
            QQ mid=(a+b)/2;
            std::vector<size_t> p_s_1;p_s_1.reserve(p_s.size());
            size_t is_mid=0;
            for (auto  i:p_s)
            {
                uint64_t v=coeffsignchanges(_upolynomial_Rtoab(G[i],a,b));
                if (!v)
                {
                    ++v0;
                    continue;
                }
                if (!assign<QQ,ZZ,QQ>(G[i],mid))
                {
                    is_mid=i+1;
                    v-=1;
                }
                switch (v)
                {
                    case 0:
                        ++v0;
                        break;
                    case 1:
                        i_=i;
                        ++v1;
                        break;
                }
                if (v)
                    p_s_1.push_back(i);
            }
            // std::cout<<"v:"<<v0<<" "<<v1<<std::endl;
            
            if (is_mid==0 && v1==1 && v0==p_s.size()-1)
            {
                uspensky_shrink_(G[i_],a,b);
                l.push_back({a,b});
                index.push_back(i_);
                continue;
            }
            if (!p_s_1.empty())
            {
                queue.push_front({{mid,b},p_s_1});
            }
            if (is_mid!=0)
            {
                queue.push_front({{mid,mid},{is_mid-1}});    
            }
            if (!p_s_1.empty())
            {
                queue.push_front({{a,mid},std::move(p_s_1)});
            }
            
        }
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
            auto B1=B;
            auto E1=E;
            uspensky_shrink(G[i_],B1,E1);
            l.push_back({B1,E1});
            index.push_back(I[i_]);
            return void();
        }
        
        QQ mid=(B+E)/2;
        std::vector<upolynomial_<ZZ>> G_1;
        std::vector<size_t> I_1;
        G_1.reserve(G.size());
        I_1.reserve(G.size());
        // std::cout<<mid<<" " <<__QQ_1_2<<std::endl;
        std::vector<std::pair<QQ,QQ>> _l;
        std::vector<size_t> _index;
        for (size_t i=0;i<G_.size();++i)
        {
            // std::cout<<G[i]<<":"<<assign<QQ,ZZ,QQ>(G[i],__QQ_1_2)<<std::endl;
            if (!assign<QQ,ZZ,QQ>(G[i],__QQ_1_2))
            {
                _l.push_back({mid,mid});
                _index.push_back(I[i]);
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
        l.insert(l.end(),_l.begin(),_l.end());
        index.insert(index.end(),_index.begin(),_index.end());
        
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
        
        // {
        //     std::vector<std::pair<QQ,QQ>> l;
        //     std::vector<size_t> index;
        //     ZZ B=0;
        //     for (auto &i:G)
        //     {
        //         auto tmp=RealRootBound(i);
        //         if (tmp>B)
        //             B=std::move(tmp);
        //     }
        //     auto t=clock();
        //     subuspensky(G,-B,B,l,index);
        //     std::cout<<"new subuspensky time="<<double(clock()-t)/CLOCKS_PER_SEC<<"s\n";
        //     for (auto& i:l)
        //     {
        //         std::cout<<"{"<<i.first<<" , "<<i.second<<"}, ";
        //     }
        //     std::cout<<std::endl;
        // }
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
        // auto t=clock();
        subuspensky(_upolynomial_Bto1(G,-B),I,l1,index1,0,B);
        subuspensky(_upolynomial_Bto1(G,B),I,l2,index2,0,B);
        // std::cout<<"subuspensky time="<<double(clock()-t)/CLOCKS_PER_SEC<<"s\n";
    
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
    
    

    void _uroot_check(uroot* r,const QQ & mid)
    {
        if (mid<r->left() || mid >r->right() ||r->left()==r->right() )
            return void();
        QQ rig_ass=assign<QQ,ZZ,QQ>(r->poly(),r->right());
        QQ mid_ass=assign<QQ,ZZ,QQ>(r->poly(),mid);
        QQ left_ass=assign<QQ,ZZ,QQ>(r->poly(),r->left());
        if (mid_ass==0)
        {
            r->left()=r->right()=mid;
            return void();
        }
        if (rig_ass==0)
        {
            if (left_ass==0) // 不应该进入
            {
                auto p=_upolynomial_Rtoab(r->poly(),r->left(),mid);
                uint64_t v=coeffsignchanges(p);
                if (v==0)
                {
                    r->left()=mid;return void();
                }
                r->right()=mid;return void();
            }
            else
            {
                if (sgn(mid_ass)==sgn(left_ass))
                {
                    r->left()=mid;return void();
                }
                r->right()=mid;return void();
            }
        }
        if (sgn(mid_ass)==sgn(rig_ass))
        {
            r->right()=mid;return void();  
        }
        r->left()=mid;return void();
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
        if (r1->isneginf())
            return int(!r2->isneginf());
        if (r1->isinf())
            return -int(!r2->isinf());
        if (r2->isinf())
            return 1;
        if (r2->isneginf())
            return -1;  
        if (r1->left()==r1->right())
        {
            return -uroot::comp(r2,r1->left());
        }
        if (r2->left()==r2->right())
        {
            return uroot::comp(r1,r2->left());
        }
        // if (r1->upolymap!=r2->upolymap && r1->upolys!=r2->upolys)
        //     return 3;
        //  std::cout<<*r1<<*r2<<std::endl;
        int status=1;
        if (r1->left()>r2->left())
        {
            std::swap(r1,r2);status*=-1;
        }
        if (r1->right()<r2->left())
            return status;
        //  std::cout<<*r1<<*r2<<std::endl;
        if (r1->left()!=r2->left() || r1->right()!=r2->right())
        {
            _uroot_check(r1,r2->left());
            _uroot_check(r1,r2->right());
            _uroot_check(r2,r1->right());
        }
        //  std::cout<<*r1<<*r2<<std::endl;
        if (r1->left()>r2->left())
        {
            std::swap(r1,r2);status*=-1;
        }
        if (r1->left()!=r2->left() || r1->right()!=r2->right())
            return status;
        if (r1->poly()==r2->poly())
            return 0;
        auto p=polynomial_GCD(r1->poly(),r2->poly());
        if (_uroot_check(p,r1->left(),r1->right()))
        {
            r1->poly()=p;
            r2->poly()=p;
            return 0;
        }
        else{
            r1->poly()=r1->poly()/p;
            r2->poly()=r2->poly()/p;
        }
        // std::cout<<r1->poly()<<" "<<r2->poly()<<std::endl;
        while (r1->left()==r2->left()  && r1->right()==r2->right())
        {
            // std::cout<<*r1<<*r2<<std::endl;
            _uroot_check(r1,(r1->left()+r1->right())/2);
            _uroot_check(r2,(r2->left()+r2->right())/2);
        }
        if (r1->left()<=r2->left())
            return status;
        return -status;
    }

    int uroot::comp(uroot * r,const QQ  & q) // 1:r<q;0:r=q;-1:r>q
    {
        if (r->isneginf())
            return 1;
        if (r->isinf())
            return -1;
        
        if (q==r->left() && r->left()==r->right() )
            return 0;
        
        if (q<=r->left())
            return -1;
        if (q>=r->right())
            return 1;
        
        if (assign<QQ,ZZ,QQ>(r->poly(),q)==0)
        {
            r->left()=r->right()=q;
            return 0;
        }
        _uroot_check(r,q);
        if (q<=r->left())
            return -1;
        return 1;
    }
}