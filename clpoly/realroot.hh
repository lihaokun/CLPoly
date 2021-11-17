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
#include <stdexcept>
namespace clpoly{
    
    uint64_t coeffsignchanges(const upolynomial_<ZZ>& G);
    upolynomial_<ZZ> _upolynomial_1toinf(const upolynomial_<ZZ>& G);
    upolynomial_<ZZ> _upolynomial_Bto1(upolynomial_<ZZ> G,const ZZ &B);
    upolynomial_<ZZ> _upolynomial_01to012(const upolynomial_<ZZ>& G);
    upolynomial_<ZZ> _upolynomial_01to121(const upolynomial_<ZZ>& G);
    
    // void subcontraction_root_interval(upolynomial_<ZZ> G,QQ& B,QQ & E,const  QQ&  width=0)
    // {
    //     assert(E>B && width>=0);
    //     if (!width || (E-B<=width))
    //         return void();
    //     upolynomial_<ZZ> G_,G1;
    //     while (E-B>width)
    //     {
    //         QQ mid=(B+E)/2;
    //         if (!assign<QQ,ZZ,QQ>(G,mid))
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
    
    void subuspensky(const upolynomial_<ZZ>& G, std::vector<std::pair<QQ,QQ>>& l,const QQ & B=0,const QQ & E=1);
    ZZ RealRootBound(const upolynomial_<ZZ>& G);
    std::vector<std::pair<QQ,QQ>> uspensky(const upolynomial_<ZZ>& G);
    std::vector<upolynomial_<ZZ>> _upolynomial_1toinf(const std::vector<upolynomial_<ZZ>> &Gs);
    std::vector<upolynomial_<ZZ>>  _upolynomial_Bto1(std::vector<upolynomial_<ZZ>> Gs,const ZZ &B);
    void subuspensky
    (const std::vector<upolynomial_<ZZ>>& G,const std::vector<size_t>& I, std::vector<std::pair<QQ,QQ>>& l,std::vector<size_t>& index,const QQ & B=0,const QQ & E=1);
    std::pair<std::vector<std::pair<QQ,QQ>>,std::vector<size_t>> uspensky(std::vector<upolynomial_<ZZ>> G);//输入无平方基

    class uroot;
    template <class comp>
    std::pair<std::vector<uroot>,
    std::vector<std::vector<std::pair<uint64_t,uint64_t>>>>  //{poly,multiple}
     realroot(const std::vector<polynomial_<ZZ,comp>> & F)
    {
        if (F.empty())
            return {{},{}};
        
        if (get_variables(F).size()>1)
            throw std::invalid_argument("realroot:不是单变量的.");
        auto L=squarefreebasis(F);
        // std::cout<<L.first<<std::endl;
        std::vector<upolynomial_<ZZ>> G;
        G.reserve(L.first.size());
        for (auto &i:L.first)
        {
            G.push_back(i);
        }
        auto root=uspensky(G);
        // for (size_t i=0;i<root.first.size();++i)
        // {
        //     std::cout<<root.second[i]<<" ";
        // }
        // std::cout<<std::endl;
        std::vector<uroot> uroots;
        std::vector<std::vector<std::pair<uint64_t,uint64_t>>> I;
        uroots.reserve(root.first.size());
        I.reserve(root.first.size());
        for (size_t i=0;i<root.first.size();++i)
        {
            uroots.push_back(uroot(G[root.second[i]],std::move(root.first[i].first),std::move(root.first[i].second)));
            I.push_back(L.second[root.second[i]]);
        }
        return {std::move(uroots),std::move(I)};
    }
    // template <class comp>
    // std::vector<std::pair<QQ,QQ>> realroot(const upolynomial_<ZZ>& f)
    // {
    //     upolynomial_<ZZ> G=f/polynomial_GCD(f,derivative(f));
    //     return uspensky(G);
    // }

    class uroot;
    void _uroot_check(uroot* r,const QQ & mid);
    class uroot
    {
        private:
            int _is_inf=0;
            upolynomial_ZZ _poly;
            QQ _left;
            QQ _right;
        public:


            
            uroot(){}
            uroot(upolynomial_ZZ p,QQ _l,QQ _r)
            :_poly(std::move(p)),_left(std::move(_l)),_right(std::move(_r)),_is_inf(0)
            {}
            static uroot inf()
            {
                uroot u;
                u._is_inf=1;
                return u;
            }
            static uroot neginf()
            {
                uroot u;
                u._is_inf=-1;
                return u;
            }
            bool isinf()const 
            {
            return _is_inf==1;
            }
            bool isneginf()const
            {
            return _is_inf==-1;
            }

            const upolynomial_ZZ &  poly() const
            {
                return this->_poly;
            }
            upolynomial_ZZ &  poly() 
            {
                return this->_poly;
            }
            const QQ &  right() const
            {
                return this->_right;
            }
            QQ &  right() 
            {
                return this->_right;
            }
             const QQ &  left() const
            {
                return this->_left;
            }
            QQ &  left() 
            {
                return this->_left;
            }
            int static comp(uroot * r1,uroot* r2);
            int static comp(uroot * r,const QQ  & q);
            
            inline bool operator==(uroot & u)
            {
                return comp(this,&u)==0;
            }
            inline bool operator<(uroot & u)
            {
                return comp(this,&u)==1;
            }
            inline bool operator>(uroot & u)
            {
                return comp(this,&u)==-1;
            }
            inline bool operator<=(uroot & u)
            {
                return comp(this,&u)>=0;
            }
            inline bool operator>=(uroot & u)
            {
                return comp(this,&u)<=0;
            } 
            inline bool operator!=(uroot & u)
            {
                return comp(this,&u)!=0;
            }   
            void shrinkinterval()
            {
                if (left()!=right())
                    _uroot_check(this,(left()+right())/2);
            }
            bool is_single() const
            {
                return left()==right();
            }
            friend std::ostream& operator<<  (std::ostream& stream, const uroot& c) 
            {
                if (c.isinf())
                {
                stream<<"inf";
                return stream;
                }
                if (c.isneginf())
                {
                stream<<"-inf";
                return stream;
                }
                
                
                stream<<"{";
                stream<<c._left;
                stream<<",";
                stream<<c._right;
                stream<<"}";
                return stream;
            }         
    };
    inline bool operator!=(uroot & u,const QQ & q)
    {
        return uroot::comp(&u,q)!=0;
    }
    inline bool operator!=(const QQ & q,uroot & u)
    {
        return uroot::comp(&u,q)!=0;
    }  
    inline bool operator==(uroot & u,const QQ & q)
    {
        return uroot::comp(&u,q)==0;
    }
    inline bool operator==(const QQ & q,uroot & u)
    {
        return uroot::comp(&u,q)==0;
    }   
    inline bool operator>(uroot & u,const QQ & q)
    {
        return uroot::comp(&u,q)==-1;
    }
    inline bool operator<(const QQ & q,uroot & u)
    {
        return uroot::comp(&u,q)==-1;
    }
    inline bool operator<(uroot & u,const QQ & q)
    {
        return uroot::comp(&u,q)==1;
    }
    inline bool operator>(const QQ & q,uroot & u)
    {
        return uroot::comp(&u,q)==1;
    }  
    inline bool operator>=(uroot & u,const QQ & q)
    {
        return uroot::comp(&u,q)!=1;
    }
    inline bool operator<=(const QQ & q,uroot & u)
    {
        return uroot::comp(&u,q)!=1;
    } 
    inline bool operator<=(uroot & u,const QQ & q)
    {
        return uroot::comp(&u,q)!=-1;
    }
    inline bool operator>=(const QQ & q,uroot & u)
    {
        return uroot::comp(&u,q)!=-1;
    } 
    upolynomial_<ZZ> _upolynomial_Rtoab(const upolynomial_<ZZ>& G,const QQ &a,const QQ& b);

}
#endif