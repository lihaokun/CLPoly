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
    const upolynomial_<ZZ> __upolynomial_x_plus_1={{1,1},{0,1}};
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


    template <class comp>
    std::pair<std::vector<std::pair<QQ,QQ>>,
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
        std::vector<std::vector<std::pair<uint64_t,uint64_t>>> I;
        I.reserve(root.first.size());
        for (size_t i=0;i<root.first.size();++i)
        {
            I.push_back(L.second[root.second[i]]);
        }
        return {std::move(root.first),std::move(I)};
    }
    // template <class comp>
    // std::vector<std::pair<QQ,QQ>> realroot(const upolynomial_<ZZ>& f)
    // {
    //     upolynomial_<ZZ> G=f/polynomial_GCD(f,derivative(f));
    //     return uspensky(G);
    // }


   
    class uroot
    {
        public:
            // std::vector<upolynomial_ZZ>* upolys;
            // std::map<upolynomial_ZZ,size_t>* upolymap;
            upolynomial_ZZ poly;
            // size_t poly_index;
            QQ left;
            QQ right;
            int is_inf=0;
            uroot(){}
            uroot(upolynomial_ZZ p,QQ _l,QQ _r)
            :poly(p),left(_l),right(_r),is_inf(0)
            {}
            uroot static inf()
            {
                uroot u;
                u.is_inf=1;
                return u;
            }
            uroot static neginf()
            {
                uroot u;
                u.is_inf=-1;
                return u;
            }
            bool isinf()const 
            {
            return is_inf==1;
            }
            bool isneginf()const
            {
            return is_inf==-1;
            }
            
            int static comp(uroot * r1,uroot* r2);
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
                stream<<c.left;
                stream<<",";
                stream<<c.right;
                stream<<"}";
                return stream;
            }         
    };
    upolynomial_<ZZ> _upolynomial_Rtoab(const upolynomial_<ZZ>& G,const QQ &a,const QQ& b);

}
#endif