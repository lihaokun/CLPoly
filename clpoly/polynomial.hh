/*
Module Name:
    polynomial.hh
Abstract:
    定义polynomial和相关的运算
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_POLYNOMIAL_HH
#define CLPOLY_POLYNOMIAL_HH

#include <clpoly/monomial.hh>
#include <clpoly/number.hh>
#include <clpoly/basic_polynomial.hh>
#include <list>
#include <string>
#include <random>

namespace clpoly{
    template <class Tc,class comp=grlex>
    using polynomial_=basic_polynomial<basic_monomial<comp>,Tc,comp>;
    using polynomial_ZZ=polynomial_<ZZ>;
    using polynomial_QQ=polynomial_<QQ>;
    
/*======================================主要函数=========================================*/
    // template <class Tc,class comp>
    // int64_t degree(const polynomial_<Tc,comp> & p);
    // template <class Tc,class comp>
    // int64_t degree(const polynomial_<Tc,comp> & p,const variable & v);
    // template<class Tc>
    // polynomial_<Tc> random_polynomial(const std::vector<variable> & v,uint64_t deg,double p,int up,int down);
    // template<class T1,class T2,class comp1,class comp2>
    // void poly_convert(const polynomial_<T1,comp1>& p_in,polynomial_<T2,comp2> & p_out);
    // template <class Tc,class comp>
    // inline polynomial_<Tc,comp> leadcoeff(const polynomial_<Tc,comp>&F);
    // template <class Tc,class comp>
    // inline polynomial_<Tc,comp> prem(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v);
    // template <class Tc,class comp>
    // inline polynomial_<Tc,comp> resultant(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v);
    // template <class Tc, class To>
    // polynomial_<Tc,To> association(const polynomial_<Tc,To>& p,const variable & v,const Tc & c);
    // template <class compare>
    // polynomial_<Zp,compare> polynomial_mod(polynomial_<ZZ,compare> && p, uint32_t prime)
/*======================================实现===============================================*/
    inline monomial pow(const variable & v,int64_t i)
    {
        return monomial({{v,i}});
    }
    template <class compare>
    inline  basic_monomial<compare> pow(const basic_monomial<compare> & v,int64_t i)
    {
        basic_monomial<compare> m1(v.comp_ptr());
        if (i)
        {
            m1=v;
            for (auto & i:m1)
            {
                i.second*=i;
            }
        }
        return m1;
    }


    template <class Tc,class compare>
    inline polynomial_<Tc,compare> pow(const polynomial_<Tc,compare> & p,int64_t i)
    {
        assert(i>=0);
        polynomial_<Tc,compare> o(p.comp_ptr());
        o={{{},1}};
        switch (i)
        {
        case 0:
            return o;
            break;
        case 1:
            return p;
            break;
        case 2:
            return p*p;
            break;
        case 3:
            return p*p*p;
            break;         
        default:
            pair_vec_power(o.data(),p.data(),i,p.comp());
            return o;
            break;
        }
    }
    
    monomial operator* (const variable & v1,const variable & v2)
    {
        if (monomial::compare_type()(v1,v2))  
            return monomial({{v1,1},{v2,1}});
        return monomial({{v2,1},{v1,1}});          
    }
    polynomial_ZZ operator -(const monomial & m)
    {
        return polynomial_ZZ({{m,-1}});
    }

    polynomial_ZZ operator+(const monomial & m,int64_t i)
    {
        return polynomial_ZZ({{m,1},{{},i}});
    }
    // polynomial_ZZ operator+(const monomial & m,const variable & v)
    // {
    //     return polynomial_ZZ({{m,1},{{{v,1}},1}});
    // }
    polynomial_ZZ operator-(const monomial & m,int64_t i)
    {
        return polynomial_ZZ({{m,1},{{},-i}});
    }

    monomial operator* (const monomial & v1,const variable & v2)
    {
        return v1*monomial(v2);          
    }
    monomial operator* (const variable & v2,const monomial & v1)
    {
        return v1*monomial(v2);          
    }

    polynomial_ZZ operator* (const ZZ & v1,const monomial & v2)
    {
        return polynomial_ZZ({{v2,v1}});          
    }

    polynomial_ZZ operator* (const ZZ & v1,const variable & v2)
    {
        return polynomial_ZZ({{monomial(v2),v1}});          
    }

    polynomial_ZZ operator* (int64_t v1,const variable & v2)
    {
        return polynomial_ZZ({{monomial(v2),ZZ(v1)}});          
    }
    
    polynomial_ZZ operator* (const polynomial_ZZ &  v1,const variable & v2)
    {
        return v1*monomial(v2); 
    }

    // template <class Tc>
    // polynomial_<Tc> operator* (const polynomial_<Tc> &  v1,const monomial & v2)
    // {
    //     polynomial_<Tc> O;
    //     O.reserve(v1.size());
    //     for (auto &i:v1)
    //         O.push_back({i.first*v2,i.second});
    //     return O;
    // }

    polynomial_ZZ operator+ (const monomial & m1,const monomial & m2)
    {
        if (polynomial_ZZ::compare_type()(m1,m2))
            return polynomial_ZZ({{m1,1},{m2,1}});
        return polynomial_ZZ({{m2,1},{m1,1}});
    }
    
    template<class Tc>
    polynomial_<Tc> operator+ (const monomial & m,const polynomial_<Tc> & p)
    {
        return p+polynomial_<Tc>({{m,1}});
    }

    template<class Tc>
    polynomial_<Tc> operator+ (const polynomial_<Tc> & p,const monomial & m)
    {
        return p+polynomial_<Tc>({{m,1}});
    }
    template<class Tc>
    polynomial_<Tc> operator+ (const polynomial_<Tc> & p,const variable & m)
    {
        return p+polynomial_<Tc>({{{{m,1}},1}});
    }
    polynomial_ZZ operator- (const monomial & m1,const monomial & m2)
    {
        if (polynomial_ZZ::compare_type()(m1,m2))
            return polynomial_ZZ({{m1,1},{m2,-1}});
        return polynomial_ZZ({{m2,-1},{m1,1}});
    }
    template<class Tc>
    polynomial_<Tc> operator- (const polynomial_<Tc> & p,const monomial & m)
    {
        return p-polynomial_<Tc>({{m,1}});
    }
    // template<class Tc>
    // polynomial_<Tc> operator+ (const monomial & m,const Tc & p)
    // {
    //     return polynomial_<Tc>({{m,1},{{},p}});
    // }

    template<class Tc>
    polynomial_<Tc> operator+ (const Tc & m,const polynomial_<Tc> & p)
    {
        return p+m;
        
    }
    template<class Tc>
    polynomial_<Tc> operator+ (const polynomial_<Tc> & p,const Tc & m)
    {
        polynomial_<Tc> O;
        O=p;
        if (p.back().first.empty())
            O.back().second+=m;
        else
            O.push_back({{},m});
        return O;
    }
    template<class Tc>
    polynomial_<Tc> operator+ (const polynomial_<Tc> & p,int64_t m)
    {
        return p+Tc(m);
    }
    polynomial_ZZ operator+ (const variable & p,int64_t m)
    {
        return {{{{p,1}},1},{{},ZZ(m)}};
    }
    template<class Tc>
    polynomial_<Tc> operator- (const polynomial_<Tc> & p,const Tc & m)
    {
        polynomial_<Tc> O;
        O=p;
        if (p.back().first.empty())
            O.back().second-=m;
        else
            O.push_back({{},-m});
        return O;
    }
    template<class Tc>
    polynomial_<Tc> operator- (const polynomial_<Tc> & p,int64_t m)
    {
        return p-Tc(m);
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
                p1.push_back({{},Tc(dis(gen))});
            }
            
        }
        return polynomial_<Tc>(p1);
    }
    template<class Tc,class comp>
    std::list<std::pair<variable,int64_t>> get_variables(const polynomial_<Tc,comp>& p)
    {
        std::list<std::pair<variable,int64_t>> l;
        __pair_vec_variables(p.data(),l);
        return l;
    }

    template <class Tc>
    int64_t get_deg(const polynomial_<Tc> & p)
    {
        if (p.empty())
            return 0;
        return p.front().first.deg(); 
    }

    template <class Tc,class comp>
    int64_t get_deg(const polynomial_<Tc,comp> & p)
    {
        if (p.empty())
            return 0;
        int64_t deg=p.front().first.deg();
        for(auto &i:p)
            deg=std::max(i.first.deg(),deg);
        return deg;
    }


    template <class Tc,class comp>
    int64_t degree(const polynomial_<Tc,comp> & p)
    {
        return p.degree();
    }
    template <class Tc,class comp>
    int64_t degree(const polynomial_<Tc,comp> & p,const variable & v)
    {
        // auto & l=p.variables();
        // for (auto &i:l)
        //     if (l.first==v)
        //         return l.second;
        // return 0;
        int64_t deg=0,tmp;
        for (auto &i:p)
            if ((tmp=i.first.deg(v))>deg)
                deg=tmp;
        return deg;
    }
    
    template <class Tc>
    constexpr const int64_t get_up_deg(const polynomial_<Tc,univariate_priority_order>& p)
    {
        return p.empty()?-1:get_up_deg(p.begin()->first);
    }
    
    template <class Tc>
    int64_t degree(const polynomial_<Tc,univariate_priority_order> & p,const variable & v)
    {
        if (p.comp().v==v)
            return get_up_deg(p);
        int64_t deg=0,tmp;
        for (auto &i:p)
            if ((tmp=i.first.deg(v))>deg)
                deg=tmp;
        return deg;
    }

    template<class T1,class T2,class comp1,class comp2>
    void poly_convert(const polynomial_<T1,comp1>& p_in,polynomial_<T2,comp2> & p_out)
    {
        assert((void*)&p_in!=(void*)&p_out);
        p_out.clear();
        basic_monomial<comp2> m(p_out.comp_ptr());
        for (auto &i:p_in)
        {
            m=i.first.data();
            p_out.push_back({std::move(m),i.second});
        }
        p_out.normalization();
    }
    template<class T1,class T2,class comp1,class comp2>
    void poly_convert(polynomial_<T1,comp1>&& p_in,polynomial_<T2,comp2> & p_out)
    {
        assert((void*)&p_in!=(void*)&p_out);
        p_out.clear();
        basic_monomial<comp2> m(p_out.comp_ptr());
        for (auto &i:p_in)
        {
            m=std::move(i.first.data());
            p_out.push_back({std::move(m),std::move(i.second)});
        }
        p_out.normalization();
        p_in.clear();
    }
    template<class T2,class comp1,class comp2>
    void poly_convert(const polynomial_<Zp,comp1>& p_in,polynomial_<T2,comp2> & p_out)
    {
        assert((void*)&p_in!=(void*)&p_out);
        p_out.clear();
        basic_monomial<comp2> m(p_out.comp_ptr());
        for (auto &i:p_in)
        {
            m=i.first.data();
            p_out.push_back({std::move(m),uint64_t(i.second)});
        }
        p_out.normalization();
    }
    template<class T2,class comp1,class comp2>
    void poly_convert(polynomial_<Zp,comp1>&& p_in,polynomial_<T2,comp2> & p_out)
    {
        assert(&p_in!=&p_out);
        p_out.clear();
        basic_monomial<comp2> m(p_out.comp_ptr());
        for (auto &i:p_in)
        {
            m=std::move(i.first.data());
            p_out.push_back({std::move(m),uint64_t(i.second)});
        }
        p_out.normalization();
        p_in.clear();
    }



    template <class Tc>
    inline int64_t leadcoeff(polynomial_<Tc,univariate_priority_order>&O,const polynomial_<Tc,univariate_priority_order>&F)
    {
        assert(&O!=&F);
        int64_t d=get_up_deg(F);
        if (d && !F.empty())
        {
            O.clear();
            O.comp(F.comp_ptr());
            for (auto & i:F)
            {
                if (get_up_deg(i.first)!=d)
                    return d;
                O.push_back(i);
                O.back().first.pop_back();
            }
        }
        else
            O=F;
        return d;
    }
    template <class Tc,class comp >
    inline int64_t leadcoeff(polynomial_<Tc,comp>&O,const polynomial_<Tc,comp>&F,const variable & v)
    {
        assert(&O!=&F);
        int64_t deg=0,tmp;
        for (auto &i:F)
            if ((tmp=i.first.deg(v))>deg)
                deg=tmp;
        if (deg)
        {
            O.clear();
            O.comp(F.comp_ptr());
            basic_monomial<comp> m(F.comp_ptr());
            typename basic_monomial<comp>::const_iterator tmp_I;
            for (auto &i:F)
                if ((tmp_I=i.first.find(v))!=i.first.end() && tmp_I->second==deg)
                {
                    m=i.first;
                    m[tmp_I-i.first.begin()].second=0;
                    O.push_back({std::move(m),i.second});
                }
            O.normalization();
        }
        else 
            O=F;
        return deg;
    }
    
    template <class Tc>
    inline polynomial_<Tc,univariate_priority_order> leadcoeff(const polynomial_<Tc,univariate_priority_order>&F)
    {
        polynomial_<Tc,univariate_priority_order> O(F.comp_ptr());
        leadcoeff(O,F);
        return O;
    }
    
    template <class Tc,class comp>
    inline polynomial_<Tc,comp> leadcoeff(const polynomial_<Tc,comp>&F)
    {
        polynomial_<Tc,comp> O(F.comp_ptr());
        leadcoeff(O,F);
        return O;
    }

    #define prem_v1 prem
    
    template <class Tc>
    inline polynomial_<Tc,univariate_priority_order> prem(const polynomial_<Tc,univariate_priority_order> &G,const polynomial_<Tc,univariate_priority_order> & F)
    {
        polynomial_<Tc,univariate_priority_order>  O1(&G.comp_ptr());
        prem(O1,G,F);
        return O1;
    }
    template <class Tc,class comp>
    inline polynomial_<Tc,comp> prem(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v)
    {
        if (F.size()==0)
        {
//#ifndef NDEBUG  
            throw std::invalid_argument("Error:polynomial prem div 0.");
//#endif            
            return polynomial_<Tc,comp>();
        }
        polynomial_<Tc,comp> O(G.comp_ptr());
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order>  G1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        polynomial_<Tc,univariate_priority_order>  O1(&comp_v);
        prem(O1,G1,F1);
        poly_convert(std::move(O1),O);
        return O;
    }

    
    template <class Tc>
    inline void prem_v2(polynomial_<Tc,univariate_priority_order> &O,const polynomial_<Tc,univariate_priority_order> &G,const polynomial_<Tc,univariate_priority_order> & F)
    {
        if (F.size()==0)
        {
//#ifndef NDEBUG  
            throw std::invalid_argument("Error:polynomial prem div 0.");
//#endif            
            return void();
        }
        assert(comp_consistent(G.comp(),F.comp()));
        assert((void*)&O!=(void*)&G && (void*)&O!=(void*)&F);
        // std::cout<<"G="<<G<<std::endl;
        // std::cout<<"F="<<F<<std::endl;
        polynomial_<Tc,univariate_priority_order>  Fl(F.comp_ptr());
        int64_t F_deg=leadcoeff(Fl,F);
        int64_t G_deg=get_up_deg(G);
        if (F_deg>G_deg) 
        {
            O=G;
            return void();
        }
        // std::cout<<"G_deg="<<G_deg<<" F_deg="<<F_deg<<std::endl;
        // std::cout<<"Fl="<<Fl<<std::endl;
        O=polynomial_rem(G*pow(Fl,G_deg-F_deg+1),F);
    }


    template <class Tc,class compare>
    void __onestep__prem_v1(
        const std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & G,
        const std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & F0,
        int64_t f_deg,
        const std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & F,
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & O,
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & F1,
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & G1,
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> & G2,
        compare comp 
    )
    {

        auto G_ptr=G.begin();
        auto G_end=G.end();
        int64_t g_deg=get_up_deg(G_ptr->first);
        G1.clear();
        G1.reserve(G.size());
        if (g_deg-f_deg)
        {
            while(G_ptr!=G_end && get_up_deg(G_ptr->first)==g_deg)
            {
                G1.push_back(*(G_ptr++));
                G1.back().first.back().second=g_deg-f_deg;
                G1.back().first.deg()-=f_deg;
            }
        }
        else
        {
            while(G_ptr!=G_end && get_up_deg(G_ptr->first)==g_deg)
            {
                G1.push_back(*(G_ptr++));
                G1.back().first.pop_back();
            }
        }

        pair_vec_multiplies(F1,G1,F,comp);
        if(G_ptr==G_end)
        {
            pair_vec_negate(F1);
            std::swap(F1,O);
        }
        else{
            G2.clear();
            G2.reserve((G_end-G_ptr));
            while(G_ptr!=G_end)
                G2.push_back(*(G_ptr++));
            pair_vec_multiplies(G1,G2,F0,comp);
            pair_vec_sub(O,G1,F1,comp);
        }
    }
    
    template <class Tc>
    void prem_v1(
            polynomial_<Tc,univariate_priority_order>&O,
            //polynomial_<Tc,univariate_priority_order>&L,
            const polynomial_<Tc,univariate_priority_order>&G,
            const polynomial_<Tc,univariate_priority_order>&F,
            bool is_L=true
            )
    {
        if (F.size()==0)
        {
//#ifndef NDEBUG  
            throw std::invalid_argument("Error:polynomial prem div 0.");
//#endif            
            return void();
        }
        O.clear();
        int64_t f_deg=get_up_deg(F);
        int64_t g_deg=get_up_deg(G);
        if (f_deg<=0)
            return void();
        if (g_deg<f_deg)
        {
            O=G;
            return void();
        }
        int64_t d=g_deg-f_deg;
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_G;
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_G1;
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_G2;
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_F0;
        tmp_F0.reserve(F.size());
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_F_;
        tmp_F_.reserve(F.size());
        for (auto & i:F)
            if (get_up_deg(i.first)!=f_deg)
                tmp_F_.push_back(i);
            else
            {
                tmp_F0.push_back(i);
                tmp_F0.back().first.pop_back();
            }
        std::vector<std::pair<basic_monomial<univariate_priority_order>,Tc>> tmp_F1;
        O.data().reserve(G.size());
        __onestep__prem_v1(G.data(),tmp_F0,f_deg,tmp_F_,O.data(),tmp_F1,tmp_G1,tmp_G2,G.comp());
        // std::cout<<"prem_v1_:"<<O<<std::endl;
        while (get_up_deg(O)>=f_deg)
        {
            std::swap(tmp_G,O.data());
            __onestep__prem_v1(tmp_G,tmp_F0,f_deg,tmp_F_,O.data(),tmp_F1,tmp_G1,tmp_G2,G.comp());    
            --d;
            // std::cout<<"prem_v1_:"<<O<<std::endl;
        }
        if (is_L && d>0)
        {
            if (d==1)
            {
                pair_vec_multiplies(tmp_G,O.data(),tmp_F0,G.comp());
                std::swap(tmp_G,O.data());
            }
            else
            {
                pair_vec_power(tmp_G1,tmp_F0,d,G.comp());
                pair_vec_multiplies(tmp_G,O.data(),tmp_G1,G.comp());
                std::swap(tmp_G,O.data());
            }
            
        }
        
        //L.data()=std::move(tmp_F0);
        // L=pow(L,d);
    }
    
    template <class Tc,class comp>
    inline polynomial_<Tc,comp> resultant(const polynomial_<Tc,comp> &G,const polynomial_<Tc,comp> & F,const variable & v)
    {
        assert(comp_consistent(G.comp(),F.comp()));
        univariate_priority_order comp_v(v);
        polynomial_<Tc,univariate_priority_order>  G1(&comp_v);
        polynomial_<Tc,univariate_priority_order>  F1(&comp_v);
        poly_convert(G,G1);poly_convert(F,F1);
        polynomial_<Tc,univariate_priority_order>  O1(&comp_v);
        resultant(O1,G1,F1);
        polynomial_<Tc,comp> O(G.comp_ptr());
        poly_convert(std::move(O1),O);
        return O;
    }
    
    
    template <class Tc>
    void resultant
        (   
            polynomial_<Tc,univariate_priority_order>&O,
            const polynomial_<Tc,univariate_priority_order>&F,
            const polynomial_<Tc,univariate_priority_order>&G
            
        )
    {
        O.clear();
        const univariate_priority_order &comp=G.comp();
        int64_t l=get_up_deg(G);
        int64_t m=get_up_deg(F);
        if (m<0 || l<0)
            return void();
        if (m<l)
        {
            resultant(O,G,F);
            if ((l&1) &&(m&1))
                pair_vec_negate(O.data());
            return void();
        }
        if (m==0)
        {
            O=pow(F,l);
            return void();
        }
        if (l==0)
        {
            O=pow(G,m);
            return  void();
        }

        polynomial_<Tc,univariate_priority_order> S_j_1(&comp);
        polynomial_<Tc,univariate_priority_order> S_j(&comp);
        polynomial_<Tc,univariate_priority_order> S_r_1(&comp);
        polynomial_<Tc,univariate_priority_order> S_r(&comp);
        
        polynomial_<Tc,univariate_priority_order> R_(&comp);
        polynomial_<Tc,univariate_priority_order> tmp1(&comp);
        polynomial_<Tc,univariate_priority_order> tmp2(&comp);
        int64_t j,r,u;
        // std::cout<<"F="<<F<<std::endl;
        // std::cout<<"G="<<G<<std::endl;
        
        if (l<m)
        {
            j=m-1;
            if(j==l)
            {
                prem(S_j,F,G);
                --j;
                //  std::cout<<j<<":"<<S_j<<std::endl;
                S_j_1=G;
            }
            else{
                leadcoeff(R_,G);
                if (j-l<2)
                    pair_vec_multiplies(S_r.data(),R_.data(),G.data(),comp);
                else 
                {
                    pair_vec_power(tmp1.data(),R_.data(),j-l,comp);
                    pair_vec_multiplies(S_r.data(),tmp1.data(),G.data(),comp);
                }
                //  std::cout<<l<<":"<<S_r<<std::endl;
                prem(S_r_1,F,G);
                if (j-l & 1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<l-1<<":"<<S_r_1<<std::endl;
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
                j=l-1;
            }

        }else
        {
            j=m;
            prem(S_j,F,G);
            pair_vec_negate(S_j.data());
            //  std::cout<<j-1<<":"<<S_j<<std::endl;
            if (!(--j))
            {
                O.data()=std::move(S_j.data());
                return void();
            }
            r=get_up_deg(S_j);
            if (r<0)
                return void();
            if (r<j)
            {
                
                if (r!=1)
                {
                    leadcoeff(R_,S_j);
                    if (j-r<2)
                        pair_vec_multiplies(S_r.data(),R_.data(),S_j.data(),comp);
                    else 
                    {
                        pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                        pair_vec_multiplies(S_r.data(),tmp1.data(),S_j.data(),comp);
                    }
                }
                
                //  std::cout<<r<<":"<<S_r<<std::endl;
                if (!r)
                {
                    O.data()=std::move(S_r.data());
                    return void();
                }
                prem(tmp1,G,S_j);
                //std::cout<<"G"<<":"<<G<<std::endl;
                //std::cout<<"S_"<<j<<":"<<S_j<<std::endl;
                
                //  std::cout<<r-1<<":"<<tmp1<<std::endl;
                leadcoeff(R_,G);
                pair_vec_div(S_r_1.data(),tmp1.data(),R_.data(),comp);
                if ((j-r) & 1==1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<r-1<<":"<<S_r_1<<std::endl;
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
                j=r-1;
            }
            else
            {
                prem(tmp1,G,S_j);
                leadcoeff(R_,G);
                pair_vec_div(S_j_1.data(),tmp1.data(),R_.data(),comp);
                //  std::cout<<j-1<<":"<<S_j_1<<std::endl;
                --j;
                swap(S_j.data(),S_j_1.data());
            }
        }
        while (j>0)
        {
            r=get_up_deg(S_j);
            if (r<0)
                return void();
            if (r<j)
            {
                
                //std::cout<<"R_j"<<":"<<tmp1<<std::endl;
                leadcoeff(R_,S_j_1);
                //std::cout<<"R_j+1"<<":"<<R_<<std::endl;
                if (r!=1)
                {
                    leadcoeff(tmp1,S_j);
                    if (j-r<2)
                    {
                        pair_vec_multiplies(tmp2.data(),tmp1.data(),S_j.data(),comp);
                        pair_vec_div(S_r.data(),tmp2.data(),R_.data(),comp);
                    }
                    else 
                    {
                        pair_vec_power(S_r.data(),tmp1.data(),j-r,comp);
                        pair_vec_multiplies(tmp2.data(),S_r.data(),S_j.data(),comp);
                        pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                        pair_vec_div(S_r.data(),tmp2.data(),tmp1.data(),comp);
                    }
                }
                else
                {
                    pair_vec_power(tmp1.data(),R_.data(),j-r,comp);
                }
                
                //  std::cout<<r<<":"<<S_r<<std::endl;
                if (!r)
                {
                    O.data()=std::move(S_r.data());
                    return void();
                }
                pair_vec_multiplies(tmp2.data(),R_.data(),R_.data(),comp);
                if (j-r<2)
                {
                    pair_vec_multiplies(tmp1.data(),R_.data(),tmp2.data(),comp);    
                    prem(tmp2,S_j_1,S_j);
                    pair_vec_div(S_r_1.data(),tmp2.data(),tmp1.data(),comp);    
                }        
                else
                {
                    pair_vec_multiplies(R_.data(),tmp2.data(),tmp1.data(),comp);
                    prem(tmp2,S_j_1,S_j);
                    pair_vec_div(S_r_1.data(),tmp2.data(),R_.data(),comp);
                }
                if ((j-r) & 1)
                    pair_vec_negate(S_r_1.data());
                //  std::cout<<r-1<<":"<<S_r_1<<std::endl;
                swap(S_j.data(),S_r_1.data());
                swap(S_j_1.data(),S_r.data());
                j=r-1;
            }
            else
            {
                prem(tmp1,S_j_1,S_j);
                // std::cout<<"tmp1:"<<tmp1<<std::endl;
                leadcoeff(R_,S_j_1);
                // std::cout<<"R_:"<<R_<<std::endl;
                pair_vec_multiplies(tmp2.data(),R_.data(),R_.data(),comp);
                // std::cout<<"tmp2:"<<tmp2<<std::endl;
                pair_vec_div(S_j_1.data(),tmp1.data(),tmp2.data(),comp);
                --j;
                //  std::cout<<j<<":"<<S_j_1<<std::endl;
                swap(S_j.data(),S_j_1.data());
            }
            
        }
        O.data()=std::move(S_j.data());
        
    }

    
    template <class compare>
    polynomial_<Zp,compare> polynomial_mod(const polynomial_<ZZ,compare> & p, uint32_t prime)
    {
        polynomial_<Zp,compare> new_p(p.comp_ptr());
        // std::cout<<p<<" mod "<<prime;
        Zp coeff(prime);
        for (auto & i:p)
        {
            coeff=i.second; 
            if (coeff)
                new_p.push_back({i.first,std::move(coeff)});
        }
        // std::cout<<" ="<<new_p<<std::endl;
        return new_p;
    }
    template <class compare>
    polynomial_<Zp,compare> polynomial_mod(polynomial_<ZZ,compare> && p, uint32_t prime)
    {
        polynomial_<Zp,compare> new_p(p.comp_ptr());
        Zp coeff(prime);
        for (auto & i:p)
        {
            coeff=i.second; 
            if (coeff)
            new_p.push_back({std::move(i.first),std::move(coeff)});
        }
        p.clear();
        return new_p;
    }

    template <class Tc, class To>
    polynomial_<Tc,To> association(const polynomial_<Tc,To>& p,const variable & v,const Tc & c)
    {
        polynomial_<Tc,To> Pout(p.comp_ptr());
        basic_monomial<To> m(p.comp_ptr());
        basic_monomial<To> m1(p.comp_ptr());
        bool f=true;
        Tc z, z1;
        for (auto &i:p)
        {

            m.clear();m.reserve(i.first.size());
            z=i.second;
            for (auto& j:i.first)
                if (j.first==v)
                {
                    z*=pow(c,j.second);
                }
                else
                    m.push_back(j);
            if (!f && m!=m1)
            {
                if (z1)
                    Pout.push_back({std::move(m1),std::move(z1)});
                m1=std::move(m);
                z1=std::move(z);
            }
            else if (f)
            {
                m1=std::move(m);
                z1=std::move(z);
                f=false;
            }
            else
            {
                z1+=z;
            }

            
        }
        if (z1)
            Pout.push_back({std::move(m1),std::move(z1)});
        return Pout;
    }
    



}
#endif