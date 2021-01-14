/**
 * @file polynomial_.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义polynomial_
 * 
 */
#ifndef CLPOLY_POLYNOMIAL__HH
#define CLPOLY_POLYNOMIAL__HH

#include <clpoly/monomial.hh>
#include <clpoly/number.hh>
#include <clpoly/basic_polynomial.hh>
namespace clpoly{
    template <class Tc,class comp=grlex>
    using polynomial_=basic_polynomial<basic_monomial<comp>,Tc,comp>;
    using polynomial_ZZ=polynomial_<ZZ>;
    using polynomial_QQ=polynomial_<QQ>;

    template<class Tc,class comp>
    void  convert(polynomial_<Tc,comp> & p,variable v)
    {
        p={{monomial(v),1}};
    }
    template<class Tc,class comp>
    void  convert(polynomial_<Tc,comp> & p,monomial m)
    {
        p={{std::move(m),1}};
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


    template <class Tc,class Tm,class compare>
    inline basic_polynomial<Tm,Tc,compare> pow(const  basic_polynomial<Tm,Tc,compare> & p,int64_t i)
    {
        assert(i>=0);
        basic_polynomial<Tm,Tc,compare> o(p.comp_ptr());
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
    monomial operator +(monomial  m)
    {
        return m;
    }
    polynomial_ZZ operator -(monomial  m)
    {
        return polynomial_ZZ({{std::move(m),-1}});
    }
    polynomial_ZZ operator+ (monomial m1,monomial m2)
    {
        if (polynomial_ZZ::compare_type()(m1,m2))
            return polynomial_ZZ({{std::move(m1),1},{std::move(m2),1}});
        return polynomial_ZZ({{std::move(m2),1},{std::move(m1),1}});
    }
    polynomial_ZZ operator- (monomial m1,monomial m2)
    {
        if (polynomial_ZZ::compare_type()(m1,m2))
            return polynomial_ZZ({{std::move(m1),1},{std::move(m2),-1}});
        return polynomial_ZZ({{std::move(m2),-1},{std::move(m1),1}});
    }

 
    
    polynomial_ZZ operator+(monomial m,int64_t i)
    {
        return polynomial_ZZ({{std::move(m),1},{{},i}});
    }
    
    polynomial_ZZ operator-(monomial  m,int64_t i)
    {
        return polynomial_ZZ({{std::move(m),1},{{},-i}});
    }

    polynomial_ZZ operator+(int64_t i,monomial  m)
    {
        return polynomial_ZZ({{std::move(m),1},{{},i}});
    }
    
    polynomial_ZZ operator-(int64_t i,monomial  m)
    {
        return polynomial_ZZ({{std::move(m),1},{{},-i}});
    }
    
    polynomial_ZZ operator*( monomial  m,int64_t i)
    {
        return polynomial_ZZ({{std::move(m),i}});
    }
    polynomial_ZZ operator*(int64_t i, monomial  m)
    {
        return polynomial_ZZ({{std::move(m),i}});
    }

    
    template<class Tc>
    polynomial_<Tc> operator- (polynomial_<Tc> O,const Tc & m)
    {

        return O+(-m);
    }
    template<class Tc>
    polynomial_<Tc> operator- (Tc  m,polynomial_<Tc> O)
    {

        return (-O)+m;
    }
    template<class Tc>
    polynomial_<Tc> operator- (polynomial_<Tc> p,int64_t m)
    {
        return p+Tc(-m);
    }
    template<class Tc>
    polynomial_<Tc> operator- (int64_t m, const polynomial_<Tc> &  p)
    {
        return (-p)+Tc(m);
    }
    template<class Tc>
    polynomial_<Tc> operator+ (const polynomial_<Tc> & p,int64_t m)
    {
        return p+Tc(m);
    }
    template<class Tc>
    polynomial_<Tc> operator+ (int64_t m,const polynomial_<Tc> & p)
    {
        return p+Tc(m);
    }
    template<class Tc>
    polynomial_<Tc> operator+ (polynomial_<Tc>  O,Tc  m)
    {
        if (O.empty() || O.back().first.empty())
            O.back().second+=m;
        else
            O.push_back({{},std::move(m)});
        return O;
    }
    template<class Tc>
    polynomial_<Tc> operator+ (Tc  m,polynomial_<Tc>  O)
    {
       return O+m;
    }
    template<class Tc>
    polynomial_<Tc> operator+ (polynomial_<Tc>  O,monomial m)
    {
        if (!O.empty() && !O.comp(O.back().first,m))
            return O+polynomial_<Tc>({{m,1}});
        else
            O.push_back({std::move(m),1});
        return O;
    }
    template<class Tc>
    polynomial_<Tc> operator+ (monomial m,polynomial_<Tc>  O)
    {
        return O+m;
    }
    template<class Tc>
    polynomial_<Tc> operator- (monomial m,const polynomial_<Tc>&  O)
    {
        return m+(-O);
    }
    template<class Tc>
    polynomial_<Tc> operator- (polynomial_<Tc>  O,monomial m)
    {
        if (!O.empty() && !O.comp(O.back().first,m))
            return O+polynomial_<Tc>({{m,-1}});
        else
            O.push_back({std::move(m),-1});
        return O;
    }
    template<class Tc,class comp>
    polynomial_<Tc,comp> operator* (polynomial_<Tc,comp>  O,Tc c)
    {
        if (!c)
            return polynomial_<Tc,comp>();
        for (auto & i:O)
            i.second*=c;
        O.normalization();
        return O;
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

}   


#endif