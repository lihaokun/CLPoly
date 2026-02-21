/**
 * @file polynomial_convert.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义poly_convert
 * 
 */
#ifndef CLPOLY_POLYNOMIAL_CONVERT__HH
#define CLPOLY_POLYNOMIAL_CONVERT__HH
#include <clpoly/monomial.hh>
#include <clpoly/number.hh>
#include <clpoly/basic_polynomial.hh>
#include <clpoly/polynomial_type.hh>
namespace clpoly{
    template<class Tc,class comp>
    Tc tonum(const polynomial_<Tc,comp> & p)
    {
        if (!is_number(p))
            throw std::invalid_argument("clpoly::tonum::not num.");
        if (p.empty())
            return 0;
        return p.front().second;
    }

    template<class Tc,class Tb>
    void  poly_convert(Tc a,Tb p)= delete;
    template<class Tc,class comp>
    void  poly_convert(int a,polynomial_<Tc,comp> & p)
    {
        basic_monomial<comp> m(p.comp_ptr());
        
        p={{std::move(m),a}};
    }
    template<class Tc,class comp>
    void  poly_convert(variable v,polynomial_<Tc,comp> & p)
    {
        basic_monomial<comp> m(p.comp_ptr());
        m.push_back({v,1});
        p={{std::move(m),1}};
    }
    template<class Tc,class comp>
    void  poly_convert(basic_monomial<comp> m,polynomial_<Tc,comp> & p)
    {
        p={{std::move(m),1}};
    }
    template<class comp>
    void  poly_convert(QQ c,polynomial_<QQ,comp> & p)
    {
        basic_monomial<comp> m(p.comp_ptr());
        p={{std::move(m),c}};
    }
    template<class comp>
    void  poly_convert(ZZ c,polynomial_<ZZ,comp> & p)
    {
        basic_monomial<comp> m(p.comp_ptr());
        p={{std::move(m),c}};
    }
    template<class comp>
    void  poly_convert(ZZ c,polynomial_<QQ,comp> & p)
    {
        basic_monomial<comp> m(p.comp_ptr());
        p={{std::move(m),c}};
    }
    template<class T1,class T2,class comp1,class comp2>
    void poly_convert(const polynomial_<T1,comp1>& p_in,polynomial_<T2,comp2> & p_out)
    {
        if ((void*)&p_in==(void*)&p_out)
        {
            polynomial_<T1,comp1> tmp(p_in);
            poly_convert(tmp, p_out);
            return;
        }
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
        if ((void*)&p_in==(void*)&p_out)
        {
            polynomial_<T1,comp1> tmp(std::move(p_in));
            poly_convert(std::move(tmp), p_out);
            return;
        }
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
        assert((void*)&p_in!=(void*)&p_out);
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
    template<class comp1,class comp2>
    void poly_convert(const polynomial_<QQ,comp1>& p_in,polynomial_<ZZ,comp2> & p_out)
    {
        assert((void*)&p_in!=(void*)&p_out);
        p_out.clear();
        ZZ den=1;
        for (auto &i:p_in)
        {
            den=lcm(den,i.second.get_den());
        }
        basic_monomial<comp2> m(p_out.comp_ptr());
        for (auto &i:p_in)
        {
            m=i.first.data();
            p_out.push_back({std::move(m),i.second.get_num()*(den/i.second.get_den())});
        }
        p_out.normalization();
    }
    template<class comp1,class comp2>
    void poly_convert(polynomial_<QQ,comp1>&& p_in,polynomial_<ZZ,comp2> & p_out)
    {
        assert((void*)&p_in!=(void*)&p_out);
        p_out.clear();
        ZZ den=1;
        for (auto &i:p_in)
        {
            den=lcm(den,i.second.get_den());
        }
        basic_monomial<comp2> m(p_out.comp_ptr());
        for (auto &i:p_in)
        {
            m=std::move(i.first.data());
            p_out.push_back({std::move(m),i.second.get_num()*(den/i.second.get_den())});
        }
        p_out.normalization();
        p_in.clear();
    }
}
#endif