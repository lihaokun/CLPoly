/**
 * @file interval.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义interval class
*/

#ifndef CLPOLY_INTERVAL_HH
#define CLPOLY_INTERVAL_HH
#include <mpria.h>
#include <sstream>
#include <iostream>
#include <clpoly/number.hh>
#include <clpoly/polynomial_.hh>
namespace clpoly{
    class interval
    {
        private:
            mpri_t _I;

            inline static std::ostream& _to_stream(std::ostream& stream, const mpq_t q)
            {
                if (mpria_mpq_is_nan(q))
                {
                    stream<<"NAN";
                }
                else if (mpria_mpq_is_infinite(q)){
                    if (mpria_mpq_sgn(q)<0)
                        stream<<"-";
                    stream<<"INF";
                }
                else{
                    stream<<q;
                }
                return stream;
            }
        public:
            inline interval(){
                mpri_init(_I);
                MPRIA_MPQ_SET_NEG_INF(mpri_lepref(_I));
                MPRIA_MPQ_SET_POS_INF(mpri_repref(_I));
            }

            inline interval(const interval& I){
                mpri_init(_I);
                mpri_set(this->_I,I._I);
            }
            inline interval(interval&& I){
                mpri_init(_I);
                mpri_swap(this->_I,I._I);
            }
            inline interval(const QQ& q){
                mpri_init(_I);
                MPRI_SET_Q(_I,q.get_mpq_t());
            }

            inline  virtual ~interval(){
                mpri_clear(_I);
            }
            inline  interval& operator=(interval I)
            {
                if (I._I!=this->_I)
                    mpri_swap(this->_I,I._I);
                return *this;
            }
            inline  interval& operator=(const QQ& q)
            {
                MPRI_SET_Q(_I,q.get_mpq_t());
                return *this;
            }
        
            inline void set_l(const QQ&q)
            {
                mpq_set(mpri_lepref(this->_I),q.get_mpq_t());
            }
            inline void set_r(const QQ&q)
            {
                mpq_set(mpri_repref(this->_I),q.get_mpq_t());
            }

            inline bool is_zero() const
            {
                return mpri_is_zero(_I);
            }
            inline bool is_notzero() const
            {
                return mpri_is_nonzero(_I);
            }
            
            inline operator bool() const
            {
                return (this->is_notzero());
            }

            inline  bool operator==(const interval& I) const
            {
                return mpri_equal(this->_I,I._I);
            }
            inline  interval operator+(const interval& I) const
            {
                interval ans;
                mpri_add(ans._I,this->_I,I._I);
                return ans;
            }
            inline  interval& operator+=(const interval& I) 
            {
                interval ans;
                mpri_add(ans._I,this->_I,I._I);
                *this=std::move(ans);
                return *this;
            }
            
            
            inline  interval operator-(const interval& I) const
            {
                interval ans;
                mpri_sub(ans._I,this->_I,I._I);
                return ans;
            }

            inline  interval& operator-=(const interval& I) 
            {
                interval ans;
                mpri_sub(ans._I,this->_I,I._I);
                *this=std::move(ans);
                return *this;
            }


            inline  interval operator*(const interval& I) const
            {
                interval ans;
                mpri_mul(ans._I,this->_I,I._I);
                return ans;
            }
            inline  interval& operator*=(const interval& I) 
            {
                interval ans;
                mpri_mul(ans._I,this->_I,I._I);
                *this=std::move(ans);
                return *this;
            }


            inline  interval operator/(const interval& I) const
            {
                interval ans;
                mpri_div(ans._I,this->_I,I._I);
                return ans;
            }
            inline  interval& operator/=(const interval& I) 
            {
                interval ans;
                mpri_div(ans._I,this->_I,I._I);
                *this=std::move(ans);
                return *this;
            }

            inline interval sqr() const
            {
                interval ans;
                mpri_sqr(ans._I,this->_I);
                return ans;
            }
            
            inline interval pow(size_t n) const
            {
                if (n==1)
                    return *this;
                if (n==2)
                    return this->sqr();
                interval ans(1);
                if (n==0)
                    return ans;
                interval tmp(*this);
                size_t index=1;
                while (1){
                    if (index&n)
                    {
                        ans=ans*tmp;
                    }
                    index<<=1;
                    if (index>n)
                        break;
                    tmp=tmp.sqr();
                    if (index==n)
                        return tmp;
                }
                return ans;
            }
            


            friend std::ostream& operator<<(std::ostream& stream, const interval& I)
            {
                _to_stream(_to_stream(stream<<"[",mpri_lepref(I._I))<<",",mpri_repref(I._I))<<"]";
                return stream;
            }

            
            
    };

    template <class comp=grlex>
    using interval_poly=polynomial_<interval,comp>;


    template <class Tc, class To>
    interval_poly<To> assign(const polynomial_<Tc,To>& p,const std::map<variable,interval> & ass_list)
    {
        interval_poly<To> Pout(p.comp_ptr());
        basic_monomial<To> m(p.comp_ptr());
        basic_monomial<To> m1(p.comp_ptr());
        bool f=true;
        interval z, z1;
        for (auto &i:p)
        {

            m.clear();m.reserve(i.first.size());
            z=i.second;
            for (auto& j:i.first)
            {
                auto ptr=ass_list.find(j.first);
                if (ptr!=ass_list.end())
                {
                    z*=ptr->second.pow(j.second);
                }
                else
                    m.push_back(j);
            }
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
        std::cout<<Pout<<std::endl;
        Pout.normalization();
        return Pout;
    }


    template <> 
    struct zore_check<interval>
    {
        inline bool operator()(const interval & I)
        {
            return I.is_zero();
        } 
    };

    template <class To>
    std::ostream& operator<<  (std::ostream& stream, const interval_poly<To>& v) 
    {
        if (v.size()==0)
            return stream<<'0';
        bool is_print=false;
        for(auto i=v.begin();i!=v.end();++i)
        {
            if (!zore_check<basic_monomial<To>>()(i->first) //&& !zore_check<Tc>()(i->second)
                ){
                if (is_print) //greater
                    stream<<"+";
                is_print=true;
                stream<< i->second<<"*";
                stream<<i->first;
                continue;
            }
            if (is_print) //greater
                stream<<"+";
            is_print=true;
            stream<<i->second;
        }
        return stream;
    }
}
#endif