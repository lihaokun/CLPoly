/*
Module Name:
    upolynomial.hh
Abstract:
    定义upolynomial
Author:
    haokun li
*/
#ifndef CLPOLY_UPOLYNOMIAL_HH
#define CLPOLY_UPOLYNOMIAL_HH
#include <clpoly/polynomial.hh>
namespace clpoly{
    class umonomial
    {
        private:
            int64_t __deg;
        public:
            umonomial():__deg(0){}
            umonomial(int64_t d):__deg(d){}
            constexpr int64_t deg() const {return this->__deg;}
            constexpr int64_t & deg() {return this->__deg;}
            constexpr bool empty() const  {return __deg==0;}
            constexpr operator int64_t () const{return this->__deg;}
            friend inline umonomial operator*  (umonomial p1,umonomial p2)
            {
                return umonomial(p1.deg()+p2.deg());
            }
            friend inline umonomial operator/  (umonomial p1,umonomial p2)
            {
                return umonomial(p1.deg()-p2.deg());
            }       
            std::string str() const
            {
                std::ostringstream ss;
                ss<<(*this);
                return ss.str();
            }
            friend inline std::ostream& operator<<  (std::ostream& stream, umonomial v) {
                if (!v)
                {
                    stream<<"1";
                }
                else
                {
                    stream<<"#";
                    if (v.deg()!=1)
                        stream<<"^"<<v.deg();
                }
                return stream;
            }
    };
    umonomial pow(umonomial m,int64_t i)
    {return umonomial(m.deg()*i);}
    struct uless
    {
        static uless init;
        constexpr bool operator()(const variable & v1,const variable & v2) const 
        {
            return  v1<v2; 
        }
        constexpr bool operator()(const umonomial & v1,const umonomial & v2) const 
        {
            return  v1.deg()>v2.deg(); 
        }
        constexpr bool operator==( uless g1)const {return true;}
        constexpr bool operator!=(uless g1)const {return false;}
    };
    uless uless::init;
    template <class coeff>
    using upolynomial_=basic_polynomial<umonomial,coeff,uless>;

    template<class Tc>
    Tc association(const upolynomial_<Tc> &P,const Tc & a)
    {
        Tc O=0;
        for (auto &i:P)
        {
            O+=i.second*pow(a,i.first.deg());
        }
        return O;
    }
    template<class Tc,class Tc2,class Tc3>
    Tc association(const upolynomial_<Tc2> &P,const Tc3 & a)
    {
        Tc O=0;
        for (auto &i:P)
        {
            O+=i.second*pow(a,i.first.deg());
        }
        return O;
    }
    
}
#endif