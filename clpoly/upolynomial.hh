/**
 * @file upolynomial.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 定义upolynomial
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
    using upolynomial_ZZ=upolynomial_<ZZ>;
    bool is_divexact(umonomial & op,const umonomial & op1,const umonomial & op2)
    {
        op=umonomial(op1.deg()-op2.deg());
        return op.deg()>=0;    
    }
    template<class Tc>
    Tc assign(const upolynomial_<Tc> &P,const Tc & a)
    {
        Tc O=0;
        for (auto &i:P)
        {
            O+=i.second*pow(a,i.first.deg());
        }
        return O;
    }
    template<class Tc,class Tc2,class Tc3>
    Tc assign(const upolynomial_<Tc2> &P,const Tc3 & a)
    {
        Tc O=0;
        for (auto &i:P)
        {
            O+=i.second*pow(a,i.first.deg());
        }
        return O;
    }
    template<class Tc>
    inline int64_t get_deg(const upolynomial_<Tc> &p)
    {
        return p.empty()?0:p.front().first.deg();
    }
    template <class Tc>
    int64_t degree(const upolynomial_<Tc> & p)
    {
        return p.degree();
    }
    upolynomial_<Zp> polynomial_mod(const upolynomial_<ZZ> & p, uint32_t prime)
    {
        upolynomial_<Zp> new_p;
        Zp coeff(prime);
        for (auto & i:p)
        {
            coeff=i.second; 
            if (coeff)
            new_p.push_back({i.first,std::move(coeff)});
        }
        return new_p;
    }
    template<class T1,class T2,class comp1>
    void poly_convert(const polynomial_<T1,comp1>& p_in,upolynomial_<T2> & p_out)
    {
        p_out.clear();
        p_out.reserve(p_in.size());
        for(auto &i:p_in)
        {
            p_out.push_back(std::pair<umonomial,T2>(i.first.deg(),i.second));
        }
        p_out.normalization();
    }
    template <class T>
    upolynomial_<T>  derivative(const upolynomial_<T> & p)
    {
        upolynomial_<T> Pout;
        int64_t b;
        for (auto &i:p)
        {
            if (i.first.deg())
            {
                Pout.push_back({umonomial(i.first.deg()-1),i.second*i.first.deg()});
            }    
        }
        return Pout;
    }
    
    


}
#endif