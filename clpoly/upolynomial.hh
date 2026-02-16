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
    inline umonomial pow(umonomial m,int64_t i)
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
    // uless uless::init;
    template <class coeff>
    using upolynomial_=basic_polynomial<umonomial,coeff,uless>;
    using upolynomial_ZZ=upolynomial_<ZZ>;
    inline bool is_divexact(umonomial & op,const umonomial & op1,const umonomial & op2)
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
    inline upolynomial_<Zp> polynomial_mod(const upolynomial_<ZZ> & p, uint32_t prime)
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
        // std::cout<<"poly_convert(const polynomial_<T1,comp1>& p_in,upolynomial_<T2> & p_out)"<<std::endl;
        // std::cout<<p_in<<std::endl;
        // std::cout<<p_out<<std::endl;
        
        p_out.clear();
        p_out.reserve(p_in.size());
        for(auto &i:p_in)
        {
            // std::cout<<i.first.deg()<<" "<<i.second<<std::endl;
        
            p_out.push_back(std::pair<umonomial,T2>(i.first.deg(),i.second));
            // std::cout<<p_out.back().first.deg()<<" "<<p_out.back().second<<std::endl;
        
        }
        // std::cout<<p_out<<std::endl;
        
        // std::cout<<"p_out:"<<p_out<<std::end;
        p_out.normalization();
    }
    void poly_convert(const upolynomial_<ZZ>& p_in,upolynomial_<QQ> & p_out);
    void poly_convert(const upolynomial_<QQ>& p_in,upolynomial_<ZZ> & p_out);
    template <class T>
    upolynomial_<T>  derivative(const upolynomial_<T> & p)
    {
        upolynomial_<T> Pout;
        for (auto &i:p)
        {
            if (i.first.deg())
            {
                T c = i.second * i.first.deg();
                if (!zore_check<T>()(c))
                    Pout.push_back(std::make_pair(umonomial(i.first.deg()-1), c));
            }
        }
        return Pout;
    }
    template <class Tc>
    constexpr bool is_number(const upolynomial_<Tc> & F)
    {
        return (F.empty() || F.size()==1 && F.front().first.empty());
    }
    // 标量乘法：upolynomial * scalar
    template<class Tc>
    upolynomial_<Tc> operator*(upolynomial_<Tc> O, const Tc& c) {
        if (!c) return upolynomial_<Tc>();
        for (auto& i : O)
            i.second *= c;
        O.normalization();
        return O;
    }

    // 标量乘法：scalar * upolynomial
    template<class Tc>
    upolynomial_<Tc> operator*(const Tc& c, upolynomial_<Tc> O) {
        return O * c;
    }

    template <class Tc>
    Tc cont(const upolynomial_<Tc> & F)
    {
        if (F.empty())
            return 1;
        auto ptr=F.begin();
        auto I=(ptr++).second;
        for (;ptr!=F.end();++ptr)
        {
            I=gcd(I,ptr->second);
        }
        return I;
    }
    


}
#endif