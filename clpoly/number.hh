/*
Module Name:
    number.hh
Abstract:
    关于高精度的定义
Author:
    haokun li
Notes:
*/
#ifndef CLPOLY_NUMBER_HH
#define CLPOLY_NUMBER_HH
#include <gmpxx.h>
#include <clpoly/basic.hh>
#include <cmath>
#include <cassert>
namespace clpoly{
    typedef mpz_class ZZ;
    typedef mpq_class QQ;
    // using ZZ=cln::cl_I;
    template<>
    inline void addmul(mpz_class &op,const mpz_class &op1,const mpz_class&op2)
    {
        mpz_addmul(op.get_mpz_t(),op1.get_mpz_t(),op2.get_mpz_t());
    }
    template<>
    inline void submul(mpz_class &op,const mpz_class &op1,const mpz_class&op2)
    {
        mpz_submul(op.get_mpz_t(),op1.get_mpz_t(),op2.get_mpz_t());
    }
    template<>
    inline void __div(mpz_class &op,const mpz_class &op1,const mpz_class&op2)
    {
        mpz_fdiv_q(op.get_mpz_t(),op1.get_mpz_t(),op2.get_mpz_t());
    }
    template<>
    inline void __div(mpz_class &op,mpz_class &op_r,const mpz_class &op1,const mpz_class&op2)
    {
        mpz_fdiv_qr(op.get_mpz_t(),op_r.get_mpz_t(),op1.get_mpz_t(),op2.get_mpz_t());
    }
    template<>
    inline void __div(mpq_class &op,const mpz_class &op1,const mpz_class&op2)
    {
        op=mpq_class(op1,op2);
    }
    template<>
    inline void __div(mpq_class &op,mpq_class &op_r,const mpz_class &op1,const mpz_class&op2)
    {
        op=mpq_class(op1,op2);
        op_r=0;
    }
    template <>
    void set_zero(mpz_class& op)
    {
        mpz_set_si(op.get_mpz_t(),0);
    }
    template<>
    struct zore_check<mpz_class>: public std::unary_function<mpz_class, bool>
    {
        bool operator()(const mpz_class & op)
        {
            return !op;
        } 
    };
    template<>
    struct zore_check<mpq_class>: public std::unary_function<mpq_class, bool>
    {
        bool operator()(const mpq_class & op)
        {
            return !op;
        } 
    };

    uint64_t inv_prime(uint64_t _i,uint32_t _p)
    {
        assert(_p!=0);
        uint64_t a=_p,b=_i,c;
        uint64_t s1=0,s2=1,s3;
        while (c=(a%b))
        {
            s3=(s1+_p-(s2*(a/b))%_p)%_p;
            a=b;b=c;
            s1=s2;s2=s3;
        }
        return s2;
    }
    class Zp
    {
    private:
        uint64_t _i;
        uint32_t _p; 
    public:
        Zp():_i(0),_p(0){}
        explicit  Zp(uint32_t p):_i(0),_p(p){}
        Zp(uint64_t i,uint32_t p):_i(i%p),_p(p){}
        Zp(int64_t i,uint32_t p):_i(i>0?i%p:p-(-i)%p),_p(p){}
        Zp(int i,uint32_t p):_i(i>0?i%p:p-(-i)%p),_p(p){}
        Zp(ZZ i,uint32_t p):_i(mpz_fdiv_ui(i.get_mpz_t(),p)),_p(p){}
        inline Zp inv() const
        {
            assert(this->_p!=0);
            Zp new_op(this->_p);
            new_op._i=inv_prime(this->_i,this->_p);
            return new_op;
        }
        constexpr Zp& operator=(int64_t i)
        {
            assert(this->_p!=0);
            this->_i=i>0?i%this->_p:this->_p-(-i)%this->_p;
            return *this;
        }
        inline Zp& operator=(const ZZ& i)
        {
            assert(this->_p!=0);
            this->_i=mpz_fdiv_ui(i.get_mpz_t(),this->_p);
            return *this;
        }
        constexpr operator std::uint64_t() const {return this->_i;}
        constexpr uint32_t prime() const {return this->_p;}
        constexpr uint32_t & prime() {return this->_p;}
        constexpr uint64_t number() const {return this->_i;}
        constexpr uint64_t & number() {return this->_i;}
        constexpr void normalization(){assert(this->_p);this->_i%=this->_p;}
        constexpr void prime(uint32_t p) {this->_p=p;}
        constexpr Zp & operator+()
        { return *this;}
        constexpr Zp & operator-()
        {
            this->_i=this->_p-this->_i;
            return *this;
        }
        friend inline Zp operator+(Zp op1,const Zp & op2)
        {
            assert(op1._p==op2._p);
            op1._i+=op2._i;
            op1._i%=op1._p;
            return op1;
        }
        friend inline Zp operator-(Zp op1,const Zp & op2)
        {
            assert(op1._p==op2._p);
            op1._i+=op1._p-op2._i;
            op1._i%=op1._p;
            return op1;
        }
        friend inline Zp operator*(Zp op1,const Zp & op2)
        {
            assert(op1._p==op2._p);
            op1._i*=op2._i;
            op1._i%=op1._p;
            return op1;
        }
        friend inline Zp operator/(Zp op1,const Zp & op2)
        {
            assert(op1._p==op2._p);
            op1._i*=inv_prime(op2._i,op1._p);
            op1._i%=op1._p;
            return op1;
        }
        friend std::ostream& operator<<  (std::ostream& stream, const Zp& v) 
        {
            stream << v._i;
            return stream;
        }
    };
    template<>
    struct zore_check<Zp>: public std::unary_function<Zp, bool>
    {
        bool operator()(const Zp & op)
        {
            return !op;
        } 
    };
    template<>
    inline void addmul(Zp &op,const Zp &op1,const Zp&op2)
    {
        assert((op.prime()==op1.prime()|| op.prime()==0) && op1.prime()==op2.prime() && op1.prime());
        op.prime()|=op1.prime();
        op.number()+=op1.number()*op2.number();
        op.number()%=op.prime();
    }
    template<>
    inline void submul(Zp &op,const Zp &op1,const Zp&op2)
    {
        assert((op.prime()==op1.prime()|| op.prime()==0) && op1.prime()==op2.prime() && op1.prime());
        op.prime()|=op1.prime();
        op.number()-=op1.number()*op2.number();
        op.number()%=op.prime();
    }
    // template<>
    // constexpr void __div(Zp &op,const Zp &op1,const Zp&op2)
    // {
    //     op=op1/op2;
    // }
    template<>
    inline void __div(Zp &op,Zp &op_r,const Zp &op1,const Zp&op2)
    {
        op=op1/op2;
        op_r=0;
    }

}
#endif