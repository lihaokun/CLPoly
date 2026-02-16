/**
 * @file number.hh
 * @author 李昊坤 (ker@pm.me)
 * @brief 关于数域的定义
 *
 */
#ifndef CLPOLY_NUMBER_HH
#define CLPOLY_NUMBER_HH
#include <clpoly/number/ZZ.hh>
#include <clpoly/number/QQ.hh>
#include <clpoly/basic.hh>
#include <cmath>
#include <cstdint>
#include <cassert>

namespace clpoly{
    // ---- ZZ template specializations ----
    template<>
    inline void addmul(ZZ &op,const ZZ &op1,const ZZ&op2)
    {
        op.addmul(op1, op2);
    }
    template<>
    inline void submul(ZZ &op,const ZZ &op1,const ZZ&op2)
    {
        op.submul(op1, op2);
    }
    template<>
    inline void __div(ZZ &op,const ZZ &op1,const ZZ&op2)
    {
        ZZ::fdiv_q(op, op1, op2);
    }
    template<>
    inline void __div(ZZ &op,ZZ &op_r,const ZZ &op1,const ZZ&op2)
    {
        ZZ::fdiv_qr(op, op_r, op1, op2);
    }
    template<>
    inline void __div(QQ &op,const ZZ &op1,const ZZ&op2)
    {
        op=QQ(op1,op2);
    }
    template<>
    inline void __div(QQ &op,QQ &op_r,const ZZ &op1,const ZZ&op2)
    {
        op=QQ(op1,op2);
        op_r=QQ(0);
    }
    template <>
    inline void set_zero(ZZ& op)
    {
        op=0LL;
    }
    template<>
    struct zore_check<ZZ>
    {
        bool operator()(const ZZ & op)
        {
            return !op;
        }
    };
    template<>
    struct zore_check<QQ>
    {
        bool operator()(const QQ & op)
        {
            return !op;
        }
    };

    // ---- Zp ----
    inline uint64_t inv_prime(uint64_t _i,uint32_t _p)
    {
        assert(_p!=0 && _i!=0);
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
        Zp(uint32_t i,uint32_t p):_i(i%p),_p(p){}
        Zp(int64_t i,uint32_t p):_i(i>=0?i%p:p-(-i)%p),_p(p){}
        Zp(int i,uint32_t p):_i(i>=0?i%p:p-(-i)%p),_p(p){}
        Zp(const ZZ& i,uint32_t p):_i(i.fdiv_ui(p)),_p(p){}
        inline Zp inv() const
        {
            //assert(this->_p!=0 && this->_i!=0);
            Zp new_op(this->_p);
            new_op._i=inv_prime(this->_i,this->_p);
            return new_op;
        }
        constexpr Zp& operator=(int64_t i)
        {
            assert(this->_p!=0);
            this->_i=i>=0?i%this->_p:this->_p-(-i)%this->_p;
            return *this;
        }
        inline Zp& operator=(const ZZ& i)
        {
            assert(this->_p!=0);
            this->_i=i.fdiv_ui(this->_p);
            return *this;
        }
        explicit constexpr operator std::uint64_t() const {return this->_i;}
        explicit constexpr operator bool() const {return this->_i != 0;}
        constexpr uint32_t prime() const {return this->_p;}
        constexpr uint32_t & prime() {return this->_p;}
        constexpr uint64_t number() const {return this->_i;}
        constexpr uint64_t & number() {return this->_i;}
        constexpr void normalization(){assert(this->_p);this->_i%=this->_p;}
        constexpr void prime(uint32_t p) {this->_p=p;}
        constexpr Zp & operator+()
        { return *this;}
        Zp operator-() const
        {
            return Zp(this->_p-this->_i,this->_p);
        }
        friend inline Zp operator+(Zp op1,const Zp &  op2)
        {
            assert(op1._p==op2._p);
            op1._i+=op2._i;
            op1._i%=op1._p;
            return op1;
        }
        inline Zp & operator+=(const Zp &  op2)
        {
            assert(this->_p==op2._p);
            this->_i+=op2._i;
            this->_i%=this->_p;
            return *this;
        }
        friend inline Zp operator-(Zp op1,const Zp &  op2)
        {
            assert(op1._p==op2._p);
            op1._i+=op1._p-op2._i;
            op1._i%=op1._p;
            return op1;
        }
        friend inline Zp operator*(Zp op1,const Zp &  op2)
        {
            assert(op1._p==op2._p);
            op1._i*=op2._i;
            op1._i%=op1._p;
            return op1;
        }
        friend inline Zp operator*(Zp op1, int64_t op2)
        {
            Zp tmp(op2, op1._p);
            op1._i*=tmp._i;
            op1._i%=op1._p;
            return op1;
        }
        friend inline Zp operator*(int64_t op1, const Zp & op2)
        {
            return op2 * op1;
        }
        inline Zp & operator*=(const Zp &  op2)
        {
            assert(this->_p==op2._p);
            this->_i*=op2._i;
            this->_i%=this->_p;
            return *this;
        }
        friend inline Zp operator/(Zp op1,const Zp & op2)
        {
            assert(op1._p==op2._p);
            op1._i*=inv_prime(op2._i,op1._p);
            op1._i%=op1._p;
            return op1;
        }
        inline Zp & operator/=(const Zp & op2)
        {
            assert(this->_p==op2._p);
            this->_i*=inv_prime(op2._i,this->_p);
            this->_i%=this->_p;
            return *this;
        }
        friend inline bool operator==(const Zp& op1, const Zp& op2)
        {
            assert(op1._p==op2._p);
            return op1._i==op2._i;
        }
        friend inline bool operator!=(const Zp& op1, const Zp& op2)
        {
            assert(op1._p==op2._p);
            return op1._i!=op2._i;
        }
        friend inline bool operator==(const Zp& op1, int64_t op2)
        {
            Zp tmp(op2, op1._p);
            return op1._i==tmp._i;
        }
        friend inline bool operator==(int64_t op1, const Zp& op2)
        {
            return op2==op1;
        }
        friend inline bool operator!=(const Zp& op1, int64_t op2)
        {
            return !(op1==op2);
        }
        friend inline bool operator!=(int64_t op1, const Zp& op2)
        {
            return !(op2==op1);
        }
        friend inline bool operator>=(const Zp& op1, int64_t op2)
        {
            if (op2<=0) return true;
            return op1._i>=(uint64_t)op2;
        }
        friend std::ostream& operator<<  (std::ostream& stream, const Zp& v)
        {
            stream << v._i;
            return stream;
        }
    };
    inline std::size_t hash_value(Zp  const& v)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, v.number());
        boost::hash_combine(seed, v.prime());
        return seed;
    }
    template<>
    struct zore_check<Zp>
    {
        bool operator()(const Zp & op)
        {
            return op.number() == 0;
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
        op.number()+=op.prime()-(op1.number()*op2.number())%op.prime();
        op.number()%=op.prime();
    }
    inline Zp pow(const Zp & z,int64_t i)
    {
        Zp o(1,z.prime());
        Zp z_=z;
        while (i>0)
        {
            if (i %2)
                o*=z_;
            i>>=1;
            if (i)
                z_*=z_;
        }
        return o;
    }
    template<>
    inline void __div(Zp &op,const Zp &op1,const Zp&op2)
    {
        op=op1/op2;
    }

    template<>
    inline void __div(Zp &op,Zp &op_r,const Zp &op1,const Zp&op2)
    {
        op=op1/op2;
        op_r=Zp(0,op.prime());
    }

}
#endif
