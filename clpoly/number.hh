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
        uint32_t _i;     // 值 ∈ [0, p)
        uint32_t _p;     // 素数 p（0 = 未初始化哨兵）
        uint64_t _ninv;  // Barrett 常数：UINT64_MAX / p（0 = 未初始化）

        inline static uint32_t _s_prime = 0;  // 类级别当前素数，默认 0（未设置）

        // 预计算 Barrett 常数
        // p==0：保持旧 Zp() 未初始化哨兵行为，addmul/submul 延迟初始化依赖此语义
        // p>=2：UINT64_MAX/p 等价于 floor(2^64/p)（product < p² < 2^64，整数误差=0）
        static uint64_t __barrett_ninv(uint32_t p)
        {
            if (p == 0) return 0;
            assert(p >= 2 && p < (1u << 31));  // p < 2^31：保证加法 _i+_i 不溢出 uint32_t
            return UINT64_MAX / p;
        }

        // Barrett 归约：将 product（< p²，即 < 2^64）归约到 [0, p)
        uint32_t __barrett_reduce(uint64_t product) const
        {
            assert(_p != 0);
            uint64_t q = (unsigned __int128)product * _ninv >> 64;
            uint64_t r = product - q * _p;
            return (uint32_t)(r >= _p ? r - _p : r);
        }

    public:
        // 类级别素数管理
        static void     set_prime(uint32_t p) { assert(p >= 2 && p < (1u << 31)); _s_prime = p; }
        static uint32_t cur_prime()           { return _s_prime; }

        Zp() : _i(0), _p(_s_prime), _ninv(__barrett_ninv(_s_prime)) {}
        explicit Zp(uint32_t p) : _i(0), _p(p), _ninv(__barrett_ninv(p)) {}
        Zp(uint64_t i, uint32_t p) : _p(p), _ninv(__barrett_ninv(p)) { _i = (uint32_t)(i % p); }
        Zp(uint32_t i, uint32_t p) : _p(p), _ninv(__barrett_ninv(p)) { _i = i % p; }
        Zp(int64_t i, uint32_t p) : _p(p), _ninv(__barrett_ninv(p))
        {
            // 修复原始 bug：i=-p 时 p - (-p)%p = p - 0 = p（越界）
            int64_t r = i % (int64_t)p;
            _i = (uint32_t)(r >= 0 ? r : r + (int64_t)p);
        }
        Zp(int i, uint32_t p) : _p(p), _ninv(__barrett_ninv(p))
        {
            int64_t r = (int64_t)i % (int64_t)p;
            _i = (uint32_t)(r >= 0 ? r : r + (int64_t)p);
        }
        Zp(const ZZ& i, uint32_t p) : _p(p), _ninv(__barrett_ninv(p)) { _i = (uint32_t)i.fdiv_ui(p); }

        inline Zp inv() const
        {
            Zp new_op(this->_p);
            new_op._i = (uint32_t)inv_prime(this->_i, this->_p);
            return new_op;
        }
        constexpr Zp& operator=(int64_t i)
        {
            if (this->_p == 0) { assert(i == 0); this->_i = 0; return *this; }
            int64_t r = i % (int64_t)this->_p;
            this->_i = (uint32_t)(r >= 0 ? r : r + (int64_t)this->_p);
            return *this;
        }
        inline Zp& operator=(const ZZ& i)
        {
            if (this->_p == 0) { assert(!i); this->_i = 0; return *this; }
            this->_i = (uint32_t)i.fdiv_ui(this->_p);
            return *this;
        }
        explicit constexpr operator std::uint64_t() const { return this->_i; }
        explicit constexpr operator bool() const { return this->_i != 0; }
        constexpr uint32_t prime() const { return this->_p; }
        // uint32_t& prime() 删除：_p 和 _ninv 必须同步，不能暴露可变引用
        constexpr uint32_t  number() const { return this->_i; }
        constexpr uint32_t& number()       { return this->_i; }
        constexpr void normalization() { assert(this->_p); this->_i %= this->_p; }
        void prime(uint32_t p) { this->_p = p; this->_ninv = __barrett_ninv(p); }

        constexpr Zp& operator+() { return *this; }
        Zp operator-() const
        {
            return Zp(_i == 0 ? 0u : _p - _i, _p);
        }

        friend inline Zp operator+(Zp op1, const Zp& op2)
        {
            assert(op1._p == op2._p);
            uint32_t r = op1._i + op2._i;
            op1._i = (r >= op1._p) ? r - op1._p : r;
            return op1;
        }
        inline Zp& operator+=(const Zp& op2)
        {
            assert(this->_p == op2._p);
            uint32_t r = this->_i + op2._i;
            this->_i = (r >= this->_p) ? r - this->_p : r;
            return *this;
        }
        friend inline Zp operator-(Zp op1, const Zp& op2)
        {
            assert(op1._p == op2._p);
            op1._i = op1._i >= op2._i ? op1._i - op2._i : op1._p - op2._i + op1._i;
            return op1;
        }
        inline Zp& operator-=(const Zp& op2)
        {
            assert(this->_p == op2._p);
            this->_i = this->_i >= op2._i ? this->_i - op2._i : this->_p - op2._i + this->_i;
            return *this;
        }
        friend inline Zp operator*(Zp op1, const Zp& op2)
        {
            assert(op1._p == op2._p);
            op1._i = op1.__barrett_reduce((uint64_t)op1._i * op2._i);
            return op1;
        }
        friend inline Zp operator*(Zp op1, int64_t op2)
        {
            op1._i = op1.__barrett_reduce((uint64_t)op1._i * Zp(op2, op1._p)._i);
            return op1;
        }
        friend inline Zp operator*(int64_t op1, const Zp& op2) { return op2 * op1; }
        inline Zp& operator*=(const Zp& op2)
        {
            assert(this->_p == op2._p);
            this->_i = this->__barrett_reduce((uint64_t)this->_i * op2._i);
            return *this;
        }
        friend inline Zp operator/(Zp op1, const Zp& op2)
        {
            assert(op1._p == op2._p);
            op1._i = op1.__barrett_reduce((uint64_t)op1._i * inv_prime(op2._i, op1._p));
            return op1;
        }
        inline Zp& operator/=(const Zp& op2)
        {
            assert(this->_p == op2._p);
            this->_i = this->__barrett_reduce((uint64_t)this->_i * inv_prime(op2._i, this->_p));
            return *this;
        }
        friend inline bool operator==(const Zp& op1, const Zp& op2)
        {
            assert(op1._p == op2._p);
            return op1._i == op2._i;
        }
        friend inline bool operator!=(const Zp& op1, const Zp& op2)
        {
            assert(op1._p == op2._p);
            return op1._i != op2._i;
        }
        friend inline bool operator==(const Zp& op1, int64_t op2)
        {
            Zp tmp(op2, op1._p);
            return op1._i == tmp._i;
        }
        friend inline bool operator==(int64_t op1, const Zp& op2)
        {
            return op2 == op1;
        }
        friend inline bool operator!=(const Zp& op1, int64_t op2)
        {
            return !(op1 == op2);
        }
        friend inline bool operator!=(int64_t op1, const Zp& op2)
        {
            return !(op2 == op1);
        }
        friend inline bool operator>=(const Zp& op1, int64_t op2)
        {
            if (op2 <= 0) return true;
            return op1._i >= (uint32_t)op2;
        }
        friend std::ostream& operator<<(std::ostream& stream, const Zp& v)
        {
            stream << v._i;
            return stream;
        }
    };
    inline std::size_t hash_value(Zp const& v)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, v.number());
        boost::hash_combine(seed, v.prime());
        return seed;
    }
    template<>
    struct zore_check<Zp>
    {
        bool operator()(const Zp& op)
        {
            return op.number() == 0;
        }
    };
    // set_zero<Zp>: 使用泛型模板（op = 0，即 operator=(int64_t)），无需特化
    template<>
    inline void addmul(Zp& op, const Zp& op1, const Zp& op2)
    {
        assert((op.prime() == op1.prime() || op.prime() == 0)
               && op1.prime() == op2.prime() && op1.prime());
        if (op.prime() == 0) op.prime(op1.prime());  // setter 同步更新 _p 和 _ninv
        op += op1 * op2;
    }
    template<>
    inline void submul(Zp& op, const Zp& op1, const Zp& op2)
    {
        assert((op.prime() == op1.prime() || op.prime() == 0)
               && op1.prime() == op2.prime() && op1.prime());
        if (op.prime() == 0) op.prime(op1.prime());  // setter 同步更新 _p 和 _ninv
        op -= op1 * op2;
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
