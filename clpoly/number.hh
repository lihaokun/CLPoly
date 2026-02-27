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
    inline uint64_t inv_prime(uint64_t _i, uint64_t _p)
    {
        assert(_p != 0 && _i != 0);
        uint64_t a = _p, b = _i, c;
        uint64_t s1 = 0, s2 = 1, s3;
        while ((c = a % b))
        {
            uint64_t q = a / b;
            uint64_t sq_mod = (unsigned __int128)s2 * q % _p;
            s3 = s1 >= sq_mod ? s1 - sq_mod : _p - sq_mod + s1;
            a = b; b = c;
            s1 = s2; s2 = s3;
        }
        return s2;
    }
    class Zp
    {
    private:
        uint64_t _i;      // 值 ∈ [0, p)
        uint64_t _p;      // 素数 p（0 = 未初始化哨兵）
        uint64_t _ninv;   // FLINT preinvert(p << _norm)
        uint32_t _norm;   // __builtin_clzll(p)（p=0 时为 0）

        inline static uint64_t _s_prime = 0;  // 类级别当前素数，默认 0（未设置）

        // FLINT 归一化逆元（pn 必须最高位为 1）
        static uint64_t __preinvert_limb(uint64_t pn)
        {
            assert(pn >> 63);
            unsigned __int128 num = ((unsigned __int128)(~pn) << 64) | ~(uint64_t)0;
            return (uint64_t)(num / pn);
        }

        // 预计算 ninv 和 norm（替换 __barrett_ninv）
        static void __precompute(uint64_t p, uint64_t& ninv, uint32_t& norm)
        {
            if (p == 0) { ninv = 0; norm = 0; return; }
            assert(p >= 2 && p < (1ULL << 63));
            norm = (uint32_t)__builtin_clzll(p);
            ninv = __preinvert_limb(p << norm);
        }

        // FLINT 归一化 Barrett 乘法归约（替换 __barrett_reduce）
        uint64_t __nmod_mul(uint64_t a, uint64_t b) const
        {
            assert(_p != 0);
            uint64_t pn = _p << _norm;
            uint64_t a_shifted = a << _norm;
            unsigned __int128 prod = (unsigned __int128)a_shifted * b;
            uint64_t hi = (uint64_t)(prod >> 64);
            uint64_t lo = (uint64_t)prod;
            unsigned __int128 qm = (unsigned __int128)hi * _ninv;
            uint64_t q1 = (uint64_t)(qm >> 64);
            uint64_t q0 = (uint64_t)qm;
            q0 += lo;
            q1 += hi + (q0 < lo ? 1 : 0);
            uint64_t r = lo - (q1 + 1) * pn;
            if (r > q0) r += pn;
            if (r >= pn) r -= pn;
            return r >> _norm;
        }

    public:
        // 类级别素数管理
        static void     set_prime(uint64_t p) { assert(p >= 2 && p < (1ULL << 63)); _s_prime = p; }
        static uint64_t cur_prime()           { return _s_prime; }

        Zp() : _i(0), _p(_s_prime) { __precompute(_p, _ninv, _norm); }
        explicit Zp(uint64_t p) : _i(0), _p(p) { __precompute(p, _ninv, _norm); }
        Zp(uint64_t i, uint64_t p) : _p(p) { __precompute(p, _ninv, _norm); _i = i % p; }
        Zp(int64_t i, uint64_t p) : _p(p)
        {
            __precompute(p, _ninv, _norm);
            int64_t r = i % (int64_t)p;
            _i = (uint64_t)(r >= 0 ? r : r + (int64_t)p);
        }
        Zp(int i, uint64_t p) : Zp((int64_t)i, p) {}
        Zp(const ZZ& i, uint64_t p) : _p(p)
        {
            __precompute(p, _ninv, _norm);
            _i = i.fdiv_ui(p);
        }

        inline Zp inv() const
        {
            Zp new_op(this->_p);
            new_op._i = inv_prime(this->_i, this->_p);
            return new_op;
        }
        constexpr Zp& operator=(int64_t i)
        {
            if (this->_p == 0) { assert(i == 0); this->_i = 0; return *this; }
            int64_t r = i % (int64_t)this->_p;
            this->_i = (uint64_t)(r >= 0 ? r : r + (int64_t)this->_p);
            return *this;
        }
        inline Zp& operator=(const ZZ& i)
        {
            if (this->_p == 0) { assert(!i); this->_i = 0; return *this; }
            this->_i = i.fdiv_ui(this->_p);
            return *this;
        }
        explicit constexpr operator std::uint64_t() const { return this->_i; }
        explicit constexpr operator bool() const { return this->_i != 0; }
        constexpr uint64_t  prime() const { return this->_p; }
        constexpr uint64_t  number() const { return this->_i; }
        constexpr uint64_t& number()       { return this->_i; }
        constexpr void normalization() { assert(this->_p); this->_i %= this->_p; }
        void prime(uint64_t p) { this->_p = p; __precompute(p, _ninv, _norm); }

        constexpr Zp& operator+() { return *this; }
        Zp operator-() const
        {
            Zp r(*this);
            r._i = (_i == 0) ? 0 : _p - _i;
            return r;
        }

        friend inline Zp operator+(Zp op1, const Zp& op2)
        {
            assert(op1._p == op2._p);
            uint64_t r = op1._i + op2._i;
            op1._i = (r >= op1._p) ? r - op1._p : r;
            return op1;
        }
        inline Zp& operator+=(const Zp& op2)
        {
            assert(this->_p == op2._p);
            uint64_t r = this->_i + op2._i;
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
            op1._i = op1.__nmod_mul(op1._i, op2._i);
            return op1;
        }
        friend inline Zp operator*(Zp op1, int64_t op2)
        {
            op1._i = op1.__nmod_mul(op1._i, Zp(op2, op1._p)._i);
            return op1;
        }
        friend inline Zp operator*(int64_t op1, const Zp& op2) { return op2 * op1; }
        inline Zp& operator*=(const Zp& op2)
        {
            assert(this->_p == op2._p);
            this->_i = this->__nmod_mul(this->_i, op2._i);
            return *this;
        }
        friend inline Zp operator/(Zp op1, const Zp& op2)
        {
            assert(op1._p == op2._p);
            op1._i = op1.__nmod_mul(op1._i, inv_prime(op2._i, op1._p));
            return op1;
        }
        inline Zp& operator/=(const Zp& op2)
        {
            assert(this->_p == op2._p);
            this->_i = this->__nmod_mul(this->_i, inv_prime(op2._i, this->_p));
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
            return op1._i >= (uint64_t)op2;
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
