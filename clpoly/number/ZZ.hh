/**
 * @file ZZ.hh
 * @brief Custom arbitrary-precision integer with small integer optimization.
 *
 * 16-byte dual-field layout: int64_t _val + mpz_ptr _mpz.
 * When _mpz == nullptr the value lives in _val (zero heap allocation).
 * When _mpz != nullptr the value lives in the pointed-to mpz_t.
 */
#ifndef CLPOLY_NUMBER_ZZ_HH
#define CLPOLY_NUMBER_ZZ_HH

#include <gmp.h>
#include <cstdint>
#include <cstdlib>
#include <climits>
#include <iostream>
#include <string>
#include <functional>
#include <boost/container_hash/hash.hpp>

namespace clpoly {

class ZZ {
private:
    int64_t _val;
    mpz_ptr _mpz;

    // ---- internal helpers ----
    bool _is_small() const { return _mpz == nullptr; }

    static mpz_ptr _mpz_new() {
        mpz_ptr p = new __mpz_struct;
        mpz_init(p);
        return p;
    }
    static void _mpz_del(mpz_ptr p) {
        mpz_clear(p);
        delete p;
    }

    void _promote() {
        _mpz = _mpz_new();
        mpz_set_si(_mpz, _val);
    }

    // Returns true if mpz value fits int64_t
    static bool _fits_si(mpz_srcptr z) {
        // mpz_fits_slong_p checks if it fits in a long.
        // On LP64 (Linux), long == int64_t so this is fine.
        // On LLP64 (Windows), long is 32-bit, so we do a manual check.
        if (sizeof(long) >= 8) {
            return mpz_fits_slong_p(z) != 0;
        }
        // Manual: check bit size <= 63
        size_t bits = mpz_sizeinbase(z, 2);
        if (bits < 63) return true;
        if (bits > 63) return false;
        // bits == 63: could be INT64_MIN .. INT64_MAX
        // Compare against limits
        if (mpz_sgn(z) >= 0) {
            // check z <= INT64_MAX
            static mpz_t max_val;
            static bool init = false;
            if (!init) {
                mpz_init(max_val);
                mpz_set_si(max_val, INT64_MAX);
                init = true;
            }
            return mpz_cmp(z, max_val) <= 0;
        } else {
            // check z >= INT64_MIN
            static mpz_t min_val;
            static bool init = false;
            if (!init) {
                mpz_init(min_val);
                mpz_set_si(min_val, INT64_MIN);
                // On LLP64 where long is 32-bit, INT64_MIN can't be set by mpz_set_si.
                // Use string fallback.
                if (sizeof(long) < 8) {
                    mpz_set_str(min_val, "-9223372036854775808", 10);
                }
                init = true;
            }
            return mpz_cmp(z, min_val) >= 0;
        }
    }

    void _demote_if_small() {
        if (_mpz && _fits_si(_mpz)) {
            _val = mpz_get_si(_mpz);
            _mpz_del(_mpz);
            _mpz = nullptr;
        }
    }

    mpz_ptr _ensure_mpz() {
        if (_is_small()) _promote();
        return _mpz;
    }

public:
    // ---- constructors ----
    ZZ() : _val(0), _mpz(nullptr) {}
    ZZ(int v) : _val(v), _mpz(nullptr) {}
    ZZ(unsigned int v) : _val(static_cast<int64_t>(v)), _mpz(nullptr) {}
    ZZ(long v) : _val(v), _mpz(nullptr) {}
    ZZ(long long v) : _val(v), _mpz(nullptr) {}

    ZZ(unsigned long v) : _mpz(nullptr) {
        if (v <= static_cast<unsigned long>(INT64_MAX)) {
            _val = static_cast<int64_t>(v);
        } else {
            _val = 0;
            _mpz = _mpz_new();
            mpz_set_ui(_mpz, v);
        }
    }

    ZZ(unsigned long long v) : _mpz(nullptr) {
        if (v <= static_cast<unsigned long long>(INT64_MAX)) {
            _val = static_cast<int64_t>(v);
        } else {
            _val = 0;
            _mpz = _mpz_new();
            if (sizeof(unsigned long) >= 8) {
                mpz_set_ui(_mpz, static_cast<unsigned long>(v));
            } else {
                // LLP64: set high and low parts
                mpz_set_ui(_mpz, static_cast<unsigned long>(v >> 32));
                mpz_mul_2exp(_mpz, _mpz, 32);
                mpz_add_ui(_mpz, _mpz, static_cast<unsigned long>(v & 0xFFFFFFFF));
            }
        }
    }

    ZZ(const char* str, int base = 10) : _val(0), _mpz(nullptr) {
        mpz_ptr tmp = _mpz_new();
        if (mpz_set_str(tmp, str, base) != 0) {
            _mpz_del(tmp);
            throw std::invalid_argument(std::string("ZZ: invalid string: ") + str);
        }
        if (_fits_si(tmp)) {
            _val = mpz_get_si(tmp);
            _mpz_del(tmp);
        } else {
            _mpz = tmp;
        }
    }

    ZZ(const std::string& str, int base = 10) : ZZ(str.c_str(), base) {}

    // copy
    ZZ(const ZZ& o) : _val(o._val), _mpz(nullptr) {
        if (o._mpz) {
            _mpz = _mpz_new();
            mpz_set(_mpz, o._mpz);
        }
    }

    // move
    ZZ(ZZ&& o) noexcept : _val(o._val), _mpz(o._mpz) {
        o._val = 0;
        o._mpz = nullptr;
    }

    ~ZZ() {
        if (_mpz) _mpz_del(_mpz);
    }

    // ---- assignment ----
    ZZ& operator=(const ZZ& o) {
        if (this == &o) return *this;
        if (o._is_small()) {
            if (_mpz) { _mpz_del(_mpz); _mpz = nullptr; }
            _val = o._val;
        } else {
            if (_mpz) {
                mpz_set(_mpz, o._mpz);
            } else {
                _mpz = _mpz_new();
                mpz_set(_mpz, o._mpz);
            }
        }
        return *this;
    }

    ZZ& operator=(ZZ&& o) noexcept {
        if (this == &o) return *this;
        if (_mpz) _mpz_del(_mpz);
        _val = o._val;
        _mpz = o._mpz;
        o._val = 0;
        o._mpz = nullptr;
        return *this;
    }

    ZZ& operator=(long long v) {
        if (_mpz) { _mpz_del(_mpz); _mpz = nullptr; }
        _val = v;
        return *this;
    }

    // ---- swap ----
    void swap(ZZ& o) noexcept {
        std::swap(_val, o._val);
        std::swap(_mpz, o._mpz);
    }

    // ---- comparison ----
    friend bool operator==(const ZZ& a, const ZZ& b) {
        if (a._is_small() && b._is_small()) return a._val == b._val;
        if (a._is_small()) {
            // a small, b big
            return mpz_cmp_si(b._mpz, static_cast<long>(a._val)) == 0
                   || (sizeof(long) < 8 && _fits_si(b._mpz) && mpz_get_si(b._mpz) == a._val);
        }
        if (b._is_small()) {
            return mpz_cmp_si(a._mpz, static_cast<long>(b._val)) == 0
                   || (sizeof(long) < 8 && _fits_si(a._mpz) && mpz_get_si(a._mpz) == b._val);
        }
        return mpz_cmp(a._mpz, b._mpz) == 0;
    }
    friend bool operator!=(const ZZ& a, const ZZ& b) { return !(a == b); }

    friend bool operator<(const ZZ& a, const ZZ& b) {
        if (a._is_small() && b._is_small()) return a._val < b._val;
        if (a._is_small() && !b._is_small()) {
            // a < b  iff  cmp(b, a) > 0
            return mpz_cmp_si(b._mpz, a._val) > 0;
        }
        if (!a._is_small() && b._is_small()) {
            return mpz_cmp_si(a._mpz, b._val) < 0;
        }
        return mpz_cmp(a._mpz, b._mpz) < 0;
    }
    friend bool operator>(const ZZ& a, const ZZ& b)  { return b < a; }
    friend bool operator<=(const ZZ& a, const ZZ& b) { return !(b < a); }
    friend bool operator>=(const ZZ& a, const ZZ& b) { return !(a < b); }

    // ---- state query ----
    explicit operator bool() const {
        if (_is_small()) return _val != 0;
        return mpz_sgn(_mpz) != 0;
    }
    bool is_zero() const {
        if (_is_small()) return _val == 0;
        return mpz_sgn(_mpz) == 0;
    }
    bool is_one() const {
        if (_is_small()) return _val == 1;
        return mpz_cmp_si(_mpz, 1) == 0;
    }
    bool is_odd() const {
        if (_is_small()) return (_val & 1) != 0;
        return mpz_odd_p(_mpz) != 0;
    }
    int sgn() const {
        if (_is_small()) return (_val > 0) ? 1 : ((_val < 0) ? -1 : 0);
        return mpz_sgn(_mpz);
    }

    // ---- unary minus ----
    friend ZZ operator-(const ZZ& a) {
        if (a._is_small()) {
            if (a._val == INT64_MIN) {
                // overflow: promote
                ZZ r;
                r._mpz = _mpz_new();
                mpz_set_si(r._mpz, INT64_MIN);
                mpz_neg(r._mpz, r._mpz);
                return r;
            }
            return ZZ(-a._val);
        }
        ZZ r;
        r._mpz = _mpz_new();
        mpz_neg(r._mpz, a._mpz);
        r._demote_if_small();
        return r;
    }

    // ---- arithmetic: addition ----
    friend ZZ operator+(const ZZ& a, const ZZ& b) {
        if (a._is_small() && b._is_small()) {
            int64_t r;
            if (!__builtin_add_overflow(a._val, b._val, &r))
                return ZZ(r);
        }
        ZZ res;
        res._mpz = _mpz_new();
        if (a._is_small() && b._is_small()) {
            // overflow case
            mpz_set_si(res._mpz, a._val);
            if (b._val >= 0)
                mpz_add_ui(res._mpz, res._mpz, static_cast<unsigned long>(b._val));
            else
                mpz_sub_ui(res._mpz, res._mpz, static_cast<unsigned long>(-b._val));
        } else if (a._is_small()) {
            if (a._val >= 0)
                mpz_add_ui(res._mpz, b._mpz, static_cast<unsigned long>(a._val));
            else {
                mpz_sub_ui(res._mpz, b._mpz, static_cast<unsigned long>(-a._val));
            }
        } else if (b._is_small()) {
            if (b._val >= 0)
                mpz_add_ui(res._mpz, a._mpz, static_cast<unsigned long>(b._val));
            else
                mpz_sub_ui(res._mpz, a._mpz, static_cast<unsigned long>(-b._val));
        } else {
            mpz_add(res._mpz, a._mpz, b._mpz);
        }
        res._demote_if_small();
        return res;
    }

    // ---- arithmetic: subtraction ----
    friend ZZ operator-(const ZZ& a, const ZZ& b) {
        if (a._is_small() && b._is_small()) {
            int64_t r;
            if (!__builtin_sub_overflow(a._val, b._val, &r))
                return ZZ(r);
        }
        ZZ res;
        res._mpz = _mpz_new();
        if (a._is_small() && b._is_small()) {
            mpz_set_si(res._mpz, a._val);
            if (b._val >= 0)
                mpz_sub_ui(res._mpz, res._mpz, static_cast<unsigned long>(b._val));
            else
                mpz_add_ui(res._mpz, res._mpz, static_cast<unsigned long>(-b._val));
        } else if (a._is_small()) {
            mpz_set_si(res._mpz, a._val);
            mpz_sub(res._mpz, res._mpz, b._mpz);
        } else if (b._is_small()) {
            if (b._val >= 0)
                mpz_sub_ui(res._mpz, a._mpz, static_cast<unsigned long>(b._val));
            else
                mpz_add_ui(res._mpz, a._mpz, static_cast<unsigned long>(-b._val));
        } else {
            mpz_sub(res._mpz, a._mpz, b._mpz);
        }
        res._demote_if_small();
        return res;
    }

    // ---- arithmetic: multiplication ----
    friend ZZ operator*(const ZZ& a, const ZZ& b) {
        if (a._is_small() && b._is_small()) {
#ifdef __SIZEOF_INT128__
            __int128 r = static_cast<__int128>(a._val) * static_cast<__int128>(b._val);
            if (r >= INT64_MIN && r <= INT64_MAX)
                return ZZ(static_cast<long long>(r));
#else
            int64_t r;
            if (!__builtin_mul_overflow(a._val, b._val, &r))
                return ZZ(r);
#endif
        }
        ZZ res;
        res._mpz = _mpz_new();
        if (a._is_small() && b._is_small()) {
            mpz_set_si(res._mpz, a._val);
            mpz_mul_si(res._mpz, res._mpz, b._val);
        } else if (a._is_small()) {
            mpz_mul_si(res._mpz, b._mpz, a._val);
        } else if (b._is_small()) {
            mpz_mul_si(res._mpz, a._mpz, b._val);
        } else {
            mpz_mul(res._mpz, a._mpz, b._mpz);
        }
        res._demote_if_small();
        return res;
    }

    // ---- arithmetic: truncated division (consistent with mpz_class) ----
    friend ZZ operator/(const ZZ& a, const ZZ& b) {
        if (a._is_small() && b._is_small()) {
            if (b._val == 0) throw std::domain_error("ZZ: division by zero");
            // handle INT64_MIN / -1 overflow
            if (a._val == INT64_MIN && b._val == -1) {
                ZZ r;
                r._mpz = _mpz_new();
                mpz_set_si(r._mpz, INT64_MIN);
                mpz_neg(r._mpz, r._mpz);
                return r;
            }
            return ZZ(a._val / b._val);
        }
        ZZ res;
        res._mpz = _mpz_new();
        if (a._is_small()) {
            mpz_t tmp;
            mpz_init_set_si(tmp, a._val);
            mpz_tdiv_q(res._mpz, tmp, b._mpz);
            mpz_clear(tmp);
        } else if (b._is_small()) {
            mpz_t tmp;
            mpz_init_set_si(tmp, b._val);
            mpz_tdiv_q(res._mpz, a._mpz, tmp);
            mpz_clear(tmp);
        } else {
            mpz_tdiv_q(res._mpz, a._mpz, b._mpz);
        }
        res._demote_if_small();
        return res;
    }

    // ---- arithmetic: modulo (truncated, consistent with mpz_class) ----
    friend ZZ operator%(const ZZ& a, const ZZ& b) {
        if (a._is_small() && b._is_small()) {
            if (b._val == 0) throw std::domain_error("ZZ: modulo by zero");
            if (a._val == INT64_MIN && b._val == -1) return ZZ(0);
            return ZZ(a._val % b._val);
        }
        ZZ res;
        res._mpz = _mpz_new();
        if (a._is_small()) {
            mpz_t tmp;
            mpz_init_set_si(tmp, a._val);
            mpz_tdiv_r(res._mpz, tmp, b._mpz);
            mpz_clear(tmp);
        } else if (b._is_small()) {
            mpz_t tmp;
            mpz_init_set_si(tmp, b._val);
            mpz_tdiv_r(res._mpz, a._mpz, tmp);
            mpz_clear(tmp);
        } else {
            mpz_tdiv_r(res._mpz, a._mpz, b._mpz);
        }
        res._demote_if_small();
        return res;
    }

    // ---- compound assignment ----
    ZZ& operator+=(const ZZ& b) {
        if (_is_small() && b._is_small()) {
            int64_t r;
            if (!__builtin_add_overflow(_val, b._val, &r)) {
                _val = r;
                return *this;
            }
        }
        _ensure_mpz();
        if (b._is_small()) {
            if (b._val >= 0)
                mpz_add_ui(_mpz, _mpz, static_cast<unsigned long>(b._val));
            else
                mpz_sub_ui(_mpz, _mpz, static_cast<unsigned long>(-b._val));
        } else {
            mpz_add(_mpz, _mpz, b._mpz);
        }
        _demote_if_small();
        return *this;
    }

    ZZ& operator-=(const ZZ& b) {
        if (_is_small() && b._is_small()) {
            int64_t r;
            if (!__builtin_sub_overflow(_val, b._val, &r)) {
                _val = r;
                return *this;
            }
        }
        _ensure_mpz();
        if (b._is_small()) {
            if (b._val >= 0)
                mpz_sub_ui(_mpz, _mpz, static_cast<unsigned long>(b._val));
            else
                mpz_add_ui(_mpz, _mpz, static_cast<unsigned long>(-b._val));
        } else {
            mpz_sub(_mpz, _mpz, b._mpz);
        }
        _demote_if_small();
        return *this;
    }

    ZZ& operator*=(const ZZ& b) {
        if (_is_small() && b._is_small()) {
#ifdef __SIZEOF_INT128__
            __int128 r = static_cast<__int128>(_val) * static_cast<__int128>(b._val);
            if (r >= INT64_MIN && r <= INT64_MAX) {
                _val = static_cast<int64_t>(r);
                return *this;
            }
#else
            int64_t r;
            if (!__builtin_mul_overflow(_val, b._val, &r)) {
                _val = r;
                return *this;
            }
#endif
        }
        _ensure_mpz();
        if (b._is_small()) {
            mpz_mul_si(_mpz, _mpz, b._val);
        } else {
            mpz_mul(_mpz, _mpz, b._mpz);
        }
        _demote_if_small();
        return *this;
    }

    ZZ& operator/=(const ZZ& b) { *this = *this / b; return *this; }
    ZZ& operator%=(const ZZ& b) { *this = *this % b; return *this; }
    ZZ& operator>>=(unsigned long n) {
        if (_is_small()) {
            if (n >= 63) { _val = (_val < 0) ? -1 : 0; }
            else { _val >>= n; }
        } else {
            mpz_fdiv_q_2exp(_mpz, _mpz, n);
            _demote_if_small();
        }
        return *this;
    }
    ZZ& operator<<=(unsigned long n) {
        _ensure_mpz();
        mpz_mul_2exp(_mpz, _mpz, n);
        _demote_if_small();
        return *this;
    }
    friend ZZ operator>>(ZZ&& a,       unsigned long n) { a >>= n; return std::move(a); }
    friend ZZ operator>>(const ZZ& a,  unsigned long n) { ZZ r(a); r >>= n; return r; }
    friend ZZ operator<<(ZZ&& a,       unsigned long n) { a <<= n; return std::move(a); }
    friend ZZ operator<<(const ZZ& a,  unsigned long n) { ZZ r(a); r <<= n; return r; }

    // ---- special: addmul / submul ----
    void addmul(const ZZ& a, const ZZ& b) {
        if (_is_small() && a._is_small() && b._is_small()) {
#ifdef __SIZEOF_INT128__
            __int128 prod = static_cast<__int128>(a._val) * static_cast<__int128>(b._val);
            __int128 sum = static_cast<__int128>(_val) + prod;
            if (sum >= INT64_MIN && sum <= INT64_MAX) {
                _val = static_cast<int64_t>(sum);
                return;
            }
#endif
        }
        _ensure_mpz();
        if (a._is_small() && b._is_small()) {
            // promote both to mpz temporaries for mpz_addmul
            mpz_t ta, tb;
            mpz_init_set_si(ta, a._val);
            mpz_init_set_si(tb, b._val);
            mpz_addmul(_mpz, ta, tb);
            mpz_clear(ta);
            mpz_clear(tb);
        } else if (a._is_small()) {
            mpz_t ta;
            mpz_init_set_si(ta, a._val);
            mpz_addmul(_mpz, ta, b._mpz);
            mpz_clear(ta);
        } else if (b._is_small()) {
            mpz_t tb;
            mpz_init_set_si(tb, b._val);
            mpz_addmul(_mpz, a._mpz, tb);
            mpz_clear(tb);
        } else {
            mpz_addmul(_mpz, a._mpz, b._mpz);
        }
        _demote_if_small();
    }

    void submul(const ZZ& a, const ZZ& b) {
        if (_is_small() && a._is_small() && b._is_small()) {
#ifdef __SIZEOF_INT128__
            __int128 prod = static_cast<__int128>(a._val) * static_cast<__int128>(b._val);
            __int128 diff = static_cast<__int128>(_val) - prod;
            if (diff >= INT64_MIN && diff <= INT64_MAX) {
                _val = static_cast<int64_t>(diff);
                return;
            }
#endif
        }
        _ensure_mpz();
        if (a._is_small() && b._is_small()) {
            mpz_t ta, tb;
            mpz_init_set_si(ta, a._val);
            mpz_init_set_si(tb, b._val);
            mpz_submul(_mpz, ta, tb);
            mpz_clear(ta);
            mpz_clear(tb);
        } else if (a._is_small()) {
            mpz_t ta;
            mpz_init_set_si(ta, a._val);
            mpz_submul(_mpz, ta, b._mpz);
            mpz_clear(ta);
        } else if (b._is_small()) {
            mpz_t tb;
            mpz_init_set_si(tb, b._val);
            mpz_submul(_mpz, a._mpz, tb);
            mpz_clear(tb);
        } else {
            mpz_submul(_mpz, a._mpz, b._mpz);
        }
        _demote_if_small();
    }

    // ---- fdiv_ui: floor-mod for Zp ----
    uint64_t fdiv_ui(uint64_t d) const {
        if (_is_small()) {
            if (_val >= 0) return static_cast<uint64_t>(_val) % d;
            // negative: floor mod = d - ((-val) % d), but 0 if divides evenly
            uint64_t r = static_cast<uint64_t>(-_val) % d;
            return r == 0 ? 0 : d - r;
        }
        return mpz_fdiv_ui(_mpz, d);
    }

    // ---- sizeinbase ----
    size_t sizeinbase(int base) const {
        if (_is_small()) {
            if (_val == 0) return 1;
            int64_t v = _val < 0 ? -_val : _val;
            // Special case: _val == INT64_MIN, -_val overflows
            if (_val == INT64_MIN) {
                mpz_t tmp;
                mpz_init_set_si(tmp, INT64_MIN);
                mpz_neg(tmp, tmp);
                size_t r = mpz_sizeinbase(tmp, base);
                mpz_clear(tmp);
                return r;
            }
            size_t count = 0;
            while (v > 0) { v /= base; ++count; }
            return count;
        }
        return mpz_sizeinbase(_mpz, base);
    }

    // ---- get_si: extract as long long ----
    long long get_si() const {
        if (_is_small()) return _val;
        return static_cast<long long>(mpz_get_si(_mpz));
    }

    // ---- floor division (static, for number.hh template specialization) ----
    static void fdiv_q(ZZ& q, const ZZ& a, const ZZ& b) {
        if (a._is_small() && b._is_small()) {
            if (b._val == 0) throw std::domain_error("ZZ: division by zero");
            if (a._val == INT64_MIN && b._val == -1) {
                q._ensure_mpz();
                mpz_set_si(q._mpz, INT64_MIN);
                mpz_neg(q._mpz, q._mpz);
                return;
            }
            // floor division for signed integers
            int64_t quot = a._val / b._val;
            int64_t rem  = a._val % b._val;
            // Adjust: if remainder != 0 and signs differ, floor towards -inf
            if (rem != 0 && ((rem ^ b._val) < 0))
                --quot;
            if (q._mpz) { _mpz_del(q._mpz); q._mpz = nullptr; }
            q._val = quot;
            return;
        }
        q._ensure_mpz();
        if (a._is_small()) {
            mpz_t ta;
            mpz_init_set_si(ta, a._val);
            mpz_fdiv_q(q._mpz, ta, b._mpz);
            mpz_clear(ta);
        } else if (b._is_small()) {
            mpz_t tb;
            mpz_init_set_si(tb, b._val);
            mpz_fdiv_q(q._mpz, a._mpz, tb);
            mpz_clear(tb);
        } else {
            mpz_fdiv_q(q._mpz, a._mpz, b._mpz);
        }
        q._demote_if_small();
    }

    static void fdiv_qr(ZZ& q, ZZ& r, const ZZ& a, const ZZ& b) {
        if (a._is_small() && b._is_small()) {
            if (b._val == 0) throw std::domain_error("ZZ: division by zero");
            if (a._val == INT64_MIN && b._val == -1) {
                q._ensure_mpz();
                mpz_set_si(q._mpz, INT64_MIN);
                mpz_neg(q._mpz, q._mpz);
                if (r._mpz) { _mpz_del(r._mpz); r._mpz = nullptr; }
                r._val = 0;
                return;
            }
            int64_t quot = a._val / b._val;
            int64_t rem  = a._val % b._val;
            if (rem != 0 && ((rem ^ b._val) < 0)) {
                --quot;
                rem += b._val;
            }
            if (q._mpz) { _mpz_del(q._mpz); q._mpz = nullptr; }
            q._val = quot;
            if (r._mpz) { _mpz_del(r._mpz); r._mpz = nullptr; }
            r._val = rem;
            return;
        }
        q._ensure_mpz();
        r._ensure_mpz();
        if (a._is_small()) {
            mpz_t ta;
            mpz_init_set_si(ta, a._val);
            mpz_fdiv_qr(q._mpz, r._mpz, ta, b._mpz);
            mpz_clear(ta);
        } else if (b._is_small()) {
            mpz_t tb;
            mpz_init_set_si(tb, b._val);
            mpz_fdiv_qr(q._mpz, r._mpz, a._mpz, tb);
            mpz_clear(tb);
        } else {
            mpz_fdiv_qr(q._mpz, r._mpz, a._mpz, b._mpz);
        }
        q._demote_if_small();
        r._demote_if_small();
    }

    // ---- floor remainder (static) ----
    static void fdiv_r(ZZ& r, const ZZ& a, const ZZ& b) {
        if (a._is_small() && b._is_small()) {
            if (b._val == 0) throw std::domain_error("ZZ: division by zero");
            int64_t rem = a._val % b._val;
            if (rem != 0 && ((rem ^ b._val) < 0))
                rem += b._val;
            if (r._mpz) { _mpz_del(r._mpz); r._mpz = nullptr; }
            r._val = rem;
            return;
        }
        r._ensure_mpz();
        if (a._is_small()) {
            mpz_t ta;
            mpz_init_set_si(ta, a._val);
            mpz_fdiv_r(r._mpz, ta, b._mpz);
            mpz_clear(ta);
        } else if (b._is_small()) {
            mpz_t tb;
            mpz_init_set_si(tb, b._val);
            mpz_fdiv_r(r._mpz, a._mpz, tb);
            mpz_clear(tb);
        } else {
            mpz_fdiv_r(r._mpz, a._mpz, b._mpz);
        }
        r._demote_if_small();
    }

    // ---- modular inverse (static) ----
    // Computes result = a^{-1} mod m. Returns true if inverse exists.
    static bool invert(ZZ& result, const ZZ& a, const ZZ& m) {
        result._ensure_mpz();
        mpz_t ta, tm;
        bool need_clear_a = a._is_small(), need_clear_m = m._is_small();
        if (need_clear_a) { mpz_init_set_si(ta, a._val); }
        if (need_clear_m) { mpz_init_set_si(tm, m._val); }
        const mpz_ptr pa = need_clear_a ? ta : a._mpz;
        const mpz_ptr pm = need_clear_m ? tm : m._mpz;
        int ok = mpz_invert(result._mpz, pa, pm);
        if (need_clear_a) mpz_clear(ta);
        if (need_clear_m) mpz_clear(tm);
        result._demote_if_small();
        return ok != 0;
    }

    // ---- IO ----
    friend std::ostream& operator<<(std::ostream& os, const ZZ& z) {
        if (z._is_small()) {
            os << z._val;
        } else {
            char* str = mpz_get_str(nullptr, 10, z._mpz);
            os << str;
            void (*freefunc)(void*, size_t);
            mp_get_memory_functions(nullptr, nullptr, &freefunc);
            freefunc(str, std::strlen(str) + 1);
        }
        return os;
    }

    friend std::istream& operator>>(std::istream& is, ZZ& z) {
        std::string s;
        is >> s;
        if (!s.empty()) {
            z = ZZ(s);
        }
        return is;
    }

    // ---- hash ----
    std::size_t hash() const {
        if (_is_small()) {
            return std::hash<int64_t>{}(_val);
        }
        // hash mpz limbs (compatible with old hash_value for mpz_class)
        std::size_t seed = boost::hash_range(
            _mpz->_mp_d,
            _mpz->_mp_d + std::abs(_mpz->_mp_size));
        boost::hash_combine(seed, _mpz->_mp_size);
        return seed;
    }

    // ---- friend free functions (declared here, defined below) ----
    friend ZZ pow(const ZZ& base, uint64_t exp);
    friend ZZ gcd(const ZZ& a, const ZZ& b);
    friend ZZ lcm(const ZZ& a, const ZZ& b);
    friend ZZ abs(const ZZ& a);
    friend int sgn(const ZZ& a);
};

// ---- free functions ----

inline ZZ pow(const ZZ& base, uint64_t exp) {
    if (exp == 0) return ZZ(1);
    if (exp == 1) return base;
    if (base.is_zero()) return ZZ(0);
    if (base.is_one()) return ZZ(1);

    // For small base and small exponent, try direct computation
    // But for safety, delegate to mpz_pow_ui
    ZZ r;
    r._mpz = ZZ::_mpz_new();
    if (base._is_small()) {
        mpz_set_si(r._mpz, base._val);
    } else {
        mpz_set(r._mpz, base._mpz);
    }
    mpz_pow_ui(r._mpz, r._mpz, static_cast<unsigned long>(exp));
    r._demote_if_small();
    return r;
}

inline ZZ gcd(const ZZ& a, const ZZ& b) {
    if (a._is_small() && b._is_small()) {
        // Simple GCD for small integers
        int64_t x = a._val < 0 ? -a._val : a._val;
        int64_t y = b._val < 0 ? -b._val : b._val;
        // Handle INT64_MIN
        if (a._val == INT64_MIN || b._val == INT64_MIN) {
            // fall through to mpz
        } else {
            while (y) { int64_t t = y; y = x % y; x = t; }
            return ZZ(x);
        }
    }
    ZZ r;
    r._mpz = ZZ::_mpz_new();
    if (a._is_small() && b._is_small()) {
        mpz_t ta, tb;
        mpz_init_set_si(ta, a._val);
        mpz_init_set_si(tb, b._val);
        mpz_gcd(r._mpz, ta, tb);
        mpz_clear(ta); mpz_clear(tb);
    } else if (a._is_small()) {
        mpz_t ta;
        mpz_init_set_si(ta, a._val);
        mpz_gcd(r._mpz, ta, b._mpz);
        mpz_clear(ta);
    } else if (b._is_small()) {
        mpz_t tb;
        mpz_init_set_si(tb, b._val);
        mpz_gcd(r._mpz, a._mpz, tb);
        mpz_clear(tb);
    } else {
        mpz_gcd(r._mpz, a._mpz, b._mpz);
    }
    r._demote_if_small();
    return r;
}

inline ZZ lcm(const ZZ& a, const ZZ& b) {
    if (a.is_zero() || b.is_zero()) return ZZ(0);
    if (a._is_small() && b._is_small() && a._val != INT64_MIN && b._val != INT64_MIN) {
        ZZ g = gcd(a, b);
        // |a| / g * |b| to avoid overflow
        int64_t av = a._val < 0 ? -a._val : a._val;
        int64_t bv = b._val < 0 ? -b._val : b._val;
        int64_t gv = g._val;
        int64_t div_result = av / gv;
#ifdef __SIZEOF_INT128__
        __int128 r = static_cast<__int128>(div_result) * static_cast<__int128>(bv);
        if (r >= 0 && r <= INT64_MAX)
            return ZZ(static_cast<long long>(r));
#endif
    }
    ZZ r;
    r._mpz = ZZ::_mpz_new();
    if (a._is_small() && b._is_small()) {
        mpz_t ta, tb;
        mpz_init_set_si(ta, a._val);
        mpz_init_set_si(tb, b._val);
        mpz_lcm(r._mpz, ta, tb);
        mpz_clear(ta); mpz_clear(tb);
    } else if (a._is_small()) {
        mpz_t ta;
        mpz_init_set_si(ta, a._val);
        mpz_lcm(r._mpz, ta, b._mpz);
        mpz_clear(ta);
    } else if (b._is_small()) {
        mpz_t tb;
        mpz_init_set_si(tb, b._val);
        mpz_lcm(r._mpz, a._mpz, tb);
        mpz_clear(tb);
    } else {
        mpz_lcm(r._mpz, a._mpz, b._mpz);
    }
    r._demote_if_small();
    return r;
}

inline ZZ abs(const ZZ& a) {
    if (a._is_small()) {
        if (a._val == INT64_MIN) {
            ZZ r;
            r._mpz = ZZ::_mpz_new();
            mpz_set_si(r._mpz, INT64_MIN);
            mpz_abs(r._mpz, r._mpz);
            return r;
        }
        return ZZ(a._val < 0 ? -a._val : a._val);
    }
    ZZ r;
    r._mpz = ZZ::_mpz_new();
    mpz_abs(r._mpz, a._mpz);
    r._demote_if_small();
    return r;
}

inline int sgn(const ZZ& a) {
    return a.sgn();
}

inline size_t sizeinbase(const ZZ& z, int base) {
    return z.sizeinbase(base);
}

inline std::size_t hash_value(const ZZ& v) {
    return v.hash();
}

inline void swap(ZZ& a, ZZ& b) noexcept {
    a.swap(b);
}

} // namespace clpoly

namespace std {
    template<> struct hash<clpoly::ZZ> {
        std::size_t operator()(const clpoly::ZZ& v) const { return v.hash(); }
    };
}

#endif // CLPOLY_NUMBER_ZZ_HH
