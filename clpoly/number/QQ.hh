/**
 * @file QQ.hh
 * @brief Custom rational number type built on ZZ.
 *
 * Stores numerator (_num) and denominator (_den) as ZZ values.
 * Denominator is always > 0. Auto-canonicalized after construction and arithmetic.
 */
#ifndef CLPOLY_NUMBER_QQ_HH
#define CLPOLY_NUMBER_QQ_HH

#include <clpoly/number/ZZ.hh>

namespace clpoly {

class QQ {
private:
    ZZ _num;
    ZZ _den; // always > 0

    void _canonicalize() {
        if (_num.is_zero()) {
            _den = ZZ(1);
            return;
        }
        if (_den.sgn() < 0) {
            _num = -_num;
            _den = -_den;
        }
        ZZ g = gcd(abs(_num), _den);
        if (!g.is_one()) {
            _num /= g;
            _den /= g;
        }
    }

public:
    // ---- constructors ----
    QQ() : _num(0), _den(1) {}
    QQ(int v) : _num(v), _den(1) {}
    QQ(long v) : _num(v), _den(1) {}
    QQ(long long v) : _num(v), _den(1) {}
    QQ(const ZZ& v) : _num(v), _den(1) {}

    QQ(const ZZ& num, const ZZ& den) : _num(num), _den(den) {
        if (_den.is_zero()) throw std::domain_error("QQ: zero denominator");
        _canonicalize();
    }

    QQ(long long num, long long den) : _num(num), _den(den) {
        if (_den.is_zero()) throw std::domain_error("QQ: zero denominator");
        _canonicalize();
    }

    QQ(const QQ& o) : _num(o._num), _den(o._den) {}
    QQ(QQ&& o) noexcept : _num(std::move(o._num)), _den(std::move(o._den)) {}
    ~QQ() = default;

    // ---- assignment ----
    QQ& operator=(const QQ& o) {
        if (this != &o) { _num = o._num; _den = o._den; }
        return *this;
    }
    QQ& operator=(QQ&& o) noexcept {
        if (this != &o) { _num = std::move(o._num); _den = std::move(o._den); }
        return *this;
    }
    QQ& operator=(long long v) {
        _num = v; _den = ZZ(1);
        return *this;
    }

    // ---- accessors (compatible with mpq_class API) ----
    const ZZ& num() const { return _num; }
    const ZZ& den() const { return _den; }
    const ZZ& get_num() const { return _num; }
    const ZZ& get_den() const { return _den; }

    void canonicalize() { _canonicalize(); }

    // ---- swap ----
    void swap(QQ& o) noexcept {
        _num.swap(o._num);
        _den.swap(o._den);
    }

    // ---- comparison ----
    friend bool operator==(const QQ& a, const QQ& b) {
        // Both canonicalized, so direct comparison works
        return a._num == b._num && a._den == b._den;
    }
    friend bool operator!=(const QQ& a, const QQ& b) { return !(a == b); }

    friend bool operator<(const QQ& a, const QQ& b) {
        // a/ad < b/bd  iff  a*bd < b*ad  (both denominators positive)
        return a._num * b._den < b._num * a._den;
    }
    friend bool operator>(const QQ& a, const QQ& b)  { return b < a; }
    friend bool operator<=(const QQ& a, const QQ& b) { return !(b < a); }
    friend bool operator>=(const QQ& a, const QQ& b) { return !(a < b); }

    // ---- state query ----
    explicit operator bool() const { return !_num.is_zero(); }
    bool is_zero() const { return _num.is_zero(); }

    // ---- unary minus ----
    friend QQ operator-(const QQ& a) {
        QQ r;
        r._num = -a._num;
        r._den = a._den;
        return r;
    }

    // ---- arithmetic ----
    friend QQ operator+(const QQ& a, const QQ& b) {
        // a/ad + b/bd = (a*bd + b*ad) / (ad*bd), then canonicalize
        QQ r;
        r._num = a._num * b._den + b._num * a._den;
        r._den = a._den * b._den;
        r._canonicalize();
        return r;
    }

    friend QQ operator-(const QQ& a, const QQ& b) {
        QQ r;
        r._num = a._num * b._den - b._num * a._den;
        r._den = a._den * b._den;
        r._canonicalize();
        return r;
    }

    friend QQ operator*(const QQ& a, const QQ& b) {
        QQ r;
        r._num = a._num * b._num;
        r._den = a._den * b._den;
        r._canonicalize();
        return r;
    }

    friend QQ operator/(const QQ& a, const QQ& b) {
        if (b._num.is_zero()) throw std::domain_error("QQ: division by zero");
        QQ r;
        r._num = a._num * b._den;
        r._den = a._den * b._num;
        r._canonicalize();
        return r;
    }

    // ---- compound assignment ----
    QQ& operator+=(const QQ& b) { *this = *this + b; return *this; }
    QQ& operator-=(const QQ& b) { *this = *this - b; return *this; }
    QQ& operator*=(const QQ& b) { *this = *this * b; return *this; }
    QQ& operator/=(const QQ& b) { *this = *this / b; return *this; }

    // ---- IO ----
    friend std::ostream& operator<<(std::ostream& os, const QQ& q) {
        os << q._num;
        if (!q._den.is_one()) os << "/" << q._den;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, QQ& q) {
        std::string s;
        is >> s;
        auto pos = s.find('/');
        if (pos == std::string::npos) {
            q._num = ZZ(s);
            q._den = ZZ(1);
        } else {
            q._num = ZZ(s.substr(0, pos));
            q._den = ZZ(s.substr(pos + 1));
            q._canonicalize();
        }
        return is;
    }

    // ---- hash ----
    std::size_t hash() const {
        std::size_t seed = 0;
        boost::hash_combine(seed, _num.hash());
        boost::hash_combine(seed, _den.hash());
        return seed;
    }

    // ---- friend free functions ----
    friend QQ pow(const QQ& base, uint64_t exp);
    friend int sgn(const QQ& q);
    friend QQ abs(const QQ& q);
};

// ---- free functions ----

inline QQ pow(const QQ& base, uint64_t exp) {
    QQ r;
    r._num = pow(base._num, exp);
    r._den = pow(base._den, exp);
    // pow preserves canonicalization if base was canonical
    return r;
}

inline int sgn(const QQ& q) {
    return q._num.sgn();
}

inline QQ abs(const QQ& q) {
    QQ r;
    r._num = abs(q._num);
    r._den = q._den;
    return r;
}

inline std::size_t hash_value(const QQ& v) {
    return v.hash();
}

inline void swap(QQ& a, QQ& b) noexcept {
    a.swap(b);
}

} // namespace clpoly

namespace std {
    template<> struct hash<clpoly::QQ> {
        std::size_t operator()(const clpoly::QQ& v) const { return v.hash(); }
    };
}

#endif // CLPOLY_NUMBER_QQ_HH
