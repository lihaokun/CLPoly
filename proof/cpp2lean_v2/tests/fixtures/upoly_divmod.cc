// fixture: __upoly_divmod
// 覆盖：void 函数、non-const ref 参数（输出参数）
// 源：clpoly/polynomial_factorize_zp.hh:49-56

#include <clpoly/polynomial_factorize.hh>

namespace clpoly {
    // §4.1.3 商和余式: f = q*g + r in Zp[x]
    inline void __upoly_divmod(
        upolynomial_<Zp>& q,
        upolynomial_<Zp>& r,
        const upolynomial_<Zp>& f,
        const upolynomial_<Zp>& g)
    {
        pair_vec_div(q.data(), r.data(), f.data(), g.data(), f.comp());
    }
}

// 强制实例化
void _instantiate_upoly_divmod() {
    clpoly::upolynomial_<clpoly::Zp> q, r, f, g;
    clpoly::__upoly_divmod(q, r, f, g);
}
