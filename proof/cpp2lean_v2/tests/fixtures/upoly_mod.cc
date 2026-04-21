// fixture: __upoly_mod
// 覆盖：多语句、VarDecl+LetStmt、CallExpr、MemberExpr
// 源：clpoly/polynomial_factorize_zp.hh:39-46

#include <clpoly/polynomial_factorize.hh>

namespace clpoly {
    // §4.1.2 多项式取模: f mod g in Zp[x]
    inline upolynomial_<Zp> __upoly_mod(
        const upolynomial_<Zp>& f,
        const upolynomial_<Zp>& g)
    {
        upolynomial_<Zp> q, r;
        pair_vec_div(q.data(), r.data(), f.data(), g.data(), f.comp());
        return r;
    }
}

// 强制实例化
void _instantiate_upoly_mod() {
    clpoly::upolynomial_<clpoly::Zp> f, g;
    auto r = clpoly::__upoly_mod(f, g);
}
