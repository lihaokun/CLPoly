// fixture: __upoly_const_term
// 覆盖：if-return 链、CXXMemberCallExpr (.empty(), .back(), .deg())、
//       MemberExpr、FieldAccess、ZZ literal
// 源：clpoly/polynomial_factorize_univar.hh:743-748

#include <clpoly/polynomial_factorize.hh>

namespace clpoly {
    // 获取 ZZ 多项式的常数项
    inline ZZ __upoly_const_term(const upolynomial_<ZZ>& f)
    {
        if (f.empty()) return ZZ(0);
        if (f.back().first.deg() == 0) return f.back().second;
        return ZZ(0);
    }
}

// 强制实例化
void _instantiate_upoly_const_term() {
    clpoly::upolynomial_<clpoly::ZZ> f;
    clpoly::__upoly_const_term(f);
}
