// fixture: __symmetric_mod
// 覆盖：if-else、BinaryOperator、CompoundAssignOperator (-=)、ZZ 构造
// 源：clpoly/polynomial_factorize_univar.hh:94-103

#include <clpoly/polynomial_factorize.hh>

namespace clpoly {
    // §4.2.1 标量对称模: 将 a 约化到 (-m/2, m/2]
    inline ZZ __symmetric_mod(const ZZ& a, const ZZ& m)
    {
        ZZ r;
        ZZ::fdiv_r(r, a, m);   // r ∈ [0, m)
        ZZ half;
        ZZ::fdiv_q(half, m, ZZ(2));
        if (r > half)
            r -= m;
        return r;
    }
}

// 强制实例化
void _instantiate_symmetric_mod() {
    clpoly::__symmetric_mod(clpoly::ZZ(0), clpoly::ZZ(1));
}
