// fixture: __make_zp (最简单的 CLPoly 函数)
// 覆盖：FunctionDecl, ReturnStmt, CXXConstructExpr, ParmVarDecl
// 源：clpoly/polynomial_factorize_zp.hh:20

#include <clpoly/polynomial_factorize.hh>

namespace clpoly {
    // Helper: 构造 Zp(int64_t, uint64_t) 避免重载歧义
    inline Zp __make_zp(int64_t val, uint64_t p) { return Zp(val, p); }
}

// 强制实例化
void _instantiate_make_zp() {
    clpoly::__make_zp(0, 1);
}
