// 强制实例化全部因式分解模板函数。
// 用途：clang++ -Xclang -ast-dump=json -fsyntax-only -std=c++17 -I../../ instantiate.cc
// 这样 Clang AST 包含所有具体化的函数定义。

#include <clpoly/polynomial_factorize.hh>

void force_instantiate() {
    using namespace clpoly;

    // grlex 是 CLPoly 默认的单项式序（多变量）
    // 调用顶层 factorize 会沿调用链实例化 wang/univar/zp 的全部模板
    polynomial_<ZZ, grlex> f_mv;
    auto r1 = factorize(f_mv);

    // lex 特化
    polynomial_<ZZ, lex> f_lex;
    auto r2 = factorize(f_lex);

    // 单变量（upolynomial = polynomial_<ZZ, lex> 单变量特化）
    upolynomial_<ZZ> f_uni;
    auto r3 = factorize(f_uni);

    // Zp 模块
    upolynomial_<Zp> f_zp;
    auto r4 = __factor_Zp(f_zp);
}
