// 强制实例化全部因式分解模板函数。
// 用途：clang++ -Xclang -ast-dump=json -fsyntax-only -std=c++17 -I../../ instantiate.cc
//
// 注意 (2026-05-19): clang ≥ 18.1.8 在 `-fsyntax-only` 下若 result 未被
// ODR-used，会跳过沿调用链的内部模板实例化。因此每个 result 都必须真正
// 被读取（这里通过 .factors.size() + volatile sink 实现），否则
// __upoly_divmod / __upoly_mod 等内部函数会退化为 primary template 声明
// （参数 qualType 变成 `int &`，丢失 body）。

#include <cstddef>
#include <clpoly/polynomial_factorize.hh>

namespace { volatile std::size_t __force_instantiate_sink = 0; }

std::size_t force_instantiate() {
    using namespace clpoly;

    // grlex 是 CLPoly 默认的单项式序（多变量）
    polynomial_<ZZ, grlex> f_mv;
    auto r1 = factorize(f_mv);

    // lex 特化
    polynomial_<ZZ, lex> f_lex;
    auto r2 = factorize(f_lex);

    // 单变量（upolynomial = polynomial_<ZZ, lex> 单变量特化）
    upolynomial_<ZZ> f_uni;
    auto r3 = factorize(f_uni);

    // Zp 模块（__factor_Zp 返回 std::pair<Zp, vector<pair<upoly, u64>>>，
    // 不是 factorization<>，所以用 r4.second 而非 r4.factors）
    upolynomial_<Zp> f_zp;
    auto r4 = __factor_Zp(f_zp);

    // 真正使用 r1..r4 的成员，强制 clang 沿调用链 ODR-instantiate
    std::size_t total = r1.factors.size() + r2.factors.size()
                      + r3.factors.size() + r4.second.size();
    __force_instantiate_sink = total;
    return total;
}
