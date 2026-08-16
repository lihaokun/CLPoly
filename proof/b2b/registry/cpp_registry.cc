// B2B C++ 函数注册表实现
//
// 每加一个函数：写一个 if (fn == "<name>") { ... return ...; } 分支

#include "cpp_registry.hh"
#include "../types/b2b_types.hh"

#include "clpoly/polynomial_factorize_zp.hh"
#include "clpoly/polynomial_factorize_univar.hh"
#include "clpoly/polynomial.hh"
#include <cstring>
#include <gmp.h>
#include <sstream>

namespace b2b {

using nlohmann::json;

json dispatch(const std::string& fn, const json& args) {
    if (!args.is_array()) {
        throw std::runtime_error("args: expected array");
    }

    // === 函数 dispatch ===
    if (fn == "__needs_zassenhaus_safety_net") {
        uint64_t result_count = parse_UInt64(args.at(0));
        uint64_t modular_count = parse_UInt64(args.at(1));
        bool at_full_precision = parse_Bool(args.at(2));
        return serialize_Bool(clpoly::__needs_zassenhaus_safety_net(
            result_count, modular_count, at_full_precision));
    }

    if (fn == "__make_zp") {
        int64_t  val = parse_Int64(args.at(0));
        uint64_t p   = parse_UInt64(args.at(1));
        return serialize_Zp(clpoly::__make_zp(val, p));
    }

    // ---- Zp 算术（B2B wrapper）：operator 重载在 C++ 没有命名函数，B2B 起 ____xx_zp 名 ----
    if (fn == "__add_zp") {
        Zp a = parse_Zp(args.at(0));
        Zp b = parse_Zp(args.at(1));
        return serialize_Zp(a + b);
    }
    if (fn == "__sub_zp") {
        Zp a = parse_Zp(args.at(0));
        Zp b = parse_Zp(args.at(1));
        return serialize_Zp(a - b);
    }
    if (fn == "__mul_zp") {
        Zp a = parse_Zp(args.at(0));
        Zp b = parse_Zp(args.at(1));
        return serialize_Zp(a * b);
    }
    if (fn == "__neg_zp") {
        Zp a = parse_Zp(args.at(0));
        return serialize_Zp(-a);
    }
    if (fn == "__inv_zp") {
        Zp a = parse_Zp(args.at(0));
        return serialize_Zp(a.inv());
    }
    if (fn == "__pow_zp") {
        Zp       a = parse_Zp(args.at(0));
        int64_t  i = parse_Int64(args.at(1));
        return serialize_Zp(clpoly::pow(a, i));
    }

    // ---- ZZ 算术（C++ 端用 gcd/lcm/static fdiv_q/r 入口）----
    if (fn == "__gcd_zz") {
        ZZ a = parse_ZZ(args.at(0));
        ZZ b = parse_ZZ(args.at(1));
        return serialize_ZZ(clpoly::gcd(a, b));
    }
    if (fn == "__lcm_zz") {
        ZZ a = parse_ZZ(args.at(0));
        ZZ b = parse_ZZ(args.at(1));
        return serialize_ZZ(clpoly::lcm(a, b));
    }
    if (fn == "__fdiv_q_zz") {
        ZZ a = parse_ZZ(args.at(0));
        ZZ b = parse_ZZ(args.at(1));
        ZZ q;
        ZZ::fdiv_q(q, a, b);
        return serialize_ZZ(q);
    }
    if (fn == "__fdiv_r_zz") {
        ZZ a = parse_ZZ(args.at(0));
        ZZ b = parse_ZZ(args.at(1));
        ZZ r;
        ZZ::fdiv_r(r, a, b);
        return serialize_ZZ(r);
    }
    if (fn == "__fdiv_ui_zz") {
        ZZ       a = parse_ZZ(args.at(0));
        uint64_t b = parse_UInt64(args.at(1));
        return serialize_UInt64(a.fdiv_ui(b));
    }
    if (fn == "__sizeinbase_zz") {
        ZZ      z    = parse_ZZ(args.at(0));
        int32_t base = parse_Int32(args.at(1));
        return serialize_UInt64(z.sizeinbase(base));
    }
    if (fn == "__invert_zz") {
        ZZ a = parse_ZZ(args.at(0));
        ZZ m = parse_ZZ(args.at(1));
        ZZ result;
        bool ok = ZZ::invert(result, a, m);
        return serialize_BoolZZ(ok, result);
    }

    // ---- SparsePolyZZ ops ----
    if (fn == "__cont_zz_poly") {
        auto p = parse_SparsePolyZZ(args.at(0));
        return serialize_ZZ(clpoly::cont(p));
    }
    if (fn == "__pp_zz_poly") {
        auto p = parse_SparsePolyZZ(args.at(0));
        ZZ c = clpoly::cont(p);
        // pp = p / cont(p)（精确除法）。空 p 或 cont=0：返原值（防 div by zero）
        if (c != ZZ(0)) {
            for (auto& term : p) term.second = term.second / c;
        }
        return serialize_SparsePolyZZ(p);
    }
    if (fn == "__derivative_zz_poly") {
        auto p = parse_SparsePolyZZ(args.at(0));
        return serialize_SparsePolyZZ(clpoly::derivative(p));
    }
    if (fn == "__normalization_zz_poly") {
        auto p = parse_SparsePolyZZ(args.at(0));
        p.normalization();
        return serialize_SparsePolyZZ(p);
    }
    if (fn == "__polynomial_mod_zz") {
        auto     p     = parse_SparsePolyZZ(args.at(0));
        uint64_t prime = parse_UInt64(args.at(1));
        return serialize_SparsePolyZp(clpoly::polynomial_mod(p, prime));
    }

    // ---- SparsePolyZp ops ----
    if (fn == "__add_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        return serialize_SparsePolyZp(f + g);
    }
    if (fn == "__sub_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        return serialize_SparsePolyZp(f - g);
    }
    if (fn == "__mul_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        return serialize_SparsePolyZp(f * g);
    }
    if (fn == "__divmod_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        upolynomial_<Zp> q, r;
        clpoly::__upoly_divmod(q, r, f, g);
        return serialize_PairSPZp(q, r);
    }
    if (fn == "__gcd_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        return serialize_SparsePolyZp(clpoly::polynomial_GCD(f, g));
    }
    if (fn == "__gcd_eea_zp_poly") {
        auto f = parse_SparsePolyZp(args.at(0));
        auto g = parse_SparsePolyZp(args.at(1));
        upolynomial_<Zp> s, t;
        auto gcd = clpoly::polynomial_GCD(f, g, s, t);
        return serialize_TripleSPZp(gcd, s, t);
    }
    if (fn == "__derivative_zp_poly") {
        auto p = parse_SparsePolyZp(args.at(0));
        return serialize_SparsePolyZp(clpoly::derivative(p));
    }
    if (fn == "__normalization_zp_poly") {
        auto p = parse_SparsePolyZp(args.at(0));
        p.normalization();
        return serialize_SparsePolyZp(p);
    }

    // ---- A5: MvPoly ops ----
    if (fn == "__assign_mv_zz") {
        auto p = parse_MvPolyZZ(args.at(0));
        auto v = parse_Variable(args.at(1));
        auto c = parse_ZZ(args.at(2));
        return serialize_MvPolyZZ(clpoly::assign(p, v, c));
    }
    if (fn == "__assign2_mv_zz") {
        auto p = parse_MvPolyZZ(args.at(0));
        auto m = parse_VarMapZZ(args.at(1));
        return serialize_MvPolyZZ(clpoly::assign(p, m));
    }
    if (fn == "__leadcoeff_mv_zz") {
        auto p = parse_MvPolyZZ(args.at(0));
        auto v = parse_Variable(args.at(1));
        return serialize_MvPolyZZ(clpoly::leadcoeff(p, v));
    }
    if (fn == "__poly_convert_spzp_to_spzz") {
        auto p_in = parse_SparsePolyZp(args.at(0));
        upolynomial_<ZZ> p_out;
        clpoly::poly_convert(p_in, p_out);
        return serialize_SparsePolyZZ(p_out);
    }
    if (fn == "__poly_convert3_spzz_to_mvzz") {
        auto p_in = parse_SparsePolyZZ(args.at(0));
        auto v    = parse_Variable(args.at(1));
        MvPolyZZ p_out;
        clpoly::poly_convert(p_in, p_out, v);
        return serialize_MvPolyZZ(p_out);
    }

    // ---- A6: Bezout extGcd ----
    if (fn == "__nat_extgcd") {
        // 输入两个非负整数 a, b，返回 (g, s, t) 满足 s*a + t*b = g (= gcd(a,b))
        // 直接调 GMP mpz_gcdext，然后用字符串往返回 ZZ（绕过私有 mpz 访问）
        ZZ a = parse_ZZ(args.at(0));
        ZZ b = parse_ZZ(args.at(1));
        mpz_t g_mp, s_mp, t_mp, a_mp, b_mp;
        mpz_inits(g_mp, s_mp, t_mp, NULL);
        // 把 ZZ 转 mpz_t（用 ZZ::operator<<）
        std::ostringstream a_os; a_os << a;
        std::ostringstream b_os; b_os << b;
        mpz_init_set_str(a_mp, a_os.str().c_str(), 10);
        mpz_init_set_str(b_mp, b_os.str().c_str(), 10);
        mpz_gcdext(g_mp, s_mp, t_mp, a_mp, b_mp);
        // 从 mpz_t 重建 ZZ
        char* g_str = mpz_get_str(nullptr, 10, g_mp);
        char* s_str = mpz_get_str(nullptr, 10, s_mp);
        char* t_str = mpz_get_str(nullptr, 10, t_mp);
        ZZ g(g_str), s(s_str), t(t_str);
        void (*freefunc)(void*, size_t);
        mp_get_memory_functions(nullptr, nullptr, &freefunc);
        freefunc(g_str, std::strlen(g_str) + 1);
        freefunc(s_str, std::strlen(s_str) + 1);
        freefunc(t_str, std::strlen(t_str) + 1);
        mpz_clears(g_mp, s_mp, t_mp, a_mp, b_mp, NULL);
        return serialize_ZZTriple(g, s, t);
    }

    throw std::runtime_error("unknown fn: " + fn);
}

}  // namespace b2b
