// B2B C++ 函数注册表实现
//
// 每加一个函数：写一个 if (fn == "<name>") { ... return ...; } 分支

#include "cpp_registry.hh"
#include "../types/b2b_types.hh"

#include "clpoly/polynomial_factorize_zp.hh"

namespace b2b {

using nlohmann::json;

json dispatch(const std::string& fn, const json& args) {
    if (!args.is_array()) {
        throw std::runtime_error("args: expected array");
    }

    // === 函数 dispatch ===
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

    throw std::runtime_error("unknown fn: " + fn);
}

}  // namespace b2b
