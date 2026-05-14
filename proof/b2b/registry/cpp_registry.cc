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

    throw std::runtime_error("unknown fn: " + fn);
}

}  // namespace b2b
