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

    throw std::runtime_error("unknown fn: " + fn);
}

}  // namespace b2b
