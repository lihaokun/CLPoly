// B2B C++ 端类型库：JSON ↔ C++ 类型双向转换
//
// 协议：每个值在 JSON 中表示为 {"type": "<类型名>", "val": <载荷>}
//   parse_X 收 json 返回 C++ 值（类型不匹配抛 std::runtime_error）
//   serialize_X 收 C++ 值返回 json
//
// Step 1：实现 6 个核心类型：Int32 / Int64 / UInt64 / ZZ / Zp / SparsePolyZp

#ifndef CLPOLY_B2B_TYPES_HH
#define CLPOLY_B2B_TYPES_HH

#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include "nlohmann/json.hpp"
#include "clpoly/number.hh"
#include "clpoly/polynomial.hh"
#include "clpoly/upolynomial.hh"

namespace b2b {

using nlohmann::json;
using clpoly::ZZ;
using clpoly::Zp;
using clpoly::umonomial;
template<class T> using upolynomial_ = clpoly::upolynomial_<T>;

// ---- 内部工具 ----
inline void _check_type(const json& j, const char* expected) {
    if (!j.is_object()) {
        throw std::runtime_error("expected object, got: " + j.dump());
    }
    auto it = j.find("type");
    if (it == j.end() || !it->is_string()) {
        throw std::runtime_error("missing or non-string 'type' field");
    }
    if (it->get<std::string>() != expected) {
        throw std::runtime_error(std::string("type mismatch: expected ") + expected +
                                 ", got " + it->get<std::string>());
    }
}

// ---- 标量 ----

int32_t  parse_Int32(const json& j);
int64_t  parse_Int64(const json& j);
uint64_t parse_UInt64(const json& j);

json serialize_Int32(int32_t v);
json serialize_Int64(int64_t v);
json serialize_UInt64(uint64_t v);

// ---- 大整数 ----
// JSON: {"type":"ZZ","val":"<十进制字符串>"}

ZZ   parse_ZZ(const json& j);
json serialize_ZZ(const ZZ& z);

// ---- Zp ----
// JSON: {"type":"Zp","val":[<val>,<prime>]}

Zp   parse_Zp(const json& j);
json serialize_Zp(const Zp& z);

// ---- BoolZZ tuple ----
// JSON: {"type":"BoolZZ","val":[<bool>,<zz_string>]}
// 用于 ZZ::invert 返回 (success, inverse)

json serialize_BoolZZ(bool ok, const ZZ& z);

// ---- SparsePolyZp ----
// JSON: {"type":"SparsePolyZp","val":[[<deg>,[<val>,<prime>]], ...]}
// 注意：内部 [val,prime] 不带 type 包装（紧凑表示）

upolynomial_<Zp> parse_SparsePolyZp(const json& j);
json             serialize_SparsePolyZp(const upolynomial_<Zp>& p);

// ---- SparsePolyZZ ----
// JSON: {"type":"SparsePolyZZ","val":[[<deg>,"<zz_decimal>"], ...]}

upolynomial_<ZZ> parse_SparsePolyZZ(const json& j);
json             serialize_SparsePolyZZ(const upolynomial_<ZZ>& p);

// ---- 复合返回类型（仅 serialize；解析端不接受这些）----
// PairSPZp: divmod 返回 (q, r)
// JSON: {"type":"PairSPZp","val":[<spzp>,<spzp>]} 内嵌完整 SparsePolyZp 对象
json serialize_PairSPZp(const upolynomial_<Zp>& q, const upolynomial_<Zp>& r);
// TripleSPZp: gcd_eea 返回 (gcd, s, t)
json serialize_TripleSPZp(
    const upolynomial_<Zp>& g,
    const upolynomial_<Zp>& s,
    const upolynomial_<Zp>& t);

// ---- A5 类型：Variable + Monomial + MvPolyZZ/Zp + VarMap ----
// Variable: clpoly::variable 通过 serial 复原（B2B 测试预声明 x/y/z/w = serial 1-4）
// JSON: {"type":"Variable","val":<uint64>}
using Variable = clpoly::variable;
using Monomial = clpoly::basic_monomial<clpoly::lex>;
using MvPolyZZ = clpoly::polynomial_<ZZ, clpoly::lex>;
using MvPolyZp = clpoly::polynomial_<Zp, clpoly::lex>;

Variable parse_Variable(const json& j);
json     serialize_Variable(const Variable& v);

// Monomial: {"type":"Monomial","val":[[<var_serial>, <exp>], ...]}
Monomial parse_Monomial(const json& j);
json     serialize_Monomial(const Monomial& m);

// MvPolyZZ: {"type":"MvPolyZZ","val":[[<mono_val>, "<zz_str>"], ...]}
// 注：mono_val 是 Monomial 的 val 字段（裸数组，无 type 包装）
MvPolyZZ parse_MvPolyZZ(const json& j);
json     serialize_MvPolyZZ(const MvPolyZZ& p);

// MvPolyZp: {"type":"MvPolyZp","val":[[<mono_val>, [<zp_val>,<zp_prime>]], ...]}
MvPolyZp parse_MvPolyZp(const json& j);
json     serialize_MvPolyZp(const MvPolyZp& p);

// VarMap (for assign2): JSON {"type":"VarMap","val":[[<var_serial>, "<zz_str>"], ...]}
// 解析为 std::map<variable, ZZ>
std::map<Variable, ZZ> parse_VarMapZZ(const json& j);

// ---- A6 类型：ZZTriple（Bezout 返回 (g, s, t)）----
// JSON: {"type":"ZZTriple","val":["<g>","<s>","<t>"]}
json serialize_ZZTriple(const ZZ& g, const ZZ& s, const ZZ& t);

}  // namespace b2b

#endif
