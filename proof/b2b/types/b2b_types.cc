// B2B C++ 端类型库实现

#include "b2b_types.hh"
#include <sstream>

namespace b2b {

// ---- 标量 ----

int32_t parse_Int32(const json& j) {
    _check_type(j, "Int32");
    return j.at("val").get<int32_t>();
}

int64_t parse_Int64(const json& j) {
    _check_type(j, "Int64");
    return j.at("val").get<int64_t>();
}

uint64_t parse_UInt64(const json& j) {
    _check_type(j, "UInt64");
    // JSON 数字范围限制：UInt64 大值用字符串避免精度丢失
    const auto& v = j.at("val");
    if (v.is_string()) {
        return std::stoull(v.get<std::string>());
    }
    return v.get<uint64_t>();
}

json serialize_Int32(int32_t v)  { return {{"type","Int32"},  {"val", v}}; }
json serialize_Int64(int64_t v)  { return {{"type","Int64"},  {"val", v}}; }
json serialize_UInt64(uint64_t v){ return {{"type","UInt64"}, {"val", v}}; }

// ---- ZZ ----

ZZ parse_ZZ(const json& j) {
    _check_type(j, "ZZ");
    const auto& v = j.at("val");
    if (v.is_string()) {
        return ZZ(v.get<std::string>());
    }
    // 容错：小整数 JSON number 也接受
    if (v.is_number_integer()) {
        return ZZ(v.get<int64_t>());
    }
    throw std::runtime_error("ZZ val: expected string or integer");
}

json serialize_ZZ(const ZZ& z) {
    std::ostringstream os;
    os << z;
    return {{"type","ZZ"}, {"val", os.str()}};
}

// ---- Zp ----

Zp parse_Zp(const json& j) {
    _check_type(j, "Zp");
    const auto& v = j.at("val");
    if (!v.is_array() || v.size() != 2) {
        throw std::runtime_error("Zp val: expected [val, prime]");
    }
    int64_t  val   = v[0].is_string() ? std::stoll(v[0].get<std::string>())
                                      : v[0].get<int64_t>();
    uint64_t prime = v[1].is_string() ? std::stoull(v[1].get<std::string>())
                                      : v[1].get<uint64_t>();
    return Zp(val, prime);
}

json serialize_Zp(const Zp& z) {
    // val/prime 都用 UInt64 → 字符串编码（与 Lean 端一致，避免 JSON number 精度）
    return {{"type","Zp"}, {"val", json::array({
        std::to_string(z.number()),
        std::to_string(z.prime())
    })}};
}

// ---- SparsePolyZp ----

upolynomial_<Zp> parse_SparsePolyZp(const json& j) {
    _check_type(j, "SparsePolyZp");
    const auto& v = j.at("val");
    if (!v.is_array()) {
        throw std::runtime_error("SparsePolyZp val: expected array of [deg, [val,prime]]");
    }
    upolynomial_<Zp> result;
    for (const auto& term : v) {
        if (!term.is_array() || term.size() != 2) {
            throw std::runtime_error("SparsePolyZp term: expected [deg, [val,prime]]");
        }
        int64_t deg = term[0].get<int64_t>();
        const auto& zp_arr = term[1];
        if (!zp_arr.is_array() || zp_arr.size() != 2) {
            throw std::runtime_error("SparsePolyZp term zp: expected [val,prime]");
        }
        int64_t  val   = zp_arr[0].get<int64_t>();
        uint64_t prime = zp_arr[1].is_string() ? std::stoull(zp_arr[1].get<std::string>())
                                               : zp_arr[1].get<uint64_t>();
        result.push_back(std::make_pair(umonomial(deg), Zp(val, prime)));
    }
    return result;
}

json serialize_SparsePolyZp(const upolynomial_<Zp>& p) {
    json terms = json::array();
    for (const auto& term : p) {
        terms.push_back(json::array({
            (int64_t)term.first.deg(),
            json::array({
                std::to_string(term.second.number()),
                std::to_string(term.second.prime())
            })
        }));
    }
    return {{"type","SparsePolyZp"}, {"val", terms}};
}

}  // namespace b2b
