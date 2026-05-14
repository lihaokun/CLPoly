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
// UInt64 用字符串避免 IEEE754 精度（>2^53），与 Lean encodeUInt64 一致
json serialize_UInt64(uint64_t v){ return {{"type","UInt64"}, {"val", std::to_string(v)}}; }

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

// ---- BoolZZ ----

json serialize_BoolZZ(bool ok, const ZZ& z) {
    std::ostringstream os;
    os << z;
    return {{"type","BoolZZ"}, {"val", json::array({ ok, os.str() })}};
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
        int64_t  val   = zp_arr[0].is_string() ? std::stoll(zp_arr[0].get<std::string>())
                                               : zp_arr[0].get<int64_t>();
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

// ---- SparsePolyZZ ----

upolynomial_<ZZ> parse_SparsePolyZZ(const json& j) {
    _check_type(j, "SparsePolyZZ");
    const auto& v = j.at("val");
    if (!v.is_array()) {
        throw std::runtime_error("SparsePolyZZ val: expected array of [deg, zz_decimal_str]");
    }
    upolynomial_<ZZ> result;
    for (const auto& term : v) {
        if (!term.is_array() || term.size() != 2) {
            throw std::runtime_error("SparsePolyZZ term: expected [deg, zz_decimal_str]");
        }
        int64_t deg = term[0].get<int64_t>();
        ZZ coef = term[1].is_string()
            ? ZZ(term[1].get<std::string>())
            : ZZ(term[1].get<int64_t>());
        result.push_back(std::make_pair(umonomial(deg), std::move(coef)));
    }
    return result;
}

json serialize_SparsePolyZZ(const upolynomial_<ZZ>& p) {
    json terms = json::array();
    for (const auto& term : p) {
        std::ostringstream os;
        os << term.second;
        terms.push_back(json::array({
            (int64_t)term.first.deg(),
            os.str()
        }));
    }
    return {{"type","SparsePolyZZ"}, {"val", terms}};
}

// ---- 复合返回类型 ----

json serialize_PairSPZp(const upolynomial_<Zp>& q, const upolynomial_<Zp>& r) {
    return {{"type","PairSPZp"}, {"val", json::array({
        serialize_SparsePolyZp(q),
        serialize_SparsePolyZp(r)
    })}};
}

json serialize_TripleSPZp(const upolynomial_<Zp>& g,
                          const upolynomial_<Zp>& s,
                          const upolynomial_<Zp>& t) {
    return {{"type","TripleSPZp"}, {"val", json::array({
        serialize_SparsePolyZp(g),
        serialize_SparsePolyZp(s),
        serialize_SparsePolyZp(t)
    })}};
}

// ---- A5: Variable ----

Variable parse_Variable(const json& j) {
    _check_type(j, "Variable");
    const auto& v = j.at("val");
    uint64_t serial = v.is_string() ? std::stoull(v.get<std::string>())
                                    : v.get<uint64_t>();
    return Variable::get_variable((std::size_t)serial);
}

json serialize_Variable(const Variable& v) {
    return {{"type","Variable"}, {"val", std::to_string((uint64_t)v.serial())}};
}

// ---- A5: Monomial ----
// 内部 helper：从 JSON array 构造 Monomial（不带 type 包装，用于 MvPoly 嵌套）
static Monomial _parse_Monomial_inner(const json& v) {
    Monomial m;
    if (!v.is_array()) throw std::runtime_error("Monomial val: expected array");
    for (const auto& entry : v) {
        if (!entry.is_array() || entry.size() != 2) {
            throw std::runtime_error("Monomial entry: expected [var_serial, exp]");
        }
        uint64_t serial = entry[0].is_string() ? std::stoull(entry[0].get<std::string>())
                                                : entry[0].get<uint64_t>();
        int64_t exp = entry[1].is_string() ? std::stoll(entry[1].get<std::string>())
                                            : entry[1].get<int64_t>();
        m.push_back({Variable::get_variable((std::size_t)serial), exp});
    }
    m.normalization();
    return m;
}

static json _serialize_Monomial_inner(const Monomial& m) {
    json entries = json::array();
    for (const auto& entry : m) {
        entries.push_back(json::array({
            std::to_string((uint64_t)entry.first.serial()),
            (int64_t)entry.second
        }));
    }
    return entries;
}

Monomial parse_Monomial(const json& j) {
    _check_type(j, "Monomial");
    return _parse_Monomial_inner(j.at("val"));
}

json serialize_Monomial(const Monomial& m) {
    return {{"type","Monomial"}, {"val", _serialize_Monomial_inner(m)}};
}

// ---- A5: MvPolyZZ / MvPolyZp ----

MvPolyZZ parse_MvPolyZZ(const json& j) {
    _check_type(j, "MvPolyZZ");
    const auto& v = j.at("val");
    if (!v.is_array()) throw std::runtime_error("MvPolyZZ val: expected array");
    MvPolyZZ p;
    for (const auto& term : v) {
        if (!term.is_array() || term.size() != 2) {
            throw std::runtime_error("MvPolyZZ term: expected [mono, zz_str]");
        }
        Monomial m = _parse_Monomial_inner(term[0]);
        ZZ coef = term[1].is_string() ? ZZ(term[1].get<std::string>())
                                       : ZZ(term[1].get<int64_t>());
        p.push_back({std::move(m), std::move(coef)});
    }
    p.normalization();
    return p;
}

json serialize_MvPolyZZ(const MvPolyZZ& p) {
    json terms = json::array();
    for (const auto& term : p) {
        std::ostringstream os;
        os << term.second;
        terms.push_back(json::array({
            _serialize_Monomial_inner(term.first),
            os.str()
        }));
    }
    return {{"type","MvPolyZZ"}, {"val", terms}};
}

MvPolyZp parse_MvPolyZp(const json& j) {
    _check_type(j, "MvPolyZp");
    const auto& v = j.at("val");
    if (!v.is_array()) throw std::runtime_error("MvPolyZp val: expected array");
    MvPolyZp p;
    for (const auto& term : v) {
        if (!term.is_array() || term.size() != 2) {
            throw std::runtime_error("MvPolyZp term: expected [mono, [val,prime]]");
        }
        Monomial m = _parse_Monomial_inner(term[0]);
        const auto& zp_arr = term[1];
        if (!zp_arr.is_array() || zp_arr.size() != 2) {
            throw std::runtime_error("MvPolyZp term zp: expected [val,prime]");
        }
        int64_t val = zp_arr[0].is_string() ? std::stoll(zp_arr[0].get<std::string>())
                                             : zp_arr[0].get<int64_t>();
        uint64_t prime = zp_arr[1].is_string() ? std::stoull(zp_arr[1].get<std::string>())
                                                : zp_arr[1].get<uint64_t>();
        p.push_back({std::move(m), Zp(val, prime)});
    }
    p.normalization();
    return p;
}

json serialize_MvPolyZp(const MvPolyZp& p) {
    json terms = json::array();
    for (const auto& term : p) {
        terms.push_back(json::array({
            _serialize_Monomial_inner(term.first),
            json::array({
                std::to_string(term.second.number()),
                std::to_string(term.second.prime())
            })
        }));
    }
    return {{"type","MvPolyZp"}, {"val", terms}};
}

// ---- A5: VarMap (Variable → ZZ) ----

std::map<Variable, ZZ> parse_VarMapZZ(const json& j) {
    _check_type(j, "VarMap");
    const auto& v = j.at("val");
    if (!v.is_array()) throw std::runtime_error("VarMap val: expected array");
    std::map<Variable, ZZ> result;
    for (const auto& entry : v) {
        if (!entry.is_array() || entry.size() != 2) {
            throw std::runtime_error("VarMap entry: expected [var_serial, zz_str]");
        }
        uint64_t serial = entry[0].is_string() ? std::stoull(entry[0].get<std::string>())
                                                : entry[0].get<uint64_t>();
        ZZ val = entry[1].is_string() ? ZZ(entry[1].get<std::string>())
                                       : ZZ(entry[1].get<int64_t>());
        result.emplace(Variable::get_variable((std::size_t)serial), std::move(val));
    }
    return result;
}

// ---- A6: ZZTriple ----

json serialize_ZZTriple(const ZZ& g, const ZZ& s, const ZZ& t) {
    std::ostringstream gs, ss, ts;
    gs << g; ss << s; ts << t;
    return {{"type","ZZTriple"}, {"val", json::array({gs.str(), ss.str(), ts.str()})}};
}

}  // namespace b2b
