/**
 * @file crosscheck_flint.hh
 * @brief CLPoly <-> FLINT conversion utilities for cross-library testing.
 *
 * Conversions use string round-trip (no performance requirement in tests).
 * Requires: libflint-dev (FLINT 3.x)
 */
#ifndef CROSSCHECK_FLINT_HH
#define CROSSCHECK_FLINT_HH

#include <clpoly/clpoly.hh>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpq_mpoly.h>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <stdexcept>

namespace crosscheck {

// ---- ZZ <-> fmpz_t conversion (string round-trip) ----

inline void clpoly_zz_to_fmpz(fmpz_t out, const clpoly::ZZ& z) {
    std::ostringstream ss;
    ss << z;
    std::string s = ss.str();
    if (fmpz_set_str(out, s.c_str(), 10) != 0)
        throw std::runtime_error("fmpz_set_str failed for: " + s);
}

inline clpoly::ZZ fmpz_to_clpoly_zz(const fmpz_t z) {
    char* str = fmpz_get_str(nullptr, 10, z);
    clpoly::ZZ result(str);
    flint_free(str);
    return result;
}

// ---- RAII wrapper for fmpz_mpoly_t + fmpz_mpoly_ctx_t ----

struct FlintPoly {
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_t poly;
    std::vector<clpoly::variable> vars;  // variable ordering
    bool _owns;  // true if this instance owns ctx/poly

    FlintPoly(const std::vector<clpoly::variable>& v)
        : vars(v), _owns(true)
    {
        fmpz_mpoly_ctx_init(ctx, static_cast<slong>(v.size()), ORD_DEGLEX);
        fmpz_mpoly_init(poly, ctx);
    }

    ~FlintPoly() {
        if (_owns) {
            fmpz_mpoly_clear(poly, ctx);
            fmpz_mpoly_ctx_clear(ctx);
        }
    }

    // Move constructor
    FlintPoly(FlintPoly&& other) noexcept
        : vars(std::move(other.vars)), _owns(other._owns)
    {
        std::memcpy(ctx, other.ctx, sizeof(fmpz_mpoly_ctx_t));
        std::memcpy(poly, other.poly, sizeof(fmpz_mpoly_t));
        other._owns = false;
    }

    FlintPoly(const FlintPoly&) = delete;
    FlintPoly& operator=(const FlintPoly&) = delete;
    FlintPoly& operator=(FlintPoly&&) = delete;

    slong length() const { return fmpz_mpoly_length(poly, ctx); }
};

// ---- polynomial_ZZ -> FlintPoly ----

inline FlintPoly clpoly_to_flint(
    const clpoly::polynomial_ZZ& p,
    const std::vector<clpoly::variable>& vars)
{
    // Build variable -> index map
    std::map<clpoly::variable, slong> var_idx;
    for (slong i = 0; i < static_cast<slong>(vars.size()); ++i)
        var_idx[vars[i]] = i;

    FlintPoly fp(vars);
    slong nvars = static_cast<slong>(vars.size());
    std::vector<ulong> exp(nvars, 0);
    fmpz_t coeff;
    fmpz_init(coeff);

    for (auto& term : p) {
        // Build exponent vector
        std::fill(exp.begin(), exp.end(), 0);
        for (auto& ve : term.first) {
            auto it = var_idx.find(ve.first);
            if (it != var_idx.end())
                exp[it->second] = static_cast<ulong>(ve.second);
        }
        // Convert coefficient
        clpoly_zz_to_fmpz(coeff, term.second);
        fmpz_mpoly_push_term_fmpz_ui(fp.poly, coeff, exp.data(), fp.ctx);
    }

    fmpz_clear(coeff);
    fmpz_mpoly_sort_terms(fp.poly, fp.ctx);
    fmpz_mpoly_combine_like_terms(fp.poly, fp.ctx);
    return fp;
}

// ---- FlintPoly -> polynomial_ZZ ----

inline clpoly::polynomial_ZZ flint_to_clpoly(const FlintPoly& fp) {
    slong nvars = static_cast<slong>(fp.vars.size());
    slong len = fp.length();

    std::vector<std::pair<clpoly::monomial, clpoly::ZZ>> terms;
    terms.reserve(len);

    std::vector<ulong> exp(nvars);
    fmpz_t coeff;
    fmpz_init(coeff);

    for (slong i = 0; i < len; ++i) {
        fmpz_mpoly_get_term_coeff_fmpz(coeff, fp.poly, i, fp.ctx);
        fmpz_mpoly_get_term_exp_ui(exp.data(), fp.poly, i, fp.ctx);

        clpoly::ZZ c = fmpz_to_clpoly_zz(coeff);
        std::vector<std::pair<clpoly::variable, int64_t>> mono_data;
        for (slong j = 0; j < nvars; ++j) {
            if (exp[j] != 0)
                mono_data.push_back({fp.vars[j], static_cast<int64_t>(exp[j])});
        }
        terms.emplace_back(clpoly::monomial(mono_data), std::move(c));
    }

    fmpz_clear(coeff);
    return clpoly::polynomial_ZZ(terms);
}

// ---- Collect variables from one or two polynomials ----

inline std::vector<clpoly::variable> collect_vars(
    const clpoly::polynomial_ZZ& f,
    const clpoly::polynomial_ZZ& g = clpoly::polynomial_ZZ())
{
    auto vlist = clpoly::get_variables(f);
    if (!g.empty()) {
        auto vlist2 = clpoly::get_variables(g);
        for (auto& v : vlist2) {
            bool found = false;
            for (auto& u : vlist)
                if (u.first == v.first) { found = true; break; }
            if (!found)
                vlist.push_back(v);
        }
    }
    // Sort by variable serial for deterministic ordering
    vlist.sort([](const std::pair<clpoly::variable, int64_t>& a,
                  const std::pair<clpoly::variable, int64_t>& b) {
        return a.first < b.first;
    });
    std::vector<clpoly::variable> result;
    for (auto& v : vlist)
        result.push_back(v.first);
    return result;
}

// ======== QQ <-> fmpq_t conversion ========

inline void clpoly_qq_to_fmpq(fmpq_t out, const clpoly::QQ& q) {
    fmpz_t num, den;
    fmpz_init(num);
    fmpz_init(den);
    clpoly_zz_to_fmpz(num, q.get_num());
    clpoly_zz_to_fmpz(den, q.get_den());
    fmpq_set_fmpz_frac(out, num, den);
    fmpz_clear(num);
    fmpz_clear(den);
}

inline clpoly::QQ fmpq_to_clpoly_qq(const fmpq_t q) {
    clpoly::ZZ num = fmpz_to_clpoly_zz(fmpq_numref(q));
    clpoly::ZZ den = fmpz_to_clpoly_zz(fmpq_denref(q));
    return clpoly::QQ(num, den);
}

// ---- RAII wrapper for fmpq_mpoly_t + fmpq_mpoly_ctx_t ----

struct FlintQPoly {
    fmpq_mpoly_ctx_t ctx;
    fmpq_mpoly_t poly;
    std::vector<clpoly::variable> vars;
    bool _owns;

    FlintQPoly(const std::vector<clpoly::variable>& v)
        : vars(v), _owns(true)
    {
        fmpq_mpoly_ctx_init(ctx, static_cast<slong>(v.size()), ORD_DEGLEX);
        fmpq_mpoly_init(poly, ctx);
    }

    ~FlintQPoly() {
        if (_owns) {
            fmpq_mpoly_clear(poly, ctx);
            fmpq_mpoly_ctx_clear(ctx);
        }
    }

    FlintQPoly(FlintQPoly&& other) noexcept
        : vars(std::move(other.vars)), _owns(other._owns)
    {
        std::memcpy(ctx, other.ctx, sizeof(fmpq_mpoly_ctx_t));
        std::memcpy(poly, other.poly, sizeof(fmpq_mpoly_t));
        other._owns = false;
    }

    FlintQPoly(const FlintQPoly&) = delete;
    FlintQPoly& operator=(const FlintQPoly&) = delete;
    FlintQPoly& operator=(FlintQPoly&&) = delete;

    slong length() const { return fmpq_mpoly_length(poly, ctx); }
};

// ---- polynomial_QQ -> FlintQPoly ----

inline FlintQPoly clpoly_qq_to_flint(
    const clpoly::polynomial_QQ& p,
    const std::vector<clpoly::variable>& vars)
{
    std::map<clpoly::variable, slong> var_idx;
    for (slong i = 0; i < static_cast<slong>(vars.size()); ++i)
        var_idx[vars[i]] = i;

    FlintQPoly fp(vars);
    slong nvars = static_cast<slong>(vars.size());
    std::vector<ulong> exp(nvars, 0);
    fmpq_t coeff;
    fmpq_init(coeff);

    for (auto& term : p) {
        std::fill(exp.begin(), exp.end(), 0);
        for (auto& ve : term.first) {
            auto it = var_idx.find(ve.first);
            if (it != var_idx.end())
                exp[it->second] = static_cast<ulong>(ve.second);
        }
        clpoly_qq_to_fmpq(coeff, term.second);
        fmpq_mpoly_push_term_fmpq_ui(fp.poly, coeff, exp.data(), fp.ctx);
    }

    fmpq_clear(coeff);
    fmpq_mpoly_sort_terms(fp.poly, fp.ctx);
    fmpq_mpoly_reduce(fp.poly, fp.ctx);
    return fp;
}

// ---- FlintQPoly -> polynomial_QQ ----

inline clpoly::polynomial_QQ flint_qq_to_clpoly(const FlintQPoly& fp) {
    slong nvars = static_cast<slong>(fp.vars.size());
    slong len = fp.length();

    std::vector<std::pair<clpoly::monomial, clpoly::QQ>> terms;
    terms.reserve(len);

    std::vector<ulong> exp(nvars);
    fmpq_t coeff;
    fmpq_init(coeff);

    for (slong i = 0; i < len; ++i) {
        fmpq_mpoly_get_term_coeff_fmpq(coeff, fp.poly, i, fp.ctx);
        fmpq_mpoly_get_term_exp_ui(exp.data(), fp.poly, i, fp.ctx);

        clpoly::QQ c = fmpq_to_clpoly_qq(coeff);
        std::vector<std::pair<clpoly::variable, int64_t>> mono_data;
        for (slong j = 0; j < nvars; ++j) {
            if (exp[j] != 0)
                mono_data.push_back({fp.vars[j], static_cast<int64_t>(exp[j])});
        }
        terms.emplace_back(clpoly::monomial(mono_data), std::move(c));
    }

    fmpq_clear(coeff);
    return clpoly::polynomial_QQ(terms);
}

// ---- collect_vars overload for polynomial_QQ ----

inline std::vector<clpoly::variable> collect_vars_qq(
    const clpoly::polynomial_QQ& f,
    const clpoly::polynomial_QQ& g = clpoly::polynomial_QQ())
{
    auto vlist = clpoly::get_variables(f);
    if (!g.empty()) {
        auto vlist2 = clpoly::get_variables(g);
        for (auto& v : vlist2) {
            bool found = false;
            for (auto& u : vlist)
                if (u.first == v.first) { found = true; break; }
            if (!found)
                vlist.push_back(v);
        }
    }
    vlist.sort([](const std::pair<clpoly::variable, int64_t>& a,
                  const std::pair<clpoly::variable, int64_t>& b) {
        return a.first < b.first;
    });
    std::vector<clpoly::variable> result;
    for (auto& v : vlist)
        result.push_back(v.first);
    return result;
}

} // namespace crosscheck

#endif // CROSSCHECK_FLINT_HH
