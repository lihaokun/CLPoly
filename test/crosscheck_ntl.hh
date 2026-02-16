/**
 * @file crosscheck_ntl.hh
 * @brief CLPoly <-> NTL conversion utilities for cross-library testing.
 *
 * Conversions use string round-trip (no performance requirement in tests).
 * Requires: libntl-dev
 *
 * NOTE: Both CLPoly and NTL define a `ZZ` class. All references use
 * fully qualified names: clpoly::ZZ vs NTL::ZZ.
 */
#ifndef CROSSCHECK_NTL_HH
#define CROSSCHECK_NTL_HH

#include <clpoly/clpoly.hh>
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <sstream>
#include <string>
#include <vector>

namespace crosscheck {

// ---- clpoly::ZZ <-> NTL::ZZ conversion (string round-trip) ----

inline NTL::ZZ clpoly_zz_to_ntl(const clpoly::ZZ& z) {
    std::ostringstream ss;
    ss << z;
    NTL::ZZ result;
    std::istringstream is(ss.str());
    is >> result;
    return result;
}

inline clpoly::ZZ ntl_zz_to_clpoly(const NTL::ZZ& z) {
    std::ostringstream ss;
    ss << z;
    return clpoly::ZZ(ss.str());
}

// ---- upolynomial_ZZ -> NTL::ZZX ----

inline NTL::ZZX clpoly_upoly_to_ntl(const clpoly::upolynomial_ZZ& p) {
    NTL::ZZX f;
    for (auto& term : p) {
        long deg = static_cast<long>(term.first.deg());
        NTL::ZZ coeff = clpoly_zz_to_ntl(term.second);
        NTL::SetCoeff(f, deg, coeff);
    }
    return f;
}

// ---- NTL::ZZX -> upolynomial_ZZ ----

inline clpoly::upolynomial_ZZ ntl_to_clpoly_upoly(const NTL::ZZX& f) {
    std::vector<std::pair<clpoly::umonomial, clpoly::ZZ>> terms;
    long d = NTL::deg(f);
    for (long i = d; i >= 0; --i) {
        const NTL::ZZ& c = NTL::coeff(f, i);
        if (!NTL::IsZero(c)) {
            terms.emplace_back(clpoly::umonomial(i), ntl_zz_to_clpoly(c));
        }
    }
    return clpoly::upolynomial_ZZ(terms);
}

// ---- NTL point evaluation: evaluate ZZX at NTL::ZZ (Horner's method) ----

inline NTL::ZZ ntl_eval(const NTL::ZZX& f, const NTL::ZZ& a) {
    long d = NTL::deg(f);
    if (d < 0) return NTL::ZZ::zero();
    NTL::ZZ result = NTL::coeff(f, d);
    for (long i = d - 1; i >= 0; --i) {
        result = result * a + NTL::coeff(f, i);
    }
    return result;
}

// ---- NTL derivative ----

inline NTL::ZZX ntl_diff(const NTL::ZZX& f) {
    return NTL::diff(f);
}

} // namespace crosscheck

#endif // CROSSCHECK_NTL_HH
