#ifndef CDFCI_ENERGY_CORRECTION_H
#define CDFCI_ENERGY_CORRECTION_H

#include "config.h"
#include "ip_entry.h"
#include <algorithm>
#include <cmath>
#include <type_traits>

namespace correction_detail {
template<class W, class = void> struct compact_ip_mode : std::false_type {};
template<class W>
struct compact_ip_mode<W, std::void_t<decltype(W::uses_compact_ip)>>
    : std::bool_constant<W::uses_compact_ip> {};

template<class W, class Entry>
NumericalType coefficient(const W &wf, const Entry &entry)
{
    if constexpr (compact_ip_mode<W>::value)
        return wf.stored_c(entry.first, entry.second);
    else return ipentry::c(entry.second);
}
}

// Energy-only projected diagonal/Olsen correction (Algorithm 1, Eqs. 15-19
// in cdfci_true_perturbation.pdf). No H application, determinant insertion,
// coefficient update, diagonal cache, or correction vector is needed.
//
// For unnormalised x, put a_i = (z_p/x_p) x_i - z_i, m_i = H_ii - E.
// With A=sum(a_i^2/m_i), B=sum(a_i*x_i/m_i), C=sum(x_i^2/m_i), i != p,
// delta_internal = (-A + m_p*B^2/(x_p^2 + m_p*C)) / (x'x).
// x_i=0 entries affect only A, so their contribution can be reported separately
// as the predecessor's stored-space external PT2. Compressed external z entries
// give an approximate external correction, not full-space EN-PT2.
template <typename H, typename W>
EnergyCorrectionResult compute_energy_correction(const H &ham, const W &wf,
                                                 bool include_external = true,
                                                 bool include_internal = true)
{
    static_assert(W::NSTATES_val == 1, "Energy correction currently supports one state");
    using Quad = QUAD_PRECISION;
    using Key = typename W::key_type;
    EnergyCorrectionResult out;
    const Quad norm = wf.get_xx();
    if (!(norm > 0) || !std::isfinite(static_cast<double>(norm))) {
        out.status = "invalid_norm";
        return out;
    }
    const Quad energy = wf.get_xz() / norm;
    out.variational_energy = static_cast<double>(energy);

    Key pivot{};
    NumericalType xp = 0, zp = 0;
    double pivot_diagonal = 0;
    bool pivot_diagonal_cached = false;
    // Select a large coefficient without evaluating any diagonals.
    wf.loop([&](const auto &entry) {
        auto x = correction_detail::coefficient(wf, entry), z = ipentry::b(entry.second);
        if (include_internal && std::abs(x) > std::abs(xp)) {
            pivot = entry.first;
            xp = x;
            zp = z;
            pivot_diagonal_cached = ipentry::has_hii(entry.second);
            if (pivot_diagonal_cached)
                pivot_diagonal = ipentry::hii(entry.second);
        }
    });

    double dp = 0;
    if (include_internal) {
        if (pivot_diagonal_cached) {
            dp = pivot_diagonal;
        } else {
            dp = ham.get_diagonal(pivot);
            ++out.diagonal_evaluations;
        }
    }
    const Quad mp = static_cast<Quad>(dp) - energy;
    const Quad ratio = include_internal ? static_cast<Quad>(zp) / xp : 0;
    Quad aa = 0, ax = 0, xx = 0, external = 0;
    Quad residual_internal = 0, residual_external = 0;
    wf.loop([&](const auto &entry) {
        const double x = correction_detail::coefficient(wf, entry), z = ipentry::b(entry.second);
        if (x == 0 && (z == 0 || !include_external)) return;
        const Quad residual = energy * x - z;
        if (x != 0) {
            if (include_internal) residual_internal += residual * residual;
        } else {
            residual_external += residual * residual;
        }
        if (x != 0 && !include_internal) return;
        if (include_internal && typename W::key_equal{}(entry.first, pivot)) return;
        auto det = entry.first;
        const bool diagonal_cached = ipentry::has_hii(entry.second);
        const double diagonal = diagonal_cached
            ? ipentry::hii(entry.second) : ham.get_diagonal(det);
        if (!diagonal_cached) ++out.diagonal_evaluations;
        const Quad gap = static_cast<Quad>(diagonal) - energy;
        if (x == 0) {
            external -= static_cast<Quad>(z) * z / gap;
        } else {
            const Quad a = ratio * x - z;
            aa += a * a / gap;
            ax += a * x / gap;
            xx += static_cast<Quad>(x) * x / gap;
        }
    });
    out.internal_residual_norm = static_cast<double>(sqrt(residual_internal / norm));
    out.stored_residual_norm = static_cast<double>(sqrt((residual_internal + residual_external) / norm));
    out.external_valid = true;
    out.external_correction = static_cast<double>(external / norm);
    if (!include_internal) {
        out.internal_correction = 0;
        out.internal_valid = true;
        out.corrected_energy = static_cast<double>(energy + external / norm);
        out.valid = true;
        out.status = "unchecked";
        return out;
    }
    const Quad pivot_square = static_cast<Quad>(xp) * xp;
    const Quad denominator = pivot_square + mp * xx;
    const Quad internal = (-aa + mp * ax * ax / denominator) / norm;
    out.internal_correction = static_cast<double>(internal);
    out.internal_valid = true;
    out.valid = true;
    out.corrected_energy = static_cast<double>(energy + internal + external / norm);
    out.status = "unchecked";
    return out;
}

// The IP solver keeps the internal rows in a compact array as they change.
// Olsen still uses the current Rayleigh energy at each report, but neither
// pivot selection nor the projected sums visit the external hash table.
template <typename W>
EnergyCorrectionResult compute_cached_internal_energy_correction(const W &wf)
{
    static_assert(W::NSTATES_val == 1, "Olsen correction supports one state");
    using Quad = QUAD_PRECISION;
    EnergyCorrectionResult out;
    const Quad norm = wf.get_xx();
    if (!(norm > 0) || !std::isfinite(static_cast<double>(norm))) {
        out.status = "invalid_norm";
        return out;
    }
    const Quad energy = wf.get_xz() / norm;
    out.variational_energy = static_cast<double>(energy);
    const auto &rows = wf.ip_internal_rows();
    size_t pivot = rows.size();
    NumericalType xp = 0;
    for (size_t i = 0; i < rows.size(); ++i) {
        if (std::abs(rows[i].c) > std::abs(xp)) {
            pivot = i;
            xp = rows[i].c;
        }
    }
    if (pivot == rows.size()) {
        out.status = "invalid_internal_cache";
        return out;
    }
    const Quad mp = static_cast<Quad>(rows[pivot].hii) - energy;
    const Quad ratio = static_cast<Quad>(rows[pivot].b) / xp;
    Quad aa = 0, ax = 0, xx = 0, residual_internal = 0;
    for (size_t i = 0; i < rows.size(); ++i) {
        const auto &row = rows[i];
        if (row.c == 0) continue;
        const Quad residual = energy * row.c - row.b;
        residual_internal += residual * residual;
        if (i == pivot) continue;
        const Quad gap = static_cast<Quad>(row.hii) - energy;
        const Quad a = ratio * row.c - row.b;
        aa += a * a / gap;
        ax += a * row.c / gap;
        xx += static_cast<Quad>(row.c) * row.c / gap;
    }
    const Quad pivot_square = static_cast<Quad>(xp) * xp;
    const Quad denominator = pivot_square + mp * xx;
    const Quad internal = (-aa + mp * ax * ax / denominator) / norm;
    out.internal_residual_norm = static_cast<double>(sqrt(residual_internal / norm));
    out.stored_residual_norm = out.internal_residual_norm;
    out.internal_correction = static_cast<double>(internal);
    out.internal_valid = true;
    out.external_valid = true;
    out.valid = true;
    out.corrected_energy = static_cast<double>(energy + internal);
    out.status = "unchecked";
    return out;
}

#endif
