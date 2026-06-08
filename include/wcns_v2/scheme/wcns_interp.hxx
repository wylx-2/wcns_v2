#pragma once

#include "wcns_interp.h"
#include "wcns_v2/scheme/interp_diff.h"
#include "wcns_v2/field/field.h"
#include "wcns_v2/parallel/local_block.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

// ============================================================================
// WENO-JS constants (anonymous namespace)
// ============================================================================

namespace {

constexpr Real weno_eps   = Real(1.0e-6);   ///< ε to avoid division by zero
constexpr Real weno_power = Real(2.0);      ///< p = 2 for WENO-JS

/// Optimal linear weights for 5th-order WENO (3 sub-stencils)
constexpr Real weno_d0 = Real(1.0 / 16.0);
constexpr Real weno_d1 = Real(10.0 / 16.0);
constexpr Real weno_d2 = Real(5.0 / 16.0);

} // anonymous namespace

// ============================================================================
// WcnsInterpBase — shared helpers
// ============================================================================

inline bool WcnsInterpBase::fill_boundary(const Real* a, Int h, Int n,
                                           Real& ql, Real& qr) {
    // For boundary faces, use robust low-order interpolation that preserves
    // positivity.  High-order one-sided formulas can produce non-physical
    // values (negative density/pressure) near strong discontinuities,
    // especially at MPI split boundaries where the ghost data comes from
    // a neighboring block with a different flow state.
    //
    // Strategy: use 2nd-order linear extrapolation (weighted average of
    // nearest two cells) which is positivity-preserving for monotone data.
    // qL uses the two cells left of the face, qR uses the two cells right.

    if (h == 0) {
        // Leftmost face: between a[-1] (nonexistent) and a[0]
        // For density: just copy nearest interior cell (safe, positivity-preserving)
        ql = a[0];
        qr = a[0];
        return true;
    }
    if (h == 1) {
        // Between a[0] and a[1]
        ql = 0.5 * (a[0] + a[1]);
        qr = 0.5 * (a[0] + a[1]);
        return true;
    }
    if (h == 2) {
        // Between a[1] and a[2]
        ql = 0.5 * (a[1] + a[2]);
        qr = 0.5 * (a[1] + a[2]);
        return true;
    }
    if (h == n) {
        // Rightmost face: between a[n-1] and a[n] (nonexistent)
        ql = 1.5 * a[n-1] - 0.5 * a[n-2];  // linear extrapolation from left
        qr = 1.5 * a[n-1] - 0.5 * a[n-2];
        if (ql <= 0.0) { ql = a[n-1]; qr = a[n-1]; }
        return true;
    }
    if (h == n - 1) {
        // Between a[n-2] and a[n-1]
        ql = 0.5 * (a[n-2] + a[n-1]);
        qr = 0.5 * (a[n-2] + a[n-1]);
        return true;
    }
    if (h == n - 2) {
        // Between a[n-3] and a[n-2]
        ql = 0.5 * (a[n-3] + a[n-2]);
        qr = 0.5 * (a[n-3] + a[n-2]);
        return true;
    }
    return false;
}

inline void WcnsInterpBase::check_interp_vars(const Config& cfg,
                                               const char* scheme_name) {
    if (cfg.interp_vars != "conservative") {
        throw std::runtime_error(
            std::string(scheme_name) + ": interp_vars='" + cfg.interp_vars +
            "' not yet implemented. Only 'conservative' is supported.");
    }
}

// ============================================================================
// WcnsInterpBase — factory
// ============================================================================

inline std::unique_ptr<WcnsInterpBase> WcnsInterpBase::create(const Config& cfg) {
    if (cfg.interp_type == "weno_js") {
        return std::make_unique<WcnsNonlinearInterp>();
    } else if (cfg.interp_type == "mdcd_linear") {
        return std::make_unique<WcnsMdcdLinearInterp>();
    } else if (cfg.interp_type == "mdcd_hybrid") {
        return std::make_unique<WcnsMdcdHybridInterp>();
    } else {
        throw std::runtime_error(
            "WcnsInterpBase::create: unknown interp_type='" +
            cfg.interp_type + "'. Supported: weno_js, mdcd_linear, mdcd_hybrid.");
    }
}

// ============================================================================
// WcnsNonlinearInterp — weno_left
// ============================================================================

inline Real WcnsNonlinearInterp::weno_left(const Real* a, Int i) {
    // ---- 3rd-order sub-stencil interpolations to i+1/2 ----
    // S0 = {i-2, i-1, i}
    Real v0 = (Real(3.0) * a[i-2] - Real(10.0) * a[i-1] + Real(15.0) * a[i]) / Real(8.0);

    // S1 = {i-1, i, i+1}
    Real v1 = (Real(-1.0) * a[i-1] + Real(6.0) * a[i] + Real(3.0) * a[i+1]) / Real(8.0);

    // S2 = {i, i+1, i+2}
    Real v2 = (Real(3.0) * a[i] + Real(6.0) * a[i+1] - Real(1.0) * a[i+2]) / Real(8.0);

    // ---- Smoothness indicators (Jiang–Shu 1996) ----
    Real b0 = Real(13.0 / 12.0) * (a[i-2] - Real(2.0) * a[i-1] + a[i])
                               * (a[i-2] - Real(2.0) * a[i-1] + a[i])
            + Real(0.25) * (a[i-2] - Real(4.0) * a[i-1] + Real(3.0) * a[i])
                         * (a[i-2] - Real(4.0) * a[i-1] + Real(3.0) * a[i]);

    Real b1 = Real(13.0 / 12.0) * (a[i-1] - Real(2.0) * a[i] + a[i+1])
                               * (a[i-1] - Real(2.0) * a[i] + a[i+1])
            + Real(0.25) * (a[i-1] - a[i+1]) * (a[i-1] - a[i+1]);

    Real b2 = Real(13.0 / 12.0) * (a[i] - Real(2.0) * a[i+1] + a[i+2])
                               * (a[i] - Real(2.0) * a[i+1] + a[i+2])
            + Real(0.25) * (Real(3.0) * a[i] - Real(4.0) * a[i+1] + a[i+2])
                         * (Real(3.0) * a[i] - Real(4.0) * a[i+1] + a[i+2]);

    // ---- Nonlinear weights ----
    Real a0 = weno_d0 / ((weno_eps + b0) * (weno_eps + b0));  // d0/(ε+β0)²
    Real a1 = weno_d1 / ((weno_eps + b1) * (weno_eps + b1));
    Real a2 = weno_d2 / ((weno_eps + b2) * (weno_eps + b2));

    Real inv_sum = Real(1.0) / (a0 + a1 + a2);
    Real w0 = a0 * inv_sum;
    Real w1 = a1 * inv_sum;
    Real w2 = a2 * inv_sum;

    return w0 * v0 + w1 * v1 + w2 * v2;
}

// ============================================================================
// WcnsNonlinearInterp — weno_right
// ============================================================================

inline Real WcnsNonlinearInterp::weno_right(const Real* a, Int i) {
    // ---- 3rd-order sub-stencil interpolations to i+1/2 (mirrored) ----
    // S0' = {i+3, i+2, i+1}
    Real v0 = (Real(3.0) * a[i+3] - Real(10.0) * a[i+2] + Real(15.0) * a[i+1]) / Real(8.0);

    // S1' = {i+2, i+1, i}
    Real v1 = (Real(-1.0) * a[i+2] + Real(6.0) * a[i+1] + Real(3.0) * a[i]) / Real(8.0);

    // S2' = {i+1, i, i-1}
    Real v2 = (Real(3.0) * a[i+1] + Real(6.0) * a[i] - Real(1.0) * a[i-1]) / Real(8.0);

    // ---- Smoothness indicators (mirror-symmetric) ----
    Real b0 = Real(13.0 / 12.0) * (a[i+3] - Real(2.0) * a[i+2] + a[i+1])
                               * (a[i+3] - Real(2.0) * a[i+2] + a[i+1])
            + Real(0.25) * (a[i+3] - Real(4.0) * a[i+2] + Real(3.0) * a[i+1])
                         * (a[i+3] - Real(4.0) * a[i+2] + Real(3.0) * a[i+1]);

    Real b1 = Real(13.0 / 12.0) * (a[i+2] - Real(2.0) * a[i+1] + a[i])
                               * (a[i+2] - Real(2.0) * a[i+1] + a[i])
            + Real(0.25) * (a[i+2] - a[i]) * (a[i+2] - a[i]);

    Real b2 = Real(13.0 / 12.0) * (a[i+1] - Real(2.0) * a[i] + a[i-1])
                               * (a[i+1] - Real(2.0) * a[i] + a[i-1])
            + Real(0.25) * (Real(3.0) * a[i+1] - Real(4.0) * a[i] + a[i-1])
                         * (Real(3.0) * a[i+1] - Real(4.0) * a[i] + a[i-1]);

    // ---- Nonlinear weights (same optimal weights) ----
    Real a0 = weno_d0 / ((weno_eps + b0) * (weno_eps + b0));
    Real a1 = weno_d1 / ((weno_eps + b1) * (weno_eps + b1));
    Real a2 = weno_d2 / ((weno_eps + b2) * (weno_eps + b2));

    Real inv_sum = Real(1.0) / (a0 + a1 + a2);
    Real w0 = a0 * inv_sum;
    Real w1 = a1 * inv_sum;
    Real w2 = a2 * inv_sum;

    return w0 * v0 + w1 * v1 + w2 * v2;
}

// ============================================================================
// WcnsNonlinearInterp — interp_1d
// ============================================================================

inline void WcnsNonlinearInterp::interp_1d(const Real* a, Real* ql, Real* qr, Int n) {
    for (Int h = 0; h <= n; ++h) {
        if (fill_boundary(a, h, n, ql[h], qr[h])) continue;

        // Interior: WENO nonlinear interpolation
        // h corresponds to half-node a_{h-1/2}, i = h-1
        Int i = h - 1;
        ql[h] = weno_left(a, i);
        qr[h] = weno_right(a, i);
    }
}

// ============================================================================
// WcnsNonlinearInterp — per-direction methods
// ============================================================================

inline void WcnsNonlinearInterp::interp_xi(LocalBlock& lb, const Config& cfg) const {
    check_interp_vars(cfg, "WcnsNonlinearInterp");

    const Field& f = lb.field;
    Int nci = f.ni();
    Int ncj = f.nj();
    Int nck = f.nk();

    for (Int k = 0; k < nck; ++k) {
    for (Int j = 0; j < ncj; ++j) {
        interp_1d(&f.cons.rho(0, j, k),  &lb.field.ql_xi.rho(0, j, k),  &lb.field.qr_xi.rho(0, j, k),  nci);
        interp_1d(&f.cons.rhou(0, j, k), &lb.field.ql_xi.rhou(0, j, k), &lb.field.qr_xi.rhou(0, j, k), nci);
        interp_1d(&f.cons.rhov(0, j, k), &lb.field.ql_xi.rhov(0, j, k), &lb.field.qr_xi.rhov(0, j, k), nci);
        interp_1d(&f.cons.rhow(0, j, k), &lb.field.ql_xi.rhow(0, j, k), &lb.field.qr_xi.rhow(0, j, k), nci);
        interp_1d(&f.cons.rhoE(0, j, k), &lb.field.ql_xi.rhoE(0, j, k), &lb.field.qr_xi.rhoE(0, j, k), nci);
    }}
}

inline void WcnsNonlinearInterp::interp_eta(LocalBlock& lb, const Config& cfg) const {
    check_interp_vars(cfg, "WcnsNonlinearInterp");

    const Field& f = lb.field;
    Int nci = f.ni();
    Int ncj = f.nj();
    Int nck = f.nk();

    for (Int k = 0; k < nck; ++k) {
    for (Int i = 0; i < nci; ++i) {
        std::vector<Real> line_rho(ncj),  ql_rho(ncj + 1),  qr_rho(ncj + 1);
        std::vector<Real> line_rhou(ncj), ql_rhou(ncj + 1), qr_rhou(ncj + 1);
        std::vector<Real> line_rhov(ncj), ql_rhov(ncj + 1), qr_rhov(ncj + 1);
        std::vector<Real> line_rhow(ncj), ql_rhow(ncj + 1), qr_rhow(ncj + 1);
        std::vector<Real> line_rhoE(ncj), ql_rhoE(ncj + 1), qr_rhoE(ncj + 1);

        for (Int j = 0; j < ncj; ++j) {
            line_rho[j]  = f.cons.rho(i, j, k);
            line_rhou[j] = f.cons.rhou(i, j, k);
            line_rhov[j] = f.cons.rhov(i, j, k);
            line_rhow[j] = f.cons.rhow(i, j, k);
            line_rhoE[j] = f.cons.rhoE(i, j, k);
        }

        interp_1d(line_rho.data(),  ql_rho.data(),  qr_rho.data(),  ncj);
        interp_1d(line_rhou.data(), ql_rhou.data(), qr_rhou.data(), ncj);
        interp_1d(line_rhov.data(), ql_rhov.data(), qr_rhov.data(), ncj);
        interp_1d(line_rhow.data(), ql_rhow.data(), qr_rhow.data(), ncj);
        interp_1d(line_rhoE.data(), ql_rhoE.data(), qr_rhoE.data(), ncj);

        for (Int j = 0; j <= ncj; ++j) {
            lb.field.ql_eta.rho(i, j, k)  = ql_rho[j];
            lb.field.qr_eta.rho(i, j, k)  = qr_rho[j];
            lb.field.ql_eta.rhou(i, j, k) = ql_rhou[j];
            lb.field.qr_eta.rhou(i, j, k) = qr_rhou[j];
            lb.field.ql_eta.rhov(i, j, k) = ql_rhov[j];
            lb.field.qr_eta.rhov(i, j, k) = qr_rhov[j];
            lb.field.ql_eta.rhow(i, j, k) = ql_rhow[j];
            lb.field.qr_eta.rhow(i, j, k) = qr_rhow[j];
            lb.field.ql_eta.rhoE(i, j, k) = ql_rhoE[j];
            lb.field.qr_eta.rhoE(i, j, k) = qr_rhoE[j];
        }
    }}
}

inline void WcnsNonlinearInterp::interp_zeta(LocalBlock& lb, const Config& cfg) const {
    check_interp_vars(cfg, "WcnsNonlinearInterp");

    const Field& f = lb.field;
    Int nci = f.ni();
    Int ncj = f.nj();
    Int nck = f.nk();

    for (Int j = 0; j < ncj; ++j) {
    for (Int i = 0; i < nci; ++i) {
        std::vector<Real> line_rho(nck),  ql_rho(nck + 1),  qr_rho(nck + 1);
        std::vector<Real> line_rhou(nck), ql_rhou(nck + 1), qr_rhou(nck + 1);
        std::vector<Real> line_rhov(nck), ql_rhov(nck + 1), qr_rhov(nck + 1);
        std::vector<Real> line_rhow(nck), ql_rhow(nck + 1), qr_rhow(nck + 1);
        std::vector<Real> line_rhoE(nck), ql_rhoE(nck + 1), qr_rhoE(nck + 1);

        for (Int k = 0; k < nck; ++k) {
            line_rho[k]  = f.cons.rho(i, j, k);
            line_rhou[k] = f.cons.rhou(i, j, k);
            line_rhov[k] = f.cons.rhov(i, j, k);
            line_rhow[k] = f.cons.rhow(i, j, k);
            line_rhoE[k] = f.cons.rhoE(i, j, k);
        }

        interp_1d(line_rho.data(),  ql_rho.data(),  qr_rho.data(),  nck);
        interp_1d(line_rhou.data(), ql_rhou.data(), qr_rhou.data(), nck);
        interp_1d(line_rhov.data(), ql_rhov.data(), qr_rhov.data(), nck);
        interp_1d(line_rhow.data(), ql_rhow.data(), qr_rhow.data(), nck);
        interp_1d(line_rhoE.data(), ql_rhoE.data(), qr_rhoE.data(), nck);

        for (Int k = 0; k <= nck; ++k) {
            lb.field.ql_zeta.rho(i, j, k)  = ql_rho[k];
            lb.field.qr_zeta.rho(i, j, k)  = qr_rho[k];
            lb.field.ql_zeta.rhou(i, j, k) = ql_rhou[k];
            lb.field.qr_zeta.rhou(i, j, k) = qr_rhou[k];
            lb.field.ql_zeta.rhov(i, j, k) = ql_rhov[k];
            lb.field.qr_zeta.rhov(i, j, k) = qr_rhov[k];
            lb.field.ql_zeta.rhow(i, j, k) = ql_rhow[k];
            lb.field.qr_zeta.rhow(i, j, k) = qr_rhow[k];
            lb.field.ql_zeta.rhoE(i, j, k) = ql_rhoE[k];
            lb.field.qr_zeta.rhoE(i, j, k) = qr_rhoE[k];
        }
    }}
}

// ============================================================================
// WcnsNonlinearInterp — interp_all
// ============================================================================

inline void WcnsNonlinearInterp::interp_all(LocalBlock& lb, const Config& cfg) {
    WcnsNonlinearInterp interp;
    interp.interp_xi(lb, cfg);
    interp.interp_eta(lb, cfg);
    interp.interp_zeta(lb, cfg);
}

// ============================================================================
// WcnsMdcdLinearInterp — mdcd_left
// ============================================================================

inline Real WcnsMdcdLinearInterp::mdcd_left(const Real* s, Real diss, Real disp) {
    // 6-point left-biased MDCD linear interpolation to i+1/2
    // stencil s[0..5] = {a[i-2], a[i-1], a[i], a[i+1], a[i+2], a[i+3]}
    return (Real(3.0) * disp + Real(3.0) * diss) / Real(8.0) * s[0]
         + (Real(-18.0) * disp - Real(30.0) * diss - Real(1.0)) / Real(16.0) * s[1]
         + (Real(12.0) * disp + Real(60.0) * diss + Real(9.0)) / Real(16.0) * s[2]
         + (Real(12.0) * disp - Real(60.0) * diss + Real(9.0)) / Real(16.0) * s[3]
         + (Real(-18.0) * disp + Real(30.0) * diss - Real(1.0)) / Real(16.0) * s[4]
         + (Real(3.0) * disp - Real(3.0) * diss) / Real(8.0) * s[5];
}

// ============================================================================
// WcnsMdcdLinearInterp — interp_1d
// ============================================================================

inline void WcnsMdcdLinearInterp::interp_1d(const Real* a, Real* ql, Real* qr,
                                             Int n, Real diss, Real disp) {
    for (Int h = 0; h <= n; ++h) {
        if (WcnsInterpBase::fill_boundary(a, h, n, ql[h], qr[h])) continue;

        // Interior: 6-point MDCD linear
        // Q_L at half-node h: stencil a[h-3 .. h+2]
        ql[h] = mdcd_left(&a[h - 3], diss, disp);

        // Q_R at half-node h: mirror-reversed stencil
        //   left stencil  {a[h-3], a[h-2], a[h-1], a[h], a[h+1], a[h+2]}
        //   right stencil {a[h+2], a[h+1], a[h], a[h-1], a[h-2], a[h-3]}
        Real rev[6] = {a[h + 2], a[h + 1], a[h], a[h - 1], a[h - 2], a[h - 3]};
        qr[h] = mdcd_left(rev, diss, disp);
    }
}

// ============================================================================
// WcnsMdcdLinearInterp — per-direction methods
// ============================================================================

inline void WcnsMdcdLinearInterp::interp_xi(LocalBlock& lb, const Config& cfg) const {
    check_interp_vars(cfg, "WcnsMdcdLinearInterp");

    const Field& f = lb.field;
    Int nci = f.ni();
    Int ncj = f.nj();
    Int nck = f.nk();
    Real diss = cfg.mdcd_diss;
    Real disp = cfg.mdcd_disp;

    for (Int k = 0; k < nck; ++k) {
    for (Int j = 0; j < ncj; ++j) {
        interp_1d(&f.cons.rho(0, j, k),  &lb.field.ql_xi.rho(0, j, k),  &lb.field.qr_xi.rho(0, j, k),  nci, diss, disp);
        interp_1d(&f.cons.rhou(0, j, k), &lb.field.ql_xi.rhou(0, j, k), &lb.field.qr_xi.rhou(0, j, k), nci, diss, disp);
        interp_1d(&f.cons.rhov(0, j, k), &lb.field.ql_xi.rhov(0, j, k), &lb.field.qr_xi.rhov(0, j, k), nci, diss, disp);
        interp_1d(&f.cons.rhow(0, j, k), &lb.field.ql_xi.rhow(0, j, k), &lb.field.qr_xi.rhow(0, j, k), nci, diss, disp);
        interp_1d(&f.cons.rhoE(0, j, k), &lb.field.ql_xi.rhoE(0, j, k), &lb.field.qr_xi.rhoE(0, j, k), nci, diss, disp);
    }}
}

inline void WcnsMdcdLinearInterp::interp_eta(LocalBlock& lb, const Config& cfg) const {
    check_interp_vars(cfg, "WcnsMdcdLinearInterp");

    const Field& f = lb.field;
    Int nci = f.ni();
    Int ncj = f.nj();
    Int nck = f.nk();
    Real diss = cfg.mdcd_diss;
    Real disp = cfg.mdcd_disp;

    for (Int k = 0; k < nck; ++k) {
    for (Int i = 0; i < nci; ++i) {
        std::vector<Real> line_rho(ncj),  ql_rho(ncj + 1),  qr_rho(ncj + 1);
        std::vector<Real> line_rhou(ncj), ql_rhou(ncj + 1), qr_rhou(ncj + 1);
        std::vector<Real> line_rhov(ncj), ql_rhov(ncj + 1), qr_rhov(ncj + 1);
        std::vector<Real> line_rhow(ncj), ql_rhow(ncj + 1), qr_rhow(ncj + 1);
        std::vector<Real> line_rhoE(ncj), ql_rhoE(ncj + 1), qr_rhoE(ncj + 1);

        for (Int j = 0; j < ncj; ++j) {
            line_rho[j]  = f.cons.rho(i, j, k);
            line_rhou[j] = f.cons.rhou(i, j, k);
            line_rhov[j] = f.cons.rhov(i, j, k);
            line_rhow[j] = f.cons.rhow(i, j, k);
            line_rhoE[j] = f.cons.rhoE(i, j, k);
        }

        interp_1d(line_rho.data(),  ql_rho.data(),  qr_rho.data(),  ncj, diss, disp);
        interp_1d(line_rhou.data(), ql_rhou.data(), qr_rhou.data(), ncj, diss, disp);
        interp_1d(line_rhov.data(), ql_rhov.data(), qr_rhov.data(), ncj, diss, disp);
        interp_1d(line_rhow.data(), ql_rhow.data(), qr_rhow.data(), ncj, diss, disp);
        interp_1d(line_rhoE.data(), ql_rhoE.data(), qr_rhoE.data(), ncj, diss, disp);

        for (Int j = 0; j <= ncj; ++j) {
            lb.field.ql_eta.rho(i, j, k)  = ql_rho[j];
            lb.field.qr_eta.rho(i, j, k)  = qr_rho[j];
            lb.field.ql_eta.rhou(i, j, k) = ql_rhou[j];
            lb.field.qr_eta.rhou(i, j, k) = qr_rhou[j];
            lb.field.ql_eta.rhov(i, j, k) = ql_rhov[j];
            lb.field.qr_eta.rhov(i, j, k) = qr_rhov[j];
            lb.field.ql_eta.rhow(i, j, k) = ql_rhow[j];
            lb.field.qr_eta.rhow(i, j, k) = qr_rhow[j];
            lb.field.ql_eta.rhoE(i, j, k) = ql_rhoE[j];
            lb.field.qr_eta.rhoE(i, j, k) = qr_rhoE[j];
        }
    }}
}

inline void WcnsMdcdLinearInterp::interp_zeta(LocalBlock& lb, const Config& cfg) const {
    check_interp_vars(cfg, "WcnsMdcdLinearInterp");

    const Field& f = lb.field;
    Int nci = f.ni();
    Int ncj = f.nj();
    Int nck = f.nk();
    Real diss = cfg.mdcd_diss;
    Real disp = cfg.mdcd_disp;

    for (Int j = 0; j < ncj; ++j) {
    for (Int i = 0; i < nci; ++i) {
        std::vector<Real> line_rho(nck),  ql_rho(nck + 1),  qr_rho(nck + 1);
        std::vector<Real> line_rhou(nck), ql_rhou(nck + 1), qr_rhou(nck + 1);
        std::vector<Real> line_rhov(nck), ql_rhov(nck + 1), qr_rhov(nck + 1);
        std::vector<Real> line_rhow(nck), ql_rhow(nck + 1), qr_rhow(nck + 1);
        std::vector<Real> line_rhoE(nck), ql_rhoE(nck + 1), qr_rhoE(nck + 1);

        for (Int k = 0; k < nck; ++k) {
            line_rho[k]  = f.cons.rho(i, j, k);
            line_rhou[k] = f.cons.rhou(i, j, k);
            line_rhov[k] = f.cons.rhov(i, j, k);
            line_rhow[k] = f.cons.rhow(i, j, k);
            line_rhoE[k] = f.cons.rhoE(i, j, k);
        }

        interp_1d(line_rho.data(),  ql_rho.data(),  qr_rho.data(),  nck, diss, disp);
        interp_1d(line_rhou.data(), ql_rhou.data(), qr_rhou.data(), nck, diss, disp);
        interp_1d(line_rhov.data(), ql_rhov.data(), qr_rhov.data(), nck, diss, disp);
        interp_1d(line_rhow.data(), ql_rhow.data(), qr_rhow.data(), nck, diss, disp);
        interp_1d(line_rhoE.data(), ql_rhoE.data(), qr_rhoE.data(), nck, diss, disp);

        for (Int k = 0; k <= nck; ++k) {
            lb.field.ql_zeta.rho(i, j, k)  = ql_rho[k];
            lb.field.qr_zeta.rho(i, j, k)  = qr_rho[k];
            lb.field.ql_zeta.rhou(i, j, k) = ql_rhou[k];
            lb.field.qr_zeta.rhou(i, j, k) = qr_rhou[k];
            lb.field.ql_zeta.rhov(i, j, k) = ql_rhov[k];
            lb.field.qr_zeta.rhov(i, j, k) = qr_rhov[k];
            lb.field.ql_zeta.rhow(i, j, k) = ql_rhow[k];
            lb.field.qr_zeta.rhow(i, j, k) = qr_rhow[k];
            lb.field.ql_zeta.rhoE(i, j, k) = ql_rhoE[k];
            lb.field.qr_zeta.rhoE(i, j, k) = qr_rhoE[k];
        }
    }}
}

// ============================================================================
// WcnsMdcdHybridInterp — mdcd_hybrid_left
// ============================================================================

inline Real WcnsMdcdHybridInterp::mdcd_hybrid_left(const Real* s, Real diss,
                                                    Real disp, Real sai_ref) {
    Real f0 = s[0], f1 = s[1], f2 = s[2];
    Real f3 = s[3], f4 = s[4], f5 = s[5];

    constexpr Real eps_small = Real(1e-40);

    // ---- Discontinuity detector (sigma) ----
    Real eps = Real(0.9) * sai_ref / (Real(1.0) - Real(0.9) * sai_ref) * Real(1e-4);

    Real a1  = std::abs(f2 - f1) + std::abs(f2 - Real(2.0) * f1 + f0);
    Real b1  = std::abs(f2 - f3) + std::abs(f2 - Real(2.0) * f3 + f4);
    Real a2  = std::abs(f3 - f2) + std::abs(f3 - Real(2.0) * f2 + f1);
    Real b2  = std::abs(f3 - f4) + std::abs(f3 - Real(2.0) * f4 + f5);
    Real sai1 = (Real(2.0) * a1 * b1 + eps) / (a1 * a1 + b1 * b1 + eps);
    Real sai2 = (Real(2.0) * a2 * b2 + eps) / (a2 * a2 + b2 * b2 + eps);
    Real sai  = std::min(sai1, sai2);
    bool sigma = (sai > sai_ref);

    // ---- Discontinuous: fall back to MDCD linear ----
    if (sigma) {
        return WcnsMdcdLinearInterp::mdcd_left(s, diss, disp);
    }

    // ---- Smooth: MDCD-WENO nonlinear reconstruction ----
    // 4 sub-stencils (3-point 4th-order interpolations to i+1/2)
    Real q0 = (Real(3.0) * f0 - Real(10.0) * f1 + Real(15.0) * f2) / Real(8.0);
    Real q1 = (Real(-1.0) * f1 + Real(6.0) * f2 + Real(3.0) * f3) / Real(8.0);
    Real q2 = (Real(3.0) * f2 + Real(6.0) * f3 - Real(1.0) * f4) / Real(8.0);
    Real q3 = (Real(15.0) * f3 - Real(10.0) * f4 + Real(3.0) * f5) / Real(8.0);

    // Optimal linear weights (depend on diss, disp)
    Real d0 = Real(1.5) * (disp + diss);
    Real d1 = Real(0.5) - Real(1.5) * (disp - Real(3.0) * diss);
    Real d2 = Real(0.5) - Real(1.5) * (disp + Real(3.0) * diss);
    Real d3 = Real(1.5) * (disp - diss);

    // Smoothness indicators β_k (k = 0,1,2 for 3-point sub-stencils)
    Real t0 = f0 - Real(2.0) * f1 + f2;
    Real t1 = f1 - Real(2.0) * f2 + f3;
    Real t2 = f2 - Real(2.0) * f3 + f4;

    Real beta0 = Real(13.0 / 12.0) * t0 * t0
               + Real(0.25) * (f0 - Real(4.0) * f1 + Real(3.0) * f2)
                           * (f0 - Real(4.0) * f1 + Real(3.0) * f2);

    Real beta1 = Real(13.0 / 12.0) * t1 * t1
               + Real(0.25) * (f1 - f3) * (f1 - f3);

    Real beta2 = Real(13.0 / 12.0) * t2 * t2
               + Real(0.25) * (Real(3.0) * f2 - Real(4.0) * f3 + f4)
                           * (Real(3.0) * f2 - Real(4.0) * f3 + f4);

    // β3: 6th-order global smoothness indicator (full 6-point stencil)
    Real beta3 = (Real(271799)               * f0 * f0
        + f0 * (Real(-2380800) * f1 + Real(4086352)  * f2
              + Real(-3462252) * f3 + Real(1458762)  * f4
              + Real(-245620)  * f5)
        + f1 * (Real(5653317)  * f1 + Real(-20427884) * f2
              + Real(17905032) * f3 + Real(-7727988)  * f4
              + Real(1325006)  * f5)
        + f2 * (Real(19510972) * f2 + Real(-35817664) * f3
              + Real(15929912) * f4 + Real(-2792660)  * f5)
        + f3 * (Real(17195652) * f3 + Real(-15880404) * f4
              + Real(2863984)  * f5)
        + f4 * (Real(3824847)  * f4 + Real(-1429976)  * f5)
        + Real(139633)         * f5 * f5) / Real(120960.0);

    // τ6: reference smoothness indicator (WENO-Z style)
    Real tau6 = std::abs(beta3 - (beta0 + Real(4.0) * beta1 + beta2) / Real(6.0));

    // WENO-Z style nonlinear weights: α_k = d_k * (C + τ6/(β_k+ε))^2
    Real alpha0 = d0 * (Real(20.0) + tau6 / (beta0 + eps_small))
                    * (Real(20.0) + tau6 / (beta0 + eps_small));
    Real alpha1 = d1 * (Real(20.0) + tau6 / (beta1 + eps_small))
                    * (Real(20.0) + tau6 / (beta1 + eps_small));
    Real alpha2 = d2 * (Real(20.0) + tau6 / (beta2 + eps_small))
                    * (Real(20.0) + tau6 / (beta2 + eps_small));
    Real alpha3 = d3 * (Real(20.0) + tau6 / (beta3 + eps_small))
                    * (Real(20.0) + tau6 / (beta3 + eps_small));

    return (alpha0 * q0 + alpha1 * q1 + alpha2 * q2 + alpha3 * q3)
         / (alpha0 + alpha1 + alpha2 + alpha3);
}

// ============================================================================
// WcnsMdcdHybridInterp — interp_1d
// ============================================================================

inline void WcnsMdcdHybridInterp::interp_1d(const Real* a, Real* ql, Real* qr,
                                             Int n, Real diss, Real disp,
                                             Real sai_ref) {
    for (Int h = 0; h <= n; ++h) {
        if (WcnsInterpBase::fill_boundary(a, h, n, ql[h], qr[h])) continue;

        // Interior: 6-point MDCD-WENO hybrid
        ql[h] = mdcd_hybrid_left(&a[h - 3], diss, disp, sai_ref);

        // Q_R: mirror-reversed stencil
        Real rev[6] = {a[h + 2], a[h + 1], a[h], a[h - 1], a[h - 2], a[h - 3]};
        qr[h] = mdcd_hybrid_left(rev, diss, disp, sai_ref);
    }
}

// ============================================================================
// WcnsMdcdHybridInterp — per-direction methods
// ============================================================================

inline void WcnsMdcdHybridInterp::interp_xi(LocalBlock& lb, const Config& cfg) const {
    check_interp_vars(cfg, "WcnsMdcdHybridInterp");

    const Field& f = lb.field;
    Int nci = f.ni();
    Int ncj = f.nj();
    Int nck = f.nk();
    Real diss    = cfg.mdcd_diss;
    Real disp    = cfg.mdcd_disp;
    Real sai_ref = cfg.mdcd_sai_ref;

    for (Int k = 0; k < nck; ++k) {
    for (Int j = 0; j < ncj; ++j) {
        interp_1d(&f.cons.rho(0, j, k),  &lb.field.ql_xi.rho(0, j, k),  &lb.field.qr_xi.rho(0, j, k),  nci, diss, disp, sai_ref);
        interp_1d(&f.cons.rhou(0, j, k), &lb.field.ql_xi.rhou(0, j, k), &lb.field.qr_xi.rhou(0, j, k), nci, diss, disp, sai_ref);
        interp_1d(&f.cons.rhov(0, j, k), &lb.field.ql_xi.rhov(0, j, k), &lb.field.qr_xi.rhov(0, j, k), nci, diss, disp, sai_ref);
        interp_1d(&f.cons.rhow(0, j, k), &lb.field.ql_xi.rhow(0, j, k), &lb.field.qr_xi.rhow(0, j, k), nci, diss, disp, sai_ref);
        interp_1d(&f.cons.rhoE(0, j, k), &lb.field.ql_xi.rhoE(0, j, k), &lb.field.qr_xi.rhoE(0, j, k), nci, diss, disp, sai_ref);
    }}
}

inline void WcnsMdcdHybridInterp::interp_eta(LocalBlock& lb, const Config& cfg) const {
    check_interp_vars(cfg, "WcnsMdcdHybridInterp");

    const Field& f = lb.field;
    Int nci = f.ni();
    Int ncj = f.nj();
    Int nck = f.nk();
    Real diss    = cfg.mdcd_diss;
    Real disp    = cfg.mdcd_disp;
    Real sai_ref = cfg.mdcd_sai_ref;

    for (Int k = 0; k < nck; ++k) {
    for (Int i = 0; i < nci; ++i) {
        std::vector<Real> line_rho(ncj),  ql_rho(ncj + 1),  qr_rho(ncj + 1);
        std::vector<Real> line_rhou(ncj), ql_rhou(ncj + 1), qr_rhou(ncj + 1);
        std::vector<Real> line_rhov(ncj), ql_rhov(ncj + 1), qr_rhov(ncj + 1);
        std::vector<Real> line_rhow(ncj), ql_rhow(ncj + 1), qr_rhow(ncj + 1);
        std::vector<Real> line_rhoE(ncj), ql_rhoE(ncj + 1), qr_rhoE(ncj + 1);

        for (Int j = 0; j < ncj; ++j) {
            line_rho[j]  = f.cons.rho(i, j, k);
            line_rhou[j] = f.cons.rhou(i, j, k);
            line_rhov[j] = f.cons.rhov(i, j, k);
            line_rhow[j] = f.cons.rhow(i, j, k);
            line_rhoE[j] = f.cons.rhoE(i, j, k);
        }

        interp_1d(line_rho.data(),  ql_rho.data(),  qr_rho.data(),  ncj, diss, disp, sai_ref);
        interp_1d(line_rhou.data(), ql_rhou.data(), qr_rhou.data(), ncj, diss, disp, sai_ref);
        interp_1d(line_rhov.data(), ql_rhov.data(), qr_rhov.data(), ncj, diss, disp, sai_ref);
        interp_1d(line_rhow.data(), ql_rhow.data(), qr_rhow.data(), ncj, diss, disp, sai_ref);
        interp_1d(line_rhoE.data(), ql_rhoE.data(), qr_rhoE.data(), ncj, diss, disp, sai_ref);

        for (Int j = 0; j <= ncj; ++j) {
            lb.field.ql_eta.rho(i, j, k)  = ql_rho[j];
            lb.field.qr_eta.rho(i, j, k)  = qr_rho[j];
            lb.field.ql_eta.rhou(i, j, k) = ql_rhou[j];
            lb.field.qr_eta.rhou(i, j, k) = qr_rhou[j];
            lb.field.ql_eta.rhov(i, j, k) = ql_rhov[j];
            lb.field.qr_eta.rhov(i, j, k) = qr_rhov[j];
            lb.field.ql_eta.rhow(i, j, k) = ql_rhow[j];
            lb.field.qr_eta.rhow(i, j, k) = qr_rhow[j];
            lb.field.ql_eta.rhoE(i, j, k) = ql_rhoE[j];
            lb.field.qr_eta.rhoE(i, j, k) = qr_rhoE[j];
        }
    }}
}

inline void WcnsMdcdHybridInterp::interp_zeta(LocalBlock& lb, const Config& cfg) const {
    check_interp_vars(cfg, "WcnsMdcdHybridInterp");

    const Field& f = lb.field;
    Int nci = f.ni();
    Int ncj = f.nj();
    Int nck = f.nk();
    Real diss    = cfg.mdcd_diss;
    Real disp    = cfg.mdcd_disp;
    Real sai_ref = cfg.mdcd_sai_ref;

    for (Int j = 0; j < ncj; ++j) {
    for (Int i = 0; i < nci; ++i) {
        std::vector<Real> line_rho(nck),  ql_rho(nck + 1),  qr_rho(nck + 1);
        std::vector<Real> line_rhou(nck), ql_rhou(nck + 1), qr_rhou(nck + 1);
        std::vector<Real> line_rhov(nck), ql_rhov(nck + 1), qr_rhov(nck + 1);
        std::vector<Real> line_rhow(nck), ql_rhow(nck + 1), qr_rhow(nck + 1);
        std::vector<Real> line_rhoE(nck), ql_rhoE(nck + 1), qr_rhoE(nck + 1);

        for (Int k = 0; k < nck; ++k) {
            line_rho[k]  = f.cons.rho(i, j, k);
            line_rhou[k] = f.cons.rhou(i, j, k);
            line_rhov[k] = f.cons.rhov(i, j, k);
            line_rhow[k] = f.cons.rhow(i, j, k);
            line_rhoE[k] = f.cons.rhoE(i, j, k);
        }

        interp_1d(line_rho.data(),  ql_rho.data(),  qr_rho.data(),  nck, diss, disp, sai_ref);
        interp_1d(line_rhou.data(), ql_rhou.data(), qr_rhou.data(), nck, diss, disp, sai_ref);
        interp_1d(line_rhov.data(), ql_rhov.data(), qr_rhov.data(), nck, diss, disp, sai_ref);
        interp_1d(line_rhow.data(), ql_rhow.data(), qr_rhow.data(), nck, diss, disp, sai_ref);
        interp_1d(line_rhoE.data(), ql_rhoE.data(), qr_rhoE.data(), nck, diss, disp, sai_ref);

        for (Int k = 0; k <= nck; ++k) {
            lb.field.ql_zeta.rho(i, j, k)  = ql_rho[k];
            lb.field.qr_zeta.rho(i, j, k)  = qr_rho[k];
            lb.field.ql_zeta.rhou(i, j, k) = ql_rhou[k];
            lb.field.qr_zeta.rhou(i, j, k) = qr_rhou[k];
            lb.field.ql_zeta.rhov(i, j, k) = ql_rhov[k];
            lb.field.qr_zeta.rhov(i, j, k) = qr_rhov[k];
            lb.field.ql_zeta.rhow(i, j, k) = ql_rhow[k];
            lb.field.qr_zeta.rhow(i, j, k) = qr_rhow[k];
            lb.field.ql_zeta.rhoE(i, j, k) = ql_rhoE[k];
            lb.field.qr_zeta.rhoE(i, j, k) = qr_rhoE[k];
        }
    }}
}
