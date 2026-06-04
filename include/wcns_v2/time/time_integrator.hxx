#pragma once

#include "time_integrator.h"
#include "wcns_v2/parallel/local_block.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// ============================================================================
// n_stages
// ============================================================================

inline int TimeIntegrator::n_stages(const Config& cfg) {
    if (cfg.time_scheme == "rk3-tvd") return 3;
    if (cfg.time_scheme == "rk4")     return 4;
    if (cfg.time_scheme == "lu-sgs")  return 1;
    throw std::runtime_error("TimeIntegrator: unknown time_scheme \""
                             + cfg.time_scheme + "\"");
}

// ============================================================================
// rk_coeffs — returns alpha, beta for the given stage (0-indexed)
// ============================================================================

inline void TimeIntegrator::rk_coeffs(const Config& cfg, int stage,
                                       Real& alpha, Real& beta) {
    if (cfg.time_scheme == "rk3-tvd") {
        // RK3-TVD (Shu-Osher)
        switch (stage) {
        case 0: alpha = 0.0;       beta = 1.0;        break;
        case 1: alpha = 3.0/4.0;   beta = 1.0/4.0;   break;
        case 2: alpha = 1.0/3.0;   beta = 2.0/3.0;   break;
        default:
            throw std::runtime_error("TimeIntegrator: stage " +
                                     std::to_string(stage) + " out of range for rk3-tvd");
        }
    } else if (cfg.time_scheme == "rk4") {
        // Classical RK4
        switch (stage) {
        case 0: alpha = 1.0;       beta = 1.0/2.0;   break;
        case 1: alpha = 1.0;       beta = 1.0/2.0;   break;
        case 2: alpha = 1.0;       beta = 1.0;        break;
        case 3: alpha = 0.0;       beta = 1.0/6.0;   break;
        default:
            throw std::runtime_error("TimeIntegrator: stage " +
                                     std::to_string(stage) + " out of range for rk4");
        }
    } else if (cfg.time_scheme == "lu-sgs") {
        alpha = 0.0; beta = 1.0;
    } else {
        throw std::runtime_error("TimeIntegrator: unknown time_scheme \""
                                 + cfg.time_scheme + "\"");
    }
}

// ============================================================================
// advance_stage — single RK stage update
// ============================================================================

inline void TimeIntegrator::advance_stage(LocalBlock& lb, const Config& cfg,
                                           Real dt, int stage,
                                           const ConservativeVars& Q0) {
    Real alpha, beta;
    rk_coeffs(cfg, stage, alpha, beta);

    auto& f     = lb.field;
    const Int nci = f.ni();
    const Int ncj = f.nj();
    const Int nck = f.nk();

    // Interior cells only (same range as RHS computation)
    const Int i0 = 3, i1 = nci - 4;
    const Int j0 = 3, j1 = ncj - 4;
    const Int k0 = 3, k1 = nck - 4;

    // Q_new = alpha * Q0 + beta * (Q_old + dt * rhs)
    //   where rhs = ∂Q/∂t = -(∂F̃/∂ξ + ∂G̃/∂η + ∂H̃/∂ζ) / J
    //   (SCMM: J·∂Q/∂t + ∂F̃/∂ξ + ... = 0 ⇒ ∂Q/∂t = -(∂F̃/∂ξ + ...)/J)
    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        // Stage updates for each conservative component
        f.cons.rho(i,j,k)  = alpha * Q0.rho(i,j,k)
            + beta * (f.cons.rho(i,j,k)  + dt * f.rhs.rho(i,j,k));
        f.cons.rhou(i,j,k) = alpha * Q0.rhou(i,j,k)
            + beta * (f.cons.rhou(i,j,k) + dt * f.rhs.rhou(i,j,k));
        f.cons.rhov(i,j,k) = alpha * Q0.rhov(i,j,k)
            + beta * (f.cons.rhov(i,j,k) + dt * f.rhs.rhov(i,j,k));
        f.cons.rhow(i,j,k) = alpha * Q0.rhow(i,j,k)
            + beta * (f.cons.rhow(i,j,k) + dt * f.rhs.rhow(i,j,k));
        f.cons.rhoE(i,j,k) = alpha * Q0.rhoE(i,j,k)
            + beta * (f.cons.rhoE(i,j,k) + dt * f.rhs.rhoE(i,j,k));
    }}}
}
