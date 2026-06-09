#pragma once

#include "residual.h"
#include "parallel/local_block.h"
#include "parallel/parallel_manager.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

// ============================================================================
// Residual::compute
// ============================================================================

inline ResidualNorms Residual::compute(const std::vector<LocalBlock>& blocks) {
    Real local_sq[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    Real local_N = 0.0;

    for (const auto& lb : blocks) {
        const auto& f = lb.field;
        const Int nci = f.ni();
        const Int ncj = f.nj();
        const Int nck = f.nk();

        // Interior cells only — same range as RHS computation
        const Int i0 = 3, i1 = nci - 4;
        const Int j0 = 3, j1 = ncj - 4;
        const Int k0 = 3, k1 = nck - 4;

        for (Int k = k0; k <= k1; ++k) {
        for (Int j = j0; j <= j1; ++j) {
        for (Int i = i0; i <= i1; ++i) {
            Real r;
            r = f.rhs.rho(i,j,k);   local_sq[0] += r * r;
            r = f.rhs.rhou(i,j,k);  local_sq[1] += r * r;
            r = f.rhs.rhov(i,j,k);  local_sq[2] += r * r;
            r = f.rhs.rhow(i,j,k);  local_sq[3] += r * r;
            r = f.rhs.rhoE(i,j,k);  local_sq[4] += r * r;
        }}}

        Int nci_core = i1 - i0 + 1;
        Int ncj_core = j1 - j0 + 1;
        Int nck_core = k1 - k0 + 1;
        local_N += static_cast<Real>(nci_core) *
                   static_cast<Real>(ncj_core) *
                   static_cast<Real>(nck_core);
    }

    // Global reduction: each component's squared sum and cell count
    Real global_N = ParallelManager::global_sum(local_N);

    ResidualNorms res;
    if (global_N > 0.0) {
        res.rho  = std::sqrt(ParallelManager::global_sum(local_sq[0]) / global_N);
        res.rhou = std::sqrt(ParallelManager::global_sum(local_sq[1]) / global_N);
        res.rhov = std::sqrt(ParallelManager::global_sum(local_sq[2]) / global_N);
        res.rhow = std::sqrt(ParallelManager::global_sum(local_sq[3]) / global_N);
        res.rhoE = std::sqrt(ParallelManager::global_sum(local_sq[4]) / global_N);
    }

    return res;
}

// ============================================================================
// Residual::monitor
// ============================================================================

inline FlowMonitor Residual::monitor(const std::vector<LocalBlock>& blocks,
                                      const Config& cfg) {
    // Pre-compute derived constants
    const Real gm2   = cfg.gamma * cfg.Mach * cfg.Mach;  // gamma * Mach²
    const Real gamma = cfg.gamma;

    FlowMonitor mon;
    // Initialise with identity elements for min/max operations
    mon.rho_min  =  std::numeric_limits<Real>::max();
    mon.rho_max  = -std::numeric_limits<Real>::max();
    mon.T_min    =  std::numeric_limits<Real>::max();
    mon.T_max    = -std::numeric_limits<Real>::max();
    mon.p_min    =  std::numeric_limits<Real>::max();
    mon.p_max    = -std::numeric_limits<Real>::max();
    mon.vel_min  =  std::numeric_limits<Real>::max();
    mon.vel_max  = -std::numeric_limits<Real>::max();
    mon.Mach_max = -std::numeric_limits<Real>::max();

    for (const auto& lb : blocks) {
        const auto& f = lb.field;
        const Int nci = f.ni();
        const Int ncj = f.nj();
        const Int nck = f.nk();

        const Int i0 = 3, i1 = nci - 4;
        const Int j0 = 3, j1 = ncj - 4;
        const Int k0 = 3, k1 = nck - 4;

        for (Int k = k0; k <= k1; ++k) {
        for (Int j = j0; j <= j1; ++j) {
        for (Int i = i0; i <= i1; ++i) {
            Real rho = f.prim.rho(i,j,k);
            Real p   = f.prim.p(i,j,k);
            Real u   = f.prim.u(i,j,k);
            Real v   = f.prim.v(i,j,k);
            Real w   = f.prim.w(i,j,k);

            // ---- NaN / Inf detection ----
            if (std::isnan(rho) || std::isnan(p) ||
                std::isnan(u) || std::isnan(v) || std::isnan(w)) {
                mon.has_nan = true;
            }
            if (std::isinf(rho) || std::isinf(p) ||
                std::isinf(u) || std::isinf(v) || std::isinf(w)) {
                mon.has_inf = true;
            }

            // ---- Density ----
            mon.rho_min = std::min(mon.rho_min, rho);
            mon.rho_max = std::max(mon.rho_max, rho);

            // ---- Pressure ----
            mon.p_min = std::min(mon.p_min, p);
            mon.p_max = std::max(mon.p_max, p);

            // ---- Temperature  (T* = p* · γ · Mach² / ρ*) ----
            Real T = p * gm2 / rho;
            mon.T_min = std::min(mon.T_min, T);
            mon.T_max = std::max(mon.T_max, T);

            // ---- Velocity magnitude ----
            Real vel = std::sqrt(u * u + v * v + w * w);
            mon.vel_min = std::min(mon.vel_min, vel);
            mon.vel_max = std::max(mon.vel_max, vel);

            // ---- Local Mach number  (M = |U| / c,  c = √(γ·p/ρ)) ----
            Real c_sound = std::sqrt(gamma * p / rho);
            Real Mach_local = vel / c_sound;
            mon.Mach_max = std::max(mon.Mach_max, Mach_local);
        }}}
    }

    // ---- Global reductions across MPI ranks ----
    mon.rho_min  = ParallelManager::global_min(mon.rho_min);
    mon.rho_max  = ParallelManager::global_max(mon.rho_max);
    mon.T_min    = ParallelManager::global_min(mon.T_min);
    mon.T_max    = ParallelManager::global_max(mon.T_max);
    mon.p_min    = ParallelManager::global_min(mon.p_min);
    mon.p_max    = ParallelManager::global_max(mon.p_max);
    mon.vel_min  = ParallelManager::global_min(mon.vel_min);
    mon.vel_max  = ParallelManager::global_max(mon.vel_max);
    mon.Mach_max = ParallelManager::global_max(mon.Mach_max);

    // NaN / Inf flags: OR across all ranks (max of 0/1 flags)
    Real nan_flag = mon.has_nan ? 1.0 : 0.0;
    Real inf_flag = mon.has_inf ? 1.0 : 0.0;
    mon.has_nan = (ParallelManager::global_max(nan_flag) > 0.5);
    mon.has_inf = (ParallelManager::global_max(inf_flag) > 0.5);

    return mon;
}

// ============================================================================
// Residual::write_header
// ============================================================================

inline void Residual::write_header(std::ostream& os) {
    os << "# iter       dt             "
       << "res_rho        res_rhou       res_rhov       res_rhow       res_rhoE        "
       << "rho_min        rho_max        T_min          T_max          "
       << "p_min          p_max          vel_max        Mach_max\n";
}

// ============================================================================
// Residual::log
// ============================================================================

inline void Residual::log(std::ostream& os, Int iter, Real dt,
                           const ResidualNorms& res, const FlowMonitor& mon) {
    constexpr int W  = 16;  // field width
    constexpr int PR = 6;   // decimal precision

    os << std::scientific << std::setprecision(PR)
       << std::setw(6)  << iter
       << std::setw(W)  << dt
       << std::setw(W)  << res.rho
       << std::setw(W)  << res.rhou
       << std::setw(W)  << res.rhov
       << std::setw(W)  << res.rhow
       << std::setw(W)  << res.rhoE
       << std::setw(W)  << mon.rho_min
       << std::setw(W)  << mon.rho_max
       << std::setw(W)  << mon.T_min
       << std::setw(W)  << mon.T_max
       << std::setw(W)  << mon.p_min
       << std::setw(W)  << mon.p_max
       << std::setw(W)  << mon.vel_max
       << std::setw(W)  << mon.Mach_max
       << std::endl;
}

// ============================================================================
// Residual::converged  /  Residual::diverged
// ============================================================================

inline bool Residual::converged(const ResidualNorms& res, Real tol) {
    return res.rho  < tol &&
           res.rhou < tol &&
           res.rhov < tol &&
           res.rhow < tol &&
           res.rhoE < tol;
}

inline bool Residual::diverged(const FlowMonitor& mon) {
    return mon.has_nan || mon.has_inf;
}
