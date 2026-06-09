#pragma once

#include "time_step.h"
#include <algorithm>
#include <cmath>

// ============================================================================
// Public interface
// ============================================================================

inline Real TimeStep::compute(const std::vector<LocalBlock>& blocks,
                               const Config& cfg) {
    // Fixed time step override
    if (cfg.fixed_dt > Real(0)) {
        return cfg.fixed_dt;
    }

    Real dt_min = Real(1e30);

    for (const auto& lb : blocks) {
        Real dt_loc = compute_block(lb, cfg);
        dt_min = std::min(dt_min, dt_loc);
    }

    // MPI: global minimum across all processes
    // (will be added when integrated with ParallelManager)
    return dt_min;
}

// ============================================================================
// Per-block computation
// ============================================================================

inline Real TimeStep::compute_block(const LocalBlock& lb, const Config& cfg) {
    Int ng = cfg.ng;

    // Interior cell range
    Int i0 = ng, i1 = ng + lb.nci_core() - 1;
    Int j0 = ng, j1 = ng + lb.ncj_core() - 1;
    Int k0 = ng, k1 = ng + lb.nck_core() - 1;

    Real cfl   = cfg.cfl;
    Real gamma = cfg.gamma;
    Real Pr    = cfg.Prandtl;
    Real Re    = cfg.Re;
    Real beta  = std::max(Real(4.0 / 3.0), gamma);

    // Speed-of-sound factor: c = sqrt(gamma * p / rho)
    // Temperature factor: T = gamma * Mach^2 * p / rho
    // (non-dimensional EOS: p = rho * T / (gamma * Mach^2))

    Real dt_min = Real(1e30);

    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        // Primitive variables
        Real rho = lb.field.prim.rho(i, j, k);
        Real u   = lb.field.prim.u(i, j, k);
        Real v   = lb.field.prim.v(i, j, k);
        Real w   = lb.field.prim.w(i, j, k);
        Real p   = lb.field.prim.p(i, j, k);

        // Speed of sound
        Real c_sound = std::sqrt(gamma * p / rho);

        // Dynamic viscosity (constant, non-dimensional: μ_nd = 1/Re)
        // TODO: Sutherland formula for temperature-dependent viscosity
        Real mu = Real(1.0) / Re;
        Real mu_rho = mu / rho;

        // Cell volume in physical space: Ω = |J|.
        // SCMM Jacobian J IS the cell volume, since Δξ=Δη=Δζ=1.
        // abs() handles reversed coordinate handedness in some zones.
        Real omega = std::abs(lb.grid.jacobian(i, j, k));

        // -- ξ direction --
        Real met_xi[3] = {lb.grid.met_xi_x(i, j, k),
                          lb.grid.met_xi_y(i, j, k),
                          lb.grid.met_xi_z(i, j, k)};
        Real lambda_c_xi = convective_sum(met_xi, u, v, w, c_sound);

        // -- η direction --
        Real met_eta[3] = {lb.grid.met_eta_x(i, j, k),
                           lb.grid.met_eta_y(i, j, k),
                           lb.grid.met_eta_z(i, j, k)};
        Real lambda_c_eta = convective_sum(met_eta, u, v, w, c_sound);

        // -- ζ direction --
        Real met_zeta[3] = {lb.grid.met_zeta_x(i, j, k),
                            lb.grid.met_zeta_y(i, j, k),
                            lb.grid.met_zeta_z(i, j, k)};
        Real lambda_c_zeta = convective_sum(met_zeta, u, v, w, c_sound);

        Real lambda_c = lambda_c_xi + lambda_c_eta + lambda_c_zeta;

        Real lambda_v = viscous_sum(met_xi, met_eta, met_zeta,
                                     omega, mu_rho, beta, Pr);

        // Combined:  Δt = CFL * Ω / (Λ_c + Λ_v)
        Real dt_cell = cfl * omega / (lambda_c + lambda_v);

        dt_min = std::min(dt_min, dt_cell);
    }}}

    return dt_min;
}

// ============================================================================
// Spectral radii helpers
// ============================================================================

inline Real TimeStep::convective_sum(const Real met[3],
                                      Real u, Real v, Real w,
                                      Real c_sound) {
    // Contravariant velocity × Jacobian:     U^d = u·met_d_x + v·met_d_y + w·met_d_z
    // Face area × Jacobian (magnitude):      S^d = sqrt(met_d_x² + met_d_y² + met_d_z²)
    // Convective spectral radius:            Λ_c^d = |U^d| + c·S^d

    Real u_contra = u * met[0] + v * met[1] + w * met[2];
    Real s_face   = std::sqrt(met[0] * met[0] + met[1] * met[1] + met[2] * met[2]);

    return std::abs(u_contra) + c_sound * s_face;
}

inline Real TimeStep::viscous_sum(const Real met_xi[3],
                                   const Real met_eta[3],
                                   const Real met_zeta[3],
                                   Real omega,
                                   Real mu_rho, Real beta, Real Pr) {
    // Viscous spectral radius in direction d:
    //   Λ_v^d = β * (μ/ρ) / Pr * (S^d)² / Ω
    //
    // where β = max(4/3, γ) is the maximum eigenvalue of the viscous
    // flux Jacobian, and S^d is the face area (×J).

    auto s2 = [](const Real m[3]) {
        return m[0] * m[0] + m[1] * m[1] + m[2] * m[2];
    };

    Real factor = beta * mu_rho / (Pr * omega);

    return factor * (s2(met_xi) + s2(met_eta) + s2(met_zeta));
}
