#pragma once

#include "wcns_v2/core/config.h"
#include "wcns_v2/parallel/local_block.h"
#include <vector>

/// @file time_step.h
/// @brief Compute the global time step from the CFL condition.
///
/// Supports two modes:
///   - CFL-based:  Δt = CFL * min_cells( Ω / Σ_d Λ_d )
///   - Fixed:      Δt = cfg.fixed_dt  (when fixed_dt > 0)
///
/// Spectral radii include convective and viscous contributions in
/// each curvilinear direction (ξ, η, ζ).

class TimeStep {
public:
    /// Compute the global minimum time step across all local blocks.
    ///
    /// If cfg.fixed_dt > 0, returns the fixed value directly.
    /// Otherwise loops over interior cells of all blocks, evaluates
    /// the convective + viscous spectral radii, and returns
    /// Δt = CFL * min( Ω / (Λ_c + Λ_v) ).
    ///
    /// @param blocks   All local blocks owned by this process
    /// @param cfg      Solver configuration (cfl, fixed_dt, gamma, Re, Pr)
    /// @return Global time step (non-dimensional)
    static Real compute(const std::vector<LocalBlock>& blocks,
                        const Config& cfg);

private:
    /// Compute the local time step for one block.
    /// Returns the minimum over interior cells.
    static Real compute_block(const LocalBlock& lb, const Config& cfg);

    /// Compute the convective spectral radius sum ΣΛ_c for one cell.
    /// Λ_c^d = |U_contra^d| + c * S^d
    /// where U_contra^d = u·met_d_x + v·met_d_y + w·met_d_z
    ///       S^d = sqrt(met_d_x^2 + met_d_y^2 + met_d_z^2)
    static Real convective_sum(const Real met[3], Real u, Real v, Real w,
                               Real c_sound);

    /// Compute the viscous spectral radius sum ΣΛ_v for one cell.
    /// Λ_v^d = β * (μ/ρ) / Pr * (S^d)^2 / Ω
    /// where β = max(4/3, gamma), Ω = 1/J.
    static Real viscous_sum(const Real met_xi[3], const Real met_eta[3],
                            const Real met_zeta[3], Real omega,
                            Real mu_rho, Real beta, Real Pr);
};

#include "time_step.hxx"
