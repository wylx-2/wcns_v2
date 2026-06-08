#pragma once

#include "wcns_v2/utils/types.h"
#include "wcns_v2/core/config.h"
#include <ostream>
#include <vector>

/// @file residual.h
/// @brief Residual computation and flow-field monitoring for convergence control.
///
/// Residual computes the L2 norm of the RHS (all 5 conservative components)
/// over interior cells and tracks min/max of key physical quantities (density,
/// temperature, pressure, velocity, Mach number).  Results are globally reduced
/// across MPI ranks so every rank sees the same values.

// ============================================================================
// Forward declarations
// ============================================================================

class LocalBlock;

// ============================================================================
// Data structures
// ============================================================================

/// L2 residual norms for all 5 conservative variables.
///
/// Computed from the cell-center RHS array over interior cells
/// [3..nci-4]×[3..ncj-4]×[3..nck-4] (same range as RHS computation):
///   res[c] = sqrt( Σ RHS_c² / N_interior )
struct ResidualNorms {
    Real rho  = 0.0;  ///< ||RHS_rho||_L2  (continuity)
    Real rhou = 0.0;  ///< ||RHS_rhou||_L2 (x-momentum)
    Real rhov = 0.0;  ///< ||RHS_rhov||_L2 (y-momentum)
    Real rhow = 0.0;  ///< ||RHS_rhow||_L2 (z-momentum)
    Real rhoE = 0.0;  ///< ||RHS_rhoE||_L2 (energy)
};

/// Global min/max of monitored physical quantities.
///
/// Collected over interior cells across all blocks and MPI ranks.
/// Temperature is derived from the non-dimensional EOS:
///   T* = p* * gamma * Mach² / rho*
struct FlowMonitor {
    Real rho_min = 0.0, rho_max = 0.0;  ///< Density range
    Real T_min   = 0.0, T_max   = 0.0;  ///< Temperature range
    Real p_min   = 0.0, p_max   = 0.0;  ///< Pressure range
    Real vel_min = 0.0, vel_max = 0.0;  ///< Velocity magnitude range
    Real Mach_max = 0.0;                 ///< Maximum local Mach number

    bool has_nan = false;  ///< True if any monitored quantity is NaN
    bool has_inf = false;  ///< True if any monitored quantity is Inf

    /// Reserved for future extensions (turbulence k/ω, species mass fractions, …)
    // std::map<std::string, Real> extras;
};

// ============================================================================
// Residual class — purely static methods
// ============================================================================

/// Convergence and flow-field diagnostics.
///
/// All methods that return a struct (compute, monitor) perform MPI global
/// reductions internally and return the same result on every rank.
class Residual {
public:
    /// Compute L2 residual norms for all 5 conservative components.
    ///
    /// Iterates over interior cells [3..nci-4]×[3..ncj-4]×[3..nck-4] of
    /// every local block, accumulates Σ(RHS²) per component, then globally
    /// reduces across MPI ranks.  Returns sqrt(global_sum / global_N).
    static ResidualNorms compute(const std::vector<LocalBlock>& blocks);

    /// Collect global min/max of physical quantities over interior cells.
    ///
    /// Monitors: density, temperature, pressure, velocity magnitude,
    /// and local Mach number.  NaN/Inf flags are OR-ed across all ranks
    /// so every process sees a divergence condition.
    static FlowMonitor monitor(const std::vector<LocalBlock>& blocks,
                               const Config& cfg);

    /// Write the header line for residual.dat (rank 0 only).
    static void write_header(std::ostream& os);

    /// Write one data row to residual.dat (rank 0 only).
    static void log(std::ostream& os, Int iter, Real dt,
                    const ResidualNorms& res, const FlowMonitor& mon);

    /// Convergence test: density residual below the given tolerance.
    static bool converged(const ResidualNorms& res, Real tol);

    /// Divergence test: any monitored quantity is NaN or Inf.
    static bool diverged(const FlowMonitor& mon);
};

#include "residual.hxx"
