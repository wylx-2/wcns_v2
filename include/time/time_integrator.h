#pragma once

#include "utils/types.h"

// Forward declarations
class LocalBlock;
struct Config;
struct ConservativeVars;

/// @file time_integrator.h
/// @brief Explicit Runge-Kutta time integrators for semi-discrete NS equations.
///
/// Supports RK3-TVD (current), RK4 (reserved), and LU-SGS (reserved).
///
/// The RHS computed by InviscidRHS/ViscidRHS stores:
///     rhs = -(∂F̂/∂ξ + ∂Ĝ/∂η + ∂Ĥ/∂ζ) = J⁻¹·∂Q/∂t
///
/// The time advance must convert this to ∂Q/∂t by multiplying by J (Jacobian):
///     Q_new = α·Q0 + β·(Q_old + Δt · J · rhs)
///
/// Only interior cells [3..nci-4]×[3..ncj-4]×[3..nck-4] are updated.
/// Ghost and boundary cells are updated by subsequent BC application and
/// halo exchange.

class TimeIntegrator {
public:
    // ========================================================================
    // Scheme parameters
    // ========================================================================

    /// Return the number of stages for the configured time scheme.
    /// "rk3-tvd" → 3, "rk4" → 4, "lu-sgs" → 1
    static int n_stages(const Config& cfg);

    /// Return alpha and beta coefficients for a given stage (0-indexed).
    /// For "rk3-tvd":
    ///   stage=0: α=0,   β=1
    ///   stage=1: α=3/4, β=1/4
    ///   stage=2: α=1/3, β=2/3
    static void rk_coeffs(const Config& cfg, int stage, Real& alpha, Real& beta);

    // ========================================================================
    // Stage advance
    // ========================================================================

    /// Advance conservative variables for one RK stage.
    ///
    /// Q_new(i,j,k) = α·Q0(i,j,k) + β·(Q_old(i,j,k) + Δt·J(i,j,k)·rhs(i,j,k))
    ///
    /// @param lb     Local block (reads cons, rhs, grid.jacobian; writes cons)
    /// @param cfg    Configuration (time scheme)
    /// @param dt     Time step size
    /// @param stage  0-indexed RK stage number
    /// @param Q0     Conserved state at the beginning of the time step
    ///
    /// Prerequisites:
    ///   - RHS (inviscid + viscous) must be pre-computed and stored in lb.field.rhs
    ///   - Q0 must contain the conserved state at the start of the time step
    ///   - lb.field.cons must be the current stage's conserved state
    ///
    /// Post-condition:
    ///   - lb.field.cons interior cells contain the updated conserved state
    ///   - Ghost/boundary cells in lb.field.cons are UNCHANGED (caller must
    ///     subsequently convert to primitives, apply BC, and exchange halos)
    static void advance_stage(LocalBlock& lb, const Config& cfg,
                              Real dt, int stage,
                              const ConservativeVars& Q0);
};

#include "time_integrator.hxx"
