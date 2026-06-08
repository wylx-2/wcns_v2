#pragma once

#include "wcns_v2/utils/types.h"

// Forward declarations
class LocalBlock;
struct Config;

/// @file lu_sgs.h
/// @brief LU-SGS (Lower-Upper Symmetric Gauss-Seidel) implicit time integration.
///
/// Solves the implicit Euler system:
///     [Ω/Δt - ∂R/∂Q] ΔQ = R(Q^n)
///
/// using the approximate LU factorization:
///     (D + L) D^{-1} (D + U) ΔQ = RHS
///
/// where D, L, U are built from scalar spectral radii of the inviscid flux
/// Jacobians in each curvilinear direction.  Viscous spectral radii are
/// included when applicable.
///
/// The algorithm consists of:
///   1. Precompute cell-center spectral radii σ_ξ, σ_η, σ_ζ
///   2. Forward sweep  (i↑, j↑, k↑):  solve (D + L) ΔQ* = Ω·RHS
///   3. Backward sweep (i↓, j↓, k↓):  solve (D + U) ΔQ  = D ΔQ*
///   4. Update: Q^{n+1} = Q^n + ΔQ
///
/// Only interior cells [ng .. nci-ng-1]×[ncj-ng-1]×[nck-ng-1] are updated.
/// Ghost/boundary cells are updated by subsequent BC application and
/// halo exchange.

class LuSgs {
public:
    /// Advance conservative variables for one time step using LU-SGS.
    ///
    /// Preconditions:
    ///   - lb.field.rhs  contains dQ/dt at the current state  (as computed
    ///                    by InviscidRHS/ViscidRHS/BodyForce pipelines)
    ///   - lb.field.Q0   contains the conserved state Q^n at the start of
    ///                    the time step
    ///
    /// Postconditions:
    ///   - lb.field.cons interior cells contain Q^{n+1} = Q^n + ΔQ
    ///   - Ghost/boundary cells in lb.field.cons are UNCHANGED
    ///
    /// @param lb   Local block (reads rhs, Q0, grid metrics, prim; writes cons)
    /// @param cfg  Solver configuration
    /// @param dt   Time step size
    static void advance(LocalBlock& lb, const Config& cfg, Real dt);
};

#include "lu_sgs.hxx"
