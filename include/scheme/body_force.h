#pragma once

#include "utils/types.h"

// Forward declarations
class LocalBlock;
struct Config;

/// @file body_force.h
/// @brief Volumetric body force source term for the Navier-Stokes equations.
///
/// The physical source term S = [0, ρ·f_x, ρ·f_y, ρ·f_z, ρ·f·v]^T is
/// converted to computational space by dividing by the Jacobian:
///
///   rhs += J⁻¹ · S
///
/// This ensures the time integrator's "multiply by J" step correctly
/// recovers the physical source:
///
///   Q_new = α·Q0 + β·(Q_old + dt·J·rhs)
///         = α·Q0 + β·(Q_old + dt·J·rhs_flux + dt·S_physical)
///
/// Only interior cells [3..nci-4]×[3..ncj-4]×[3..nck-4] are updated.
/// Supported modes (via cfg.body_force_type):
///   - "none":     no body force (no-op)
///   - "constant": uniform acceleration vector f = (body_force_x, body_force_y, body_force_z)
///   - (reserved) "gravity", "centrifugal", "coriolis"

class BodyForce {
public:
    /// Compute body force source term and add to RHS.
    ///
    /// @param lb  Local block (reads prim.{rho,u,v,w}, grid.jacobian;
    ///             adds to field.rhs)
    /// @param cfg Configuration (body_force_type, body_force_{x,y,z})
    ///
    /// If body_force_type == "none", returns immediately (no-op).
    /// Only interior cells are updated.
    static void add_to_rhs(LocalBlock& lb, const Config& cfg);
};

#include "body_force.hxx"
