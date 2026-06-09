#pragma once

#include "utils/types.h"

// Forward declarations
class LocalBlock;

/// @file inviscid_rhs.h
/// @brief Inviscid RHS computation — chains interpolation → Riemann → differentiation.
///
/// Computes the inviscid contribution to the right-hand side of the Euler equations:
///
///   RHS = -(∂F_inv/∂ξ + ∂G_inv/∂η + ∂H_inv/∂ζ)
///
/// where F_inv, G_inv, H_inv are the inviscid numerical fluxes at faces (ξ, η, ζ),
/// pre-computed by the Riemann solver and stored in Field::inv_xi/eta/zeta.
///
/// Differentiation uses the 6th-order centered difference on the face half-nodes.
/// Only interior cells (3..nci-4, 3..ncj-4, 3..nck-4) are updated — boundary
/// ghost cells are left untouched (they will be filled by BC / halo exchange).

class InviscidRHS {
public:
    /// Compute inviscid RHS contribution at interior cells.
    ///
    /// Reads face fluxes inv_xi/eta/zeta from lb.field (must be pre-computed
    /// by Riemann solver).  Writes to lb.field.rhs (conservative RHS).
    ///
    /// RHS = -(∂F_inv/∂ξ + ∂G_inv/∂η + ∂H_inv/∂ζ)
    /// Computed with 6th-order centered differences, interior cells only.
    static void compute(LocalBlock& lb);
};

#include "inviscid_rhs.hxx"
