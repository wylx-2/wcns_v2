#pragma once

#include "body_force.h"
#include "wcns_v2/parallel/local_block.h"
#include <string>

// ============================================================================
// BodyForce::add_to_rhs
// ============================================================================

inline void BodyForce::add_to_rhs(LocalBlock& lb, const Config& cfg) {
    // No-op for "none" mode
    if (cfg.body_force_type == "none") return;

    Real fx, fy, fz;
    if (cfg.body_force_type == "constant") {
        fx = cfg.body_force_x;
        fy = cfg.body_force_y;
        fz = cfg.body_force_z;
    } else {
        // Reserved for future types: "gravity", "centrifugal", "coriolis"
        throw std::runtime_error("BodyForce: unknown body_force_type \""
                                 + cfg.body_force_type + "\"");
    }

    auto& f = lb.field;

    const Int nci = f.ni();
    const Int ncj = f.nj();
    const Int nck = f.nk();

    // Interior cells only (same range as RHS computation)
    const Int i0 = 3, i1 = nci - 4;
    const Int j0 = 3, j1 = ncj - 4;
    const Int k0 = 3, k1 = nck - 4;

    // Physical source: S = [0, ρ·f_x, ρ·f_y, ρ·f_z, ρ·f·v]^T
    // Note: rhs stores ∂Q/∂t (physical time derivative), so the source term
    // is added directly WITHOUT the 1/J factor.  (The inviscid and viscous
    // RHS use 1/J because their flux divergences are in computational space.)
    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        Real rho  = f.prim.rho(i,j,k);
        Real u    = f.prim.u(i,j,k);
        Real v    = f.prim.v(i,j,k);
        Real w    = f.prim.w(i,j,k);

        // S_rho  = 0
        f.rhs.rhou(i,j,k) += rho * fx;
        f.rhs.rhov(i,j,k) += rho * fy;
        f.rhs.rhow(i,j,k) += rho * fz;
        if (cfg.body_force_energy) {
            f.rhs.rhoE(i,j,k) += rho * (fx * u + fy * v + fz * w);
        }
    }}}
}
