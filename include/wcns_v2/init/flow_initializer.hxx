#pragma once

#include "flow_initializer.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// ============================================================================
// Interior range helper
// ============================================================================

inline void FlowInitializer::interior_range(const LocalBlock& lb,
                                              Int& i0, Int& i1,
                                              Int& j0, Int& j1,
                                              Int& k0, Int& k1) {
    Int ng = lb.grid.ng;
    i0 = ng;
    i1 = ng + lb.nci_core() - 1;
    j0 = ng;
    j1 = ng + lb.ncj_core() - 1;
    k0 = ng;
    k1 = ng + lb.nck_core() - 1;
}

// ============================================================================
// Main dispatch
// ============================================================================

inline void FlowInitializer::initialize(LocalBlock& lb, const Config& cfg) {
    if (cfg.init_type == "uniform") {
        init_uniform(lb, cfg);
    } else if (cfg.init_type == "poiseuille") {
        init_poiseuille(lb, cfg);
    } else {
        throw std::runtime_error("FlowInitializer: unknown init_type \"" +
                                 cfg.init_type + "\"");
    }

    // Synchronize conservative variables from primitive
    lb.field.prim_to_cons(cfg.gamma);
}

// ============================================================================
// Uniform free-stream
// ============================================================================

inline void FlowInitializer::init_uniform(LocalBlock& lb, const Config& cfg) {
    Int i0, i1, j0, j1, k0, k1;
    interior_range(lb, i0, i1, j0, j1, k0, k1);

    Real u_inf, v_inf, w_inf;
    cfg.free_stream_velocity(u_inf, v_inf, w_inf);

    // Non-dimensional free-stream pressure
    Real p_inf = 1.0 / (cfg.gamma * cfg.Mach * cfg.Mach);

    for (Int k = k0; k <= k1; ++k)
    for (Int j = j0; j <= j1; ++j)
    for (Int i = i0; i <= i1; ++i) {
        lb.field.prim.rho(i,j,k) = 1.0;
        lb.field.prim.u(i,j,k)   = u_inf;
        lb.field.prim.v(i,j,k)   = v_inf;
        lb.field.prim.w(i,j,k)   = w_inf;
        lb.field.prim.p(i,j,k)   = p_inf;
    }
}

// ============================================================================
// Poiseuille channel flow
// ============================================================================

inline void FlowInitializer::init_poiseuille(LocalBlock& lb, const Config& cfg) {
    Int i0, i1, j0, j1, k0, k1;
    interior_range(lb, i0, i1, j0, j1, k0, k1);

    // Determine flow direction from body force sign
    // body_force_x > 0 ⇒ flow in +x direction (parabola in y)
    // body_force_y > 0 ⇒ flow in +y direction (parabola in x or z)
    Real bfx = cfg.body_force_x;
    Real bfy = cfg.body_force_y;
    Real bfz = cfg.body_force_z;

    Real y_min = cfg.poiseuille_y_min;
    Real y_max = cfg.poiseuille_y_max;
    Real umax  = cfg.poiseuille_umax;

    // Non-dimensional free-stream pressure
    Real p_inf = 1.0 / (cfg.gamma * cfg.Mach * cfg.Mach);

    for (Int k = k0; k <= k1; ++k)
    for (Int j = j0; j <= j1; ++j)
    for (Int i = i0; i <= i1; ++i) {
        lb.field.prim.rho(i,j,k) = 1.0;
        lb.field.prim.p(i,j,k)   = p_inf;

        // Determine transverse coordinate and velocity component
        Real u_flow = 0.0, v_flow = 0.0, w_flow = 0.0;

        if (std::abs(bfx) > std::abs(bfy) && std::abs(bfx) > std::abs(bfz)) {
            // Flow along x, parabolic profile in y
            Real yc = lb.grid.cell_y(i,j,k);
            Real yn = (yc - y_min) / (y_max - y_min);  // [0, 1]
            u_flow = 4.0 * umax * yn * (1.0 - yn);
        } else if (std::abs(bfy) > std::abs(bfx) && std::abs(bfy) > std::abs(bfz)) {
            // Flow along y, parabolic profile in x
            Real xc = lb.grid.cell_x(i,j,k);
            // Use y_min/y_max as bounds for the transverse direction
            Real xn = (xc - y_min) / (y_max - y_min);  // [0, 1]
            v_flow = 4.0 * umax * xn * (1.0 - xn);
        } else if (std::abs(bfz) > 0) {
            // Flow along z, parabolic profile in x
            Real xc = lb.grid.cell_x(i,j,k);
            Real xn = (xc - y_min) / (y_max - y_min);  // [0, 1]
            w_flow = 4.0 * umax * xn * (1.0 - xn);
        } else {
            // Default: flow along x
            Real yc = lb.grid.cell_y(i,j,k);
            Real yn = (yc - y_min) / (y_max - y_min);
            u_flow = 4.0 * umax * yn * (1.0 - yn);
        }

        lb.field.prim.u(i,j,k) = u_flow;
        lb.field.prim.v(i,j,k) = v_flow;
        lb.field.prim.w(i,j,k) = w_flow;
    }
}
