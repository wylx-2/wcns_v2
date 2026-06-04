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
    } else if (cfg.init_type == "riemann_2d") {
        init_riemann_2d(lb, cfg);
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

// ============================================================================
// 2D Riemann problem — 4-quadrant discontinuous initial states
// ============================================================================
//
// Classic configurations from Lax & Liu (1998), "Solution of Two-Dimensional
// Riemann Problems of Gas Dynamics by Positive Schemes".
//
// Domain: [0,1] × [0,1], split into 4 quadrants at (x_split, y_split).
//
// Configuration 3 (riemann_config=3):
//   A single shock wave in each quadrant, all pointing inward.  Produces
//   a complex interaction with a Mach stem at late time.
//
//         NW           NE
//   ρ=1.5           ρ=0.5323
//   u=0, v=0        u=1.206, v=0
//   p=1.5           p=0.3
//         SW           SE
//   ρ=0.138         ρ=0.5323
//   u=1.206, v=1.206 u=0, v=1.206
//   p=0.029         p=0.3
//
// Configuration 12 (riemann_config=12):
//   Four initial planar shock waves interacting.
//
//         NW           NE
//   ρ=1.0           ρ=2.0
//   u=0.0, v=0.3    u=-0.3, v=0.0
//   p=1.0           p=1.0
//         SW           SE
//   ρ=1.0625        ρ=0.5197
//   u=0.0, v=0.8145 u=0.0, v=0.3
//   p=0.4           p=0.4

inline void FlowInitializer::init_riemann_2d(LocalBlock& lb, const Config& cfg) {
    Int i0, i1, j0, j1, k0, k1;
    interior_range(lb, i0, i1, j0, j1, k0, k1);

    Real xs = cfg.riemann_x_split;
    Real ys = cfg.riemann_y_split;

    // Quadrant primitive states: {rho, u, v, p}
    // Order: NW, NE, SW, SE
    Real q[4][4];

    if (cfg.riemann_config == 12) {
        // Configuration 12
        // NW (x<xs, y>ys)
        q[0][0]=1.0;     q[0][1]=0.0;    q[0][2]=0.3;    q[0][3]=1.0;
        // NE (x>xs, y>ys)
        q[1][0]=2.0;     q[1][1]=-0.3;   q[1][2]=0.0;    q[1][3]=1.0;
        // SW (x<xs, y<ys)
        q[2][0]=1.0625;  q[2][1]=0.0;    q[2][2]=0.8145; q[2][3]=0.4;
        // SE (x>xs, y<ys)
        q[3][0]=0.5197;  q[3][1]=0.0;    q[3][2]=0.3;    q[3][3]=0.4;
    } else {
        // Configuration 3 (default)
        // NW (x<xs, y>ys)
        q[0][0]=1.5;     q[0][1]=0.0;    q[0][2]=0.0;    q[0][3]=1.5;
        // NE (x>xs, y>ys)
        q[1][0]=0.5323;  q[1][1]=1.206;  q[1][2]=0.0;    q[1][3]=0.3;
        // SW (x<xs, y<ys)
        q[2][0]=0.138;   q[2][1]=1.206;  q[2][2]=1.206;  q[2][3]=0.029;
        // SE (x>xs, y<ys)
        q[3][0]=0.5323;  q[3][1]=0.0;    q[3][2]=1.206;  q[3][3]=0.3;
    }

    for (Int k = k0; k <= k1; ++k)
    for (Int j = j0; j <= j1; ++j)
    for (Int i = i0; i <= i1; ++i) {
        Real xc = lb.grid.cell_x(i,j,k);
        Real yc = lb.grid.cell_y(i,j,k);

        int quad = (xc > xs ? 1 : 0)  // NE/SE vs NW/SW
                 | (yc > ys ? 0 : 2); // NW/NE vs SW/SE
        // quad: 0=NW, 1=NE, 2=SW, 3=SE

        lb.field.prim.rho(i,j,k) = q[quad][0];
        lb.field.prim.u(i,j,k)   = q[quad][1];
        lb.field.prim.v(i,j,k)   = q[quad][2];
        lb.field.prim.w(i,j,k)   = 0.0;  // 2D: w = 0
        lb.field.prim.p(i,j,k)   = q[quad][3];
    }
}
