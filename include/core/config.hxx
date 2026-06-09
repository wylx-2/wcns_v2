#pragma once

#include "config.h"
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline void Config::finalize() {
    // Speed of sound: c_ref = sqrt(gamma * R_gas * T_ref)   [dimensional]
    Real c_ref = std::sqrt(gamma * R_gas * T_ref);

    // Reference velocity: U_ref = Mach * c_ref
    U_ref = Mach * c_ref;

    // Reference pressure (based on dynamic pressure scaling)
    p_ref = rho_ref * U_ref * U_ref;

    // Reference viscosity from Reynolds number
    mu_ref = rho_ref * U_ref * L_ref / Re;

    // Sutherland constant (non-dimensional)
    sutherland_S_nd = sutherland_S / T_ref;
}

inline Real Config::eos_factor() const {
    // p* = rho* * T* / (gamma * Mach^2)
    return 1.0 / (gamma * Mach * Mach);
}

inline void Config::free_stream_velocity(Real& u, Real& v, Real& w) const {
    Real aoa_rad  = AoA  * M_PI / 180.0;
    Real beta_rad = beta * M_PI / 180.0;

    Real ca = std::cos(aoa_rad);
    Real sa = std::sin(aoa_rad);
    Real cb = std::cos(beta_rad);
    Real sb = std::sin(beta_rad);

    // Non-dimensional: velocity magnitude = 1.0 (normalized by U_ref)
    u = ca * cb;
    v = sb;
    w = sa * cb;
}
