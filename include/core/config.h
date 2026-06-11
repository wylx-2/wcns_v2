#pragma once

#include "utils/types.h"
#include <string>
#include <cmath>

/// @file config.h
/// @brief Physical and computation control parameters for the WCNS solver.
///
/// All parameters are stored in non-dimensional form where applicable.
/// Reference values are used to non-dimensionalize inputs (grid, initial/boundary
/// conditions) and to re-dimensionalize outputs for post-processing.
///
/// Non-dimensionalization convention:
///   rho* = rho / rho_ref      u* = u / U_ref        T* = T / T_ref
///   p*   = p   / (rho_ref * U_ref^2)
///   L*   = L   / L_ref        t* = t * U_ref / L_ref
///   mu*  = mu  / mu_ref        (mu_ref = rho_ref * U_ref * L_ref / Re)
///
/// Equation of state (non-dimensional):
///   p* = rho* * T* / (gamma * Mach^2)
///
/// Derived reference quantities are computed by finalize().

struct Config {

    // =========================================================================
    // Physical parameters
    // =========================================================================

    Real gamma   = 1.4;   ///< Specific heat ratio
    Real Prandtl = 0.72;  ///< Prandtl number
    Real Re      = 1.0e6; ///< Reynolds number  (Re = rho_ref * U_ref * L_ref / mu_ref)
    Real Mach    = 0.5;   ///< Free-stream Mach number

    Real AoA  = 0.0;      ///< Angle of attack  [degrees]
    Real beta = 0.0;      ///< Sideslip angle    [degrees]

    // Free-stream velocity direction in Cartesian coordinates (computed from AoA, beta)
    // u_inf = (cos(AoA)*cos(beta), sin(beta), sin(AoA)*cos(beta))

    // =========================================================================
    // Reference values (dimensional) — for non-dimensionalization
    // =========================================================================

    Real L_ref   = 1.0;   ///< Reference length (e.g. chord, body length)
    Real rho_ref = 1.0;   ///< Reference density (e.g. free-stream density)
    Real T_ref   = 1.0;   ///< Reference temperature (e.g. free-stream T)

    // =========================================================================
    // Derived reference quantities (set by finalize())
    // =========================================================================

    Real U_ref  = 0.0;    ///< Reference velocity = Mach * sqrt(gamma * R_gas * T_ref)
    Real p_ref  = 0.0;    ///< Reference pressure = rho_ref * U_ref^2
    Real mu_ref = 0.0;    ///< Reference viscosity = rho_ref * U_ref * L_ref / Re
    Real R_gas  = 287.0;  ///< Dimensional gas constant [J/(kg·K)] (air)
                          ///< Used only to derive U_ref; not needed in ND solver

    // =========================================================================
    // Computation control
    // =========================================================================

    Real    cfl           = 0.5;        ///< CFL number
    Real    fixed_dt      = 0.0;        ///< Fixed time step (0 = compute from CFL)
    Int     max_iter      = 1000;       ///< Maximum number of iterations
    Real    max_time      = 0.0;        ///< Maximum physical time (0 = disabled, step-based only)
    Int     output_freq   = 100;        ///< Write solution every N iterations
    Int     restart_freq  = 500;        ///< Write restart files every N iterations
    Int     residual_freq = 1;          ///< Compute residuals & monitor every N iterations
    std::string restart_file = "";      ///< Restart file basename for resuming ("" = start from scratch)
    Real    converge_tol  = 1.0e-6;     ///< Convergence tolerance (L2 residual norm)
    std::string time_scheme = "rk3-tvd";///< Time integration: "rk3-tvd", "rk4", "lu-sgs"
    Real    lu_sgs_kappa  = 4.0;        ///< LU-SGS over-relaxation factor (>1 for stability)
    Int     ng            = 3;          ///< Ghost layer count (authoritative source)

    // =========================================================================
    // Scheme options
    // =========================================================================

    /// Variable space for WCNS interpolation
    /// "conservative"  — interpolate conservative vars component-wise (default)
    /// "characteristic" — project to characteristic space, interpolate, project back (reserved)
    std::string interp_vars = "conservative";

    /// WCNS interpolation scheme type
    /// "weno_js"       — WENO-JS 5th-order nonlinear interpolation (default)
    /// "mdcd_linear"   — MDCD linear interpolation with tunable diss/disp
    /// "mdcd_hybrid"   — MDCD-WENO hybrid with discontinuity detector
    std::string interp_type = "weno_js";

    // ---- MDCD interpolation parameters (used by mdcd_linear / mdcd_hybrid) ----

    Real mdcd_diss    = 0.02;   ///< MDCD dissipation coefficient
    Real mdcd_disp    = 0.046;   ///< MDCD dispersion coefficient
    Real mdcd_sai_ref = 0.4;   ///< Discontinuity detector threshold for MDCD hybrid

    // ---- Riemann solver ----

    /// Riemann solver type: "roe", "rusanov", "hll", "hllc"
    std::string riemann_type = "roe";

    /// Entropy fix coefficient (Harten correction). 0 = no entropy fix.
    Real entropy_fix_eps = 0.1;

    // ---- Viscous flux ----

    /// Viscous computation mode:
    ///   "none"       — inviscid only, RHS_vis = 0
    ///   "constant"   — fixed viscosity: μ* = mu_const / Re
    ///   "sutherland" — Sutherland's law: μ* = T^(1.5)·(1+S*)/(T+S*) / Re
    std::string viscous_type = "constant";

    /// Constant viscosity multiplier (viscous_type="constant" only, typically 1.0).
    Real mu_const = 1.0;

    /// Sutherland constant [dimensional, K]. S = 110.4 K for air.
    /// Non-dimensionalized in finalize(): sutherland_S_nd = S / T_ref.
    Real sutherland_S      = 110.4;   ///< Dimensional Sutherland constant [K]
    Real sutherland_S_nd   = 0.0;     ///< Non-dimensional S* = S / T_ref (set by finalize())

    // =========================================================================
    // Initialization
    // =========================================================================

    std::string init_type = "uniform";  ///< Initial condition: "uniform" | "poiseuille" | "riemann_2d" | "isentropic_vortex"

    // ---- Poiseuille flow parameters ----
    Real poiseuille_umax  = 1.0;   ///< Theoretical max velocity (parabolic profile peak)
    Real poiseuille_y_min = -1.0;  ///< Lower wall y-coordinate
    Real poiseuille_y_max = 1.0;   ///< Upper wall y-coordinate

    // ---- 2D Riemann problem parameters ----
    Int  riemann_config  = 3;      ///< Riemann configuration (3 = classic 4-shock, 12 = alternate)
    Real riemann_x_split = 0.5;    ///< Quadrant split x-coordinate (domain [0,1])
    Real riemann_y_split = 0.5;    ///< Quadrant split y-coordinate (domain [0,1])

    // ---- Isentropic vortex parameters ----
    Real isentropic_vortex_strength = 5.0; ///< Vortex strength β (perturbation amplitude)
    Real isentropic_vortex_radius   = 0.05;///< Vortex core radius R
    Real isentropic_vortex_xc       = 0.25;///< Initial vortex center x-coordinate
    Real isentropic_vortex_yc       = 0.25;///< Initial vortex center y-coordinate
    Real isentropic_vortex_u_inf    = 1.0; ///< Mean flow x-velocity (overrides AoA)
    Real isentropic_vortex_v_inf    = 0.0; ///< Mean flow y-velocity (overrides beta)

    // =========================================================================
    // Body force (source term in momentum equations)
    // =========================================================================

    /// Body force type: "none" (no source), "constant" (uniform acceleration vector)
    std::string body_force_type = "none";

    Real body_force_x = 0.0;   ///< Non-dimensional body force x-component (acceleration)
    Real body_force_y = 0.0;   ///< Non-dimensional body force y-component (acceleration)
    Real body_force_z = 0.0;   ///< Non-dimensional body force z-component (acceleration)

    /// Include body-force work in the energy equation (ρ f·u).
    /// Set false for low-Mach / incompressible-like flows (e.g. Poiseuille)
    /// where the energy equation is decoupled from momentum and adiabatic walls
    /// would otherwise cause thermal runaway.
    bool body_force_energy = true;

    // =========================================================================
    // Output control
    // =========================================================================

    /// Grid metrics computation: "auto" (SCMM from grid) | "uniform" (analytical Cartesian)
    std::string grid_metrics = "auto";

    /// Solution output format: "tecplot" (implemented) | "cgns" (reserved) | "vtk" (reserved)
    std::string output_format = "tecplot";

    /// Output directory ("" = current working directory)
    std::string output_dir = "";

    /// Physical time interval for output (0 = disabled, use step-based only)
    Real output_time_interval = 0.0;

    /// History monitor: x-locations for u-velocity cross-section averaging.
    /// Comma-separated values parsed from config, e.g. "0.0, 3.14159"
    /// Default empty = no history monitoring.
    std::string history_x_locations = "";

    // =========================================================================
    // Boundary condition
    // =========================================================================

    std::string wall_type = "adiabatic_noslip"; ///< Wall BC: "adiabatic_noslip" | "isothermal_noslip" | "slip"

    /// Wall temperature for isothermal wall (non-dimensional, T* = T / T_ref).
    /// Default 1.0 = same as free-stream static temperature.
    Real wall_temperature = 1.0;

    // =========================================================================
    // Derived non-dimensional constants
    // =========================================================================

    /// Compute all derived quantities.
    /// Must be called after setting gamma, Mach, Re, L_ref, rho_ref, T_ref.
    /// Computes U_ref, p_ref, mu_ref, and the non-dimensional EOS factor.
    void finalize();

    /// Non-dimensional equation-of-state factor: p* = rho* * T* * eos_factor
    Real eos_factor() const;

    /// Free-stream velocity vector in Cartesian coordinates (non-dimensional).
    void free_stream_velocity(Real& u, Real& v, Real& w) const;
};

#include "config.hxx"
