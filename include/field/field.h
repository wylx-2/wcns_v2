#pragma once

#include "utils/types.h"
#include <map>
#include <string>

/// @file field.h
/// @brief Field class — manages all computed quantities on a structured block.
///
/// Stores primitive variables, conservative variables, inviscid and viscous
/// fluxes in each coordinate direction.  Also provides an extension mechanism
/// for additional variables (turbulence model, passive scalars, etc.).
///
/// All arrays are cell-centered with size (nci × ncj × nck), matching the
/// primary computational domain defined by Grid.

// ============================================================================
// Sub-structs: 5-component fluid state vectors
// ============================================================================

/// Primitive variables:  (rho, u, v, w, p)
struct PrimitiveVars {
    MultiArray3D<Real> rho;  ///< Density
    MultiArray3D<Real> u;    ///< x-velocity
    MultiArray3D<Real> v;    ///< y-velocity
    MultiArray3D<Real> w;    ///< z-velocity
    MultiArray3D<Real> p;    ///< Pressure

    void allocate(Int ni, Int nj, Int nk);
    void fill(Real val);
};

/// Conservative variables:  (rho, rho*u, rho*v, rho*w, rho*E)
struct ConservativeVars {
    MultiArray3D<Real> rho;  ///< Density (same as primitive)
    MultiArray3D<Real> rhou; ///< x-momentum
    MultiArray3D<Real> rhov; ///< y-momentum
    MultiArray3D<Real> rhow; ///< z-momentum
    MultiArray3D<Real> rhoE; ///< Total energy per unit volume

    void allocate(Int ni, Int nj, Int nk);
    void fill(Real val);
};

/// Inviscid or viscous flux vector in one coordinate direction.
/// Components:  (F_rho, F_rhou, F_rhov, F_rhow, F_rhoE)
struct FluxVars {
    MultiArray3D<Real> f1;  ///< Density flux  (continuity)
    MultiArray3D<Real> f2;  ///< x-momentum flux
    MultiArray3D<Real> f3;  ///< y-momentum flux
    MultiArray3D<Real> f4;  ///< z-momentum flux
    MultiArray3D<Real> f5;  ///< Energy flux

    void allocate(Int ni, Int nj, Int nk);
    void fill(Real val);
};

// ============================================================================
// Field class
// ============================================================================

/// Holds all flow-field quantities for one structured block.
///
/// Core arrays:
///   - prim  — primitive variables     (cell-centered)
///   - cons  — conservative variables  (cell-centered)
///   - inv_xi, inv_eta, inv_zeta  — inviscid fluxes at faces (3 directions)
///   - vis_xi, vis_eta, vis_zeta  — viscous fluxes  at faces (3 directions)
///
/// Extension fields (turbulence k/ω, passive scalars, etc.) can be added
/// dynamically via add_extension() / ext().
///
/// Usage:
///   Field f;
///   f.allocate(nci, ncj, nck);
///   f.add_extension("k_turb");
///   f.ext("k_turb").fill(0.0);
class Field {
public:
    // ---- Core state vectors ----
    PrimitiveVars    prim;       ///< Primitive variables at cell centers
    ConservativeVars cons;       ///< Conservative variables at cell centers

    // ---- Inviscid (Euler) fluxes at face half-nodes ----
    FluxVars inv_xi;             ///< Inviscid flux through ξ-face (i+1/2)
    FluxVars inv_eta;            ///< Inviscid flux through η-face (j+1/2)
    FluxVars inv_zeta;           ///< Inviscid flux through ζ-face (k+1/2)

    // ---- Interpolated left/right states at face half-nodes (for Riemann solver) ----
    ConservativeVars ql_xi;      ///< Left state at ξ-face (i+1/2)
    ConservativeVars qr_xi;      ///< Right state at ξ-face (i+1/2)
    ConservativeVars ql_eta;     ///< Left state at η-face (j+1/2)
    ConservativeVars qr_eta;     ///< Right state at η-face (j+1/2)
    ConservativeVars ql_zeta;    ///< Left state at ζ-face (k+1/2)
    ConservativeVars qr_zeta;    ///< Right state at ζ-face (k+1/2)

    // ---- Viscous fluxes at face half-nodes ----
    FluxVars vis_xi;             ///< Viscous flux through ξ-face (i+1/2)
    FluxVars vis_eta;            ///< Viscous flux through η-face (j+1/2)
    FluxVars vis_zeta;           ///< Viscous flux through ζ-face (k+1/2)

    // ---- Cell-center Cartesian viscous flux vectors (5c output) ----
    // F_vis (x-direction), G_vis (y-direction), H_vis (z-direction)
    // Each is a FluxVars (f1..f5) at cell centers (nci×ncj×nck).
    // f1 (density flux) ≡ 0 for viscous terms.
    // These are interpolated to faces in 5d to assemble vis_xi/eta/zeta.
    FluxVars vis_x;              ///< F_vis — Cartesian x-direction viscous flux
    FluxVars vis_y;              ///< G_vis — Cartesian y-direction viscous flux
    FluxVars vis_z;              ///< H_vis — Cartesian z-direction viscous flux

    // ---- Face-interpolated physical quantities (for viscous flux chain rule) ----
    // ξ-face (i+1/2): size (nci+1)×ncj×nck
    MultiArray3D<Real> u_face_xi, v_face_xi, w_face_xi, T_face_xi;
    // η-face (j+1/2): size nci×(ncj+1)×nck
    MultiArray3D<Real> u_face_eta, v_face_eta, w_face_eta, T_face_eta;
    // ζ-face (k+1/2): size nci×ncj×(nck+1)
    MultiArray3D<Real> u_face_zeta, v_face_zeta, w_face_zeta, T_face_zeta;

    // ---- Velocity and temperature gradients (for viscous flux, cell-center) ----
    MultiArray3D<Real> du_dx, du_dy, du_dz;
    MultiArray3D<Real> dv_dx, dv_dy, dv_dz;
    MultiArray3D<Real> dw_dx, dw_dy, dw_dz;
    MultiArray3D<Real> dT_dx, dT_dy, dT_dz;

    // ---- 5d temporary: face-interpolated Cartesian viscous flux components ----
    // Allocated per direction by interp_cart_flux_to_faces, freed after assembly.
    // Peak memory ~12 face arrays for the current direction.
    MultiArray3D<Real> vis_x_f2_face, vis_x_f3_face, vis_x_f4_face, vis_x_f5_face;
    MultiArray3D<Real> vis_y_f2_face, vis_y_f3_face, vis_y_f4_face, vis_y_f5_face;
    MultiArray3D<Real> vis_z_f2_face, vis_z_f3_face, vis_z_f4_face, vis_z_f5_face;

    // ---- Conservative variable right-hand side (for time stepping) ----
    ConservativeVars rhs;        ///< RHS ∂(Q/J)/∂t = -(∂F̂/∂ξ+∂Ĝ/∂η+∂Ĥ/∂ζ) in conservative form

    // ---- Initial conservative state snapshot (for RK multi-stage) ----
    ConservativeVars Q0;         ///< Q^(0) — conserved state at the start of each time step

    // ========================================================================
    // Allocation
    // ========================================================================

    /// Allocate all core arrays with given cell-center dimensions.
    /// Face-flux arrays are allocated at (nci+1, ncj, nck) etc.
    void allocate(Int nci, Int ncj, Int nck);

    // ========================================================================
    // Primitive ↔ Conservative conversion
    // ========================================================================

    /// Convert primitive variables to conservative variables.
    /// EOS: p = (gamma-1) * (rhoE - 0.5*rho*|u|^2)  — used to validate.
    /// The conversion is the forward direction:
    ///   cons.rho  = prim.rho
    ///   cons.rhou = prim.rho * prim.u
    ///   cons.rhov = prim.rho * prim.v
    ///   cons.rhow = prim.rho * prim.w
    ///   cons.rhoE = prim.p/(gamma-1) + 0.5*prim.rho*(u^2+v^2+w^2)
    void prim_to_cons(Real gamma);

    /// Convert conservative variables to primitive variables.
    ///   prim.rho = cons.rho
    ///   prim.u   = cons.rhou / cons.rho
    ///   prim.v   = cons.rhov / cons.rho
    ///   prim.w   = cons.rhow / cons.rho
    ///   prim.p   = (gamma-1)*(cons.rhoE - 0.5*(rhou^2+rhov^2+rhow^2)/cons.rho)
    void cons_to_prim(Real gamma);

    // ========================================================================
    // Extension fields (for future additions: turbulence, scalars, etc.)
    // ========================================================================

    /// Add a named extension field (allocated with same size as core arrays).
    /// Does nothing if the field already exists.
    void add_extension(const std::string& name);

    /// Check if an extension field exists.
    bool has_extension(const std::string& name) const;

    /// Access an extension field.  Throws std::runtime_error if not found.
    MultiArray3D<Real>&       ext(const std::string& name);
    const MultiArray3D<Real>& ext(const std::string& name) const;

    // ========================================================================
    // Query
    // ========================================================================

    Int ni() const;  ///< Number of cells in i-direction
    Int nj() const;  ///< Number of cells in j-direction
    Int nk() const;  ///< Number of cells in k-direction

private:
    Int nci_, ncj_, nck_;   // cell-center dimensions
    bool allocated_ = false;

    std::map<std::string, MultiArray3D<Real>> ext_fields_;
};

#include "field.hxx"
