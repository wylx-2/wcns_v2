#pragma once

#include "wcns_v2/utils/types.h"

// Forward declarations
class LocalBlock;
struct Config;
class FluxHaloExchange;
class HaloExchange;

/// @file viscid_rhs.h
/// @brief Viscous RHS computation — velocity/temperature gradients and viscous flux.
///
/// Subtask 5a: Interpolate primitive (u,v,w,T) from cell centers to face half-nodes.
/// Subtask 5b: Compute d(u,v,w,T)/d(x,y,z) at cell centers via chain rule:
///
///   dφ/dx = 1/J * [ d(φ·ξ̂_x)/dξ + d(φ·η̂_x)/dη + d(φ·ζ̂_x)/dζ ]
///
/// using SCMM face metric coefficients (ξ̂_x, ξ̂_y, ξ̂_z, …) and the 6th-order
/// centered differentiation from half-node to cell center.
///
/// Subtask 5c: Compute Cartesian viscous flux vectors (F_vis, G_vis, H_vis) at
///   cell centers from the gradients computed in 5b.  Stress tensor τ and heat
///   flux q are intermediate quantities — only the final flux vectors vis_x,
///   vis_y, vis_z (FluxVars, 5 components each) are stored in Field.
/// Subtask 5d-5e (face flux assembly, differentiation, RHS assembly) to follow.

class ViscidRHS {
public:
    // ========================================================================
    // 5a — Interpolate primitive variables to face half-nodes
    // ========================================================================

    /// Interpolate (u,v,w,T) from cell centers to faces in all 3 directions.
    ///
    /// Reads prim.u, prim.v, prim.w from lb.field, computes T from p and rho,
    /// and writes to u_face_xi, v_face_xi, ..., T_face_zeta in lb.field.
    ///
    /// T = gamma * Mach^2 * p / rho  (non-dimensional EOS)
    static void interp_to_faces(LocalBlock& lb, const Config& cfg);

    // ========================================================================
    // 5b — Compute velocity/temperature gradients at cell centers
    // ========================================================================

    /// Compute d(u,v,w,T)/d(x,y,z) at all cell centers (including ghost cells).
    ///
    /// Prerequisites:
    ///   1. interp_to_faces() must have been called (fills u_face_xi etc.)
    ///   2. Face-interpolated values at connectivity boundaries must be exchanged
    ///      (via FluxHaloExchange::exchange_face_arrays, handled internally).
    ///
    /// The 12 gradient arrays are written to lb.field:
    ///   du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz,
    ///   dw_dx, dw_dy, dw_dz, dT_dx, dT_dy, dT_dz
    ///
    /// After this call, gradients at connectivity-boundary ghost cells must be
    /// exchanged via HaloExchange::exchange_multi() before using them in 5c.
    static void compute_gradients(LocalBlock& lb,
                                   const std::vector<LocalBlock>& all_blocks,
                                   FluxHaloExchange& flux_ex,
                                   const Config& cfg);

    // ========================================================================
    // 5c — Compute cell-center Cartesian viscous flux vectors
    // ========================================================================

    /// Compute F_vis, G_vis, H_vis (Cartesian viscous flux vectors) at all
    /// cell centers (including ghost cells).
    ///
    /// Prerequisites:
    ///   1. compute_gradients() must have been called.
    ///   2. Gradient halos must have been exchanged (exchange_gradient_halos).
    ///
    /// Reads gradients (du_dx, du_dy, ..., dT_dz) and prim.{u,v,w,p,rho}
    /// from Field.  Computes μ via the configured viscous_type model, then
    /// stress tensor τ and heat flux q as intermediate quantities, and
    /// assembles the final flux vectors.
    ///
    /// Writes to lb.field:
    ///   vis_x (F_vis), vis_y (G_vis), vis_z (H_vis) — FluxVars at nci×ncj×nck
    ///
    /// If cfg.viscous_type == "none", sets all components to zero.
    ///
    /// After this call, vis_x/y/z ghost cells must be exchanged via
    /// exchange_viscous_flux_halos() before 5d interpolation.
    static void compute_cell_viscous_flux(LocalBlock& lb, const Config& cfg);

    // ========================================================================
    // 5d — Assemble viscous face fluxes from cell-center Cartesian flux vectors
    // ========================================================================

    /// 5d-step1 — Interpolate Cartesian viscous flux components to faces.
    ///
    /// Interpolates vis_x.{f2,f3,f4,f5}, vis_y.{f2,f3,f4,f5}, vis_z.{f2,f3,f4,f5}
    /// (12 scalar arrays) from cell centers to face half-nodes in direction @p dir.
    /// Results are stored in lb.field temporary face arrays (vis_x_f2_face, ...).
    ///
    /// Must be called on ALL blocks for the same @p dir before
    /// exchange_and_assemble_face_flux().
    ///
    /// @param dir  0=ξ, 1=η, 2=ζ
    static void interp_cart_flux_to_faces(LocalBlock& lb, int dir, const Config& cfg);

    /// 5d-step2+3 — Exchange face-interpolated Cartesian fluxes and assemble vis_xi/eta/zeta.
    ///
    /// Prerequisites: interp_cart_flux_to_faces() must have been called on ALL blocks
    /// for the same @p dir (the 12 temporary face arrays in Field must be populated).
    ///
    /// 1. MPI exchange via FluxHaloExchange::exchange_face_arrays (12 arrays)
    /// 2. Same-process local copy via copy_local_face_array (12 arrays)
    /// 3. Assemble vis_xi/eta/zeta: f1=0, f2..f5 = Σ vis_*_face_comp × face_metric_comp
    ///
    /// The 12 temporary face arrays are deallocated after assembly.
    ///
    /// @param dir  0=ξ, 1=η, 2=ζ
    static void exchange_and_assemble_face_flux(LocalBlock& lb, int dir,
                                                 const std::vector<LocalBlock>& all_blocks,
                                                 FluxHaloExchange& flux_ex,
                                                 const Config& cfg);

    // ========================================================================
    // 5e — Compute viscous RHS contribution from face fluxes
    // ========================================================================

    /// Compute viscous RHS contribution at interior cells.
    ///
    /// Prerequisites: exchange_and_assemble_face_flux() must have been called
    /// for all 3 directions (fills vis_xi, vis_eta, vis_zeta in lb.field).
    ///
    /// Differentiates vis_xi, vis_eta, vis_zeta at interior cells using 6th-order
    /// centered differences.  Subtracts each flux derivative from lb.field.rhs
    /// (accumulating with the inviscid RHS from InviscidRHS::compute).
    ///
    /// Only interior cells [3..nci-4]×[3..ncj-4]×[3..nck-4] are updated.
    static void compute_rhs(LocalBlock& lb);

private:
    /// Build the face_is_periodic mask from the block's neighbor info.
    static void build_periodic_mask(const LocalBlock& lb, bool fp[6]);

    /// Compute T = gamma * Mach^2 * p / rho at cell centers.
    /// Fills a temporary cell-center array; caller provides it pre-allocated.
    static void compute_temperature(const LocalBlock& lb, const Config& cfg,
                                     MultiArray3D<Real>& T_cell);

    /// Form face half-node products for one physical quantity in one direction.
    /// phi_face = face-interpolated value (u_face_xi etc.)
    /// Writes 3 product arrays: phi*face_?_x, phi*face_?_y, phi*face_?_z
    static void form_face_products(
            const MultiArray3D<Real>& phi_face,
            const MultiArray3D<Real>& face_m_x,
            const MultiArray3D<Real>& face_m_y,
            const MultiArray3D<Real>& face_m_z,
            MultiArray3D<Real>& prod_x,
            MultiArray3D<Real>& prod_y,
            MultiArray3D<Real>& prod_z);

    /// Copy face array slices between local blocks on the same process.
    /// Reads from neighbor's face array and writes to our ghost face region.
    /// Used for connectivity-boundary face exchange when neighbor is local.
    static void copy_local_face_array(
            MultiArray3D<Real>& my_face_arr,
            const MultiArray3D<Real>& nbr_face_arr,
            int dir,
            const LocalBlock& block, const LocalBlock& neighbor,
            int face);

    /// Compute gradient contribution from one direction.
    /// dφ/dx += J^-1 * d(φ·ξ̂_x)/dξ, etc.
    static void accumulate_gradient(
            const MultiArray3D<Real>& dxi,   ///< d(φ·ξ̂_*)/dξ
            const MultiArray3D<Real>& deta,  ///< d(φ·η̂_*)/dη
            const MultiArray3D<Real>& dzeta, ///< d(φ·ζ̂_*)/dζ
            const MultiArray3D<Real>& ja_inv,
            MultiArray3D<Real>& grad);       ///< dφ/d* (accumulated)
};

#include "viscid_rhs.hxx"
