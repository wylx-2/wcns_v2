#pragma once

#include "viscid_rhs.h"
#include "core/config.h"
#include "field/field.h"
#include "grid/grid.h"
#include "parallel/local_block.h"
#include "parallel/flux_halo_exchange.h"
#include "parallel/halo_exchange.h"
#include "scheme/interp_diff.h"
#include <cmath>
#include <vector>

// ============================================================================
// Helper: build periodic mask
// ============================================================================

inline void ViscidRHS::build_periodic_mask(const LocalBlock& lb, bool fp[6]) {
    for (int f = 0; f < 6; ++f) {
        fp[f] = lb.neighbors[f].is_periodic;
    }
}

// ============================================================================
// Helper: compute T = gamma * Mach^2 * p / rho
// ============================================================================

inline void ViscidRHS::compute_temperature(const LocalBlock& lb, const Config& cfg,
                                            MultiArray3D<Real>& T_cell) {
    Int nci = lb.field.ni();
    Int ncj = lb.field.nj();
    Int nck = lb.field.nk();
    Real factor = cfg.gamma * cfg.Mach * cfg.Mach;

    for (Int k = 0; k < nck; ++k) {
    for (Int j = 0; j < ncj; ++j) {
    for (Int i = 0; i < nci; ++i) {
        T_cell(i,j,k) = factor * lb.field.prim.p(i,j,k) / lb.field.prim.rho(i,j,k);
    }}}
}

// ============================================================================
// Helper: form face half-node products
// ============================================================================

inline void ViscidRHS::form_face_products(
        const MultiArray3D<Real>& phi_face,
        const MultiArray3D<Real>& face_m_x,
        const MultiArray3D<Real>& face_m_y,
        const MultiArray3D<Real>& face_m_z,
        MultiArray3D<Real>& prod_x,
        MultiArray3D<Real>& prod_y,
        MultiArray3D<Real>& prod_z) {

    Int ni = phi_face.ni();
    Int nj = phi_face.nj();
    Int nk = phi_face.nk();

    for (Int k = 0; k < nk; ++k) {
    for (Int j = 0; j < nj; ++j) {
    for (Int i = 0; i < ni; ++i) {
        Real phi = phi_face(i,j,k);
        prod_x(i,j,k) = phi * face_m_x(i,j,k);
        prod_y(i,j,k) = phi * face_m_y(i,j,k);
        prod_z(i,j,k) = phi * face_m_z(i,j,k);
    }}}
}

// ============================================================================
// Helper: local face array copy between same-process blocks
// ============================================================================

inline void ViscidRHS::copy_local_face_array(
        MultiArray3D<Real>& my_face_arr,
        const MultiArray3D<Real>& nbr_face_arr,
        int dir,
        const LocalBlock& block, const LocalBlock& neighbor,
        int face) {

    // This is a minimal implementation that copies face-interpolated values
    // directly between same-process blocks.  The FluxHaloExchange info
    // (send/recv ranges) is private, so we recompute the ranges here.
    // For a production code we would make FluxFaceInfo accessible instead.

    Int ng = block.grid.ng;
    Int ndim, dim1, dim2;

    switch (dir) {
    case 0: ndim = block.field.ni();     dim1 = block.field.nj(); dim2 = block.field.nk(); break;
    case 1: ndim = block.field.nj();     dim1 = block.field.ni(); dim2 = block.field.nk(); break;
    case 2: ndim = block.field.nk();     dim1 = block.field.ni(); dim2 = block.field.nj(); break;
    default: return;
    }

    Int n_faces = ng + 1;
    Int nbr_ndim;
    switch (dir) {
    case 0: nbr_ndim = neighbor.field.ni(); break;
    case 1: nbr_ndim = neighbor.field.nj(); break;
    case 2: nbr_ndim = neighbor.field.nk(); break;
    default: nbr_ndim = 0; break;
    }

    Int r0, nbr_s0;
    if (face == 1 || face == 3 || face == 5) {
        // MAX face
        r0     = ndim - ng;
        nbr_s0 = ng;
    } else {
        // MIN face
        r0     = 0;
        nbr_s0 = nbr_ndim - 2 * ng;
    }

    for (Int d = 0; d < n_faces; ++d) {
        Int my_idx  = r0 + d;
        Int nbr_idx = nbr_s0 + d;
        for (Int b = 0; b < dim2; ++b) {
        for (Int a = 0; a < dim1; ++a) {
            Real val = 0.0;
            switch (dir) {
            case 0: val = nbr_face_arr(nbr_idx, a, b); break;
            case 1: val = nbr_face_arr(a, nbr_idx, b); break;
            case 2: val = nbr_face_arr(a, b, nbr_idx); break;
            }
            switch (dir) {
            case 0: my_face_arr(my_idx, a, b) = val; break;
            case 1: my_face_arr(a, my_idx, b) = val; break;
            case 2: my_face_arr(a, b, my_idx) = val; break;
            }
        }}
    }
}

// ============================================================================
// Helper: accumulate gradient contribution from one direction
// ============================================================================

inline void ViscidRHS::accumulate_gradient(
        const MultiArray3D<Real>& dxi,
        const MultiArray3D<Real>& deta,
        const MultiArray3D<Real>& dzeta,
        const MultiArray3D<Real>& ja_inv,
        MultiArray3D<Real>& grad) {

    Int ni = grad.ni();
    Int nj = grad.nj();
    Int nk = grad.nk();

    for (Int k = 0; k < nk; ++k) {
    for (Int j = 0; j < nj; ++j) {
    for (Int i = 0; i < ni; ++i) {
        Real sum = dxi(i,j,k) + deta(i,j,k) + dzeta(i,j,k);
        grad(i,j,k) += sum * ja_inv(i,j,k);
    }}}
}

// ============================================================================
// 5a — Interpolate primitive variables to face half-nodes
// ============================================================================

inline void ViscidRHS::interp_to_faces(LocalBlock& lb, const Config& cfg) {
    Int nci = lb.field.ni();
    Int ncj = lb.field.nj();
    Int nck = lb.field.nk();
    Int ng  = lb.grid.ng;

    bool fp[6];
    build_periodic_mask(lb, fp);

    // Compute T at cell centers
    MultiArray3D<Real> T_cell;
    T_cell.allocate(nci, ncj, nck);
    compute_temperature(lb, cfg, T_cell);

    // Interpolate u, v, w, T to faces in all 3 directions

    // ---- ξ direction ----
    InterpDiff::interp_to_faces(lb.field.prim.u, lb.field.u_face_xi, 0, ng, fp);
    InterpDiff::interp_to_faces(lb.field.prim.v, lb.field.v_face_xi, 0, ng, fp);
    InterpDiff::interp_to_faces(lb.field.prim.w, lb.field.w_face_xi, 0, ng, fp);
    InterpDiff::interp_to_faces(T_cell,           lb.field.T_face_xi, 0, ng, fp);

    // ---- η direction ----
    InterpDiff::interp_to_faces(lb.field.prim.u, lb.field.u_face_eta, 1, ng, fp);
    InterpDiff::interp_to_faces(lb.field.prim.v, lb.field.v_face_eta, 1, ng, fp);
    InterpDiff::interp_to_faces(lb.field.prim.w, lb.field.w_face_eta, 1, ng, fp);
    InterpDiff::interp_to_faces(T_cell,           lb.field.T_face_eta, 1, ng, fp);

    // ---- ζ direction ----
    InterpDiff::interp_to_faces(lb.field.prim.u, lb.field.u_face_zeta, 2, ng, fp);
    InterpDiff::interp_to_faces(lb.field.prim.v, lb.field.v_face_zeta, 2, ng, fp);
    InterpDiff::interp_to_faces(lb.field.prim.w, lb.field.w_face_zeta, 2, ng, fp);
    InterpDiff::interp_to_faces(T_cell,           lb.field.T_face_zeta, 2, ng, fp);
}

// ============================================================================
// 5b — Compute velocity/temperature gradients at cell centers
// ============================================================================

inline void ViscidRHS::compute_gradients(LocalBlock& lb,
                                          const std::vector<LocalBlock>& all_blocks,
                                          FluxHaloExchange& flux_ex,
                                          const Config& /*cfg*/) {

    Int nci = lb.field.ni();
    Int ncj = lb.field.nj();
    Int nck = lb.field.nk();
    Int ng  = lb.grid.ng;

    bool fp[6];
    build_periodic_mask(lb, fp);
    Real dh = 1.0;

    // Compute ja_inv = 1 / jacobian
    MultiArray3D<Real> ja_inv;
    ja_inv.allocate(nci, ncj, nck);
    for (Int k = 0; k < nck; ++k) {
    for (Int j = 0; j < ncj; ++j) {
    for (Int i = 0; i < nci; ++i) {
        ja_inv(i,j,k) = 1.0 / lb.grid.jacobian(i,j,k);
    }}}

    // Zero out gradient arrays
    lb.field.du_dx.fill(0); lb.field.du_dy.fill(0); lb.field.du_dz.fill(0);
    lb.field.dv_dx.fill(0); lb.field.dv_dy.fill(0); lb.field.dv_dz.fill(0);
    lb.field.dw_dx.fill(0); lb.field.dw_dy.fill(0); lb.field.dw_dz.fill(0);
    lb.field.dT_dx.fill(0); lb.field.dT_dy.fill(0); lb.field.dT_dz.fill(0);

    // Temporary arrays for face products and their derivatives.
    // Reused across directions (reallocated when direction changes size).

    // ---- Process each direction ----
    // The product arrays are formed, exchanged, differentiated, and accumulated
    // one direction at a time to minimize memory footprint.

    // ξ direction ==============================================================
    {
        Int f_ni = nci + 1, f_nj = ncj, f_nk = nck;

        // Temporary product arrays for ξ-direction (reused per quantity)
        MultiArray3D<Real> px, py, pz, dpx, dpy, dpz;
        px.allocate(f_ni, f_nj, f_nk);
        py.allocate(f_ni, f_nj, f_nk);
        pz.allocate(f_ni, f_nj, f_nk);
        dpx.allocate(nci, ncj, nck);
        dpy.allocate(nci, ncj, nck);
        dpz.allocate(nci, ncj, nck);

        // ---- u ----
        form_face_products(lb.field.u_face_xi,
                           lb.grid.face_xi_x, lb.grid.face_xi_y, lb.grid.face_xi_z,
                           px, py, pz);
        // Exchange products
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 0, lb);
            // Local copy for same-process neighbors
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 0) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                // Find the neighbor block
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        // Compute neighbor's products on-the-fly and copy
                        // Form products from neighbor's data for the send region,
                        // then copy to our recv region.
                        // For simplicity, we copy the face-interpolated values
                        // from the neighbor, which makes our product correct.
                        copy_local_face_array(lb.field.u_face_xi, nbr_block.field.u_face_xi, 0, lb, nbr_block, face);
                        // Now re-form products (our ghost face values are updated)
                        form_face_products(lb.field.u_face_xi,
                                           lb.grid.face_xi_x, lb.grid.face_xi_y, lb.grid.face_xi_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        // Differentiate
        InterpDiff::derivative_from_faces(px, dpx, 0, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 0, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 0, dh, ng, fp);
        // Accumulate (init with dxi part, deta/dzeta added in later directions)
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.du_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.du_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.du_dz(i,j,k) += dpz(i,j,k) * invJ;
        }

        // ---- v ----
        form_face_products(lb.field.v_face_xi,
                           lb.grid.face_xi_x, lb.grid.face_xi_y, lb.grid.face_xi_z,
                           px, py, pz);
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 0, lb);
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 0) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        copy_local_face_array(lb.field.v_face_xi, nbr_block.field.v_face_xi, 0, lb, nbr_block, face);
                        form_face_products(lb.field.v_face_xi,
                                           lb.grid.face_xi_x, lb.grid.face_xi_y, lb.grid.face_xi_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        InterpDiff::derivative_from_faces(px, dpx, 0, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 0, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 0, dh, ng, fp);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.dv_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.dv_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.dv_dz(i,j,k) += dpz(i,j,k) * invJ;
        }

        // ---- w ----
        form_face_products(lb.field.w_face_xi,
                           lb.grid.face_xi_x, lb.grid.face_xi_y, lb.grid.face_xi_z,
                           px, py, pz);
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 0, lb);
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 0) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        copy_local_face_array(lb.field.w_face_xi, nbr_block.field.w_face_xi, 0, lb, nbr_block, face);
                        form_face_products(lb.field.w_face_xi,
                                           lb.grid.face_xi_x, lb.grid.face_xi_y, lb.grid.face_xi_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        InterpDiff::derivative_from_faces(px, dpx, 0, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 0, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 0, dh, ng, fp);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.dw_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.dw_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.dw_dz(i,j,k) += dpz(i,j,k) * invJ;
        }

        // ---- T ----
        form_face_products(lb.field.T_face_xi,
                           lb.grid.face_xi_x, lb.grid.face_xi_y, lb.grid.face_xi_z,
                           px, py, pz);
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 0, lb);
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 0) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        copy_local_face_array(lb.field.T_face_xi, nbr_block.field.T_face_xi, 0, lb, nbr_block, face);
                        form_face_products(lb.field.T_face_xi,
                                           lb.grid.face_xi_x, lb.grid.face_xi_y, lb.grid.face_xi_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        InterpDiff::derivative_from_faces(px, dpx, 0, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 0, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 0, dh, ng, fp);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.dT_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.dT_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.dT_dz(i,j,k) += dpz(i,j,k) * invJ;
        }
    } // ξ direction

    // η direction ==============================================================
    {
        Int f_ni = nci, f_nj = ncj + 1, f_nk = nck;

        MultiArray3D<Real> px, py, pz, dpx, dpy, dpz;
        px.allocate(f_ni, f_nj, f_nk);
        py.allocate(f_ni, f_nj, f_nk);
        pz.allocate(f_ni, f_nj, f_nk);
        dpx.allocate(nci, ncj, nck);
        dpy.allocate(nci, ncj, nck);
        dpz.allocate(nci, ncj, nck);

        // ---- u ----
        form_face_products(lb.field.u_face_eta,
                           lb.grid.face_eta_x, lb.grid.face_eta_y, lb.grid.face_eta_z,
                           px, py, pz);
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 1, lb);
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 1) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        copy_local_face_array(lb.field.u_face_eta, nbr_block.field.u_face_eta, 1, lb, nbr_block, face);
                        form_face_products(lb.field.u_face_eta,
                                           lb.grid.face_eta_x, lb.grid.face_eta_y, lb.grid.face_eta_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        InterpDiff::derivative_from_faces(px, dpx, 1, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 1, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 1, dh, ng, fp);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.du_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.du_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.du_dz(i,j,k) += dpz(i,j,k) * invJ;
        }

        // ---- v ----
        form_face_products(lb.field.v_face_eta,
                           lb.grid.face_eta_x, lb.grid.face_eta_y, lb.grid.face_eta_z,
                           px, py, pz);
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 1, lb);
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 1) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        copy_local_face_array(lb.field.v_face_eta, nbr_block.field.v_face_eta, 1, lb, nbr_block, face);
                        form_face_products(lb.field.v_face_eta,
                                           lb.grid.face_eta_x, lb.grid.face_eta_y, lb.grid.face_eta_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        InterpDiff::derivative_from_faces(px, dpx, 1, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 1, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 1, dh, ng, fp);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.dv_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.dv_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.dv_dz(i,j,k) += dpz(i,j,k) * invJ;
        }

        // ---- w ----
        form_face_products(lb.field.w_face_eta,
                           lb.grid.face_eta_x, lb.grid.face_eta_y, lb.grid.face_eta_z,
                           px, py, pz);
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 1, lb);
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 1) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        copy_local_face_array(lb.field.w_face_eta, nbr_block.field.w_face_eta, 1, lb, nbr_block, face);
                        form_face_products(lb.field.w_face_eta,
                                           lb.grid.face_eta_x, lb.grid.face_eta_y, lb.grid.face_eta_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        InterpDiff::derivative_from_faces(px, dpx, 1, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 1, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 1, dh, ng, fp);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.dw_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.dw_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.dw_dz(i,j,k) += dpz(i,j,k) * invJ;
        }

        // ---- T ----
        form_face_products(lb.field.T_face_eta,
                           lb.grid.face_eta_x, lb.grid.face_eta_y, lb.grid.face_eta_z,
                           px, py, pz);
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 1, lb);
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 1) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        copy_local_face_array(lb.field.T_face_eta, nbr_block.field.T_face_eta, 1, lb, nbr_block, face);
                        form_face_products(lb.field.T_face_eta,
                                           lb.grid.face_eta_x, lb.grid.face_eta_y, lb.grid.face_eta_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        InterpDiff::derivative_from_faces(px, dpx, 1, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 1, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 1, dh, ng, fp);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.dT_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.dT_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.dT_dz(i,j,k) += dpz(i,j,k) * invJ;
        }
    } // η direction

    // ζ direction ==============================================================
    {
        Int f_ni = nci, f_nj = ncj, f_nk = nck + 1;

        MultiArray3D<Real> px, py, pz, dpx, dpy, dpz;
        px.allocate(f_ni, f_nj, f_nk);
        py.allocate(f_ni, f_nj, f_nk);
        pz.allocate(f_ni, f_nj, f_nk);
        dpx.allocate(nci, ncj, nck);
        dpy.allocate(nci, ncj, nck);
        dpz.allocate(nci, ncj, nck);

        // ---- u ----
        form_face_products(lb.field.u_face_zeta,
                           lb.grid.face_zeta_x, lb.grid.face_zeta_y, lb.grid.face_zeta_z,
                           px, py, pz);
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 2, lb);
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 2) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        copy_local_face_array(lb.field.u_face_zeta, nbr_block.field.u_face_zeta, 2, lb, nbr_block, face);
                        form_face_products(lb.field.u_face_zeta,
                                           lb.grid.face_zeta_x, lb.grid.face_zeta_y, lb.grid.face_zeta_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        InterpDiff::derivative_from_faces(px, dpx, 2, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 2, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 2, dh, ng, fp);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.du_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.du_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.du_dz(i,j,k) += dpz(i,j,k) * invJ;
        }

        // ---- v ----
        form_face_products(lb.field.v_face_zeta,
                           lb.grid.face_zeta_x, lb.grid.face_zeta_y, lb.grid.face_zeta_z,
                           px, py, pz);
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 2, lb);
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 2) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        copy_local_face_array(lb.field.v_face_zeta, nbr_block.field.v_face_zeta, 2, lb, nbr_block, face);
                        form_face_products(lb.field.v_face_zeta,
                                           lb.grid.face_zeta_x, lb.grid.face_zeta_y, lb.grid.face_zeta_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        InterpDiff::derivative_from_faces(px, dpx, 2, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 2, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 2, dh, ng, fp);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.dv_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.dv_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.dv_dz(i,j,k) += dpz(i,j,k) * invJ;
        }

        // ---- w ----
        form_face_products(lb.field.w_face_zeta,
                           lb.grid.face_zeta_x, lb.grid.face_zeta_y, lb.grid.face_zeta_z,
                           px, py, pz);
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 2, lb);
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 2) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        copy_local_face_array(lb.field.w_face_zeta, nbr_block.field.w_face_zeta, 2, lb, nbr_block, face);
                        form_face_products(lb.field.w_face_zeta,
                                           lb.grid.face_zeta_x, lb.grid.face_zeta_y, lb.grid.face_zeta_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        InterpDiff::derivative_from_faces(px, dpx, 2, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 2, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 2, dh, ng, fp);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.dw_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.dw_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.dw_dz(i,j,k) += dpz(i,j,k) * invJ;
        }

        // ---- T ----
        form_face_products(lb.field.T_face_zeta,
                           lb.grid.face_zeta_x, lb.grid.face_zeta_y, lb.grid.face_zeta_z,
                           px, py, pz);
        {
            std::vector<MultiArray3D<Real>*> arrs = {&px, &py, &pz};
            flux_ex.exchange_face_arrays(arrs, 2, lb);
            for (int face = 0; face < 6; ++face) {
                if (face / 2 != 2) continue;
                const auto& nbr_info = lb.neighbors[face];
                if (!nbr_info.active || nbr_info.is_periodic) continue;
                for (const auto& nbr_block : all_blocks) {
                    if (nbr_block.block_id == nbr_info.target_block) {
                        copy_local_face_array(lb.field.T_face_zeta, nbr_block.field.T_face_zeta, 2, lb, nbr_block, face);
                        form_face_products(lb.field.T_face_zeta,
                                           lb.grid.face_zeta_x, lb.grid.face_zeta_y, lb.grid.face_zeta_z,
                                           px, py, pz);
                        break;
                    }
                }
            }
        }
        InterpDiff::derivative_from_faces(px, dpx, 2, dh, ng, fp);
        InterpDiff::derivative_from_faces(py, dpy, 2, dh, ng, fp);
        InterpDiff::derivative_from_faces(pz, dpz, 2, dh, ng, fp);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            Real invJ = ja_inv(i,j,k);
            lb.field.dT_dx(i,j,k) += dpx(i,j,k) * invJ;
            lb.field.dT_dy(i,j,k) += dpy(i,j,k) * invJ;
            lb.field.dT_dz(i,j,k) += dpz(i,j,k) * invJ;
        }
    } // ζ direction
}

// ============================================================================
// 5c — Compute cell-center Cartesian viscous flux vectors
// ============================================================================

inline void ViscidRHS::compute_cell_viscous_flux(LocalBlock& lb, const Config& cfg) {

    Int nci = lb.field.ni();
    Int ncj = lb.field.nj();
    Int nck = lb.field.nk();

    // If inviscid mode: set all components to zero and return
    if (cfg.viscous_type == "none") {
        lb.field.vis_x.fill(0);
        lb.field.vis_y.fill(0);
        lb.field.vis_z.fill(0);
        return;
    }

    // ---- Pre-compute temperature at cell centers ----
    MultiArray3D<Real> T_cell;
    T_cell.allocate(nci, ncj, nck);
    compute_temperature(lb, cfg, T_cell);

    // Pre-compute thermal conductivity factor
    // k_nd = mu_nd / (Pr * (gamma-1) * Mach^2)
    Real k_factor = 1.0 / (cfg.Prandtl * (cfg.gamma - 1.0) * cfg.Mach * cfg.Mach);
    Real inv_Re  = 1.0 / cfg.Re;

    for (Int k = 0; k < nck; ++k) {
    for (Int j = 0; j < ncj; ++j) {
    for (Int i = 0; i < nci; ++i) {

        // ---- Viscosity coefficient ----
        Real mu_star;  // dimensionless mu (before Re scaling)
        if (cfg.viscous_type == "constant") {
            mu_star = cfg.mu_const;
        } else {
            // Sutherland
            Real T   = T_cell(i,j,k);
            Real Snd = cfg.sutherland_S_nd;
            mu_star = std::pow(T, 1.5) * (1.0 + Snd) / (T + Snd);
        }
        Real mu_nd = mu_star * inv_Re;

        // ---- Thermal conductivity ----
        Real k_nd = mu_nd * k_factor;

        // ---- Velocity gradients ----
        Real dudx = lb.field.du_dx(i,j,k);
        Real dudy = lb.field.du_dy(i,j,k);
        Real dudz = lb.field.du_dz(i,j,k);
        Real dvdx = lb.field.dv_dx(i,j,k);
        Real dvdy = lb.field.dv_dy(i,j,k);
        Real dvdz = lb.field.dv_dz(i,j,k);
        Real dwdx = lb.field.dw_dx(i,j,k);
        Real dwdy = lb.field.dw_dy(i,j,k);
        Real dwdz = lb.field.dw_dz(i,j,k);

        Real dTdx = lb.field.dT_dx(i,j,k);
        Real dTdy = lb.field.dT_dy(i,j,k);
        Real dTdz = lb.field.dT_dz(i,j,k);

        // ---- Divergence and stress tensor ----
        Real divV = dudx + dvdy + dwdz;
        Real lam_div = (-2.0 / 3.0) * mu_nd * divV;

        Real tau_xx = mu_nd * 2.0 * dudx + lam_div;
        Real tau_yy = mu_nd * 2.0 * dvdy + lam_div;
        Real tau_zz = mu_nd * 2.0 * dwdz + lam_div;
        Real tau_xy = mu_nd * (dudy + dvdx);
        Real tau_xz = mu_nd * (dudz + dwdx);
        Real tau_yz = mu_nd * (dvdz + dwdy);

        // ---- Heat flux ----
        Real q_x = -k_nd * dTdx;
        Real q_y = -k_nd * dTdy;
        Real q_z = -k_nd * dTdz;

        // ---- Primitive velocities ----
        Real u = lb.field.prim.u(i,j,k);
        Real v = lb.field.prim.v(i,j,k);
        Real w = lb.field.prim.w(i,j,k);

        // ---- F_vis (x-direction Cartesian viscous flux) ----
        lb.field.vis_x.f1(i,j,k) = 0.0;
        lb.field.vis_x.f2(i,j,k) = tau_xx;
        lb.field.vis_x.f3(i,j,k) = tau_xy;
        lb.field.vis_x.f4(i,j,k) = tau_xz;
        lb.field.vis_x.f5(i,j,k) = u * tau_xx + v * tau_xy + w * tau_xz - q_x;

        // ---- G_vis (y-direction Cartesian viscous flux) ----
        lb.field.vis_y.f1(i,j,k) = 0.0;
        lb.field.vis_y.f2(i,j,k) = tau_xy;  // τ_yx = τ_xy (symmetric)
        lb.field.vis_y.f3(i,j,k) = tau_yy;
        lb.field.vis_y.f4(i,j,k) = tau_yz;
        lb.field.vis_y.f5(i,j,k) = u * tau_xy + v * tau_yy + w * tau_yz - q_y;

        // ---- H_vis (z-direction Cartesian viscous flux) ----
        lb.field.vis_z.f1(i,j,k) = 0.0;
        lb.field.vis_z.f2(i,j,k) = tau_xz;  // τ_zx = τ_xz (symmetric)
        lb.field.vis_z.f3(i,j,k) = tau_yz;  // τ_zy = τ_yz (symmetric)
        lb.field.vis_z.f4(i,j,k) = tau_zz;
        lb.field.vis_z.f5(i,j,k) = u * tau_xz + v * tau_yz + w * tau_zz - q_z;

    }}}
}

// ============================================================================
// 5d-step1 — Interpolate Cartesian viscous flux components to faces
// ============================================================================

inline void ViscidRHS::interp_cart_flux_to_faces(LocalBlock& lb, int dir,
                                                  const Config& /*cfg*/) {

    Int nci = lb.field.ni();
    Int ncj = lb.field.nj();
    Int nck = lb.field.nk();
    Int ng  = lb.grid.ng;

    bool fp[6];
    build_periodic_mask(lb, fp);

    // Determine face array dimensions for this direction
    Int f_ni, f_nj, f_nk;
    switch (dir) {
    case 0: f_ni = nci + 1; f_nj = ncj;     f_nk = nck;     break;  // ξ-face
    case 1: f_ni = nci;     f_nj = ncj + 1; f_nk = nck;     break;  // η-face
    case 2: f_ni = nci;     f_nj = ncj;     f_nk = nck + 1; break;  // ζ-face
    default: return;
    }

    // Allocate 12 temp face arrays
    lb.field.vis_x_f2_face.allocate(f_ni, f_nj, f_nk);
    lb.field.vis_x_f3_face.allocate(f_ni, f_nj, f_nk);
    lb.field.vis_x_f4_face.allocate(f_ni, f_nj, f_nk);
    lb.field.vis_x_f5_face.allocate(f_ni, f_nj, f_nk);
    lb.field.vis_y_f2_face.allocate(f_ni, f_nj, f_nk);
    lb.field.vis_y_f3_face.allocate(f_ni, f_nj, f_nk);
    lb.field.vis_y_f4_face.allocate(f_ni, f_nj, f_nk);
    lb.field.vis_y_f5_face.allocate(f_ni, f_nj, f_nk);
    lb.field.vis_z_f2_face.allocate(f_ni, f_nj, f_nk);
    lb.field.vis_z_f3_face.allocate(f_ni, f_nj, f_nk);
    lb.field.vis_z_f4_face.allocate(f_ni, f_nj, f_nk);
    lb.field.vis_z_f5_face.allocate(f_ni, f_nj, f_nk);

    // Interpolate vis_x.{f2,f3,f4,f5} to faces
    InterpDiff::interp_to_faces(lb.field.vis_x.f2, lb.field.vis_x_f2_face, dir, ng, fp);
    InterpDiff::interp_to_faces(lb.field.vis_x.f3, lb.field.vis_x_f3_face, dir, ng, fp);
    InterpDiff::interp_to_faces(lb.field.vis_x.f4, lb.field.vis_x_f4_face, dir, ng, fp);
    InterpDiff::interp_to_faces(lb.field.vis_x.f5, lb.field.vis_x_f5_face, dir, ng, fp);

    // Interpolate vis_y.{f2,f3,f4,f5} to faces
    InterpDiff::interp_to_faces(lb.field.vis_y.f2, lb.field.vis_y_f2_face, dir, ng, fp);
    InterpDiff::interp_to_faces(lb.field.vis_y.f3, lb.field.vis_y_f3_face, dir, ng, fp);
    InterpDiff::interp_to_faces(lb.field.vis_y.f4, lb.field.vis_y_f4_face, dir, ng, fp);
    InterpDiff::interp_to_faces(lb.field.vis_y.f5, lb.field.vis_y_f5_face, dir, ng, fp);

    // Interpolate vis_z.{f2,f3,f4,f5} to faces
    InterpDiff::interp_to_faces(lb.field.vis_z.f2, lb.field.vis_z_f2_face, dir, ng, fp);
    InterpDiff::interp_to_faces(lb.field.vis_z.f3, lb.field.vis_z_f3_face, dir, ng, fp);
    InterpDiff::interp_to_faces(lb.field.vis_z.f4, lb.field.vis_z_f4_face, dir, ng, fp);
    InterpDiff::interp_to_faces(lb.field.vis_z.f5, lb.field.vis_z_f5_face, dir, ng, fp);
}

// ============================================================================
// 5d-step2+3 — Exchange face-interpolated Cartesian fluxes and assemble
// ============================================================================

inline void ViscidRHS::exchange_and_assemble_face_flux(
        LocalBlock& lb, int dir,
        const std::vector<LocalBlock>& all_blocks,
        FluxHaloExchange& flux_ex,
        const Config& /*cfg*/) {

    // Gather 12 temp face arrays for exchange
    std::vector<MultiArray3D<Real>*> cart_face = {
        &lb.field.vis_x_f2_face, &lb.field.vis_x_f3_face,
        &lb.field.vis_x_f4_face, &lb.field.vis_x_f5_face,
        &lb.field.vis_y_f2_face, &lb.field.vis_y_f3_face,
        &lb.field.vis_y_f4_face, &lb.field.vis_y_f5_face,
        &lb.field.vis_z_f2_face, &lb.field.vis_z_f3_face,
        &lb.field.vis_z_f4_face, &lb.field.vis_z_f5_face
    };

    // MPI exchange
    flux_ex.exchange_face_arrays(cart_face, dir, lb);

    // Same-process local copy
    for (int face = 0; face < 6; ++face) {
        if (face / 2 != dir) continue;  // only faces in the exchange direction
        const auto& nbr_info = lb.neighbors[face];
        if (!nbr_info.active || nbr_info.is_periodic) continue;

        for (const auto& nbr_block : all_blocks) {
            if (nbr_block.block_id == nbr_info.target_block) {
                copy_local_face_array(lb.field.vis_x_f2_face, nbr_block.field.vis_x_f2_face, dir, lb, nbr_block, face);
                copy_local_face_array(lb.field.vis_x_f3_face, nbr_block.field.vis_x_f3_face, dir, lb, nbr_block, face);
                copy_local_face_array(lb.field.vis_x_f4_face, nbr_block.field.vis_x_f4_face, dir, lb, nbr_block, face);
                copy_local_face_array(lb.field.vis_x_f5_face, nbr_block.field.vis_x_f5_face, dir, lb, nbr_block, face);
                copy_local_face_array(lb.field.vis_y_f2_face, nbr_block.field.vis_y_f2_face, dir, lb, nbr_block, face);
                copy_local_face_array(lb.field.vis_y_f3_face, nbr_block.field.vis_y_f3_face, dir, lb, nbr_block, face);
                copy_local_face_array(lb.field.vis_y_f4_face, nbr_block.field.vis_y_f4_face, dir, lb, nbr_block, face);
                copy_local_face_array(lb.field.vis_y_f5_face, nbr_block.field.vis_y_f5_face, dir, lb, nbr_block, face);
                copy_local_face_array(lb.field.vis_z_f2_face, nbr_block.field.vis_z_f2_face, dir, lb, nbr_block, face);
                copy_local_face_array(lb.field.vis_z_f3_face, nbr_block.field.vis_z_f3_face, dir, lb, nbr_block, face);
                copy_local_face_array(lb.field.vis_z_f4_face, nbr_block.field.vis_z_f4_face, dir, lb, nbr_block, face);
                copy_local_face_array(lb.field.vis_z_f5_face, nbr_block.field.vis_z_f5_face, dir, lb, nbr_block, face);
                break;
            }
        }
    }

    // Assemble vis_xi/eta/zeta using face metrics
    FluxVars* out = nullptr;
    const MultiArray3D<Real>* fm_x = nullptr;
    const MultiArray3D<Real>* fm_y = nullptr;
    const MultiArray3D<Real>* fm_z = nullptr;

    switch (dir) {
    case 0: out = &lb.field.vis_xi;   fm_x = &lb.grid.face_xi_x;   fm_y = &lb.grid.face_xi_y;   fm_z = &lb.grid.face_xi_z;   break;
    case 1: out = &lb.field.vis_eta;  fm_x = &lb.grid.face_eta_x;  fm_y = &lb.grid.face_eta_y;  fm_z = &lb.grid.face_eta_z;  break;
    case 2: out = &lb.field.vis_zeta; fm_x = &lb.grid.face_zeta_x; fm_y = &lb.grid.face_zeta_y; fm_z = &lb.grid.face_zeta_z; break;
    default: return;
    }

    Int ni = out->f1.ni();
    Int nj = out->f1.nj();
    Int nk = out->f1.nk();

    for (Int k = 0; k < nk; ++k) {
    for (Int j = 0; j < nj; ++j) {
    for (Int i = 0; i < ni; ++i) {
        out->f1(i,j,k) = 0.0;
        out->f2(i,j,k) = lb.field.vis_x_f2_face(i,j,k) * (*fm_x)(i,j,k)
                       + lb.field.vis_y_f2_face(i,j,k) * (*fm_y)(i,j,k)
                       + lb.field.vis_z_f2_face(i,j,k) * (*fm_z)(i,j,k);
        out->f3(i,j,k) = lb.field.vis_x_f3_face(i,j,k) * (*fm_x)(i,j,k)
                       + lb.field.vis_y_f3_face(i,j,k) * (*fm_y)(i,j,k)
                       + lb.field.vis_z_f3_face(i,j,k) * (*fm_z)(i,j,k);
        out->f4(i,j,k) = lb.field.vis_x_f4_face(i,j,k) * (*fm_x)(i,j,k)
                       + lb.field.vis_y_f4_face(i,j,k) * (*fm_y)(i,j,k)
                       + lb.field.vis_z_f4_face(i,j,k) * (*fm_z)(i,j,k);
        out->f5(i,j,k) = lb.field.vis_x_f5_face(i,j,k) * (*fm_x)(i,j,k)
                       + lb.field.vis_y_f5_face(i,j,k) * (*fm_y)(i,j,k)
                       + lb.field.vis_z_f5_face(i,j,k) * (*fm_z)(i,j,k);
    }}}
}

// ============================================================================
// 5e — Compute viscous RHS from face fluxes
// ============================================================================

// 6th-order centered difference coefficients (same as InviscidRHS)
namespace {
    constexpr Real c0_vis =  75.0 / 64.0;    //  1.171875
    constexpr Real c1_vis = -25.0 / 384.0;   // -0.065104166...
    constexpr Real c2_vis =   3.0 / 640.0;   //  0.0046875
}

inline void ViscidRHS::compute_rhs(LocalBlock& lb) {

    auto& f = lb.field;
    const Int nci = f.ni();
    const Int ncj = f.nj();
    const Int nck = f.nk();

    // Interior cell range (same as InviscidRHS)
    const Int i0 = 3, i1 = nci - 4;
    const Int j0 = 3, j1 = ncj - 4;
    const Int k0 = 3, k1 = nck - 4;

    // ---- ξ-direction: differentiate vis_xi along i ----
    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        Real dF1 = c0_vis * (f.vis_xi.f1(i+1,j,k) - f.vis_xi.f1(i,j,k))
                 + c1_vis * (f.vis_xi.f1(i+2,j,k) - f.vis_xi.f1(i-1,j,k))
                 + c2_vis * (f.vis_xi.f1(i+3,j,k) - f.vis_xi.f1(i-2,j,k));
        Real dF2 = c0_vis * (f.vis_xi.f2(i+1,j,k) - f.vis_xi.f2(i,j,k))
                 + c1_vis * (f.vis_xi.f2(i+2,j,k) - f.vis_xi.f2(i-1,j,k))
                 + c2_vis * (f.vis_xi.f2(i+3,j,k) - f.vis_xi.f2(i-2,j,k));
        Real dF3 = c0_vis * (f.vis_xi.f3(i+1,j,k) - f.vis_xi.f3(i,j,k))
                 + c1_vis * (f.vis_xi.f3(i+2,j,k) - f.vis_xi.f3(i-1,j,k))
                 + c2_vis * (f.vis_xi.f3(i+3,j,k) - f.vis_xi.f3(i-2,j,k));
        Real dF4 = c0_vis * (f.vis_xi.f4(i+1,j,k) - f.vis_xi.f4(i,j,k))
                 + c1_vis * (f.vis_xi.f4(i+2,j,k) - f.vis_xi.f4(i-1,j,k))
                 + c2_vis * (f.vis_xi.f4(i+3,j,k) - f.vis_xi.f4(i-2,j,k));
        Real dF5 = c0_vis * (f.vis_xi.f5(i+1,j,k) - f.vis_xi.f5(i,j,k))
                 + c1_vis * (f.vis_xi.f5(i+2,j,k) - f.vis_xi.f5(i-1,j,k))
                 + c2_vis * (f.vis_xi.f5(i+3,j,k) - f.vis_xi.f5(i-2,j,k));

        Real inv_J = Real(1.0) / std::abs(lb.grid.jacobian(i, j, k));
        f.rhs.rho(i,j,k)  += dF1 * inv_J;
        f.rhs.rhou(i,j,k) += dF2 * inv_J;
        f.rhs.rhov(i,j,k) += dF3 * inv_J;
        f.rhs.rhow(i,j,k) += dF4 * inv_J;
        f.rhs.rhoE(i,j,k) += dF5 * inv_J;
    }}}

    // ---- η-direction: differentiate vis_eta along j ----
    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        Real inv_J = Real(1.0) / std::abs(lb.grid.jacobian(i, j, k));
        Real dG1 = c0_vis * (f.vis_eta.f1(i,j+1,k) - f.vis_eta.f1(i,j,k))
                 + c1_vis * (f.vis_eta.f1(i,j+2,k) - f.vis_eta.f1(i,j-1,k))
                 + c2_vis * (f.vis_eta.f1(i,j+3,k) - f.vis_eta.f1(i,j-2,k));
        Real dG2 = c0_vis * (f.vis_eta.f2(i,j+1,k) - f.vis_eta.f2(i,j,k))
                 + c1_vis * (f.vis_eta.f2(i,j+2,k) - f.vis_eta.f2(i,j-1,k))
                 + c2_vis * (f.vis_eta.f2(i,j+3,k) - f.vis_eta.f2(i,j-2,k));
        Real dG3 = c0_vis * (f.vis_eta.f3(i,j+1,k) - f.vis_eta.f3(i,j,k))
                 + c1_vis * (f.vis_eta.f3(i,j+2,k) - f.vis_eta.f3(i,j-1,k))
                 + c2_vis * (f.vis_eta.f3(i,j+3,k) - f.vis_eta.f3(i,j-2,k));
        Real dG4 = c0_vis * (f.vis_eta.f4(i,j+1,k) - f.vis_eta.f4(i,j,k))
                 + c1_vis * (f.vis_eta.f4(i,j+2,k) - f.vis_eta.f4(i,j-1,k))
                 + c2_vis * (f.vis_eta.f4(i,j+3,k) - f.vis_eta.f4(i,j-2,k));
        Real dG5 = c0_vis * (f.vis_eta.f5(i,j+1,k) - f.vis_eta.f5(i,j,k))
                 + c1_vis * (f.vis_eta.f5(i,j+2,k) - f.vis_eta.f5(i,j-1,k))
                 + c2_vis * (f.vis_eta.f5(i,j+3,k) - f.vis_eta.f5(i,j-2,k));

        f.rhs.rho(i,j,k)  += dG1 * inv_J;
        f.rhs.rhou(i,j,k) += dG2 * inv_J;
        f.rhs.rhov(i,j,k) += dG3 * inv_J;
        f.rhs.rhow(i,j,k) += dG4 * inv_J;
        f.rhs.rhoE(i,j,k) += dG5 * inv_J;
    }}}

    // ---- ζ-direction: differentiate vis_zeta along k ----
    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        Real inv_J = Real(1.0) / std::abs(lb.grid.jacobian(i, j, k));
        Real dH1 = c0_vis * (f.vis_zeta.f1(i,j,k+1) - f.vis_zeta.f1(i,j,k))
                 + c1_vis * (f.vis_zeta.f1(i,j,k+2) - f.vis_zeta.f1(i,j,k-1))
                 + c2_vis * (f.vis_zeta.f1(i,j,k+3) - f.vis_zeta.f1(i,j,k-2));
        Real dH2 = c0_vis * (f.vis_zeta.f2(i,j,k+1) - f.vis_zeta.f2(i,j,k))
                 + c1_vis * (f.vis_zeta.f2(i,j,k+2) - f.vis_zeta.f2(i,j,k-1))
                 + c2_vis * (f.vis_zeta.f2(i,j,k+3) - f.vis_zeta.f2(i,j,k-2));
        Real dH3 = c0_vis * (f.vis_zeta.f3(i,j,k+1) - f.vis_zeta.f3(i,j,k))
                 + c1_vis * (f.vis_zeta.f3(i,j,k+2) - f.vis_zeta.f3(i,j,k-1))
                 + c2_vis * (f.vis_zeta.f3(i,j,k+3) - f.vis_zeta.f3(i,j,k-2));
        Real dH4 = c0_vis * (f.vis_zeta.f4(i,j,k+1) - f.vis_zeta.f4(i,j,k))
                 + c1_vis * (f.vis_zeta.f4(i,j,k+2) - f.vis_zeta.f4(i,j,k-1))
                 + c2_vis * (f.vis_zeta.f4(i,j,k+3) - f.vis_zeta.f4(i,j,k-2));
        Real dH5 = c0_vis * (f.vis_zeta.f5(i,j,k+1) - f.vis_zeta.f5(i,j,k))
                 + c1_vis * (f.vis_zeta.f5(i,j,k+2) - f.vis_zeta.f5(i,j,k-1))
                 + c2_vis * (f.vis_zeta.f5(i,j,k+3) - f.vis_zeta.f5(i,j,k-2));

        f.rhs.rho(i,j,k)  += dH1 * inv_J;
        f.rhs.rhou(i,j,k) += dH2 * inv_J;
        f.rhs.rhov(i,j,k) += dH3 * inv_J;
        f.rhs.rhow(i,j,k) += dH4 * inv_J;
        f.rhs.rhoE(i,j,k) += dH5 * inv_J;
    }}}
}
