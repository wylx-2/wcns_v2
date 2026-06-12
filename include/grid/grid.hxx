#pragma once

#include "grid.h"
#include <cmath>
#include <iostream>

inline Grid::Grid() : ni(0), nj(0), nk(0), nci(0), ncj(0), nck(0), ng(0),
                      ni_core(0), nj_core(0), nk_core(0) {}

inline Int Grid::node_count() const { return ni * nj * nk; }
inline Int Grid::cell_count() const { return nci * ncj * nck; }

inline void Grid::allocate(Int ni_, Int nj_, Int nk_, Int ng_) {
    ni = ni_; nj = nj_; nk = nk_; ng = ng_;
    nci = ni - 1; ncj = nj - 1; nck = nk - 1;
    ni_core = ni; nj_core = nj; nk_core = nk;

    node_x.allocate(ni, nj, nk);
    node_y.allocate(ni, nj, nk);
    node_z.allocate(ni, nj, nk);

    cell_x.allocate(nci, ncj, nck);
    cell_y.allocate(nci, ncj, nck);
    cell_z.allocate(nci, ncj, nck);

    cell_vol.allocate(nci, ncj, nck);
}

inline void Grid::compute_cell_centers() {
    for (Int k = 0; k < nck; ++k)
    for (Int j = 0; j < ncj; ++j)
    for (Int i = 0; i < nci; ++i) {
        // Cell (i,j,k) is bounded by vertices (i,j,k) through (i+1,j+1,k+1)
        // Center = average of 8 corner vertices
        Real cx = 0.0, cy = 0.0, cz = 0.0;
        for (Int dk = 0; dk <= 1; ++dk)
        for (Int dj = 0; dj <= 1; ++dj)
        for (Int di = 0; di <= 1; ++di) {
            cx += node_x(i+di, j+dj, k+dk);
            cy += node_y(i+di, j+dj, k+dk);
            cz += node_z(i+di, j+dj, k+dk);
        }
        cell_x(i,j,k) = cx * 0.125;
        cell_y(i,j,k) = cy * 0.125;
        cell_z(i,j,k) = cz * 0.125;
    }
}

inline void Grid::compute_cell_volumes() {
    // Uniform Cartesian approximation: cell volume = dx * dy * dz
    // For structured curvilinear grids, this is a rough approximation.
    // Proper metric-based volumes will be added when SCMM is implemented.
    for (Int k = 0; k < nck; ++k)
    for (Int j = 0; j < ncj; ++j)
    for (Int i = 0; i < nci; ++i) {
        Real dx = node_x(i+1,j,k) - node_x(i,j,k);
        Real dy = node_y(i,j+1,k) - node_y(i,j,k);
        Real dz = node_z(i,j,k+1) - node_z(i,j,k);
        cell_vol(i,j,k) = std::abs(dx * dy * dz);
    }
}

inline void Grid::print_summary() const {
    std::cout << "Grid \"" << name << "\"\n"
              << "  Core vertices (ni,nj,nk): " << ni_core << " x " << nj_core << " x " << nk_core
              << " = " << (ni_core * nj_core * nk_core) << "\n"
              << "  Total vertices (ni,nj,nk): " << ni << " x " << nj << " x " << nk
              << " = " << node_count() << "\n"
              << "  Total cells (nci,ncj,nck): " << nci << " x " << ncj << " x " << nck
              << " = " << cell_count() << "\n"
              << "  Ghost layers: " << ng << "\n"
              << "  BC patches: " << bc.num_patches() << "\n"
              << "  1-to-1 connections: " << connections.count() << "\n";
}

// ============================================================================
// Ghost layer extension
// ============================================================================

inline void Grid::extend_ghost_layers() {
    // Save original core dimensions (before extension)
    Int ni_src = ni_core;
    Int nj_src = nj_core;
    Int nk_src = nk_core;

    // --- Save old node data ---
    MultiArray3D<Real> old_nx = std::move(node_x);
    MultiArray3D<Real> old_ny = std::move(node_y);
    MultiArray3D<Real> old_nz = std::move(node_z);

    // --- Allocate extended node arrays ---
    Int nie = ni_src + 2 * ng;
    Int nje = nj_src + 2 * ng;
    Int nke = nk_src + 2 * ng;

    node_x.allocate(nie, nje, nke);
    node_y.allocate(nie, nje, nke);
    node_z.allocate(nie, nje, nke);

    // --- Copy core node data to interior (offset by ng) ---
    for (Int k = 0; k < nk_src; ++k)
    for (Int j = 0; j < nj_src; ++j)
    for (Int i = 0; i < ni_src; ++i) {
        node_x(i + ng, j + ng, k + ng) = old_nx(i, j, k);
        node_y(i + ng, j + ng, k + ng) = old_ny(i, j, k);
        node_z(i + ng, j + ng, k + ng) = old_nz(i, j, k);
    }

    // --- Update dimensions (now total = core + 2*ng) ---
    // Must be done BEFORE fill_ghost_nodes, since it uses ni/nj/nk as loop bounds
    ni = nie;  nj = nje;  nk = nke;
    nci = ni - 1;  ncj = nj - 1;  nck = nk - 1;

    // --- Fill ghost layers ---
    fill_ghost_nodes();

    // --- Allocate extended cell arrays and recompute ---
    cell_x.allocate(nci, ncj, nck);
    cell_y.allocate(nci, ncj, nck);
    cell_z.allocate(nci, ncj, nck);
    cell_vol.allocate(nci, ncj, nck);

    compute_cell_centers();
    compute_cell_volumes();

    std::cout << "  Ghost layers extended: ng=" << ng
              << ", total nodes " << ni << "x" << nj << "x" << nk
              << ", total cells " << nci << "x" << ncj << "x" << nck << "\n";
}

inline void Grid::fill_ghost_nodes() {
    // All faces start with linear extrapolation.
    // Faces with 1-to-1 connectivity are later overwritten by
    // fill_ghost_face_from_donor() (called from ParallelManager after
    // all zones have been extended, so donor interior data is available).
    for (int face = 0; face < 6; ++face) {
        fill_ghost_face_extrapolate(face);
    }
}

inline const Connectivity* Grid::find_face_connection(int face) const {
    // Use the pre-computed conn.face when available.
    for (Int c = 0; c < connections.count(); ++c) {
        const Connectivity& conn = connections[c];
        if (conn.face == face) return &conn;
    }
    return nullptr;
}

inline void Grid::fill_ghost_face_from_donor(int face, const Grid& donor,
                                              const Connectivity& conn) {
    // ========================================================================
    // Unified ghost node coordinate filling for all 1-to-1 connections.
    //
    // CGNS 1-to-1 point correspondence:
    //   Id2 = T · (Id1 - Begin1) + Begin2
    // where Id1 is the 1-based index on the current zone, Id2 on the donor,
    // Begin1 = (imin, jmin, kmin), Begin2 = (donor_imin, donor_jmin, donor_kmin),
    // and T is the signed permutation matrix encoded by transform[3].
    //
    // Ghost coordinate = donor_interior_coord - translation
    // ========================================================================

    // ---- Identify the donor face from the connection donor range ----
    int donor_face = -1;
    if (conn.donor_imin == conn.donor_imax) {
        donor_face = (conn.donor_imin == 1) ? 0 : 1;
    } else if (conn.donor_jmin == conn.donor_jmax) {
        donor_face = (conn.donor_jmin == 1) ? 2 : 3;
    } else if (conn.donor_kmin == conn.donor_kmax) {
        donor_face = (conn.donor_kmin == 1) ? 4 : 5;
    }

    int cur_dim      = face / 2;
    bool cur_is_high = (face % 2 == 1);
    int donor_dim    = donor_face / 2;
    bool donor_is_high = (donor_face % 2 == 1);

    // Direction-mismatched connections (e.g. rotated blocks) — TODO
    if (cur_dim != donor_dim) return;

    Real tx = conn.translation[0];
    Real ty = conn.translation[1];
    Real tz = conn.translation[2];

    // ---- CGNS transform helper ----
    // Maps a face point from current to donor using the CGNS 1-to-1 formula.
    //
    // Input:  (I1, J1, K1) — 1-based core indices on the current zone.
    //         These should lie within [imin:imax, jmin:jmax, kmin:kmax];
    //         values outside are clamped before the mapping.
    //
    // Output: (I2, J2, K2) — 1-based core indices on the donor zone.
    //
    // For each current-direction d with offset = Id1[d] - Begin1[d]:
    //   If transform[d] > 0: output[|t|-1] = donor_begin[|t|-1] + offset
    //   If transform[d] < 0: output[|t|-1] = donor_end[|t|-1]   - offset
    // where |t| ∈ {1,2,3} identifies the donor index direction.
    auto map_cgns_point = [&](Int I1, Int J1, Int K1,
                               Int& I2, Int& J2, Int& K2) {
        // Clamp to connection range (corner ghost regions may fall outside)
        Int i1c = std::clamp(I1, conn.imin, conn.imax);
        Int j1c = std::clamp(J1, conn.jmin, conn.jmax);
        Int k1c = std::clamp(K1, conn.kmin, conn.kmax);

        Int offset[3] = {i1c - conn.imin, j1c - conn.jmin, k1c - conn.kmin};
        Int dbeg[3]   = {conn.donor_imin, conn.donor_jmin, conn.donor_kmin};
        Int dend[3]   = {conn.donor_imax, conn.donor_jmax, conn.donor_kmax};
        Int result[3] = {0, 0, 0};

        for (int d = 0; d < 3; ++d) {
            Int t = conn.transform[d];
            int tgt = std::abs(static_cast<int>(t)) - 1;  // 0=i, 1=j, 2=k
            if (t > 0)
                result[tgt] = dbeg[tgt] + offset[d];
            else
                result[tgt] = dend[tgt] - offset[d];
        }
        I2 = result[0];  J2 = result[1];  K2 = result[2];
    };

    // ---- Helper: convert 1-based core index to 0-based extended index ----
    auto to_extended = [](Int idx_1b, Int ng_d) -> Int {
        return idx_1b - 1 + ng_d;
    };

    // ---- Helper: step inward from a donor face by d ghost layers ----
    // For a MIN donor face, interior node is at ng + d (inward).
    // For a MAX donor face, interior node is at ng + core - 1 - d (inward).
    auto step_inward = [&](Int d, Int core, Int ng_d) -> Int {
        if (donor_is_high)
            return ng_d + core - 1 - d;
        else
            return ng_d + d;
    };

    // ---- Ghost loop ----
    // For each ghost node on the current face, we:
    //  1. Project it to the face boundary → get 1-based face point
    //  2. Map the face point to the donor face via CGNS transform
    //  3. Step inward from the donor face by the ghost layer distance
    //  4. Copy the coordinate (minus translation)

    if (cur_dim == 0) {
        // ---- I-face (IMIN or IMAX) ----
        Int i_start, i_end;
        Int i_face_1b;  // 1-based i-index of the face boundary
        if (!cur_is_high) {
            i_start = 0;        i_end = ng;            // IMIN ghost
            i_face_1b = 1;                             // boundary at i=1
        } else {
            i_start = ng + ni_core;  i_end = ni;       // IMAX ghost
            i_face_1b = ni_core;                       // boundary at i=ni_core
        }

        for (Int k = 0; k < nk; ++k) {
        for (Int j = 0; j < nj; ++j) {
        for (Int i = i_start; i < i_end; ++i) {
            Int d = (!cur_is_high) ? (ng - i) : (i - (ng + ni_core) + 1);

            // Current face point in 1-based core
            Int J1 = j - ng + 1;
            Int K1 = k - ng + 1;

            // Map to donor face
            Int I2, J2, K2;
            map_cgns_point(i_face_1b, J1, K1, I2, J2, K2);

            // Donor interior: normal = step inward, tangential = mapped face
            Int di = step_inward(d, donor.ni_core, donor.ng);
            Int dj = to_extended(J2, donor.ng);
            Int dk = to_extended(K2, donor.ng);

            node_x(i, j, k) = donor.node_x(di, dj, dk) - tx;
            node_y(i, j, k) = donor.node_y(di, dj, dk) - ty;
            node_z(i, j, k) = donor.node_z(di, dj, dk) - tz;
        }}}
    } else if (cur_dim == 1) {
        // ---- J-face (JMIN or JMAX) ----
        Int j_start, j_end;
        Int j_face_1b;
        if (!cur_is_high) {
            j_start = 0;        j_end = ng;
            j_face_1b = 1;
        } else {
            j_start = ng + nj_core;  j_end = nj;
            j_face_1b = nj_core;
        }

        for (Int k = 0; k < nk; ++k) {
        for (Int j = j_start; j < j_end; ++j) {
        for (Int i = 0; i < ni; ++i) {
            Int d = (!cur_is_high) ? (ng - j) : (j - (ng + nj_core) + 1);

            Int I1 = i - ng + 1;
            Int K1 = k - ng + 1;

            Int I2, J2, K2;
            map_cgns_point(I1, j_face_1b, K1, I2, J2, K2);

            Int di = to_extended(I2, donor.ng);
            Int dj = step_inward(d, donor.nj_core, donor.ng);
            Int dk = to_extended(K2, donor.ng);

            node_x(i, j, k) = donor.node_x(di, dj, dk) - tx;
            node_y(i, j, k) = donor.node_y(di, dj, dk) - ty;
            node_z(i, j, k) = donor.node_z(di, dj, dk) - tz;
        }}}
    } else { // cur_dim == 2
        // ---- K-face (KMIN or KMAX) ----
        Int k_start, k_end;
        Int k_face_1b;
        if (!cur_is_high) {
            k_start = 0;        k_end = ng;
            k_face_1b = 1;
        } else {
            k_start = ng + nk_core;  k_end = nk;
            k_face_1b = nk_core;
        }

        for (Int k = k_start; k < k_end; ++k) {
        for (Int j = 0; j < nj; ++j) {
        for (Int i = 0; i < ni; ++i) {
            Int d = (!cur_is_high) ? (ng - k) : (k - (ng + nk_core) + 1);

            Int I1 = i - ng + 1;
            Int J1 = j - ng + 1;

            Int I2, J2, K2;
            map_cgns_point(I1, J1, k_face_1b, I2, J2, K2);

            Int di = to_extended(I2, donor.ng);
            Int dj = to_extended(J2, donor.ng);
            Int dk = step_inward(d, donor.nk_core, donor.ng);

            node_x(i, j, k) = donor.node_x(di, dj, dk) - tx;
            node_y(i, j, k) = donor.node_y(di, dj, dk) - ty;
            node_z(i, j, k) = donor.node_z(di, dj, dk) - tz;
        }}}
    }
}

inline void Grid::fill_ghost_face_extrapolate(int face) {
    // Linear extrapolation: x(-d) = 2*x(-(d-1)) - x(-(d-2))
    // Applied layer-by-layer outward from the physical face.

    int dim = face / 2;            // 0=i, 1=j, 2=k
    bool is_high = (face % 2 == 1); // MAX side

    for (Int layer = 1; layer <= ng; ++layer) {
        if (dim == 0) {
            Int ghost_i, base_i, prev_i;
            if (!is_high) {
                // IMIN: ghost at ng - layer, base at ng (face), prev at ng + 1
                ghost_i = ng - layer;
                base_i  = ng - layer + 1;   // one layer closer to face
                prev_i  = ng - layer + 2;   // two layers closer to face
            } else {
                // IMAX: ghost at ng + ni_core - 1 + layer
                ghost_i = ng + ni_core - 1 + layer;
                base_i  = ng + ni_core - 1 + layer - 1;
                prev_i  = ng + ni_core - 1 + layer - 2;
            }
            for (Int k = 0; k < nk; ++k) {
            for (Int j = 0; j < nj; ++j) {
                node_x(ghost_i, j, k) = 2.0 * node_x(base_i, j, k) - node_x(prev_i, j, k);
                node_y(ghost_i, j, k) = 2.0 * node_y(base_i, j, k) - node_y(prev_i, j, k);
                node_z(ghost_i, j, k) = 2.0 * node_z(base_i, j, k) - node_z(prev_i, j, k);
            }}
        } else if (dim == 1) {
            Int ghost_j, base_j, prev_j;
            if (!is_high) {
                ghost_j = ng - layer;
                base_j  = ng - layer + 1;
                prev_j  = ng - layer + 2;
            } else {
                ghost_j = ng + nj_core - 1 + layer;
                base_j  = ng + nj_core - 1 + layer - 1;
                prev_j  = ng + nj_core - 1 + layer - 2;
            }
            // ghost_j is fixed; iterate over i,k plane
            for (Int k = 0; k < nk; ++k) {
            for (Int i = 0; i < ni; ++i) {
                node_x(i, ghost_j, k) = 2.0 * node_x(i, base_j, k) - node_x(i, prev_j, k);
                node_y(i, ghost_j, k) = 2.0 * node_y(i, base_j, k) - node_y(i, prev_j, k);
                node_z(i, ghost_j, k) = 2.0 * node_z(i, base_j, k) - node_z(i, prev_j, k);
            }}
        } else {
            Int ghost_k, base_k, prev_k;
            if (!is_high) {
                ghost_k = ng - layer;
                base_k  = ng - layer + 1;
                prev_k  = ng - layer + 2;
            } else {
                ghost_k = ng + nk_core - 1 + layer;
                base_k  = ng + nk_core - 1 + layer - 1;
                prev_k  = ng + nk_core - 1 + layer - 2;
            }
            // ghost_k is fixed; iterate over i,j plane
            for (Int j = 0; j < nj; ++j) {
            for (Int i = 0; i < ni; ++i) {
                node_x(i, j, ghost_k) = 2.0 * node_x(i, j, base_k) - node_x(i, j, prev_k);
                node_y(i, j, ghost_k) = 2.0 * node_y(i, j, base_k) - node_y(i, j, prev_k);
                node_z(i, j, ghost_k) = 2.0 * node_z(i, j, base_k) - node_z(i, j, prev_k);
            }}
        }
    }
}

// ============================================================================
// Metric coefficient computation (SCMM)
// ============================================================================

inline void Grid::build_face_connected(bool face_connected[6]) const {
    for (int f = 0; f < 6; ++f) {
        face_connected[f] = (find_face_connection(f) != nullptr);
    }
}

inline void Grid::compute_coord_deriv(int dir, Real dh, const bool fp[6],
                                       MultiArray3D<Real>& dx_d,
                                       MultiArray3D<Real>& dy_d,
                                       MultiArray3D<Real>& dz_d) {
    InterpDiff::derivative(cell_x, dx_d, dir, dh, ng, fp);
    InterpDiff::derivative(cell_y, dy_d, dir, dh, ng, fp);
    InterpDiff::derivative(cell_z, dz_d, dir, dh, ng, fp);
}

inline void Grid::compute_metrics() {
    if (metrics_type == "uniform") {
        compute_metrics_uniform();
        return;
    }

    // ---- 1. Determine periodic faces ----
    bool fp[6];
    build_face_connected(fp);

    const Real dh = 1.0;  // computational-space grid spacing

    // ---- 2. Allocate cell-center metric arrays ----
    met_xi_x.allocate(nci, ncj, nck);
    met_xi_y.allocate(nci, ncj, nck);
    met_xi_z.allocate(nci, ncj, nck);
    met_eta_x.allocate(nci, ncj, nck);
    met_eta_y.allocate(nci, ncj, nck);
    met_eta_z.allocate(nci, ncj, nck);
    met_zeta_x.allocate(nci, ncj, nck);
    met_zeta_y.allocate(nci, ncj, nck);
    met_zeta_z.allocate(nci, ncj, nck);
    jacobian.allocate(nci, ncj, nck);

    // ---- 3. Compute coordinate derivatives at cell centers ----
    // x_ξ, x_η, x_ζ,  y_ξ, y_η, y_ζ,  z_ξ, z_η, z_ζ
    MultiArray3D<Real> x_xi(nci,ncj,nck), x_eta(nci,ncj,nck), x_zeta(nci,ncj,nck);
    MultiArray3D<Real> y_xi(nci,ncj,nck), y_eta(nci,ncj,nck), y_zeta(nci,ncj,nck);
    MultiArray3D<Real> z_xi(nci,ncj,nck), z_eta(nci,ncj,nck), z_zeta(nci,ncj,nck);

    compute_coord_deriv(0, dh, fp, x_xi,  y_xi,  z_xi);   // ξ derivatives
    compute_coord_deriv(1, dh, fp, x_eta, y_eta, z_eta);  // η derivatives
    compute_coord_deriv(2, dh, fp, x_zeta,y_zeta,z_zeta); // ζ derivatives

    // ---- 4. Temporary arrays for products and derivatives ----
    MultiArray3D<Real> t_prod(nci, ncj, nck);
    MultiArray3D<Real> t_deriv(nci, ncj, nck);

    // Lambda: pointwise product  t = a * b
    auto mul = [&](const MultiArray3D<Real>& a, const MultiArray3D<Real>& b) {
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i)
            t_prod(i,j,k) = a(i,j,k) * b(i,j,k);
    };

    // Lambda: accumulate  met += coef * deriv(t_prod, dir)
    auto accum_deriv = [&](MultiArray3D<Real>& met, Real coef, int d) {
        InterpDiff::derivative(t_prod, t_deriv, d, dh, ng, fp);
        if (coef == 1.0) {
            for (Int k = 0; k < nck; ++k)
            for (Int j = 0; j < ncj; ++j)
            for (Int i = 0; i < nci; ++i)
                met(i,j,k) += t_deriv(i,j,k);
        } else if (coef == -1.0) {
            for (Int k = 0; k < nck; ++k)
            for (Int j = 0; j < ncj; ++j)
            for (Int i = 0; i < nci; ++i)
                met(i,j,k) -= t_deriv(i,j,k);
        } else {
            for (Int k = 0; k < nck; ++k)
            for (Int j = 0; j < ncj; ++j)
            for (Int i = 0; i < nci; ++i)
                met(i,j,k) += coef * t_deriv(i,j,k);
        }
    };

    // ================================================================
    // 5. Compute ξ̂ metrics
    // ================================================================

    // --- ξ̂_x = 0.5*[(z·y_η)_ζ + (y·z_ζ)_η - (z·y_ζ)_η - (y·z_η)_ζ] ---
    met_xi_x.fill(0.0);
    mul(cell_z, y_eta);   accum_deriv(met_xi_x,  0.5, 2);  // +(z·y_η)_ζ
    mul(cell_y, z_zeta);  accum_deriv(met_xi_x,  0.5, 1);  // +(y·z_ζ)_η
    mul(cell_z, y_zeta);  accum_deriv(met_xi_x, -0.5, 1);  // -(z·y_ζ)_η
    mul(cell_y, z_eta);   accum_deriv(met_xi_x, -0.5, 2);  // -(y·z_η)_ζ

    // --- ξ̂_y = 0.5*[(x·z_η)_ζ + (z·x_ζ)_η - (x·z_ζ)_η - (z·x_η)_ζ] ---
    met_xi_y.fill(0.0);
    mul(cell_x, z_eta);   accum_deriv(met_xi_y,  0.5, 2);  // +(x·z_η)_ζ
    mul(cell_z, x_zeta);  accum_deriv(met_xi_y,  0.5, 1);  // +(z·x_ζ)_η
    mul(cell_x, z_zeta);  accum_deriv(met_xi_y, -0.5, 1);  // -(x·z_ζ)_η
    mul(cell_z, x_eta);   accum_deriv(met_xi_y, -0.5, 2);  // -(z·x_η)_ζ

    // --- ξ̂_z = 0.5*[(y·x_η)_ζ + (x·y_ζ)_η - (y·x_ζ)_η - (x·y_η)_ζ] ---
    met_xi_z.fill(0.0);
    mul(cell_y, x_eta);   accum_deriv(met_xi_z,  0.5, 2);  // +(y·x_η)_ζ
    mul(cell_x, y_zeta);  accum_deriv(met_xi_z,  0.5, 1);  // +(x·y_ζ)_η
    mul(cell_y, x_zeta);  accum_deriv(met_xi_z, -0.5, 1);  // -(y·x_ζ)_η
    mul(cell_x, y_eta);   accum_deriv(met_xi_z, -0.5, 2);  // -(x·y_η)_ζ

    // ================================================================
    // 6. Compute η̂ metrics
    // ================================================================

    // --- η̂_x = 0.5*[(z·y_ζ)_ξ + (y·z_ξ)_ζ - (z·y_ξ)_ζ - (y·z_ζ)_ξ] ---
    met_eta_x.fill(0.0);
    mul(cell_z, y_zeta);  accum_deriv(met_eta_x,  0.5, 0);  // +(z·y_ζ)_ξ
    mul(cell_y, z_xi);    accum_deriv(met_eta_x,  0.5, 2);  // +(y·z_ξ)_ζ
    mul(cell_z, y_xi);    accum_deriv(met_eta_x, -0.5, 2);  // -(z·y_ξ)_ζ
    mul(cell_y, z_zeta);  accum_deriv(met_eta_x, -0.5, 0);  // -(y·z_ζ)_ξ

    // --- η̂_y = 0.5*[(x·z_ζ)_ξ + (z·x_ξ)_ζ - (x·z_ξ)_ζ - (z·x_ζ)_ξ] ---
    met_eta_y.fill(0.0);
    mul(cell_x, z_zeta);  accum_deriv(met_eta_y,  0.5, 0);  // +(x·z_ζ)_ξ
    mul(cell_z, x_xi);    accum_deriv(met_eta_y,  0.5, 2);  // +(z·x_ξ)_ζ
    mul(cell_x, z_xi);    accum_deriv(met_eta_y, -0.5, 2);  // -(x·z_ξ)_ζ
    mul(cell_z, x_zeta);  accum_deriv(met_eta_y, -0.5, 0);  // -(z·x_ζ)_ξ

    // --- η̂_z = 0.5*[(y·x_ζ)_ξ + (x·y_ξ)_ζ - (y·x_ξ)_ζ - (x·y_ζ)_ξ] ---
    met_eta_z.fill(0.0);
    mul(cell_y, x_zeta);  accum_deriv(met_eta_z,  0.5, 0);  // +(y·x_ζ)_ξ
    mul(cell_x, y_xi);    accum_deriv(met_eta_z,  0.5, 2);  // +(x·y_ξ)_ζ
    mul(cell_y, x_xi);    accum_deriv(met_eta_z, -0.5, 2);  // -(y·x_ξ)_ζ
    mul(cell_x, y_zeta);  accum_deriv(met_eta_z, -0.5, 0);  // -(x·y_ζ)_ξ

    // ================================================================
    // 7. Compute ζ̂ metrics
    // ================================================================

    // --- ζ̂_x = 0.5*[(z·y_ξ)_η + (y·z_η)_ξ - (z·y_η)_ξ - (y·z_ξ)_η] ---
    met_zeta_x.fill(0.0);
    mul(cell_z, y_xi);    accum_deriv(met_zeta_x,  0.5, 1);  // +(z·y_ξ)_η
    mul(cell_y, z_eta);   accum_deriv(met_zeta_x,  0.5, 0);  // +(y·z_η)_ξ
    mul(cell_z, y_eta);   accum_deriv(met_zeta_x, -0.5, 0);  // -(z·y_η)_ξ
    mul(cell_y, z_xi);    accum_deriv(met_zeta_x, -0.5, 1);  // -(y·z_ξ)_η

    // --- ζ̂_y = 0.5*[(x·z_ξ)_η + (z·x_η)_ξ - (x·z_η)_ξ - (z·x_ξ)_η] ---
    met_zeta_y.fill(0.0);
    mul(cell_x, z_xi);    accum_deriv(met_zeta_y,  0.5, 1);  // +(x·z_ξ)_η
    mul(cell_z, x_eta);   accum_deriv(met_zeta_y,  0.5, 0);  // +(z·x_η)_ξ
    mul(cell_x, z_eta);   accum_deriv(met_zeta_y, -0.5, 0);  // -(x·z_η)_ξ
    mul(cell_z, x_xi);    accum_deriv(met_zeta_y, -0.5, 1);  // -(z·x_ξ)_η

    // --- ζ̂_z = 0.5*[(y·x_ξ)_η + (x·y_η)_ξ - (y·x_η)_ξ - (x·y_ξ)_η] ---
    met_zeta_z.fill(0.0);
    mul(cell_y, x_xi);    accum_deriv(met_zeta_z,  0.5, 1);  // +(y·x_ξ)_η
    mul(cell_x, y_eta);   accum_deriv(met_zeta_z,  0.5, 0);  // +(x·y_η)_ξ
    mul(cell_y, x_eta);   accum_deriv(met_zeta_z, -0.5, 0);  // -(y·x_η)_ξ
    mul(cell_x, y_xi);    accum_deriv(met_zeta_z, -0.5, 1);  // -(x·y_ξ)_η

    // ================================================================
    // 8. Compute Jacobian (SCMM form)
    //    J = (1/3) * [ (S_ξ)_ξ + (S_η)_η + (S_ζ)_ζ ]
    //    where S_dir = x·met_dir_x + y·met_dir_y + z·met_dir_z
    // ================================================================

    // S_ξ and (S_ξ)_ξ
    for (Int k = 0; k < nck; ++k)
    for (Int j = 0; j < ncj; ++j)
    for (Int i = 0; i < nci; ++i)
        t_prod(i,j,k) = cell_x(i,j,k) * met_xi_x(i,j,k)
                      + cell_y(i,j,k) * met_xi_y(i,j,k)
                      + cell_z(i,j,k) * met_xi_z(i,j,k);
    InterpDiff::derivative(t_prod, jacobian, 0, dh, ng, fp);  // jacobian = (S_ξ)_ξ

    // S_η and (S_η)_η
    for (Int k = 0; k < nck; ++k)
    for (Int j = 0; j < ncj; ++j)
    for (Int i = 0; i < nci; ++i)
        t_prod(i,j,k) = cell_x(i,j,k) * met_eta_x(i,j,k)
                      + cell_y(i,j,k) * met_eta_y(i,j,k)
                      + cell_z(i,j,k) * met_eta_z(i,j,k);
    InterpDiff::derivative(t_prod, t_deriv, 1, dh, ng, fp);
    for (Int k = 0; k < nck; ++k)
    for (Int j = 0; j < ncj; ++j)
    for (Int i = 0; i < nci; ++i)
        jacobian(i,j,k) += t_deriv(i,j,k);

    // S_ζ and (S_ζ)_ζ
    for (Int k = 0; k < nck; ++k)
    for (Int j = 0; j < ncj; ++j)
    for (Int i = 0; i < nci; ++i)
        t_prod(i,j,k) = cell_x(i,j,k) * met_zeta_x(i,j,k)
                      + cell_y(i,j,k) * met_zeta_y(i,j,k)
                      + cell_z(i,j,k) * met_zeta_z(i,j,k);
    InterpDiff::derivative(t_prod, t_deriv, 2, dh, ng, fp);
    for (Int k = 0; k < nck; ++k)
    for (Int j = 0; j < ncj; ++j)
    for (Int i = 0; i < nci; ++i)
        jacobian(i,j,k) += t_deriv(i,j,k);

    // J = sum/3  (sum = (S_ξ)_ξ + (S_η)_η + (S_ζ)_ζ = 3J)
    for (Int k = 0; k < nck; ++k)
    for (Int j = 0; j < ncj; ++j)
    for (Int i = 0; i < nci; ++i)
        jacobian(i,j,k) = jacobian(i,j,k) / 3.0;

    std::cout << "  SCMM metrics computed: 9 terms + Jacobian at cell centers\n";
}

// ============================================================================
// Analytical uniform Cartesian metrics
// ============================================================================

inline void Grid::compute_metrics_uniform() {
    // Infer dx, dy, dz from interior cell centers.
    // Use cells away from boundaries (interior region) for robustness.
    Int i0 = ng + nci / 4;
    Int j0 = ng + ncj / 4;
    Int k0 = ng + nck / 4;

    Real dx = cell_x(i0+1, j0, k0) - cell_x(i0, j0, k0);
    Real dy = cell_y(i0, j0+1, k0) - cell_y(i0, j0, k0);
    Real dz = cell_z(i0, j0, k0+1) - cell_z(i0, j0, k0);

    std::cout << "  Analytical uniform metrics: dx=" << dx
              << ", dy=" << dy << ", dz=" << dz << "\n";

    // Allocate cell-center metric arrays
    met_xi_x.allocate(nci, ncj, nck);
    met_xi_y.allocate(nci, ncj, nck);
    met_xi_z.allocate(nci, ncj, nck);
    met_eta_x.allocate(nci, ncj, nck);
    met_eta_y.allocate(nci, ncj, nck);
    met_eta_z.allocate(nci, ncj, nck);
    met_zeta_x.allocate(nci, ncj, nck);
    met_zeta_y.allocate(nci, ncj, nck);
    met_zeta_z.allocate(nci, ncj, nck);
    jacobian.allocate(nci, ncj, nck);

    // For uniform Cartesian grid, metrics are constant:
    //   S_ξ = (dy*dz, 0, 0)     — face area in x-direction
    //   S_η = (0, dx*dz, 0)     — face area in y-direction
    //   S_ζ = (0, 0, dx*dy)     — face area in z-direction
    //   J   = dx*dy*dz          — cell volume
    Real xi_area  = dy * dz;   // |S_ξ|
    Real eta_area = dx * dz;   // |S_η|
    Real zeta_area = dx * dy;  // |S_ζ|
    Real J = dx * dy * dz;

    for (Int k = 0; k < nck; ++k) {
    for (Int j = 0; j < ncj; ++j) {
    for (Int i = 0; i < nci; ++i) {
        // ξ-direction metrics (x-normal faces)
        met_xi_x(i,j,k) = xi_area;
        met_xi_y(i,j,k) = 0.0;
        met_xi_z(i,j,k) = 0.0;

        // η-direction metrics (y-normal faces)
        met_eta_x(i,j,k) = 0.0;
        met_eta_y(i,j,k) = eta_area;
        met_eta_z(i,j,k) = 0.0;

        // ζ-direction metrics (z-normal faces)
        met_zeta_x(i,j,k) = 0.0;
        met_zeta_y(i,j,k) = 0.0;
        met_zeta_z(i,j,k) = zeta_area;

        // Jacobian
        jacobian(i,j,k) = J;
    }}}

    std::cout << "  Uniform metrics: |S_xi|=" << xi_area
              << ", |S_eta|=" << eta_area
              << ", |S_zeta|=" << zeta_area
              << ", J=" << J << "\n";
}

// ============================================================================
// Face metric interpolation
// ============================================================================

inline void Grid::compute_face_metrics() {
    bool fp[6];
    build_face_connected(fp);

    if (metrics_type == "uniform") {
        // Uniform grid: face metrics equal cell-center metrics (all constant).
        // Infer dx,dy,dz from interior cells (same as compute_metrics_uniform)
        Int i0 = ng + nci / 4;
        Int j0 = ng + ncj / 4;
        Int k0 = ng + nck / 4;
        Real dx = cell_x(i0+1, j0, k0) - cell_x(i0, j0, k0);
        Real dy = cell_y(i0, j0+1, k0) - cell_y(i0, j0, k0);
        Real dz = cell_z(i0, j0, k0+1) - cell_z(i0, j0, k0);
        Real xi_area  = dy * dz;
        Real eta_area = dx * dz;
        Real zeta_area = dx * dy;

        // X-faces (i+1/2): size (nci+1) x ncj x nck
        face_xi_x.allocate(nci + 1, ncj, nck);
        face_xi_y.allocate(nci + 1, ncj, nck);
        face_xi_z.allocate(nci + 1, ncj, nck);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci + 1; ++i) {
            face_xi_x(i,j,k) = xi_area;
            face_xi_y(i,j,k) = 0.0;
            face_xi_z(i,j,k) = 0.0;
        }

        // Y-faces (j+1/2): size nci x (ncj+1) x nck
        face_eta_x.allocate(nci, ncj + 1, nck);
        face_eta_y.allocate(nci, ncj + 1, nck);
        face_eta_z.allocate(nci, ncj + 1, nck);
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj + 1; ++j)
        for (Int i = 0; i < nci; ++i) {
            face_eta_x(i,j,k) = 0.0;
            face_eta_y(i,j,k) = eta_area;
            face_eta_z(i,j,k) = 0.0;
        }

        // Z-faces (k+1/2): size nci x ncj x (nck+1)
        face_zeta_x.allocate(nci, ncj, nck + 1);
        face_zeta_y.allocate(nci, ncj, nck + 1);
        face_zeta_z.allocate(nci, ncj, nck + 1);
        for (Int k = 0; k < nck + 1; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i) {
            face_zeta_x(i,j,k) = 0.0;
            face_zeta_y(i,j,k) = 0.0;
            face_zeta_z(i,j,k) = zeta_area;
        }

        std::cout << "  Uniform face metrics set directly (no interpolation)\n";
        return;
    }

    // --- X-faces (i+1/2): size (nci+1) x ncj x nck ---
    face_xi_x.allocate(nci + 1, ncj, nck);
    face_xi_y.allocate(nci + 1, ncj, nck);
    face_xi_z.allocate(nci + 1, ncj, nck);

    InterpDiff::interp_to_faces(met_xi_x, face_xi_x, 0, ng, fp);
    InterpDiff::interp_to_faces(met_xi_y, face_xi_y, 0, ng, fp);
    InterpDiff::interp_to_faces(met_xi_z, face_xi_z, 0, ng, fp);

    // --- Y-faces (j+1/2): size nci x (ncj+1) x nck ---
    face_eta_x.allocate(nci, ncj + 1, nck);
    face_eta_y.allocate(nci, ncj + 1, nck);
    face_eta_z.allocate(nci, ncj + 1, nck);

    InterpDiff::interp_to_faces(met_eta_x, face_eta_x, 1, ng, fp);
    InterpDiff::interp_to_faces(met_eta_y, face_eta_y, 1, ng, fp);
    InterpDiff::interp_to_faces(met_eta_z, face_eta_z, 1, ng, fp);

    // --- Z-faces (k+1/2): size nci x ncj x (nck+1) ---
    face_zeta_x.allocate(nci, ncj, nck + 1);
    face_zeta_y.allocate(nci, ncj, nck + 1);
    face_zeta_z.allocate(nci, ncj, nck + 1);

    InterpDiff::interp_to_faces(met_zeta_x, face_zeta_x, 2, ng, fp);
    InterpDiff::interp_to_faces(met_zeta_y, face_zeta_y, 2, ng, fp);
    InterpDiff::interp_to_faces(met_zeta_z, face_zeta_z, 2, ng, fp);

    std::cout << "  Face metrics interpolated: "
              << "X-faces (" << face_xi_x.ni() << "x" << face_xi_x.nj() << "x" << face_xi_x.nk() << "), "
              << "Y-faces (" << face_eta_x.ni() << "x" << face_eta_x.nj() << "x" << face_eta_x.nk() << "), "
              << "Z-faces (" << face_zeta_x.ni() << "x" << face_zeta_x.nj() << "x" << face_zeta_x.nk() << ")\n";
}

// ============================================================================
// Extract metrics from donor grid
// ============================================================================

inline void Grid::extract_metrics_from(const Grid& donor, Int ci0, Int cj0, Int ck0) {
    // ---- Cell-center metrics (nci x ncj x nck) ----
    met_xi_x.allocate(nci, ncj, nck);
    met_xi_y.allocate(nci, ncj, nck);
    met_xi_z.allocate(nci, ncj, nck);
    met_eta_x.allocate(nci, ncj, nck);
    met_eta_y.allocate(nci, ncj, nck);
    met_eta_z.allocate(nci, ncj, nck);
    met_zeta_x.allocate(nci, ncj, nck);
    met_zeta_y.allocate(nci, ncj, nck);
    met_zeta_z.allocate(nci, ncj, nck);
    jacobian.allocate(nci, ncj, nck);

    for (Int k = 0; k < nck; ++k) {
    for (Int j = 0; j < ncj; ++j) {
    for (Int i = 0; i < nci; ++i) {
        Int di = ci0 + i;
        Int dj = cj0 + j;
        Int dk = ck0 + k;

        met_xi_x(i,j,k)   = donor.met_xi_x(di, dj, dk);
        met_xi_y(i,j,k)   = donor.met_xi_y(di, dj, dk);
        met_xi_z(i,j,k)   = donor.met_xi_z(di, dj, dk);
        met_eta_x(i,j,k)  = donor.met_eta_x(di, dj, dk);
        met_eta_y(i,j,k)  = donor.met_eta_y(di, dj, dk);
        met_eta_z(i,j,k)  = donor.met_eta_z(di, dj, dk);
        met_zeta_x(i,j,k) = donor.met_zeta_x(di, dj, dk);
        met_zeta_y(i,j,k) = donor.met_zeta_y(di, dj, dk);
        met_zeta_z(i,j,k) = donor.met_zeta_z(di, dj, dk);
        jacobian(i,j,k)   = donor.jacobian(di, dj, dk);
    }}}

    // ---- Face metrics ----
    // X-faces: (nci+1) x ncj x nck
    face_xi_x.allocate(nci + 1, ncj, nck);
    face_xi_y.allocate(nci + 1, ncj, nck);
    face_xi_z.allocate(nci + 1, ncj, nck);
    for (Int k = 0; k < nck; ++k) {
    for (Int j = 0; j < ncj; ++j) {
    for (Int i = 0; i < nci + 1; ++i) {
        Int di = ci0 + i;
        Int dj = cj0 + j;
        Int dk = ck0 + k;
        face_xi_x(i,j,k) = donor.face_xi_x(di, dj, dk);
        face_xi_y(i,j,k) = donor.face_xi_y(di, dj, dk);
        face_xi_z(i,j,k) = donor.face_xi_z(di, dj, dk);
    }}}

    // Y-faces: nci x (ncj+1) x nck
    face_eta_x.allocate(nci, ncj + 1, nck);
    face_eta_y.allocate(nci, ncj + 1, nck);
    face_eta_z.allocate(nci, ncj + 1, nck);
    for (Int k = 0; k < nck; ++k) {
    for (Int j = 0; j < ncj + 1; ++j) {
    for (Int i = 0; i < nci; ++i) {
        Int di = ci0 + i;
        Int dj = cj0 + j;
        Int dk = ck0 + k;
        face_eta_x(i,j,k) = donor.face_eta_x(di, dj, dk);
        face_eta_y(i,j,k) = donor.face_eta_y(di, dj, dk);
        face_eta_z(i,j,k) = donor.face_eta_z(di, dj, dk);
    }}}

    // Z-faces: nci x ncj x (nck+1)
    face_zeta_x.allocate(nci, ncj, nck + 1);
    face_zeta_y.allocate(nci, ncj, nck + 1);
    face_zeta_z.allocate(nci, ncj, nck + 1);
    for (Int k = 0; k < nck + 1; ++k) {
    for (Int j = 0; j < ncj; ++j) {
    for (Int i = 0; i < nci; ++i) {
        Int di = ci0 + i;
        Int dj = cj0 + j;
        Int dk = ck0 + k;
        face_zeta_x(i,j,k) = donor.face_zeta_x(di, dj, dk);
        face_zeta_y(i,j,k) = donor.face_zeta_y(di, dj, dk);
        face_zeta_z(i,j,k) = donor.face_zeta_z(di, dj, dk);
    }}}

    // Also copy cell_vol (used by history monitor)
    cell_vol.allocate(nci, ncj, nck);
    for (Int k = 0; k < nck; ++k) {
    for (Int j = 0; j < ncj; ++j) {
    for (Int i = 0; i < nci; ++i) {
        cell_vol(i,j,k) = donor.cell_vol(ci0 + i, cj0 + j, ck0 + k);
    }}}

    std::cout << "  Metrics extracted from donor at offset ("
              << ci0 << "," << cj0 << "," << ck0 << ")\n";
}
