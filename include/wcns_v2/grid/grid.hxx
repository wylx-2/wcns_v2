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
    for (int face = 0; face < 6; ++face) {
        const Connectivity* conn = find_periodic_connection(face);
        if (conn) {
            fill_ghost_face_periodic(face, *conn);
        } else {
            fill_ghost_face_extrapolate(face);
        }
    }
}

inline const Connectivity* Grid::find_periodic_connection(int face) const {
    // face: 0=IMIN, 1=IMAX, 2=JMIN, 3=JMAX, 4=KMIN, 5=KMAX
    // A connection covers a face when the face's normal dimension is
    // collapsed in the point range and the constant index lies at the
    // corresponding boundary (1 for MIN, ni_core for MAX, etc.).
    // CGNS range indices are 1-based.

    for (Int c = 0; c < connections.count(); ++c) {
        const Connectivity& conn = connections[c];
        if (!conn.is_periodic) continue;

        switch (face) {
        case 0: // IMIN
            if (conn.imin == conn.imax && conn.imin == 1) return &conn;
            break;
        case 1: // IMAX
            if (conn.imin == conn.imax && conn.imin == ni_core) return &conn;
            break;
        case 2: // JMIN
            if (conn.jmin == conn.jmax && conn.jmin == 1) return &conn;
            break;
        case 3: // JMAX
            if (conn.jmin == conn.jmax && conn.jmin == nj_core) return &conn;
            break;
        case 4: // KMIN
            if (conn.kmin == conn.kmax && conn.kmin == 1) return &conn;
            break;
        case 5: // KMAX
            if (conn.kmin == conn.kmax && conn.kmin == nk_core) return &conn;
            break;
        }
    }
    return nullptr;
}

inline const Connectivity* Grid::find_face_connection(int face) const {
    // Same logic as find_periodic_connection, but accepts ALL 1-to-1
    // connections (both periodic and inter-zone interfaces).
    for (Int c = 0; c < connections.count(); ++c) {
        const Connectivity& conn = connections[c];

        switch (face) {
        case 0: // IMIN
            if (conn.imin == conn.imax && conn.imin == 1) return &conn;
            break;
        case 1: // IMAX
            if (conn.imin == conn.imax && conn.imin == ni_core) return &conn;
            break;
        case 2: // JMIN
            if (conn.jmin == conn.jmax && conn.jmin == 1) return &conn;
            break;
        case 3: // JMAX
            if (conn.jmin == conn.jmax && conn.jmin == nj_core) return &conn;
            break;
        case 4: // KMIN
            if (conn.kmin == conn.kmax && conn.kmin == 1) return &conn;
            break;
        case 5: // KMAX
            if (conn.kmin == conn.kmax && conn.kmin == nk_core) return &conn;
            break;
        }
    }
    return nullptr;
}

inline void Grid::fix_interface_ghost(int face, const Grid& donor,
                                       const Connectivity& conn) {
    // Identify the donor face from the connection donor range.
    // The collapsed dimension gives the face direction; the constant
    // index (1 for MIN, donor.ni_core for MAX) gives the side.
    int donor_face = -1;
    if (conn.donor_imin == conn.donor_imax) {
        donor_face = (conn.donor_imin == 1) ? 0 : 1;  // IMIN or IMAX
    } else if (conn.donor_jmin == conn.donor_jmax) {
        donor_face = (conn.donor_jmin == 1) ? 2 : 3;  // JMIN or JMAX
    } else if (conn.donor_kmin == conn.donor_kmax) {
        donor_face = (conn.donor_kmin == 1) ? 4 : 5;  // KMIN or KMAX
    }

    int cur_dim   = face / 2;
    bool cur_is_high  = (face % 2 == 1);
    int donor_dim = donor_face / 2;
    bool donor_is_high = (donor_face % 2 == 1);

    // Periodic translation — for inter-zone periodic connections,
    // ghost = donor_node - translation (same convention as fill_ghost_face_periodic)
    Real tx = conn.translation[0];
    Real ty = conn.translation[1];
    Real tz = conn.translation[2];

    if (cur_dim != donor_dim) {
        // Direction-mismatched connections (e.g. rotated blocks)
        // are not yet supported — fall back to extrapolation value.
        return;
    }

    // Ghost node copying: the ghost on this zone's face corresponds to
    // interior nodes of the donor zone.  For a MIN-side face, donor
    // interior is at donor.ng + d (inward).  For a MAX-side face,
    // donor interior is at donor.ng + donor.ni_core - 1 - d (inward).

    if (cur_dim == 0) {
        Int i_start, i_end;
        if (!cur_is_high) {
            // IMIN ghost: i = 0 .. ng-1
            i_start = 0;
            i_end   = ng;
        } else {
            // IMAX ghost: i = ng+ni_core .. ni-1
            i_start = ng + ni_core;
            i_end   = ni;
        }

        for (Int k = 0; k < nk; ++k) {
        for (Int j = 0; j < nj; ++j) {
        for (Int i = i_start; i < i_end; ++i) {
            Int d;
            if (!cur_is_high) {
                d = ng - i;                          // i=0→d=ng, i=ng-1→d=1
            } else {
                d = i - (ng + ni_core) + 1;          // i=ng+ni_core→d=1
            }

            Int donor_i;
            if (!donor_is_high) {
                donor_i = donor.ng + d;              // inward from donor MIN face
            } else {
                donor_i = donor.ng + donor.ni_core - 1 - d;  // inward from donor MAX face
            }

            node_x(i, j, k) = donor.node_x(donor_i, j, k) - tx;
            node_y(i, j, k) = donor.node_y(donor_i, j, k) - ty;
            node_z(i, j, k) = donor.node_z(donor_i, j, k) - tz;
        }}}
    } else if (cur_dim == 1) {
        Int j_start, j_end;
        if (!cur_is_high) {
            j_start = 0;
            j_end   = ng;
        } else {
            j_start = ng + nj_core;
            j_end   = nj;
        }

        for (Int k = 0; k < nk; ++k) {
        for (Int j = j_start; j < j_end; ++j) {
        for (Int i = 0; i < ni; ++i) {
            Int d;
            if (!cur_is_high) {
                d = ng - j;
            } else {
                d = j - (ng + nj_core) + 1;
            }

            Int donor_j;
            if (!donor_is_high) {
                donor_j = donor.ng + d;
            } else {
                donor_j = donor.ng + donor.nj_core - 1 - d;
            }

            node_x(i, j, k) = donor.node_x(i, donor_j, k) - tx;
            node_y(i, j, k) = donor.node_y(i, donor_j, k) - ty;
            node_z(i, j, k) = donor.node_z(i, donor_j, k) - tz;
        }}}
    } else { // cur_dim == 2
        Int k_start, k_end;
        if (!cur_is_high) {
            k_start = 0;
            k_end   = ng;
        } else {
            k_start = ng + nk_core;
            k_end   = nk;
        }

        for (Int k = k_start; k < k_end; ++k) {
        for (Int j = 0; j < nj; ++j) {
        for (Int i = 0; i < ni; ++i) {
            Int d;
            if (!cur_is_high) {
                d = ng - k;
            } else {
                d = k - (ng + nk_core) + 1;
            }

            Int donor_k;
            if (!donor_is_high) {
                donor_k = donor.ng + d;
            } else {
                donor_k = donor.ng + donor.nk_core - 1 - d;
            }

            node_x(i, j, k) = donor.node_x(i, j, donor_k) - tx;
            node_y(i, j, k) = donor.node_y(i, j, donor_k) - ty;
            node_z(i, j, k) = donor.node_z(i, j, donor_k) - tz;
        }}}
    }
}

inline void Grid::fill_ghost_face_periodic(int face, const Connectivity& conn) {
    // This implementation handles the simple case where transform = {1,2,3}
    // (no direction reversal).  Reversed transforms are TODO.

    // Determine donor face side (MIN or MAX).
    // For transform={1,2,3} the donor dimension matches cur_dim, so we
    // only need the side (high/low).  Direction-reversed transforms TODO.
    bool donor_is_high = false;
    if (conn.donor_imin == conn.donor_imax) {
        donor_is_high = (conn.donor_imin == ni_core);  // 1-based donor index
    } else if (conn.donor_jmin == conn.donor_jmax) {
        donor_is_high = (conn.donor_jmin == nj_core);
    } else if (conn.donor_kmin == conn.donor_kmax) {
        donor_is_high = (conn.donor_kmin == nk_core);
    }

    // Current face information
    // face: 0=IMIN, 1=IMAX, 2=JMIN, 3=JMAX, 4=KMIN, 5=KMAX
    int cur_dim = face / 2;          // 0=i, 1=j, 2=k
    bool cur_is_high = (face % 2 == 1); // odd = MAX side

    // The ghost sign: for the LOW side, ghost nodes are at -d direction;
    // donor interior lies inward from the donor face in the same direction.
    // Ghost coord = donor_interior_coord - translation
    // (derived from: ghost is the periodic image of donor interior,
    //  where periodic mapping gives donor_periodic = current + translation
    //  → ghost = donor_interior - translation for all faces)
    Real tx = conn.translation[0];
    Real ty = conn.translation[1];
    Real tz = conn.translation[2];

    // Iterate over the ghost region for this face
    if (cur_dim == 0) {
        // IMIN or IMAX: ghost i varies, j,k span the full extended range
        Int i_start, i_end, i_step;
        if (!cur_is_high) {
            // IMIN: ghost i = 0 .. ng-1 (outward in -i direction)
            i_start = 0;
            i_end   = ng;
            i_step  = 1;
        } else {
            // IMAX: ghost i = ng+ni_core .. ni-1 (outward in +i direction)
            i_start = ng + ni_core;
            i_end   = ni;
            i_step  = 1;
        }

        for (Int k = 0; k < nk; ++k) {
        for (Int j = 0; j < nj; ++j) {
        for (Int i = i_start; i < i_end; i += i_step) {
            // Distance from the face (1-based): ghost layer d
            Int d;
            if (!cur_is_high) {
                d = ng - i;  // i=0 → d=ng (farthest), i=ng-1 → d=1 (closest)
            } else {
                d = i - (ng + ni_core) + 1;  // i=ng+ni_core → d=1, i=ni-1 → d=ng
            }

            // Donor interior index (for the normal dimension)
            Int donor_index;
            if (donor_is_high) {
                // Donor face is a MAX side → interior is inward at face - d
                donor_index = ng + ni_core - 1 - d;  // last core node - d
            } else {
                // Donor face is a MIN side → interior is inward at face + d
                donor_index = ng + d;
            }

            // TODO: handle direction reversals in transform
            // For now assume transform = {1,2,3}: j→j, k→k
            Real xn = node_x(donor_index, j, k);
            Real yn = node_y(donor_index, j, k);
            Real zn = node_z(donor_index, j, k);

            node_x(i, j, k) = xn - tx;
            node_y(i, j, k) = yn - ty;
            node_z(i, j, k) = zn - tz;
        }}}
    } else if (cur_dim == 1) {
        // JMIN or JMAX
        Int j_start, j_end, j_step;
        if (!cur_is_high) {
            j_start = 0;
            j_end   = ng;
            j_step  = 1;
        } else {
            j_start = ng + nj_core;
            j_end   = nj;
            j_step  = 1;
        }

        for (Int k = 0; k < nk; ++k) {
        for (Int j = j_start; j < j_end; j += j_step) {
        for (Int i = 0; i < ni; ++i) {
            Int d;
            if (!cur_is_high) {
                d = ng - j;
            } else {
                d = j - (ng + nj_core) + 1;
            }

            Int donor_index;
            if (donor_is_high) {
                donor_index = ng + nj_core - 1 - d;
            } else {
                donor_index = ng + d;
            }

            Real xn = node_x(i, donor_index, k);
            Real yn = node_y(i, donor_index, k);
            Real zn = node_z(i, donor_index, k);

            node_x(i, j, k) = xn - tx;
            node_y(i, j, k) = yn - ty;
            node_z(i, j, k) = zn - tz;
        }}}
    } else {
        // KMIN or KMAX (cur_dim == 2)
        Int k_start, k_end, k_step;
        if (!cur_is_high) {
            k_start = 0;
            k_end   = ng;
            k_step  = 1;
        } else {
            k_start = ng + nk_core;
            k_end   = nk;
            k_step  = 1;
        }

        for (Int k = k_start; k < k_end; k += k_step) {
        for (Int j = 0; j < nj; ++j) {
        for (Int i = 0; i < ni; ++i) {
            Int d;
            if (!cur_is_high) {
                d = ng - k;
            } else {
                d = k - (ng + nk_core) + 1;
            }

            Int donor_index;
            if (donor_is_high) {
                donor_index = ng + nk_core - 1 - d;
            } else {
                donor_index = ng + d;
            }

            Real xn = node_x(i, j, donor_index);
            Real yn = node_y(i, j, donor_index);
            Real zn = node_z(i, j, donor_index);

            node_x(i, j, k) = xn - tx;
            node_y(i, j, k) = yn - ty;
            node_z(i, j, k) = zn - tz;
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

inline void Grid::build_face_periodic(bool face_periodic[6]) const {
    for (int f = 0; f < 6; ++f) {
        face_periodic[f] = (find_periodic_connection(f) != nullptr);
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
    // ---- 1. Determine periodic faces ----
    bool fp[6];
    build_face_periodic(fp);

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
// Face metric interpolation
// ============================================================================

inline void Grid::compute_face_metrics() {
    bool fp[6];
    build_face_periodic(fp);

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
