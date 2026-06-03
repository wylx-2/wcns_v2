#pragma once

#include "bc_applier.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// ============================================================================
// Main entry
// ============================================================================

inline void BoundaryConditionApplier::apply_all(LocalBlock& lb,
                                                  const Config& cfg) {
    // Stage 1: fill face ghost cells (BC applied to full face region,
    //          including edge/corner overlaps)
    apply_face_ghost(lb, cfg);

    // Stage 2: correct edge ghost cells (intersection of 2 faces)
    apply_edge_ghost(lb);

    // Stage 3: correct corner ghost cells (intersection of 3 faces)
    apply_corner_ghost(lb);
}

// ============================================================================
// is_bc_face
// ============================================================================

inline bool BoundaryConditionApplier::is_bc_face(const LocalBlock& lb, int face) {
    return !lb.neighbors[face].active;
}

// ============================================================================
// patch_cell_range — convert BCPatch node range to cell range for ghost filling
// ============================================================================

inline void BoundaryConditionApplier::patch_cell_range(
        int face, const BCPatch& patch, Int ng,
        Int nci_core, Int ncj_core, Int nck_core,
        Int& r0, Int& r1, Int& s0, Int& s1) {

    // BCPatch stores 1-based node indices in the sub-block's core coordinate
    // system.  After ghost extension, core node index p (1-based) maps to
    // array node index (p - 1 + ng).
    //
    // Ghost cells sit between ghost nodes.  For a given face, we need the
    // cell range in the two tangential directions.
    //
    // Example — IMIN face (face=0):
    //   The patch covers nodes with j in [jmin, jmax] and k in [kmin, kmax].
    //   The ghost cells adjacent to these nodes have cell indices:
    //     j: [jmin-1+ng, jmax-2+ng]     (the last node jmax is the right node
    //                                     of the last cell jmax-2)
    //     k: [kmin-1+ng, kmax-2+ng]
    //
    // Note: patch.imin/patch.imax are typically both 1 for IMIN face
    //       (they indicate the i-node layer of the face, not a range in i).

    switch (face) {
    case 0: // IMIN — tangential: j, k
        r0 = patch.jmin - 1 + ng;
        r1 = patch.jmax - 2 + ng;
        s0 = patch.kmin - 1 + ng;
        s1 = patch.kmax - 2 + ng;
        break;
    case 1: // IMAX — tangential: j, k
        r0 = patch.jmin - 1 + ng;
        r1 = patch.jmax - 2 + ng;
        s0 = patch.kmin - 1 + ng;
        s1 = patch.kmax - 2 + ng;
        break;
    case 2: // JMIN — tangential: i, k
        r0 = patch.imin - 1 + ng;
        r1 = patch.imax - 2 + ng;
        s0 = patch.kmin - 1 + ng;
        s1 = patch.kmax - 2 + ng;
        break;
    case 3: // JMAX — tangential: i, k
        r0 = patch.imin - 1 + ng;
        r1 = patch.imax - 2 + ng;
        s0 = patch.kmin - 1 + ng;
        s1 = patch.kmax - 2 + ng;
        break;
    case 4: // KMIN — tangential: i, j
        r0 = patch.imin - 1 + ng;
        r1 = patch.imax - 2 + ng;
        s0 = patch.jmin - 1 + ng;
        s1 = patch.jmax - 2 + ng;
        break;
    case 5: // KMAX — tangential: i, j
        r0 = patch.imin - 1 + ng;
        r1 = patch.imax - 2 + ng;
        s0 = patch.jmin - 1 + ng;
        s1 = patch.jmax - 2 + ng;
        break;
    default:
        r0 = r1 = s0 = s1 = 0;
        return;
    }

    // Clamp to valid cell range
    r0 = std::max(r0, Int(0));
    s0 = std::max(s0, Int(0));

    switch (face) {
    case 0: case 1:
        r1 = std::min(r1, ncj_core + 2 * ng - 1);
        s1 = std::min(s1, nck_core + 2 * ng - 1);
        break;
    case 2: case 3:
        r1 = std::min(r1, nci_core + 2 * ng - 1);
        s1 = std::min(s1, nck_core + 2 * ng - 1);
        break;
    case 4: case 5:
        r1 = std::min(r1, nci_core + 2 * ng - 1);
        s1 = std::min(s1, ncj_core + 2 * ng - 1);
        break;
    }
}

// ============================================================================
// mirror_index — even reflection about a face
// ============================================================================

inline Int BoundaryConditionApplier::mirror_index(int face, Int ghost_idx, Int ng) {
    // For ghost cells immediately adjacent to the interior:
    //   ghost(ng-1) reflects interior(ng) for IMIN
    //   i.e., ghost_idx and its mirror add to 2*ng - 1 for IMIN
    // For IMAX:
    //   ghost(nci - ng) reflects interior(nci - ng - 1)
    //   i.e., ghost_idx + mirror = 2*(nci-ng) - 1
    //
    // But we don't have nci here, so we return the offset.
    // The caller uses this to compute the mirror index.
    //
    // For a ghost cell at index g (0-based from face), the corresponding
    // interior cell is at:
    //   IMIN: interior = ng + (ng - 1 - g)      (g ∈ [0, ng-1])
    //   IMAX: interior = (nci - 1 - ng) - (g - (nci - ng))
    //                  = 2*(nci - ng) - 1 - g
    //
    // Since we don't have total dimensions here, we return a "mirror offset"
    // that the caller uses as: mirror = 2*base - 1 - ghost_idx
    // where base = ng for IMIN, base = nci-ng for IMAX.
    (void)face;
    (void)ghost_idx;
    (void)ng;
    return 0;  // not used directly; see per-face logic in apply_wall_*
}

// ============================================================================
// Stage 1: face ghost
// ============================================================================

inline void BoundaryConditionApplier::apply_face_ghost(LocalBlock& lb,
                                                         const Config& cfg) {
    Int nci = lb.grid.nci;
    Int ncj = lb.grid.ncj;
    Int nck = lb.grid.nck;

    // Determine the BC type for each face from the first matching patch.
    // Apply BC to the FULL face ghost region (including edges/corners),
    // since edge/corner corrections need face ghost values everywhere.
    BCType face_bc_type[6];
    bool    face_has_bc[6] = {false};

    for (Int p = 0; p < lb.grid.bc.num_patches(); ++p) {
        const BCPatch& patch = lb.grid.bc.patch(p);
        int face = patch.face;
        if (face >= 0 && face < 6 && !face_has_bc[face] && is_bc_face(lb, face)) {
            face_bc_type[face] = patch.type;
            face_has_bc[face]   = true;
        }
    }

    // Apply BC to each face, covering the ENTIRE tangential extent
    for (int face = 0; face < 6; ++face) {
        if (!is_bc_face(lb, face)) continue;

        // Full tangential range
        Int r0 = 0, r1 = 0, s0 = 0, s1 = 0;
        switch (face) {
        case 0: case 1: r1 = ncj - 1; s1 = nck - 1; break;
        case 2: case 3: r1 = nci - 1; s1 = nck - 1; break;
        case 4: case 5: r1 = nci - 1; s1 = ncj - 1; break;
        }

        if (face_has_bc[face]) {
            apply_bc_on_face(face, face_bc_type[face], lb, cfg, r0, r1, s0, s1);
        } else {
            // No BC patch — extrapolate from interior as fallback
            apply_outflow(face, lb, r0, r1, s0, s1);
        }
    }

    // Fill self-periodic faces (same-process, target_rank == -1):
    // copy from opposite face's interior cells.
    for (int face = 0; face < 6; ++face) {
        const NeighborInfo& ni = lb.neighbors[face];
        if (!ni.active) continue;
        if (ni.target_rank != -1) continue; // remote neighbor, handled by MPI

        // Self-periodic: donor is this block's opposite face.
        // Copy from the opposite face's interior layers outward.

        // Determine ranges for the full face ghost
        Int r0 = 0, r1 = 0, s0 = 0, s1 = 0;
        switch (face) {
        case 0: case 1: r1 = ncj - 1; s1 = nck - 1; break;
        case 2: case 3: r1 = nci - 1; s1 = nck - 1; break;
        case 4: case 5: r1 = nci - 1; s1 = ncj - 1; break;
        }

        Int ng = lb.grid.ng;
        for (Int d = 0; d < ng; ++d) {
            Int gi; // ghost index
            Int si; // source (donor interior) index
            switch (face) {
            case 0: gi = d;              si = nci - 1 - ng - d; break; // IMIN ← IMAX donor
            case 1: gi = nci - 1 - d;    si = ng + d;           break; // IMAX ← IMIN donor
            case 2: gi = d;              si = ncj - 1 - ng - d; break;
            case 3: gi = ncj - 1 - d;    si = ng + d;           break;
            case 4: gi = d;              si = nck - 1 - ng - d; break;
            case 5: gi = nck - 1 - d;    si = ng + d;           break;
            default: continue;
            }

            // r0..r1 = first tangential range (j for I-face, i for J/K-faces)
            // s0..s1 = second tangential range (k for I/J-face, j for K-face)
            for (Int t2 = s0; t2 <= s1; ++t2)   // second tangential index
            for (Int t1 = r0; t1 <= r1; ++t1) {  // first tangential index
                Real r, u, v, w, p;
                switch (face) {
                case 0: case 1:  // t1=j, t2=k
                    r = lb.field.prim.rho(si, t1, t2);
                    u = lb.field.prim.u(si, t1, t2);
                    v = lb.field.prim.v(si, t1, t2);
                    w = lb.field.prim.w(si, t1, t2);
                    p = lb.field.prim.p(si, t1, t2);
                    lb.field.prim.rho(gi, t1, t2) = r;
                    lb.field.prim.u(gi, t1, t2) = u;
                    lb.field.prim.v(gi, t1, t2) = v;
                    lb.field.prim.w(gi, t1, t2) = w;
                    lb.field.prim.p(gi, t1, t2) = p;
                    break;
                case 2: case 3:  // t1=i, t2=k
                    r = lb.field.prim.rho(t1, si, t2);
                    u = lb.field.prim.u(t1, si, t2);
                    v = lb.field.prim.v(t1, si, t2);
                    w = lb.field.prim.w(t1, si, t2);
                    p = lb.field.prim.p(t1, si, t2);
                    lb.field.prim.rho(t1, gi, t2) = r;
                    lb.field.prim.u(t1, gi, t2) = u;
                    lb.field.prim.v(t1, gi, t2) = v;
                    lb.field.prim.w(t1, gi, t2) = w;
                    lb.field.prim.p(t1, gi, t2) = p;
                    break;
                case 4: case 5:  // t1=i, t2=j
                    r = lb.field.prim.rho(t1, t2, si);
                    u = lb.field.prim.u(t1, t2, si);
                    v = lb.field.prim.v(t1, t2, si);
                    w = lb.field.prim.w(t1, t2, si);
                    p = lb.field.prim.p(t1, t2, si);
                    lb.field.prim.rho(t1, t2, gi) = r;
                    lb.field.prim.u(t1, t2, gi) = u;
                    lb.field.prim.v(t1, t2, gi) = v;
                    lb.field.prim.w(t1, t2, gi) = w;
                    lb.field.prim.p(t1, t2, gi) = p;
                    break;
                }
            }
        }
    }
}

// ============================================================================
// apply_bc_on_face — dispatch to specific BC type
// ============================================================================

inline void BoundaryConditionApplier::apply_bc_on_face(
        int face, BCType bc_type, LocalBlock& lb, const Config& cfg,
        Int j0, Int j1, Int k0, Int k1) {

    switch (bc_type) {
    case BCType::Farfield:
        apply_farfield(face, lb, cfg, j0, j1, k0, k1);
        break;
    case BCType::Wall:
        if (cfg.wall_type == "slip") {
            apply_wall_slip(face, lb, j0, j1, k0, k1);
        } else {
            apply_wall_noslip(face, lb, j0, j1, k0, k1);
        }
        break;
    case BCType::Symmetry:
        apply_symmetry(face, lb, j0, j1, k0, k1);
        break;
    case BCType::Inflow:
        apply_inflow(face, lb, cfg, j0, j1, k0, k1);
        break;
    case BCType::Outflow:
        apply_outflow(face, lb, j0, j1, k0, k1);
        break;
    default:
        // Unknown or Periodic — extrapolate as fallback
        apply_outflow(face, lb, j0, j1, k0, k1);
        break;
    }
}

// ============================================================================
// Farfield — 1D Riemann invariants
// ============================================================================

inline void BoundaryConditionApplier::apply_farfield(
        int face, LocalBlock& lb, const Config& cfg,
        Int j0, Int j1, Int k0, Int k1) {

    Int ng  = lb.grid.ng;
    Int nci = lb.grid.nci;
    Int ncj = lb.grid.ncj;
    Int nck = lb.grid.nck;
    Real gamma = cfg.gamma;
    Real gm1   = gamma - 1.0;

    // Free-stream values (non-dimensional)
    Real rho_inf = 1.0;
    Real u_inf, v_inf, w_inf;
    cfg.free_stream_velocity(u_inf, v_inf, w_inf);
    Real p_inf = 1.0 / (gamma * cfg.Mach * cfg.Mach);
    Real c_inf = std::sqrt(gamma * p_inf / rho_inf);

    // Determine normal velocity index and sign
    int    vel_idx = 0;     // 0=u, 1=v, 2=w  — which component is face-normal
    Real   sign    = 1.0;   // +1 for IMIN/JMIN/KMIN (ghost on low-index side)
                            // -1 for IMAX/JMAX/KMAX (ghost on high-index side)

    switch (face) {
    case 0: vel_idx = 0; sign =  1.0; break; // IMIN: +x normal
    case 1: vel_idx = 0; sign = -1.0; break; // IMAX: -x normal
    case 2: vel_idx = 1; sign =  1.0; break; // JMIN: +y normal
    case 3: vel_idx = 1; sign = -1.0; break; // JMAX: -y normal
    case 4: vel_idx = 2; sign =  1.0; break; // KMIN: +z normal
    case 5: vel_idx = 2; sign = -1.0; break; // KMAX: -z normal
    }

    // Loop over ghost layers (from inside out)
    for (Int d = 0; d < ng; ++d) {
        Int gi, ii;  // ghost index, interior reference index
        switch (face) {
        case 0: gi = d;         ii = ng + d;            break;
        case 1: gi = nci - 1 - d; ii = nci - 1 - ng - d; break;
        case 2: gi = d;         ii = ng + d;            break;
        case 3: gi = ncj - 1 - d; ii = ncj - 1 - ng - d; break;
        case 4: gi = d;         ii = ng + d;            break;
        case 5: gi = nck - 1 - d; ii = nck - 1 - ng - d; break;
        default: continue;
        }

        for (Int k = k0; k <= k1; ++k)
        for (Int j = j0; j <= j1; ++j) {

            // Interior values
            Real rho_int = 0, u_int = 0, v_int = 0, w_int = 0, p_int = 0;
            switch (face) {
            case 0: case 1:
                rho_int = lb.field.prim.rho(ii, j, k);
                u_int   = lb.field.prim.u(ii, j, k);
                v_int   = lb.field.prim.v(ii, j, k);
                w_int   = lb.field.prim.w(ii, j, k);
                p_int   = lb.field.prim.p(ii, j, k);
                break;
            case 2: case 3:
                rho_int = lb.field.prim.rho(j, ii, k);
                u_int   = lb.field.prim.u(j, ii, k);
                v_int   = lb.field.prim.v(j, ii, k);
                w_int   = lb.field.prim.w(j, ii, k);
                p_int   = lb.field.prim.p(j, ii, k);
                break;
            case 4: case 5:
                rho_int = lb.field.prim.rho(j, k, ii);
                u_int   = lb.field.prim.u(j, k, ii);
                v_int   = lb.field.prim.v(j, k, ii);
                w_int   = lb.field.prim.w(j, k, ii);
                p_int   = lb.field.prim.p(j, k, ii);
                break;
            }

            Real c_int = std::sqrt(gamma * p_int / rho_int);

            // Normal velocity (positive = outward from domain)
            Real u_n_int, u_n_inf;
            Real v_t1_int, v_t2_int, v_t1_inf, v_t2_inf;
            if (vel_idx == 0) {
                u_n_int   = sign * u_int;
                v_t1_int  = v_int;
                v_t2_int  = w_int;
                u_n_inf   = sign * u_inf;
                v_t1_inf  = v_inf;
                v_t2_inf  = w_inf;
            } else if (vel_idx == 1) {
                u_n_int   = sign * v_int;
                v_t1_int  = u_int;
                v_t2_int  = w_int;
                u_n_inf   = sign * v_inf;
                v_t1_inf  = u_inf;
                v_t2_inf  = w_inf;
            } else {
                u_n_int   = sign * w_int;
                v_t1_int  = u_int;
                v_t2_int  = v_int;
                u_n_inf   = sign * w_inf;
                v_t1_inf  = u_inf;
                v_t2_inf  = v_inf;
            }

            // Riemann invariants
            Real R_plus  = u_n_int + 2.0 * c_int / gm1;     // outgoing
            Real R_minus = u_n_inf  - 2.0 * c_inf / gm1;    // incoming

            Real u_n_bnd = 0.5 * (R_plus + R_minus);
            Real c_bnd   = 0.25 * gm1 * (R_plus - R_minus);

            // Entropy: based on flow direction at boundary
            Real s_bnd;
            Real v_t1_bnd, v_t2_bnd;
            if (u_n_bnd > 0.0) {
                // Outflow: use interior values
                s_bnd = p_int / std::pow(rho_int, gamma);
                v_t1_bnd = v_t1_int;
                v_t2_bnd = v_t2_int;
            } else {
                // Inflow: use free-stream values
                s_bnd = p_inf / std::pow(rho_inf, gamma);
                v_t1_bnd = v_t1_inf;
                v_t2_bnd = v_t2_inf;
            }

            Real rho_bnd = std::pow(c_bnd * c_bnd / (gamma * s_bnd), 1.0 / gm1);
            Real p_bnd   = s_bnd * std::pow(rho_bnd, gamma);

            // Convert normal/tangential back to Cartesian
            Real u_bnd, v_bnd, w_bnd;
            if (vel_idx == 0) {
                u_bnd = sign * u_n_bnd;
                v_bnd = v_t1_bnd;
                w_bnd = v_t2_bnd;
            } else if (vel_idx == 1) {
                u_bnd = v_t1_bnd;
                v_bnd = sign * u_n_bnd;
                w_bnd = v_t2_bnd;
            } else {
                u_bnd = v_t1_bnd;
                v_bnd = v_t2_bnd;
                w_bnd = sign * u_n_bnd;
            }

            // Write ghost cell
            switch (face) {
            case 0: case 1:
                lb.field.prim.rho(gi, j, k) = rho_bnd;
                lb.field.prim.u(gi, j, k)   = u_bnd;
                lb.field.prim.v(gi, j, k)   = v_bnd;
                lb.field.prim.w(gi, j, k)   = w_bnd;
                lb.field.prim.p(gi, j, k)   = p_bnd;
                break;
            case 2: case 3:
                lb.field.prim.rho(j, gi, k) = rho_bnd;
                lb.field.prim.u(j, gi, k)   = u_bnd;
                lb.field.prim.v(j, gi, k)   = v_bnd;
                lb.field.prim.w(j, gi, k)   = w_bnd;
                lb.field.prim.p(j, gi, k)   = p_bnd;
                break;
            case 4: case 5:
                lb.field.prim.rho(j, k, gi) = rho_bnd;
                lb.field.prim.u(j, k, gi)   = u_bnd;
                lb.field.prim.v(j, k, gi)   = v_bnd;
                lb.field.prim.w(j, k, gi)   = w_bnd;
                lb.field.prim.p(j, k, gi)   = p_bnd;
                break;
            }
        }
    }
}

// ============================================================================
// Wall — no-slip, adiabatic (even reflection)
// ============================================================================

inline void BoundaryConditionApplier::apply_wall_noslip(
        int face, LocalBlock& lb,
        Int j0, Int j1, Int k0, Int k1) {

    Int ng  = lb.grid.ng;
    Int nci = lb.grid.nci;
    Int ncj = lb.grid.ncj;
    Int nck = lb.grid.nck;

    for (Int d = 0; d < ng; ++d) {
        Int gi, mi;  // ghost index, mirror interior index
        switch (face) {
        case 0: gi = d;              mi = ng + (ng - 1 - d);    break;
        case 1: gi = nci - 1 - d;    mi = nci - 1 - ng - (ng - 1 - d); break;
        case 2: gi = d;              mi = ng + (ng - 1 - d);    break;
        case 3: gi = ncj - 1 - d;    mi = ncj - 1 - ng - (ng - 1 - d); break;
        case 4: gi = d;              mi = ng + (ng - 1 - d);    break;
        case 5: gi = nck - 1 - d;    mi = nck - 1 - ng - (ng - 1 - d); break;
        default: continue;
        }

        for (Int k = k0; k <= k1; ++k)
        for (Int j = j0; j <= j1; ++j) {

            Real rho_mir = 0, u_mir = 0, v_mir = 0, w_mir = 0, p_mir = 0;
            switch (face) {
            case 0: case 1:
                rho_mir = lb.field.prim.rho(mi, j, k);
                u_mir   = lb.field.prim.u(mi, j, k);
                v_mir   = lb.field.prim.v(mi, j, k);
                w_mir   = lb.field.prim.w(mi, j, k);
                p_mir   = lb.field.prim.p(mi, j, k);
                break;
            case 2: case 3:
                rho_mir = lb.field.prim.rho(j, mi, k);
                u_mir   = lb.field.prim.u(j, mi, k);
                v_mir   = lb.field.prim.v(j, mi, k);
                w_mir   = lb.field.prim.w(j, mi, k);
                p_mir   = lb.field.prim.p(j, mi, k);
                break;
            case 4: case 5:
            default:
                rho_mir = lb.field.prim.rho(j, k, mi);
                u_mir   = lb.field.prim.u(j, k, mi);
                v_mir   = lb.field.prim.v(j, k, mi);
                w_mir   = lb.field.prim.w(j, k, mi);
                p_mir   = lb.field.prim.p(j, k, mi);
                break;
            }

            // No-slip: all velocity components reflected (negated)
            // Adiabatic: density and pressure zero-gradient (copy)
            switch (face) {
            case 0: case 1:
                lb.field.prim.rho(gi, j, k) = rho_mir;
                lb.field.prim.u(gi, j, k)   = -u_mir;
                lb.field.prim.v(gi, j, k)   = -v_mir;
                lb.field.prim.w(gi, j, k)   = -w_mir;
                lb.field.prim.p(gi, j, k)   = p_mir;
                break;
            case 2: case 3:
                lb.field.prim.rho(j, gi, k) = rho_mir;
                lb.field.prim.u(j, gi, k)   = -u_mir;
                lb.field.prim.v(j, gi, k)   = -v_mir;
                lb.field.prim.w(j, gi, k)   = -w_mir;
                lb.field.prim.p(j, gi, k)   = p_mir;
                break;
            case 4: case 5:
                lb.field.prim.rho(j, k, gi) = rho_mir;
                lb.field.prim.u(j, k, gi)   = -u_mir;
                lb.field.prim.v(j, k, gi)   = -v_mir;
                lb.field.prim.w(j, k, gi)   = -w_mir;
                lb.field.prim.p(j, k, gi)   = p_mir;
                break;
            }
        }
    }
}

// ============================================================================
// Wall — slip / inviscid (normal velocity reflected only)
// ============================================================================

inline void BoundaryConditionApplier::apply_wall_slip(
        int face, LocalBlock& lb,
        Int j0, Int j1, Int k0, Int k1) {
    // For slip wall, treat same as symmetry: normal velocity reflected,
    // tangential and scalars zero-gradient.
    apply_symmetry(face, lb, j0, j1, k0, k1);
}

// ============================================================================
// Symmetry — normal velocity reflected, others zero-gradient
// ============================================================================

inline void BoundaryConditionApplier::apply_symmetry(
        int face, LocalBlock& lb,
        Int j0, Int j1, Int k0, Int k1) {

    Int ng  = lb.grid.ng;
    Int nci = lb.grid.nci;
    Int ncj = lb.grid.ncj;
    Int nck = lb.grid.nck;

    // Determine which velocity component is face-normal
    int vel_idx = 0; // 0=u, 1=v, 2=w
    switch (face) {
    case 0: case 1: vel_idx = 0; break;
    case 2: case 3: vel_idx = 1; break;
    case 4: case 5: vel_idx = 2; break;
    }

    for (Int d = 0; d < ng; ++d) {
        Int gi, mi;
        switch (face) {
        case 0: gi = d;              mi = ng + (ng - 1 - d);    break;
        case 1: gi = nci - 1 - d;    mi = nci - 1 - ng - (ng - 1 - d); break;
        case 2: gi = d;              mi = ng + (ng - 1 - d);    break;
        case 3: gi = ncj - 1 - d;    mi = ncj - 1 - ng - (ng - 1 - d); break;
        case 4: gi = d;              mi = ng + (ng - 1 - d);    break;
        case 5: gi = nck - 1 - d;    mi = nck - 1 - ng - (ng - 1 - d); break;
        default: continue;
        }

        for (Int k = k0; k <= k1; ++k)
        for (Int j = j0; j <= j1; ++j) {

            Real rho_mir = 0, u_mir = 0, v_mir = 0, w_mir = 0, p_mir = 0;
            switch (face) {
            case 0: case 1:
                rho_mir = lb.field.prim.rho(mi, j, k);
                u_mir   = lb.field.prim.u(mi, j, k);
                v_mir   = lb.field.prim.v(mi, j, k);
                w_mir   = lb.field.prim.w(mi, j, k);
                p_mir   = lb.field.prim.p(mi, j, k);
                break;
            case 2: case 3:
                rho_mir = lb.field.prim.rho(j, mi, k);
                u_mir   = lb.field.prim.u(j, mi, k);
                v_mir   = lb.field.prim.v(j, mi, k);
                w_mir   = lb.field.prim.w(j, mi, k);
                p_mir   = lb.field.prim.p(j, mi, k);
                break;
            case 4: case 5:
            default:
                rho_mir = lb.field.prim.rho(j, k, mi);
                u_mir   = lb.field.prim.u(j, k, mi);
                v_mir   = lb.field.prim.v(j, k, mi);
                w_mir   = lb.field.prim.w(j, k, mi);
                p_mir   = lb.field.prim.p(j, k, mi);
                break;
            }

            Real u_g = u_mir, v_g = v_mir, w_g = w_mir;
            if (vel_idx == 0)      u_g = -u_mir;
            else if (vel_idx == 1) v_g = -v_mir;
            else                   w_g = -w_mir;

            switch (face) {
            case 0: case 1:
                lb.field.prim.rho(gi, j, k) = rho_mir;
                lb.field.prim.u(gi, j, k)   = u_g;
                lb.field.prim.v(gi, j, k)   = v_g;
                lb.field.prim.w(gi, j, k)   = w_g;
                lb.field.prim.p(gi, j, k)   = p_mir;
                break;
            case 2: case 3:
                lb.field.prim.rho(j, gi, k) = rho_mir;
                lb.field.prim.u(j, gi, k)   = u_g;
                lb.field.prim.v(j, gi, k)   = v_g;
                lb.field.prim.w(j, gi, k)   = w_g;
                lb.field.prim.p(j, gi, k)   = p_mir;
                break;
            case 4: case 5:
                lb.field.prim.rho(j, k, gi) = rho_mir;
                lb.field.prim.u(j, k, gi)   = u_g;
                lb.field.prim.v(j, k, gi)   = v_g;
                lb.field.prim.w(j, k, gi)   = w_g;
                lb.field.prim.p(j, k, gi)   = p_mir;
                break;
            }
        }
    }
}

// ============================================================================
// Inflow — all variables specified from free-stream
// ============================================================================

inline void BoundaryConditionApplier::apply_inflow(
        int face, LocalBlock& lb, const Config& cfg,
        Int j0, Int j1, Int k0, Int k1) {

    Int ng  = lb.grid.ng;
    Int nci = lb.grid.nci;
    Int ncj = lb.grid.ncj;
    Int nck = lb.grid.nck;

    Real rho_inf = 1.0;
    Real u_inf, v_inf, w_inf;
    cfg.free_stream_velocity(u_inf, v_inf, w_inf);
    Real p_inf = 1.0 / (cfg.gamma * cfg.Mach * cfg.Mach);

    for (Int d = 0; d < ng; ++d) {
        Int gi;
        switch (face) {
        case 0: gi = d;              break;
        case 1: gi = nci - 1 - d;    break;
        case 2: gi = d;              break;
        case 3: gi = ncj - 1 - d;    break;
        case 4: gi = d;              break;
        case 5: gi = nck - 1 - d;    break;
        default: continue;
        }

        for (Int k = k0; k <= k1; ++k)
        for (Int j = j0; j <= j1; ++j) {
            switch (face) {
            case 0: case 1:
                lb.field.prim.rho(gi, j, k) = rho_inf;
                lb.field.prim.u(gi, j, k)   = u_inf;
                lb.field.prim.v(gi, j, k)   = v_inf;
                lb.field.prim.w(gi, j, k)   = w_inf;
                lb.field.prim.p(gi, j, k)   = p_inf;
                break;
            case 2: case 3:
                lb.field.prim.rho(j, gi, k) = rho_inf;
                lb.field.prim.u(j, gi, k)   = u_inf;
                lb.field.prim.v(j, gi, k)   = v_inf;
                lb.field.prim.w(j, gi, k)   = w_inf;
                lb.field.prim.p(j, gi, k)   = p_inf;
                break;
            case 4: case 5:
                lb.field.prim.rho(j, k, gi) = rho_inf;
                lb.field.prim.u(j, k, gi)   = u_inf;
                lb.field.prim.v(j, k, gi)   = v_inf;
                lb.field.prim.w(j, k, gi)   = w_inf;
                lb.field.prim.p(j, k, gi)   = p_inf;
                break;
            }
        }
    }
}

// ============================================================================
// Outflow — zero-order extrapolation from interior
// ============================================================================

inline void BoundaryConditionApplier::apply_outflow(
        int face, LocalBlock& lb,
        Int j0, Int j1, Int k0, Int k1) {

    Int ng  = lb.grid.ng;
    Int nci = lb.grid.nci;
    Int ncj = lb.grid.ncj;
    Int nck = lb.grid.nck;

    // Reference interior index: the first interior cell adjacent to the face
    Int int_i = 0;
    switch (face) {
    case 0: int_i = ng;               break;
    case 1: int_i = nci - 1 - ng;     break;
    case 2: int_i = ng;               break;
    case 3: int_i = ncj - 1 - ng;     break;
    case 4: int_i = ng;               break;
    case 5: int_i = nck - 1 - ng;     break;
    }

    for (Int d = 0; d < ng; ++d) {
        Int gi;
        switch (face) {
        case 0: gi = d;              break;
        case 1: gi = nci - 1 - d;    break;
        case 2: gi = d;              break;
        case 3: gi = ncj - 1 - d;    break;
        case 4: gi = d;              break;
        case 5: gi = nck - 1 - d;    break;
        default: continue;
        }

        for (Int k = k0; k <= k1; ++k)
        for (Int j = j0; j <= j1; ++j) {
            Real rho_ref, u_ref, v_ref, w_ref, p_ref;
            switch (face) {
            case 0: case 1:
                rho_ref = lb.field.prim.rho(int_i, j, k);
                u_ref   = lb.field.prim.u(int_i, j, k);
                v_ref   = lb.field.prim.v(int_i, j, k);
                w_ref   = lb.field.prim.w(int_i, j, k);
                p_ref   = lb.field.prim.p(int_i, j, k);
                break;
            case 2: case 3:
                rho_ref = lb.field.prim.rho(j, int_i, k);
                u_ref   = lb.field.prim.u(j, int_i, k);
                v_ref   = lb.field.prim.v(j, int_i, k);
                w_ref   = lb.field.prim.w(j, int_i, k);
                p_ref   = lb.field.prim.p(j, int_i, k);
                break;
            case 4: case 5:
                rho_ref = lb.field.prim.rho(j, k, int_i);
                u_ref   = lb.field.prim.u(j, k, int_i);
                v_ref   = lb.field.prim.v(j, k, int_i);
                w_ref   = lb.field.prim.w(j, k, int_i);
                p_ref   = lb.field.prim.p(j, k, int_i);
                break;
            }

            switch (face) {
            case 0: case 1:
                lb.field.prim.rho(gi, j, k) = rho_ref;
                lb.field.prim.u(gi, j, k)   = u_ref;
                lb.field.prim.v(gi, j, k)   = v_ref;
                lb.field.prim.w(gi, j, k)   = w_ref;
                lb.field.prim.p(gi, j, k)   = p_ref;
                break;
            case 2: case 3:
                lb.field.prim.rho(j, gi, k) = rho_ref;
                lb.field.prim.u(j, gi, k)   = u_ref;
                lb.field.prim.v(j, gi, k)   = v_ref;
                lb.field.prim.w(j, gi, k)   = w_ref;
                lb.field.prim.p(j, gi, k)   = p_ref;
                break;
            case 4: case 5:
                lb.field.prim.rho(j, k, gi) = rho_ref;
                lb.field.prim.u(j, k, gi)   = u_ref;
                lb.field.prim.v(j, k, gi)   = v_ref;
                lb.field.prim.w(j, k, gi)   = w_ref;
                lb.field.prim.p(j, k, gi)   = p_ref;
                break;
            }
        }
    }
}

// ============================================================================
// Stage 2: edge ghost (intersection of two faces)
// ============================================================================

inline void BoundaryConditionApplier::apply_edge_ghost(LocalBlock& lb) {
    Int ng  = lb.grid.ng;
    Int nci = lb.grid.nci;
    Int ncj = lb.grid.ncj;
    Int nck = lb.grid.nck;

    // Helper lambda: average primitives from two cells
    auto avg_prim = [](Real r1, Real u1, Real v1, Real w1, Real p1,
                       Real r2, Real u2, Real v2, Real w2, Real p2,
                       Real& ra, Real& ua, Real& va, Real& wa, Real& pa) {
        ra = 0.5 * (r1 + r2);
        ua = 0.5 * (u1 + u2);
        va = 0.5 * (v1 + v2);
        wa = 0.5 * (w1 + w2);
        pa = 0.5 * (p1 + p2);
    };

    // Lambda: get value from face ghost cell at offset from edge
    // dir == 0: value from face A (first face of the edge)
    // dir == 1: value from face B (second face of the edge)
    auto get_face_val = [&](int face, Int gi, Int gj, Int gk,
                            Real& r, Real& u, Real& v, Real& w, Real& p) {
        switch (face) {
        case 0: case 1:
            r = lb.field.prim.rho(gi, gj, gk);
            u = lb.field.prim.u(gi, gj, gk);
            v = lb.field.prim.v(gi, gj, gk);
            w = lb.field.prim.w(gi, gj, gk);
            p = lb.field.prim.p(gi, gj, gk);
            break;
        case 2: case 3:
            r = lb.field.prim.rho(gj, gi, gk);
            u = lb.field.prim.u(gj, gi, gk);
            v = lb.field.prim.v(gj, gi, gk);
            w = lb.field.prim.w(gj, gi, gk);
            p = lb.field.prim.p(gj, gi, gk);
            break;
        case 4: case 5:
        default:
            r = lb.field.prim.rho(gj, gk, gi);
            u = lb.field.prim.u(gj, gk, gi);
            v = lb.field.prim.v(gj, gk, gi);
            w = lb.field.prim.w(gj, gk, gi);
            p = lb.field.prim.p(gj, gk, gi);
            break;
        }
    };

    auto set_edge_val = [&](int face_a, int /*face_b*/,
                            Int gi, Int gj, Int gk,
                            Real r, Real u, Real v, Real w, Real p) {
        // Write using the coordinate convention of face_a
        // (both faces index the same physical cell, just with different
        //  index ordering — the cell is unambiguous from (gi,gj,gk) context)
        // We write via face_a's convention.
        switch (face_a) {
        case 0: case 1:
            lb.field.prim.rho(gi, gj, gk) = r;
            lb.field.prim.u(gi, gj, gk) = u;
            lb.field.prim.v(gi, gj, gk) = v;
            lb.field.prim.w(gi, gj, gk) = w;
            lb.field.prim.p(gi, gj, gk) = p;
            break;
        case 2: case 3:
            lb.field.prim.rho(gj, gi, gk) = r;
            lb.field.prim.u(gj, gi, gk) = u;
            lb.field.prim.v(gj, gi, gk) = v;
            lb.field.prim.w(gj, gi, gk) = w;
            lb.field.prim.p(gj, gi, gk) = p;
            break;
        case 4: case 5:
            lb.field.prim.rho(gj, gk, gi) = r;
            lb.field.prim.u(gj, gk, gi) = u;
            lb.field.prim.v(gj, gk, gi) = v;
            lb.field.prim.w(gj, gk, gi) = w;
            lb.field.prim.p(gj, gk, gi) = p;
            break;
        }
    };

    // Define all 12 edges: {{face_a, face_b}, ...}
    // Edge cells: indices i in [0,ng-1] or [end-ng, end-1] for each direction
    const int edges[12][2] = {
        {0,2}, {0,3}, {0,4}, {0,5},   // IMIN + each J/K face
        {1,2}, {1,3}, {1,4}, {1,5},   // IMAX + each J/K face
        {2,4}, {2,5}, {3,4}, {3,5}    // JMIN/JMAX + each K face
    };

    for (auto& e : edges) {
        int fa = e[0], fb = e[1];
        bool bc_a = is_bc_face(lb, fa);
        bool bc_b = is_bc_face(lb, fb);

        if (!bc_a && !bc_b) continue; // both are parallel neighbors

        // Determine index ranges for this edge
        Int ia0 = 0, ia1 = -1, ja0 = 0, ja1 = -1, ka0 = 0, ka1 = -1;
        Int ib0 = 0, ib1 = -1, jb0 = 0, jb1 = -1, kb0 = 0, kb1 = -1;

        // face_a ghost range
        switch (fa) {
        case 0: ia0 = 0;        ia1 = ng-1;      ja0 = 0; ja1 = ncj-1; ka0 = 0; ka1 = nck-1; break;
        case 1: ia0 = nci-ng;   ia1 = nci-1;     ja0 = 0; ja1 = ncj-1; ka0 = 0; ka1 = nck-1; break;
        case 2: ja0 = 0;        ja1 = ng-1;      ia0 = 0; ia1 = nci-1; ka0 = 0; ka1 = nck-1; break;
        case 3: ja0 = ncj-ng;   ja1 = ncj-1;     ia0 = 0; ia1 = nci-1; ka0 = 0; ka1 = nck-1; break;
        case 4: ka0 = 0;        ka1 = ng-1;      ia0 = 0; ia1 = nci-1; ja0 = 0; ja1 = ncj-1; break;
        case 5: ka0 = nck-ng;   ka1 = nck-1;     ia0 = 0; ia1 = nci-1; ja0 = 0; ja1 = ncj-1; break;
        }

        // face_b ghost range
        switch (fb) {
        case 0: ib0 = 0;        ib1 = ng-1;      jb0 = 0; jb1 = ncj-1; kb0 = 0; kb1 = nck-1; break;
        case 1: ib0 = nci-ng;   ib1 = nci-1;     jb0 = 0; jb1 = ncj-1; kb0 = 0; kb1 = nck-1; break;
        case 2: jb0 = 0;        jb1 = ng-1;      ib0 = 0; ib1 = nci-1; kb0 = 0; kb1 = nck-1; break;
        case 3: jb0 = ncj-ng;   jb1 = ncj-1;     ib0 = 0; ib1 = nci-1; kb0 = 0; kb1 = nck-1; break;
        case 4: kb0 = 0;        kb1 = ng-1;      ib0 = 0; ib1 = nci-1; jb0 = 0; jb1 = ncj-1; break;
        case 5: kb0 = nck-ng;   kb1 = nck-1;     ib0 = 0; ib1 = nci-1; jb0 = 0; jb1 = ncj-1; break;
        }

        // Intersection of the two ghost ranges = edge region
        Int i0 = std::max(ia0, ib0), i1 = std::min(ia1, ib1);
        Int j0 = std::max(ja0, jb0), j1 = std::min(ja1, jb1);
        Int k0 = std::max(ka0, kb0), k1 = std::min(ka1, kb1);

        if (i0 > i1 || j0 > j1 || k0 > k1) continue;

        for (Int k = k0; k <= k1; ++k)
        for (Int j = j0; j <= j1; ++j)
        for (Int i = i0; i <= i1; ++i) {

            Real r_a = 0, u_a = 0, v_a = 0, w_a = 0, p_a = 0;
            Real r_b = 0, u_b = 0, v_b = 0, w_b = 0, p_b = 0;

            get_face_val(fa, i, j, k, r_a, u_a, v_a, w_a, p_a);
            get_face_val(fb, i, j, k, r_b, u_b, v_b, w_b, p_b);

            Real ra, ua, va, wa, pa;
            if (bc_a && bc_b) {
                avg_prim(r_a, u_a, v_a, w_a, p_a,
                         r_b, u_b, v_b, w_b, p_b,
                         ra, ua, va, wa, pa);
            } else if (bc_a) {
                ra = r_a; ua = u_a; va = v_a; wa = w_a; pa = p_a;
            } else {
                ra = r_b; ua = u_b; va = v_b; wa = w_b; pa = p_b;
            }

            set_edge_val(fa, fb, i, j, k, ra, ua, va, wa, pa);
        }
    }
}

// ============================================================================
// Stage 3: corner ghost (intersection of three faces)
// ============================================================================

inline void BoundaryConditionApplier::apply_corner_ghost(LocalBlock& lb) {
    Int ng  = lb.grid.ng;
    Int nci = lb.grid.nci;
    Int ncj = lb.grid.ncj;
    Int nck = lb.grid.nck;

    // 8 corners: I-face × J-face × K-face
    const int corners[8][3] = {
        {0,2,4}, {0,2,5}, {0,3,4}, {0,3,5},
        {1,2,4}, {1,2,5}, {1,3,4}, {1,3,5}
    };

    for (auto& c : corners) {
        int f0 = c[0], f1 = c[1], f2 = c[2];

        int bc_count = 0;
        if (is_bc_face(lb, f0)) ++bc_count;
        if (is_bc_face(lb, f1)) ++bc_count;
        if (is_bc_face(lb, f2)) ++bc_count;
        if (bc_count == 0) continue; // all parallel

        // Determine i,j,k range for this corner
        Int i0 = 0, i1 = -1, j0 = 0, j1 = -1, k0 = 0, k1 = -1;

        // I-face determines i range
        switch (f0) {
        case 0: i0 = 0; i1 = ng-1; break;
        case 1: i0 = nci-ng; i1 = nci-1; break;
        }
        // J-face determines j range
        switch (f1) {
        case 2: j0 = 0; j1 = ng-1; break;
        case 3: j0 = ncj-ng; j1 = ncj-1; break;
        }
        // K-face determines k range
        switch (f2) {
        case 4: k0 = 0; k1 = ng-1; break;
        case 5: k0 = nck-ng; k1 = nck-1; break;
        }

        for (Int k = k0; k <= k1; ++k)
        for (Int j = j0; j <= j1; ++j)
        for (Int i = i0; i <= i1; ++i) {

            // Collect values from BC faces only
            Real r_sum = 0, u_sum = 0, v_sum = 0, w_sum = 0, p_sum = 0;
            int count = 0;

            auto accum = [&](int face) {
                if (!is_bc_face(lb, face)) return;
                Real r, u, v, w, pv;
                // Use the same indexing as the main grid (i,j,k)
                r  = lb.field.prim.rho(i, j, k);
                u  = lb.field.prim.u(i, j, k);
                v  = lb.field.prim.v(i, j, k);
                w  = lb.field.prim.w(i, j, k);
                pv = lb.field.prim.p(i, j, k);
                r_sum += r; u_sum += u; v_sum += v; w_sum += w; p_sum += pv;
                ++count;
            };

            accum(f0);
            accum(f1);
            accum(f2);

            if (count > 0) {
                Real inv = 1.0 / static_cast<Real>(count);
                lb.field.prim.rho(i, j, k) = r_sum * inv;
                lb.field.prim.u(i, j, k)   = u_sum * inv;
                lb.field.prim.v(i, j, k)   = v_sum * inv;
                lb.field.prim.w(i, j, k)   = w_sum * inv;
                lb.field.prim.p(i, j, k)   = p_sum * inv;
            }
        }
    }
}
