#pragma once

#include "local_block.h"
#include <algorithm>
#include <iostream>

// ============================================================================
// NeighborInfo
// ============================================================================

inline NeighborInfo::NeighborInfo()
    : target_rank(-1), target_block(-1), target_face(-1), active(false), is_periodic(false) {
    transform[0] = 1; transform[1] = 2; transform[2] = 3;
}

// ============================================================================
// LocalBlock
// ============================================================================

inline LocalBlock::LocalBlock() : block_id(-1), zone_id(-1), sub_id(-1) {}

inline bool LocalBlock::neighbor_is_remote(int face) const {
    return neighbors[face].active && neighbors[face].target_rank >= 0
        && neighbors[face].target_rank != -1;  // -1 means same-process internal
}

inline int LocalBlock::neighbor_rank(int face) const {
    return neighbors[face].target_rank;
}

inline Int LocalBlock::nci_core() const { return grid.ni_core - 1; }
inline Int LocalBlock::ncj_core() const { return grid.nj_core - 1; }
inline Int LocalBlock::nck_core() const { return grid.nk_core - 1; }

// ============================================================================
// Construction — full zone (no splitting)
// ============================================================================

inline LocalBlock LocalBlock::from_full_zone(
        const Grid& full_zone, int block_id, int zone_id, int my_rank,
        const std::vector<SubBlock>& all_decomp) {

    LocalBlock lb;
    lb.block_id = block_id;
    lb.zone_id  = zone_id;
    lb.sub_id   = 0;

    // Deep-copy the grid (node arrays, cell arrays, metrics)
    lb.grid = full_zone;

    // Allocate field matching grid dimensions
    lb.field.allocate(lb.grid.nci, lb.grid.ncj, lb.grid.nck);

    // Build neighbor info
    SubBlock whole_zone;
    whole_zone.zone_id       = zone_id;
    whole_zone.sub_id        = 0;
    whole_zone.assigned_rank = my_rank;
    whole_zone.ci_min = 0;  whole_zone.ci_max = full_zone.ni_core - 2;
    whole_zone.cj_min = 0;  whole_zone.cj_max = full_zone.nj_core - 2;
    whole_zone.ck_min = 0;  whole_zone.ck_max = full_zone.nk_core - 2;

    build_neighbors(lb, full_zone, whole_zone, my_rank, all_decomp);

    return lb;
}

// ============================================================================
// Construction — sub-zone (after Cartesian splitting)
// ============================================================================

inline LocalBlock LocalBlock::from_sub_zone(
        const Grid& full_zone, const SubBlock& sub, int block_id,
        int my_rank, const std::vector<SubBlock>& all_decomp) {

    LocalBlock lb;
    lb.block_id = block_id;
    lb.zone_id  = sub.zone_id;
    lb.sub_id   = sub.sub_id;

    Int ng = full_zone.ng;

    // ---- Compute sub-block dimensions ----
    Int nci_sub_core = sub.ci_max - sub.ci_min + 1;
    Int ncj_sub_core = sub.cj_max - sub.cj_min + 1;
    Int nck_sub_core = sub.ck_max - sub.ck_min + 1;

    Int ni_sub_core = nci_sub_core + 1;
    Int nj_sub_core = ncj_sub_core + 1;
    Int nk_sub_core = nck_sub_core + 1;

    Int ni_tot = ni_sub_core + 2 * ng;
    Int nj_tot = nj_sub_core + 2 * ng;
    Int nk_tot = nk_sub_core + 2 * ng;

    // ---- Allocate sub-grid node arrays ----
    lb.grid.ni = ni_tot;
    lb.grid.nj = nj_tot;
    lb.grid.nk = nk_tot;
    lb.grid.nci = ni_tot - 1;
    lb.grid.ncj = nj_tot - 1;
    lb.grid.nck = nk_tot - 1;
    lb.grid.ng = ng;
    lb.grid.ni_core = ni_sub_core;
    lb.grid.nj_core = nj_sub_core;
    lb.grid.nk_core = nk_sub_core;
    lb.grid.name = full_zone.name + "_sub" + std::to_string(sub.sub_id);

    lb.grid.node_x.allocate(ni_tot, nj_tot, nk_tot);
    lb.grid.node_y.allocate(ni_tot, nj_tot, nk_tot);
    lb.grid.node_z.allocate(ni_tot, nj_tot, nk_tot);

    // Allocate cell-center arrays (needed by compute_cell_centers and metric routines)
    lb.grid.cell_x.allocate(lb.grid.nci, lb.grid.ncj, lb.grid.nck);
    lb.grid.cell_y.allocate(lb.grid.nci, lb.grid.ncj, lb.grid.nck);
    lb.grid.cell_z.allocate(lb.grid.nci, lb.grid.ncj, lb.grid.nck);
    lb.grid.cell_vol.allocate(lb.grid.nci, lb.grid.ncj, lb.grid.nck);

    // ---- Copy node coordinates from full zone ----
    // Map: sub-block node (i,j,k) -> full zone node (i + offset_i, j + offset_j, k + offset_k)
    // sub-block node 0 corresponds to ghost region left of cell ci_min
    // In full zone, cell ci_min's left node is at node index ci_min + ng
    // So sub-block node 0 = full zone node (ci_min + ng) - ng = ci_min
    Int offset_i = sub.ci_min;        // full zone node index for sub-block ghost start
    Int offset_j = sub.cj_min;
    Int offset_k = sub.ck_min;

    for (Int k = 0; k < nk_tot; ++k)
    for (Int j = 0; j < nj_tot; ++j)
    for (Int i = 0; i < ni_tot; ++i) {
        lb.grid.node_x(i,j,k) = full_zone.node_x(i + offset_i, j + offset_j, k + offset_k);
        lb.grid.node_y(i,j,k) = full_zone.node_y(i + offset_i, j + offset_j, k + offset_k);
        lb.grid.node_z(i,j,k) = full_zone.node_z(i + offset_i, j + offset_j, k + offset_k);
    }

    // ---- Compute cell centers and volumes ----
    lb.grid.compute_cell_centers();
    lb.grid.compute_cell_volumes();

    // ---- Compute metrics ----
    lb.grid.compute_metrics();
    lb.grid.compute_face_metrics();

    // ---- Copy boundary conditions (only patches that overlap this sub-block) ----
    for (Int b = 0; b < full_zone.bc.num_patches(); ++b) {
        const BCPatch& src = full_zone.bc.patch(b);

        // Determine if this BC overlaps with the sub-block's face
        // src.range is in node indices (1-based) of original core
        // Map to 0-based cell range for comparison
        Int src_ci_min = src.imin - 1;
        Int src_ci_max = src.imax - 1;
        Int src_cj_min = src.jmin - 1;
        Int src_cj_max = src.jmax - 1;
        Int src_ck_min = src.kmin - 1;
        Int src_ck_max = src.kmax - 1;

        // Check overlap with sub-block's cell range (also 0-based in original space)
        bool overlap =
            !(src_ci_max < sub.ci_min || src_ci_min > sub.ci_max ||
              src_cj_max < sub.cj_min || src_cj_min > sub.cj_max ||
              src_ck_max < sub.ck_min || src_ck_min > sub.ck_max);

        if (overlap) {
            BCPatch patch = src;
            // Remap ranges to sub-block local coordinates (1-based)
            patch.imin = std::max(src.imin, static_cast<Int>(sub.ci_min + 1))
                       - static_cast<Int>(sub.ci_min);
            patch.imax = std::min(src.imax, static_cast<Int>(sub.ci_max + 2))
                       - static_cast<Int>(sub.ci_min);
            patch.jmin = std::max(src.jmin, static_cast<Int>(sub.cj_min + 1))
                       - static_cast<Int>(sub.cj_min);
            patch.jmax = std::min(src.jmax, static_cast<Int>(sub.cj_max + 2))
                       - static_cast<Int>(sub.cj_min);
            patch.kmin = std::max(src.kmin, static_cast<Int>(sub.ck_min + 1))
                       - static_cast<Int>(sub.ck_min);
            patch.kmax = std::min(src.kmax, static_cast<Int>(sub.ck_max + 2))
                       - static_cast<Int>(sub.ck_min);
            lb.grid.bc.add_patch(patch);
        }
    }

    // ---- Periodic connectivity is handled in build_neighbors via the full zone ----
    // No need to copy here — neighbor info is built separately.

    // ---- Allocate field ----
    lb.field.allocate(lb.grid.nci, lb.grid.ncj, lb.grid.nck);

    // ---- Build neighbor info ----
    build_neighbors(lb, full_zone, sub, my_rank, all_decomp);

    return lb;
}

// ============================================================================
// Neighbor construction
// ============================================================================

inline void LocalBlock::build_neighbors(
        LocalBlock& lb, const Grid& full_zone, const SubBlock& sub,
        int /*my_rank*/, const std::vector<SubBlock>& all_decomp) {

    for (int face = 0; face < 6; ++face) {
        NeighborInfo& ni = lb.neighbors[face];
        ni = NeighborInfo();  // reset to inactive

        // Check if this face is external (touches original zone boundary)
        if (is_external_face(face, sub, full_zone)) {
            // Check for periodic connection
            const Connectivity* pc = find_periodic_at_face(full_zone, face, sub);

            if (pc && pc->is_periodic) {
                ni.active      = true;
                ni.is_periodic = true;
                ni.transform[0] = pc->transform[0];
                ni.transform[1] = pc->transform[1];
                ni.transform[2] = pc->transform[2];
                // Periodic target: same block or donor block on another process
                // For single-zone periodic, target is self
                ni.target_rank  = -1;  // same process
                ni.target_block = -1;  // will be resolved later
                ni.target_face  = face; // will be resolved from connectivity
            } else {
                // BC face — no halo exchange needed
                ni.active = false;
            }
        } else {
            // Internal split face — find the adjacent sub-block
            ni.active = true;
            ni.is_periodic = false;
            ni.transform[0] = 1;
            ni.transform[1] = 2;
            ni.transform[2] = 3;

            int neighbor_idx = find_split_neighbor(face, sub, all_decomp);
            if (neighbor_idx >= 0) {
                const SubBlock& ns = all_decomp[static_cast<std::size_t>(neighbor_idx)];
                ni.target_rank  = ns.assigned_rank;
                ni.target_block = neighbor_idx;
                // The neighbor's connecting face is the opposite side
                ni.target_face  = (face % 2 == 0) ? face + 1 : face - 1;
            } else {
                // Shouldn't happen — internal split should always have a neighbor
                ni.active = false;
            }
        }
    }
}

inline bool LocalBlock::is_external_face(int face, const SubBlock& sub,
                                          const Grid& full_zone) {
    Int nci_core = full_zone.ni_core - 1;
    Int ncj_core = full_zone.nj_core - 1;
    Int nck_core = full_zone.nk_core - 1;

    switch (face) {
    case 0: return sub.ci_min == 0;                // IMIN
    case 1: return sub.ci_max == nci_core - 1;     // IMAX
    case 2: return sub.cj_min == 0;                // JMIN
    case 3: return sub.cj_max == ncj_core - 1;     // JMAX
    case 4: return sub.ck_min == 0;                // KMIN
    case 5: return sub.ck_max == nck_core - 1;     // KMAX
    default: return false;
    }
}

inline const Connectivity* LocalBlock::find_periodic_at_face(
        const Grid& full_zone, int face, const SubBlock& /*sub*/) {
    // Use the existing Grid method
    return full_zone.find_periodic_connection(face);
}

inline int LocalBlock::find_split_neighbor(
        int face, const SubBlock& sub,
        const std::vector<SubBlock>& all_decomp) {

    for (std::size_t i = 0; i < all_decomp.size(); ++i) {
        const SubBlock& other = all_decomp[i];
        if (other.zone_id != sub.zone_id) continue;  // must be same zone

        switch (face) {
        case 0: // IMIN — neighbor must have ci_max == ci_min - 1
            if (other.ci_max == sub.ci_min - 1
                && other.cj_min == sub.cj_min && other.cj_max == sub.cj_max
                && other.ck_min == sub.ck_min && other.ck_max == sub.ck_max)
                return static_cast<int>(i);
            break;
        case 1: // IMAX — neighbor must have ci_min == ci_max + 1
            if (other.ci_min == sub.ci_max + 1
                && other.cj_min == sub.cj_min && other.cj_max == sub.cj_max
                && other.ck_min == sub.ck_min && other.ck_max == sub.ck_max)
                return static_cast<int>(i);
            break;
        case 2: // JMIN
            if (other.cj_max == sub.cj_min - 1
                && other.ci_min == sub.ci_min && other.ci_max == sub.ci_max
                && other.ck_min == sub.ck_min && other.ck_max == sub.ck_max)
                return static_cast<int>(i);
            break;
        case 3: // JMAX
            if (other.cj_min == sub.cj_max + 1
                && other.ci_min == sub.ci_min && other.ci_max == sub.ci_max
                && other.ck_min == sub.ck_min && other.ck_max == sub.ck_max)
                return static_cast<int>(i);
            break;
        case 4: // KMIN
            if (other.ck_max == sub.ck_min - 1
                && other.ci_min == sub.ci_min && other.ci_max == sub.ci_max
                && other.cj_min == sub.cj_min && other.cj_max == sub.cj_max)
                return static_cast<int>(i);
            break;
        case 5: // KMAX
            if (other.ck_min == sub.ck_max + 1
                && other.ci_min == sub.ci_min && other.ci_max == sub.ci_max
                && other.cj_min == sub.cj_min && other.cj_max == sub.cj_max)
                return static_cast<int>(i);
            break;
        }
    }

    return -1;
}

// ============================================================================
// Summary
// ============================================================================

inline void LocalBlock::print_summary() const {
    std::cout << "LocalBlock[" << block_id << "] zone=" << zone_id
              << " sub=" << sub_id
              << " dims=" << grid.nci << "x" << grid.ncj << "x" << grid.nck
              << " (core=" << nci_core() << "x" << ncj_core() << "x" << nck_core() << ")"
              << " ng=" << grid.ng << "\n";
    const char* face_names[] = {"IMIN","IMAX","JMIN","JMAX","KMIN","KMAX"};
    for (int f = 0; f < 6; ++f) {
        const auto& n = neighbors[f];
        if (n.active) {
            std::cout << "  " << face_names[f]
                      << " → block=" << n.target_block
                      << " face=" << n.target_face
                      << " rank=" << n.target_rank
                      << (n.is_periodic ? " [periodic]" : " [internal]") << "\n";
        } else {
            std::cout << "  " << face_names[f] << " → BC (no exchange)\n";
        }
    }
}
