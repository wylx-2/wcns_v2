#pragma once

#include "wcns_v2/core/config.h"
#include "wcns_v2/parallel/local_block.h"
#include "wcns_v2/grid/boundary_condition.h"

/// @file bc_applier.h
/// @brief Apply boundary conditions to ghost cells.
///
/// Three-stage approach:
///   1. Face ghost — fill cells where exactly one index is in the ghost range,
///      based on the BC patch type for that face.
///   2. Edge ghost — fill cells where two indices are in the ghost range
///      (intersection of two faces), by averaging the two face values.
///   3. Corner ghost — fill cells where all three indices are in the ghost
///      range (intersection of three faces), by averaging the three face values.
///
/// For faces with active parallel neighbors, ghost data comes from
/// HaloExchange instead — those faces are skipped during BC application.

class BoundaryConditionApplier {
public:
    /// Apply all boundary conditions: face → edge → corner ghost filling.
    /// Only processes faces where neighbors[f].active == false (BC faces).
    static void apply_all(LocalBlock& lb, const Config& cfg);

private:
    // ---- Stage 1: face ghost ----
    /// Fill face ghost cells (including edge/corner overlaps — these will be
    /// corrected in stages 2 and 3).
    static void apply_face_ghost(LocalBlock& lb, const Config& cfg);

    /// Apply a specific BC type to a rectangular patch on one face.
    /// (j0,j1,k0,k1) are the cell index ranges in the tangential directions.
    static void apply_bc_on_face(int face, BCType bc_type,
                                  LocalBlock& lb, const Config& cfg,
                                  Int j0, Int j1, Int k0, Int k1);

    // ---- Individual BC type implementations ----
    static void apply_farfield(int face, LocalBlock& lb, const Config& cfg,
                               Int j0, Int j1, Int k0, Int k1);
    static void apply_wall_noslip(int face, LocalBlock& lb,
                                  Int j0, Int j1, Int k0, Int k1);
    static void apply_wall_slip(int face, LocalBlock& lb,
                                Int j0, Int j1, Int k0, Int k1);
    static void apply_symmetry(int face, LocalBlock& lb,
                               Int j0, Int j1, Int k0, Int k1);
    static void apply_inflow(int face, LocalBlock& lb, const Config& cfg,
                             Int j0, Int j1, Int k0, Int k1);
    static void apply_outflow(int face, LocalBlock& lb,
                              Int j0, Int j1, Int k0, Int k1);

    // ---- Stage 2 & 3: edge and corner ----
    static void apply_edge_ghost(LocalBlock& lb);
    static void apply_corner_ghost(LocalBlock& lb);

    // ---- Helpers ----
    /// Convert a BCPatch node range to cell index range in the face-tangential
    /// directions.  Returns (r0,r1,s0,s1) which are interpreted as (j0,j1,k0,k1)
    /// for I-faces, or (i0,i1,k0,k1) for J-faces, etc.
    static void patch_cell_range(int face, const BCPatch& patch, Int ng,
                                 Int nci_core, Int ncj_core, Int nck_core,
                                 Int& r0, Int& r1, Int& s0, Int& s1);

    /// Mirror index for wall BC: ghost layer d maps to interior layer
    /// (even reflection about the face).
    static Int mirror_index(int face, Int ghost_idx, Int ng);

    /// Check if a face is a BC face (not a parallel neighbor).
    static bool is_bc_face(const LocalBlock& lb, int face);
};

#include "bc_applier.hxx"
