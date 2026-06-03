#pragma once

#include "wcns_v2/utils/types.h"
#include "wcns_v2/grid/grid.h"
#include "wcns_v2/field/field.h"
#include "wcns_v2/parallel/decomposer.h"
#include <vector>

/// @file local_block.h
/// @brief LocalBlock — the fundamental runtime unit combining Grid + Field + adjacency.
///
/// Each MPI process owns one or more LocalBlocks.  A LocalBlock holds:
///   - grid  — geometry, metrics, ghost layers for this sub-domain
///   - field — flow variables on this sub-domain
///   - neighbors[6] — adjacency info for halo exchange
///
/// Construction supports two paths:
///   1. from_full_zone()  — the entire zone maps to this block (no splitting)
///   2. from_sub_zone()  — extract a sub-region from a zone (after decomposition)

/// Information about the neighbor on one face, used for halo exchange.
struct NeighborInfo {
    int  target_rank;     ///< MPI rank of the neighbor process (-1 = same process)
    int  target_block;    ///< Global block ID of the neighbor
    int  target_face;     ///< Which face of the neighbor connects to us (0..5)
    int  transform[3];    ///< Direction mapping (same as CGNS transform)
    bool active;          ///< false = physical boundary (no halo exchange needed)
    bool is_periodic;     ///< true = periodic connection (even if same block)

    NeighborInfo();
};

/// A single local sub-domain with geometry, flow field, and connectivity.
class LocalBlock {
public:
    // ---- Identity ----
    int block_id;              ///< Global unique block ID
    int zone_id;               ///< Original zone index
    int sub_id;                ///< Sub-block index within zone (0 if not split)

    // ---- Core data ----
    Grid  grid;                ///< Geometry and metrics (includes ghost layers)
    Field field;               ///< Flow variables (allocated to match grid)

    // ---- Adjacency (6 faces: 0=IMIN, 1=IMAX, 2=JMIN, 3=JMAX, 4=KMIN, 5=KMAX) ----
    NeighborInfo neighbors[6];

    LocalBlock();

    // ========================================================================
    // Construction from a full zone (no splitting)
    // ========================================================================

    /// Build a LocalBlock covering the entire zone.
    /// The zone must already have ghost layers extended and metrics computed.
    /// BC and 1-to-1 connectivity from the zone are mapped to NeighborInfo.
    static LocalBlock from_full_zone(const Grid& full_zone, int block_id,
                                      int zone_id, int my_rank,
                                      const std::vector<SubBlock>& all_decomp);

    // ========================================================================
    // Construction from a sub-zone (after Cartesian splitting)
    // ========================================================================

    /// Build a LocalBlock from a sub-region of a full zone.
    /// Extracts node coordinates and re-computes cell centers and metrics.
    /// Builds NeighborInfo for both external boundaries and internal splits.
    static LocalBlock from_sub_zone(const Grid& full_zone,
                                     const SubBlock& sub, int block_id,
                                     int my_rank,
                                     const std::vector<SubBlock>& all_decomp);

    // ========================================================================
    // Queries
    // ========================================================================

    /// True if the given face neighbors a block on a different MPI rank.
    bool neighbor_is_remote(int face) const;

    /// MPI rank of the neighbor on the given face (-1 if same process or BC).
    int  neighbor_rank(int face) const;

    /// Number of core cells in each direction.
    Int nci_core() const;
    Int ncj_core() const;
    Int nck_core() const;

    /// Print a summary of this block.
    void print_summary() const;

private:
    /// Build NeighborInfo for all 6 faces from existing connectivity data.
    /// Handles BC faces, periodic connections, and internal split interfaces.
    static void build_neighbors(LocalBlock& lb, const Grid& full_zone,
                                 const SubBlock& sub, int my_rank,
                                 const std::vector<SubBlock>& all_decomp);

    /// Determine if a face of the sub-block touches the boundary of the
    /// original zone (i.e., is an external face, not an internal split).
    static bool is_external_face(int face, const SubBlock& sub,
                                  const Grid& full_zone);

    /// Find a periodic connection on the original zone at the given face.
    static const Connectivity* find_periodic_at_face(const Grid& full_zone,
                                                      int face,
                                                      const SubBlock& sub);

    /// Find an adjacent sub-block on the given face (internal split interface).
    /// Returns index into all_decomp, or -1 if not found.
    static int find_split_neighbor(int face, const SubBlock& sub,
                                    const std::vector<SubBlock>& all_decomp);
};

#include "local_block.hxx"
