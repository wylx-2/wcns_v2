#pragma once

#include "wcns_v2/utils/types.h"
#include "wcns_v2/grid/grid.h"
#include <vector>

/// @file decomposer.h
/// @brief Domain decomposition for structured multi-block grids.
///
/// Distributes zones (or sub-zones) across MPI processes with load balancing.
/// If the number of processes exceeds the number of zones, the largest zones
/// are split using Cartesian decomposition while respecting a minimum block
/// size constraint (>= 2*ng in each direction for valid ghost layers).

/// Describes a sub-block produced by domain decomposition.
/// All index ranges refer to the *original zone's core cells* (0-based,
/// BEFORE ghost extension).  A sub-block that is the entire zone has
/// ci_min=0, ci_max=nci_core-1, etc.
struct SubBlock {
    Int zone_id;             ///< Index into the original zone list
    Int sub_id;              ///< Sub-block number within this zone (0 if not split)
    Int assigned_rank;       ///< MPI process that owns this sub-block

    /// Cell-index range in the original zone (0-based, inclusive).
    /// The zone's core cell count is nci_core = ni_core-1, so valid
    /// ranges are [0 .. ni_core-2].
    Int ci_min, ci_max;      ///< i-direction cell range
    Int cj_min, cj_max;      ///< j-direction cell range
    Int ck_min, ck_max;      ///< k-direction cell range
};

/// Decompose a list of zones across processes.
class BlockDecomposer {
public:
    /// Compute the decomposition.
    ///
    /// @param zones    Input zones (core dimensions only, before ghost extension)
    /// @param nprocs   Number of MPI processes
    /// @param ng       Ghost layer count (minimum sub-block dimension >= 2*ng+1 cells)
    /// @return         Vector of SubBlock describing the assignment
    static std::vector<SubBlock> decompose(
        const std::vector<Grid>& zones, int nprocs, Int ng);

private:
    /// Greedy assignment: assign zones to processes by descending cell count.
    static void assign_greedy(const std::vector<Int>& cell_counts,
                               int nprocs,
                               std::vector<int>& proc_loads,
                               std::vector<int>& assignments);

    /// Find optimal integer split factors for a grid of size (ni,nj,nk) cells,
    /// splitting into approximately `target` sub-blocks.
    /// Constraint: each sub-block must have >= min_cells cells in each direction.
    static void find_split_factors(Int nc_i, Int nc_j, Int nc_k,
                                    int target, Int min_cells,
                                    int& si, int& sj, int& sk);

    /// Compute the total cell count of a zone (core, no ghost).
    static Int cell_count_core(const Grid& g);

    /// Compute the split factors and add sub-blocks to the result.
    /// @param n_blocks  How many sub-blocks to split this zone into
    /// @param nprocs    Total number of MPI processes (for round-robin assignment)
    /// @param next_rank Next available rank slot (incremented for each sub-block)
    static void split_zone(const Grid& zone, int zone_id,
                            int n_blocks, int nprocs, Int ng,
                            std::vector<SubBlock>& result,
                            int& next_rank);
};

#include "decomposer.hxx"
