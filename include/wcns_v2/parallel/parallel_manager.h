#pragma once

#include "wcns_v2/utils/types.h"
#include "wcns_v2/core/config.h"
#include "wcns_v2/grid/grid.h"
#include "wcns_v2/parallel/parallel_env.h"
#include "wcns_v2/parallel/decomposer.h"
#include "wcns_v2/parallel/local_block.h"
#include "wcns_v2/parallel/halo_exchange.h"
#include "wcns_v2/parallel/flux_halo_exchange.h"
#include <map>
#include <string>
#include <vector>

/// @file parallel_manager.h
/// @brief Top-level parallel manager — drives domain decomposition, LocalBlock
///        construction, halo exchange, and global reductions.

class ParallelManager {
public:
    /// Initialize the parallel environment and set up domain decomposition.
    ///
    /// All processes read the same CGNS file independently, extend ghost layers,
    /// compute metrics, and run the same decomposition algorithm.  Each process
    /// then builds only the LocalBlocks assigned to it.
    ///
    /// @param grid_file   Path to CGNS grid file
    /// @param cfg         Solver configuration (provides ng)
    /// @param blocks      [out] LocalBlocks owned by this process
    void initialize(const std::string& grid_file,
                    const Config& cfg,
                    std::vector<LocalBlock>& blocks);

    /// Exchange halo (ghost) data for all arrays of all local blocks.
    /// This is the primary communication call during time stepping.
    void exchange_all_halos(std::vector<LocalBlock>& blocks);

    /// Exchange inviscid face-flux halos for all local blocks.
    ///
    /// After the Riemann solver computes face fluxes, this exchanges ng+1
    /// face-flux slices at inter-zone connectivity boundaries so that
    /// ghost-face fluxes are consistent with the neighbor block's interior
    /// face fluxes.  Must be called after Riemann solver and before the
    /// 6-point centered difference in InviscidRHS.
    void exchange_flux_halos(std::vector<LocalBlock>& blocks);

    // ========================================================================
    // Global reductions (for convergence checks, CFL computation, etc.)
    // ========================================================================

    /// Global maximum across all MPI processes.
    static Real global_max(Real local_val);

    /// Global minimum across all MPI processes.
    static Real global_min(Real local_val);

    /// Global sum across all MPI processes.
    static Real global_sum(Real local_val);

    /// Broadcast a value from the master rank to all processes.
    static Real broadcast(Real val);

    /// Number of sub-blocks in the global decomposition.
    int total_blocks() const { return total_blocks_; }

private:
    int rank_;
    int nprocs_;
    int total_blocks_ = 0;

    /// One HaloExchange per local block.
    std::vector<HaloExchange> halo_ex_;

    /// One FluxHaloExchange per local block.
    std::vector<FluxHaloExchange> flux_halo_ex_;

    /// Map from original zone name → (rank, global_block_id).
    /// Built during initialize() for resolving inter-zone connectivity.
    std::map<std::string, std::pair<int,int>> zone_to_block_;
};

#include "parallel_manager.hxx"
