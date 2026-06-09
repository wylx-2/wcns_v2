#pragma once

#include "utils/types.h"
#include "core/config.h"
#include "grid/grid.h"
#include "parallel/parallel_env.h"
#include "parallel/decomposer.h"
#include "parallel/local_block.h"
#include "parallel/halo_exchange.h"
#include "parallel/flux_halo_exchange.h"
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

    /// Access the FluxHaloExchange object for a given local block index.
    /// Used by ViscidRHS for face-product exchange.
    FluxHaloExchange& flux_halo_ex(Int idx);

    /// Exchange velocity/temperature gradient ghost cells for all local blocks.
    /// Exchanges the 12 gradient arrays (du_dx, du_dy, ..., dT_dz) stored in Field.
    void exchange_gradient_halos(std::vector<LocalBlock>& blocks);

    /// Exchange cell-center viscous flux ghost cells for all local blocks.
    /// Exchanges vis_x, vis_y, vis_z (15 arrays total) stored in Field.
    void exchange_viscous_flux_halos(std::vector<LocalBlock>& blocks);

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
