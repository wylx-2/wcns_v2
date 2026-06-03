#pragma once

#include "parallel_manager.h"
#include "wcns_v2/io/cgns_reader.h"
#include <iostream>

// ============================================================================
// Initialize
// ============================================================================

inline void ParallelManager::initialize(
        const std::string& grid_file, const Config& cfg,
        std::vector<LocalBlock>& blocks) {

    rank_   = ParallelEnv::rank();
    nprocs_ = ParallelEnv::size();
    blocks.clear();

    // ---- Step 1: All processes read the CGNS file ----
    std::vector<Grid> zones;
    {
        CGNSReader reader;
        reader.open(grid_file);
        reader.read_all(zones, /*ng=*/cfg.ng);
        reader.close();
    }
    Int nz = static_cast<Int>(zones.size());

    // ---- Step 2: Extend ghost layers on all zones ----
    for (auto& z : zones) {
        z.extend_ghost_layers();
    }

    // ---- Step 2b: Fix interface ghost nodes using 1-to-1 connectivity ----
    // For inter-zone connections, ghost node coordinates should come from
    // the donor zone's interior nodes, not from linear extrapolation.
    for (auto& z : zones) {
        for (int face = 0; face < 6; ++face) {
            const Connectivity* conn = z.find_face_connection(face);
            if (!conn) continue;
            // Within-zone connections are already handled correctly
            // by fill_ghost_face_periodic during extend_ghost_layers.
            if (conn->donor_name == z.name) continue;

            // Find donor zone by name
            Grid* donor = nullptr;
            for (auto& d : zones) {
                if (d.name == conn->donor_name) { donor = &d; break; }
            }
            if (donor) {
                z.fix_interface_ghost(face, *donor, *conn);
            }
        }
    }

    // ---- Step 2c: Recompute cell centers (ghost nodes changed in 2b) ----
    // Cell centers were originally computed inside extend_ghost_layers()
    // from extrapolated ghost nodes.  After fix_interface_ghost corrected
    // the ghost node coordinates, the cell-center coordinates and volumes
    // must be recomputed.
    for (auto& z : zones) {
        z.compute_cell_centers();
        z.compute_cell_volumes();
    }

    // ---- Step 2d: Compute metrics on all zones ----
    for (auto& z : zones) {
        z.compute_metrics();
        z.compute_face_metrics();
    }

    // ---- Step 3: Domain decomposition (deterministic, all procs same result) ----
    std::vector<SubBlock> decomp = BlockDecomposer::decompose(zones, nprocs_, cfg.ng);
    total_blocks_ = static_cast<int>(decomp.size());

    // ---- Build zone name → (rank, block_id) mapping for flux halo exchange ----
    zone_to_block_.clear();
    for (std::size_t i = 0; i < decomp.size(); ++i) {
        const SubBlock& sb = decomp[i];
        // For full-zone blocks (sub_id == 0, not split), use the original zone name.
        // The zone_id indexes into the original zones vector.
        if (sb.sub_id == 0) {
            // Check if this zone was split
            bool was_split = false;
            for (const auto& other : decomp) {
                if (other.zone_id == sb.zone_id && other.sub_id != sb.sub_id) {
                    was_split = true;
                    break;
                }
            }
            if (!was_split) {
                const std::string& zname = zones[static_cast<std::size_t>(sb.zone_id)].name;
                zone_to_block_[zname] = {sb.assigned_rank, static_cast<int>(i)};
            }
        }
    }

    // ---- Step 4: Build LocalBlocks assigned to this process ----
    for (std::size_t i = 0; i < decomp.size(); ++i) {
        const SubBlock& sb = decomp[i];
        if (sb.assigned_rank != rank_) continue;

        int block_id = static_cast<int>(i);

        // Check if any other SubBlock shares the same zone_id (meaning this zone was split)
        bool zone_was_split = false;
        for (const auto& other : decomp) {
            if (other.zone_id == sb.zone_id && other.sub_id != sb.sub_id) {
                zone_was_split = true;
                break;
            }
        }

        LocalBlock lb;
        if (!zone_was_split) {
            // Full zone — just deep-copy the already-constructed zone
            lb = LocalBlock::from_full_zone(
                zones[static_cast<std::size_t>(sb.zone_id)],
                block_id, sb.zone_id, rank_, decomp);
        } else {
            // Sub-zone — extract from the full zone
            lb = LocalBlock::from_sub_zone(
                zones[static_cast<std::size_t>(sb.zone_id)],
                sb, block_id, rank_, decomp);
        }

        blocks.push_back(std::move(lb));
    }

    // ---- Step 5: Setup HaloExchange for each local block ----
    halo_ex_.resize(blocks.size());
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        halo_ex_[i].setup(blocks[i]);
    }

    // ---- Step 5b: Setup FluxHaloExchange for each local block ----
    flux_halo_ex_.resize(blocks.size());
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        flux_halo_ex_[i].setup(blocks[i], zone_to_block_);
    }

    // ---- Print summary ----
    if (rank_ == 0) {
        std::cout << "\n=== Parallel Setup Complete ===\n"
                  << "  Total zones:    " << nz << "\n"
                  << "  Total sub-blocks: " << total_blocks_ << "\n"
                  << "  Processes:      " << nprocs_ << "\n"
                  << "  Ghost layers:   " << cfg.ng << "\n"
                  << "================================\n" << std::endl;
    }

    for (int r = 0; r < nprocs_; ++r) {
        ParallelEnv::barrier();
        if (r == rank_) {
            std::cout << "[Rank " << rank_ << "] owns " << blocks.size()
                      << " block(s):\n";
            for (const auto& lb : blocks) {
                lb.print_summary();
            }
            std::cout << std::flush;
        }
        ParallelEnv::barrier();
    }
}

// ============================================================================
// Exchange all halos
// ============================================================================

inline void ParallelManager::exchange_all_halos(std::vector<LocalBlock>& blocks) {
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        // Exchange primitive variables
        std::vector<MultiArray3D<Real>*> arrays = {
            &blocks[i].field.prim.rho,
            &blocks[i].field.prim.u,
            &blocks[i].field.prim.v,
            &blocks[i].field.prim.w,
            &blocks[i].field.prim.p
        };
        halo_ex_[i].exchange_multi(arrays, blocks[i]);
    }
}

// ============================================================================
// Exchange flux halos
// ============================================================================

inline void ParallelManager::exchange_flux_halos(std::vector<LocalBlock>& blocks) {
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        flux_halo_ex_[i].exchange(blocks[i], blocks);
    }
}

// ============================================================================
// Global reductions
// ============================================================================

inline Real ParallelManager::global_max(Real local_val) {
    if (!ParallelEnv::is_parallel()) return local_val;
    Real global = local_val;
    MPI_Allreduce(&local_val, &global, 1, MPI_DOUBLE, MPI_MAX,
                  ParallelEnv::communicator());
    return global;
}

inline Real ParallelManager::global_min(Real local_val) {
    if (!ParallelEnv::is_parallel()) return local_val;
    Real global = local_val;
    MPI_Allreduce(&local_val, &global, 1, MPI_DOUBLE, MPI_MIN,
                  ParallelEnv::communicator());
    return global;
}

inline Real ParallelManager::global_sum(Real local_val) {
    if (!ParallelEnv::is_parallel()) return local_val;
    Real global = local_val;
    MPI_Allreduce(&local_val, &global, 1, MPI_DOUBLE, MPI_SUM,
                  ParallelEnv::communicator());
    return global;
}

inline Real ParallelManager::broadcast(Real val) {
    if (!ParallelEnv::is_parallel()) return val;
    // Only rank 0's value is meaningful; all ranks participate in the broadcast
    MPI_Bcast(&val, 1, MPI_DOUBLE, 0, ParallelEnv::communicator());
    return val;
}
