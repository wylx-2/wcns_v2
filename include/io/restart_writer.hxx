#pragma once

#include "restart_writer.h"
#include "parallel/local_block.h"
#include "parallel/parallel_env.h"
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

// ============================================================================
// Helper: restart filename generation
// ============================================================================

namespace {

inline std::string restart_basename(const std::string& dir,
                                     Int iter) {
    std::ostringstream oss;
    if (!dir.empty()) oss << dir << "/";
    oss << "restart_" << std::setw(6) << std::setfill('0') << iter;
    return oss.str();
}

inline std::string restart_rank_filename(const std::string& dir,
                                          Int iter, int rank) {
    std::ostringstream oss;
    if (!dir.empty()) oss << dir << "/";
    oss << "restart_" << std::setw(6) << std::setfill('0') << iter
        << "_r" << rank << ".bin";
    return oss.str();
}

}  // anonymous namespace

// ============================================================================
// RestartWriter::write
// ============================================================================

inline void RestartWriter::write(const std::vector<LocalBlock>& blocks,
                                  const Config& cfg, Int iter, Real time) {
    int my_rank = ParallelEnv::rank();
    int nprocs  = ParallelEnv::size();

    // ---- Per-rank data file ----
    {
        auto fname = restart_rank_filename(cfg.output_dir, iter, my_rank);
        std::ofstream os(fname, std::ios::binary);
        if (!os.is_open()) {
            throw std::runtime_error("RestartWriter: cannot open \""
                                     + fname + "\" for writing");
        }

        Int n_local = static_cast<Int>(blocks.size());
        os.write(reinterpret_cast<const char*>(&n_local), sizeof(Int));

        for (const auto& lb : blocks) {
            Int nci = lb.field.ni();
            Int ncj = lb.field.nj();
            Int nck = lb.field.nk();
            Int ng  = lb.grid.ng;
            Int bid = lb.block_id;
            Int zid = lb.zone_id;

            // Metadata
            os.write(reinterpret_cast<const char*>(&bid), sizeof(Int));
            os.write(reinterpret_cast<const char*>(&zid), sizeof(Int));
            os.write(reinterpret_cast<const char*>(&nci), sizeof(Int));
            os.write(reinterpret_cast<const char*>(&ncj), sizeof(Int));
            os.write(reinterpret_cast<const char*>(&nck), sizeof(Int));
            os.write(reinterpret_cast<const char*>(&ng),  sizeof(Int));

            // Conservative arrays (flat binary dump in memory order)
            os.write(reinterpret_cast<const char*>(lb.field.cons.rho.data()),
                     lb.field.cons.rho.size()   * sizeof(Real));
            os.write(reinterpret_cast<const char*>(lb.field.cons.rhou.data()),
                     lb.field.cons.rhou.size()  * sizeof(Real));
            os.write(reinterpret_cast<const char*>(lb.field.cons.rhov.data()),
                     lb.field.cons.rhov.size()  * sizeof(Real));
            os.write(reinterpret_cast<const char*>(lb.field.cons.rhow.data()),
                     lb.field.cons.rhow.size()  * sizeof(Real));
            os.write(reinterpret_cast<const char*>(lb.field.cons.rhoE.data()),
                     lb.field.cons.rhoE.size()  * sizeof(Real));
        }

        os.close();
    }

    // ---- Index file ----
    // MPI collectives must be called by ALL ranks, not just rank 0.
    Int n_local = static_cast<Int>(blocks.size());

    // Pack local block IDs
    std::vector<Int> local_bids;
    local_bids.reserve(static_cast<std::size_t>(n_local));
    for (const auto& lb : blocks)
        local_bids.push_back(lb.block_id);

    // Gather counts and block IDs from all ranks to rank 0
    std::vector<Int> all_counts;
    std::vector<Int> all_block_ids;

    if (nprocs > 1) {
        all_counts.resize(static_cast<std::size_t>(nprocs));
        MPI_Gather(&n_local, 1, MPI_INT, all_counts.data(), 1,
                   MPI_INT, 0, ParallelEnv::communicator());

        std::vector<Int> recv_counts(static_cast<std::size_t>(nprocs));
        std::vector<Int> recv_displs(static_cast<std::size_t>(nprocs));
        Int total_bids = 0;
        if (my_rank == 0) {
            for (int r = 0; r < nprocs; ++r) {
                recv_counts[static_cast<std::size_t>(r)] =
                    all_counts[static_cast<std::size_t>(r)];
                recv_displs[static_cast<std::size_t>(r)] = total_bids;
                total_bids += recv_counts[static_cast<std::size_t>(r)];
            }
            all_block_ids.resize(static_cast<std::size_t>(total_bids));
        }

        MPI_Gatherv(local_bids.data(), n_local, MPI_INT,
                    all_block_ids.data(), recv_counts.data(),
                    recv_displs.data(), MPI_INT, 0,
                    ParallelEnv::communicator());
    } else {
        all_counts.push_back(n_local);
        all_block_ids = local_bids;
    }

    // Only rank 0 writes the index file
    if (my_rank == 0) {
        auto fname = restart_basename(cfg.output_dir, iter) + ".bin";
        std::ofstream os(fname, std::ios::binary);
        if (!os.is_open()) {
            throw std::runtime_error("RestartWriter: cannot open index \""
                                     + fname + "\" for writing");
        }

        os.write(reinterpret_cast<const char*>(&kMagic),   sizeof(std::uint32_t));
        os.write(reinterpret_cast<const char*>(&kVersion), sizeof(std::uint32_t));
        os.write(reinterpret_cast<const char*>(&iter),     sizeof(Int));
        os.write(reinterpret_cast<const char*>(&time),     sizeof(Real));

        Int nprocs_val = nprocs;
        os.write(reinterpret_cast<const char*>(&nprocs_val), sizeof(Int));

        Int n_blocks_total = static_cast<Int>(all_block_ids.size());
        os.write(reinterpret_cast<const char*>(&n_blocks_total), sizeof(Int));

        Int rank_offset = 0;
        for (int r = 0; r < nprocs; ++r) {
            Int nb = all_counts[static_cast<std::size_t>(r)];
            for (Int b = 0; b < nb; ++b) {
                Int bid = all_block_ids[static_cast<std::size_t>(rank_offset + b)];
                Int rnk = r;
                os.write(reinterpret_cast<const char*>(&bid), sizeof(Int));
                os.write(reinterpret_cast<const char*>(&rnk), sizeof(Int));
            }
            rank_offset += nb;
        }

        os.close();

        std::cout << "Restart written: iter=" << iter << ", time=" << time
                  << ", blocks=" << n_blocks_total << "\n";
    }
}

// ============================================================================
// RestartWriter::read
// ============================================================================

inline std::pair<Int, Real> RestartWriter::read(
        const std::string& basename,
        std::vector<LocalBlock>& blocks,
        const Config& cfg) {

    int my_rank = ParallelEnv::rank();

    auto fname = basename + "_r" + std::to_string(my_rank) + ".bin";
    std::ifstream is(fname, std::ios::binary);
    if (!is.is_open()) {
        throw std::runtime_error("RestartWriter: cannot open \""
                                 + fname + "\" for reading");
    }

    Int n_local;
    is.read(reinterpret_cast<char*>(&n_local), sizeof(Int));
    if (!is) {
        throw std::runtime_error("RestartWriter: failed to read header from \""
                                 + fname + "\"");
    }

    if (n_local != static_cast<Int>(blocks.size())) {
        throw std::runtime_error(
            "RestartWriter: block count mismatch in \""
            + fname + "\" (expected " + std::to_string(blocks.size())
            + ", got " + std::to_string(n_local) + ")");
    }

    // Read iteration/time from index file (rank 0)
    // For simplicity, each rank reads the index too
    Int saved_iter = 0;
    Real saved_time = 0.0;

    {
        std::string idx_name = basename + ".bin";
        std::ifstream idx(idx_name, std::ios::binary);
        if (idx.is_open()) {
            std::uint32_t magic, version;
            idx.read(reinterpret_cast<char*>(&magic),   sizeof(std::uint32_t));
            idx.read(reinterpret_cast<char*>(&version), sizeof(std::uint32_t));
            if (magic == kMagic) {
                idx.read(reinterpret_cast<char*>(&saved_iter), sizeof(Int));
                idx.read(reinterpret_cast<char*>(&saved_time), sizeof(Real));
            }
        }
    }

    // Broadcast iter/time from rank 0 in case other ranks couldn't read index
    if (ParallelEnv::is_parallel()) {
        MPI_Bcast(&saved_iter, 1, MPI_INT, 0, ParallelEnv::communicator());
        MPI_Bcast(&saved_time, 1, MPI_DOUBLE, 0, ParallelEnv::communicator());
    }

    for (Int ib = 0; ib < n_local; ++ib) {
        auto& lb = blocks[static_cast<std::size_t>(ib)];

        Int bid, zid, nci, ncj, nck, ng;
        is.read(reinterpret_cast<char*>(&bid), sizeof(Int));
        is.read(reinterpret_cast<char*>(&zid), sizeof(Int));
        is.read(reinterpret_cast<char*>(&nci), sizeof(Int));
        is.read(reinterpret_cast<char*>(&ncj), sizeof(Int));
        is.read(reinterpret_cast<char*>(&nck), sizeof(Int));
        is.read(reinterpret_cast<char*>(&ng),  sizeof(Int));

        if (!is) {
            throw std::runtime_error("RestartWriter: failed to read block "
                                     + std::to_string(ib) + " metadata from \""
                                     + fname + "\"");
        }

        // Verify dimensions match
        if (nci != lb.field.ni() || ncj != lb.field.nj() || nck != lb.field.nk()) {
            throw std::runtime_error(
                "RestartWriter: dimension mismatch for block " + std::to_string(ib)
                + " (field: " + std::to_string(lb.field.ni()) + "x"
                + std::to_string(lb.field.nj()) + "x"
                + std::to_string(lb.field.nk())
                + ", restart: " + std::to_string(nci) + "x"
                + std::to_string(ncj) + "x" + std::to_string(nck) + ")");
        }

        // Read conservative arrays
        is.read(reinterpret_cast<char*>(lb.field.cons.rho.data()),
                lb.field.cons.rho.size()   * sizeof(Real));
        is.read(reinterpret_cast<char*>(lb.field.cons.rhou.data()),
                lb.field.cons.rhou.size()  * sizeof(Real));
        is.read(reinterpret_cast<char*>(lb.field.cons.rhov.data()),
                lb.field.cons.rhov.size()  * sizeof(Real));
        is.read(reinterpret_cast<char*>(lb.field.cons.rhow.data()),
                lb.field.cons.rhow.size()  * sizeof(Real));
        is.read(reinterpret_cast<char*>(lb.field.cons.rhoE.data()),
                lb.field.cons.rhoE.size()  * sizeof(Real));

        if (!is) {
            throw std::runtime_error("RestartWriter: failed to read block "
                                     + std::to_string(ib) + " data from \""
                                     + fname + "\"");
        }

        // Recompute primitive variables from conservative
        lb.field.cons_to_prim(cfg.gamma);
    }

    is.close();

    if (my_rank == 0) {
        std::cout << "Restart loaded: iter=" << saved_iter
                  << ", time=" << saved_time << "\n";
    }

    return {saved_iter, saved_time};
}
