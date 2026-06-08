#pragma once

#include "solution_writer.h"
#include "wcns_v2/parallel/local_block.h"
#include "wcns_v2/parallel/parallel_env.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

// ============================================================================
// Helper: make_filename
// ============================================================================

inline std::string make_filename(const std::string& dir,
                                  const std::string& prefix,
                                  Int iter,
                                  const std::string& suffix) {
    std::ostringstream oss;
    if (!dir.empty()) oss << dir << "/";
    oss << prefix << "_" << std::setw(6) << std::setfill('0') << iter << suffix;
    return oss.str();
}

// ============================================================================
// Per-block metadata for gather (used internally by TecplotWriter)
// ============================================================================

struct BlockMeta {
    Int block_id;
    Int ni_core;    // interior cells in i
    Int nj_core;    // interior cells in j
    Int nk_core;    // interior cells in k
    Int data_size;  // number of Reals = ni_core * nj_core * nk_core * 9
};

// ============================================================================
// TecplotWriter::write
// ============================================================================

inline void TecplotWriter::write(const std::vector<LocalBlock>& blocks,
                                  const Config& cfg,
                                  const std::string& filename,
                                  Int iter, Real time) {

    const Int ng = (blocks.empty()) ? 0 : blocks[0].grid.ng;

    // ---- Step 1: build metadata and pack local data ----
    std::vector<BlockMeta> local_meta;
    local_meta.reserve(blocks.size());

    // First pass: compute sizes and build metadata
    Int local_total_size = 0;
    for (const auto& lb : blocks) {
        BlockMeta m;
        m.block_id = lb.block_id;
        m.ni_core  = lb.field.ni() - 2 * ng;
        m.nj_core  = lb.field.nj() - 2 * ng;
        m.nk_core  = lb.field.nk() - 2 * ng;
        m.data_size = m.ni_core * m.nj_core * m.nk_core * 9;
        local_meta.push_back(m);
        local_total_size += m.data_size;
    }

    Int n_local_blocks = static_cast<Int>(local_meta.size());

    // Second pass: pack floating-point data
    // Order: x, y, z, rho, u, v, w, p per interior cell, block by block
    std::vector<Real> local_data;
    local_data.reserve(static_cast<std::size_t>(local_total_size));

    for (std::size_t ib = 0; ib < blocks.size(); ++ib) {
        const auto& lb = blocks[ib];
        const auto& f = lb.field;
        const auto& g = lb.grid;

        const Int i0 = ng;
        const Int i1 = f.ni() - 1 - ng;
        const Int j0 = ng;
        const Int j1 = f.nj() - 1 - ng;
        const Int k0 = ng;
        const Int k1 = f.nk() - 1 - ng;

        Real gm2 = cfg.gamma * cfg.Mach * cfg.Mach;  // 1 / eos_factor

        for (Int k = k0; k <= k1; ++k) {
        for (Int j = j0; j <= j1; ++j) {
        for (Int i = i0; i <= i1; ++i) {
            local_data.push_back(g.cell_x(i,j,k));
            local_data.push_back(g.cell_y(i,j,k));
            local_data.push_back(g.cell_z(i,j,k));
            local_data.push_back(f.prim.rho(i,j,k));
            local_data.push_back(f.prim.u(i,j,k));
            local_data.push_back(f.prim.v(i,j,k));
            local_data.push_back(f.prim.w(i,j,k));
            local_data.push_back(f.prim.p(i,j,k));
            local_data.push_back(f.prim.p(i,j,k) * gm2 / f.prim.rho(i,j,k));
        }}}
    }

    // ---- Step 2: gather to rank 0 (or direct write in serial) ----

    if (!ParallelEnv::is_parallel()) {
        // Serial: rank 0 writes directly
        std::ofstream os(filename);
        if (!os.is_open()) {
            throw std::runtime_error("TecplotWriter: cannot open \"" +
                                     filename + "\" for writing");
        }

        os << std::scientific << std::setprecision(6);
        os << "TITLE = \"WCNS v2 Solution — iter=" << iter
           << ", time=" << time << "\"\n";
        os << "VARIABLES = \"X\", \"Y\", \"Z\", "
           << "\"rho\", \"u\", \"v\", \"w\", \"p\", \"T\"\n";

        const Real* ptr = local_data.data();
        for (std::size_t ib = 0; ib < blocks.size(); ++ib) {
            const BlockMeta& m = local_meta[ib];
            os << "ZONE T=\"block_" << m.block_id
               << "\", I=" << m.ni_core
               << ", J=" << m.nj_core
               << ", K=" << m.nk_core
               << ", DATAPACKING=POINT\n";

            Int n_vals = m.data_size;
            for (Int v = 0; v < n_vals; v += 9) {
                os << " "  << ptr[0] << " " << ptr[1] << " " << ptr[2]
                   << " "  << ptr[3] << " " << ptr[4] << " " << ptr[5]
                   << " "  << ptr[6] << " " << ptr[7] << " " << ptr[8] << "\n";
                ptr += 9;
            }
        }
        os.close();
        return;
    }

    // ---- Parallel: MPI gather to rank 0 ----

    int my_rank = ParallelEnv::rank();
    int nprocs  = ParallelEnv::size();
    MPI_Comm comm = ParallelEnv::communicator();

    // --- Gather block counts ---
    std::vector<Int> nblocks_per_rank(static_cast<std::size_t>(nprocs));
    MPI_Gather(&n_local_blocks, 1, MPI_INT,
               nblocks_per_rank.data(), 1, MPI_INT, 0, comm);

    // --- Gather metadata ---
    // Pack local metadata into Int array: {block_id, ni_core, nj_core, nk_core, data_size} per block
    std::vector<Int> local_meta_buf;
    local_meta_buf.reserve(static_cast<std::size_t>(n_local_blocks) * 5);
    for (const auto& m : local_meta) {
        local_meta_buf.push_back(m.block_id);
        local_meta_buf.push_back(m.ni_core);
        local_meta_buf.push_back(m.nj_core);
        local_meta_buf.push_back(m.nk_core);
        local_meta_buf.push_back(m.data_size);
    }

    std::vector<Int> recv_counts_meta(static_cast<std::size_t>(nprocs));
    std::vector<Int> recv_displs_meta(static_cast<std::size_t>(nprocs));
    Int total_meta_entries = 0;

    if (my_rank == 0) {
        for (int r = 0; r < nprocs; ++r) {
            recv_counts_meta[static_cast<std::size_t>(r)] =
                nblocks_per_rank[static_cast<std::size_t>(r)] * 5;
            recv_displs_meta[static_cast<std::size_t>(r)] = total_meta_entries;
            total_meta_entries += recv_counts_meta[static_cast<std::size_t>(r)];
        }
    }

    std::vector<Int> global_meta_buf(static_cast<std::size_t>(total_meta_entries));
    MPI_Gatherv(local_meta_buf.data(), n_local_blocks * 5, MPI_INT,
                global_meta_buf.data(), recv_counts_meta.data(),
                recv_displs_meta.data(), MPI_INT, 0, comm);

    // --- Gather data sizes ---
    std::vector<Int> data_sizes_per_rank(static_cast<std::size_t>(nprocs));
    MPI_Gather(&local_total_size, 1, MPI_INT,
               data_sizes_per_rank.data(), 1, MPI_INT, 0, comm);

    // --- Gather floating-point data ---
    std::vector<Int> recv_counts_data(static_cast<std::size_t>(nprocs));
    std::vector<Int> recv_displs_data(static_cast<std::size_t>(nprocs));
    Int total_data_entries = 0;

    if (my_rank == 0) {
        for (int r = 0; r < nprocs; ++r) {
            recv_counts_data[static_cast<std::size_t>(r)] =
                data_sizes_per_rank[static_cast<std::size_t>(r)];
            recv_displs_data[static_cast<std::size_t>(r)] = total_data_entries;
            total_data_entries += recv_counts_data[static_cast<std::size_t>(r)];
        }
    }

    std::vector<Real> global_data(static_cast<std::size_t>(total_data_entries));
    MPI_Gatherv(local_data.data(), local_total_size, MPI_DOUBLE,
                global_data.data(), recv_counts_data.data(),
                recv_displs_data.data(), MPI_DOUBLE, 0, comm);

    // ---- Step 3: rank 0 writes the Tecplot file ----
    if (my_rank != 0) return;

    // Reconstruct per-block metadata and offsets
    struct GlobalBlockInfo {
        Int block_id, ni_core, nj_core, nk_core, data_size;
        Int data_offset;  // offset into global_data (in Reals)
    };
    std::vector<GlobalBlockInfo> all_blocks;

    {
        const Int* meta_ptr = global_meta_buf.data();
        Int data_offset = 0;
        for (int r = 0; r < nprocs; ++r) {
            Int nb = nblocks_per_rank[static_cast<std::size_t>(r)];
            for (Int b = 0; b < nb; ++b) {
                GlobalBlockInfo info;
                info.block_id   = meta_ptr[0];
                info.ni_core    = meta_ptr[1];
                info.nj_core    = meta_ptr[2];
                info.nk_core    = meta_ptr[3];
                info.data_size  = meta_ptr[4];
                info.data_offset = data_offset;
                all_blocks.push_back(info);
                data_offset += info.data_size;
                meta_ptr += 5;
            }
        }
    }

    // Sort by block_id for deterministic output
    std::sort(all_blocks.begin(), all_blocks.end(),
              [](const GlobalBlockInfo& a, const GlobalBlockInfo& b) {
                  return a.block_id < b.block_id;
              });

    // Write
    std::ofstream os(filename);
    if (!os.is_open()) {
        throw std::runtime_error("TecplotWriter: cannot open \"" +
                                 filename + "\" for writing");
    }

    os << std::scientific << std::setprecision(6);
    os << "TITLE = \"WCNS v2 Solution — iter=" << iter
       << ", time=" << time << "\"\n";
    os << "VARIABLES = \"X\", \"Y\", \"Z\", "
       << "\"rho\", \"u\", \"v\", \"w\", \"p\", \"T\"\n";

    for (const auto& blk : all_blocks) {
        os << "ZONE T=\"block_" << blk.block_id
           << "\", I=" << blk.ni_core
           << ", J=" << blk.nj_core
           << ", K=" << blk.nk_core
           << ", DATAPACKING=POINT\n";

        const Real* ptr = global_data.data() + blk.data_offset;
        Int n_vals = blk.data_size;
        for (Int v = 0; v < n_vals; v += 9) {
            os << " "  << ptr[0] << " " << ptr[1] << " " << ptr[2]
               << " "  << ptr[3] << " " << ptr[4] << " " << ptr[5]
               << " "  << ptr[6] << " " << ptr[7] << " " << ptr[8] << "\n";
            ptr += 9;
        }
    }

    os.close();
}

// ============================================================================
// CgnsWriter::write (reserved stub)
// ============================================================================

inline void CgnsWriter::write(const std::vector<LocalBlock>& /*blocks*/,
                               const std::string& /*filename*/,
                               Int /*iter*/, Real /*time*/) {
    throw std::runtime_error("CGNS output not yet implemented");
}

// ============================================================================
// VtkWriter::write (reserved stub)
// ============================================================================

inline void VtkWriter::write(const std::vector<LocalBlock>& /*blocks*/,
                              const std::string& /*filename*/,
                              Int /*iter*/, Real /*time*/) {
    throw std::runtime_error("VTK output not yet implemented");
}

// ============================================================================
// SolutionWriter::write (top-level dispatch)
// ============================================================================

inline void SolutionWriter::write(const std::vector<LocalBlock>& blocks,
                                   const Config& cfg, Int iter, Real time) {
    auto filename = make_filename(cfg.output_dir, "sol", iter, ".plt");

    if (cfg.output_format == "tecplot") {
        TecplotWriter::write(blocks, cfg, filename, iter, time);
    } else if (cfg.output_format == "cgns") {
        CgnsWriter::write(blocks, filename, iter, time);
    } else if (cfg.output_format == "vtk") {
        VtkWriter::write(blocks, filename, iter, time);
    } else {
        throw std::runtime_error("SolutionWriter: unknown output_format \""
                                 + cfg.output_format + "\"");
    }
}
