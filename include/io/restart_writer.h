#pragma once

#include "utils/types.h"
#include "core/config.h"
#include <string>
#include <utility>
#include <vector>

/// @file restart_writer.h
/// @brief Binary restart I/O for checkpoint / resume.
///
/// RestartWriter saves the complete conservative variable state (including
/// ghost cells) so that the simulation can be resumed exactly.  Each MPI rank
/// writes its own data file; rank 0 additionally writes an index file that
/// maps global_block_id → rank for redistribution on a different process count.
///
/// File layout
/// -----------
/// Index file (restart_XXXXXX.bin, rank 0 only):
///   [uint32] magic      = 0x57434E53 ("WCNS")
///   [uint32] version    = 1
///   [Int]    iter
///   [Real]   time
///   [Int]    nprocs
///   [Int]    n_blocks_total
///   for each block: [Int] block_id, [Int] rank
///
/// Per-rank data file (restart_XXXXXX_rN.bin):
///   [Int]    n_blocks_local
///   for each local block:
///     [Int×6] block_id, zone_id, nci, ncj, nck, ng
///     [Real × nci*ncj*nck] cons.rho, cons.rhou, cons.rhov, cons.rhow, cons.rhoE
///     (each array stored as flat binary in (k,j,i)-major order)

// ============================================================================
// Forward declarations
// ============================================================================

class LocalBlock;

// ============================================================================
// RestartWriter
// ============================================================================

class RestartWriter {
public:
    /// Write restart files (binary).
    ///
    /// Each rank writes its own data to restart_XXXXXX_rN.bin.
    /// Rank 0 additionally writes the index file restart_XXXXXX.bin.
    ///
    /// @param blocks   All local blocks on this rank
    /// @param cfg      Configuration (output_dir)
    /// @param iter     Current iteration number
    /// @param time     Current physical time
    static void write(const std::vector<LocalBlock>& blocks,
                      const Config& cfg, Int iter, Real time);

    /// Read restart files and populate fields.
    ///
    /// Each rank reads its own restart_XXXXXX_rN.bin, loads conservative
    /// variables, then computes primitive variables via
    /// Field::cons_to_prim(cfg.gamma).
    ///
    /// @param basename  e.g. "restart_000500" (without rank suffix)
    /// @param blocks    [in/out] Local blocks — fields will be overwritten
    /// @param cfg       Configuration (provides gamma for cons→prim conversion)
    /// @return          {iteration, time} from the restart point
    static std::pair<Int, Real> read(const std::string& basename,
                                      std::vector<LocalBlock>& blocks,
                                      const Config& cfg);

private:
    static constexpr std::uint32_t kMagic   = 0x57434E53u;  // "WCNS"
    static constexpr std::uint32_t kVersion = 1;
};

#include "restart_writer.hxx"
