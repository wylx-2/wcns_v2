#pragma once

#include "utils/types.h"
#include "core/config.h"
#include <string>
#include <vector>

/// @file solution_writer.h
/// @brief Solution output in Tecplot / CGNS / VTK format.
///
/// SolutionWriter dispatches to the format-specific writer based on
/// cfg.output_format.  Currently only Tecplot (ASCII .plt) is implemented;
/// CGNS and VTK are reserved and throw at runtime.
///
/// Output range: interior cells only [ng .. nci-1-ng], i.e. [3 .. nci-4].
/// Output variables: cell-center coordinates (x,y,z) + primitive vars (rho,u,v,w,p,T).

// ============================================================================
// Forward declarations
// ============================================================================

class LocalBlock;

// ============================================================================
// Helper: filename generation
// ============================================================================

/// Build a filename with zero-padded iteration number.
///   make_filename("output", "sol", 100, ".plt") → "output/sol_000100.plt"
inline std::string make_filename(const std::string& dir,
                                  const std::string& prefix,
                                  Int iter,
                                  const std::string& suffix);

// ============================================================================
// TecplotWriter — multi-zone ASCII Tecplot format
// ============================================================================

/// Writes a single .plt file with one ZONE per block.
///
/// Gathers interior primitive vars and cell-center coordinates from all
/// blocks (across MPI ranks) to rank 0, then writes sequentially.
class TecplotWriter {
public:
    /// Write Tecplot multi-zone solution file.
    ///
    /// @param blocks    All local blocks on this rank
    /// @param cfg       Configuration (for eos_factor to compute temperature)
    /// @param filename  Full path (e.g. "output/sol_000100.plt")
    /// @param iter      Current iteration (for title annotation)
    /// @param time      Current physical time (for title annotation)
    static void write(const std::vector<LocalBlock>& blocks,
                      const Config& cfg,
                      const std::string& filename,
                      Int iter, Real time);
};

// ============================================================================
// CgnsWriter — CGNS output (reserved)
// ============================================================================

class CgnsWriter {
public:
    /// Write CGNS solution file (NOT YET IMPLEMENTED).
    /// @throws std::runtime_error always
    static void write(const std::vector<LocalBlock>& blocks,
                      const std::string& filename,
                      Int iter, Real time);
};

// ============================================================================
// VtkWriter — VTK structured grid output (reserved)
// ============================================================================

class VtkWriter {
public:
    /// Write VTK solution file (NOT YET IMPLEMENTED).
    /// @throws std::runtime_error always
    static void write(const std::vector<LocalBlock>& blocks,
                      const std::string& filename,
                      Int iter, Real time);
};

// ============================================================================
// SolutionWriter — top-level dispatch
// ============================================================================

/// Top-level solution output interface.
///
/// Usage (main loop):
///   if (iter % cfg.output_freq == 0)
///       SolutionWriter::write(blocks, cfg, iter, time);
class SolutionWriter {
public:
    /// Write solution file.  Dispatches to format-specific writer.
    ///
    /// @param blocks   All local blocks on this rank
    /// @param cfg      Configuration (output_format, output_dir, output_freq)
    /// @param iter     Current iteration number
    /// @param time     Current physical time
    static void write(const std::vector<LocalBlock>& blocks,
                      const Config& cfg, Int iter, Real time);
};

#include "solution_writer.hxx"
