#pragma once

#include "parallel/local_block.h"
#include "utils/types.h"
#include <iosfwd>
#include <vector>

/// @file history_monitor.h
/// @brief Cross-section average velocity monitoring for channel/pipe flows.
///
/// Computes area-weighted (or arithmetic) average of u-velocity at specified
/// x-coordinate cross-sections.  Works with MPI — each rank contributes its
/// local cells, then results are reduced globally.

class HistoryMonitor {
public:
    /// Write CSV header to a history file.
    /// Columns: iter, time, dt, Uavg_x0, Uavg_x1, ...
    static void write_header(std::ostream& os,
                             const std::vector<Real>& x_locations);

    /// Compute cross-section-average u-velocity at the given x-locations.
    ///
    /// For each target x, cells whose cell-center x-coordinate lies within
    /// `tolerance` of the target are included.  If tolerance <= 0, it is
    /// auto-computed from the local grid spacing.
    ///
    /// Returns one average per target x-location (same order).
    /// If no cell matches a target, the corresponding entry is NaN.
    static std::vector<Real> compute_averages(
        const std::vector<LocalBlock>& blocks,
        const std::vector<Real>& x_locations,
        Real tolerance = 0.0);

    /// Write one data line to a history file.
    static void log(std::ostream& os, Int iter, Real time, Real dt,
                    const std::vector<Real>& averages);
};

#include "history_monitor.hxx"
