#pragma once

#include "history_monitor.h"
#include "wcns_v2/parallel/parallel_env.h"
#include "wcns_v2/parallel/parallel_manager.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <mpi.h>

// ============================================================================
// Helper: estimate local grid spacing in x-direction
// ============================================================================
namespace {

inline Real estimate_dx(const std::vector<LocalBlock>& blocks) {
    Real dx_min = std::numeric_limits<Real>::max();

    for (const auto& lb : blocks) {
        Int ng  = lb.grid.ng;
        Int i0  = ng;
        Int i1  = ng + lb.nci_core() - 1;
        Int j0  = ng;
        Int k0  = ng;

        if (i1 <= i0) continue;

        // Sample at a few j,k locations
        Int j_sample = j0;
        Int k_sample = k0;

        for (Int i = i0; i < i1; ++i) {
            Real dx_i = std::abs(lb.grid.cell_x(i+1, j_sample, k_sample)
                                - lb.grid.cell_x(i,   j_sample, k_sample));
            if (dx_i > 0) dx_min = std::min(dx_min, dx_i);
        }
    }

    return ParallelManager::global_min(dx_min);
}

}  // anonymous namespace

// ============================================================================
// HistoryMonitor::write_header
// ============================================================================

inline void HistoryMonitor::write_header(std::ostream& os,
                                          const std::vector<Real>& x_locations) {
    os << "# iter       time        dt         ";
    for (std::size_t n = 0; n < x_locations.size(); ++n) {
        os << "Uavg_x" << std::setw(2) << n << "(" << std::fixed
           << std::setprecision(3) << x_locations[n] << ")";
        if (n + 1 < x_locations.size()) os << "  ";
    }
    os << "\n";
}

// ============================================================================
// HistoryMonitor::compute_averages
// ============================================================================

inline std::vector<Real> HistoryMonitor::compute_averages(
        const std::vector<LocalBlock>& blocks,
        const std::vector<Real>& x_locations,
        Real tolerance) {

    std::size_t n_targets = x_locations.size();

    // Local accumulators: volume-weighted sum_u and total_vol per target
    std::vector<Real> local_sum_u(n_targets, 0.0);
    std::vector<Real> local_vol(n_targets, 0.0);

    // Auto-compute tolerance if not provided
    Real tol = tolerance;
    if (tol <= 0.0) {
        tol = 0.6 * estimate_dx(blocks);
    }

    for (const auto& lb : blocks) {
        Int ng  = lb.grid.ng;
        Int i0  = ng;
        Int i1  = ng + lb.nci_core() - 1;
        Int j0  = ng;
        Int j1  = ng + lb.ncj_core() - 1;
        Int k0  = ng;
        Int k1  = ng + lb.nck_core() - 1;

        for (Int k = k0; k <= k1; ++k) {
            for (Int j = j0; j <= j1; ++j) {
                for (Int i = i0; i <= i1; ++i) {
                    Real xc = lb.grid.cell_x(i,j,k);

                    for (std::size_t n = 0; n < n_targets; ++n) {
                        if (std::abs(xc - x_locations[n]) <= tol) {
                            Real vol = lb.grid.cell_vol(i,j,k);
                            local_sum_u[n] += lb.field.prim.u(i,j,k) * vol;
                            local_vol[n]  += vol;
                        }
                    }
                }
            }
        }
    }

    // MPI reduction
    std::vector<Real> global_sum_u(n_targets, 0.0);
    std::vector<Real> global_vol(n_targets, 0.0);

    if (ParallelEnv::is_parallel()) {
        MPI_Allreduce(local_sum_u.data(), global_sum_u.data(),
                      static_cast<int>(n_targets), MPI_DOUBLE, MPI_SUM,
                      ParallelEnv::communicator());
        MPI_Allreduce(local_vol.data(), global_vol.data(),
                      static_cast<int>(n_targets), MPI_DOUBLE, MPI_SUM,
                      ParallelEnv::communicator());
    } else {
        global_sum_u = local_sum_u;
        global_vol   = local_vol;
    }

    // Compute volume-weighted averages
    std::vector<Real> averages(n_targets);
    for (std::size_t n = 0; n < n_targets; ++n) {
        if (global_vol[n] > 0.0) {
            averages[n] = global_sum_u[n] / global_vol[n];
        } else {
            averages[n] = std::numeric_limits<Real>::quiet_NaN();
        }
    }

    return averages;
}

// ============================================================================
// HistoryMonitor::log
// ============================================================================

inline void HistoryMonitor::log(std::ostream& os, Int iter, Real time, Real dt,
                                 const std::vector<Real>& averages) {
    os << std::scientific << std::setprecision(8)
       << std::setw(8)  << iter  << " "
       << std::setw(15) << time  << " "
       << std::setw(15) << dt;

    for (const auto& avg : averages) {
        os << " " << std::setw(15) << avg;
    }
    os << std::endl;
}
