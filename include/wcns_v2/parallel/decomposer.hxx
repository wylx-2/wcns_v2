#pragma once

#include "decomposer.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

// ============================================================================
// Public interface
// ============================================================================

inline std::vector<SubBlock> BlockDecomposer::decompose(
        const std::vector<Grid>& zones, int nprocs, Int ng) {
    Int nz = static_cast<Int>(zones.size());
    Int min_cells = 2 * ng;  // minimum cells per direction after splitting

    // ---- Collect cell counts ----
    std::vector<Int> counts(nz);
    for (Int z = 0; z < nz; ++z) counts[z] = cell_count_core(zones[z]);

    // ---- Phase 1: Greedy assignment ----
    std::vector<int> proc_loads(nprocs, 0);
    std::vector<int> assignments(nz, 0);
    assign_greedy(counts, nprocs, proc_loads, assignments);

    // ---- Phase 2: If nprocs > nz, split the largest blocks ----
    std::vector<SubBlock> result;

    if (nprocs <= nz) {
        // No splitting needed — one sub-block per zone
        for (Int z = 0; z < nz; ++z) {
            SubBlock sb;
            sb.zone_id       = z;
            sb.sub_id        = 0;
            sb.assigned_rank = assignments[z];
            sb.ci_min = 0;  sb.ci_max = zones[z].ni_core - 2;  // nci_core = ni_core-1
            sb.cj_min = 0;  sb.cj_max = zones[z].nj_core - 2;
            sb.ck_min = 0;  sb.ck_max = zones[z].nk_core - 2;
            result.push_back(sb);
        }
    } else {
        // Need to split: assign one sub-block per zone first, then split rest
        int n_extra = nprocs - nz;  // how many more sub-blocks we need

        // Build list of (zone_index, cell_count) sorted by cell count descending
        std::vector<std::pair<Int, Int>> idx_by_size;  // (cell_count, zone_index)
        for (Int z = 0; z < nz; ++z)
            idx_by_size.emplace_back(counts[z], z);
        std::sort(idx_by_size.begin(), idx_by_size.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        // Determine splits per zone
        std::vector<int> zone_splits(nz, 1);  // how many sub-blocks per zone
        int remaining = n_extra;

        for (auto& [cnt, zid] : idx_by_size) {
            if (remaining <= 0) break;
            Int nc_i = zones[zid].ni_core - 1;
            Int nc_j = zones[zid].nj_core - 1;
            Int nc_k = zones[zid].nk_core - 1;

            // Max possible subdivisions given min_cell constraint
            int max_possible = static_cast<int>(nc_i / min_cells)
                             * static_cast<int>(nc_j / min_cells)
                             * static_cast<int>(nc_k / min_cells);
            max_possible = std::max(1, max_possible);

            // How many additional splits for this zone?
            int extra_for_zone = std::min(remaining, max_possible - 1);
            zone_splits[zid] = 1 + extra_for_zone;
            remaining -= extra_for_zone;
        }

        // Now do the actual splitting, round-robin across all processes
        int next_rank = 0;
        for (Int z = 0; z < nz; ++z) {
            split_zone(zones[z], z, zone_splits[z], nprocs, ng, result, next_rank);
        }
    }

    // ---- Print summary on rank 0 ----
    if (result.size() >= static_cast<std::size_t>(nz)) {
        std::cout << "Domain decomposition: " << nz << " zone(s) → "
                  << result.size() << " sub-block(s) on " << nprocs
                  << " process(es)\n";
        for (const auto& sb : result) {
            Int nc_i = sb.ci_max - sb.ci_min + 1;
            Int nc_j = sb.cj_max - sb.cj_min + 1;
            Int nc_k = sb.ck_max - sb.ck_min + 1;
            std::cout << "  zone=" << sb.zone_id << " sub=" << sb.sub_id
                      << " rank=" << sb.assigned_rank
                      << " cells=[" << nc_i << "x" << nc_j << "x" << nc_k << "]"
                      << " range_i=[" << sb.ci_min << ":" << sb.ci_max << "]"
                      << " range_j=[" << sb.cj_min << ":" << sb.cj_max << "]"
                      << " range_k=[" << sb.ck_min << ":" << sb.ck_max << "]\n";
        }
    }

    return result;
}

// ============================================================================
// Private
// ============================================================================

inline void BlockDecomposer::assign_greedy(
        const std::vector<Int>& cell_counts, int nprocs,
        std::vector<int>& proc_loads, std::vector<int>& assignments) {
    Int n = static_cast<Int>(cell_counts.size());

    // Index list sorted by cell count descending
    std::vector<Int> idx_by_size(n);
    std::iota(idx_by_size.begin(), idx_by_size.end(), 0);
    std::sort(idx_by_size.begin(), idx_by_size.end(),
              [&](Int a, Int b) { return cell_counts[a] > cell_counts[b]; });

    for (Int z_idx : idx_by_size) {
        // Find least-loaded process
        int best_proc = 0;
        for (int p = 1; p < nprocs; ++p)
            if (proc_loads[p] < proc_loads[best_proc]) best_proc = p;

        assignments[z_idx] = best_proc;
        proc_loads[best_proc] += static_cast<int>(cell_counts[z_idx]);
    }
}

inline void BlockDecomposer::find_split_factors(
        Int nc_i, Int nc_j, Int nc_k, int target, Int min_cells,
        int& si, int& sj, int& sk) {

    si = sj = sk = 1;

    if (target <= 1) return;

    // Maximum splits per direction
    int max_i = std::max(1, static_cast<int>(nc_i / min_cells));
    int max_j = std::max(1, static_cast<int>(nc_j / min_cells));
    int max_k = std::max(1, static_cast<int>(nc_k / min_cells));

    // Try all factorizations of target, pick the one with the most balanced
    // split ratios (closest to 1:1:1 in terms of resulting sub-block shape)
    int best = 1;
    Real best_score = 1e30;

    for (int i = 1; i <= max_i && i <= target; ++i) {
        if (target % i != 0) continue;
        int rem = target / i;
        for (int j = 1; j <= max_j && j <= rem; ++j) {
            if (rem % j != 0) continue;
            int k = rem / j;
            if (k < 1 || k > max_k) continue;
            if (i * j * k != target) continue;

            // Score: how close the sub-block shape is to a cube
            Real rx = static_cast<Real>(nc_i) / i;
            Real ry = static_cast<Real>(nc_j) / j;
            Real rz = static_cast<Real>(nc_k) / k;
            Real r_avg = (rx + ry + rz) / 3.0;
            Real score = std::abs(rx - r_avg) + std::abs(ry - r_avg) + std::abs(rz - r_avg);

            if (score < best_score) {
                best_score = score;
                si = i; sj = j; sk = k;
                best = i * j * k;
            }
        }
    }

    // If target factorization didn't work, try nearest values
    if (best < target) {
        // Try factors that give product >= target
        for (int i = 1; i <= max_i; ++i) {
            for (int j = 1; j <= max_j; ++j) {
                for (int k = 1; k <= max_k; ++k) {
                    int prod = i * j * k;
                    if (prod < target) continue;
                    Real rx = static_cast<Real>(nc_i) / i;
                    Real ry = static_cast<Real>(nc_j) / j;
                    Real rz = static_cast<Real>(nc_k) / k;
                    Real r_avg = (rx + ry + rz) / 3.0;
                    Real score = std::abs(rx - r_avg) + std::abs(ry - r_avg) + std::abs(rz - r_avg)
                               + 100.0 * (prod - target);  // penalty for over-splitting
                    if (score < best_score) {
                        best_score = score;
                        si = i; sj = j; sk = k;
                    }
                }
            }
        }
    }
}

inline Int BlockDecomposer::cell_count_core(const Grid& g) {
    return (g.ni_core - 1) * (g.nj_core - 1) * (g.nk_core - 1);
}

inline void BlockDecomposer::split_zone(
        const Grid& zone, int zone_id, int n_blocks, int nprocs,
        Int ng, std::vector<SubBlock>& result, int& next_rank) {

    Int nc_i = zone.ni_core - 1;
    Int nc_j = zone.nj_core - 1;
    Int nc_k = zone.nk_core - 1;

    Int min_cells = 2 * ng;

    int si, sj, sk;
    find_split_factors(nc_i, nc_j, nc_k, n_blocks, min_cells, si, sj, sk);

    // Compute chunk sizes as evenly as possible
    Int chunk_i = nc_i / si;
    Int rem_i   = nc_i % si;
    Int chunk_j = nc_j / sj;
    Int rem_j   = nc_j % sj;
    Int chunk_k = nc_k / sk;
    Int rem_k   = nc_k % sk;

    int sub_id = 0;
    Int ck_start = 0;
    for (int kk = 0; kk < sk; ++kk) {
        Int ck_size = chunk_k + (kk < rem_k ? 1 : 0);
        Int cj_start = 0;
        for (int jj = 0; jj < sj; ++jj) {
            Int cj_size = chunk_j + (jj < rem_j ? 1 : 0);
            Int ci_start = 0;
            for (int ii = 0; ii < si; ++ii) {
                Int ci_size = chunk_i + (ii < rem_i ? 1 : 0);

                SubBlock sb;
                sb.zone_id       = zone_id;
                sb.sub_id        = sub_id;
                sb.assigned_rank = (next_rank++) % nprocs;

                sb.ci_min = ci_start;
                sb.ci_max = ci_start + ci_size - 1;
                sb.cj_min = cj_start;
                sb.cj_max = cj_start + cj_size - 1;
                sb.ck_min = ck_start;
                sb.ck_max = ck_start + ck_size - 1;

                result.push_back(sb);
                ++sub_id;
                ci_start += ci_size;
            }
            cj_start += cj_size;
        }
        ck_start += ck_size;
    }
}
