#pragma once

#include "halo_exchange.h"
#include <algorithm>
#include <cstring>
#include <iostream>

// ============================================================================
// Constructor
// ============================================================================

inline HaloExchange::HaloExchange() {
    for (int f = 0; f < 6; ++f) {
        bufs_[f].active = false;
        bufs_[f].is_remote = false;
    }
}

// ============================================================================
// Setup
// ============================================================================

inline void HaloExchange::setup(const LocalBlock& block) {
    for (int face = 0; face < 6; ++face) {
        FaceBuffer& fb = bufs_[face];
        const NeighborInfo& ni = block.neighbors[face];

        fb.active    = ni.active;
        fb.is_remote = ni.active && ni.target_rank >= 0 &&
                        ni.target_rank != -1 &&                 // -1 means same-block periodic
                        ni.target_rank != ParallelEnv::rank();  // same rank → local copy

        if (!fb.active) continue;

        fb.target_rank  = ni.target_rank;
        fb.target_block = ni.target_block;

        Int buf_sz = buffer_size(face, block);
        fb.send_buf.resize(buf_sz);
        fb.recv_buf.resize(buf_sz);

        // Set sub-volume dimensions for pack/unpack
        Int nci = block.grid.nci;
        Int ncj = block.grid.ncj;
        Int nck = block.grid.nck;
        Int ng  = block.grid.ng;

        switch (face) {
        case 0: case 1: // I-face
            fb.ni_send = ng; fb.nj_send = ncj; fb.nk_send = nck; break;
        case 2: case 3: // J-face
            fb.ni_send = nci; fb.nj_send = ng; fb.nk_send = nck; break;
        case 4: case 5: // K-face
            fb.ni_send = nci; fb.nj_send = ncj; fb.nk_send = ng; break;
        }
    }
}

// ============================================================================
// Buffer size
// ============================================================================

inline Int HaloExchange::buffer_size(int face, const LocalBlock& block) const {
    Int nci = block.grid.nci;
    Int ncj = block.grid.ncj;
    Int nck = block.grid.nck;
    Int ng  = block.grid.ng;

    switch (face) {
    case 0: case 1: return ng  * ncj * nck;  // I-face
    case 2: case 3: return nci * ng  * nck;  // J-face
    case 4: case 5: return nci * ncj * ng;   // K-face
    default: return 0;
    }
}

// ============================================================================
// Pack — copy interior data into send buffer
// ============================================================================

inline void HaloExchange::pack_face(const MultiArray3D<Real>& arr, int face,
                                     const LocalBlock& block, Real* buf, Int offset) {
    Int nci = block.grid.nci;
    Int ncj = block.grid.ncj;
    Int nck = block.grid.nck;
    Int ng  = block.grid.ng;

    switch (face) {
    case 0: { // IMIN: send interior ng..2*ng-1
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int d = 0; d < ng; ++d)
            buf[offset + d + ng * (j + ncj * k)] = arr(ng + d, j, k);
        break;
    }
    case 1: { // IMAX: send interior nci-2*ng .. nci-ng-1
        Int i0 = nci - 2 * ng;
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int d = 0; d < ng; ++d)
            buf[offset + d + ng * (j + ncj * k)] = arr(i0 + d, j, k);
        break;
    }
    case 2: { // JMIN
        for (Int k = 0; k < nck; ++k)
        for (Int d = 0; d < ng; ++d)
        for (Int i = 0; i < nci; ++i)
            buf[offset + i + nci * (d + ng * k)] = arr(i, ng + d, k);
        break;
    }
    case 3: { // JMAX
        Int j0 = ncj - 2 * ng;
        for (Int k = 0; k < nck; ++k)
        for (Int d = 0; d < ng; ++d)
        for (Int i = 0; i < nci; ++i)
            buf[offset + i + nci * (d + ng * k)] = arr(i, j0 + d, k);
        break;
    }
    case 4: { // KMIN
        for (Int d = 0; d < ng; ++d)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i)
            buf[offset + i + nci * (j + ncj * d)] = arr(i, j, ng + d);
        break;
    }
    case 5: { // KMAX
        Int k0 = nck - 2 * ng;
        for (Int d = 0; d < ng; ++d)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i)
            buf[offset + i + nci * (j + ncj * d)] = arr(i, j, k0 + d);
        break;
    }
    }
}

// ============================================================================
// Unpack — copy from recv buffer into ghost cells
// ============================================================================

inline void HaloExchange::unpack_face(MultiArray3D<Real>& arr, int face,
                                       const LocalBlock& block,
                                       const Real* buf, Int offset) {
    Int nci = block.grid.nci;
    Int ncj = block.grid.ncj;
    Int nck = block.grid.nck;
    Int ng  = block.grid.ng;

    switch (face) {
    case 0: { // IMIN: ghost 0..ng-1
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int d = 0; d < ng; ++d)
            arr(d, j, k) = buf[offset + d + ng * (j + ncj * k)];
        break;
    }
    case 1: { // IMAX: ghost nci-ng .. nci-1
        Int i0 = nci - ng;
        for (Int k = 0; k < nck; ++k)
        for (Int j = 0; j < ncj; ++j)
        for (Int d = 0; d < ng; ++d)
            arr(i0 + d, j, k) = buf[offset + d + ng * (j + ncj * k)];
        break;
    }
    case 2: { // JMIN: ghost 0..ng-1
        for (Int k = 0; k < nck; ++k)
        for (Int d = 0; d < ng; ++d)
        for (Int i = 0; i < nci; ++i)
            arr(i, d, k) = buf[offset + i + nci * (d + ng * k)];
        break;
    }
    case 3: { // JMAX: ghost ncj-ng .. ncj-1
        Int j0 = ncj - ng;
        for (Int k = 0; k < nck; ++k)
        for (Int d = 0; d < ng; ++d)
        for (Int i = 0; i < nci; ++i)
            arr(i, j0 + d, k) = buf[offset + i + nci * (d + ng * k)];
        break;
    }
    case 4: { // KMIN: ghost 0..ng-1
        for (Int d = 0; d < ng; ++d)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i)
            arr(i, j, d) = buf[offset + i + nci * (j + ncj * d)];
        break;
    }
    case 5: { // KMAX: ghost nck-ng .. nck-1
        Int k0 = nck - ng;
        for (Int d = 0; d < ng; ++d)
        for (Int j = 0; j < ncj; ++j)
        for (Int i = 0; i < nci; ++i)
            arr(i, j, k0 + d) = buf[offset + i + nci * (j + ncj * d)];
        break;
    }
    }
}

// ============================================================================
// Local copy (same process, no MPI)
// ============================================================================

inline void HaloExchange::copy_local(MultiArray3D<Real>& arr, int face,
                                      const LocalBlock& block,
                                      const LocalBlock& neighbor) {
    // Direct (same-process) copy from neighbor's interior to block's ghost.
    //
    // Example: Block A (IMAX, face=1) ↔ Block B (IMIN, face=0)
    //   pack_face at IMIN from Block B → Block B's interior [ng, 2*ng-1]
    //   unpack_face at IMAX to Block A → Block A's ghost [nci-ng, nci-1]
    //
    // The neighbor's connecting face is ni.target_face (the opposite face
    // of the current block's face).  Data is read from neighbor's interior
    // cells adjacent to that face and written to the current block's ghost.

    const NeighborInfo& ni = block.neighbors[face];
    Int single_sz = buffer_size(face, block);

    // Use a temporary buffer for the pack→unpack sequence
    std::vector<Real> tmp_buf(static_cast<std::size_t>(single_sz));

    // Pack from neighbor's interior at the neighbor's connecting face
    pack_face(arr, ni.target_face, neighbor, tmp_buf.data(), 0);

    // Unpack to the current block's ghost at the current face
    unpack_face(arr, face, block, tmp_buf.data(), 0);
}

// ============================================================================
// Single-array exchange (blocking)
// ============================================================================

inline void HaloExchange::exchange(MultiArray3D<Real>& arr,
                                    const LocalBlock& block,
                                    const std::vector<LocalBlock>& all_blocks) {
    std::vector<MultiArray3D<Real>*> arrays = {&arr};
    exchange_multi(arrays, block, all_blocks);
}

// ============================================================================
// Multi-array exchange (blocking)
// ============================================================================

inline void HaloExchange::exchange_multi(
        const std::vector<MultiArray3D<Real>*>& arrays,
        const LocalBlock& block,
        const std::vector<LocalBlock>& all_blocks) {

    Int n_arrays = static_cast<Int>(arrays.size());
    MPI_Comm comm = MPI_COMM_WORLD;

    // We use a 3-phase approach to guarantee deadlock-free communication
    // even when periodic faces connect different face indices across ranks
    // (e.g. rank 0 IMIN ↔ rank 1 IMAX).
    //
    // Phase 1: Post all Irecvs for all remote faces
    // Phase 2: Pack + Isend for all remote faces  (same-process copies inline)
    // Phase 3: Waitall + Unpack

    MPI_Request all_reqs[12];  // max 6 sends + 6 recvs
    int n_reqs = 0;
    int face_recv_idx[6] = {-1, -1, -1, -1, -1, -1};
    int face_send_idx[6] = {-1, -1, -1, -1, -1, -1};

    // ---- Phase 1: Post all remote Irecvs ----
    for (int face = 0; face < 6; ++face) {
        FaceBuffer& fb = bufs_[face];
        if (!fb.active || !fb.is_remote) continue;

        const NeighborInfo& ni = block.neighbors[face];
        Int single_sz = buffer_size(face, block);
        Int total_sz  = single_sz * n_arrays;

        if (static_cast<Int>(fb.recv_buf.size()) < total_sz) {
            fb.recv_buf.assign(total_sz, 0.0);
        }

        int tag = std::min(block.block_id, ni.target_block) * 1000
                + std::min(face, ni.target_face) * 10
                + (ni.is_periodic ? 5 : 0)
                + 1;

        face_recv_idx[face] = n_reqs;
        MPI_Irecv(fb.recv_buf.data(), static_cast<int>(total_sz), MPI_DOUBLE,
                  ni.target_rank, tag, comm, &all_reqs[n_reqs++]);
    }

    // ---- Phase 2: Same-process copies + pack & Isend for remote faces ----
    for (int face = 0; face < 6; ++face) {
        FaceBuffer& fb = bufs_[face];
        const NeighborInfo& ni = block.neighbors[face];
        if (!fb.active) continue;

        Int single_sz = buffer_size(face, block);
        Int total_sz  = single_sz * n_arrays;

        if (!fb.is_remote) {
            // Same-process copy.  Two cases:
            //   1) Same-block periodic (target_block == -1 or self):
            //      pack from self at target_face, unpack to self at face.
            //      This copies interior data from the opposite periodic side.
            //   2) Same-rank different block (internal split boundary):
            //      pack from neighbor block's interior at target_face,
            //      unpack to current block's ghost at face.
            if (ni.target_face >= 0 && ni.target_face < 6) {
                if (ni.target_block == -1 ||
                    ni.target_block == block.block_id) {
                    // Case 1: Same-block periodic — copy from self
                    Int needed_sz = single_sz * n_arrays;
                    if (static_cast<Int>(fb.send_buf.size()) < needed_sz) {
                        fb.send_buf.resize(static_cast<std::size_t>(needed_sz));
                    }
                    for (Int a = 0; a < n_arrays; ++a) {
                        pack_face(*arrays[static_cast<std::size_t>(a)], ni.target_face,
                                  block, fb.send_buf.data(), a * single_sz);
                        unpack_face(*arrays[static_cast<std::size_t>(a)], face,
                                    block, fb.send_buf.data(), a * single_sz);
                    }
                } else {
                    // Case 2: Same-rank different block — copy from neighbor
                    const LocalBlock* neighbor = nullptr;
                    for (const auto& nbr : all_blocks) {
                        if (nbr.block_id == ni.target_block) {
                            neighbor = &nbr;
                            break;
                        }
                    }
                    if (neighbor) {
                        for (Int a = 0; a < n_arrays; ++a) {
                            copy_local(*arrays[static_cast<std::size_t>(a)],
                                       face, block, *neighbor);
                        }
                    }
                }
            }
            continue;
        }

        // Remote face: pack and Isend
        if (static_cast<Int>(fb.send_buf.size()) < total_sz) {
            fb.send_buf.assign(total_sz, 0.0);
        }

        for (Int a = 0; a < n_arrays; ++a) {
            pack_face(*arrays[static_cast<std::size_t>(a)], face, block,
                      fb.send_buf.data(), a * single_sz);
        }

        int tag = std::min(block.block_id, ni.target_block) * 1000
                + std::min(face, ni.target_face) * 10
                + (ni.is_periodic ? 5 : 0)
                + 1;

        face_send_idx[face] = n_reqs;
        MPI_Isend(fb.send_buf.data(), static_cast<int>(total_sz), MPI_DOUBLE,
                  ni.target_rank, tag, comm, &all_reqs[n_reqs++]);
    }

    // ---- Phase 3: Wait for all communications + unpack ----
    if (n_reqs > 0) {
        MPI_Waitall(n_reqs, all_reqs, MPI_STATUSES_IGNORE);
    }

    for (int face = 0; face < 6; ++face) {
        FaceBuffer& fb = bufs_[face];
        if (!fb.active || !fb.is_remote) continue;

        Int single_sz = buffer_size(face, block);

        for (Int a = 0; a < n_arrays; ++a) {
            unpack_face(*arrays[static_cast<std::size_t>(a)], face, block,
                        fb.recv_buf.data(), a * single_sz);
        }
    }
}

// ============================================================================
// Non-blocking exchange: start
// ============================================================================

inline void HaloExchange::start_exchange(
        const std::vector<MultiArray3D<Real>*>& arrays,
        const LocalBlock& block,
        const std::vector<LocalBlock>& /*all_blocks*/) {

    Int n_arrays = static_cast<Int>(arrays.size());
    n_arrays_packed_ = n_arrays;
    MPI_Comm comm = MPI_COMM_WORLD;

    for (int face = 0; face < 6; ++face) {
        FaceBuffer& fb = bufs_[face];
        if (!fb.active) continue;

        Int single_sz = buffer_size(face, block);
        Int total_sz  = single_sz * n_arrays;

        fb.send_buf.assign(total_sz, 0.0);
        fb.recv_buf.assign(total_sz, 0.0);

        if (!fb.is_remote) {
            // Same-process copy — handled in start_exchange for consistency.
            // Same-block periodic: pack from self at target_face.
            // Same-rank multi-block: cannot do direct copy here (no access
            // to neighbor's arrays); skip and rely on blocking exchange_multi
            // for correct behavior.
            // For now, same-rank multi-block is handled by exchange_multi only.
            const NeighborInfo& ni = block.neighbors[face];
            if (ni.target_block == -1 || ni.target_block == block.block_id) {
                // Same-block periodic: do the pack→unpack inline
                for (Int a = 0; a < n_arrays; ++a) {
                    pack_face(*arrays[static_cast<std::size_t>(a)], ni.target_face,
                              block, fb.send_buf.data(), a * single_sz);
                    unpack_face(*arrays[static_cast<std::size_t>(a)], face,
                                block, fb.send_buf.data(), a * single_sz);
                }
            }
            continue;
        }

        // Pack all arrays into send buffer
        for (Int a = 0; a < n_arrays; ++a) {
            pack_face(*arrays[static_cast<std::size_t>(a)], face, block,
                      fb.send_buf.data(), a * single_sz);
        }

        // Start non-blocking send and receive
        // Tag must be symmetric so sender and receiver agree.
        int tag = std::min(block.block_id, fb.target_block) * 1000
                + std::min(face, block.neighbors[face].target_face) * 10
                + (block.neighbors[face].is_periodic ? 5 : 0)
                + 1;

        MPI_Isend(fb.send_buf.data(), static_cast<int>(total_sz), MPI_DOUBLE,
                  fb.target_rank, tag, comm, &fb.send_req);

        MPI_Irecv(fb.recv_buf.data(), static_cast<int>(total_sz), MPI_DOUBLE,
                  fb.target_rank, tag, comm, &fb.recv_req);
    }
}

// ============================================================================
// Non-blocking exchange: wait
// ============================================================================

inline void HaloExchange::wait_exchange(const LocalBlock& /*block*/) {
    for (int face = 0; face < 6; ++face) {
        FaceBuffer& fb = bufs_[face];
        if (!fb.active || !fb.is_remote) continue;

        MPI_Wait(&fb.send_req, MPI_STATUS_IGNORE);
        MPI_Wait(&fb.recv_req, MPI_STATUS_IGNORE);
    }

    // NOTE: Unpack requires knowing which arrays were packed. Since start_exchange
    // only has the array list in scope, the caller should unpack explicitly,
    // or the arrays should be stored. For now, use exchange_multi (blocking) which
    // does pack → send → recv → unpack in one shot.
}
