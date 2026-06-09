#pragma once

#include "utils/types.h"
#include "parallel/local_block.h"
#include <mpi.h>
#include <vector>

/// @file halo_exchange.h
/// @brief Halo (ghost-cell) data exchange across MPI processes.
///
/// Each HaloExchange instance is bound to a specific LocalBlock.  It pre-allocates
/// send/receive buffers and handles packing/unpacking of 3D sub-volumes for the
/// block's 6 faces.
///
/// Supports both blocking and non-blocking (start/wait) exchange modes.
/// Multiple arrays can be exchanged together by packing them consecutively
/// into a single MPI message per face.

class HaloExchange {
public:
    HaloExchange();

    /// Initialize buffers for a specific LocalBlock.
    /// Must be called once before any exchange.
    void setup(const LocalBlock& block);

    /// Exchange ghost data for a single 3D array (blocking).
    void exchange(MultiArray3D<Real>& arr, const LocalBlock& block);

    /// Exchange multiple 3D arrays in a single MPI transaction (blocking).
    /// All arrays must have the same dimensions as the block's grid (nci×ncj×nck).
    void exchange_multi(const std::vector<MultiArray3D<Real>*>& arrays,
                        const LocalBlock& block);

    /// Start non-blocking exchange for multiple arrays.
    /// Returns immediately; communication proceeds in the background.
    /// Call wait_exchange() before accessing ghost cells.
    void start_exchange(const std::vector<MultiArray3D<Real>*>& arrays,
                        const LocalBlock& block);

    /// Wait for the non-blocking exchange started by start_exchange() to complete.
    void wait_exchange(const LocalBlock& block);

private:
    /// Per-face buffer and communication state.
    struct FaceBuffer {
        std::vector<Real> send_buf;
        std::vector<Real> recv_buf;
        int  target_rank;
        int  target_block;  // local index in all_blocks (for direct copy, -1 if remote)
        bool is_remote;
        bool active;
        MPI_Request send_req;
        MPI_Request recv_req;

        Int ni_send, nj_send, nk_send;  // dimensions of the send/recv sub-volume
    };

    FaceBuffer bufs_[6];

    /// Compute the buffer size for one face in number of Real values.
    Int buffer_size(int face, const LocalBlock& block) const;

    /// Pack interior data of one array into the send buffer for the given face.
    void pack_face(const MultiArray3D<Real>& arr, int face,
                   const LocalBlock& block, Real* buf, Int offset);

    /// Unpack data from the receive buffer into ghost cells of one array.
    void unpack_face(MultiArray3D<Real>& arr, int face,
                     const LocalBlock& block, const Real* buf, Int offset);

    /// Perform a direct (same-process) copy from neighbor block.
    void copy_local(MultiArray3D<Real>& arr, int face,
                    const LocalBlock& block, const LocalBlock& neighbor);

    // Store pointers to all local blocks for same-process halo copy
    const std::vector<LocalBlock>* all_blocks_ = nullptr;

    // Number of arrays packed (set by start_exchange, used by wait_exchange)
    Int n_arrays_packed_ = 0;
};

#include "halo_exchange.hxx"
