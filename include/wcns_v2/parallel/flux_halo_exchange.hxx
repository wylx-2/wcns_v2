#pragma once

#include "flux_halo_exchange.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <mpi.h>

// ============================================================================
// Constructor
// ============================================================================

inline FluxHaloExchange::FluxHaloExchange() {
    for (int f = 0; f < 6; ++f) {
        faces_[f].active = false;
        faces_[f].is_remote = false;
    }
}

// ============================================================================
// Setup — resolve inter-zone connectivity for one block
// ============================================================================

inline void FluxHaloExchange::setup(
        const LocalBlock& block,
        const std::map<std::string, std::pair<int,int>>& zone_to_block) {

    Int nci = block.grid.nci;
    Int ncj = block.grid.ncj;
    Int nck = block.grid.nck;
    Int ng  = block.grid.ng;

    for (int face = 0; face < 6; ++face) {
        FluxFaceInfo& info = faces_[face];
        info.active = false;

        // Check for 1-to-1 connectivity on this face
        const Connectivity* conn = block.grid.find_face_connection(face);
        if (!conn) continue;

        // Look up the donor zone in the zone→block map
        auto it = zone_to_block.find(conn->donor_name);
        if (it == zone_to_block.end()) {
            // Donor zone not found — might be on another process that we
            // don't track?  Skip for now (should not happen in valid config).
            continue;
        }

        info.active   = true;
        info.target_rank  = it->second.first;
        info.target_block = it->second.second;

        // Determine if remote (same-process blocks have target_rank == rank,
        // but we check by comparing with the block's own neighbor logic).
        // For now: same-process if target_rank == ParallelEnv::rank().
        // We'll resolve this in exchange() using the all_blocks_ list.
        info.is_remote = true;  // default; will check at exchange time

        // Compute send/recv face index ranges based on face direction
        Int ndim;  // number of faces in the varying dimension
        switch (face) {
        case 0: case 1: // I-face → inv_xi
            ndim = nci;  // inv_xi has nci+1 faces (0..nci)
            info.dim1 = ncj;
            info.dim2 = nck;
            break;
        case 2: case 3: // J-face → inv_eta
            ndim = ncj;  // inv_eta has ncj+1 faces (0..ncj)
            info.dim1 = nci;
            info.dim2 = nck;
            break;
        case 4: case 5: // K-face → inv_zeta
            ndim = nck;  // inv_zeta has nck+1 faces (0..nck)
            info.dim1 = nci;
            info.dim2 = ncj;
            break;
        default:
            info.active = false;
            continue;
        }

        info.n_faces = ng + 1;

        // IMIN=0, IMAX=1, JMIN=2, JMAX=3, KMIN=4, KMAX=5
        // MAX faces send interior near end, recv ghost at end
        // MIN faces send interior near start, recv ghost at start
        if (face == 1 || face == 3 || face == 5) {
            // MAX face
            info.send_begin = ndim - 2 * ng;
            info.send_end   = ndim - ng;
            info.recv_begin = ndim - ng;
            info.recv_end   = ndim;
        } else {
            // MIN face
            info.send_begin = ng;
            info.send_end   = 2 * ng;
            info.recv_begin = 0;
            info.recv_end   = ng;
        }

        // Allocate buffers: 5 components × n_faces × dim1 × dim2
        Int single_comp_sz = info.n_faces * info.dim1 * info.dim2;
        info.send_buf.resize(5 * single_comp_sz);
        info.recv_buf.resize(5 * single_comp_sz);
    }
}

// ============================================================================
// Pack — copy interior face slices into send buffer
// ============================================================================

inline void FluxHaloExchange::pack_flux(const FluxVars& fv, int dir,
                                         FluxFaceInfo& info) {
    Int n_faces = info.n_faces;
    Int dim1    = info.dim1;
    Int dim2    = info.dim2;
    Int s0      = info.send_begin;

    // One component's data size: n_faces × dim1 × dim2
    Int comp_sz = n_faces * dim1 * dim2;

    // Pack 5 components sequentially
    // f1 (density flux)
    {
        Real* buf = info.send_buf.data();
        for (Int d = 0; d < n_faces; ++d) {
            Int face_idx = s0 + d;
            for (Int b = 0; b < dim2; ++b) {
            for (Int a = 0; a < dim1; ++a) {
                Int buf_idx = d * dim1 * dim2 + a + dim1 * b;
                Real val = 0.0;
                switch (dir) {
                case 0: val = fv.f1(face_idx, a, b); break;
                case 1: val = fv.f1(a, face_idx, b); break;
                case 2: val = fv.f1(a, b, face_idx); break;
                }
                buf[buf_idx] = val;
            }}
        }
    }
    // f2 (x-momentum flux)
    {
        Real* buf = info.send_buf.data() + comp_sz;
        for (Int d = 0; d < n_faces; ++d) {
            Int face_idx = s0 + d;
            for (Int b = 0; b < dim2; ++b) {
            for (Int a = 0; a < dim1; ++a) {
                Int buf_idx = d * dim1 * dim2 + a + dim1 * b;
                Real val = 0.0;
                switch (dir) {
                case 0: val = fv.f2(face_idx, a, b); break;
                case 1: val = fv.f2(a, face_idx, b); break;
                case 2: val = fv.f2(a, b, face_idx); break;
                }
                buf[buf_idx] = val;
            }}
        }
    }
    // f3 (y-momentum flux)
    {
        Real* buf = info.send_buf.data() + 2 * comp_sz;
        for (Int d = 0; d < n_faces; ++d) {
            Int face_idx = s0 + d;
            for (Int b = 0; b < dim2; ++b) {
            for (Int a = 0; a < dim1; ++a) {
                Int buf_idx = d * dim1 * dim2 + a + dim1 * b;
                Real val = 0.0;
                switch (dir) {
                case 0: val = fv.f3(face_idx, a, b); break;
                case 1: val = fv.f3(a, face_idx, b); break;
                case 2: val = fv.f3(a, b, face_idx); break;
                }
                buf[buf_idx] = val;
            }}
        }
    }
    // f4 (z-momentum flux)
    {
        Real* buf = info.send_buf.data() + 3 * comp_sz;
        for (Int d = 0; d < n_faces; ++d) {
            Int face_idx = s0 + d;
            for (Int b = 0; b < dim2; ++b) {
            for (Int a = 0; a < dim1; ++a) {
                Int buf_idx = d * dim1 * dim2 + a + dim1 * b;
                Real val = 0.0;
                switch (dir) {
                case 0: val = fv.f4(face_idx, a, b); break;
                case 1: val = fv.f4(a, face_idx, b); break;
                case 2: val = fv.f4(a, b, face_idx); break;
                }
                buf[buf_idx] = val;
            }}
        }
    }
    // f5 (energy flux)
    {
        Real* buf = info.send_buf.data() + 4 * comp_sz;
        for (Int d = 0; d < n_faces; ++d) {
            Int face_idx = s0 + d;
            for (Int b = 0; b < dim2; ++b) {
            for (Int a = 0; a < dim1; ++a) {
                Int buf_idx = d * dim1 * dim2 + a + dim1 * b;
                Real val = 0.0;
                switch (dir) {
                case 0: val = fv.f5(face_idx, a, b); break;
                case 1: val = fv.f5(a, face_idx, b); break;
                case 2: val = fv.f5(a, b, face_idx); break;
                }
                buf[buf_idx] = val;
            }}
        }
    }
}

// ============================================================================
// Unpack — copy received data into ghost face positions
// ============================================================================

inline void FluxHaloExchange::unpack_flux(FluxVars& fv, int dir,
                                           const FluxFaceInfo& info) {
    Int n_faces = info.n_faces;
    Int dim1    = info.dim1;
    Int dim2    = info.dim2;
    Int r0      = info.recv_begin;

    Int comp_sz = n_faces * dim1 * dim2;

    // f1
    {
        const Real* buf = info.recv_buf.data();
        for (Int d = 0; d < n_faces; ++d) {
            Int face_idx = r0 + d;
            for (Int b = 0; b < dim2; ++b) {
            for (Int a = 0; a < dim1; ++a) {
                Int buf_idx = d * dim1 * dim2 + a + dim1 * b;
                Real val = buf[buf_idx];
                switch (dir) {
                case 0: fv.f1(face_idx, a, b) = val; break;
                case 1: fv.f1(a, face_idx, b) = val; break;
                case 2: fv.f1(a, b, face_idx) = val; break;
                }
            }}
        }
    }
    // f2
    {
        const Real* buf = info.recv_buf.data() + comp_sz;
        for (Int d = 0; d < n_faces; ++d) {
            Int face_idx = r0 + d;
            for (Int b = 0; b < dim2; ++b) {
            for (Int a = 0; a < dim1; ++a) {
                Int buf_idx = d * dim1 * dim2 + a + dim1 * b;
                Real val = buf[buf_idx];
                switch (dir) {
                case 0: fv.f2(face_idx, a, b) = val; break;
                case 1: fv.f2(a, face_idx, b) = val; break;
                case 2: fv.f2(a, b, face_idx) = val; break;
                }
            }}
        }
    }
    // f3
    {
        const Real* buf = info.recv_buf.data() + 2 * comp_sz;
        for (Int d = 0; d < n_faces; ++d) {
            Int face_idx = r0 + d;
            for (Int b = 0; b < dim2; ++b) {
            for (Int a = 0; a < dim1; ++a) {
                Int buf_idx = d * dim1 * dim2 + a + dim1 * b;
                Real val = buf[buf_idx];
                switch (dir) {
                case 0: fv.f3(face_idx, a, b) = val; break;
                case 1: fv.f3(a, face_idx, b) = val; break;
                case 2: fv.f3(a, b, face_idx) = val; break;
                }
            }}
        }
    }
    // f4
    {
        const Real* buf = info.recv_buf.data() + 3 * comp_sz;
        for (Int d = 0; d < n_faces; ++d) {
            Int face_idx = r0 + d;
            for (Int b = 0; b < dim2; ++b) {
            for (Int a = 0; a < dim1; ++a) {
                Int buf_idx = d * dim1 * dim2 + a + dim1 * b;
                Real val = buf[buf_idx];
                switch (dir) {
                case 0: fv.f4(face_idx, a, b) = val; break;
                case 1: fv.f4(a, face_idx, b) = val; break;
                case 2: fv.f4(a, b, face_idx) = val; break;
                }
            }}
        }
    }
    // f5
    {
        const Real* buf = info.recv_buf.data() + 4 * comp_sz;
        for (Int d = 0; d < n_faces; ++d) {
            Int face_idx = r0 + d;
            for (Int b = 0; b < dim2; ++b) {
            for (Int a = 0; a < dim1; ++a) {
                Int buf_idx = d * dim1 * dim2 + a + dim1 * b;
                Real val = buf[buf_idx];
                switch (dir) {
                case 0: fv.f5(face_idx, a, b) = val; break;
                case 1: fv.f5(a, face_idx, b) = val; break;
                case 2: fv.f5(a, b, face_idx) = val; break;
                }
            }}
        }
    }
}

// ============================================================================
// Local copy (same process)
// ============================================================================

inline void FluxHaloExchange::copy_local_flux(
        LocalBlock& block, int face, int dir,
        const FluxFaceInfo& info,
        const LocalBlock& neighbor) {

    Int n_faces = info.n_faces;
    Int dim1    = info.dim1;
    Int dim2    = info.dim2;
    Int r0      = info.recv_begin;

    // Determine the neighbor's send range
    // If we're receiving at our ghost faces (recv_begin..recv_end),
    // the neighbor sends from its send range.
    // For a MIN face (0,2,4), neighbor is on MAX face → sends from [ndim_nbr-2*ng, ndim_nbr-ng]
    // For a MAX face (1,3,5), neighbor is on MIN face → sends from [ng, 2*ng]

    Int nbr_ndim;
    Int nbr_s0;  // neighbor send start

    switch (dir) {
    case 0: nbr_ndim = neighbor.field.ni(); break;
    case 1: nbr_ndim = neighbor.field.nj(); break;
    case 2: nbr_ndim = neighbor.field.nk(); break;
    default: return;
    }

    Int ng = block.grid.ng;
    if (face == 1 || face == 3 || face == 5) {
        // We're MAX → neighbor is MIN → neighbor sends from [ng, 2*ng]
        nbr_s0 = ng;
    } else {
        // We're MIN → neighbor is MAX → neighbor sends from [ndim-2*ng, ndim-ng]
        nbr_s0 = nbr_ndim - 2 * ng;
    }

    // Access neighbor's flux arrays
    const FluxVars* nbr_fv = nullptr;
    switch (dir) {
    case 0: nbr_fv = &neighbor.field.inv_xi; break;
    case 1: nbr_fv = &neighbor.field.inv_eta; break;
    case 2: nbr_fv = &neighbor.field.inv_zeta; break;
    }
    if (!nbr_fv) return;

    // Access our flux arrays
    FluxVars* my_fv = nullptr;
    switch (dir) {
    case 0: my_fv = &block.field.inv_xi; break;
    case 1: my_fv = &block.field.inv_eta; break;
    case 2: my_fv = &block.field.inv_zeta; break;
    }
    if (!my_fv) return;

    // Copy 5 components
    // We need a lambda to do the copy for each component accessor pattern
    auto copy_comp = [&](auto& my_access, auto& nbr_access) {
        for (Int d = 0; d < n_faces; ++d) {
            Int my_idx  = r0 + d;
            Int nbr_idx = nbr_s0 + d;
            for (Int b = 0; b < dim2; ++b) {
            for (Int a = 0; a < dim1; ++a) {
                Real val = 0.0;
                switch (dir) {
                case 0: val = nbr_access(nbr_idx, a, b); break;
                case 1: val = nbr_access(a, nbr_idx, b); break;
                case 2: val = nbr_access(a, b, nbr_idx); break;
                }
                switch (dir) {
                case 0: my_access(my_idx, a, b) = val; break;
                case 1: my_access(a, my_idx, b) = val; break;
                case 2: my_access(a, b, my_idx) = val; break;
                }
            }}
        }
    };

    copy_comp(my_fv->f1, nbr_fv->f1);
    copy_comp(my_fv->f2, nbr_fv->f2);
    copy_comp(my_fv->f3, nbr_fv->f3);
    copy_comp(my_fv->f4, nbr_fv->f4);
    copy_comp(my_fv->f5, nbr_fv->f5);
}

// ============================================================================
// Exchange — main entry point
// ============================================================================

inline void FluxHaloExchange::exchange(LocalBlock& block,
                                         const std::vector<LocalBlock>& all_blocks) {

    MPI_Comm comm = MPI_COMM_WORLD;

    for (int face = 0; face < 6; ++face) {
        FluxFaceInfo& info = faces_[face];
        if (!info.active) continue;

        // Determine direction from face index
        int dir = face / 2;  // 0=ξ, 1=η, 2=ζ

        // Determine if the neighbor is local (same process)
        bool is_local = false;
        const LocalBlock* neighbor_ptr = nullptr;
        for (const auto& nbr : all_blocks) {
            if (nbr.block_id == info.target_block) {
                is_local = true;
                neighbor_ptr = &nbr;
                break;
            }
        }
        info.is_remote = !is_local;

        if (is_local && neighbor_ptr) {
            // Same-process direct copy
            copy_local_flux(block, face, dir, info, *neighbor_ptr);
        } else if (info.is_remote) {
            // MPI exchange
            Int total_sz = static_cast<Int>(info.send_buf.size());

            // Pack interior face fluxes into send buffer
            const FluxVars* fv = nullptr;
            switch (dir) {
            case 0: fv = &block.field.inv_xi; break;
            case 1: fv = &block.field.inv_eta; break;
            case 2: fv = &block.field.inv_zeta; break;
            }
            if (fv) pack_flux(*fv, dir, info);

            // Post non-blocking send and receive
            int tag = std::min(block.block_id, info.target_block) * 100
                    + face * 10 + dir + 200;  // +200 offset vs cell-centered exchange

            MPI_Request sreq, rreq;
            MPI_Isend(info.send_buf.data(), total_sz, MPI_DOUBLE,
                      info.target_rank, tag, comm, &sreq);
            MPI_Irecv(info.recv_buf.data(), total_sz, MPI_DOUBLE,
                      info.target_rank, tag, comm, &rreq);

            MPI_Wait(&sreq, MPI_STATUS_IGNORE);
            MPI_Wait(&rreq, MPI_STATUS_IGNORE);

            // Unpack received data into ghost face positions
            FluxVars* fv_out = nullptr;
            switch (dir) {
            case 0: fv_out = &block.field.inv_xi; break;
            case 1: fv_out = &block.field.inv_eta; break;
            case 2: fv_out = &block.field.inv_zeta; break;
            }
            if (fv_out) unpack_flux(*fv_out, dir, info);
        }
    }
}
