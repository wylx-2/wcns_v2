#pragma once

#include "wcns_v2/utils/types.h"
#include "wcns_v2/parallel/local_block.h"
#include <mpi.h>
#include <map>
#include <string>
#include <vector>

/// @file flux_halo_exchange.h
/// @brief Face-flux halo exchange for inviscid fluxes at connectivity boundaries.
///
/// After the Riemann solver computes inviscid face fluxes (inv_xi, inv_eta,
/// inv_zeta), blocks that share a 1-to-1 connectivity (inter-zone or periodic)
/// must exchange face flux values so that the ghost-face fluxes match the
/// corresponding interior-face fluxes of the neighbor block.
///
/// This is necessary because the cell-centered halo exchange may not correctly
/// populate ghost cells at interface boundaries (particularly for same-process
/// multi-block), which degrades the WCNS interpolation at faces near the
/// interface.  By exchanging the face fluxes directly, we ensure that the
/// 6-point centered difference in InviscidRHS sees consistent flux values
/// across the interface.
///
/// Exchange count: ng+1 faces per connectivity boundary per direction.
///
/// Face index mapping (JMAX ↔ JMIN example, ncj cells, ng ghost layers):
///
///   Block A (JMAX, face=3):
///     Send: faces [ncj-2*ng, ncj-ng]  →  Block B receives at [0, ng]
///     Recv: faces [ncj-ng, ncj]       ←  Block B sends    [ng, 2*ng]
///
///   Block B (JMIN, face=2):
///     Send: faces [ng, 2*ng]          →  Block A receives at [ncj-ng, ncj]
///     Recv: faces [0, ng]             ←  Block A sends    [ncj-2*ng, ncj-ng]

class FluxHaloExchange {
public:
    FluxHaloExchange();

    /// Resolve inter-zone connectivity for a local block.
    ///
    /// Examines the block's grid 1-to-1 connections and the zone→block mapping
    /// to determine which faces require flux exchange.  For each such face,
    /// pre-computes buffer sizes, send/recv face ranges, and target info.
    ///
    /// @param block           The local block to set up
    /// @param zone_to_block   Map from original zone name → (rank, global_block_id)
    void setup(const LocalBlock& block,
               const std::map<std::string, std::pair<int,int>>& zone_to_block);

    /// Exchange inviscid face fluxes for all three directions.
    ///
    /// Reads inv_xi, inv_eta, inv_zeta from block.field and overwrites
    /// ghost-face entries with values received from neighbor blocks.
    ///
    /// @param block        The local block (field.flux arrays are read & written)
    /// @param all_blocks   All local blocks on this process (for same-process copy)
    void exchange(LocalBlock& block, const std::vector<LocalBlock>& all_blocks);

    /// Exchange an arbitrary set of single-component face arrays in one direction.
    ///
    /// All arrays in @p face_arrs must have the same dimensions (matching the
    /// face direction) and are packed into a single MPI buffer for efficiency.
    /// Only handles MPI-remote neighbors; same-process local copy must be done
    /// separately by the caller (the arrays are typically temporaries not
    /// accessible via Field/Grid).
    ///
    /// @param face_arrs    Vector of face arrays to exchange (all in the same dir)
    /// @param dir          Direction: 0=ξ, 1=η, 2=ζ
    /// @param block        The local block
    void exchange_face_arrays(std::vector<MultiArray3D<Real>*>& face_arrs,
                              int dir, const LocalBlock& block);

private:
    /// Per-face exchange descriptor for one connectivity boundary.
    struct FluxFaceInfo {
        bool active = false;       ///< true if this face has flux exchange
        bool is_remote = false;    ///< true if target is on a different MPI rank
        int  target_rank = -1;     ///< MPI rank of neighbor (-1 = same process)
        int  target_block = -1;    ///< Global block ID of neighbor

        /// Face index ranges for send and receive (in the varying dimension).
        Int send_begin, send_end;  ///< [send_begin, send_end] inclusive faces to send
        Int recv_begin, recv_end;  ///< [recv_begin, recv_end] inclusive faces to receive
        Int n_faces;               ///< number of faces = send_end - send_begin + 1

        /// Dimensions of one face slice (the two non-varying dimensions).
        Int dim1, dim2;            ///< size of one face in the two transverse directions

        /// Pre-allocated buffers (5 components × n_faces × dim1 × dim2 values).
        std::vector<Real> send_buf;
        std::vector<Real> recv_buf;
    };

    FluxFaceInfo faces_[6];  ///< one per block face

    /// Pack ng+1 face slices of a FluxVars array into a send buffer.
    /// @param fv      The 5-component flux array (inv_xi, inv_eta, or inv_zeta)
    /// @param dir     0=ξ (varying i), 1=η (varying j), 2=ζ (varying k)
    /// @param info    Pre-computed face exchange info (send_buf is written)
    void pack_flux(const FluxVars& fv, int dir, FluxFaceInfo& info);

    /// Unpack received face slices into a FluxVars array.
    void unpack_flux(FluxVars& fv, int dir, const FluxFaceInfo& info);

    /// Direct memory copy from a neighbor block on the same process.
    void copy_local_flux(LocalBlock& block, int face, int dir,
                         const FluxFaceInfo& info,
                         const LocalBlock& neighbor);
};

#include "flux_halo_exchange.hxx"
