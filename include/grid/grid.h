#pragma once

#include "utils/types.h"
#include "grid/boundary_condition.h"
#include "grid/connectivity.h"
#include "scheme/interp_diff.h"
#include <string>
#include <vector>

/// Structured grid block (cell-centered FVM / FD formulation).
///
/// Node coordinates are stored at vertices: ni x nj x nk.
/// Cell-centered quantities (coordinates, flow variables, fluxes, metrics)
/// are stored at cell centers: nci x ncj x nck, where nci = ni-1, etc.
///
/// The primary computational domain is the cell-centered space.
/// Boundary conditions are applied to cell faces, and ghost cells
/// extend the cell-centered array for stencil support.
///
/// After extend_ghost_layers() is called, ni/nj/nk become the TOTAL
/// dimensions (core + 2*ng ghost vertices) and ni_core/nj_core/nk_core
/// store the original physical-domain vertex counts.
class Grid {
public:
    // ---- Node-space dimensions (total, including ghosts after extension) ----
    Int ni, nj, nk;          // number of vertices in each direction
    Int node_count() const;  // total vertices = ni * nj * nk

    // ---- Cell-space dimensions (total, including ghosts after extension) ----
    Int nci, ncj, nck;       // number of cells = ni-1, nj-1, nk-1
    Int cell_count() const;  // total cells  = nci * ncj * nck

    // ---- Ghost layers ----
    Int ng;                   // number of ghost cell layers (set by solver)

    // ---- Core (physical) dimensions before ghost extension ----
    Int ni_core, nj_core, nk_core;

    // ---- Node coordinates (vertex-centered, ni*nj*nk) ----
    /// After extension, x(i,j,k) uses 0-based index into the full array.
    /// The physical core is at [ng .. ng+ni_core-1] in each direction.
    MultiArray3D<Real> node_x, node_y, node_z;

    // ---- Cell-center coordinates (cell-centered, nci*ncj*nck) ----
    /// Computed as the average of 8 surrounding vertices.
    MultiArray3D<Real> cell_x, cell_y, cell_z;

    // ---- Cell volume (placeholder for metric computation) ----
    MultiArray3D<Real> cell_vol;

    // ---- SCMM metric coefficients at cell centers (nci x ncj x nck) ----
    MultiArray3D<Real> met_xi_x, met_xi_y, met_xi_z;
    MultiArray3D<Real> met_eta_x, met_eta_y, met_eta_z;
    MultiArray3D<Real> met_zeta_x, met_zeta_y, met_zeta_z;

    // ---- Jacobian (J, cell volume) at cell centers (nci x ncj x nck) ----
    MultiArray3D<Real> jacobian;

    // ---- Face metric coefficients (primary normal-direction on each face) ----
    /// X-face (i+1/2): size (nci+1) x ncj x nck
    MultiArray3D<Real> face_xi_x, face_xi_y, face_xi_z;
    /// Y-face (j+1/2): size nci x (ncj+1) x nck
    MultiArray3D<Real> face_eta_x, face_eta_y, face_eta_z;
    /// Z-face (k+1/2): size nci x ncj x (nck+1)
    MultiArray3D<Real> face_zeta_x, face_zeta_y, face_zeta_z;

    // ---- Boundary conditions ----
    BoundaryCondition bc;

    // ---- 1-to-1 connectivity (block interfaces / periodic BCs) ----
    ConnectivityList connections;

    // ---- Metadata ----
    std::string name;   // zone name from CGNS
    std::string metrics_type = "auto";  ///< "auto" = SCMM from grid, "uniform" = analytical Cartesian

    Grid();

    /// Allocate node and cell arrays for given vertex dimensions (core only).
    void allocate(Int ni_, Int nj_, Int nk_, Int ng_);

    /// Extend node and cell arrays by ng ghost layers on each side.
    /// Ghost node coordinates are filled by periodic copy (if a 1-to-1
    /// connection exists for the face) or by linear extrapolation.
    /// After this call, ni/nj/nk are the total (extended) dimensions,
    /// and ni_core/nj_core/nk_core hold the original sizes.
    void extend_ghost_layers();

    /// Compute cell-center coordinates as average of 8 vertex corners.
    void compute_cell_centers();

    /// Compute cell volumes (uniform Cartesian approximation for now).
    void compute_cell_volumes();

    /// Compute SCMM metric coefficients and Jacobian at cell centers.
    ///
    /// Requires ghost layers to be already extended and cell-center
    /// coordinates to be valid.  Metric terms are stored in the
    /// met_xi_*, met_eta_*, met_zeta_* arrays; Jacobian in jacobian.
    ///
    /// If metrics_type == "uniform", uses analytical uniform Cartesian
    /// metrics directly (ignoring the grid geometry).
    void compute_metrics();

    /// Compute analytical metrics for a uniform Cartesian grid.
    /// dx, dy, dz are inferred from interior cell centers.
    void compute_metrics_uniform();

    /// Interpolate cell-center metrics to face centers (half-nodes).
    /// Populates face_xi_*, face_eta_*, face_zeta_* arrays.
    /// Must be called after compute_metrics().
    void compute_face_metrics();

    /// Extract metrics from a donor grid at a given cell offset.
    ///
    /// Copies cell-center metrics (met_xi_*, met_eta_*, met_zeta_*, jacobian,
    /// cell_vol) and face metrics (face_xi_*, face_eta_*, face_zeta_*) from a
    /// donor (typically the full zone before decomposition) into this grid.
    ///
    /// This avoids recomputing metrics on sub-blocks where ghost cells at
    /// internal split boundaries contain extrapolated coordinates rather than
    /// the correct neighbour data — a problem for curvilinear grids.
    ///
    /// The donor must have its metrics already computed.
    /// @param donor  Source grid with pre-computed metrics and face metrics
    /// @param ci0    Starting cell i-index in the donor for this block
    /// @param cj0    Starting cell j-index in the donor for this block
    /// @param ck0    Starting cell k-index in the donor for this block
    void extract_metrics_from(const Grid& donor, Int ci0, Int cj0, Int ck0);

    /// Find a 1-to-1 connection covering the given face.
    /// face: 0=IMIN, 1=IMAX, 2=JMIN, 3=JMAX, 4=KMIN, 5=KMAX
    /// Returns nullptr if no connection exists for this face.
    const Connectivity* find_face_connection(int face) const;

    /// Fill ghost node coordinates on one face by copying from a donor zone's
    /// interior nodes, applying the connection's translation.
    ///
    /// Handles both within-zone periodic connections (donor == *this) and
    /// inter-zone block interfaces (donor is a different zone) uniformly.
    ///
    /// Ghost coordinate = donor_interior_coordinate - translation
    ///
    /// The donor zone must already have its core data in place (post extend).
    void fill_ghost_face_from_donor(int face, const Grid& donor,
                                    const Connectivity& conn);

    /// Print grid summary to stdout.
    void print_summary() const;

private:
    /// Fill ghost node coordinates on all 6 sides by linear extrapolation.
    /// Connection-face ghosts are overwritten later by fill_ghost_face_from_donor
    /// (called from ParallelManager after all zones have been extended).
    void fill_ghost_nodes();

    /// Fill ghost nodes on one face by linear extrapolation.
    void fill_ghost_face_extrapolate(int face);

    /// Build the face_connected[6] mask from the connectivity list.
    /// A face is "connected" if it has any 1-to-1 connection (periodic or interface).
    void build_face_connected(bool face_connected[6]) const;

    /// Compute the 3 coordinate derivatives in one direction (x_d, y_d, z_d)
    /// from cell-center coordinates.  Updates the three output arrays.
    void compute_coord_deriv(int dir, Real dh, const bool fp[6],
                             MultiArray3D<Real>& dx_d, MultiArray3D<Real>& dy_d,
                             MultiArray3D<Real>& dz_d);
};

#include "grid.hxx"
