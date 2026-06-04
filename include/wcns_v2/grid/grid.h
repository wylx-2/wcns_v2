#pragma once

#include "wcns_v2/utils/types.h"
#include "wcns_v2/grid/boundary_condition.h"
#include "wcns_v2/grid/connectivity.h"
#include "wcns_v2/scheme/interp_diff.h"
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
    void compute_metrics();

    /// Interpolate cell-center metrics to face centers (half-nodes).
    /// Populates face_xi_*, face_eta_*, face_zeta_* arrays.
    /// Must be called after compute_metrics().
    void compute_face_metrics();

    /// Find a periodic 1-to-1 connection covering the given face.
    /// face: 0=IMIN, 1=IMAX, 2=JMIN, 3=JMAX, 4=KMIN, 5=KMAX
    /// Returns nullptr if no periodic connection exists for this face.
    const Connectivity* find_periodic_connection(int face) const;

    /// Find any 1-to-1 connection covering the given face (periodic or interface).
    const Connectivity* find_face_connection(int face) const;

    /// Fix ghost node coordinates on one face by copying from a donor zone.
    /// Used for inter-zone block interfaces (non-periodic 1-to-1 connections).
    /// The donor zone must already have its ghost layers extended.
    void fix_interface_ghost(int face, const Grid& donor,
                             const Connectivity& conn);

    /// Print grid summary to stdout.
    void print_summary() const;

private:
    /// Fill ghost node coordinates on all 6 sides.
    void fill_ghost_nodes();

    /// Fill ghost nodes on one face by linear extrapolation.
    void fill_ghost_face_extrapolate(int face);

    /// Fill ghost nodes on one face by copying from a periodic donor.
    void fill_ghost_face_periodic(int face, const Connectivity& conn);

    /// Build the face_periodic[6] mask from the connectivity list.
    void build_face_periodic(bool face_periodic[6]) const;

    /// Compute the 3 coordinate derivatives in one direction (x_d, y_d, z_d)
    /// from cell-center coordinates.  Updates the three output arrays.
    void compute_coord_deriv(int dir, Real dh, const bool fp[6],
                             MultiArray3D<Real>& dx_d, MultiArray3D<Real>& dy_d,
                             MultiArray3D<Real>& dz_d);
};

#include "grid.hxx"
