#pragma once

#include "wcns_v2/utils/types.h"

/// @file interp_diff.h
/// High-order interpolation and differentiation operators for WCNS.
///
/// All derivatives follow the Interpolation → Differentiation workflow:
///   1. Interpolate cell-centered values to half-nodes (faces)
///   2. Differentiate the half-node values back to cell centers
///
/// Interior: 6th-order centered stencils (7-point effective width).
/// Non-periodic boundaries: 5th-order one-sided stencils.
/// Periodic boundaries: centered stencils extend into ghost region
///   (ghost data is assumed valid from periodic copy / halo exchange).

// Forward declarations for friend access to one-sided interpolation stencils
class WcnsInterpBase;

class InterpDiff {
    friend class WcnsInterpBase;
public:
    // ========================================================================
    // 3D array operations
    // ========================================================================

    /// Compute the derivative of a cell-centered 3D array in a given direction.
    ///
    /// Internally allocates a half-node array, interpolates cell centers to
    /// half-nodes along @p dir , then differentiates back to cell centers.
    ///
    /// @param a     Input cell-centered array (ni x nj x nk)
    /// @param da    Output derivative array (ni x nj x nk), must be pre-allocated
    /// @param dir   Direction: 0 = ξ(i), 1 = η(j), 2 = ζ(k)
    /// @param dh    Grid spacing in computational space (Δξ, Δη, or Δζ)
    /// @param ng    Number of ghost cell layers
    /// @param face_is_periodic  [IMIN,IMAX,JMIN,JMAX,KMIN,KMAX] — true when the
    ///               face has a periodic 1-to-1 connection
    static void derivative(const MultiArray3D<Real>& a,
                           MultiArray3D<Real>& da,
                           int dir, Real dh, Int ng,
                           const bool face_is_periodic[6]);

    /// Interpolate a cell-centered 3D array to half-nodes (faces) in a given
    /// direction.  Used for face-metric computation.
    ///
    /// @param a     Input cell-centered array (ni x nj x nk)
    /// @param af    Output face array.  Size depends on @p dir :
    ///              dir=0 → (ni+1)×nj×nk,  dir=1 → ni×(nj+1)×nk,
    ///              dir=2 → ni×nj×(nk+1).  Must be pre-allocated.
    /// @param dir   Direction: 0 = ξ(i), 1 = η(j), 2 = ζ(k)
    /// @param ng    Number of ghost cell layers
    /// @param face_is_periodic  as above
    static void interp_to_faces(const MultiArray3D<Real>& a,
                                MultiArray3D<Real>& af,
                                int dir, Int ng,
                                const bool face_is_periodic[6]);

private:
    // ========================================================================
    // 1D interpolation to half-nodes  (cell center → i+1/2)
    // ========================================================================

    /// Centered 6-point interpolation.
    /// a[i] are cell-center values; returns value at i+1/2.
    /// Uses cells [i-2 .. i+3].  Caller guarantees these indices are valid.
    static Real interp_center_6pt(const Real* a, Int i);

    /// One-sided interpolation at left boundary — half-node a_{-1/2}
    /// (before the first cell).  Uses cells a[0..4].
    static Real interp_left_1st(const Real* a);

    /// One-sided interpolation at left boundary — half-node a_{1/2}
    /// (between cell 0 and cell 1).  Uses cells a[0..4].
    static Real interp_left_2nd(const Real* a);

    /// One-sided interpolation at left boundary — half-node a_{3/2}
    /// (between cell 1 and cell 2).  Uses cells a[0..4].
    static Real interp_left_3rd(const Real* a);

    /// One-sided interpolation at right boundary — half-node a_{n-1/2}
    /// (after the last cell).  Uses the last five cells a[n-5..n-1].
    static Real interp_right_1st(const Real* a);

    /// One-sided interpolation at right boundary — half-node a_{n-3/2}
    /// (between cell n-2 and cell n-1).  Uses cells a[n-5..n-1].
    static Real interp_right_2nd(const Real* a);

    /// One-sided interpolation at right boundary — half-node a_{n-5/2}
    /// (between cell n-3 and cell n-2).  Uses cells a[n-5..n-1].
    static Real interp_right_3rd(const Real* a);

    // ========================================================================
    // 1D differentiation from half-nodes  (i+1/2 → cell center)
    // ========================================================================

    /// Centered 6-point differentiation at cell i.
    /// ah[h] is the half-node between cells h and h+1.
    /// Uses ah[i-3 .. i+2].  Caller guarantees these indices are valid.
    static Real diff_center_6pt(const Real* ah, Int i, Real dh);

    /// One-sided differentiation at absolute left end — cell 0.
    /// ah has n+1 entries with convention ah[i] = a_{i-1/2}; i0=0 is the
    /// first half-node (a_{-1/2}).
    static Real diff_left_1st(const Real* ah, Int i0, Real dh);

    /// One-sided differentiation at absolute left end — cell 1.
    /// i0=0 is the first half-node (a_{-1/2}).
    static Real diff_left_2nd(const Real* ah, Int i0, Real dh);

    /// One-sided differentiation at absolute right end — last cell (n-1).
    /// i1=n is the last half-node (a_{n-1/2}).
    static Real diff_right_1st(const Real* ah, Int i1, Real dh);

    /// One-sided differentiation at absolute right end — cell n-2.
    /// i1=n is the last half-node (a_{n-1/2}).
    static Real diff_right_2nd(const Real* ah, Int i1, Real dh);

    // ========================================================================
    // 1D line helpers
    // ========================================================================

    /// Fill half-node array ah[0..n] from cell-center array a[0..n-1].
    /// Convention: ah[i] = a_{i-1/2}.  There are n+1 half-nodes for n cells.
    /// @param n     number of cell centers in this line
    /// @param left_periodic, right_periodic — boundary type switches
    static void interp_line(const Real* a, Real* ah, Int n,
                            bool left_periodic, bool right_periodic);

    /// Fill derivative array da[0..n-1] from half-node array ah[0..n].
    /// One-sided stencils at absolute ends (cells 0,1 and n-1,n-2);
    /// centered 6-point for cells 2..n-3.
    static void diff_line(const Real* ah, Real* da, Int n, Real dh,
                          bool left_periodic, bool right_periodic);
};

#include "interp_diff.hxx"
