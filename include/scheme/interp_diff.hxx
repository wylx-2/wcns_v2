#pragma once

#include "interp_diff.h"
#include <algorithm>

// ============================================================================
// 1D interpolation to half-nodes
// ============================================================================

// Centered 6-point:  a_{i+1/2} = 75/128*(a_i + a_{i+1})
//                                - 25/256*(a_{i-1} + a_{i+2})
//                                +  3/256*(a_{i-2} + a_{i+3})
inline Real InterpDiff::interp_center_6pt(const Real* a, Int i) {
    return (75.0 / 128.0) * (a[i] + a[i+1])
         - (25.0 / 256.0) * (a[i-1] + a[i+2])
         + ( 3.0 / 256.0) * (a[i-2] + a[i+3]);
}

// Left 1st half-node: a_{-1/2} (before first cell)
// Uses cells a[0..4], extrapolated to position x = -1/2
inline Real InterpDiff::interp_left_1st(const Real* a) {
    return (315.0 * a[0] - 420.0 * a[1] + 378.0 * a[2]
          - 180.0 * a[3] +  35.0 * a[4]) / 128.0;
}

// Left 2nd half-node: a_{1/2} (between cell 0 and cell 1)
inline Real InterpDiff::interp_left_2nd(const Real* a) {
    return (35.0 * a[0] + 140.0 * a[1] - 70.0 * a[2]
          + 28.0 * a[3] -   5.0 * a[4]) / 128.0;
}

// Left 3rd half-node: a_{3/2} (between cell 1 and cell 2)
inline Real InterpDiff::interp_left_3rd(const Real* a) {
    return (-5.0 * a[0] +  60.0 * a[1] + 90.0 * a[2]
           - 20.0 * a[3] +   3.0 * a[4]) / 128.0;
}

// Right 1st half-node: a_{n-1/2} (after last cell)
// a[0..4] = last five cells a[n-5..n-1]; symmetric with left 1st
inline Real InterpDiff::interp_right_1st(const Real* a) {
    return (35.0 * a[0] - 180.0 * a[1] + 378.0 * a[2]
          - 420.0 * a[3] + 315.0 * a[4]) / 128.0;
}

// Right 2nd half-node: a_{n-3/2} (between cell n-2 and cell n-1)
// Symmetric with left 2nd
inline Real InterpDiff::interp_right_2nd(const Real* a) {
    return (-5.0 * a[0] +  28.0 * a[1] - 70.0 * a[2]
           + 140.0 * a[3] +  35.0 * a[4]) / 128.0;
}

// Right 3rd half-node: a_{n-5/2} (between cell n-3 and cell n-2)
// Symmetric with left 3rd
inline Real InterpDiff::interp_right_3rd(const Real* a) {
    return ( 3.0 * a[0] -  20.0 * a[1] + 90.0 * a[2]
           + 60.0 * a[3] -   5.0 * a[4]) / 128.0;
}

// ============================================================================
// 1D differentiation from half-nodes
// ============================================================================

// Centered 6-point:
// Convention: ah[i] = a_{i-1/2}  (half-node at position i-1/2)
// a_ζ(i) = 75/(64Δζ)*(a_{i+1/2} - a_{i-1/2})
//         - 25/(384Δζ)*(a_{i+3/2} - a_{i-3/2})
//         +  3/(640Δζ)*(a_{i+5/2} - a_{i-5/2})
inline Real InterpDiff::diff_center_6pt(const Real* ah, Int i, Real dh) {
    const Real c0 =  75.0 / ( 64.0 * dh);
    const Real c1 = -25.0 / (384.0 * dh);
    const Real c2 =   3.0 / (640.0 * dh);
    return c0 * (ah[i+1] - ah[i])
         + c1 * (ah[i+2] - ah[i-1])
         + c2 * (ah[i+3] - ah[i-2]);
}

// Left 1st cell (i=0): uses half-nodes a_{-1/2}..a_{7/2} = ah[0..4]
inline Real InterpDiff::diff_left_1st(const Real* ah, Int i0, Real dh) {
    return (-22.0 * ah[i0] + 17.0 * ah[i0+1] + 9.0 * ah[i0+2]
            - 5.0 * ah[i0+3] +        ah[i0+4]) / (24.0 * dh);
}

// Left 2nd cell (i=1): uses half-nodes a_{-1/2}..a_{5/2} = ah[0..3]
inline Real InterpDiff::diff_left_2nd(const Real* ah, Int i0, Real dh) {
    return (ah[i0] - 27.0 * ah[i0+1] + 27.0 * ah[i0+2] - ah[i0+3]) / (24.0 * dh);
}

// Right 1st cell (last cell i=n-1):
// Uses half-nodes a_{n-1/2}..a_{n-9/2} = ah[i1..i1-4] where i1=n
inline Real InterpDiff::diff_right_1st(const Real* ah, Int i1, Real dh) {
    return (22.0 * ah[i1] - 17.0 * ah[i1-1] - 9.0 * ah[i1-2]
            + 5.0 * ah[i1-3] -        ah[i1-4]) / (24.0 * dh);
}

// Right 2nd cell (second-to-last cell i=n-2):
// Uses half-nodes a_{n-1/2}..a_{n-7/2} = ah[i1..i1-3] where i1=n
inline Real InterpDiff::diff_right_2nd(const Real* ah, Int i1, Real dh) {
    return (-ah[i1] + 27.0 * ah[i1-1] - 27.0 * ah[i1-2] + ah[i1-3]) / (24.0 * dh);
}

// ============================================================================
// 1D line helpers
// ============================================================================

inline void InterpDiff::interp_line(const Real* a, Real* ah, Int n,
                                     bool left_periodic, bool right_periodic) {
    (void)left_periodic;   // one-sided closure at absolute ends is independent
    (void)right_periodic;  // of periodicity; ghost data must be valid externally

    // Produce n+1 half-nodes ah[0..n] with the convention ah[i] = a_{i-1/2}.
    // One-sided formulas at the absolute ends (h=0,1,2 on the left,
    // h=n,n-1,n-2 on the right); centered 6-point for h=3..n-3.
    //
    // One-sided closure is applied regardless of periodicity — the array
    // ends always lack stencil support on one side.  Periodicity only
    // guarantees that ghost-cell data is valid for the centered stencil.

    for (Int h = 0; h <= n; ++h) {
        // --- Absolute left end: first three half-nodes ---
        //     ah[0]=a_{-1/2}, ah[1]=a_{1/2}, ah[2]=a_{3/2}
        //     Stencil uses cells a[0..4]; ghost cells must be valid.
        if (h == 0)  { ah[h] = interp_left_1st(&a[0]);  continue; }
        if (h == 1)  { ah[h] = interp_left_2nd(&a[0]);  continue; }
        if (h == 2)  { ah[h] = interp_left_3rd(&a[0]);  continue; }

        // --- Absolute right end: last three half-nodes ---
        //     ah[n]=a_{n-1/2}, ah[n-1]=a_{n-3/2}, ah[n-2]=a_{n-5/2}
        //     Stencil uses cells a[n-5..n-1].
        if (h == n)      { ah[h] = interp_right_1st(&a[n - 5]); continue; }
        if (h == n - 1)  { ah[h] = interp_right_2nd(&a[n - 5]); continue; }
        if (h == n - 2)  { ah[h] = interp_right_3rd(&a[n - 5]); continue; }

        // --- Centered stencil: needs cells (h-3) .. (h+2) ---
        //     interp_center_6pt(a, h-1) returns a_{(h-1)+1/2} = a_{h-1/2}
        if (h >= 3 && h + 2 < n) {
            ah[h] = interp_center_6pt(a, h - 1);
            continue;
        }

        // --- Fallback (should only be reached for tiny arrays n < 6) ---
        if (h <= 1) {
            ah[h] = interp_left_1st(&a[0]);
        } else {
            Int base = (n >= 5) ? n - 5 : 0;
            ah[h] = interp_right_1st(&a[base]);
        }
    }
}

inline void InterpDiff::diff_line(const Real* ah, Real* da, Int n, Real dh,
                                   bool left_periodic, bool right_periodic) {
    (void)left_periodic;   // one-sided closure at absolute ends is independent
    (void)right_periodic;  // of periodicity; ghost data must be valid externally

    // Compute derivatives at n cell centers from n+1 half-nodes ah[0..n].
    // Convention: ah[i] = a_{i-1/2}.
    //
    // One-sided stencils at the absolute ends: cells i=0,1 (left) and
    // i=n-1,n-2 (right).  Centered 6-point for cells i=2..n-3.
    //
    // As with interp_line, one-sided closure is applied regardless of
    // periodicity — it addresses stencil availability at the array ends.
    //
    // diff_center_6pt needs half-nodes i-2 .. i+3 → valid when i>=2, i+3<=n.

    for (Int i = 0; i < n; ++i) {
        // --- Absolute left end: first two cells ---
        if (i == 0)  { da[i] = diff_left_1st(ah, 0, dh);  continue; }
        if (i == 1)  { da[i] = diff_left_2nd(ah, 0, dh);  continue; }

        // --- Absolute right end: last two cells ---
        //     i1=n is the index of the last half-node a_{n-1/2}
        if (i == n - 1)  { da[i] = diff_right_1st(ah, n, dh); continue; }
        if (i == n - 2)  { da[i] = diff_right_2nd(ah, n, dh); continue; }

        // --- Centered stencil: needs half-nodes i-2 .. i+3 ---
        if (i >= 2 && i + 3 <= n) {
            da[i] = diff_center_6pt(ah, i, dh);
            continue;
        }

        // --- Fallback (should only be reached for tiny arrays n < 5) ---
        if (i <= 1) {
            da[i] = diff_left_1st(ah, 0, dh);
        } else {
            da[i] = diff_right_1st(ah, n, dh);
        }
    }
}

// ============================================================================
// 3D array operations
// ============================================================================

inline void InterpDiff::derivative(const MultiArray3D<Real>& a,
                                     MultiArray3D<Real>& da,
                                     int dir, Real dh, Int ng,
                                     const bool face_is_periodic[6]) {
    (void)ng;  // no longer used internally; kept for API compatibility
    Int ni = a.ni(), nj = a.nj(), nk = a.nk();

    if (dir == 0) {
        // ξ direction: half-node array (ni+1, nj, nk)
        MultiArray3D<Real> ah;
        ah.allocate(ni + 1, nj, nk);

        bool lp = face_is_periodic[0];  // IMIN
        bool rp = face_is_periodic[1];  // IMAX

        for (Int k = 0; k < nk; ++k) {
        for (Int j = 0; j < nj; ++j) {
            // Extract 1D line
            const Real* a_line = &a(0, j, k);
            Real*      ah_line = &ah(0, j, k);
            Real*      da_line = &da(0, j, k);

            interp_line(a_line, ah_line, ni, lp, rp);
            diff_line(ah_line, da_line, ni, dh, lp, rp);
        }}
    } else if (dir == 1) {
        // η direction: half-node array (ni, nj+1, nk)
        MultiArray3D<Real> ah;
        ah.allocate(ni, nj + 1, nk);

        bool lp = face_is_periodic[2];  // JMIN
        bool rp = face_is_periodic[3];  // JMAX

        // Temporary 1D arrays on stack (max size ~71 for our grids)
        for (Int k = 0; k < nk; ++k) {
        for (Int i = 0; i < ni; ++i) {
            // Gather 1D line along j
            // We need temporary buffers since MultiArray3D is contiguous in i
            // and we can't get a contiguous line along j directly
            // Use stack allocation for small lines
            // nj is at most ~71 for our test grids → fine on stack
            std::vector<Real> a_line(nj), ah_line(nj + 1), da_line_buf(nj);
            for (Int j = 0; j < nj; ++j) a_line[j] = a(i, j, k);

            interp_line(a_line.data(), ah_line.data(), nj, lp, rp);

            // Scatter half-nodes back to 3D array
            for (Int j = 0; j < nj + 1; ++j) ah(i, j, k) = ah_line[j];

            diff_line(ah_line.data(), da_line_buf.data(), nj, dh, lp, rp);

            for (Int j = 0; j < nj; ++j) da(i, j, k) = da_line_buf[j];
        }}
    } else {
        // ζ direction: half-node array (ni, nj, nk+1)
        MultiArray3D<Real> ah;
        ah.allocate(ni, nj, nk + 1);

        bool lp = face_is_periodic[4];  // KMIN
        bool rp = face_is_periodic[5];  // KMAX

        for (Int j = 0; j < nj; ++j) {
        for (Int i = 0; i < ni; ++i) {
            // Gather 1D line along k
            std::vector<Real> a_line(nk), ah_line(nk + 1), da_line_buf(nk);
            for (Int k = 0; k < nk; ++k) a_line[k] = a(i, j, k);

            interp_line(a_line.data(), ah_line.data(), nk, lp, rp);

            for (Int k = 0; k < nk + 1; ++k) ah(i, j, k) = ah_line[k];

            diff_line(ah_line.data(), da_line_buf.data(), nk, dh, lp, rp);

            for (Int k = 0; k < nk; ++k) da(i, j, k) = da_line_buf[k];
        }}
    }
}

inline void InterpDiff::interp_to_faces(const MultiArray3D<Real>& a,
                                          MultiArray3D<Real>& af,
                                          int dir, Int ng,
                                          const bool face_is_periodic[6]) {
    (void)ng;  // no longer used internally; kept for API compatibility
    Int ni = a.ni(), nj = a.nj(), nk = a.nk();

    if (dir == 0) {
        // ξ faces: af size = (ni+1, nj, nk)
        bool lp = face_is_periodic[0];
        bool rp = face_is_periodic[1];

        for (Int k = 0; k < nk; ++k) {
        for (Int j = 0; j < nj; ++j) {
            const Real* a_line = &a(0, j, k);
            Real*      af_line = &af(0, j, k);

            interp_line(a_line, af_line, ni, lp, rp);
        }}
    } else if (dir == 1) {
        // η faces: af size = (ni, nj+1, nk)
        bool lp = face_is_periodic[2];
        bool rp = face_is_periodic[3];

        for (Int k = 0; k < nk; ++k) {
        for (Int i = 0; i < ni; ++i) {
            std::vector<Real> a_line(nj), af_line(nj + 1);
            for (Int j = 0; j < nj; ++j) a_line[j] = a(i, j, k);

            interp_line(a_line.data(), af_line.data(), nj, lp, rp);

            for (Int j = 0; j < nj + 1; ++j) af(i, j, k) = af_line[j];
        }}
    } else {
        // ζ faces: af size = (ni, nj, nk+1)
        bool lp = face_is_periodic[4];
        bool rp = face_is_periodic[5];

        for (Int j = 0; j < nj; ++j) {
        for (Int i = 0; i < ni; ++i) {
            std::vector<Real> a_line(nk), af_line(nk + 1);
            for (Int k = 0; k < nk; ++k) a_line[k] = a(i, j, k);

            interp_line(a_line.data(), af_line.data(), nk, lp, rp);

            for (Int k = 0; k < nk + 1; ++k) af(i, j, k) = af_line[k];
        }}
    }
}

// ============================================================================
// derivative_from_faces — 3D diff from pre-computed face values
// ============================================================================

inline void InterpDiff::derivative_from_faces(const MultiArray3D<Real>& af,
                                               MultiArray3D<Real>& da,
                                               int dir, Real dh,
                                               const bool face_is_periodic[6]) {
    // af is already a face array (half-node values).  We only need diff_line.
    // For dir=0: af is (ni+1, nj, nk) → da is (ni, nj, nk)
    // For dir=1: af is (ni, nj+1, nk) → da is (ni, nj, nk)
    // For dir=2: af is (ni, nj, nk+1) → da is (ni, nj, nk)

    Int ni = da.ni(), nj = da.nj(), nk = da.nk();

    if (dir == 0) {
        // ξ direction: af has ni+1 faces (i+1/2), da has ni cells
        bool lp = face_is_periodic[0];  // IMIN
        bool rp = face_is_periodic[1];  // IMAX

        for (Int k = 0; k < nk; ++k) {
        for (Int j = 0; j < nj; ++j) {
            const Real* af_line = &af(0, j, k);
            Real*       da_line = &da(0, j, k);

            diff_line(af_line, da_line, ni, dh, lp, rp);
        }}
    } else if (dir == 1) {
        // η direction: af has nj+1 faces (j+1/2), da has nj cells
        bool lp = face_is_periodic[2];  // JMIN
        bool rp = face_is_periodic[3];  // JMAX

        for (Int k = 0; k < nk; ++k) {
        for (Int i = 0; i < ni; ++i) {
            std::vector<Real> af_line(nj + 1), da_line_buf(nj);
            for (Int j = 0; j < nj + 1; ++j) af_line[j] = af(i, j, k);

            diff_line(af_line.data(), da_line_buf.data(), nj, dh, lp, rp);

            for (Int j = 0; j < nj; ++j) da(i, j, k) = da_line_buf[j];
        }}
    } else {
        // ζ direction: af has nk+1 faces (k+1/2), da has nk cells
        bool lp = face_is_periodic[4];  // KMIN
        bool rp = face_is_periodic[5];  // KMAX

        for (Int j = 0; j < nj; ++j) {
        for (Int i = 0; i < ni; ++i) {
            std::vector<Real> af_line(nk + 1), da_line_buf(nk);
            for (Int k = 0; k < nk + 1; ++k) af_line[k] = af(i, j, k);

            diff_line(af_line.data(), da_line_buf.data(), nk, dh, lp, rp);

            for (Int k = 0; k < nk; ++k) da(i, j, k) = da_line_buf[k];
        }}
    }
}
