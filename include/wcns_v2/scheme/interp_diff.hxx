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

// Left 1st half-node: a_{1/2} = (315*a1 - 420*a2 + 378*a3 - 180*a4 + 35*a5) / 128
// a[0..4] = first five physical cells
inline Real InterpDiff::interp_left_1st(const Real* a) {
    return (315.0 * a[0] - 420.0 * a[1] + 378.0 * a[2]
          - 180.0 * a[3] +  35.0 * a[4]) / 128.0;
}

// Left 2nd half-node: a_{3/2} = (35*a1 + 140*a2 - 70*a3 + 28*a4 - 5*a5) / 128
inline Real InterpDiff::interp_left_2nd(const Real* a) {
    return (35.0 * a[0] + 140.0 * a[1] - 70.0 * a[2]
          + 28.0 * a[3] -   5.0 * a[4]) / 128.0;
}

// Right 1st half-node (last physical half-node): symmetric with left 1st
// a[0..4] = last five physical cells, a[4] is the last physical cell
inline Real InterpDiff::interp_right_1st(const Real* a) {
    return (35.0 * a[0] - 180.0 * a[1] + 378.0 * a[2]
          - 420.0 * a[3] + 315.0 * a[4]) / 128.0;
}

// Right 2nd half-node (second-to-last physical half-node): symmetric with left 2nd
inline Real InterpDiff::interp_right_2nd(const Real* a) {
    return (-5.0 * a[0] +  28.0 * a[1] - 70.0 * a[2]
           + 140.0 * a[3] +  35.0 * a[4]) / 128.0;
}

// ============================================================================
// 1D differentiation from half-nodes
// ============================================================================

// Centered 6-point:
// a_ζ(i) = 75/(64Δζ)*(a_{i+1/2} - a_{i-1/2})
//         - 25/(384Δζ)*(a_{i+3/2} - a_{i-3/2})
//         +  3/(640Δζ)*(a_{i+5/2} - a_{i-5/2})
// ah[h] is half-node at h+1/2
inline Real InterpDiff::diff_center_6pt(const Real* ah, Int i, Real dh) {
    const Real c0 =  75.0 / ( 64.0 * dh);
    const Real c1 = -25.0 / (384.0 * dh);
    const Real c2 =   3.0 / (640.0 * dh);
    return c0 * (ah[i]   - ah[i-1])
         + c1 * (ah[i+1] - ah[i-2])
         + c2 * (ah[i+2] - ah[i-3]);
}

// Left 1st cell: (-22 a_{1/2} + 17 a_{3/2} + 9 a_{5/2} - 5 a_{7/2} + a_{9/2}) / (24 Δξ)
// i0 = index of first physical cell (= ng); ah[i0] = a_{1/2}
inline Real InterpDiff::diff_left_1st(const Real* ah, Int i0, Real dh) {
    return (-22.0 * ah[i0] + 17.0 * ah[i0+1] + 9.0 * ah[i0+2]
            - 5.0 * ah[i0+3] +        ah[i0+4]) / (24.0 * dh);
}

// Left 2nd cell: (a_{1/2} - 27 a_{3/2} + 27 a_{5/2} - a_{7/2}) / (24 Δξ)
inline Real InterpDiff::diff_left_2nd(const Real* ah, Int i0, Real dh) {
    return (ah[i0] - 27.0 * ah[i0+1] + 27.0 * ah[i0+2] - ah[i0+3]) / (24.0 * dh);
}

// Right 1st cell (last physical cell i1):
// (22 a_{N-1/2} - 17 a_{N-3/2} - 9 a_{N-5/2} + 5 a_{N-7/2} - a_{N-9/2}) / (24 Δξ)
// i1 = index of last physical cell; ah[i1-1] = a_{N-1/2}
inline Real InterpDiff::diff_right_1st(const Real* ah, Int i1, Real dh) {
    return (22.0 * ah[i1-1] - 17.0 * ah[i1-2] - 9.0 * ah[i1-3]
            + 5.0 * ah[i1-4] -        ah[i1-5]) / (24.0 * dh);
}

// Right 2nd cell (second-to-last physical cell i1-1):
// (-a_{N-1/2} + 27 a_{N-3/2} - 27 a_{N-5/2} + a_{N-7/2}) / (24 Δξ)
inline Real InterpDiff::diff_right_2nd(const Real* ah, Int i1, Real dh) {
    return (-ah[i1-1] + 27.0 * ah[i1-2] - 27.0 * ah[i1-3] + ah[i1-4]) / (24.0 * dh);
}

// ============================================================================
// 1D line helpers
// ============================================================================

inline void InterpDiff::interp_line(const Real* a, Real* ah, Int n, Int ng,
                                     bool left_periodic, bool right_periodic) {
    // Indices of first and last physical cells
    Int i0 = ng;
    Int i1 = n - 1 - ng;

    for (Int h = 0; h < n - 1; ++h) {
        // --- Non-periodic left boundary: first two physical half-nodes ---
        if (!left_periodic) {
            if (h == i0)      { ah[h] = interp_left_1st(&a[i0]);   continue; }
            if (h == i0 + 1)  { ah[h] = interp_left_2nd(&a[i0]);   continue; }
        }

        // --- Non-periodic right boundary: last two physical half-nodes ---
        if (!right_periodic) {
            if (h == i1 - 1)  { ah[h] = interp_right_1st(&a[i1 - 4]); continue; }
            if (h == i1 - 2)  { ah[h] = interp_right_2nd(&a[i1 - 4]); continue; }
        }

        // --- Centered stencil: needs cells h-2 .. h+3 ---
        if (h >= 2 && h + 3 < n) {
            ah[h] = interp_center_6pt(a, h);
            continue;
        }

        // --- Out-of-bounds: use nearest one-sided stencil ---
        if (h < 2) {
            // Left edge of array — use left one-sided formula at the leftmost
            // valid position.  With periodic BCs ghost data is valid; with
            // non-periodic BCs this only affects outer ghost half-nodes.
            // Use a[0] as base when h==0 (needs a[-2],a[-1] out of bounds),
            // and a[1] when h==1.
            Int base = (h == 0) ? 0 : 1;
            ah[h] = interp_left_1st(&a[base]);
        } else {
            // Right edge of array — symmetric logic
            // h+3 >= n, i.e. h >= n-3
            Int base = (h == n - 2) ? n - 5 : n - 6;
            if (base < 0) base = 0;
            ah[h] = interp_right_1st(&a[base]);
        }
    }
}

inline void InterpDiff::diff_line(const Real* ah, Real* da, Int n, Int ng, Real dh,
                                   bool left_periodic, bool right_periodic) {
    Int i0 = ng;
    Int i1 = n - 1 - ng;

    for (Int i = 0; i < n; ++i) {
        // --- Non-periodic left boundary: first two physical cells ---
        if (!left_periodic) {
            if (i == i0)      { da[i] = diff_left_1st(ah, i0, dh);  continue; }
            if (i == i0 + 1)  { da[i] = diff_left_2nd(ah, i0, dh);  continue; }
        }

        // --- Non-periodic right boundary: last two physical cells ---
        if (!right_periodic) {
            if (i == i1)      { da[i] = diff_right_1st(ah, i1, dh); continue; }
            if (i == i1 - 1)  { da[i] = diff_right_2nd(ah, i1, dh); continue; }
        }

        // --- Centered stencil: needs half-nodes i-3 .. i+2 ---
        if (i >= 3 && i + 2 <= n - 2) {  // i+2 <= n-2 means i+2 < n-1 (half-node count)
            da[i] = diff_center_6pt(ah, i, dh);
            continue;
        }

        // --- Out-of-bounds: use nearest one-sided stencil ---
        if (i < 3) {
            // Left edge — use left one-sided at the leftmost valid index.
            // For i=0,1,2 the centered stencil overflows the half-node array.
            // diff_left_1st uses ah[i0..i0+4] — shift to use the leftmost
            // available half-nodes by using i as the "boundary".
            Int base = i;
            if (base + 4 < n - 1) {
                da[i] = diff_left_1st(ah, base, dh);
            } else {
                // Fallback: use boundary value
                da[i] = diff_left_1st(ah, i0, dh);
            }
        } else {
            // Right edge — use right one-sided
            Int last = i;
            if (last >= 5) {
                da[i] = diff_right_1st(ah, last, dh);
            } else {
                da[i] = diff_right_1st(ah, i1, dh);
            }
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

            interp_line(a_line, ah_line, ni, ng, lp, rp);
            diff_line(ah_line, da_line, ni, ng, dh, lp, rp);
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

            interp_line(a_line.data(), ah_line.data(), nj, ng, lp, rp);

            // Scatter half-nodes back to 3D array
            for (Int j = 0; j < nj + 1; ++j) ah(i, j, k) = ah_line[j];

            diff_line(ah_line.data(), da_line_buf.data(), nj, ng, dh, lp, rp);

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

            interp_line(a_line.data(), ah_line.data(), nk, ng, lp, rp);

            for (Int k = 0; k < nk + 1; ++k) ah(i, j, k) = ah_line[k];

            diff_line(ah_line.data(), da_line_buf.data(), nk, ng, dh, lp, rp);

            for (Int k = 0; k < nk; ++k) da(i, j, k) = da_line_buf[k];
        }}
    }
}

inline void InterpDiff::interp_to_faces(const MultiArray3D<Real>& a,
                                          MultiArray3D<Real>& af,
                                          int dir, Int ng,
                                          const bool face_is_periodic[6]) {
    Int ni = a.ni(), nj = a.nj(), nk = a.nk();

    if (dir == 0) {
        // ξ faces: af size = (ni+1, nj, nk)
        bool lp = face_is_periodic[0];
        bool rp = face_is_periodic[1];

        for (Int k = 0; k < nk; ++k) {
        for (Int j = 0; j < nj; ++j) {
            const Real* a_line = &a(0, j, k);
            Real*      af_line = &af(0, j, k);

            interp_line(a_line, af_line, ni, ng, lp, rp);
        }}
    } else if (dir == 1) {
        // η faces: af size = (ni, nj+1, nk)
        bool lp = face_is_periodic[2];
        bool rp = face_is_periodic[3];

        for (Int k = 0; k < nk; ++k) {
        for (Int i = 0; i < ni; ++i) {
            std::vector<Real> a_line(nj), af_line(nj + 1);
            for (Int j = 0; j < nj; ++j) a_line[j] = a(i, j, k);

            interp_line(a_line.data(), af_line.data(), nj, ng, lp, rp);

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

            interp_line(a_line.data(), af_line.data(), nk, ng, lp, rp);

            for (Int k = 0; k < nk + 1; ++k) af(i, j, k) = af_line[k];
        }}
    }
}
