#pragma once

#include "inviscid_rhs.h"
#include "wcns_v2/parallel/local_block.h"

// ============================================================================
// 6th-order centered difference coefficients (dh = 1.0 in computational space)
// ============================================================================
// From InterpDiff::diff_center_6pt:
//   da[i] =  75/(64*dh)*(ah[i+1] - ah[i])
//          - 25/(384*dh)*(ah[i+2] - ah[i-1])
//          +  3/(640*dh)*(ah[i+3] - ah[i-2])

namespace {
    constexpr Real c0 =  75.0 / 64.0;    //  1.171875
    constexpr Real c1 = -25.0 / 384.0;   // -0.065104166...
    constexpr Real c2 =   3.0 / 640.0;   //  0.0046875
}

// ============================================================================
// InviscidRHS::compute
// ============================================================================

inline void InviscidRHS::compute(LocalBlock& lb) {
    auto& f = lb.field;  // non-const — we write to f.rhs
    const Int nci = f.ni();
    const Int ncj = f.nj();
    const Int nck = f.nk();

    // Interior cell range (RHS only needed for core cells that time stepping updates)
    const Int i0 = 3, i1 = nci - 4;
    const Int j0 = 3, j1 = ncj - 4;
    const Int k0 = 3, k1 = nck - 4;

    // Initialize RHS to zero (only interior cells matter, but zero all for safety)
    f.rhs.fill(0.0);

    // ---- ξ-direction: differentiate inv_xi along i ----
    // inv_xi has size (nci+1)×ncj×nck; inv_xi(h,j,k) = face at h-1/2
    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        Real inv_J = Real(1.0) / std::abs(lb.grid.jacobian(i, j, k));

        // ∂F1/∂ξ
        Real dF1 = c0 * (f.inv_xi.f1(i+1,j,k) - f.inv_xi.f1(i,j,k))
                 + c1 * (f.inv_xi.f1(i+2,j,k) - f.inv_xi.f1(i-1,j,k))
                 + c2 * (f.inv_xi.f1(i+3,j,k) - f.inv_xi.f1(i-2,j,k));
        Real dF2 = c0 * (f.inv_xi.f2(i+1,j,k) - f.inv_xi.f2(i,j,k))
                 + c1 * (f.inv_xi.f2(i+2,j,k) - f.inv_xi.f2(i-1,j,k))
                 + c2 * (f.inv_xi.f2(i+3,j,k) - f.inv_xi.f2(i-2,j,k));
        Real dF3 = c0 * (f.inv_xi.f3(i+1,j,k) - f.inv_xi.f3(i,j,k))
                 + c1 * (f.inv_xi.f3(i+2,j,k) - f.inv_xi.f3(i-1,j,k))
                 + c2 * (f.inv_xi.f3(i+3,j,k) - f.inv_xi.f3(i-2,j,k));
        Real dF4 = c0 * (f.inv_xi.f4(i+1,j,k) - f.inv_xi.f4(i,j,k))
                 + c1 * (f.inv_xi.f4(i+2,j,k) - f.inv_xi.f4(i-1,j,k))
                 + c2 * (f.inv_xi.f4(i+3,j,k) - f.inv_xi.f4(i-2,j,k));
        Real dF5 = c0 * (f.inv_xi.f5(i+1,j,k) - f.inv_xi.f5(i,j,k))
                 + c1 * (f.inv_xi.f5(i+2,j,k) - f.inv_xi.f5(i-1,j,k))
                 + c2 * (f.inv_xi.f5(i+3,j,k) - f.inv_xi.f5(i-2,j,k));

        f.rhs.rho(i,j,k)  -= dF1 * inv_J;
        f.rhs.rhou(i,j,k) -= dF2 * inv_J;
        f.rhs.rhov(i,j,k) -= dF3 * inv_J;
        f.rhs.rhow(i,j,k) -= dF4 * inv_J;
        f.rhs.rhoE(i,j,k) -= dF5 * inv_J;
    }}}

    // ---- η-direction: differentiate inv_eta along j ----
    // inv_eta has size nci×(ncj+1)×nck; inv_eta(i,h,k) = face at h-1/2
    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        Real inv_J = Real(1.0) / std::abs(lb.grid.jacobian(i, j, k));

        Real dG1 = c0 * (f.inv_eta.f1(i,j+1,k) - f.inv_eta.f1(i,j,k))
                 + c1 * (f.inv_eta.f1(i,j+2,k) - f.inv_eta.f1(i,j-1,k))
                 + c2 * (f.inv_eta.f1(i,j+3,k) - f.inv_eta.f1(i,j-2,k));
        Real dG2 = c0 * (f.inv_eta.f2(i,j+1,k) - f.inv_eta.f2(i,j,k))
                 + c1 * (f.inv_eta.f2(i,j+2,k) - f.inv_eta.f2(i,j-1,k))
                 + c2 * (f.inv_eta.f2(i,j+3,k) - f.inv_eta.f2(i,j-2,k));
        Real dG3 = c0 * (f.inv_eta.f3(i,j+1,k) - f.inv_eta.f3(i,j,k))
                 + c1 * (f.inv_eta.f3(i,j+2,k) - f.inv_eta.f3(i,j-1,k))
                 + c2 * (f.inv_eta.f3(i,j+3,k) - f.inv_eta.f3(i,j-2,k));
        Real dG4 = c0 * (f.inv_eta.f4(i,j+1,k) - f.inv_eta.f4(i,j,k))
                 + c1 * (f.inv_eta.f4(i,j+2,k) - f.inv_eta.f4(i,j-1,k))
                 + c2 * (f.inv_eta.f4(i,j+3,k) - f.inv_eta.f4(i,j-2,k));
        Real dG5 = c0 * (f.inv_eta.f5(i,j+1,k) - f.inv_eta.f5(i,j,k))
                 + c1 * (f.inv_eta.f5(i,j+2,k) - f.inv_eta.f5(i,j-1,k))
                 + c2 * (f.inv_eta.f5(i,j+3,k) - f.inv_eta.f5(i,j-2,k));

        f.rhs.rho(i,j,k)  -= dG1 * inv_J;
        f.rhs.rhou(i,j,k) -= dG2 * inv_J;
        f.rhs.rhov(i,j,k) -= dG3 * inv_J;
        f.rhs.rhow(i,j,k) -= dG4 * inv_J;
        f.rhs.rhoE(i,j,k) -= dG5 * inv_J;
    }}}

    // ---- ζ-direction: differentiate inv_zeta along k ----
    // inv_zeta has size nci×ncj×(nck+1); inv_zeta(i,j,h) = face at h-1/2
    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        Real inv_J = Real(1.0) / std::abs(lb.grid.jacobian(i, j, k));

        Real dH1 = c0 * (f.inv_zeta.f1(i,j,k+1) - f.inv_zeta.f1(i,j,k))
                 + c1 * (f.inv_zeta.f1(i,j,k+2) - f.inv_zeta.f1(i,j,k-1))
                 + c2 * (f.inv_zeta.f1(i,j,k+3) - f.inv_zeta.f1(i,j,k-2));
        Real dH2 = c0 * (f.inv_zeta.f2(i,j,k+1) - f.inv_zeta.f2(i,j,k))
                 + c1 * (f.inv_zeta.f2(i,j,k+2) - f.inv_zeta.f2(i,j,k-1))
                 + c2 * (f.inv_zeta.f2(i,j,k+3) - f.inv_zeta.f2(i,j,k-2));
        Real dH3 = c0 * (f.inv_zeta.f3(i,j,k+1) - f.inv_zeta.f3(i,j,k))
                 + c1 * (f.inv_zeta.f3(i,j,k+2) - f.inv_zeta.f3(i,j,k-1))
                 + c2 * (f.inv_zeta.f3(i,j,k+3) - f.inv_zeta.f3(i,j,k-2));
        Real dH4 = c0 * (f.inv_zeta.f4(i,j,k+1) - f.inv_zeta.f4(i,j,k))
                 + c1 * (f.inv_zeta.f4(i,j,k+2) - f.inv_zeta.f4(i,j,k-1))
                 + c2 * (f.inv_zeta.f4(i,j,k+3) - f.inv_zeta.f4(i,j,k-2));
        Real dH5 = c0 * (f.inv_zeta.f5(i,j,k+1) - f.inv_zeta.f5(i,j,k))
                 + c1 * (f.inv_zeta.f5(i,j,k+2) - f.inv_zeta.f5(i,j,k-1))
                 + c2 * (f.inv_zeta.f5(i,j,k+3) - f.inv_zeta.f5(i,j,k-2));

        f.rhs.rho(i,j,k)  -= dH1 * inv_J;
        f.rhs.rhou(i,j,k) -= dH2 * inv_J;
        f.rhs.rhov(i,j,k) -= dH3 * inv_J;
        f.rhs.rhow(i,j,k) -= dH4 * inv_J;
        f.rhs.rhoE(i,j,k) -= dH5 * inv_J;
    }}}
}
