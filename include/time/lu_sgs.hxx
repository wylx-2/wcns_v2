#pragma once

#include "lu_sgs.h"
#include "parallel/local_block.h"
#include <algorithm>
#include <cmath>

// ============================================================================
// LuSgs::advance
// ============================================================================

inline void LuSgs::advance(LocalBlock& lb, const Config& cfg, Real dt) {
    auto& f = lb.field;
    auto& g = lb.grid;
    const Int nci = f.ni();
    const Int ncj = f.nj();
    const Int nck = f.nk();

    // ---- Interior cell range (matches RK advance_stage) ----
    const Int i0 = 3, i1 = nci - 4;
    const Int j0 = 3, j1 = ncj - 4;
    const Int k0 = 3, k1 = nck - 4;

    const Real gamma   = cfg.gamma;
    const Real Pr      = cfg.Prandtl;
    const Real Re      = cfg.Re;
    const Real inv_dt  = 1.0 / dt;
    const bool has_viscous = (cfg.viscous_type != "none");

    // ---- Over-relaxation factor for LU-SGS diagonal ----
    // κ > 1 increases diagonal dominance; κ = 1.0 is the "standard" value.
    // For high-order spatial schemes (6th-order central), κ ≥ 4 is recommended.
    const Real kappa = cfg.lu_sgs_kappa;

    // ---- Allocate temporary work arrays ----
    MultiArray3D<Real> sigma_xi, sigma_eta, sigma_zeta;
    sigma_xi.allocate(nci, ncj, nck);
    sigma_eta.allocate(nci, ncj, nck);
    sigma_zeta.allocate(nci, ncj, nck);

    MultiArray3D<Real> diag;
    diag.allocate(nci, ncj, nck);

    MultiArray3D<Real> dQ_rho, dQ_rhou, dQ_rhov, dQ_rhow, dQ_rhoE;
    dQ_rho.allocate(nci, ncj, nck);
    dQ_rhou.allocate(nci, ncj, nck);
    dQ_rhov.allocate(nci, ncj, nck);
    dQ_rhow.allocate(nci, ncj, nck);
    dQ_rhoE.allocate(nci, ncj, nck);

    // =========================================================================
    // Step 1: Precompute spectral radii and diagonal at each interior cell
    // =========================================================================

    Real beta_visc = std::max(Real(4.0 / 3.0), gamma);

    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        Real rho = f.prim.rho(i, j, k);
        Real u   = f.prim.u(i, j, k);
        Real v   = f.prim.v(i, j, k);
        Real w   = f.prim.w(i, j, k);
        Real p   = f.prim.p(i, j, k);

        Real c_sound = std::sqrt(gamma * p / rho);
        Real omega = std::abs(g.jacobian(i, j, k));

        // ξ-direction inviscid spectral radius: ρ_ξ = |U_contra| + c · |S_ξ|
        Real u_contra = u * g.met_xi_x(i,j,k) + v * g.met_xi_y(i,j,k) + w * g.met_xi_z(i,j,k);
        Real s_face   = std::sqrt(g.met_xi_x(i,j,k)*g.met_xi_x(i,j,k)
                                + g.met_xi_y(i,j,k)*g.met_xi_y(i,j,k)
                                + g.met_xi_z(i,j,k)*g.met_xi_z(i,j,k));
        Real s_xi = std::abs(u_contra) + c_sound * s_face;
        sigma_xi(i, j, k) = s_xi;

        // η-direction
        u_contra = u * g.met_eta_x(i,j,k) + v * g.met_eta_y(i,j,k) + w * g.met_eta_z(i,j,k);
        s_face   = std::sqrt(g.met_eta_x(i,j,k)*g.met_eta_x(i,j,k)
                           + g.met_eta_y(i,j,k)*g.met_eta_y(i,j,k)
                           + g.met_eta_z(i,j,k)*g.met_eta_z(i,j,k));
        Real s_eta = std::abs(u_contra) + c_sound * s_face;
        sigma_eta(i, j, k) = s_eta;

        // ζ-direction
        u_contra = u * g.met_zeta_x(i,j,k) + v * g.met_zeta_y(i,j,k) + w * g.met_zeta_z(i,j,k);
        s_face   = std::sqrt(g.met_zeta_x(i,j,k)*g.met_zeta_x(i,j,k)
                           + g.met_zeta_y(i,j,k)*g.met_zeta_y(i,j,k)
                           + g.met_zeta_z(i,j,k)*g.met_zeta_z(i,j,k));
        Real s_zeta = std::abs(u_contra) + c_sound * s_face;
        sigma_zeta(i, j, k) = s_zeta;

        // Diagonal: D = Ω/Δt + κ·(ρ_ξ + ρ_η + ρ_ζ) + viscous
        Real d_val = omega * inv_dt + kappa * (s_xi + s_eta + s_zeta);

        if (has_viscous) {
            Real mu    = Real(1.0) / Re;
            Real mu_rho = mu / rho;
            Real factor = beta_visc * mu_rho / (Pr * omega);
            auto s2 = [](Real mx, Real my, Real mz) { return mx*mx + my*my + mz*mz; };
            Real rho_v = factor * (s2(g.met_xi_x(i,j,k),  g.met_xi_y(i,j,k),  g.met_xi_z(i,j,k))
                                 + s2(g.met_eta_x(i,j,k), g.met_eta_y(i,j,k), g.met_eta_z(i,j,k))
                                 + s2(g.met_zeta_x(i,j,k),g.met_zeta_y(i,j,k),g.met_zeta_z(i,j,k)));
            d_val += Real(2.0) * rho_v;
        }

        diag(i, j, k) = d_val;
    }}}

    // =========================================================================
    // Step 2: Forward sweep — i↑, j↑, k↑
    //   (D + L) ΔQ* = Ω·RHS
    //   → ΔQ*_i = D^{-1}[Ω·RHS_i + A^+_{i-1}ΔQ*_{i-1} + B^+_{j-1}ΔQ*_{j-1} + C^+_{k-1}ΔQ*_{k-1}]
    //   With scalar LU-SGS: A^+_{neighbor} ΔQ ≈ 0.5·max(σ_L, σ_R)·ΔQ
    // =========================================================================

    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        Real inv_diag = Real(1.0) / diag(i, j, k);

        // RHS term: Ω · rhs  (convert dQ/dt to flux-balance RHS)
        Real omega = std::abs(g.jacobian(i, j, k));
        Real d_rho  = omega * f.rhs.rho(i, j, k);
        Real d_rhou = omega * f.rhs.rhou(i, j, k);
        Real d_rhov = omega * f.rhs.rhov(i, j, k);
        Real d_rhow = omega * f.rhs.rhow(i, j, k);
        Real d_rhoE = omega * f.rhs.rhoE(i, j, k);

        // ---- L operator: contributions from i-1, j-1, k-1 ----
        // Use max of adjacent sigma for better diagonal dominance
        if (i > i0) {
            Real s_face = Real(0.5) * std::max(sigma_xi(i - 1, j, k), sigma_xi(i, j, k));
            d_rho  += s_face * dQ_rho(i - 1, j, k);
            d_rhou += s_face * dQ_rhou(i - 1, j, k);
            d_rhov += s_face * dQ_rhov(i - 1, j, k);
            d_rhow += s_face * dQ_rhow(i - 1, j, k);
            d_rhoE += s_face * dQ_rhoE(i - 1, j, k);
        }
        if (j > j0) {
            Real s_face = Real(0.5) * std::max(sigma_eta(i, j - 1, k), sigma_eta(i, j, k));
            d_rho  += s_face * dQ_rho(i, j - 1, k);
            d_rhou += s_face * dQ_rhou(i, j - 1, k);
            d_rhov += s_face * dQ_rhov(i, j - 1, k);
            d_rhow += s_face * dQ_rhow(i, j - 1, k);
            d_rhoE += s_face * dQ_rhoE(i, j - 1, k);
        }
        if (k > k0) {
            Real s_face = Real(0.5) * std::max(sigma_zeta(i, j, k - 1), sigma_zeta(i, j, k));
            d_rho  += s_face * dQ_rho(i, j, k - 1);
            d_rhou += s_face * dQ_rhou(i, j, k - 1);
            d_rhov += s_face * dQ_rhov(i, j, k - 1);
            d_rhow += s_face * dQ_rhow(i, j, k - 1);
            d_rhoE += s_face * dQ_rhoE(i, j, k - 1);
        }

        dQ_rho(i, j, k)  = d_rho  * inv_diag;
        dQ_rhou(i, j, k) = d_rhou * inv_diag;
        dQ_rhov(i, j, k) = d_rhov * inv_diag;
        dQ_rhow(i, j, k) = d_rhow * inv_diag;
        dQ_rhoE(i, j, k) = d_rhoE * inv_diag;
    }}}

    // =========================================================================
    // Step 3: Backward sweep — i↓, j↓, k↓
    //   (D + U) ΔQ = D ΔQ*
    //   → ΔQ_i = ΔQ*_i - D^{-1}[A^-_{i+1}ΔQ_{i+1} + B^-_{j+1}ΔQ_{j+1} + C^-_{k+1}ΔQ_{k+1}]
    //   With scalar LU-SGS: A^-_{neighbor} ΔQ ≈ 0.5·max(σ_L, σ_R)·ΔQ
    // =========================================================================

    for (Int k = k1; k >= k0; --k) {
    for (Int j = j1; j >= j0; --j) {
    for (Int i = i1; i >= i0; --i) {
        Real inv_diag = Real(1.0) / diag(i, j, k);

        Real corr_rho  = Real(0.0);
        Real corr_rhou = Real(0.0);
        Real corr_rhov = Real(0.0);
        Real corr_rhow = Real(0.0);
        Real corr_rhoE = Real(0.0);

        // ---- U operator: contributions from i+1, j+1, k+1 ----
        if (i < i1) {
            Real s_face = Real(0.5) * std::max(sigma_xi(i, j, k), sigma_xi(i + 1, j, k));
            corr_rho  += s_face * dQ_rho(i + 1, j, k);
            corr_rhou += s_face * dQ_rhou(i + 1, j, k);
            corr_rhov += s_face * dQ_rhov(i + 1, j, k);
            corr_rhow += s_face * dQ_rhow(i + 1, j, k);
            corr_rhoE += s_face * dQ_rhoE(i + 1, j, k);
        }
        if (j < j1) {
            Real s_face = Real(0.5) * std::max(sigma_eta(i, j, k), sigma_eta(i, j + 1, k));
            corr_rho  += s_face * dQ_rho(i, j + 1, k);
            corr_rhou += s_face * dQ_rhou(i, j + 1, k);
            corr_rhov += s_face * dQ_rhov(i, j + 1, k);
            corr_rhow += s_face * dQ_rhow(i, j + 1, k);
            corr_rhoE += s_face * dQ_rhoE(i, j + 1, k);
        }
        if (k < k1) {
            Real s_face = Real(0.5) * std::max(sigma_zeta(i, j, k), sigma_zeta(i, j, k + 1));
            corr_rho  += s_face * dQ_rho(i, j, k + 1);
            corr_rhou += s_face * dQ_rhou(i, j, k + 1);
            corr_rhov += s_face * dQ_rhov(i, j, k + 1);
            corr_rhow += s_face * dQ_rhow(i, j, k + 1);
            corr_rhoE += s_face * dQ_rhoE(i, j, k + 1);
        }

        // ΔQ = ΔQ* - D^{-1} · U_correction
        dQ_rho(i, j, k)  -= corr_rho  * inv_diag;
        dQ_rhou(i, j, k) -= corr_rhou * inv_diag;
        dQ_rhov(i, j, k) -= corr_rhov * inv_diag;
        dQ_rhow(i, j, k) -= corr_rhow * inv_diag;
        dQ_rhoE(i, j, k) -= corr_rhoE * inv_diag;
    }}}

    // =========================================================================
    // Step 4: Update conservative variables — Q^{n+1} = Q^n + ΔQ
    // =========================================================================

    for (Int k = k0; k <= k1; ++k) {
    for (Int j = j0; j <= j1; ++j) {
    for (Int i = i0; i <= i1; ++i) {
        f.cons.rho(i, j, k)  = f.Q0.rho(i, j, k)  + dQ_rho(i, j, k);
        f.cons.rhou(i, j, k) = f.Q0.rhou(i, j, k) + dQ_rhou(i, j, k);
        f.cons.rhov(i, j, k) = f.Q0.rhov(i, j, k) + dQ_rhov(i, j, k);
        f.cons.rhow(i, j, k) = f.Q0.rhow(i, j, k) + dQ_rhow(i, j, k);
        f.cons.rhoE(i, j, k) = f.Q0.rhoE(i, j, k) + dQ_rhoE(i, j, k);
    }}}
}
