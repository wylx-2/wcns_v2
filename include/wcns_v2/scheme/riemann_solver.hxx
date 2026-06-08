#pragma once

#include "riemann_solver.h"
#include "wcns_v2/parallel/local_block.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

// ============================================================================
// detail: directional loop helpers — eliminate boilerplate across solvers
// ============================================================================

namespace detail {

/// Signature of a per-face flux kernel.
/// qL, qR: left/right conservative states (5 components)
/// Sx, Sy, Sz: face area vector
/// gamma: specific heat ratio
/// eps: entropy fix coefficient (ignored by non-Roe solvers)
/// flux: [out] 5-component numerical flux through the face
using FluxKernel = void (*)(const Real qL[5], const Real qR[5],
                             Real Sx, Real Sy, Real Sz,
                             Real gamma, Real eps,
                             Real flux[5]);

inline void solve_xi_impl(LocalBlock& lb, Real gamma, Real eps, FluxKernel kernel) {
    const auto& grid = lb.grid;
    auto& f = lb.field;

    const Int ni = grid.nci + 1;
    const Int nj = grid.ncj;
    const Int nk = grid.nck;

    for (Int k = 0; k < nk; ++k) {
    for (Int j = 0; j < nj; ++j) {
    for (Int i = 0; i < ni; ++i) {
        Real qL[5] = { f.ql_xi.rho(i,j,k),  f.ql_xi.rhou(i,j,k),
                       f.ql_xi.rhov(i,j,k), f.ql_xi.rhow(i,j,k),
                       f.ql_xi.rhoE(i,j,k) };
        Real qR[5] = { f.qr_xi.rho(i,j,k),  f.qr_xi.rhou(i,j,k),
                       f.qr_xi.rhov(i,j,k), f.qr_xi.rhow(i,j,k),
                       f.qr_xi.rhoE(i,j,k) };

        Real flux[5];
        kernel(qL, qR,
               grid.face_xi_x(i,j,k), grid.face_xi_y(i,j,k), grid.face_xi_z(i,j,k),
               gamma, eps, flux);

        f.inv_xi.f1(i,j,k) = flux[0];
        f.inv_xi.f2(i,j,k) = flux[1];
        f.inv_xi.f3(i,j,k) = flux[2];
        f.inv_xi.f4(i,j,k) = flux[3];
        f.inv_xi.f5(i,j,k) = flux[4];
    }}}
}

inline void solve_eta_impl(LocalBlock& lb, Real gamma, Real eps, FluxKernel kernel) {
    const auto& grid = lb.grid;
    auto& f = lb.field;

    const Int ni = grid.nci;
    const Int nj = grid.ncj + 1;
    const Int nk = grid.nck;

    for (Int k = 0; k < nk; ++k) {
    for (Int j = 0; j < nj; ++j) {
    for (Int i = 0; i < ni; ++i) {
        Real qL[5] = { f.ql_eta.rho(i,j,k),  f.ql_eta.rhou(i,j,k),
                       f.ql_eta.rhov(i,j,k), f.ql_eta.rhow(i,j,k),
                       f.ql_eta.rhoE(i,j,k) };
        Real qR[5] = { f.qr_eta.rho(i,j,k),  f.qr_eta.rhou(i,j,k),
                       f.qr_eta.rhov(i,j,k), f.qr_eta.rhow(i,j,k),
                       f.qr_eta.rhoE(i,j,k) };

        Real flux[5];
        kernel(qL, qR,
               grid.face_eta_x(i,j,k), grid.face_eta_y(i,j,k), grid.face_eta_z(i,j,k),
               gamma, eps, flux);

        f.inv_eta.f1(i,j,k) = flux[0];
        f.inv_eta.f2(i,j,k) = flux[1];
        f.inv_eta.f3(i,j,k) = flux[2];
        f.inv_eta.f4(i,j,k) = flux[3];
        f.inv_eta.f5(i,j,k) = flux[4];
    }}}
}

inline void solve_zeta_impl(LocalBlock& lb, Real gamma, Real eps, FluxKernel kernel) {
    const auto& grid = lb.grid;
    auto& f = lb.field;

    const Int ni = grid.nci;
    const Int nj = grid.ncj;
    const Int nk = grid.nck + 1;

    for (Int k = 0; k < nk; ++k) {
    for (Int j = 0; j < nj; ++j) {
    for (Int i = 0; i < ni; ++i) {
        Real qL[5] = { f.ql_zeta.rho(i,j,k),  f.ql_zeta.rhou(i,j,k),
                       f.ql_zeta.rhov(i,j,k), f.ql_zeta.rhow(i,j,k),
                       f.ql_zeta.rhoE(i,j,k) };
        Real qR[5] = { f.qr_zeta.rho(i,j,k),  f.qr_zeta.rhou(i,j,k),
                       f.qr_zeta.rhov(i,j,k), f.qr_zeta.rhow(i,j,k),
                       f.qr_zeta.rhoE(i,j,k) };

        Real flux[5];
        kernel(qL, qR,
               grid.face_zeta_x(i,j,k), grid.face_zeta_y(i,j,k), grid.face_zeta_z(i,j,k),
               gamma, eps, flux);

        f.inv_zeta.f1(i,j,k) = flux[0];
        f.inv_zeta.f2(i,j,k) = flux[1];
        f.inv_zeta.f3(i,j,k) = flux[2];
        f.inv_zeta.f4(i,j,k) = flux[3];
        f.inv_zeta.f5(i,j,k) = flux[4];
    }}}
}

} // namespace detail

// ============================================================================
// Factory
// ============================================================================

inline std::unique_ptr<RiemannSolverBase>
RiemannSolverBase::create(const Config& cfg) {
    if (cfg.riemann_type == "roe") {
        return std::make_unique<RiemannSolverRoe>();
    }
    if (cfg.riemann_type == "rusanov") {
        return std::make_unique<RiemannSolverRusanov>();
    }
    if (cfg.riemann_type == "hll") {
        return std::make_unique<RiemannSolverHLL>();
    }
    if (cfg.riemann_type == "hllc") {
        return std::make_unique<RiemannSolverHLLC>();
    }
    throw std::runtime_error("RiemannSolverBase::create: unknown type \""
                             + cfg.riemann_type + "\"");
}

// ============================================================================
// Roe flux kernel — one face
// ============================================================================

inline void RiemannSolverRoe::roe_flux(
    const Real qL[5], const Real qR[5],
    Real Sx, Real Sy, Real Sz,
    Real gamma, Real eps,
    Real flux[5])
{
    // ---- 1. Extract primitive variables from conservative states ----

    // Left state
    Real rhoL = qL[0];
    Real uL   = qL[1] / rhoL;
    Real vL   = qL[2] / rhoL;
    Real wL   = qL[3] / rhoL;
    Real keL  = 0.5 * (qL[1]*uL + qL[2]*vL + qL[3]*wL);  // 0.5 * rho*|V|²
    Real pL   = (gamma - 1.0) * (qL[4] - keL);
    Real HL   = (qL[4] + pL) / rhoL;  // total specific enthalpy

    // Right state
    Real rhoR = qR[0];
    Real uR   = qR[1] / rhoR;
    Real vR   = qR[2] / rhoR;
    Real wR   = qR[3] / rhoR;
    Real keR  = 0.5 * (qR[1]*uR + qR[2]*vR + qR[3]*wR);
    Real pR   = (gamma - 1.0) * (qR[4] - keR);
    Real HR   = (qR[4] + pR) / rhoR;

    // ---- 2. Physical flux through the face for each state ----

    // Contravariant velocity:  U = Sx*u + Sy*v + Sz*w
    Real UL = Sx*uL + Sy*vL + Sz*wL;
    Real UR = Sx*uR + Sy*vR + Sz*wR;

    // F_face = [rho*U, rhou*U + p*Sx, rhov*U + p*Sy, rhow*U + p*Sz, (rhoE+p)*U]
    Real FL[5], FR[5];
    FL[0] = qL[0] * UL;
    FL[1] = qL[1] * UL + pL * Sx;
    FL[2] = qL[2] * UL + pL * Sy;
    FL[3] = qL[3] * UL + pL * Sz;
    FL[4] = (qL[4] + pL) * UL;

    FR[0] = qR[0] * UR;
    FR[1] = qR[1] * UR + pR * Sx;
    FR[2] = qR[2] * UR + pR * Sy;
    FR[3] = qR[3] * UR + pR * Sz;
    FR[4] = (qR[4] + pR) * UR;

    // ---- 3. Roe average state ----

    Real R  = std::sqrt(rhoR / rhoL);         // √(ρ_R/ρ_L)
    Real Rp1 = R + 1.0;                       // R + 1
    Real rho_tilde = R * rhoL;                // √(ρ_L·ρ_R)

    Real u_tilde = (R * uR + uL) / Rp1;
    Real v_tilde = (R * vR + vL) / Rp1;
    Real w_tilde = (R * wR + wL) / Rp1;
    Real H_tilde = (R * HR + HL) / Rp1;

    Real q2_tilde = u_tilde*u_tilde + v_tilde*v_tilde + w_tilde*w_tilde;
    Real c_tilde  = std::sqrt(std::max(0.0,
                        (gamma - 1.0) * (H_tilde - 0.5 * q2_tilde)));

    // Face area magnitude and Roe contravariant velocity
    Real absS = std::sqrt(Sx*Sx + Sy*Sy + Sz*Sz);
    Real U_tilde = Sx * u_tilde + Sy * v_tilde + Sz * w_tilde;  // c̃·|S| → λ

    // ---- 4. Eigenvalues of face-normal flux Jacobian ----

    Real lam1 = U_tilde - c_tilde * absS;  // left acoustic
    Real lam2 = U_tilde;                    // entropy + 2 shear (all = Ũ)
    Real lam5 = U_tilde + c_tilde * absS;  // right acoustic

    // ---- 5. Entropy fix (Harten correction) ----

    Real delta = eps * c_tilde * absS;      // threshold

    auto entropy_fix = [delta](Real lam) -> Real {
        Real al = std::abs(lam);
        if (al >= delta) return al;
        return 0.5 * (lam * lam + delta * delta) / delta;
    };

    Real a_lam1 = entropy_fix(lam1);
    Real a_lam2 = entropy_fix(lam2);        // same for λ₂, λ₃, λ₄
    Real a_lam5 = entropy_fix(lam5);

    // ---- 6. Primitive variable jumps ----

    Real drho = rhoR - rhoL;
    Real dp   = pR   - pL;
    Real du   = uR   - uL;
    Real dv   = vR   - vL;
    Real dw   = wR   - wL;

    // ---- 7. Face-normal and tangential velocity jumps ----

    Real inv_absS = 1.0 / absS;
    Real nx = Sx * inv_absS;
    Real ny = Sy * inv_absS;
    Real nz = Sz * inv_absS;

    // Δu_n = ΔV · n̂  (face-normal velocity jump)
    Real du_n = du * nx + dv * ny + dw * nz;

    // ΔV_t = ΔV - (ΔV·n̂) n̂  (tangential velocity jump vector)
    Real du_tx = du - du_n * nx;
    Real du_ty = dv - du_n * ny;
    Real du_tz = dw - du_n * nz;

    // ---- 8. Wave strengths (α = R⁻¹ · ΔQ) ----

    Real c2 = c_tilde * c_tilde;
    Real inv_c2 = 1.0 / c2;

    Real alpha1 = (dp - rho_tilde * c_tilde * du_n) * 0.5 * inv_c2;  // left acoustic
    Real alpha2 = drho - dp * inv_c2;                                  // entropy
    Real alpha3x = rho_tilde * du_tx;                                   // shear (combined)
    Real alpha3y = rho_tilde * du_ty;
    Real alpha3z = rho_tilde * du_tz;
    Real alpha5 = (dp + rho_tilde * c_tilde * du_n) * 0.5 * inv_c2;  // right acoustic

    // ---- 9. Compact dissipation vector d = R·|Λ|·α ----

    Real d[5];
    // Continuity
    d[0] = a_lam2 * alpha2 + a_lam1 * alpha1 + a_lam5 * alpha5;

    // x-momentum
    d[1] = a_lam2 * (alpha2 * u_tilde + alpha3x)
         + a_lam1 * alpha1 * (u_tilde - c_tilde * nx)
         + a_lam5 * alpha5 * (u_tilde + c_tilde * nx);

    // y-momentum
    d[2] = a_lam2 * (alpha2 * v_tilde + alpha3y)
         + a_lam1 * alpha1 * (v_tilde - c_tilde * ny)
         + a_lam5 * alpha5 * (v_tilde + c_tilde * ny);

    // z-momentum
    d[3] = a_lam2 * (alpha2 * w_tilde + alpha3z)
         + a_lam1 * alpha1 * (w_tilde - c_tilde * nz)
         + a_lam5 * alpha5 * (w_tilde + c_tilde * nz);

    // Energy: shear contribution = ρ̃ * (ũ·ΔV_t), acoustic uses H̃ ± c̃*Ũ/|S|
    Real alpha3_dot_u = rho_tilde * (u_tilde*du_tx + v_tilde*du_ty + w_tilde*du_tz);
    Real un_tilde = U_tilde * inv_absS;  // face-normal velocity at Roe state

    d[4] = a_lam2 * (alpha2 * 0.5 * q2_tilde + alpha3_dot_u)
         + a_lam1 * alpha1 * (H_tilde - c_tilde * un_tilde)
         + a_lam5 * alpha5 * (H_tilde + c_tilde * un_tilde);

    // ---- 10. Final Roe flux: F = 0.5*(F_L + F_R) - 0.5*d ----

    for (int k = 0; k < 5; ++k) {
        flux[k] = 0.5 * (FL[k] + FR[k]) - 0.5 * d[k];
    }
}

// ============================================================================
// Directional solve methods
// ============================================================================

inline void RiemannSolverRoe::solve_xi(LocalBlock& lb, const Config& cfg) const {
    const auto& grid = lb.grid;
    auto& f = lb.field;

    const Int ni = grid.nci + 1;   // number of ξ-faces in i-direction
    const Int nj = grid.ncj;
    const Int nk = grid.nck;

    const Real gamma = cfg.gamma;
    const Real eps   = cfg.entropy_fix_eps;

    for (Int k = 0; k < nk; ++k) {
    for (Int j = 0; j < nj; ++j) {
    for (Int i = 0; i < ni; ++i) {
        Real qL[5] = { f.ql_xi.rho(i,j,k),  f.ql_xi.rhou(i,j,k),
                       f.ql_xi.rhov(i,j,k), f.ql_xi.rhow(i,j,k),
                       f.ql_xi.rhoE(i,j,k) };
        Real qR[5] = { f.qr_xi.rho(i,j,k),  f.qr_xi.rhou(i,j,k),
                       f.qr_xi.rhov(i,j,k), f.qr_xi.rhow(i,j,k),
                       f.qr_xi.rhoE(i,j,k) };

        Real flux[5];
        roe_flux(qL, qR,
                 grid.face_xi_x(i,j,k), grid.face_xi_y(i,j,k), grid.face_xi_z(i,j,k),
                 gamma, eps, flux);

        f.inv_xi.f1(i,j,k) = flux[0];
        f.inv_xi.f2(i,j,k) = flux[1];
        f.inv_xi.f3(i,j,k) = flux[2];
        f.inv_xi.f4(i,j,k) = flux[3];
        f.inv_xi.f5(i,j,k) = flux[4];
    }}}
}

inline void RiemannSolverRoe::solve_eta(LocalBlock& lb, const Config& cfg) const {
    const auto& grid = lb.grid;
    auto& f = lb.field;

    const Int ni = grid.nci;
    const Int nj = grid.ncj + 1;   // number of η-faces in j-direction
    const Int nk = grid.nck;

    const Real gamma = cfg.gamma;
    const Real eps   = cfg.entropy_fix_eps;

    for (Int k = 0; k < nk; ++k) {
    for (Int j = 0; j < nj; ++j) {
    for (Int i = 0; i < ni; ++i) {
        Real qL[5] = { f.ql_eta.rho(i,j,k),  f.ql_eta.rhou(i,j,k),
                       f.ql_eta.rhov(i,j,k), f.ql_eta.rhow(i,j,k),
                       f.ql_eta.rhoE(i,j,k) };
        Real qR[5] = { f.qr_eta.rho(i,j,k),  f.qr_eta.rhou(i,j,k),
                       f.qr_eta.rhov(i,j,k), f.qr_eta.rhow(i,j,k),
                       f.qr_eta.rhoE(i,j,k) };

        Real flux[5];
        roe_flux(qL, qR,
                 grid.face_eta_x(i,j,k), grid.face_eta_y(i,j,k), grid.face_eta_z(i,j,k),
                 gamma, eps, flux);

        f.inv_eta.f1(i,j,k) = flux[0];
        f.inv_eta.f2(i,j,k) = flux[1];
        f.inv_eta.f3(i,j,k) = flux[2];
        f.inv_eta.f4(i,j,k) = flux[3];
        f.inv_eta.f5(i,j,k) = flux[4];
    }}}
}

inline void RiemannSolverRoe::solve_zeta(LocalBlock& lb, const Config& cfg) const {
    const auto& grid = lb.grid;
    auto& f = lb.field;

    const Int ni = grid.nci;
    const Int nj = grid.ncj;
    const Int nk = grid.nck + 1;   // number of ζ-faces in k-direction

    const Real gamma = cfg.gamma;
    const Real eps   = cfg.entropy_fix_eps;

    for (Int k = 0; k < nk; ++k) {
    for (Int j = 0; j < nj; ++j) {
    for (Int i = 0; i < ni; ++i) {
        Real qL[5] = { f.ql_zeta.rho(i,j,k),  f.ql_zeta.rhou(i,j,k),
                       f.ql_zeta.rhov(i,j,k), f.ql_zeta.rhow(i,j,k),
                       f.ql_zeta.rhoE(i,j,k) };
        Real qR[5] = { f.qr_zeta.rho(i,j,k),  f.qr_zeta.rhou(i,j,k),
                       f.qr_zeta.rhov(i,j,k), f.qr_zeta.rhow(i,j,k),
                       f.qr_zeta.rhoE(i,j,k) };

        Real flux[5];
        roe_flux(qL, qR,
                 grid.face_zeta_x(i,j,k), grid.face_zeta_y(i,j,k), grid.face_zeta_z(i,j,k),
                 gamma, eps, flux);

        f.inv_zeta.f1(i,j,k) = flux[0];
        f.inv_zeta.f2(i,j,k) = flux[1];
        f.inv_zeta.f3(i,j,k) = flux[2];
        f.inv_zeta.f4(i,j,k) = flux[3];
        f.inv_zeta.f5(i,j,k) = flux[4];
    }}}
}

// ============================================================================
// Rusanov (Local Lax-Friedrichs) flux kernel
// ============================================================================

namespace detail {

inline void rusanov_flux(const Real qL[5], const Real qR[5],
                          Real Sx, Real Sy, Real Sz,
                          Real gamma, Real /*eps*/,
                          Real flux[5])
{
    // ---- 1. Extract primitive variables ----

    Real rhoL = qL[0];
    Real uL   = qL[1] / rhoL;
    Real vL   = qL[2] / rhoL;
    Real wL   = qL[3] / rhoL;
    Real keL  = 0.5 * (qL[1]*uL + qL[2]*vL + qL[3]*wL);
    Real pL   = (gamma - 1.0) * (qL[4] - keL);

    Real rhoR = qR[0];
    Real uR   = qR[1] / rhoR;
    Real vR   = qR[2] / rhoR;
    Real wR   = qR[3] / rhoR;
    Real keR  = 0.5 * (qR[1]*uR + qR[2]*vR + qR[3]*wR);
    Real pR   = (gamma - 1.0) * (qR[4] - keR);

    // ---- 2. Physical flux through the face ----

    Real UL = Sx*uL + Sy*vL + Sz*wL;
    Real UR = Sx*uR + Sy*vR + Sz*wR;

    Real FL[5], FR[5];
    FL[0] = qL[0] * UL;
    FL[1] = qL[1] * UL + pL * Sx;
    FL[2] = qL[2] * UL + pL * Sy;
    FL[3] = qL[3] * UL + pL * Sz;
    FL[4] = (qL[4] + pL) * UL;

    FR[0] = qR[0] * UR;
    FR[1] = qR[1] * UR + pR * Sx;
    FR[2] = qR[2] * UR + pR * Sy;
    FR[3] = qR[3] * UR + pR * Sz;
    FR[4] = (qR[4] + pR) * UR;

    // ---- 3. Maximum local wave speed (spectral radius) ----

    Real absS  = std::sqrt(Sx*Sx + Sy*Sy + Sz*Sz);
    Real cL    = std::sqrt(std::max(0.0, gamma * pL / rhoL));
    Real cR    = std::sqrt(std::max(0.0, gamma * pR / rhoR));
    Real SmaxL = std::abs(UL) + cL * absS;
    Real SmaxR = std::abs(UR) + cR * absS;
    Real Smax  = std::max(SmaxL, SmaxR);

    // ---- 4. Rusanov flux: F = 0.5*(F_L+F_R) - 0.5*Smax*(Q_R-Q_L) ----

    for (int k = 0; k < 5; ++k) {
        flux[k] = 0.5 * (FL[k] + FR[k]) - 0.5 * Smax * (qR[k] - qL[k]);
    }
}

} // namespace detail

// --- Directional solve methods ---

inline void RiemannSolverRusanov::solve_xi  (LocalBlock& lb, const Config& cfg) const { detail::solve_xi_impl  (lb, cfg.gamma, cfg.entropy_fix_eps, detail::rusanov_flux); }
inline void RiemannSolverRusanov::solve_eta (LocalBlock& lb, const Config& cfg) const { detail::solve_eta_impl (lb, cfg.gamma, cfg.entropy_fix_eps, detail::rusanov_flux); }
inline void RiemannSolverRusanov::solve_zeta(LocalBlock& lb, const Config& cfg) const { detail::solve_zeta_impl(lb, cfg.gamma, cfg.entropy_fix_eps, detail::rusanov_flux); }

// ============================================================================
// HLL flux kernel
// ============================================================================

namespace detail {

inline void hll_flux(const Real qL[5], const Real qR[5],
                      Real Sx, Real Sy, Real Sz,
                      Real gamma, Real /*eps*/,
                      Real flux[5])
{
    // ---- 1. Extract primitive variables ----

    Real rhoL = qL[0];
    Real uL   = qL[1] / rhoL;
    Real vL   = qL[2] / rhoL;
    Real wL   = qL[3] / rhoL;
    Real keL  = 0.5 * (qL[1]*uL + qL[2]*vL + qL[3]*wL);
    Real pL   = (gamma - 1.0) * (qL[4] - keL);

    Real rhoR = qR[0];
    Real uR   = qR[1] / rhoR;
    Real vR   = qR[2] / rhoR;
    Real wR   = qR[3] / rhoR;
    Real keR  = 0.5 * (qR[1]*uR + qR[2]*vR + qR[3]*wR);
    Real pR   = (gamma - 1.0) * (qR[4] - keR);

    // ---- 2. Physical flux through the face ----

    Real UL = Sx*uL + Sy*vL + Sz*wL;
    Real UR = Sx*uR + Sy*vR + Sz*wR;

    Real FL[5], FR[5];
    FL[0] = qL[0] * UL;
    FL[1] = qL[1] * UL + pL * Sx;
    FL[2] = qL[2] * UL + pL * Sy;
    FL[3] = qL[3] * UL + pL * Sz;
    FL[4] = (qL[4] + pL) * UL;

    FR[0] = qR[0] * UR;
    FR[1] = qR[1] * UR + pR * Sx;
    FR[2] = qR[2] * UR + pR * Sy;
    FR[3] = qR[3] * UR + pR * Sz;
    FR[4] = (qR[4] + pR) * UR;

    // ---- 3. Wave speed estimates (Davis) ----

    Real absS = std::sqrt(Sx*Sx + Sy*Sy + Sz*Sz);
    Real cL   = std::sqrt(std::max(0.0, gamma * pL / rhoL));
    Real cR   = std::sqrt(std::max(0.0, gamma * pR / rhoR));

    Real SL = std::min(UL - cL * absS, UR - cR * absS);
    Real SR = std::max(UL + cL * absS, UR + cR * absS);

    // ---- 4. HLL flux ----

    if (SL >= 0.0) {
        for (int k = 0; k < 5; ++k) flux[k] = FL[k];
    } else if (SR <= 0.0) {
        for (int k = 0; k < 5; ++k) flux[k] = FR[k];
    } else {
        Real inv_dS = 1.0 / (SR - SL);
        for (int k = 0; k < 5; ++k) {
            flux[k] = (SR * FL[k] - SL * FR[k] + SL * SR * (qR[k] - qL[k])) * inv_dS;
        }
    }
}

} // namespace detail

// --- Directional solve methods ---

inline void RiemannSolverHLL::solve_xi  (LocalBlock& lb, const Config& cfg) const { detail::solve_xi_impl  (lb, cfg.gamma, cfg.entropy_fix_eps, detail::hll_flux); }
inline void RiemannSolverHLL::solve_eta (LocalBlock& lb, const Config& cfg) const { detail::solve_eta_impl (lb, cfg.gamma, cfg.entropy_fix_eps, detail::hll_flux); }
inline void RiemannSolverHLL::solve_zeta(LocalBlock& lb, const Config& cfg) const { detail::solve_zeta_impl(lb, cfg.gamma, cfg.entropy_fix_eps, detail::hll_flux); }

// ============================================================================
// HLLC flux kernel
// ============================================================================

namespace detail {

inline void hllc_flux(const Real qL[5], const Real qR[5],
                       Real Sx, Real Sy, Real Sz,
                       Real gamma, Real /*eps*/,
                       Real flux[5])
{
    // ---- 1. Extract primitive variables ----

    Real rhoL = qL[0];
    Real uL   = qL[1] / rhoL;
    Real vL   = qL[2] / rhoL;
    Real wL   = qL[3] / rhoL;
    Real keL  = 0.5 * (qL[1]*uL + qL[2]*vL + qL[3]*wL);
    Real pL   = (gamma - 1.0) * (qL[4] - keL);

    Real rhoR = qR[0];
    Real uR   = qR[1] / rhoR;
    Real vR   = qR[2] / rhoR;
    Real wR   = qR[3] / rhoR;
    Real keR  = 0.5 * (qR[1]*uR + qR[2]*vR + qR[3]*wR);
    Real pR   = (gamma - 1.0) * (qR[4] - keR);

    // ---- 2. Physical flux through the face ----

    Real UL = Sx*uL + Sy*vL + Sz*wL;
    Real UR = Sx*uR + Sy*vR + Sz*wR;

    Real FL[5], FR[5];
    FL[0] = qL[0] * UL;
    FL[1] = qL[1] * UL + pL * Sx;
    FL[2] = qL[2] * UL + pL * Sy;
    FL[3] = qL[3] * UL + pL * Sz;
    FL[4] = (qL[4] + pL) * UL;

    FR[0] = qR[0] * UR;
    FR[1] = qR[1] * UR + pR * Sx;
    FR[2] = qR[2] * UR + pR * Sy;
    FR[3] = qR[3] * UR + pR * Sz;
    FR[4] = (qR[4] + pR) * UR;

    // ---- 3. Wave speed estimates (Davis) ----

    Real absS     = std::sqrt(Sx*Sx + Sy*Sy + Sz*Sz);
    Real inv_absS = 1.0 / std::max(absS, Real(1e-30));
    Real inv_absS2 = inv_absS * inv_absS;

    Real cL = std::sqrt(std::max(0.0, gamma * pL / rhoL));
    Real cR = std::sqrt(std::max(0.0, gamma * pR / rhoR));

    Real SL = std::min(UL - cL * absS, UR - cR * absS);
    Real SR = std::max(UL + cL * absS, UR + cR * absS);

    // ---- 4. HLLC flux branching ----

    if (SL >= 0.0) {
        for (int k = 0; k < 5; ++k) flux[k] = FL[k];
        return;
    }
    if (SR <= 0.0) {
        for (int k = 0; k < 5; ++k) flux[k] = FR[k];
        return;
    }

    // ---- 5. Contact wave speed S* ----
    // Solve p*_L = p*_R for S* (in contravariant units, includes |S| factor)

    Real dL = rhoL * (SL - UL);
    Real dR = rhoR * (SR - UR);

    Real denom = dL - dR;
    if (std::abs(denom) < Real(1e-30)) {
        // Degenerate — fall back to HLL
        Real inv_dS = 1.0 / (SR - SL);
        for (int k = 0; k < 5; ++k) {
            flux[k] = (SR * FL[k] - SL * FR[k] + SL * SR * (qR[k] - qL[k])) * inv_dS;
        }
        return;
    }

    // S* = [ (p_R-p_L)*|S|² + dL*U_L - dR*U_R ] / (dL - dR)
    Real S_star = ((pR - pL) * absS * absS + dL * UL - dR * UR) / denom;

    // ---- 5a. Robustness: clamp S* to [SL, SR] ----
    S_star = std::max(SL, std::min(SR, S_star));

    // ---- 6. Star-state flux ----

    auto fallback_hll = [&]() {
        Real inv_dS = 1.0 / (SR - SL);
        for (int k = 0; k < 5; ++k) {
            flux[k] = (SR * FL[k] - SL * FR[k] + SL * SR * (qR[k] - qL[k])) * inv_dS;
        }
    };

    if (S_star >= 0.0) {
        // Left star region (S_L < 0 =< S_star)
        Real rho_star = dL / (SL - S_star);
        Real p_star   = pL + dL * (S_star - UL) * inv_absS2;

        // Fall back to HLL if star state is non-physical
        if (rho_star <= 0.0 || p_star <= 0.0) {
            fallback_hll();
            return;
        }

        Real u_star = uL + (S_star - UL) * Sx * inv_absS2;
        Real v_star = vL + (S_star - UL) * Sy * inv_absS2;
        Real w_star = wL + (S_star - UL) * Sz * inv_absS2;

        Real ke_star = 0.5 * (u_star*u_star + v_star*v_star + w_star*w_star);
        Real E_star  = p_star / (rho_star * (gamma - 1.0)) + ke_star;

        flux[0] = rho_star * S_star;
        flux[1] = rho_star * u_star * S_star + p_star * Sx;
        flux[2] = rho_star * v_star * S_star + p_star * Sy;
        flux[3] = rho_star * w_star * S_star + p_star * Sz;
        flux[4] = (rho_star * E_star + p_star) * S_star;
    } else {
        // Right star region (S_star < 0 < S_R)
        Real rho_star = dR / (SR - S_star);
        Real p_star   = pR + dR * (S_star - UR) * inv_absS2;

        // Fall back to HLL if star state is non-physical
        if (rho_star <= 0.0 || p_star <= 0.0) {
            fallback_hll();
            return;
        }

        Real u_star = uR + (S_star - UR) * Sx * inv_absS2;
        Real v_star = vR + (S_star - UR) * Sy * inv_absS2;
        Real w_star = wR + (S_star - UR) * Sz * inv_absS2;

        Real ke_star = 0.5 * (u_star*u_star + v_star*v_star + w_star*w_star);
        Real E_star  = p_star / (rho_star * (gamma - 1.0)) + ke_star;

        flux[0] = rho_star * S_star;
        flux[1] = rho_star * u_star * S_star + p_star * Sx;
        flux[2] = rho_star * v_star * S_star + p_star * Sy;
        flux[3] = rho_star * w_star * S_star + p_star * Sz;
        flux[4] = (rho_star * E_star + p_star) * S_star;
    }
}

} // namespace detail

// --- Directional solve methods ---

inline void RiemannSolverHLLC::solve_xi  (LocalBlock& lb, const Config& cfg) const { detail::solve_xi_impl  (lb, cfg.gamma, cfg.entropy_fix_eps, detail::hllc_flux); }
inline void RiemannSolverHLLC::solve_eta (LocalBlock& lb, const Config& cfg) const { detail::solve_eta_impl (lb, cfg.gamma, cfg.entropy_fix_eps, detail::hllc_flux); }
inline void RiemannSolverHLLC::solve_zeta(LocalBlock& lb, const Config& cfg) const { detail::solve_zeta_impl(lb, cfg.gamma, cfg.entropy_fix_eps, detail::hllc_flux); }
