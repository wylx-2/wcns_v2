#pragma once

#include "wcns_v2/utils/types.h"
#include "wcns_v2/core/config.h"
#include <memory>

// Forward declarations
class LocalBlock;

/// @file riemann_solver.h
/// @brief Riemann solvers for computing inviscid numerical fluxes at cell faces.
///
/// Each solver reads the left/right states (Q_L, Q_R) interpolated at half-nodes
/// and the face metric coefficients (area vectors), then computes the numerical
/// inviscid flux through each face.  Results are stored in Field::inv_xi/eta/zeta.
///
/// Face flux convention:
///   F_face = Sx*F(Q) + Sy*G(Q) + Sz*H(Q)
/// where (Sx,Sy,Sz) is the face area vector (not unit normal).
///
/// Supported schemes (controlled by cfg.riemann_type):
///   - "roe"  — Roe approximate Riemann solver (with entropy fix)


/// Abstract base for Riemann solvers.
/// Derived classes implement specific flux schemes.
class RiemannSolverBase {
public:
    virtual ~RiemannSolverBase() = default;

    /// Compute inviscid flux at ξ-faces (i+1/2).
    /// Reads ql_xi/qr_xi and face_xi_* from lb; writes inv_xi to lb.field.
    virtual void solve_xi(LocalBlock& lb, const Config& cfg) const = 0;

    /// Compute inviscid flux at η-faces (j+1/2).
    /// Reads ql_eta/qr_eta and face_eta_* from lb; writes inv_eta to lb.field.
    virtual void solve_eta(LocalBlock& lb, const Config& cfg) const = 0;

    /// Compute inviscid flux at ζ-faces (k+1/2).
    /// Reads ql_zeta/qr_zeta and face_zeta_* from lb; writes inv_zeta to lb.field.
    virtual void solve_zeta(LocalBlock& lb, const Config& cfg) const = 0;

    /// Factory: create the solver specified by cfg.riemann_type.
    /// Supported: "roe".
    static std::unique_ptr<RiemannSolverBase> create(const Config& cfg);
};


/// Roe approximate Riemann solver.
///
/// For each face, left and right conservative states are combined via
/// Roe-averaging to produce the numerical flux:
///
///   F_roe = 0.5*(F_L + F_R) - 0.5 * |A| * (Q_R - Q_L)
///
/// where |A| = R·|Λ|·R^{-1} is the absolute flux Jacobian in face-normal
/// direction evaluated at the Roe-average state.
///
/// Entropy fix (Harten correction) is applied to eigenvalues near zero,
/// controlled by cfg.entropy_fix_eps.
class RiemannSolverRoe : public RiemannSolverBase {
public:
    void solve_xi(LocalBlock& lb, const Config& cfg) const override;
    void solve_eta(LocalBlock& lb, const Config& cfg) const override;
    void solve_zeta(LocalBlock& lb, const Config& cfg) const override;

private:
    /// Core Roe flux for one face.
    ///
    /// @param qL    Left conservative state  [rho, rhou, rhov, rhow, rhoE]
    /// @param qR    Right conservative state [rho, rhou, rhov, rhow, rhoE]
    /// @param Sx,Sy,Sz  Face area vector (face metric at this half-node)
    /// @param gamma Specific heat ratio
    /// @param eps   Entropy fix coefficient (0 = no fix)
    /// @param flux  [out] Numerical flux [F1..F5]
    static void roe_flux(const Real qL[5], const Real qR[5],
                         Real Sx, Real Sy, Real Sz,
                         Real gamma, Real eps,
                         Real flux[5]);
};

#include "riemann_solver.hxx"
