#pragma once

#include "utils/types.h"
#include "core/config.h"
#include <memory>

// Forward declarations
class LocalBlock;

/// @file wcns_interp.h
/// @brief WCNS interpolation — cell-center values to half-node (face) left/right states.
///
/// Abstract base WcnsInterpBase defines the interface for all interpolation methods.
/// Derived classes implement specific schemes (WENO-JS, MDCD linear, MDCD hybrid).
///
/// Convention: half-node ah[i] = a_{i-1/2} (half-node between cell i-1 and cell i).
/// For n cells there are n+1 half-nodes indexed 0..n.
///
/// Face array sizes:
///   ξ-direction (i-face): (nci+1) × ncj × nck
///   η-direction (j-face): nci × (ncj+1) × nck
///   ζ-direction (k-face): nci × ncj × (nck+1)


/// Abstract base for WCNS interpolation schemes.
/// Derived classes implement specific methods (nonlinear WENO, MDCD linear, MDCD hybrid).
class WcnsInterpBase {
public:
    virtual ~WcnsInterpBase() = default;

    /// Interpolate cell-center conservative vars to ξ-faces (i+1/2).
    /// Results stored in lb.field.ql_xi / lb.field.qr_xi.
    virtual void interp_xi(LocalBlock& lb, const Config& cfg) const = 0;

    /// Interpolate cell-center conservative vars to η-faces (j+1/2).
    /// Results stored in lb.field.ql_eta / lb.field.qr_eta.
    virtual void interp_eta(LocalBlock& lb, const Config& cfg) const = 0;

    /// Interpolate cell-center conservative vars to ζ-faces (k+1/2).
    /// Results stored in lb.field.ql_zeta / lb.field.qr_zeta.
    virtual void interp_zeta(LocalBlock& lb, const Config& cfg) const = 0;

    /// Factory: create the interpolator specified by cfg.interp_type.
    /// Supported types: "weno_js", "mdcd_linear", "mdcd_hybrid".
    static std::unique_ptr<WcnsInterpBase> create(const Config& cfg);

protected:
    // ========================================================================
    // Shared helpers
    // ========================================================================

    /// Fill boundary half-nodes (h = 0,1,2 and h = n-2,n-1,n) using one-sided
    /// 5th-order interpolation from InterpDiff.  Both Q_L and Q_R receive the
    /// same value to ensure the Riemann solver degenerates correctly at boundaries.
    ///
    /// @param a   Cell-center array (starting at index 0)
    /// @param h   Half-node index [0..n]
    /// @param n   Number of cell centers in this line
    /// @param ql  [out] Left state at half-node h
    /// @param qr  [out] Right state at half-node h
    /// @return    true if h is a boundary half-node (ql,qr have been filled)
    static bool fill_boundary(const Real* a, Int h, Int n, Real& ql, Real& qr);

    /// Check that interp_vars is supported by this interpolator.
    static void check_interp_vars(const Config& cfg, const char* scheme_name);
};


/// Original WCNS nonlinear interpolation (5th-order WENO-JS).
///
/// For each half-node i+1/2:
///   - Q_L (left-biased):  3 sub-stencils from [i-2 .. i+2], WENO-combined
///   - Q_R (right-biased): 3 sub-stencils from [i-1 .. i+3], WENO-combined
///
/// At boundaries (h = 0,1,2 and h = n-2,n-1,n) where the full 5-point
/// stencil is unavailable, one-sided 5th-order interpolation from InterpDiff
/// is used and the same value is assigned to both Q_L and Q_R.
class WcnsNonlinearInterp : public WcnsInterpBase {
public:
    void interp_xi(LocalBlock& lb, const Config& cfg) const override;
    void interp_eta(LocalBlock& lb, const Config& cfg) const override;
    void interp_zeta(LocalBlock& lb, const Config& cfg) const override;

    /// Run interpolation for all three directions at once (convenience).
    static void interp_all(LocalBlock& lb, const Config& cfg);

private:
    /// 1D WENO-JS interpolation kernel.
    static void interp_1d(const Real* a, Real* ql, Real* qr, Int n);

    /// WENO left-biased interpolation to half-node i+1/2.
    /// Uses stencil a[i-2 .. i+2] (5 cells).
    static Real weno_left(const Real* a, Int i);

    /// WENO right-biased interpolation to half-node i+1/2.
    /// Uses stencil a[i-1 .. i+3] (5 cells).
    static Real weno_right(const Real* a, Int i);
};


/// MDCD linear interpolation — 6-point stencil with tunable dissipation/dispersion.
///
/// For each half-node i+1/2:
///   - Q_L uses 6-point left-biased stencil [i-2 .. i+3]
///   - Q_R uses mirror-reversed stencil (equivalent to right-biased)
///
/// Parameters (from Config):
///   - mdcd_diss : dissipation coefficient (default 0)
///   - mdcd_disp : dispersion coefficient (default 0)
///
/// When diss = disp = 0, the scheme reduces to the standard 6th-order linear
/// interpolation.
class WcnsMdcdLinearInterp : public WcnsInterpBase {
    friend class WcnsMdcdHybridInterp;  // access mdcd_left for fallback
public:
    void interp_xi(LocalBlock& lb, const Config& cfg) const override;
    void interp_eta(LocalBlock& lb, const Config& cfg) const override;
    void interp_zeta(LocalBlock& lb, const Config& cfg) const override;

private:
    /// 1D MDCD linear interpolation kernel.
    static void interp_1d(const Real* a, Real* ql, Real* qr, Int n,
                          Real diss, Real disp);

    /// 6-point left-biased MDCD linear interpolation to i+1/2.
    /// stencil[0..5] = {a[i-2], a[i-1], a[i], a[i+1], a[i+2], a[i+3]}.
    static Real mdcd_left(const Real* stencil, Real diss, Real disp);
};


/// MDCD-WENO hybrid interpolation — 6-point stencil with discontinuity detector.
///
/// At each half-node, a discontinuity detector (sai) is evaluated:
///   - If discontinuous (sai > sai_ref): falls back to MDCD linear interpolation
///     (dissipative to suppress oscillations)
///   - If smooth: uses WENO-Z style nonlinear reconstruction with 4 sub-stencils
///     and 6th-order global smoothness indicator τ6
///
/// Parameters (from Config):
///   - mdcd_diss, mdcd_disp : MDCD coefficients for optimal weights & fallback
///   - mdcd_sai_ref : discontinuity detector threshold (default 0.4)
class WcnsMdcdHybridInterp : public WcnsInterpBase {
public:
    void interp_xi(LocalBlock& lb, const Config& cfg) const override;
    void interp_eta(LocalBlock& lb, const Config& cfg) const override;
    void interp_zeta(LocalBlock& lb, const Config& cfg) const override;

private:
    /// 1D MDCD-WENO hybrid interpolation kernel.
    static void interp_1d(const Real* a, Real* ql, Real* qr, Int n,
                          Real diss, Real disp, Real sai_ref);

    /// 6-point left-biased MDCD-WENO hybrid interpolation to i+1/2.
    /// stencil[0..5] = {a[i-2], a[i-1], a[i], a[i+1], a[i+2], a[i+3]}.
    static Real mdcd_hybrid_left(const Real* stencil, Real diss, Real disp,
                                 Real sai_ref);
};


/// Convenience: run interpolation for all 3 directions using the scheme
/// specified by cfg.interp_type.
inline void wcns_interp_all(LocalBlock& lb, const Config& cfg) {
    auto interp = WcnsInterpBase::create(cfg);
    interp->interp_xi(lb, cfg);
    interp->interp_eta(lb, cfg);
    interp->interp_zeta(lb, cfg);
}

#include "wcns_interp.hxx"
