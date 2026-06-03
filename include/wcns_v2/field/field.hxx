#pragma once

#include "field.h"
#include <cmath>
#include <stdexcept>

// ============================================================================
// PrimitiveVars
// ============================================================================

inline void PrimitiveVars::allocate(Int ni, Int nj, Int nk) {
    rho.allocate(ni, nj, nk);
    u.allocate(ni, nj, nk);
    v.allocate(ni, nj, nk);
    w.allocate(ni, nj, nk);
    p.allocate(ni, nj, nk);
}

inline void PrimitiveVars::fill(Real val) {
    rho.fill(val);
    u.fill(val);
    v.fill(val);
    w.fill(val);
    p.fill(val);
}

// ============================================================================
// ConservativeVars
// ============================================================================

inline void ConservativeVars::allocate(Int ni, Int nj, Int nk) {
    rho.allocate(ni, nj, nk);
    rhou.allocate(ni, nj, nk);
    rhov.allocate(ni, nj, nk);
    rhow.allocate(ni, nj, nk);
    rhoE.allocate(ni, nj, nk);
}

inline void ConservativeVars::fill(Real val) {
    rho.fill(val);
    rhou.fill(val);
    rhov.fill(val);
    rhow.fill(val);
    rhoE.fill(val);
}

// ============================================================================
// FluxVars
// ============================================================================

inline void FluxVars::allocate(Int ni, Int nj, Int nk) {
    f1.allocate(ni, nj, nk);
    f2.allocate(ni, nj, nk);
    f3.allocate(ni, nj, nk);
    f4.allocate(ni, nj, nk);
    f5.allocate(ni, nj, nk);
}

inline void FluxVars::fill(Real val) {
    f1.fill(val);
    f2.fill(val);
    f3.fill(val);
    f4.fill(val);
    f5.fill(val);
}

// ============================================================================
// Field — allocation
// ============================================================================

inline void Field::allocate(Int nci, Int ncj, Int nck) {
    nci_ = nci;
    ncj_ = ncj;
    nck_ = nck;

    // Core state — cell-centered
    prim.allocate(nci, ncj, nck);
    cons.allocate(nci, ncj, nck);
    rhs.allocate(nci, ncj, nck);

    // Face fluxes — size depends on face direction
    // ξ-face (i+1/2): (nci+1) × ncj × nck
    inv_xi.allocate(nci + 1, ncj, nck);
    vis_xi.allocate(nci + 1, ncj, nck);

    // Interpolated left/right states at ξ-face (i+1/2)
    ql_xi.allocate(nci + 1, ncj, nck);
    qr_xi.allocate(nci + 1, ncj, nck);

    // η-face (j+1/2): nci × (ncj+1) × nck
    inv_eta.allocate(nci, ncj + 1, nck);
    vis_eta.allocate(nci, ncj + 1, nck);

    // Interpolated left/right states at η-face (j+1/2)
    ql_eta.allocate(nci, ncj + 1, nck);
    qr_eta.allocate(nci, ncj + 1, nck);

    // ζ-face (k+1/2): nci × ncj × (nck+1)
    inv_zeta.allocate(nci, ncj, nck + 1);
    vis_zeta.allocate(nci, ncj, nck + 1);

    // Interpolated left/right states at ζ-face (k+1/2)
    ql_zeta.allocate(nci, ncj, nck + 1);
    qr_zeta.allocate(nci, ncj, nck + 1);

    allocated_ = true;
}

// ============================================================================
// Field — primitive ↔ conservative conversion
// ============================================================================

inline void Field::prim_to_cons(Real gamma) {
    Real gm1 = gamma - 1.0;

    for (Int k = 0; k < nck_; ++k) {
    for (Int j = 0; j < ncj_; ++j) {
    for (Int i = 0; i < nci_; ++i) {
        Real rho = prim.rho(i,j,k);
        Real u   = prim.u(i,j,k);
        Real v   = prim.v(i,j,k);
        Real w   = prim.w(i,j,k);
        Real p   = prim.p(i,j,k);

        Real ke = 0.5 * rho * (u*u + v*v + w*w);

        cons.rho(i,j,k)  = rho;
        cons.rhou(i,j,k) = rho * u;
        cons.rhov(i,j,k) = rho * v;
        cons.rhow(i,j,k) = rho * w;
        cons.rhoE(i,j,k) = p / gm1 + ke;
    }}}
}

inline void Field::cons_to_prim(Real gamma) {
    Real gm1 = gamma - 1.0;

    for (Int k = 0; k < nck_; ++k) {
    for (Int j = 0; j < ncj_; ++j) {
    for (Int i = 0; i < nci_; ++i) {
        Real rho  = cons.rho(i,j,k);
        Real rhou = cons.rhou(i,j,k);
        Real rhov = cons.rhov(i,j,k);
        Real rhow = cons.rhow(i,j,k);
        Real rhoE = cons.rhoE(i,j,k);

        Real inv_rho = 1.0 / rho;
        Real ke = 0.5 * inv_rho * (rhou*rhou + rhov*rhov + rhow*rhow);

        prim.rho(i,j,k) = rho;
        prim.u(i,j,k)   = rhou * inv_rho;
        prim.v(i,j,k)   = rhov * inv_rho;
        prim.w(i,j,k)   = rhow * inv_rho;
        prim.p(i,j,k)   = gm1 * (rhoE - ke);
    }}}
}

// ============================================================================
// Field — extension fields
// ============================================================================

inline void Field::add_extension(const std::string& name) {
    if (ext_fields_.find(name) == ext_fields_.end()) {
        ext_fields_[name].allocate(nci_, ncj_, nck_);
    }
}

inline bool Field::has_extension(const std::string& name) const {
    return ext_fields_.find(name) != ext_fields_.end();
}

inline MultiArray3D<Real>& Field::ext(const std::string& name) {
    auto it = ext_fields_.find(name);
    if (it == ext_fields_.end()) {
        throw std::runtime_error("Field: extension \"" + name + "\" not found");
    }
    return it->second;
}

inline const MultiArray3D<Real>& Field::ext(const std::string& name) const {
    auto it = ext_fields_.find(name);
    if (it == ext_fields_.end()) {
        throw std::runtime_error("Field: extension \"" + name + "\" not found");
    }
    return it->second;
}

// ============================================================================
// Field — query
// ============================================================================

inline Int Field::ni() const { return nci_; }
inline Int Field::nj() const { return ncj_; }
inline Int Field::nk() const { return nck_; }
