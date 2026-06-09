#pragma once

#include "boundary_condition.h"

// CGNS BC type constants (matching cgnslib.h BCType_t enum)
namespace {
    constexpr int CG_BCInflow            =  9;
    constexpr int CG_BCInflowSubsonic    = 10;
    constexpr int CG_BCInflowSupersonic  = 11;
    constexpr int CG_BCOutflow           = 13;
    constexpr int CG_BCOutflowSubsonic   = 14;
    constexpr int CG_BCOutflowSupersonic = 15;
    constexpr int CG_BCTunnelInflow      = 16;
    constexpr int CG_BCTunnelOutflow     = 17;
    constexpr int CG_BCWall              = 20;
    constexpr int CG_BCWallInviscid      = 21;
    constexpr int CG_BCWallViscous       = 22;
    constexpr int CG_BCWallViscousHeatFlux = 23;
    constexpr int CG_BCWallViscousIsothermal = 24;
    constexpr int CG_BCSymmetryPlane     = 25;
    constexpr int CG_BCFarfield          = 26;
    constexpr int CG_BCFamilySpecified   = 38;
}

inline BoundaryCondition::BoundaryCondition() {}

inline void BoundaryCondition::add_patch(const BCPatch& patch) {
    patches_.push_back(patch);
}

inline Int BoundaryCondition::num_patches() const {
    return static_cast<Int>(patches_.size());
}

inline const BCPatch& BoundaryCondition::patch(Int i) const {
    return patches_[static_cast<std::size_t>(i)];
}

inline BCPatch& BoundaryCondition::patch(Int i) {
    return patches_[static_cast<std::size_t>(i)];
}

inline BCType BoundaryCondition::from_cgns(int cgns_bc_type) {
    switch (cgns_bc_type) {
    case CG_BCInflow:
    case CG_BCInflowSubsonic:
    case CG_BCInflowSupersonic:
    case CG_BCTunnelInflow:        return BCType::Inflow;
    case CG_BCOutflow:
    case CG_BCOutflowSubsonic:
    case CG_BCOutflowSupersonic:
    case CG_BCTunnelOutflow:       return BCType::Outflow;
    case CG_BCWall:
    case CG_BCWallInviscid:
    case CG_BCWallViscous:
    case CG_BCWallViscousHeatFlux:
    case CG_BCWallViscousIsothermal: return BCType::Wall;
    case CG_BCSymmetryPlane:       return BCType::Symmetry;
    case CG_BCFarfield:            return BCType::Farfield;
    case CG_BCFamilySpecified:     return BCType::Unknown;
    default:                       return BCType::Unknown;
    }
}
