#pragma once

#include "wcns_v2/utils/types.h"
#include <string>
#include <vector>

/// Boundary condition type enumeration.
/// CGNS BC types are mapped to these canonical types.
enum class BCType {
    Inflow,
    Outflow,
    Wall,
    Farfield,
    Symmetry,
    Periodic,
    Unknown
};

/// A single boundary-condition patch.
/// Describes a subset of faces on one side of the structured block.
struct BCPatch {
    std::string name;       // from CGNS BC name
    BCType      type;

    // Face: 0=IMIN, 1=IMAX, 2=JMIN, 3=JMAX, 4=KMIN, 5=KMAX
    int  face;

    // Range in NODE indices (CGNS convention: 1-based, inclusive)
    Int  imin, imax, jmin, jmax, kmin, kmax;
};

/// Collection of boundary-condition patches for a block.
class BoundaryCondition {
public:
    BoundaryCondition();

    void add_patch(const BCPatch& patch);

    Int  num_patches() const;
    const BCPatch& patch(Int i) const;
    BCPatch&       patch(Int i);

    /// Map a CGNS BC type enum value to our BCType.
    static BCType from_cgns(int cgns_bc_type);

private:
    std::vector<BCPatch> patches_;
};

#include "boundary_condition.hxx"
