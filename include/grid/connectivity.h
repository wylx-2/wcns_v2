#pragma once

#include "utils/types.h"
#include <string>
#include <vector>

/// A single 1-to-1 grid connectivity (block interface / periodic boundary).
///
/// In CGNS, structured 1-to-1 connections map a face range on one zone
/// (the "current" zone) to a face range on another zone (or the same
/// zone for periodic BCs).  All ranges use node indices (1-based,
/// inclusive), matching the CGNS convention.
struct Connectivity {
    std::string name;        // connection name from CGNS
    std::string donor_name;  // donor zone name

    // Range on the CURRENT zone (1-based, inclusive node indices)
    Int imin, imax, jmin, jmax, kmin, kmax;

    // Range on the DONOR zone (1-based, inclusive node indices)
    Int donor_imin, donor_imax, donor_jmin, donor_jmax, donor_kmin, donor_kmax;

    // Transform: how the current zone's (i,j,k) directions map to the
    // donor zone's index directions.  E.g. {1,2,3} = i→i, j→j, k→k.
    // A negative value means the direction is reversed.
    Int transform[3];

    // ---- Face orientation ----
    /// Which face of the current zone this connection lies on.
    /// 0=IMIN, 1=IMAX, 2=JMIN, 3=JMAX, 4=KMIN, 5=KMAX.  -1 = undetected.
    int face = -1;

    // ---- Periodic / transform properties ----
    /// For both periodic and general 1-to-1 connections, these describe
    /// the spatial relationship between the current zone and the donor.
    /// Ghost coordinates are filled as: ghost = donor_interior - translation.
    bool is_periodic;
    Real translation[3];       // translation vector
    Real rotation_center[3];   // rotation centre
    Real rotation_angle[3];    // rotation angle
};

/// Collection of 1-to-1 connectivity entries for a block.
///
/// Holds both block-interface connections (multi-block) and periodic
/// connections (single-block wrap-around or multi-block periodicity).
class ConnectivityList {
public:
    ConnectivityList() = default;

    void add(const Connectivity& conn);

    Int  count() const;
    const Connectivity& operator[](Int i) const;
    Connectivity&       operator[](Int i);

private:
    std::vector<Connectivity> connections_;
};

#include "connectivity.hxx"
