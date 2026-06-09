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

    // ---- Periodic properties ----
    bool is_periodic;
    Real translation[3];       // periodic translation vector
    Real rotation_center[3];   // rotation centre for periodic transform
    Real rotation_angle[3];    // rotation angle for periodic transform
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
