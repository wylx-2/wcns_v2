#pragma once

#include "utils/types.h"
#include "grid/grid.h"
#include <string>

/// Reader for CGNS structured-grid files.
///
/// Reads a single-zone or multi-zone CGNS file.  Only structured zones
/// (ZoneType == Structured) are supported.  Node coordinates and boundary
/// condition patches are extracted and stored into a Grid object.
///
/// Usage:
///   CGNSReader reader;
///   reader.open("grid.cgns");
///   reader.read_zone(1, grid);      // read first zone
///   // or  reader.read_all(grids);  // read all zones
///   reader.close();
class CGNSReader {
public:
    CGNSReader();
    ~CGNSReader();

    /// Open a CGNS file for reading.
    void open(const std::string& filename);

    /// Close the file.
    void close();

    /// Number of bases in the file.
    Int num_bases() const;

    /// Number of zones in the given base (1-indexed).
    Int num_zones(Int base = 1) const;

    /// Read a single structured zone into a Grid object.
    /// @param base  1-indexed base index
    /// @param zone  1-indexed zone index
    /// @param grid  Output grid (allocated and filled)
    /// @param ng    Number of ghost layers to allocate
    void read_zone(Int base, Int zone, Grid& grid, Int ng = 2);

    /// Read all zones from base 1 into a vector of Grid objects.
    void read_all(std::vector<Grid>& grids, Int ng = 2);

private:
    int file_id_;
    bool is_open_;
    Int num_bases_;

    /// Read boundary conditions for a zone.
    void read_bc(Int base, Int zone, Grid& grid);

    /// Read 1-to-1 connectivity for a zone (block interfaces / periodic BCs).
    void read_1to1(Int base, Int zone, Grid& grid);

    /// Read node coordinates for a zone.
    void read_coords(Int base, Int zone, Grid& grid);

    /// Determine which face a BC patch is on, from its point range.
    /// For a structured block, a BC on a face has one dimension collapsed.
    static int detect_face(const BCPatch& patch, Int ni, Int nj, Int nk);
};
