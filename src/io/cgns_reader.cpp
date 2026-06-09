#include "io/cgns_reader.h"
#include <cgnslib.h>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

CGNSReader::CGNSReader() : file_id_(-1), is_open_(false), num_bases_(0) {}

CGNSReader::~CGNSReader() {
    if (is_open_) close();
}

void CGNSReader::open(const std::string& filename) {
    if (is_open_) close();

    int err = cg_open(filename.c_str(), CG_MODE_READ, &file_id_);
    if (err != CG_OK) {
        throw std::runtime_error("CGNS: cannot open file \"" + filename + "\": " +
                                 std::string(cg_get_error()));
    }
    is_open_ = true;

    // Count bases
    err = cg_nbases(file_id_, &num_bases_);
    if (err != CG_OK || num_bases_ <= 0) {
        close();
        throw std::runtime_error("CGNS: no bases found in file \"" + filename + "\"");
    }

    std::cout << "CGNS file \"" << filename << "\" opened: "
              << num_bases_ << " base(s)\n";
}

void CGNSReader::close() {
    if (is_open_) {
        cg_close(file_id_);
        file_id_ = -1;
        is_open_ = false;
        num_bases_ = 0;
    }
}

Int CGNSReader::num_bases() const { return num_bases_; }

Int CGNSReader::num_zones(Int base) const {
    int nz = 0;
    if (cg_nzones(file_id_, static_cast<int>(base), &nz) != CG_OK) return 0;
    return static_cast<Int>(nz);
}

void CGNSReader::read_zone(Int base, Int zone_idx, Grid& grid, Int ng) {
    int B = static_cast<int>(base);
    int Z = static_cast<int>(zone_idx);

    // --- Zone info ---
    char zone_name[128];
    cgsize_t sizes[9];  // sizes[0]=ni, sizes[1]=nj, sizes[2]=nk for structured
    if (cg_zone_read(file_id_, B, Z, zone_name, sizes) != CG_OK) {
        throw std::runtime_error("CGNS: failed to read zone " + std::to_string(zone_idx));
    }

    // For structured: sizes[0..2] are the vertex counts
    Int ni = static_cast<Int>(sizes[0]);
    Int nj = static_cast<Int>(sizes[1]);
    Int nk = static_cast<Int>(sizes[2]);

    if (ni < 2 || nj < 2 || nk < 2) {
        throw std::runtime_error("CGNS: zone has degenerate dimensions (need >= 2 in each dir)");
    }

    grid.name = std::string(zone_name);
    grid.allocate(ni, nj, nk, ng);

    std::cout << "  Zone \"" << zone_name << "\": " << ni << " x " << nj << " x " << nk
              << " vertices, " << (ni-1) << " x " << (nj-1) << " x " << (nk-1) << " cells\n";

    // --- Coordinates ---
    read_coords(B, Z, grid);

    // --- Cell centers ---
    grid.compute_cell_centers();
    grid.compute_cell_volumes();

    // --- Boundary conditions ---
    read_bc(B, Z, grid);

    // --- 1-to-1 connectivity (periodic / block interfaces) ---
    read_1to1(B, Z, grid);
}

void CGNSReader::read_coords(int B, int Z, Grid& grid) {
    // CGNS coordinate range: [1,ni], [1,nj], [1,nk]
    cgsize_t rmin[3] = {1, 1, 1};
    cgsize_t rmax[3] = {static_cast<cgsize_t>(grid.ni),
                        static_cast<cgsize_t>(grid.nj),
                        static_cast<cgsize_t>(grid.nk)};

    // Check if CoordinateX/Y/Z exist (some files use "CoordinateX" etc.)
    const char* coord_names[3] = {"CoordinateX", "CoordinateY", "CoordinateZ"};
    MultiArray3D<Real>* coord_arrays[3] = {&grid.node_x, &grid.node_y, &grid.node_z};

    for (int d = 0; d < 3; ++d) {
        // CGNS data type: try RealDouble first
        CGNS_ENUMT(DataType_t) dt = RealDouble;
        int err = cg_coord_read(file_id_, B, Z, coord_names[d], dt, rmin, rmax,
                                 coord_arrays[d]->data());
        if (err != CG_OK) {
            // Try RealSingle
            dt = RealSingle;
            err = cg_coord_read(file_id_, B, Z, coord_names[d], dt, rmin, rmax,
                                 coord_arrays[d]->data());
        }
        if (err != CG_OK) {
            throw std::runtime_error("CGNS: failed to read coordinate \"" +
                                     std::string(coord_names[d]) + "\": " +
                                     std::string(cg_get_error()));
        }
    }
}

void CGNSReader::read_bc(int B, int Z, Grid& grid) {
    int nbocos = 0;
    if (cg_nbocos(file_id_, B, Z, &nbocos) != CG_OK) {
        std::cerr << "  Warning: could not read BC count for zone\n";
        return;
    }

    for (int bc_idx = 1; bc_idx <= nbocos; ++bc_idx) {
        char boco_name[128];
        CGNS_ENUMT(BCType_t) boco_type;
        CGNS_ENUMT(PointSetType_t) ptset_type;
        cgsize_t npnts;
        int normal_index;
        cgsize_t normal_list_size;
        CGNS_ENUMT(DataType_t) normal_data_type;
        int ndataset;

        if (cg_boco_info(file_id_, B, Z, bc_idx, boco_name, &boco_type, &ptset_type,
                         &npnts, &normal_index, &normal_list_size,
                         &normal_data_type, &ndataset) != CG_OK) {
            std::cerr << "  Warning: could not read BC info for index " << bc_idx << "\n";
            continue;
        }

        // Read the point range (for structured: 2*index_dim values)
        // index_dim = 3 for 3D
        std::vector<cgsize_t> point_range(6);
        if (cg_boco_read(file_id_, B, Z, bc_idx, point_range.data(), nullptr) != CG_OK) {
            std::cerr << "  Warning: could not read BC point range for \"" << boco_name << "\"\n";
            continue;
        }

        // Read grid location (Vertex or CellCenter)
        CGNS_ENUMT(GridLocation_t) grid_loc;
        if (cg_boco_gridlocation_read(file_id_, B, Z, bc_idx, &grid_loc) != CG_OK) {
            grid_loc = GridLocationNull;  // default
        }

        // Build BCPatch
        BCPatch patch;
        patch.name = std::string(boco_name);
        patch.type = BoundaryCondition::from_cgns(static_cast<int>(boco_type));

        // CGNS point_range for structured grids:
        //   [imin, jmin, kmin, imax, jmax, kmax]  (1-based, inclusive)
        patch.imin = static_cast<Int>(point_range[0]);
        patch.jmin = static_cast<Int>(point_range[1]);
        patch.kmin = static_cast<Int>(point_range[2]);
        patch.imax = static_cast<Int>(point_range[3]);
        patch.jmax = static_cast<Int>(point_range[4]);
        patch.kmax = static_cast<Int>(point_range[5]);

        patch.face = detect_face(patch, grid.ni, grid.nj, grid.nk);

        grid.bc.add_patch(patch);

        std::cout << "    BC \"" << patch.name << "\": type=" << static_cast<int>(boco_type)
                  << " face=" << patch.face
                  << " range=[" << patch.imin << ":" << patch.imax
                  << ", " << patch.jmin << ":" << patch.jmax
                  << ", " << patch.kmin << ":" << patch.kmax << "]"
                  << " gridloc=" << static_cast<int>(grid_loc) << "\n";
    }

    std::cout << "  Read " << nbocos << " BC patches\n";
}

void CGNSReader::read_1to1(int B, int Z, Grid& grid) {
    int n1to1 = 0;
    if (cg_n1to1(file_id_, B, Z, &n1to1) != CG_OK) {
        // n1to1 may not exist in the file — that's fine
        return;
    }

    for (int ci = 1; ci <= n1to1; ++ci) {
        char connectname[128];
        char donorname[128];
        cgsize_t range[6];
        cgsize_t donor_range[6];
        int transform[3];

        if (cg_1to1_read(file_id_, B, Z, ci, connectname, donorname,
                         range, donor_range, transform) != CG_OK) {
            std::cerr << "  Warning: could not read 1-to-1 connection " << ci << "\n";
            continue;
        }

        Connectivity conn;
        conn.name       = std::string(connectname);
        conn.donor_name = std::string(donorname);

        // CGNS range: [imin, jmin, kmin, imax, jmax, kmax] (1-based, inclusive)
        conn.imin = static_cast<Int>(range[0]);
        conn.jmin = static_cast<Int>(range[1]);
        conn.kmin = static_cast<Int>(range[2]);
        conn.imax = static_cast<Int>(range[3]);
        conn.jmax = static_cast<Int>(range[4]);
        conn.kmax = static_cast<Int>(range[5]);

        conn.donor_imin = static_cast<Int>(donor_range[0]);
        conn.donor_jmin = static_cast<Int>(donor_range[1]);
        conn.donor_kmin = static_cast<Int>(donor_range[2]);
        conn.donor_imax = static_cast<Int>(donor_range[3]);
        conn.donor_jmax = static_cast<Int>(donor_range[4]);
        conn.donor_kmax = static_cast<Int>(donor_range[5]);

        conn.transform[0] = transform[0];
        conn.transform[1] = transform[1];
        conn.transform[2] = transform[2];

        // Try to read periodic properties
        float rot_ctr[3] = {0, 0, 0};
        float rot_ang[3] = {0, 0, 0};
        float trans[3]   = {0, 0, 0};
        int perr = cg_1to1_periodic_read(file_id_, B, Z, ci,
                                         rot_ctr, rot_ang, trans);
        conn.is_periodic = (perr == CG_OK);
        if (conn.is_periodic) {
            conn.translation[0]    = static_cast<Real>(trans[0]);
            conn.translation[1]    = static_cast<Real>(trans[1]);
            conn.translation[2]    = static_cast<Real>(trans[2]);
            conn.rotation_center[0] = static_cast<Real>(rot_ctr[0]);
            conn.rotation_center[1] = static_cast<Real>(rot_ctr[1]);
            conn.rotation_center[2] = static_cast<Real>(rot_ctr[2]);
            conn.rotation_angle[0]  = static_cast<Real>(rot_ang[0]);
            conn.rotation_angle[1]  = static_cast<Real>(rot_ang[1]);
            conn.rotation_angle[2]  = static_cast<Real>(rot_ang[2]);
        } else {
            conn.translation[0] = conn.translation[1] = conn.translation[2] = 0;
            conn.rotation_center[0] = conn.rotation_center[1] = conn.rotation_center[2] = 0;
            conn.rotation_angle[0]  = conn.rotation_angle[1]  = conn.rotation_angle[2]  = 0;
        }

        grid.connections.add(conn);

        std::cout << "  1-to-1[" << ci << "]: \"" << conn.name
                  << "\" donor=\"" << conn.donor_name << "\""
                  << " range=[" << conn.imin << ":" << conn.imax
                  << "," << conn.jmin << ":" << conn.jmax
                  << "," << conn.kmin << ":" << conn.kmax << "]"
                  << " donor_range=[" << conn.donor_imin << ":" << conn.donor_imax
                  << "," << conn.donor_jmin << ":" << conn.donor_jmax
                  << "," << conn.donor_kmin << ":" << conn.donor_kmax << "]"
                  << " transform=[" << conn.transform[0]
                  << "," << conn.transform[1]
                  << "," << conn.transform[2] << "]";
        if (conn.is_periodic) {
            std::cout << " periodic=["
                      << conn.translation[0] << ","
                      << conn.translation[1] << ","
                      << conn.translation[2] << "]";
        }
        std::cout << "\n";
    }

    std::cout << "  Read " << n1to1 << " 1-to-1 connection(s)\n";
}

int CGNSReader::detect_face(const BCPatch& patch, Int ni, Int nj, Int nk) {
    // CGNS point_range: [imin, jmin, kmin, imax, jmax, kmax]
    // A face has one dimension collapsed.
    if (patch.imin == patch.imax) {
        return (patch.imin == 1) ? 0 : 1;  // IMIN or IMAX
    }
    if (patch.jmin == patch.jmax) {
        return (patch.jmin == 1) ? 2 : 3;  // JMIN or JMAX
    }
    if (patch.kmin == patch.kmax) {
        return (patch.kmin == 1) ? 4 : 5;  // KMIN or KMAX
    }
    (void)ni; (void)nj; (void)nk;
    return -1;
}

void CGNSReader::read_all(std::vector<Grid>& grids, Int ng) {
    grids.clear();
    for (Int z = 1; z <= num_zones(1); ++z) {
        Grid g;
        read_zone(1, z, g, ng);
        grids.push_back(std::move(g));
    }
}
