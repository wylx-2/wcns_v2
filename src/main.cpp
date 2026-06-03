#include "wcns_v2/core/config.h"
#include "wcns_v2/field/field.h"
#include "wcns_v2/io/config_reader.h"
#include "wcns_v2/init/flow_initializer.h"
#include "wcns_v2/bc/bc_applier.h"
#include "wcns_v2/parallel/parallel_env.h"
#include "wcns_v2/parallel/parallel_manager.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <stdexcept>

// ============================================================================
// Helper: print config summary
// ============================================================================
static void print_config(const Config& cfg) {
    if (ParallelEnv::rank() != 0) return;
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "\n========== Configuration Summary ==========\n"
              << "  Physical:\n"
              << "    gamma   = " << cfg.gamma << "\n"
              << "    Prandtl = " << cfg.Prandtl << "\n"
              << "    Re      = " << cfg.Re << "\n"
              << "    Mach    = " << cfg.Mach << "\n"
              << "    AoA     = " << cfg.AoA << " deg\n"
              << "    beta    = " << cfg.beta << " deg\n"
              << "    eos_factor = " << cfg.eos_factor()
              << "  (p* = rho* * T* * eos_factor)\n"
              << "  Control:\n"
              << "    CFL          = " << cfg.cfl << "\n"
              << "    max_iter     = " << cfg.max_iter << "\n"
              << "    output_freq  = " << cfg.output_freq << "\n"
              << "    time_scheme  = " << cfg.time_scheme << "\n"
              << "    ghost layers = " << cfg.ng << "\n"
              << "  Initialization:\n"
              << "    init_type    = " << cfg.init_type << "\n"
              << "    wall_type    = " << cfg.wall_type << "\n"
              << "    body_force   = (" << cfg.body_force_x
              << ", " << cfg.body_force_y
              << ", " << cfg.body_force_z << ")\n"
              << "============================================\n\n";
}

// ============================================================================
// Helper: test Field conversion
// ============================================================================
static void test_field_conversion(const Config& cfg) {
    Int nci = 8, ncj = 6, nck = 4;

    Field f;
    f.allocate(nci, ncj, nck);

    Real rho_inf = 1.0;
    Real u_inf, v_inf, w_inf;
    cfg.free_stream_velocity(u_inf, v_inf, w_inf);
    Real p_inf = 1.0 / (cfg.gamma * cfg.Mach * cfg.Mach);

    for (Int k = 0; k < nck; ++k)
    for (Int j = 0; j < ncj; ++j)
    for (Int i = 0; i < nci; ++i) {
        f.prim.rho(i,j,k) = rho_inf;
        f.prim.u(i,j,k)   = u_inf;
        f.prim.v(i,j,k)   = v_inf;
        f.prim.w(i,j,k)   = w_inf;
        f.prim.p(i,j,k)   = p_inf;
    }

    f.prim_to_cons(cfg.gamma);
    f.prim.fill(0.0);
    f.cons_to_prim(cfg.gamma);

    Real max_err = 0.0;
    for (Int k = 0; k < nck; ++k)
    for (Int j = 0; j < ncj; ++j)
    for (Int i = 0; i < nci; ++i) {
        max_err = std::max(max_err, std::abs(f.prim.rho(i,j,k) - rho_inf));
        max_err = std::max(max_err, std::abs(f.prim.u(i,j,k)   - u_inf));
        max_err = std::max(max_err, std::abs(f.prim.v(i,j,k)   - v_inf));
        max_err = std::max(max_err, std::abs(f.prim.w(i,j,k)   - w_inf));
        max_err = std::max(max_err, std::abs(f.prim.p(i,j,k)   - p_inf));
    }

    if (ParallelEnv::rank() == 0) {
        std::cout << "  Field round-trip max error: " << max_err << "\n";
    }
}

// ============================================================================
// Helper: verify initialization on each local block
// ============================================================================
static void verify_initialization(const std::vector<LocalBlock>& blocks,
                                   const Config& cfg) {
    Real p_inf = 1.0 / (cfg.gamma * cfg.Mach * cfg.Mach);
    Real u_inf, v_inf, w_inf;
    cfg.free_stream_velocity(u_inf, v_inf, w_inf);

    for (const auto& lb : blocks) {
        Int ng = lb.grid.ng;
        Int i0 = ng, i1 = ng + lb.nci_core() - 1;
        Int j0 = ng, j1 = ng + lb.ncj_core() - 1;
        Int k0 = ng, k1 = ng + lb.nck_core() - 1;

        // --- Check interior values ---
        Real max_rho_err = 0, max_u_err = 0, max_p_err = 0;
        Real min_rho = 1e30, max_rho = -1e30;

        for (Int k = k0; k <= k1; ++k)
        for (Int j = j0; j <= j1; ++j)
        for (Int i = i0; i <= i1; ++i) {
            Real rho = lb.field.prim.rho(i,j,k);
            min_rho = std::min(min_rho, rho);
            max_rho = std::max(max_rho, rho);

            if (cfg.init_type == "uniform") {
                max_rho_err = std::max(max_rho_err, std::abs(rho - 1.0));
                max_u_err   = std::max(max_u_err,   std::abs(lb.field.prim.u(i,j,k) - u_inf));
                max_p_err   = std::max(max_p_err,   std::abs(lb.field.prim.p(i,j,k) - p_inf));
            }
        }

        std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                  << " interior: rho∈[" << min_rho << ", " << max_rho << "]";
        if (cfg.init_type == "uniform") {
            std::cout << " | max_err: rho=" << max_rho_err
                      << " u=" << max_u_err << " p=" << max_p_err;
        }
        std::cout << "\n";

        // --- Check ghost cells are non-zero (BC applied) ---
        // IMIN ghost: i=0..ng-1
        Real g_min = 1e30, g_max = -1e30;
        for (Int k = 0; k < lb.grid.nck; ++k)
        for (Int j = 0; j < lb.grid.ncj; ++j)
        for (Int d = 0; d < ng; ++d) {
            Real v = lb.field.prim.rho(d, j, k);
            g_min = std::min(g_min, v);
            g_max = std::max(g_max, v);
        }
        std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                  << " IMIN ghost rho∈[" << g_min << ", " << g_max << "]\n";

        // --- Corner check ---
        Real c_val = lb.field.prim.rho(0, 0, 0);
        std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                  << " corner (0,0,0) rho=" << c_val << "\n";
    }
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char* argv[]) {
    // ---- Initialize MPI ----
    ParallelEnv::init(argc, argv);

    if (argc < 2) {
        if (ParallelEnv::is_master()) {
            std::cerr << "Usage: " << argv[0] << " <grid.cgns> [config.ini]\n";
        }
        ParallelEnv::finalize();
        return 1;
    }

    try {
        // ---- Read configuration ----
        Config cfg;
        if (argc >= 3) {
            cfg = ConfigReader::read(argv[2]);
        } else {
            if (ParallelEnv::is_master()) {
                std::cout << "No config file provided — using default parameters.\n";
            }
            cfg.finalize();
        }
        print_config(cfg);

        // ---- Field conversion test (rank 0 only) ----
        if (ParallelEnv::is_master()) {
            std::cout << "--- Field Conversion Test ---\n";
        }
        test_field_conversion(cfg);

        // ---- Parallel initialization ----
        ParallelManager pm;
        std::vector<LocalBlock> local_blocks;
        pm.initialize(argv[1], cfg, local_blocks);

        // ---- Flow field initialization ----
        if (ParallelEnv::is_master()) {
            std::cout << "\n--- Flow Initialization: " << cfg.init_type << " ---\n";
        }
        for (auto& lb : local_blocks) {
            FlowInitializer::initialize(lb, cfg);
        }

        // ---- Boundary condition application ----
        if (ParallelEnv::is_master()) {
            std::cout << "\n--- Boundary Condition Application ---\n";
        }
        for (auto& lb : local_blocks) {
            BoundaryConditionApplier::apply_all(lb, cfg);
        }

        // ---- Halo exchange (parallel ghost cells) ----
        pm.exchange_all_halos(local_blocks);
        if (ParallelEnv::is_master()) {
            std::cout << "Halo exchange completed.\n";
        }

        // ---- Verify initialization ----
        verify_initialization(local_blocks, cfg);

        // ---- Verify metrics (GCL) ----
        if (ParallelEnv::is_master()) {
            std::cout << "\n--- Metric Verification (GCL) ---\n";
        }
        for (auto& lb : local_blocks) {
            Int i0 = lb.grid.ng;
            Int i1 = i0 + lb.nci_core() - 1;
            Int j0 = i0;
            Int j1 = j0 + lb.ncj_core() - 1;
            Int k0 = i0;
            Int k1 = k0 + lb.nck_core() - 1;

            bool fp[6];
            for (int f = 0; f < 6; ++f) fp[f] = lb.neighbors[f].is_periodic;

            MultiArray3D<Real> gx(lb.grid.nci, lb.grid.ncj, lb.grid.nck);
            MultiArray3D<Real> gy(lb.grid.nci, lb.grid.ncj, lb.grid.nck);
            MultiArray3D<Real> gz(lb.grid.nci, lb.grid.ncj, lb.grid.nck);

            Real dh = 1.0;

            InterpDiff::derivative(lb.grid.met_xi_x, gx, 0, dh, lb.grid.ng, fp);
            InterpDiff::derivative(lb.grid.met_eta_x, gy, 1, dh, lb.grid.ng, fp);
            InterpDiff::derivative(lb.grid.met_zeta_x, gz, 2, dh, lb.grid.ng, fp);

            Real gcl_x = 0;
            for (Int k = k0; k <= k1; ++k)
            for (Int j = j0; j <= j1; ++j)
            for (Int i = i0; i <= i1; ++i)
                gcl_x = std::max(gcl_x, std::abs(gx(i,j,k) + gy(i,j,k) + gz(i,j,k)));

            InterpDiff::derivative(lb.grid.met_xi_y, gx, 0, dh, lb.grid.ng, fp);
            InterpDiff::derivative(lb.grid.met_eta_y, gy, 1, dh, lb.grid.ng, fp);
            InterpDiff::derivative(lb.grid.met_zeta_y, gz, 2, dh, lb.grid.ng, fp);

            Real gcl_y = 0;
            for (Int k = k0; k <= k1; ++k)
            for (Int j = j0; j <= j1; ++j)
            for (Int i = i0; i <= i1; ++i)
                gcl_y = std::max(gcl_y, std::abs(gx(i,j,k) + gy(i,j,k) + gz(i,j,k)));

            InterpDiff::derivative(lb.grid.met_xi_z, gx, 0, dh, lb.grid.ng, fp);
            InterpDiff::derivative(lb.grid.met_eta_z, gy, 1, dh, lb.grid.ng, fp);
            InterpDiff::derivative(lb.grid.met_zeta_z, gz, 2, dh, lb.grid.ng, fp);

            Real gcl_z = 0;
            for (Int k = k0; k <= k1; ++k)
            for (Int j = j0; j <= j1; ++j)
            for (Int i = i0; i <= i1; ++i)
                gcl_z = std::max(gcl_z, std::abs(gx(i,j,k) + gy(i,j,k) + gz(i,j,k)));

            std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                      << " GCL: x=" << gcl_x << " y=" << gcl_y << " z=" << gcl_z << "\n";

            // Jacobian range
            Real jmin = 1e30, jmax = -1e30;
            for (Int k = k0; k <= k1; ++k)
            for (Int j = j0; j <= j1; ++j)
            for (Int i = i0; i <= i1; ++i) {
                jmin = std::min(jmin, lb.grid.jacobian(i,j,k));
                jmax = std::max(jmax, lb.grid.jacobian(i,j,k));
            }
            std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                      << " Jacobian range: [" << jmin << ", " << jmax << "]\n";
        }

        // ---- Global reduction test ----
        Real loc_max = local_blocks.empty() ? -1.0 : 1.0;
        Real glob_max = ParallelManager::global_max(loc_max);
        if (ParallelEnv::is_master()) {
            std::cout << "\nGlobal max test: " << glob_max << " (expected: 1.0)\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Error [Rank " << ParallelEnv::rank() << "]: " << e.what() << "\n";
        ParallelEnv::finalize();
        return 1;
    }

    ParallelEnv::finalize();
    return 0;
}
