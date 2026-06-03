#include "wcns_v2/core/config.h"
#include "wcns_v2/field/field.h"
#include "wcns_v2/io/config_reader.h"
#include "wcns_v2/init/flow_initializer.h"
#include "wcns_v2/bc/bc_applier.h"
#include "wcns_v2/parallel/parallel_env.h"
#include "wcns_v2/parallel/parallel_manager.h"
#include "wcns_v2/time/time_step.h"
#include "wcns_v2/scheme/wcns_interp.h"
#include "wcns_v2/scheme/riemann_solver.h"
#include "wcns_v2/scheme/inviscid_rhs.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
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
              << "    fixed_dt     = " << cfg.fixed_dt
              << (cfg.fixed_dt > 0 ? " (fixed)" : " (CFL-based)") << "\n"
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

        // ---- Time step test ----
        if (ParallelEnv::is_master()) {
            std::cout << "\n--- Time Step Test ---\n";
        }
        Real dt_local = TimeStep::compute(local_blocks, cfg);
        Real dt_global = ParallelManager::global_min(dt_local);
        if (ParallelEnv::is_master()) {
            std::cout << "  Local  Δt_min = " << dt_local << "\n";
            std::cout << "  Global Δt_min = " << dt_global << "\n";
            if (cfg.fixed_dt > 0) {
                std::cout << "  (using fixed_dt = " << cfg.fixed_dt << ")\n";
            }
        }

        // ---- WCNS interpolation test ----
        if (ParallelEnv::is_master()) {
            std::cout << "\n--- WCNS Interpolation Test ---\n";
        }
        {
            // Test each interpolation type
            const char* interp_types[] = {"weno_js", "mdcd_linear", "mdcd_hybrid"};
            // Use the default config for weno_js, and set MDCD params for the others
            Config cfg_weno = cfg;  // interp_type = "weno_js" (default in input)
            Config cfg_mdcd_lin = cfg;  cfg_mdcd_lin.interp_type = "mdcd_linear";
            Config cfg_mdcd_hyb = cfg;  cfg_mdcd_hyb.interp_type = "mdcd_hybrid";

            for (const char* type_name : interp_types) {
                Config* pcfg = &cfg_weno;
                if (std::string(type_name) == "mdcd_linear") pcfg = &cfg_mdcd_lin;
                if (std::string(type_name) == "mdcd_hybrid") pcfg = &cfg_mdcd_hyb;

                if (ParallelEnv::is_master()) {
                    std::cout << "\n  Testing interp_type = " << type_name
                              << " (diss=" << pcfg->mdcd_diss
                              << ", disp=" << pcfg->mdcd_disp
                              << ", sai_ref=" << pcfg->mdcd_sai_ref << ")\n";
                }

                auto interp = WcnsInterpBase::create(*pcfg);

                for (auto& lb : local_blocks) {
                    // Run interpolation for all 3 directions
                    interp->interp_xi(lb, *pcfg);

                    // Verify ξ-direction interpolation
                    Int ng = lb.grid.ng;
                    Int nci = lb.field.ni();
                    Int ncj = lb.field.nj();
                    Int nck = lb.field.nk();

                    // --- Boundary check: Q_L == Q_R at first 3 and last 3 faces ---
                    Real max_boundary_diff = 0.0;
                    for (Int k = 0; k < nck; ++k) {
                    for (Int j = 0; j < ncj; ++j) {
                        for (int h = 0; h <= 2; ++h) {
                            Real dl = std::abs(lb.field.ql_xi.rho(h,j,k) -
                                               lb.field.qr_xi.rho(h,j,k));
                            max_boundary_diff = std::max(max_boundary_diff, dl);
                        }
                        for (int h = nci; h >= nci - 2; --h) {
                            Real dl = std::abs(lb.field.ql_xi.rho(h,j,k) -
                                               lb.field.qr_xi.rho(h,j,k));
                            max_boundary_diff = std::max(max_boundary_diff, dl);
                        }
                    }}

                    std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                              << " [" << type_name << "] xi-boundary Q_L==Q_R max diff: "
                              << max_boundary_diff << "\n";

                    // --- Interior check: uniform flow → Q_L ≈ Q_R ≈ free-stream ---
                    Real max_interior_err = 0.0;
                    Int h0 = 4, h1 = nci - 4;
                    for (Int k = 4; k < nck - 4; ++k) {
                    for (Int j = 4; j < ncj - 4; ++j) {
                    for (Int h = h0; h <= h1; ++h) {
                        Real err_rho = std::abs(lb.field.ql_xi.rho(h,j,k) - 1.0);
                        max_interior_err = std::max(max_interior_err, err_rho);

                        Real diff = std::abs(lb.field.ql_xi.rho(h,j,k) -
                                             lb.field.qr_xi.rho(h,j,k));
                        max_interior_err = std::max(max_interior_err, diff);
                    }}}

                    std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                              << " [" << type_name << "] xi-interior max deviation: "
                              << max_interior_err << "\n";

                    // --- Spot-check a few interior values ---
                    Int ci = ng + 1, cj = ng + 1, ck = ng + 1;
                    Int face_i = ci + 1;
                    std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                              << " [" << type_name << "] Spot cell(" << ci << "," << cj
                              << "," << ck << "):"
                              << " cons.rho=" << lb.field.cons.rho(ci,cj,ck)
                              << " Q_L(" << face_i << ")=" << lb.field.ql_xi.rho(face_i,cj,ck)
                              << " Q_R(" << face_i << ")=" << lb.field.qr_xi.rho(face_i,cj,ck)
                              << "\n";

                    // Also interpolate eta and zeta directions (light test)
                    interp->interp_eta(lb, *pcfg);
                    interp->interp_zeta(lb, *pcfg);
                    std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                              << " [" << type_name << "] eta/zeta interpolation done."
                              << " Q_L_eta(" << ci << "," << cj+1 << "," << ck
                              << ").rho=" << lb.field.ql_eta.rho(ci,cj+1,ck)
                              << " Q_L_zeta(" << ci << "," << cj << "," << ck+1
                              << ").rho=" << lb.field.ql_zeta.rho(ci,cj,ck+1)
                              << "\n";
                }
            }
        }

        // ---- Riemann solver test ----
        if (ParallelEnv::is_master()) {
            std::cout << "\n--- Riemann Solver Test (Roe) ---\n";
        }
        {
            auto riemann = RiemannSolverBase::create(cfg);

            for (auto& lb : local_blocks) {
                Int ng = lb.grid.ng;
                Int nci = lb.field.ni();
                Int ncj = lb.field.nj();
                Int nck = lb.field.nk();

                // Re-run interpolation to ensure Q_L/Q_R are fresh
                auto interp = WcnsInterpBase::create(cfg);
                interp->interp_xi(lb, cfg);
                interp->interp_eta(lb, cfg);
                interp->interp_zeta(lb, cfg);

                // Solve Riemann problem in all 3 directions
                riemann->solve_xi(lb, cfg);
                riemann->solve_eta(lb, cfg);
                riemann->solve_zeta(lb, cfg);

                // ---- Verify: on uniform flow, dissipation must vanish ----
                Real u_inf, v_inf, w_inf;
                cfg.free_stream_velocity(u_inf, v_inf, w_inf);
                Real rho_inf = 1.0;
                Real p_inf = cfg.eos_factor();

                // ξ-direction verification at interior faces
                Int h0 = ng + 2, h1 = nci - ng - 2;
                Real max_flux_err = 0.0;

                for (Int k = 4; k < nck - 4; ++k) {
                for (Int j = 4; j < ncj - 4; ++j) {
                for (Int i = h0; i <= h1; ++i) {
                    Real Sx = lb.grid.face_xi_x(i,j,k);
                    Real Sy = lb.grid.face_xi_y(i,j,k);
                    Real Sz = lb.grid.face_xi_z(i,j,k);
                    Real U_inf = Sx*u_inf + Sy*v_inf + Sz*w_inf;

                    // Physical flux for uniform flow
                    Real F_phys[5];
                    F_phys[0] = rho_inf * U_inf;
                    F_phys[1] = rho_inf * u_inf * U_inf + p_inf * Sx;
                    F_phys[2] = rho_inf * v_inf * U_inf + p_inf * Sy;
                    F_phys[3] = rho_inf * w_inf * U_inf + p_inf * Sz;

                    Real rhoE_inf = p_inf / (cfg.gamma - 1.0)
                        + 0.5 * rho_inf * (u_inf*u_inf+v_inf*v_inf+w_inf*w_inf);
                    F_phys[4] = (rhoE_inf + p_inf) * U_inf;

                    Real F_roe[5] = {
                        lb.field.inv_xi.f1(i,j,k),
                        lb.field.inv_xi.f2(i,j,k),
                        lb.field.inv_xi.f3(i,j,k),
                        lb.field.inv_xi.f4(i,j,k),
                        lb.field.inv_xi.f5(i,j,k)
                    };

                    for (int c = 0; c < 5; ++c) {
                        Real err = std::abs(F_roe[c] - F_phys[c]);
                        Real ref = std::max(std::abs(F_phys[c]), 1.0e-12);
                        max_flux_err = std::max(max_flux_err, err / ref);
                    }
                }}}

                std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                          << " Roe xi-flux max relative error: " << max_flux_err
                          << " (expected: ~0 on uniform flow)\n";

                // Spot-check a few interior face values
                Int ci = ng + 2, cj = ng + 2, ck = ng + 2;
                Int fi = ci + 1;
                std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                          << " Spot xi-face(" << fi << "," << cj << "," << ck
                          << "): inv_xi.f1=" << lb.field.inv_xi.f1(fi,cj,ck)
                          << " f2=" << lb.field.inv_xi.f2(fi,cj,ck)
                          << " f5=" << lb.field.inv_xi.f5(fi,cj,ck)
                          << "\n";

                // η-direction spot check
                Int fj = cj + 1;
                std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                          << " Spot eta-face(" << ci << "," << fj << "," << ck
                          << "): inv_eta.f1=" << lb.field.inv_eta.f1(ci,fj,ck)
                          << " f2=" << lb.field.inv_eta.f2(ci,fj,ck)
                          << " f5=" << lb.field.inv_eta.f5(ci,fj,ck)
                          << "\n";

                // ζ-direction spot check
                Int fk = ck + 1;
                std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                          << " Spot zeta-face(" << ci << "," << cj << "," << fk
                          << "): inv_zeta.f1=" << lb.field.inv_zeta.f1(ci,cj,fk)
                          << " f2=" << lb.field.inv_zeta.f2(ci,cj,fk)
                          << " f5=" << lb.field.inv_zeta.f5(ci,cj,fk)
                          << "\n";
            }
        }

        // ---- Inviscid RHS test ----
        if (ParallelEnv::is_master()) {
            std::cout << "\n--- Inviscid RHS Test ---\n";
        }
        {
            auto interp = WcnsInterpBase::create(cfg);
            auto riemann = RiemannSolverBase::create(cfg);

            // Step 1: Interpolation + Riemann solver for all blocks
            for (auto& lb : local_blocks) {
                interp->interp_xi(lb, cfg);
                interp->interp_eta(lb, cfg);
                interp->interp_zeta(lb, cfg);

                riemann->solve_xi(lb, cfg);
                riemann->solve_eta(lb, cfg);
                riemann->solve_zeta(lb, cfg);
            }

            // Step 2: Exchange face fluxes at connectivity boundaries
            pm.exchange_flux_halos(local_blocks);

            // Step 3: Compute inviscid RHS for all blocks
            for (auto& lb : local_blocks) {
                Int nci = lb.field.ni();
                Int ncj = lb.field.nj();
                Int nck = lb.field.nk();

                // Compute inviscid RHS
                InviscidRHS::compute(lb);

                // ---- Verify: uniform flow → RHS ≡ 0 ----
                Int i0 = 3, i1 = nci - 4;
                Int j0 = 3, j1 = ncj - 4;
                Int k0 = 3, k1 = nck - 4;

                Real max_rhs = 0.0, max_rhs_far = 0.0;
                bool has_nan = false, has_inf = false;

                // far-field region: away from walls (j in middle third)
                Int jf0 = ncj / 3, jf1 = 2 * ncj / 3;

                for (Int k = k0; k <= k1; ++k) {
                for (Int j = j0; j <= j1; ++j) {
                for (Int i = i0; i <= i1; ++i) {
                    Real v = std::abs(lb.field.rhs.rho(i,j,k));
                    if (std::isnan(lb.field.rhs.rho(i,j,k))) has_nan = true;
                    if (std::isinf(lb.field.rhs.rho(i,j,k))) has_inf = true;
                    max_rhs = std::max(max_rhs, v);
                    if (j >= jf0 && j <= jf1) max_rhs_far = std::max(max_rhs_far, v);

                    v = std::abs(lb.field.rhs.rhou(i,j,k));
                    if (std::isnan(lb.field.rhs.rhou(i,j,k))) has_nan = true;
                    if (std::isinf(lb.field.rhs.rhou(i,j,k))) has_inf = true;
                    max_rhs = std::max(max_rhs, v);
                    if (j >= jf0 && j <= jf1) max_rhs_far = std::max(max_rhs_far, v);

                    v = std::abs(lb.field.rhs.rhov(i,j,k));
                    if (std::isnan(lb.field.rhs.rhov(i,j,k))) has_nan = true;
                    if (std::isinf(lb.field.rhs.rhov(i,j,k))) has_inf = true;
                    max_rhs = std::max(max_rhs, v);
                    if (j >= jf0 && j <= jf1) max_rhs_far = std::max(max_rhs_far, v);

                    v = std::abs(lb.field.rhs.rhow(i,j,k));
                    if (std::isnan(lb.field.rhs.rhow(i,j,k))) has_nan = true;
                    if (std::isinf(lb.field.rhs.rhow(i,j,k))) has_inf = true;
                    max_rhs = std::max(max_rhs, v);
                    if (j >= jf0 && j <= jf1) max_rhs_far = std::max(max_rhs_far, v);

                    v = std::abs(lb.field.rhs.rhoE(i,j,k));
                    if (std::isnan(lb.field.rhs.rhoE(i,j,k))) has_nan = true;
                    if (std::isinf(lb.field.rhs.rhoE(i,j,k))) has_inf = true;
                    max_rhs = std::max(max_rhs, v);
                    if (j >= jf0 && j <= jf1) max_rhs_far = std::max(max_rhs_far, v);

                }}}


                std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                          << " Inviscid RHS max |value|: "
                          << max_rhs << " (all interior)"
                          << " | far-field: " << max_rhs_far
                          << " | NaN=" << (has_nan ? "YES!" : "no")
                          << " Inf=" << (has_inf ? "YES!" : "no")
                          << "\n";

                // Spot-check interior cell values
                Int ci = 4, cj = 4, ck = 4;
                std::cout << "[Rank " << ParallelEnv::rank() << "] Block " << lb.block_id
                          << " Spot cell(" << ci << "," << cj << "," << ck
                          << "): rhs.rho=" << lb.field.rhs.rho(ci,cj,ck)
                          << " rhs.rhou=" << lb.field.rhs.rhou(ci,cj,ck)
                          << " rhs.rhoE=" << lb.field.rhs.rhoE(ci,cj,ck)
                          << "\n";
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error [Rank " << ParallelEnv::rank() << "]: " << e.what() << "\n";
        ParallelEnv::finalize();
        return 1;
    }

    ParallelEnv::finalize();
    return 0;
}
