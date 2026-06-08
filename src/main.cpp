#include "wcns_v2/core/config.h"
#include "wcns_v2/core/residual.h"
#include "wcns_v2/core/history_monitor.h"
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
#include "wcns_v2/scheme/viscid_rhs.h"
#include "wcns_v2/scheme/body_force.h"
#include "wcns_v2/time/time_integrator.h"
#include "wcns_v2/time/lu_sgs.h"
#include "wcns_v2/io/solution_writer.h"
#include "wcns_v2/io/restart_writer.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>

// ============================================================================
// Helper: print config summary (rank 0 only)
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
              << "    max_time     = " << cfg.max_time
              << (cfg.max_time > 0 ? "" : " (disabled)") << "\n"
              << "    converge_tol = " << cfg.converge_tol << "\n"
              << "    output_freq  = " << cfg.output_freq << "\n"
              << "    restart_freq = " << cfg.restart_freq << "\n"
              << "    residual_freq= " << cfg.residual_freq << "\n"
              << "    time_scheme  = " << cfg.time_scheme << "\n"
              << "    ghost layers = " << cfg.ng << "\n"
              << "  Scheme:\n"
              << "    interp_type   = " << cfg.interp_type << "\n"
              << "    interp_vars   = " << cfg.interp_vars << "\n"
              << "    riemann_type  = " << cfg.riemann_type << "\n"
              << "    viscous_type  = " << cfg.viscous_type << "\n"
              << "    body_force    = (" << cfg.body_force_x
              << ", " << cfg.body_force_y
              << ", " << cfg.body_force_z
              << ") type=" << cfg.body_force_type << "\n"
              << "  Initialization:\n"
              << "    init_type    = " << cfg.init_type << "\n"
              << "    wall_type    = " << cfg.wall_type << "\n"
              << "============================================\n\n";
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
        // ================================================================
        // 1. Read configuration
        // ================================================================
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

        // ================================================================
        // 2. Parallel initialization (grid load, decomposition, metrics)
        // ================================================================
        ParallelManager pm;
        std::vector<LocalBlock> blocks;
        pm.initialize(argv[1], cfg, blocks);

        if (ParallelEnv::is_master()) {
            std::cout << "Grid loaded: " << pm.total_blocks()
                      << " total blocks on " << ParallelEnv::size()
                      << " rank(s)\n";
        }

        // ================================================================
        // 3. Flow field initialization
        // ================================================================
        if (ParallelEnv::is_master()) {
            std::cout << "\n--- Flow Initialization: " << cfg.init_type << " ---\n";
        }
        for (auto& lb : blocks) {
            FlowInitializer::initialize(lb, cfg);
        }

        // ================================================================
        // 4. Boundary conditions + halo exchange (primitive variables)
        // ================================================================
        // Correct ordering: face ghost → halo exchange → edge/corner ghost.
        // Edge/corner ghost fixup needs valid MPI-split ghost data, so the
        // halo exchange must happen BEFORE edge/corner (but AFTER face ghost).
        if (ParallelEnv::is_master()) {
            std::cout << "Applying boundary conditions...\n";
        }
        for (auto& lb : blocks) {
            BoundaryConditionApplier::apply_face_ghost(lb, cfg);
        }
        pm.exchange_all_halos(blocks);
        for (auto& lb : blocks) {
            BoundaryConditionApplier::apply_edge_ghost(lb);
            BoundaryConditionApplier::apply_corner_ghost(lb);
        }

        // ★ Critical: convert prim → cons so WCNS interpolation sees
        //   correct conservative ghost-cell values from the start.
        for (auto& lb : blocks) {
            lb.field.prim_to_cons(cfg.gamma);
        }

        if (ParallelEnv::is_master()) {
            std::cout << "Initialization complete.\n";
        }

        // ---- Write initial flow field (iter=0) ----
        {
            constexpr Int iter0 = 0;
            constexpr Real time0 = 0.0;
            if (ParallelEnv::is_master()) {
                std::cout << "Writing initial solution (iter=0)...\n";
            }
            SolutionWriter::write(blocks, cfg, iter0, time0);
        }

        // ================================================================
        // 5. Create interpolation and Riemann solver objects (reused)
        // ================================================================
        auto interp   = WcnsInterpBase::create(cfg);
        auto riemann  = RiemannSolverBase::create(cfg);
        bool has_viscous = (cfg.viscous_type != "none");

        int  n_stages = TimeIntegrator::n_stages(cfg);
        bool is_local = !ParallelEnv::is_parallel();

        // ================================================================
        // 6. Open residual log file (rank 0 only)
        // ================================================================
        std::string res_path = "residual.dat";
        std::string hist_path = "history.dat";
        if (!cfg.output_dir.empty()) {
            res_path  = cfg.output_dir + "/residual.dat";
            hist_path = cfg.output_dir + "/history.dat";
        }

        std::ofstream res_file;
        if (ParallelEnv::is_master()) {
            res_file.open(res_path);
            if (!res_file.is_open()) {
                throw std::runtime_error("Cannot open " + res_path + " for writing");
            }
            Residual::write_header(res_file);
        }

        // ---- History monitor: cross-section average velocity ----
        std::vector<Real> monitor_x_locations = {0.0, M_PI};
        std::ofstream hist_file;
        if (ParallelEnv::is_master()) {
            hist_file.open(hist_path);
            if (!hist_file.is_open()) {
                throw std::runtime_error("Cannot open " + hist_path + " for writing");
            }
            HistoryMonitor::write_header(hist_file, monitor_x_locations);
        }

        // ================================================================
        // 7. Main time loop
        // ================================================================
        Real time            = 0.0;
        Real last_output_time = 0.0;
        bool converged       = false;
        Int iter             = 0;

        // ---- Log initial history at iter=0 ----
        {
            auto section_avgs = HistoryMonitor::compute_averages(blocks, monitor_x_locations);
            if (ParallelEnv::is_master()) {
                HistoryMonitor::log(hist_file, 0, 0.0, 0.0, section_avgs);
                std::cout << "Initial U_avg(x=0) = " << section_avgs[0]
                          << ", U_avg(x=pi) = " << section_avgs[1] << "\n";
            }
        }

        if (ParallelEnv::is_master()) {
            std::cout << "\n========== Starting Time Loop ==========\n"
                      << "  Scheme: " << cfg.time_scheme
                      << " (" << n_stages << " stage(s))\n"
                      << "  Max iter: " << cfg.max_iter;
            if (cfg.max_time > 0.0) {
                std::cout << ", Max time: " << cfg.max_time;
            }
            std::cout << "\n=========================================\n\n";
        }

        for (iter = 1; iter <= cfg.max_iter; ++iter) {

            // ---- 7a. Compute time step ----
            Real dt = TimeStep::compute(blocks, cfg);
            dt = ParallelManager::global_min(dt);
            time += dt;

            // ---- 7b. Physical time termination ----
            if (cfg.max_time > 0.0 && time >= cfg.max_time) {
                if (ParallelEnv::is_master()) {
                    std::cout << "Reached max_time=" << cfg.max_time
                              << " at iter=" << iter << ", time=" << time << "\n";
                }
                break;
            }

            // ============================================================
            // RK stage loop
            // ============================================================
            for (int stage = 0; stage < n_stages; ++stage) {

                // ---- Save Q^(0) at the beginning of each time step ----
                if (stage == 0) {
                    for (auto& lb : blocks) {
                        // Use Field's built-in Q0 snapshot
                        Int nci = lb.field.ni();
                        Int ncj = lb.field.nj();
                        Int nck = lb.field.nk();
                        for (Int k = 0; k < nck; ++k)
                        for (Int j = 0; j < ncj; ++j)
                        for (Int i = 0; i < nci; ++i) {
                            lb.field.Q0.rho(i,j,k)  = lb.field.cons.rho(i,j,k);
                            lb.field.Q0.rhou(i,j,k) = lb.field.cons.rhou(i,j,k);
                            lb.field.Q0.rhov(i,j,k) = lb.field.cons.rhov(i,j,k);
                            lb.field.Q0.rhow(i,j,k) = lb.field.cons.rhow(i,j,k);
                            lb.field.Q0.rhoE(i,j,k) = lb.field.cons.rhoE(i,j,k);
                        }
                    }
                }

                // ========================================================
                // Inviscid RHS
                // ========================================================

                // WCNS interpolation: cell-center cons → face Q_L/Q_R
                for (auto& lb : blocks) {
                    interp->interp_xi(lb, cfg);
                    interp->interp_eta(lb, cfg);
                    interp->interp_zeta(lb, cfg);
                }

                // Riemann solver: Q_L/Q_R → inviscid face fluxes inv_xi/eta/zeta
                for (auto& lb : blocks) {
                    riemann->solve_xi(lb, cfg);
                    riemann->solve_eta(lb, cfg);
                    riemann->solve_zeta(lb, cfg);
                }

                // Exchange inviscid face-flux halos at connectivity boundaries
                // (also handles same-block periodic connections)
                pm.exchange_flux_halos(blocks);

                // 6th-order centered difference → inviscid RHS (writes to rhs)
                for (auto& lb : blocks) {
                    InviscidRHS::compute(lb);
                }

                // ========================================================
                // Viscous RHS (full pipeline: 5a → 5b → 5c → 5d → 5e)
                // ========================================================
                if (has_viscous) {

                    // 5a: Interpolate (u,v,w,T) from cell centers to faces
                    for (auto& lb : blocks) {
                        ViscidRHS::interp_to_faces(lb, cfg);
                    }

                    // 5b: Compute velocity/temperature gradients at cell centers
                    for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
                        ViscidRHS::compute_gradients(blocks[bi], blocks,
                                                      pm.flux_halo_ex(
                                                          static_cast<Int>(bi)),
                                                      cfg);
                    }

                    // 5b-exchange: gradient ghost cells
                    pm.exchange_gradient_halos(blocks);

                    // 5c: Compute cell-center Cartesian viscous flux vectors
                    for (auto& lb : blocks) {
                        ViscidRHS::compute_cell_viscous_flux(lb, cfg);
                    }

                    // 5c-exchange: viscous flux ghost cells
                    pm.exchange_viscous_flux_halos(blocks);

                    // 5d: Assemble viscous face fluxes (3 directions)
                    for (int dir = 0; dir < 3; ++dir) {
                        // Step 1: Interpolate Cartesian fluxes to faces (all blocks)
                        for (auto& lb : blocks) {
                            ViscidRHS::interp_cart_flux_to_faces(lb, dir, cfg);
                        }
                        // Step 2+3: Exchange face-interpolated values + assemble
                        for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
                            ViscidRHS::exchange_and_assemble_face_flux(
                                blocks[bi], dir, blocks,
                                pm.flux_halo_ex(static_cast<Int>(bi)), cfg);
                        }
                    }

                    // 5e: 6th-order diff → viscous RHS (accumulates into rhs)
                    for (auto& lb : blocks) {
                        ViscidRHS::compute_rhs(lb);
                    }
                }

                // ========================================================
                // Body force source term (accumulates into rhs)
                // ========================================================
                for (auto& lb : blocks) {
                    BodyForce::add_to_rhs(lb, cfg);
                }

                // ========================================================
                // Time integration — advance conservative variables
                // ========================================================
                if (cfg.time_scheme == "lu-sgs") {
                    for (auto& lb : blocks) {
                        LuSgs::advance(lb, cfg, dt);
                    }
                } else {
                    for (auto& lb : blocks) {
                        TimeIntegrator::advance_stage(lb, cfg, dt, stage,
                                                       lb.field.Q0);
                    }
                }

                // ========================================================
                // Post-stage: cons → prim, face-BC, halo, edge/corner-BC, prim → cons
                // ========================================================
                // Important: halo exchange must happen BETWEEN face-BC and
                // edge/corner-BC so that MPI-split ghost data is available
                // for the edge/corner ghost fixup.

                // Convert updated conservative variables to primitive
                for (auto& lb : blocks) {
                    lb.field.cons_to_prim(cfg.gamma);
                }

                // Apply face boundary conditions (fills prim ghost cells on BC faces)
                for (auto& lb : blocks) {
                    BoundaryConditionApplier::apply_face_ghost(lb, cfg);
                }

                // Halo exchange for primitive variables (fills MPI-split ghost cells)
                // Exchange primitive-variable halos (handles same-block periodic too)
                pm.exchange_all_halos(blocks);

                // Apply edge + corner ghost fixup (now has valid MPI ghost data)
                for (auto& lb : blocks) {
                    BoundaryConditionApplier::apply_edge_ghost(lb);
                    BoundaryConditionApplier::apply_corner_ghost(lb);
                }

                // ★ Critical: convert prim → cons so that conservative ghost
                //   cells are up-to-date for the next stage's WCNS interpolation.
                for (auto& lb : blocks) {
                    lb.field.prim_to_cons(cfg.gamma);
                }

                // ========================================================
                // Final stage: residual & flow monitor
                // ========================================================
                if (stage == n_stages - 1) {

                    if (iter % cfg.residual_freq == 0) {
                        auto res = Residual::compute(blocks);
                        auto mon = Residual::monitor(blocks, cfg);

                        // Divergence check
                        if (Residual::diverged(mon)) {
                            if (ParallelEnv::is_master()) {
                                std::cerr << "\nFATAL: Divergence detected at iter="
                                          << iter << ", time=" << time << "\n";
                                Residual::log(std::cerr, iter, dt, res, mon);
                            }
                            // Flush output before aborting
                            if (res_file.is_open()) {
                                Residual::log(res_file, iter, dt, res, mon);
                                res_file.close();
                            }
                            ParallelEnv::finalize();
                            return 1;
                        }

                        // Write to residual log
                        if (ParallelEnv::is_master()) {
                            Residual::log(res_file, iter, dt, res, mon);
                        }

                        // Cross-section average velocity monitoring
                        {
                            auto section_avgs = HistoryMonitor::compute_averages(
                                blocks, monitor_x_locations);
                            if (ParallelEnv::is_master()) {
                                HistoryMonitor::log(hist_file, iter, time, dt,
                                                     section_avgs);
                            }
                        }

                        // Convergence check
                        if (Residual::converged(res, cfg.converge_tol)) {
                            if (ParallelEnv::is_master()) {
                                std::cout << "\nConverged at iter=" << iter
                                          << ", time=" << time
                                          << " (res_rho=" << res.rho << ")\n";
                            }
                            converged = true;
                        }

                        // Progress to stdout
                        if (ParallelEnv::is_master() && iter % 10 == 0) {
                            std::cout << "  iter=" << std::setw(6) << iter
                                      << "  dt=" << std::scientific
                                      << std::setprecision(4) << dt
                                      << "  res_rho=" << res.rho
                                      << "  T_max=" << std::fixed
                                      << std::setprecision(4) << mon.T_max
                                      << "\n";
                        }
                    }
                }
            } // end RK stage loop

            if (converged) break;

            // ============================================================
            // 8. Solution output (step-based or time-based)
            // ============================================================
            bool step_output = (cfg.output_freq > 0 &&
                                iter % cfg.output_freq == 0);
            bool time_output = (cfg.output_time_interval > 0.0 &&
                                time - last_output_time >=
                                    cfg.output_time_interval - 1.0e-12);

            if (step_output || time_output) {
                if (ParallelEnv::is_master()) {
                    std::cout << "Writing solution at iter=" << iter
                              << ", time=" << time << "\n";
                }
                SolutionWriter::write(blocks, cfg, iter, time);
                if (time_output) last_output_time = time;
            }

            // ============================================================
            // 9. Restart output
            // ============================================================
            if (cfg.restart_freq > 0 && iter % cfg.restart_freq == 0) {
                if (ParallelEnv::is_master()) {
                    std::cout << "Writing restart at iter=" << iter
                              << ", time=" << time << "\n";
                }
                RestartWriter::write(blocks, cfg, iter, time);
            }

        } // end main time loop

        // ================================================================
        // 10. Final output
        // ================================================================
        if (ParallelEnv::is_master()) {
            std::cout << "\n========== Time Loop Finished ==========\n"
                      << "  Iterations: " << (iter > cfg.max_iter ? cfg.max_iter : iter) << "\n"
                      << "  Final time: " << time << "\n"
                      << "  Converged:  " << (converged ? "yes" : "no") << "\n"
                      << "=========================================\n";
        }

        // Final solution and restart write
        SolutionWriter::write(blocks, cfg, iter, time);
        RestartWriter::write(blocks, cfg, iter, time);

        // Close residual log and history log
        if (res_file.is_open()) {
            res_file.close();
        }
        if (hist_file.is_open()) {
            hist_file.close();
        }

    } catch (const std::exception& e) {
        std::cerr << "Error [Rank " << ParallelEnv::rank() << "]: "
                  << e.what() << "\n";
        ParallelEnv::finalize();
        return 1;
    }

    ParallelEnv::finalize();
    return 0;
}
