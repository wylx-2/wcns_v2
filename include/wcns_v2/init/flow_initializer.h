#pragma once

#include "wcns_v2/core/config.h"
#include "wcns_v2/parallel/local_block.h"

/// @file flow_initializer.h
/// @brief Flow field initialization for different test cases.
///
/// Initialization functions fill only interior cells
/// (indices [ng, ng+nci_core-1] in each direction).
/// Ghost cells are filled later by BoundaryConditionApplier and HaloExchange.

class FlowInitializer {
public:
    /// Main entry: dispatch to the correct initializer based on cfg.init_type.
    /// Fills prim arrays and then converts to conservative variables.
    static void initialize(LocalBlock& lb, const Config& cfg);

    /// Uniform free-stream flow.
    static void init_uniform(LocalBlock& lb, const Config& cfg);

    /// Poiseuille channel flow: parabolic velocity profile with uniform density/pressure.
    /// The velocity profile is oriented along the direction of the body force.
    static void init_poiseuille(LocalBlock& lb, const Config& cfg);

    /// 2D Riemann problem: 4-quadrant discontinuous initial states.
    /// Each quadrant has uniform (rho, u, v, p). The configuration is selected
    /// by cfg.riemann_config (default 3).  Split at cfg.riemann_x/y_split.
    static void init_riemann_2d(LocalBlock& lb, const Config& cfg);

private:
    /// Get the interior cell loop range (inclusive).
    /// Interior cells are at [ng, ng + nci_core - 1] etc.
    static void interior_range(const LocalBlock& lb,
                               Int& i0, Int& i1,
                               Int& j0, Int& j1,
                               Int& k0, Int& k1);
};

#include "flow_initializer.hxx"
