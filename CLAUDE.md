# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Run

```bash
# Configure
cmake -S . -B build

# Build
cmake --build build

# Run
./build/wcns_v2
```

- C++17, compiled with `-Wall -Wextra -Wpedantic`.
- `compile_commands.json` is generated at `build/compile_commands.json` for IDE language server support.

## Architecture

WCNS-based 3D compressible Navier-Stokes parallel solver on structured multi-block grids with MPI.

### Source layout

- `src/main.cpp` — entry point
- `src/core/` — solver driver, config, timer
- `src/grid/` — multi-block structured grid, metrics, I/O
- `src/field/` — conservative/primitive variables, thermodynamics
- `src/scheme/` — WCNS interpolation, compact differencing, Riemann solvers, viscous flux
- `src/bc/` — boundary conditions (farfield, wall, symmetry, inflow, outflow)
- `src/time/` — time integrators (RK3-TVD, RK4, LU-SGS)
- `src/parallel/` — MPI env, domain decomposition, halo exchange, parallel tridiagonal
- `src/io/` — config reader, solution output (Tecplot/VTK/Plot3D), restart
- `src/utils/` — type aliases, constants, math helpers, logger
- `include/` — public headers (for potential library use)
- `build/` — out-of-tree CMake build directory
