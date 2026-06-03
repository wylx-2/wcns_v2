#pragma once

#include <mpi.h>

/// @file parallel_env.h
/// @brief MPI environment wrapper for the WCNS parallel solver.
///
/// Static class that manages MPI initialization, finalization, and provides
/// rank/size/barrier queries.  All methods are safe to call from any context;
/// in single-process mode (no MPI or mpirun -np 1) the communication methods
/// become no-ops.

class ParallelEnv {
public:
    /// Initialize MPI.  Must be called exactly once before any other method.
    /// On subsequent calls, returns immediately (idempotent).
    static void init(int& argc, char**& argv);

    /// Finalize MPI.  Idempotent.
    static void finalize();

    /// Current process rank (0..size-1).
    static int rank();

    /// Total number of MPI processes.
    static int size();

    /// True if this is the master (rank 0) process.
    static bool is_master();

    /// True if more than one process is in use.
    static bool is_parallel();

    /// Barrier synchronization across all processes.
    static void barrier();

    /// The MPI communicator used by the solver.
    static MPI_Comm communicator();

private:
    inline static int  rank_    = 0;
    inline static int  size_    = 1;
    inline static bool init_    = false;
    inline static bool finalize_ = false;
};

#include "parallel_env.hxx"
