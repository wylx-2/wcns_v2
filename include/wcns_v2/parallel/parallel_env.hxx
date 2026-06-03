#pragma once

#include "parallel_env.h"
#include <iostream>

inline void ParallelEnv::init(int& argc, char**& argv) {
    if (init_) return;   // already initialized

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);

    init_ = true;

    if (rank_ == 0) {
        std::cout << "MPI initialized: " << size_ << " process(es)"
                  << " (thread support: " << provided << "/"
                  << MPI_THREAD_MULTIPLE << ")\n";
    }
}

inline void ParallelEnv::finalize() {
    if (!init_ || finalize_) return;
    finalize_ = true;
    MPI_Finalize();
}

inline int  ParallelEnv::rank()            { return rank_; }
inline int  ParallelEnv::size()            { return size_; }
inline bool ParallelEnv::is_master()       { return rank_ == 0; }
inline bool ParallelEnv::is_parallel()     { return size_ > 1; }

inline void ParallelEnv::barrier() {
    if (init_ && !finalize_) MPI_Barrier(MPI_COMM_WORLD);
}

inline MPI_Comm ParallelEnv::communicator() { return MPI_COMM_WORLD; }
