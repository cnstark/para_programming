#ifndef _MPI_ARG_HPP_
#define _MPI_ARG_HPP_
#include "mpi.h"



class MPIArg {
    public:
        int rank_size, rank;
        MPI_Comm comm;
        MPIArg(int *argc, char*** argv) {
            MPI_Init(argc, argv);
            MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            comm = MPI_COMM_WORLD;
        };
};

#endif
