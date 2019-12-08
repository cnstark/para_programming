#include <iostream>
#include "mpi_arg.hpp"

#define NODE_NUM 3
#define NODE_RANK_NUM 4

#define TAG 10

void mpi_bcast(int rank_size, int rank, MPI_Comm comm) {
    MPI_Status status_r;
    int *content = (int *)malloc(sizeof(int));

    if (rank == 0) {
        *content = 5555;
    }

    if (rank % NODE_RANK_NUM == 0) {
        // root节点广播
        if (rank == 0) {
            for (int i = 1; i < NODE_NUM; i++) {
                MPI_Send(content, 1, MPI_INT, i * NODE_RANK_NUM, TAG, comm);
            }
        } else {
            MPI_Recv(content, 1, MPI_INT, 0, TAG, comm, &status_r);
        }
    }

    if (rank % NODE_RANK_NUM == 0) {
        // send to sub
        for (int i = 1; i < NODE_RANK_NUM; i++) {
            MPI_Send(content, 1, MPI_INT, rank + i, TAG, comm);
        }
    } else {
        MPI_Recv(content, 1, MPI_INT, (rank / NODE_RANK_NUM) * NODE_RANK_NUM, TAG, comm, &status_r);
    }
    std::cout << "Rank: " << rank << ", content: " << *content << std::endl;
}

int main(int argc, char* argv[]) {
    MPIArg *mpi_arg = new MPIArg(&argc, &argv);

    mpi_bcast(mpi_arg -> rank_size, mpi_arg -> rank, mpi_arg -> comm);

    MPI_Finalize();
    return 0;
}
