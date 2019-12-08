#include <iostream>
#include "mpi_arg.hpp"

#define TAG 10

void mpi_all2all(int rank_size, int rank, MPI_Comm comm) {
    MPI_Status status;

    int *send_buf = (int *)malloc(rank_size * sizeof(int));
    int *recv_buf = (int *)malloc(rank_size * sizeof(int));

    // init send_buf
    for (int i = 0; i < rank_size; i++) {
        send_buf[i] = 100 + rank * 10 + i;
    }

    // send
    for (int i = 0; i < rank_size; i++) {
        if (i == rank) {
            recv_buf[i] = send_buf[i];
        } else {
            MPI_Send(send_buf + i, 1, MPI_INT, i, TAG, comm);
        }
    }

    // log send_buf
    std::cout << "Rank: " << rank << ", send_buf: ";
    for (int j = 0; j < rank_size; j++) {
        std::cout << send_buf[j] << " ";
    }
    std::cout << std::endl;

    // recv
    for (int i = 0; i < rank_size; i++) {
        if (i != rank) {
            MPI_Recv(recv_buf + i, 1, MPI_INT, i, TAG, comm, &status);
        }
    }

    // log recv_buf
    std::cout << "Rank: " << rank << ", recv_buf: ";
    for (int j = 0; j < rank_size; j++) {
        std::cout << recv_buf[j] << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    MPIArg *mpi_arg = new MPIArg(&argc, &argv);

    mpi_all2all(mpi_arg -> rank_size, mpi_arg -> rank, mpi_arg -> comm);

    MPI_Finalize();
    return 0;
}

