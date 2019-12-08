#include <iostream>
#include <cassert>
#include "mpi_arg.hpp"

#define TAG 10

#define N 34

int get_step(int n) {
    int step = 0;
    if ((n & (n - 1)) == 0) {
        while (n != 1) {
            step++;
            n = n >> 1;
        }
    } else {
        while (n != 0) {
            step++;
            n = n >> 1;
        }
    }
    return step;
}

void mpi_sum(MPIArg arg, int *arr, int n) {
    MPI_Status status;
    int rank = arg.rank;
    int rank_size = arg.rank_size;
    MPI_Comm comm = arg.comm;
    assert(rank_size == n);
    int step = get_step(n);
    int value = arr[rank];
    int recv_value;
    int i;
    for (i = 0; i < step; i++) {
        if (rank % (1 << (i + 1)) == 0) {
            int send_to = rank + (1 << i);
            if (send_to < rank_size) {
                MPI_Recv(&recv_value, 1, MPI_INT, send_to, TAG, comm, &status);
                value += recv_value;
            }
        } else {
            MPI_Send(&value, 1, MPI_INT, rank - (1 << i), TAG, comm);
            break;
        }
    }
    for (int j = step - 1; j >= 0; j--) {
        if (j < i) {
            int send_to = rank + (1 << j);
            if (send_to < rank_size) {
                MPI_Send(&value, 1, MPI_INT, send_to, TAG, comm);
            }
        } else if (j == i){
            MPI_Recv(&value, 1, MPI_INT, rank - (1 << j), TAG, comm, &status);
        }
    }
    std::cout << "Rank: " << rank << ", sum: " << value << std::endl;

}

int main(int argc, char* argv[]) {
    MPIArg mpi_arg = MPIArg(&argc, &argv);

    int *a = (int *)malloc(N * sizeof(int));
    for (int i = 0; i < N; i++) {
        a[i] = i;
    }
    if (mpi_arg.rank == 0) {
        std::cout << "arr: ";
        for (int i = 0; i < N; i++) {
        std::cout << a[i] << " ";
        }
        std::cout << std::endl;
    }
    mpi_sum(mpi_arg, a, N);

    MPI_Finalize();
    return 0;
}
