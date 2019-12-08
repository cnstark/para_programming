#include <iostream>
#include <cassert>
#include <random>
#include <ctime>
#include "mpi_arg.hpp"

#define TAG 10

#define N 50
#define P 6
#define Q 44

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

double mpi_sum(MPIArg arg, double value) {
    MPI_Status status;
    int rank = arg.rank;
    MPI_Comm comm = arg.comm;
    int step = get_step(P);
    double recv_value;
    int i;
    for (i = 0; i < step; i++) {
        if (rank % (1 << (i + 1)) == 0) {
            int send_to = rank + (1 << i);
            if (send_to < P) {
                MPI_Recv(&recv_value, 1, MPI_DOUBLE, send_to, TAG, comm, &status);
                value += recv_value;
            }
        } else {
            MPI_Send(&value, 1, MPI_DOUBLE, rank - (1 << i), TAG, comm);
            break;
        }
    }
    for (int j = step - 1; j >= 0; j--) {
        if (j < i) {
            int send_to = rank + (1 << j);
            if (send_to < P) {
                MPI_Send(&value, 1, MPI_DOUBLE, send_to, TAG, comm);
            }
        } else if (j == i){
            MPI_Recv(&value, 1, MPI_DOUBLE, rank - (1 << j), TAG, comm, &status);
        }
    }
    return value;
}

void param_service(MPIArg arg) {
    MPI_Status status;
    double value = 0, sum = 0;
    for (int i = arg.rank + P; i < arg.rank_size; i += P) {
        MPI_Recv(&value, 1, MPI_DOUBLE, i, TAG, arg.comm, &status);
        sum += value;
    }
    sum = mpi_sum(arg, sum);
    double ave = sum / Q;
    for (int i = arg.rank + P; i < arg.rank_size; i += P) {
        MPI_Send(&ave, 1, MPI_DOUBLE, i, TAG, arg.comm);
    }
}

void work(MPIArg arg) {
    std::default_random_engine e(arg.rank + time(NULL));
    double value = (e() % 100)* 1.0;
    // double value = arg.rank;
    double ave = 0;
    MPI_Status status;
    int rank = arg.rank;
    MPI_Comm comm = arg.comm;
    MPI_Send(&value, 1, MPI_DOUBLE, rank % P, TAG, comm);
    MPI_Recv(&ave, 1, MPI_DOUBLE, rank % P, TAG, comm, &status);
    std::cout << "Rank: " << rank << ", value: " << value << ", ave: " << ave << std::endl;
}


int main(int argc, char* argv[]) {
    MPIArg mpi_arg = MPIArg(&argc, &argv);
    assert(mpi_arg.rank_size == N);

    if (mpi_arg.rank < P) {
        param_service(mpi_arg);
    } else {
        work(mpi_arg);
    }

    MPI_Finalize();
    return 0;
}
