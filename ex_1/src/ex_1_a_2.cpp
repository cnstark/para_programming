#include <iostream>
#include <cassert>
#include <ctime>
#include "mpi_arg.hpp"

#define DEBUG 0

#define FREE(x) if (x != nullptr) {\
        free(x); \
        x = nullptr; \
    }

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

int mpi_sum(MPIArg arg, int value) {
    MPI_Status status;
    int step = get_step(arg.rank_size);
    int recv_value;
    int i;
    for (i = 0; i < step; i++) {
        if (arg.rank % (1 << (i + 1)) == 0) {
            int recv_from = arg.rank + (1 << i);
            if (recv_from < arg.rank_size) {
                MPI_Recv(&recv_value, 1, MPI_INT, recv_from, arg.rank, arg.comm, &status);
                value += recv_value;
            }
        } else {
            MPI_Send(&value, 1, MPI_INT, arg.rank - (1 << i), arg.rank - (1 << i), arg.comm);
            break;
        }
    }
    for (int j = step - 1; j >= 0; j--) {
        if (j < i) {
            int send_to = arg.rank + (1 << j);
            if (send_to < arg.rank_size) {
                MPI_Send(&value, 1, MPI_INT, send_to, send_to, arg.comm);
            }
        } else if (j == i){
            MPI_Recv(&value, 1, MPI_INT, arg.rank - (1 << j), arg.rank, arg.comm, &status);
        }
    }
    #if DEBUG
    std::cout << "Rank: " << arg.rank << ", sum: " << value << std::endl;
    #endif
    return value;
}

int main(int argc, char* argv[]) {
    MPIArg mpi_arg = MPIArg(&argc, &argv);

    // 数据生成（可替换为随机数）
    int *a = (int *)malloc(mpi_arg.rank_size * sizeof(int));
    for (int i = 0; i < mpi_arg.rank_size; i++) {
        a[i] = i;
    }

    #if DEBUG
    // 0号进程打印原始数组
    if (mpi_arg.rank == 0) {
        std::cout << "arr: ";
        for (int i = 0; i < mpi_arg.rank_size; i++) {
            std::cout << a[i] << " ";
        }
        std::cout << std::endl;
    }
    #endif

    // 二叉树求和
    clock_t start = clock();
    int sum = mpi_sum(mpi_arg, a[mpi_arg.rank]);
    clock_t end = clock();
    if (mpi_arg.rank == 0) {
        std::cout << "parallel: " << std::endl;
        std::cout << "sum: " << sum << std::endl;
        std::cout << "time: " << end - start << std::endl;
        int sum = 0;
        start = clock();
        for (int i = 0; i < mpi_arg.rank_size; i++) {
            sum += a[i];
        }
        end = clock();
        std::cout << "sequential: " << std::endl;
        std::cout << "sum: " << sum << std::endl;
        std::cout << "time: " << end - start << std::endl;
    }

    FREE(a)

    MPI_Finalize();
    return 0;
}
