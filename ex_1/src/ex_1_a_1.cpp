#include <iostream>
#include <cassert>
#include <ctime>
#include "mpi_arg.hpp"

#define DEBUG 0

#define FREE(x) if (x != nullptr) {\
        free(x); \
        x = nullptr; \
    }

/**
 * 获取步数
 */
int get_step(int n) {
    int step = 0;
    while (n != 1) {
        step++;
        n = n >> 1;
    }
    return step;
}

int mpi_sum(MPIArg arg, int value) {
    MPI_Status status;
    int step = get_step(arg.rank_size);
    int recv_value;
    for (int i = 0; i < step; i++) {
        if ((arg.rank % (1 << (i + 1))) < (1 << i)) {
            // 每组中的前半部分处理器
            MPI_Sendrecv(
                &value, 1, MPI_INT, arg.rank + (1 << i), arg.rank + (1 << i),
                &recv_value, 1, MPI_INT, arg.rank + (1 << i), arg.rank,
                arg.comm, &status);
        } else {
            // 每组中的后半部分处理器
            MPI_Sendrecv(
                &value, 1, MPI_INT, arg.rank - (1 << i), arg.rank - (1 << i),
                &recv_value, 1, MPI_INT, arg.rank - (1 << i), arg.rank,
                arg.comm, &status);
        }
        value += recv_value;
    }
    #if DEBUG
    std::cout << "Rank: " << arg.rank << ", sum: " << value << std::endl;
    #endif
    return value;
}

int main(int argc, char* argv[]) {
    MPIArg mpi_arg = MPIArg(&argc, &argv);

    // 确保rank_size为2的整数次幂
    assert((mpi_arg.rank_size & (mpi_arg.rank_size - 1)) == 0);

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

    // 蝶形求和
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
