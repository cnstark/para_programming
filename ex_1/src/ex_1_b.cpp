#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <ctime>
#include "mpi_arg.hpp"

#define DEBUG 0

#define VALUE(m, n, x, y) m[(x) * (n) + (y)]

#define FREE(x) if (x != nullptr) {\
        free(x); \
        x = nullptr; \
    }

#define TAG 10

// 根据i, j, n计算实际rank值
#define RANK(i, j, n) ((i) * (n) + (j))

// 输出矩阵
#define LOG_MATRIX(m, n) \
    std::cout << #m << ": " << std::endl; \
    for(int i = 0; i < n; i++) { \
        for(int j = 0; j < n; j++) \
            std::cout << std::fixed << std::setprecision(2) << VALUE(m, n, i, j) << "\t"; \
        std::cout << std::endl; \
    }

/**
 * 初始化矩阵为随机值
 * n为边长
 */
void init_matrix(double *matrix, int n) {
    for (int i = 0; i < n * n; i++) {
        matrix[i] = (rand() % 20) * 0.9;
    }
}

/**
 * 分割矩阵
 * 将origin_n * origin_n的原始矩阵origin分割成(origin_n / n_p) * (origin_n / n_p)的子矩阵
 * sub为(p_i, p_j)位置的子矩阵
 */
void split_matrix(double *sub, const double *origin, int origin_n, int p_i, int p_j, int n_p) {
    assert(origin_n % n_p == 0);
    int sub_n = origin_n / n_p;
    for (int i = p_i * sub_n, ii = 0; ii < sub_n; i++, ii++) {
        for (int j = p_j * sub_n, jj = 0; jj < sub_n; j++, jj++) {
            VALUE(sub, sub_n, ii, jj) = VALUE(origin, origin_n, i, j);
        }
    }
}

/**
 * 矩阵乘法
 */
void matrix_mul(double *c ,const double *a, const double *b, int n) {
    for (int i = 0; i < n * n; i++) {
        c[i] = 0;
    }
    for (int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
                VALUE(c, n, i, j) += VALUE(a, n, i, k) * VALUE(b, n, k, j);
            }
        }
    }
}

/**
 * 矩阵加法
 */
void matrix_add(double *c ,const double *a, const double *b, int n) {
    for (int i = 0; i < n * n; i++) {
        c[i] = a[i] + b[i];
    }
}

void mpi_fox(double *c, const double *a, const double *b, int sub_n, int p_i, int p_j, int n_p, MPI_Comm comm) {
    MPI_Status status;
    MPI_Comm split_world;
    // 将同行分为一个通信域
    MPI_Comm_split(comm, p_i, p_j, &split_world);

    // 申请临时变量空间
    auto *temp_a = (double *)malloc(sub_n * sub_n * sizeof(double));
    auto *temp_b = (double *)malloc(sub_n * sub_n * sizeof(double));
    auto *temp_c = (double *)malloc(sub_n * sub_n * sizeof(double));

    // 初始化结果c
    memset(c, 0, sub_n * sub_n * sizeof(double));

    // 拷贝b值temp_b
    memcpy(temp_b, b, sub_n * sub_n * sizeof(double));

    for (int k = 0; k < n_p; k++) {
        // 找到横向的根进程
        int i_root = (k + p_i) % n_p;
        // 根进程拷贝a值temp_a
        if (p_j == i_root) {
            memcpy(temp_a, a, sub_n * sub_n * sizeof(double));
        }
        // 根进程将temp_a广播值同行其他进程
        MPI_Bcast(temp_a, sub_n * sub_n, MPI_DOUBLE, i_root, split_world);

        // 计算矩阵乘法并与原来的c相加
        matrix_mul(temp_c, temp_a, temp_b, sub_n);
        matrix_add(c, c, temp_c, sub_n);

        // b块向上循环移动
        int sent_b_to = RANK(p_i == 0 ? (n_p - 1) : (p_i - 1), p_j, n_p);
        int recv_b_from = RANK(p_i == n_p - 1 ? 0 : (p_i + 1), p_j, n_p);
        MPI_Send(temp_b, sub_n * sub_n, MPI_DOUBLE, sent_b_to, TAG, comm);
        MPI_Recv(temp_b, sub_n * sub_n, MPI_DOUBLE, recv_b_from, TAG, comm, &status);
    }
    FREE(temp_a)
    FREE(temp_b)
    FREE(temp_c)
}

int main(int argc, char* argv[]) {
    const int n = 9; // 矩阵边长
    const int n_p = 3; // 一个维度上的处理器个数，即处理器个数开方

    int sub_n = n / n_p; // 每个子矩阵的边长

    MPIArg mpi_arg = MPIArg(&argc, &argv);

    if (mpi_arg.rank_size == 1) {
        auto *matrix_a = (double *)malloc(n * n * sizeof(double));
        auto *matrix_b = (double *)malloc(n * n * sizeof(double));
        auto *matrix_c = (double *)malloc(n * n * sizeof(double));

        // 初始化矩阵
        init_matrix(matrix_a, n);
        init_matrix(matrix_b, n);

        clock_t start = clock();
        // 矩阵乘
        matrix_mul(matrix_c, matrix_a, matrix_b, n);
        clock_t end = clock();
        std::cout << "sequential: " << std::endl;
        std::cout << "time: " << end - start << std::endl;
        #if DEBUG
        LOG_MATRIX(matrix_a, n)
        LOG_MATRIX(matrix_b, n)
        LOG_MATRIX(matrix_c, n)
        #endif
        FREE(matrix_a)
        FREE(matrix_b)
        FREE(matrix_c)
    } else {
        assert(n_p * n_p == mpi_arg.rank_size);

        // 申请空间
        auto *matrix_a = (double *)malloc(n * n * sizeof(double));
        auto *matrix_b = (double *)malloc(n * n * sizeof(double));
        auto *matrix_c = (double *)malloc(n * n * sizeof(double));

        auto *sub_matrix_a = (double *)malloc(sub_n * sub_n * sizeof(double));
        auto *sub_matrix_b = (double *)malloc(sub_n * sub_n * sizeof(double));
        auto *sub_matrix_c = (double *)malloc(sub_n * sub_n * sizeof(double));

        // 初始化矩阵
        init_matrix(matrix_a, n);
        init_matrix(matrix_b, n);

        // 分割矩阵
        split_matrix(sub_matrix_a, matrix_a, n, mpi_arg.rank / n_p, mpi_arg.rank % n_p, n_p);
        split_matrix(sub_matrix_b, matrix_b, n, mpi_arg.rank / n_p, mpi_arg.rank % n_p, n_p);

        clock_t start = clock();
        // fox算法
        mpi_fox(sub_matrix_c, sub_matrix_a, sub_matrix_b, sub_n, mpi_arg.rank / n_p, mpi_arg.rank % n_p, n_p, mpi_arg.comm);
        clock_t end = clock();
        if (mpi_arg.rank == 0) {
            std::cout << "parallel: " << std::endl;
            std::cout << "time: " << end - start << std::endl;
        }
        #if DEBUG
        // 0号进程输出原始矩阵及计算结果
        if (mpi_arg.rank == 0) {
            matrix_mul(matrix_c, matrix_a, matrix_b, n);
            LOG_MATRIX(matrix_a, n)
            LOG_MATRIX(matrix_b, n)
            LOG_MATRIX(matrix_c, n)
        }

        // 输出MPI计算结果
        // 为保证输出结果，从0号进程开始依次输出
        MPI_Status status;
        int msg = 0;
        if (mpi_arg.rank == 0) {
            std::cout << "Rank: " << mpi_arg.rank << std::endl;
            LOG_MATRIX(sub_matrix_c, sub_n)
            MPI_Send(&msg, 1, MPI_INT, 1, TAG, mpi_arg.comm);
        } else if (mpi_arg.rank == n_p * n_p - 1) {
            MPI_Recv(&msg, 1, MPI_INT, n_p * n_p - 2, TAG, mpi_arg.comm, &status);
            LOG_MATRIX(sub_matrix_c, sub_n)
        } else {
            MPI_Recv(&msg, 1, MPI_INT, mpi_arg.rank - 1, TAG, mpi_arg.comm, &status);
            std::cout << "Rank: " << mpi_arg.rank << std::endl;
            LOG_MATRIX(sub_matrix_c, sub_n)
            MPI_Send(&msg, 1, MPI_INT, mpi_arg.rank + 1, TAG, mpi_arg.comm);
        }
        #endif
        FREE(matrix_a)
        FREE(matrix_b)
        FREE(matrix_c)
        FREE(sub_matrix_a)
        FREE(sub_matrix_b)
        FREE(sub_matrix_c)
    }

    MPI_Finalize();
    return 0;
}
