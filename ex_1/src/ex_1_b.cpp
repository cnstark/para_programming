#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include "mpi_arg.hpp"

#define VALUE(m, n, x, y) m[(x) * (n) + (y)]

#define TAG 10

// 根据i, j, n计算实际rank值
#define RANK(i, j, n) ((i) * (n) + (j))

#define LOG_MATRIX(m, n) \
    std::cout << #m << ": " << std::endl; \
    for(int i = 0; i < n; i++) { \
        for(int j = 0; j < n; j++) \
            std::cout << std::fixed << std::setprecision(2) << VALUE(m, n, i, j) << "\t"; \
        std::cout << std::endl; \
    }

void init_matrix(double *matrix, int n) {
    for (int i = 0; i < n * n; i++) {
        matrix[i] = (rand() % 20) * 0.9;
    }
}

void split_matrix(double *sub, const double *origin, int origin_n, int p_i, int p_j, int n_p) {
    assert(origin_n % n_p == 0);
    int sub_n = origin_n / n_p;
    for (int i = p_i * sub_n, ii = 0; ii < sub_n; i++, ii++) {
        for (int j = p_j * sub_n, jj = 0; jj < sub_n; j++, jj++) {
            VALUE(sub, sub_n, ii, jj) = VALUE(origin, origin_n, i, j);
        }
    }
}

void matrix_mul(double *answer ,const double *a, const double *b, int n) {
    for (int i = 0; i < n * n; i++) {
        answer[i] = 0;
    }
    for (int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
                VALUE(answer, n, i, j) += VALUE(a, n, i, k) * VALUE(b, n, k, j);
            }
        }
    }
}

void matrix_add(double *answer ,const double *a, const double *b, int n) {
    for (int i = 0; i < n * n; i++) {
        answer[i] = a[i] + b[i];
    }
}

void mpi_fox(double *c, const double *a, const double *b, int sub_n, int p_i, int p_j, int n_p, MPI_Comm comm) {
    MPI_Status status;
    MPI_Comm split_world;
    MPI_Comm_split(comm, p_i, p_j, &split_world);

    auto *temp_a = (double *)malloc(sub_n * sub_n * sizeof(double));
    auto *temp_b = (double *)malloc(sub_n * sub_n * sizeof(double));
    auto *temp_c = (double *)malloc(sub_n * sub_n * sizeof(double));

    memset(c, 0, sub_n * sub_n * sizeof(double));

    memcpy(temp_b, b, sub_n * sub_n * sizeof(double));

    for (int k = 0; k < n_p; k++) {
        int i_root = (k + p_i) % n_p;
        if (p_j == i_root) {
            memcpy(temp_a, a, sub_n * sub_n * sizeof(double));
        }
        MPI_Bcast(temp_a, sub_n * sub_n, MPI_DOUBLE, i_root, split_world);

        matrix_mul(temp_c, temp_a, temp_b, sub_n);
        matrix_add(c, c, temp_c, sub_n);

        int sent_b_to = RANK(p_i == 0 ? (n_p - 1) : (p_i - 1), p_j, n_p);
        int recv_b_from = RANK(p_i == n_p - 1 ? 0 : (p_i + 1), p_j, n_p);
        MPI_Send(temp_b, sub_n * sub_n, MPI_DOUBLE, sent_b_to, TAG, comm);
        MPI_Recv(temp_b, sub_n * sub_n, MPI_DOUBLE, recv_b_from, TAG, comm, &status);
    }
}

int main(int argc, char* argv[]) {
    const int n = 9; // 矩阵边长
    const int n_p = 3; // 一个维度上的处理器个数，即处理器个数开方

    int sub_n = n / n_p; // 每个子矩阵的边长

    MPIArg mpi_arg = MPIArg(&argc, &argv);
    assert(n_p * n_p == mpi_arg.rank_size);

    // 申请空间
    auto *matrix_a = (double *)malloc(n * n * sizeof(double));
    auto *sub_matrix_a = (double *)malloc(sub_n * sub_n * sizeof(double));
    auto *matrix_b = (double *)malloc(n * n * sizeof(double));
    auto *sub_matrix_b = (double *)malloc(sub_n * sub_n * sizeof(double));
    auto *sub_matrix_c = (double *)malloc(sub_n * sub_n * sizeof(double));

    auto *matrix_c = (double *)malloc(n * n * sizeof(double));

    // 初始化矩阵
    init_matrix(matrix_a, n);
    init_matrix(matrix_b, n);

    // 分割矩阵
    split_matrix(sub_matrix_a, matrix_a, n, mpi_arg.rank / n_p, mpi_arg.rank % n_p, n_p);
    split_matrix(sub_matrix_b, matrix_b, n, mpi_arg.rank / n_p, mpi_arg.rank % n_p, n_p);

    // fox算法
    mpi_fox(sub_matrix_c, sub_matrix_a, sub_matrix_b, sub_n, mpi_arg.rank / n_p, mpi_arg.rank % n_p, n_p, mpi_arg.comm);

    // 输出原始矩阵及计算结果
    if (mpi_arg.rank == 0) {
        matrix_mul(matrix_c, matrix_a, matrix_b, n);
        LOG_MATRIX(matrix_a, n)
        LOG_MATRIX(matrix_b, n)
        LOG_MATRIX(matrix_c, n)
    }

    // 输出MPI计算结果
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

    MPI_Finalize();
    return 0;
}
