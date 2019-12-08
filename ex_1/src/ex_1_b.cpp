#include <iostream>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include "mpi_arg.hpp"

#define VALUE(m, n, x, y) m[(x) * (n) + (y)]

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
        matrix[i] = rand() % 20;
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

void mpi_fox(const double *a, const double *b, int sub_n, int p_i, int p_j, int n_p, MPI_Comm comm) {

}

int main(int argc, char* argv[]) {
    const int n = 9; // 矩阵边长
    const int n_p = 3; // 一个维度上的处理器个数，即处理器个数开方

    int sub_n = origin_n / n_p; // 每个子矩阵的边长

    MPIArg mpi_arg = MPIArg(&argc, &argv);
    assert(n_p * n_p == mpi_arg.rank_size);

    // 申请空间
    auto *matrix_a = (double *)malloc(n * n * sizeof(double));
    auto *sub_matrix_a = (double *)malloc(sub_n * sub_n * sizeof(double));
    auto *matrix_b = (double *)malloc(n * n * sizeof(double));
    auto *sub_matrix_b = (double *)malloc(sub_n * sub_n * sizeof(double));

    // 初始化矩阵
    init_matrix(matrix_a, n);
    init_matrix(matrix_b, n);

    // 分割矩阵
    split_matrix(sub_matrix_a, origin_matrix_a, n, mpi_arg.rank / n_p, mpi_arg.rank % n_p, n_p);
    split_matrix(sub_matrix_b, origin_matrix_b, n, mpi_arg.rank / n_p, mpi_arg.rank % n_p, n_p);

    // fox算法
    mpi_fox(sub_matrix_a, sub_matrix_b, sub_n, mpi_arg.rank / n_p, mpi_arg.rank % n_p, mpi_arg.comm);

    MPI_Finalize();
    return 0;
}
