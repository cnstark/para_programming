#include <iostream>
#include <iomanip>
#include "mpi_arg.hpp"

#define M 4

#define VALUE(m, x, y) m[(x) * M + (y)]
#define a(x, y) VALUE(a, x, y)
#define A(x, y) VALUE(A, x, y)
#define L(x, y) VALUE(L, x, y)
#define U(x, y) VALUE(U, x, y)
#define recv(x, y) VALUE(recv, x, y)

#define LOG_MATRIX(matrix, x_size, y_size) \
    for(int i = 0; i < x_size; i++) { \
        for(int j = 0; j < y_size; j++) \
            std::cout << std::fixed << std::setprecision(4) << matrix(i, j) << " "; \
        std::cout << std::endl; \
    }


void mpi_lu_decomposion(float *A, int rank_size, int rank, MPI_Comm comm) {
    MPI_Status status;
    // 各进程处理的行数
    int m = M / rank_size;
    if (M % rank_size != 0) m++;

    float *a = (float *)malloc(m * M * sizeof(float));
    if (rank == 0) {
        for(int i = 0; i < m; i++) {
            for(int j = 0; j < M; j++) {
                a(i, j) = A(i, j);
            }
        }
        for(int i = m; i < M; i++) {
            MPI_Send(&A(i, 0), M, MPI_FLOAT, i / m, i % m, comm);
        }
    } else if (rank * m < M) {
        for(int i = 0; i < m; i++) {
	        if(rank * m + i >= M) break;
            MPI_Recv(&a(i, 0), M, MPI_FLOAT, 0, i, comm, &status);
    	}
    }

    float *f = (float *)malloc(M * sizeof(float));
    for(int i = 0; i < rank_size; i++) {
	    if(i * m >= M) break;
        for(int j = 0; j < m; j++) {
            int v;
            if(i * m + j >= M) break;
            // 行变换
            if (rank == i) {
                v = i * m + j;
	            for (int k = v; k < M; k++) {
                    f[k] = a(j, k);
                }
                MPI_Bcast(f, M, MPI_FLOAT, rank, comm);
            } else {
                v = i * m + j;
                MPI_Bcast(f, M, MPI_FLOAT, i, comm);
            }

            if (rank == i) {
                for(int k = j + 1; k < m; k++) {
                    a(k, v) = a(k, v) / f[v];
                    for(int w = v+1; w < M; w++)
                        a(k, w) = a(k, w) - f[w] * a(k, v);
                }
            }
        
            if(rank > i) {
                for(int k = 0; k < m; k++) {
                    a(k, v) = a(k, v) / f[v];
                    for(int w = v + 1; w < M; w++)
                        a(k, w) = a(k, w) - f[w] * a(k,v);
                }
            }
        }
    }

    if (rank != 0 && rank * m < M) {
        for(int i = 0; i < m; i++) {
	        if(rank * m + i >= M) break;
            MPI_Send(&a(i, 0), M, MPI_FLOAT, 0, m * rank + i, comm);
        }
    }
    
    if (rank == 0) {
        float *recv = (float *)malloc(M * M * sizeof(float));
        for(int i = 0; i < m; i++) {
            for(int j = 0; j < M; j++) {
                recv(i, j) = a(i, j);
            }
        }

        for(int i = 1; i < rank_size; i++) {
            for(int j = 0; j < m; j++) {
	            if(i * m + j >= M) break;
                MPI_Recv(&a(j, 0), M, MPI_FLOAT, i, m * i + j, comm, &status);
                for(int k = 0; k < M; k++) {
                    recv(i * m + j, k) = a(j, k);
                }
            }
        }
        
        // 初始化U
        float *L = (float *)malloc(M * M * sizeof(float));
        float *U = (float *)malloc(M * M * sizeof(float));
        for(int i = 0; i < M; i++) {
            for(int j = 0; j < M; j++) {
                U(i, j) = 0.0;
            }
        }
        // 初始化L
        for(int i = 0; i < M; i++) {
            for(int j = 0; j < M; j++) {
                if (i == j)
                    L(i, j) = 1.0;
                else
                    L(i, j) = 0.0;
            }
        }
        // 为L, U赋值
        for(int i = 0; i < M; i++) {
            for(int j = 0; j < M; j++) {
                if (i > j)
                    L(i, j) = recv(i, j);
                else
                    U(i, j) = recv(i, j);
            }
        }
        
        // 输出
        std::cout << "A:" << std::endl;
        LOG_MATRIX(A, M, M)
        std::cout << "L:" << std::endl;
        LOG_MATRIX(L, M, M)
        std::cout << "U:" << std::endl;
        LOG_MATRIX(U, M, M)
    }
}

int main(int argc, char* argv[]) {
    MPIArg *mpi_arg = new MPIArg(&argc, &argv);

    float A[M * M] = {4, 2, 1, 5,
                      8, 7, 2, 10,
                      4, 8, 3, 6,
                      6, 8, 4, 9};

    mpi_lu_decomposion(A, mpi_arg -> rank_size, mpi_arg -> rank, mpi_arg -> comm);

    MPI_Finalize();
    return 0;
}