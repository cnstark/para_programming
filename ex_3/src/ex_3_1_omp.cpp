#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include "omp.h"
#include "mpi_arg.hpp"

#define VALUE(m, n, x, y) m[(x) * (n) + (y)]

#define EXIT(code, msg) \
    std::cout << msg << std::endl; \
    exit(code);

#define LOG_MATRIX(m, n) \
    std::cout << #m << ": " << std::endl; \
    for(int i = 0; i < n; i++) { \
        for(int j = 0; j < n; j++) \
            std::cout << std::fixed << std::setprecision(2) << VALUE(m, n, i, j) << "\t"; \
        std::cout << std::endl; \
    }

#define FREE(x) if (x != nullptr) {\
        free(x); \
        x = nullptr; \
    }

#define FILE_NAME "1.txt"
#define NODE_NUM 50

#define OMP_NUM 3

double get_next_num(FILE *fp) {
    char ch;

    while (!isdigit(ch)) {
        fread(&ch, sizeof(char), 1, fp);
    }
    if(isdigit(ch)) {
        char num_s[10];
        char *p = num_s;
        for (; isdigit(ch); p++) {
            *p = ch;
            fread(&ch, sizeof(char), 1, fp);
        }
        *p = 0;
        return atof(num_s);
    } else {
        return 99999;
        fread(&ch, sizeof(char), 1, fp);
    }
}

void mpi_omp_dijkstra(MPIArg arg) {
    int n, s;
    double *matrix;
    double *part;
    // 读取矩阵并广播n和s
    if (arg.rank == 0) {
        s = NODE_NUM;
        FILE *fp = fopen(FILE_NAME, "r");
        if (fp == nullptr) {
            EXIT(-1, "File name error!")
        }
        n = get_next_num(fp);
        matrix = (double *)malloc(sizeof(double) * n * n);
        
        for (int i = 0; i < n * n; i++) {
            matrix[i] = get_next_num(fp);
        }
        // LOG_MATRIX(matrix, n)

        fclose(fp);

        if(n < 0 || n > 10000 || n <= s) {
		    EXIT(-1, "n error!")
        }
        MPI_Bcast(&n, 1, MPI_INT, 0, arg.comm);
        MPI_Bcast(&s, 1, MPI_INT, 0, arg.comm);
    } else {
        MPI_Bcast(&n, 1, MPI_INT, 0, arg.comm);
        MPI_Bcast(&s, 1, MPI_INT, 0, arg.comm);
    }

    // 计算每个处理器负责的行数
    int ep = (int)ceil(n / (double)arg.rank_size);
    int mynum = ep;
    if (arg.rank == arg.rank_size - 1) {
        mynum = n - ep * (arg.rank_size - 1);
    }

    // 分配矩阵
    part = (double *)malloc(sizeof(double) * ep * n);
    if (arg.rank == 0) {
        MPI_Scatter(matrix, ep * n, MPI_DOUBLE, part, ep * n, MPI_DOUBLE, 0, arg.comm);
        FREE(matrix)
    } else {
        MPI_Scatter(part, ep * n, MPI_DOUBLE, part, ep * n, MPI_DOUBLE, 0, arg.comm);
    }

    MPI_Status status;

    // 初始化dist和bdist
    double *dist = (double *)malloc(sizeof(double) * ep);
    bool *bdist = (bool *)malloc(sizeof(bool) * ep);
    for(int i = 0; i < ep; i++) {
        if(i + arg.rank * ep == s) {
            dist[i] = 0;
            bdist[i] = true;
        } else {
            dist[i] = VALUE(part, n, i, s);
            bdist[i] = false;
        }
    }

    // 算法主循环
    for (int i = 0; i < n; i++) {
        int index = 0, index2;
        double num = 99999, num2;

        double omp_num[OMP_NUM];
        int omp_index[OMP_NUM];

        omp_set_num_threads(OMP_NUM);
        #pragma omp parallel for
        for (int a = 0; a < OMP_NUM; a++) {
            int k = omp_get_thread_num();
            omp_num[k] = 99999;
            omp_index[k] = 0;
            int et = (int)ceil(mynum / (double)OMP_NUM);
            int my_omp_num = et;
            if (k == OMP_NUM - 1) {
                my_omp_num = mynum - (OMP_NUM - 1) * my_omp_num;
            }
            for (int j = 0; j < my_omp_num; j++) {
                int p = j + k * et;
                if (dist[p] < omp_num[k] && bdist[p] == false) {
                    omp_num[k] = dist[p];
                    omp_index[k] = ep * arg.rank + p;
                }
            }
        }

        num =  omp_num[0];
        index = omp_index[0];
        for (int k = 1; k < OMP_NUM; k++) {
            if (omp_num[k] < num) {
                num =  omp_num[k];
                index = omp_index[k];
            }
        }


        /**步骤(3.1)**/
        // for (int j = 0; j < mynum; j++) {
        //     if (dist[j] < num && bdist[j] == false) {
        //         num = dist[j];
        //         index = ep * arg.rank + j;
        //     }
        // }

        MPI_Barrier(arg.comm);

        /**步骤(3.2)**/
        int calnum = arg.rank_size;
        while (calnum > 1) {
            if (calnum % 2 == 0) {
                if (arg.rank < calnum) {
                    calnum = calnum / 2;
                    if((arg.rank + 1) > calnum) {
                        if (arg.rank + 1 - calnum <= calnum) {
                            MPI_Send(&index, 1, MPI_INT, arg.rank - calnum, arg.rank - calnum, arg.comm);
                            MPI_Send(&num, 1, MPI_DOUBLE, arg.rank - calnum, arg.rank - calnum, arg.comm);
                        }
                    } else {
                        MPI_Recv(&index2, 1, MPI_INT, arg.rank + calnum, arg.rank, arg.comm, &status);
                        MPI_Recv(&num2, 1, MPI_DOUBLE, arg.rank + calnum, arg.rank, arg.comm, &status);
                        if(num2 < num) {
                            num = num2;
                            index = index2;
                        }
                    }
                } else {
                    calnum = calnum / 2;
                }
            } else {
                if (arg.rank < calnum) {
                    calnum = (calnum + 1) / 2;
                    if((arg.rank + 1) > calnum) {
                        MPI_Send(&index, 1, MPI_INT, arg.rank - calnum, arg.rank - calnum, arg.comm);
                        MPI_Send(&num, 1, MPI_DOUBLE, arg.rank - calnum, arg.rank - calnum, arg.comm);
                    } else if((arg.rank + 1) < calnum) {
                        MPI_Recv(&index2, 1, MPI_INT, arg.rank + calnum, arg.rank, arg.comm, &status);
                        MPI_Recv(&num2, 1, MPI_DOUBLE, arg.rank + calnum, arg.rank, arg.comm, &status);
                        if(num2 < num) {
                            num = num2;
                            index = index2;
                        }
                    }
                } else {
                    calnum = (calnum + 1) / 2;
                }
            }
            MPI_Barrier(arg.comm);
        }

        /**步骤(3.3)**/
        MPI_Bcast(&index, 1, MPI_INT, 0, arg.comm);
        MPI_Bcast(&num, 1, MPI_DOUBLE, 0, arg.comm);

        /**步骤(3.4)**/
        for(int j = 0; j < mynum; j++) {
            if((bdist[j] == false) && (num + VALUE(part, n, j, index) < dist[j]))
                dist[j] = num + VALUE(part, n, j, index);
        }

        /**步骤(3.5)**/
        if(arg.rank == index / ep) {
            bdist[index % ep] = true;
        }

        MPI_Barrier(arg.comm);
    }

    if (arg.rank == 0) {
        double *all_dist = (double *)malloc(sizeof(double) * (ep * arg.rank_size));
        MPI_Gather(dist, ep, MPI_DOUBLE, all_dist, ep, MPI_DOUBLE, 0, arg.comm);
        for(int i = 0; i < n; i++) {
            std::cout << "node  " << i << ":\t" << (int)all_dist[i] << std::endl;
        }
        FREE(all_dist)
    } else {
        MPI_Gather(dist, ep, MPI_DOUBLE, dist, ep, MPI_DOUBLE, 0, arg.comm);
    }

    FREE(bdist)
    FREE(dist)
    FREE(part)
}

int main(int argc, char *argv[]) {
    MPIArg mpi_arg = MPIArg(&argc, &argv);

    clock_t start = clock();
    mpi_omp_dijkstra(mpi_arg);
    clock_t end = clock();
    std::cout << "time: " << end - start << std::endl;

    MPI_Finalize();
    return 0;
}
