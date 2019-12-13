#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <iomanip>
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

#define FILE_NAME "1.txt"
#define NODE_NUM 50

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

void mpi_dijkstra(MPIArg arg) {
    int n, s;
    double *matrix;
    double *part;
    if (arg.rank == 0) {
        s = NODE_NUM;
        FILE *fp = fopen(FILE_NAME, "r");
        if (fp == nullptr) {
            std::cout << "File name error!" << std::endl;
            return;
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

    int ep = (int)ceil(n / (double)arg.rank_size);
    int mynum = ep;
    if (arg.rank == arg.rank_size - 1) {
        mynum = n - ep * (arg.rank_size - 1);
    }
    // std::cout << "Rank: " << arg.rank << ", ep: " << ep << std::endl;
    // std::cout << "Rank: " << arg.rank << ", mynum: " << mynum << std::endl;
    // std::cout << "Rank: " << arg.rank << ", n: " << n << std::endl;
    // std::cout << "Rank: " << arg.rank << ", s: " << s << std::endl;
    part = (double *)malloc(sizeof(double) * ep * n);
    if (arg.rank == 0) {
        MPI_Scatter(matrix, ep * n, MPI_DOUBLE, part, ep * n, MPI_DOUBLE, 0, arg.comm);
        free(matrix);
    } else {
        MPI_Scatter(part, ep * n, MPI_DOUBLE, part, ep * n, MPI_DOUBLE, 0, arg.comm);
    }

    // for(int i = 0; i < mynum; i++) {
    //     std::cout << "Rank " << arg.rank << ":\t";
    //     for(int j = 0; j < n; j++) {
    //         std::cout << VALUE(part, n, i, j) << "\t";
    //     }
    //     std::cout << std::endl;
    // }

    MPI_Status status;
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


    for (int i = 0; i < n; i++) {
        int index = 0, index2;
        double num = 99999, num2;

        /**步骤(3.1)**/
        for (int j = 0; j < mynum; j++) {
            if (dist[j] < num && bdist[j] == false) {
                num = dist[j];
                index = ep * arg.rank + j;
            }
        }
        
        // std::cout << "Rank " << arg.rank << " debug" << std::endl;

        MPI_Barrier(arg.comm);
        // if (arg.rank == 1)
        // std::cout << "Rank: " << arg.rank << ", num: " << num << std::endl;

        /**步骤(3.2)**/
        int calnum = arg.rank_size;

        // while (calnum > 1) {
        //     if (calnum % 2 == 0) {
        //         calnum /= 2;
        //         if (arg.rank < calnum) {
        //             MPI_Recv(&index2, 1, MPI_INT, arg.rank + calnum, arg.rank, arg.comm, &status);
        //             MPI_Recv(&num2, 1, MPI_DOUBLE, arg.rank + calnum, arg.rank, arg.comm, &status);
        //             if(num2 < num) {
        //                 num = num2;
        //                 index = index2;
        //             }
        //         } else {
        //             MPI_Send(&index, 1, MPI_INT, arg.rank - calnum, arg.rank - calnum, arg.comm);
        //             MPI_Send(&num, 1, MPI_DOUBLE, arg.rank - calnum, arg.rank - calnum, arg.comm);
        //         }
        //     } else {
        //         if (arg.rank < calnum - 1) {
        //             MPI_Recv(&index2, 1, MPI_INT, arg.rank + calnum, arg.rank, arg.comm, &status);
        //             MPI_Recv(&num2, 1, MPI_DOUBLE, arg.rank + calnum, arg.rank, arg.comm, &status);
        //             if(num2 < num) {
        //                 num = num2;
        //                 index = index2;
        //             }
        //         } else if ()
        //     }
        // }



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


        // if (arg.rank == 1)
        // std::cout << "Rank: " << arg.rank << ", num: " << num << std::endl;

        /**步骤(3.3)**/
        MPI_Bcast(&index, 1, MPI_INT, 0, arg.comm);
        MPI_Bcast(&num, 1, MPI_DOUBLE, 0, arg.comm);
        // std::cout << "Rank: " << arg.rank << ", num: " << num << std::endl;

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
    for(int i = 0; i < mynum; i++) {
        std::cout << "node  " << arg.rank * ep + i << ":\t" << (int)dist[i] << std::endl;
    }

    // if (arg.rank == 0) {
    //     // MPI_Comm new_comm1;
    //     // MPI_Comm_split(arg.comm, 0, 0, &new_comm1);/
    //     for (int i = 1; i < arg.rank_size; i++) {
    //         MPI_Send(matrix + ep * (i - 1), ep * n, MPI_DOUBLE, i, i, arg.comm);
    //     }
    //     for (int i = 0; i < n; i++) {
    //         int index;
    //         double num;
    //         MPI_Bcast(&index, 1, MPI_INT, 1, arg.comm);
	// 	    MPI_Bcast(&num, 1, MPI_DOUBLE, 1, arg.comm);
    //         // std::cout << "Rank: " << arg.rank << ", index: " << index << std::endl;
    //         // std::cout << "Rank: " << arg.rank << ", num: " << num << std::endl;
    //         MPI_Barrier(arg.comm);
    //     }
    // } else {
    //     MPI_Status status;
    //     double *dist = (double *)malloc(sizeof(double) * ep);
	// 	bool *bdist = (bool *)malloc(sizeof(bool) * ep);
    //     matrix = (double *)malloc(sizeof(double) * ep * n);

    //     // MPI_Comm new_comm;
    //     // MPI_Comm_split(arg.comm, 1, arg.rank - 1, &new_comm);

    //     MPI_Recv(matrix, ep * n, MPI_DOUBLE, 0, arg.rank, arg.comm, &status);
        
    //     for(int i = 0; i < ep; i++) {
	// 		if(i + (arg.rank - 1) * ep == s) {
	// 			dist[i] = 0;
	// 			bdist[i] = true;
	// 		} else {
	// 			dist[i] = VALUE(matrix, n, i, s);
	// 			bdist[i] = false;
	// 		}
	// 	}

    //     // for(int i = 0; i < mynum; i++) {
	// 	// 	std::cout << "Rank " << arg.rank << ":\t";
	// 	// 	for(int j = 0; j < n; j++) {
    //     //         std::cout << VALUE(matrix, n, i, j) << "\t";
	// 	// 	}
	// 	// 	std::cout << std::endl;
	// 	// }
    //     for (int i = 0; i < n; i++) {
    //         int index = 0, index2;
    //         double num = 99999, num2;

    //         /**步骤(3.1)**/
    //         for (int j = 0; j < mynum; j++) {
    //             if (dist[j] < num && bdist[j] == false) {
    //                 num = dist[j];
    //                 index = ep * (arg.rank - 1) + j;
    //             }
    //         }
            
    //         // std::cout << "Rank " << arg.rank << " debug" << std::endl;

    //         // MPI_Barrier(new_comm);
    //         if (arg.rank == 1)
    //         std::cout << "Rank: " << arg.rank << ", num: " << num << std::endl;

    //         /**步骤(3.2)**/
    //         int calnum = arg.rank_size - 1;
    //         while (calnum > 1) {
    //             if (calnum % 2 == 0) {
    //                 calnum = calnum / 2;
    //                 if(arg.rank > calnum) {
    //                     MPI_Send(&index, 1, MPI_INT, arg.rank - calnum, arg.rank - calnum, arg.comm);
    //                     MPI_Send(&num, 1, MPI_DOUBLE, arg.rank - calnum, arg.rank - calnum, arg.comm);
    //                 } else {
    //                     MPI_Recv(&index2, 1, MPI_INT, arg.rank + calnum, arg.rank, arg.comm, &status);
    //                     MPI_Recv(&num2, 1, MPI_DOUBLE, arg.rank + calnum, arg.rank, arg.comm, &status);
    //                     if(num2 < num) {
    //                         num = num2;
    //                         index = index2;
    //                     }
    //                 }
    //             } else {
    //                 calnum = (calnum + 1) / 2;
    //                 if(arg.rank > calnum) {
    //                     MPI_Send(&index, 1, MPI_INT, arg.rank - calnum, arg.rank - calnum, arg.comm);
    //                     MPI_Send(&num, 1, MPI_DOUBLE, arg.rank - calnum, arg.rank - calnum, arg.comm);
    //                 } else if(arg.rank < calnum) {
    //                     MPI_Recv(&index2, 1, MPI_INT, arg.rank + calnum, arg.rank, arg.comm, &status);
    //                     MPI_Recv(&num2, 1, MPI_DOUBLE, arg.rank + calnum, arg.rank, arg.comm, &status);
    //                     if(num2 < num) {
    //                         num = num2;
    //                         index = index2;
    //                     }
    //                 }
    //             }
    //             // MPI_Barrier(new_comm);
    //         }
    //         // if (arg.rank == 1)
    //         // std::cout << "Rank: " << arg.rank << ", num: " << num << std::endl;

    //         /**步骤(3.3)**/
    //         MPI_Bcast(&index, 1, MPI_INT, 1, arg.comm);
	// 	    MPI_Bcast(&num, 1, MPI_DOUBLE, 1, arg.comm);

    //         /**步骤(3.4)**/
    //         for(int j = 0; j < mynum; j++) {
    //             if((bdist[j] == false) && (num + VALUE(matrix, n, j, index) < dist[j]))
    //                 dist[j] = num + VALUE(matrix, n, j, index);
    //         }

    //         /**步骤(3.5)**/
    //         if(arg.rank == index / ep + 1) {
    //             bdist[index % ep] = true;
    //         }

    //         MPI_Barrier(arg.comm);
    //     }

    //     for(int i = 0; i < mynum; i++) {
    //         std::cout << "node  " << (arg.rank - 1) * ep + i << ":\t" << (int)dist[i] << std::endl;
    //     }
    // }
}

int main(int argc, char *argv[]) {
    MPIArg mpi_arg = MPIArg(&argc, &argv);

    mpi_dijkstra(mpi_arg);

    MPI_Finalize();
    return 0;
}
