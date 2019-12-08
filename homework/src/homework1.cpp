#include <iostream>
#include "mpi_arg.hpp"

#define COUNT 10
#define INVALID -999
#define END -9999

#define TAG 10

#define PRE_RANK (rank - 1)
#define NEXT_RANK (rank + 1) 

void start_rank(MPI_Comm comm) {
    MPI_Request request_s;
    MPI_Status status_s;
    int *value = (int *)malloc((COUNT + 1) * sizeof(int));
    for (int i = 0; i < COUNT; i++) {
        value[i] = i;
    }
    value[COUNT] = END;
    for (int i = 0; i < COUNT + 1; i++) {
        MPI_Isend(value + i, 1, MPI_INT, 1, TAG, comm, &request_s);
        MPI_Wait(&request_s, &status_s);
        std::cout << "Start Rank"
            << ", value: " << value[i] << std::endl;
    }
}

void end_rank(int rank, MPI_Comm comm) {
    MPI_Request request_r;
    MPI_Status status_r;
    int *recv = (int *)malloc(COUNT * sizeof(int));
    int *in_buf = (int *)malloc(sizeof(int));
    *in_buf = INVALID;
    int i = 0;
    do {
        if (*in_buf != INVALID) {
            std::cout << "End Rank"
                << ", value: " << std::dec << (*in_buf) << std::endl;
            recv[i++] = *in_buf;
        }
        MPI_Irecv(in_buf, 1, MPI_INT, PRE_RANK, TAG, comm, &request_r);
        MPI_Wait(&request_r, &status_r);
    } while (*in_buf != END);
    std::cout << "END: ";
    for (int j = 0; j < COUNT; j++) {
        std::cout << recv[j] << " ";
    }
    std::cout << std::endl;
}

void pipeline_rank(int rank, MPI_Comm comm) {
    MPI_Request request_r, request_s;
    MPI_Status status_r, status_s;

    int *Xbuf0, *Xbuf1, *Ybuf0, *Ybuf1;
    int *X, *Y, *Xin, *Yout;
    Xbuf0 = (int *)malloc(sizeof(int));
    *Xbuf0 = INVALID;
    Xbuf1 = (int *)malloc(sizeof(int));
    *Xbuf1 = INVALID;
    Ybuf0 = (int *)malloc(sizeof(int));
    *Ybuf0 = INVALID;
    Ybuf1 = (int *)malloc(sizeof(int));
    *Ybuf1 = INVALID;

    while(true) {
        if (X == Xbuf0) {
            X = Xbuf1;
            Y = Ybuf1;
            Xin = Xbuf0;
            Yout = Ybuf0;
        } else {
            X = Xbuf0;
            Y = Ybuf0;
            Xin = Xbuf1;
            Yout = Ybuf1;
        }

        if (*X != END) {
            MPI_Irecv(Xin, 1, MPI_INT, PRE_RANK, TAG, comm, &request_r);
        }

        if (*Yout != INVALID) {
            MPI_Isend(Yout, 1, MPI_INT, NEXT_RANK, TAG, comm, &request_s);
        }

        if (*X == END) {
            *Y = END;
            *Xin = END;
        } else if (*X != INVALID) {
            switch (rank) {
                case 1: // pipeline1: y = x * 2
                    *Y = (*X) * 2;
                    break;
                case 2: // pipeline2: y = x ^ 2
                    *Y = (*X) * (*X);
                    break;
                case 3: // pipeline3: y = x + 1
                    *Y = (*X) + 1;
                    break;
                default:
                    *Y = (*X);
                    break;
            }
            // std::cout " Rank: " << rank
            //     << ", value: " << std::dec << (*X) << std::endl;
        }

        if (*Yout != INVALID) {
            MPI_Wait(&request_s, &status_s);
        }
        if (*X != END) {
            MPI_Wait(&request_r, &status_r);
        }

        if (*Yout == END) {
            break;
        }
        
    }
    std::cout << "Rank: " << rank << " over" << std::endl;
}

void pipeline(int rank_size, int rank, MPI_Comm comm) {
    if (rank == 0) { // 头节点：产生数据流
        start_rank(comm);
    } else if (rank == rank_size - 1) { // 尾节点：接收数据流并输出
        end_rank(rank, comm);
    } else { // 中间节点：三级流水线
        pipeline_rank(rank, comm);
    }
}

int main(int argc, char* argv[]) {
    MPIArg *mpi_arg = new MPIArg(&argc, &argv);

    pipeline(mpi_arg -> rank_size, mpi_arg -> rank, mpi_arg -> comm);

    MPI_Finalize();
    return 0;
}
