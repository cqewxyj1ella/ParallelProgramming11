#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

const int n = 2048; // A, B, C are size nxn
// const int n = 16; // A, B, C are size nxn
int block_n, block_size;
int rank, nproc, sp;
int row_i, col_i;

void init_matrix(int** A, int** B);
void scatter_matrix(int** A, int** B, int* block_A, int* block_B);
void calculate(int* block_A, int* block_B, long int* block_C);

int main(int argc, char* argv[]) {
    double start_time, end_time;
    start_time = MPI_Wtime();

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // check if input is valid

    sp = (int)sqrt(nproc);
    if (sp * sp != nproc) {
        printf("the nproc must be quadratic number!\n");
        MPI_Finalize();
        exit(-1);
    }
    if ((n / sp) * sp != n) {
        printf("the 2048*4 must be divided by nproc!\n");
        MPI_Finalize();
        exit(-1);
    }

    int** A = (int**)malloc(sizeof(int*) * n);
    int** B = (int**)malloc(sizeof(int*) * n);
    long int** C = (long int**)malloc(sizeof(long int*) * n);
    for (int i = 0; i < n;i++) {
        A[i] = (int*)malloc(sizeof(int) * n);
        B[i] = (int*)malloc(sizeof(int) * n);
        C[i] = (long int*)malloc(sizeof(long int) * n);
    }
    block_n = n / sp; // side length of block
    block_size = block_n * block_n; // total elements in one block
    int* block_A = (int*)malloc(sizeof(int) * block_size);
    int* block_B = (int*)malloc(sizeof(int) * block_size);
    long int* block_C = (long int*)malloc(sizeof(long int) * block_size);
    memset(block_C, 0, sizeof(long int) * block_size);
    if (rank == 0) { // generate and scatter matrix to other ranks
        init_matrix(A, B);
        scatter_matrix(A, B, block_A, block_B);
        printf("%d %d\n", A[10][10], B[10][10]);
    }
    else { // recv matrix from rank 0
        MPI_Status status;
        MPI_Recv(block_A, block_size, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Recv(block_B, block_size, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
    }
    row_i = rank / sp;
    col_i = rank % sp;
    calculate(block_A, block_B, block_C);
    // now everyone has a block_C
    // rank 0 recv block_C from all other ranks
    // all other ranks send block_C to rank 0
    MPI_Status status;
    if (rank == 0) {
        int* recv_C = (int*)malloc(sizeof(int) * block_size);
        memcpy(recv_C, block_C, sizeof(int) * block_size);
        for (int i = 0; i < block_n; i++)
            for (int j = 0; j < block_n; j++)
                C[i][j] = recv_C[i * block_n + j];
        for (int i = 1; i < nproc; i++) {
            MPI_Recv(recv_C, block_size, MPI_INT, i, 5, MPI_COMM_WORLD, &status);
            // add recv_C to C
            int row_min = (i / sp) * block_n;
            int col_min = (i % sp) * block_n;
            for (int row = row_min; row < row_min + block_n; row++) {
                for (int col = col_min; col < col_min + block_n; col++) {
                    C[row][col] = recv_C[(row - row_min) * block_n + (col - col_min)];
                }
            }
        }
        // // print C
        // for (int i = 0; i < n; i++) {
        //     for (int j = 0; j < n;j++) {
        //         printf("%ld ", C[i][j]);
        //     }
        //     printf("\n");
        // }
    }
    else {
        MPI_Send(block_C, block_size, MPI_INT, 0, 5, MPI_COMM_WORLD);
    }
    // calculate total time
    end_time = MPI_Wtime();
    if (rank == 0) {
        printf("total time of parallel fox is: %f\n", end_time - start_time);
    }

    // let rank 0 calculate matrix multiply in serial
    if (rank == 0) {
        start_time = MPI_Wtime();
        int** serial_C = (int**)malloc(sizeof(int*) * n);
        for (int i = 0; i < n; i++) {
            serial_C[i] = (int*)malloc(sizeof(int) * n);
        }
        for (int i = 0; i < n;i++) {
            for (int j = 0; j < n;j++) {
                serial_C[i][j] = 0;
                for (int k = 0; k < n;k++) {
                    serial_C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        end_time = MPI_Wtime();
        printf("total time of serial is: %f\n", end_time - start_time);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(block_A);
    free(block_B);
    free(block_C);
    free(A);
    free(B);
    free(C);
    MPI_Finalize();
    return 0;
}

void init_matrix(int** A, int** B) {
    // A and B are matrix of size nxn
    srand((unsigned int)time(0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = rand() % 10;
            B[i][j] = rand() % 10;
        }
    }
}

void scatter_matrix(int** A, int** B, int* block_A, int* block_B) {
    // rank 0 will send block matrix in A, B to others
    // for index in original matrix from min to max
    int row_min, row_max, col_min, col_max;
    for (int i = 0; i < nproc; i++) {
        row_min = (i / sp) * block_n;
        row_max = row_min + block_n;
        col_min = (i % sp) * block_n;
        col_max = col_min + block_n;
        // copy items to sent to send_A and send_B
        int* send_A = (int*)malloc(sizeof(int) * block_size);
        int* send_B = (int*)malloc(sizeof(int) * block_size);
        int index = 0;
        for (int row = row_min; row < row_max; row++) {
            for (int col = col_min; col < col_max; col++) {
                send_A[index] = A[row][col];
                send_B[index] = B[row][col];
                index++;
            }
        }
        // send to rank i
        if (i == 0) {
            memcpy(block_A, send_A, sizeof(int) * block_size);
            memcpy(block_B, send_B, sizeof(int) * block_size);
        }
        else {
            MPI_Send(send_A, block_size, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(send_B, block_size, MPI_INT, i, 2, MPI_COMM_WORLD);
        }

    }
}

int get_proc_rank(int row, int col) {
    // get rank of proc in row and col
    return ((row + sp) % sp) * sp + ((col + sp) % sp);
}

void calculate(int* block_A, int* block_B, long int* block_C) {
    // use FOX algorithm to calculate matrix C
    // block_A, block_B, block_C are all 1D array of size = block_size
    MPI_Status status;
    int send_col_indx = row_i;
    int* recv_A = (int*)malloc(sizeof(int) * block_size);
    int* recv_B = (int*)malloc(sizeof(int) * block_size);
    for (int iter = 0; iter < sp; iter++) {
        if (col_i == send_col_indx) {
            // send block_A to same row procs
            int begin = get_proc_rank(row_i, 0);
            int end = get_proc_rank(row_i, sp - 1);
            for (int i = begin; i <= end; i++) {
                if (i != rank) {
                    MPI_Send(block_A, block_size, MPI_INT, i, 3, MPI_COMM_WORLD);
                }
            }
        }
        else {
            // recv block_A from same row procs
            MPI_Recv(recv_A, block_size, MPI_INT, get_proc_rank(row_i, send_col_indx), 3, MPI_COMM_WORLD, &status);
        }
        // update send_col_indx
        send_col_indx = (send_col_indx + 1) % sp;
        // calculate block_C
        for (int i = 0; i < block_n; i++) {
            for (int j = 0; j < block_n; j++) {
                int sum = 0;
                for (int k = 0; k < block_n; k++) {
                    sum += recv_A[i * block_n + k] * block_B[k * block_n + j];
                }
                block_C[i * block_n + j] += sum;
            }
        }
        // block_A stays the same, while block_B move to upper line
        MPI_Sendrecv(block_B, block_size, MPI_INT, get_proc_rank(row_i - 1, col_i), 4, \
            recv_B, block_size, MPI_INT, get_proc_rank(row_i + 1, col_i), 4, MPI_COMM_WORLD, &status);
        memcpy(block_B, recv_B, sizeof(int) * block_size);
    }
}