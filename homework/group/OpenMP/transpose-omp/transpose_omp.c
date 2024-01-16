#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"
#include "math.h"
#include <omp.h>

#define E 0.0001
#define a(x,y) a[x*m+y]
#define b(x,y) b[x*m+y]
#define A(x,y) A[x*size+y]
#define B(x,y) B[x*size+y]
#define intsize sizeof(int)
#define floatsize sizeof(float)

#define charsize sizeof(char)

int size, N;                                       /* size:矩阵的阶数; N:矩阵的阶数 */
int m;                                            /* 每个子矩阵的大小 */
int t;                                            /* 子矩阵的分块数 */
float* A, * B;                                     /* A:原始矩阵; B:转置后的矩阵 */
double starttime;                                 /* 程序开始时间 */
double time1;                                      /* 读取文件的时间 */
double time2;                                     /* 并行计算的时间 */
int my_rank;                                      /* 当前进程的排名 */
int p;                                           /* 进程数 */
MPI_Status status;                                /* MPI状态 */
FILE* fdA;                                       /* 文件指针 */

/* 释放环境资源,释放动态分配的内存 */
void Environment_Finalize(float* a, float* b)
{
    free(a);
    free(b);
}

int main(int argc, char** argv)
{
    int i, j, k, my_rank, group_size;
    float* a, * b;
    int u, v;
    float temp;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &group_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    p = group_size;

    /* 读取输入矩阵(rank=0的进程负责)，并分配内存存储矩阵A */
    if (my_rank == 0)
    {
        starttime = MPI_Wtime();
        fdA = fopen("dataIn.txt", "rb");
        /* 读取矩阵的阶数和阶数 */
        fread(&size, sizeof(int), 1, fdA);
        fread(&N, sizeof(int), 1, fdA);
        /* 判断是否为方阵，不是则报错并退出程序 */
        if (size != N)
        {
            puts("The input is error!");
            exit(0);
        }
        A = (float*)malloc(floatsize * size * size);
        B = (float*)malloc(floatsize * size * size);
        // Read the entire matrix at once
        fread(A, sizeof(int), (N) * (N), fdA);
        fclose(fdA);
    }
    /* 广播矩阵的阶数 */
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* 计算分块数，使得t*t=p */
    t = (int)sqrt(p);
    if (t > size)
        t = size;
    if (size % t != 0)
        for (;;)
        {
            t--;
            if (size % t == 0)
                break;
        }
    /* 重新计算进程数 */
    p = t * t;
    /* 计算每个子矩阵的大小 */
    m = size / t;

    /* 分配存储空间，a存储子矩阵，b用于临时存储 */
    a = (float*)malloc(floatsize * m * m);
    b = (float*)malloc(floatsize * m * m);

    if (a == NULL || b == NULL)
        printf("allocate space  fail!");

    /* 将数据分发到各个进程，每个进程存储一个子矩阵 */
    if (my_rank == 0)
    {
#pragma omp parallel for private(i, j) shared(m, a, A)
        for (i = 0;i < m;i++)
            for (j = 0;j < m;j++)
                a(i, j) = A(i, j);
    }

    /* 广播数据到其他进程，每个进程存储一个子矩阵 */
    if (my_rank == 0)
    {
        for (i = 1;i < p;i++)
        {
            v = i / t;                               /* 行块号 */
            u = i % t;                               /* 列块号 */

#pragma omp parallel for private(j, k) shared(v, u, m, b, A)
            for (j = v * m;j < (v + 1) * m;j++)
                for (k = u * m;k < (u + 1) * m;k++)
                    b((j % m), (k % m)) = A(j, k);

            /* 将子矩阵发送给其他进程 */
            MPI_Send(b, m * m, MPI_FLOAT, i, i, MPI_COMM_WORLD);
        }
    }
    else if (my_rank < p)                           /* 如果当前进程所在行数大于所在列数，则进行数据交换 */
        MPI_Recv(a, m * m, MPI_FLOAT, 0, my_rank, MPI_COMM_WORLD, &status);

    time1 = MPI_Wtime();

    /* 如果当前进程所在行数大于所在列数，则进行数据交换 */
    if ((my_rank / t) > (my_rank % t) && my_rank < p)
    {
        v = my_rank / t;                              /* �к� */
        u = my_rank % t;                              /* �к� */

        /* 发送当前子矩阵到对应的位置 */
        MPI_Send(a, m * m, MPI_FLOAT, (u * t + v), (u * t + v), MPI_COMM_WORLD);
        /* 接收对应位置的子矩阵 */
        MPI_Recv(a, m * m, MPI_FLOAT, (u * t + v), my_rank, MPI_COMM_WORLD, &status);
    }

    /* If the current process is in a row greater than its column and within the process count */
    if ((my_rank / t) < (my_rank % t) && my_rank < p)
    {
        v = my_rank / t;                              /* Row block index */
        u = my_rank % t;                              /* Column block index */
        /* Copy elements of the matrix 'a' to 'b' */
#pragma omp parallel for private(i, j) shared(m, a, b)
        for (i = 0;i < m;i++)
            for (j = 0;j < m;j++)
                b(i, j) = a(i, j);

        /* Receive the corresponding submatrix from the specified location */
        MPI_Recv(a, m * m, MPI_FLOAT, (u * t + v), my_rank, MPI_COMM_WORLD, &status);
        /* Send the copied submatrix to the specified location */
        MPI_Send(b, m * m, MPI_FLOAT, (u * t + v), (u * t + v), MPI_COMM_WORLD);
    }

    /* Transpose the submatrix 'a' */
#pragma omp parallel for private(i, j, temp)
    for (i = 1;i < m;i++)
        for (j = 0;j < i;j++)
        {
#pragma omp critical
            temp = a(i, j);
            a(i, j) = a(j, i);
            a(j, i) = temp;
        }

    /* If the current process is rank 0,
    store the transposed submatrix in matrix 'B' */
    if (my_rank == 0)
    {
#pragma omp parallel for private(i, j) shared(m, a, B)
        for (i = 0;i < m;i++)
            for (j = 0;j < m;j++)
                B(i, j) = a(i, j);
    }
    /* If the current process is rank 0, receive transposed submatrices
    from other processes and assemble them into matrix 'B' */
    if (my_rank == 0)
    {
        for (i = 1;i < p;i++)
        {
            /* Receive transposed submatrix */
            MPI_Recv(a, m * m, MPI_FLOAT, i, i, MPI_COMM_WORLD, &status);

            v = i / t;                                /* Row block index */
            u = i % t;                                /* Column block index */
#pragma omp parallel for private(j, k) shared(v, u, m, a, B)
            for (j = v * m;j < (v + 1) * m;j++)
                for (k = u * m;k < (u + 1) * m;k++)
                    B(j, k) = a((j % m), (k % m));        /* Assemble into matrix 'B' */
        }
    }
    else if (my_rank < p) /* For other processes, send the transposed submatrix to rank 0 */
        MPI_Send(a, m * m, MPI_FLOAT, 0, my_rank, MPI_COMM_WORLD);

    time2 = MPI_Wtime();
    /* If the current process is rank 0, print timing information */
    if (my_rank == 0)
    {
        printf("\n");
        printf("Whole running time    = %f seconds\n", time2 - starttime);
        printf("Distribute data time  = %f seconds\n", time1 - starttime);
        printf("Parallel compute time = %f seconds\n", time2 - time1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    Environment_Finalize(a, b);
    return(0);
}
