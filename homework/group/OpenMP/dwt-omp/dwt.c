#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <omp.h>

#define  MAX_PROCESSOR_NUM       12
#define  MAX_N                 2048*2048
#define  MAX_L                  100*100

int main(int argc, char* argv[])
{
    int c[MAX_N], d[MAX_N];
    int h[MAX_L], g[MAX_L];
    int pred[MAX_N], curr[MAX_N];
    int N, L;
    FILE* fin;
    int n;
    int offsetAddress;
    int t;
    int i, j, k;
    int size, rank;
    double  time1, time2, beginTime;
    beginTime = MPI_Wtime();
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        fin = fopen("dataIn.txt", "r");
        if (fin == NULL)
        {
            puts("Not find input data file");
            puts("Please create a file \"dataIn.txt\"");
            puts("<example for dataIn.txt> ");
            puts("10");
            puts("1 9 40 5 67 2 6 8 49 5");
            puts("1 3  4 5  6 7 8 8  0 4");
            puts("16");
            puts("1 0 4 5  5 6 7 8  5 4 0 3  4 6 1 2");
            MPI_Finalize();
            exit(-1);
        }
        else
        {
            /* read array h,g from dataIn.txt */
            // use fread to read array h, g from dataIn.txt
            fread(&L, sizeof(int), 1, fin);
            if ((L < 1) || (L > MAX_L))
            {
                puts("Input the length of array h is out of range!");
                exit(-1);
            }
            fread(h, sizeof(int), L, fin);
            fread(g, sizeof(int), L, fin);

            /* read array c from dataIn.txt */
            /* if  2^(t-1) < N < 2^t , then attach zero to expend N to 2^t */
            // use fread to read array pred from dataIn.txt
            fread(&N, sizeof(int), 1, fin);
            if ((N < 1) || (N > MAX_N))
            {
                puts("Input the length of array c is out of range!");
                exit(-1);
            }
            fread(pred, sizeof(int), N, fin);

            puts("One dimensional wavelet transform\'s input data from file dataIn.txt");
            printf("h[] and g[]'s length: %d\n", L);
            printf("c[]'s length: %d\n", N);
            printf("\n");
        }
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(h, L, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g, L, MPI_INT, 0, MPI_COMM_WORLD);
    time1 = MPI_Wtime();

    n = N;
    offsetAddress = 0;
    t = 0;
    while (n > 1)
    {
        MPI_Bcast(pred, n, MPI_INT, 0, MPI_COMM_WORLD);
        n /= 2;

        int jLength = n / size;
        int jStart = n % size + rank * jLength;

        if (rank == 0)
        {
            for (j = 0; j < jStart + jLength; j++)
            {
                curr[j] = 0;
                d[offsetAddress + j] = 0;
                for (k = 0; k < L; k++)
                {
                    curr[j] += h[k] * pred[(2 * j + k) % (2 * n)];
                    d[offsetAddress + j] += g[k] * pred[(2 * j + k) % (2 * n)];
                }
            }
            /* curr[0 .. jStart-1] ==> pred[0 .. jStart-1] */
#pragma omp parallel for private(j) shared(jStart,curr,pred)
            for (j = 0; j < jStart; j++)
                pred[j] = curr[j];
        }
        else
        {
            for (j = jStart; j < jStart + jLength; j++)
            {
                curr[j] = 0;
                d[offsetAddress + j] = 0;
                for (k = 0; k < L; k++)
                {
                    curr[j] += h[k] * pred[(2 * j + k) % (2 * n)];
                    d[offsetAddress + j] += g[k] * pred[(2 * j + k) % (2 * n)];
                }
            }
        }

        /* gather into &(pred+jStart) */
        MPI_Gather((curr + jStart), jLength, MPI_INT, (pred + n % size), jLength, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather((d + offsetAddress + jStart), jLength, MPI_INT, (d + offsetAddress + jStart), jLength, MPI_INT, 0, MPI_COMM_WORLD);

        /* store in c[] */
        if (rank == 0)
        {
            k = offsetAddress;
#pragma omp parallel for private(j,k) shared(offsetAddress,n,c,pred)
            for (j = 0; j < n; j++) {
                c[k] = pred[j];
                k++;
            }
        }

        if (rank == 0)
        {
            t++;
            printf("Running after loop  %d\n", t);
        }
        offsetAddress += n;

    }
    time2 = MPI_Wtime();
    /* If the current process is rank 0, print timing information */
    if (rank == 0)
    {
        printf("\n");
        printf("Whole running time    = %f seconds\n", time2 - beginTime);
        printf("Distribute data time  = %f seconds\n", time1 - beginTime);
        printf("Parallel compute time = %f seconds\n", time2 - time1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
