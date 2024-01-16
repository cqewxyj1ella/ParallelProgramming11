#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define  MAX_PROCESSOR_NUM       12
#define  MAX_N                 1024
#define  MAX_L                  100

void print(int b[], int start, int len)
{
    int end = start+len;
    int i;

    for(i = start; i < end; i ++)
    {
        printf("%d  ",b[i]);
    }
    printf("\n");

}

int main(int argc, char *argv[])
{
    int c[MAX_N], d[MAX_N];
    int h[MAX_L], g[MAX_L];
    int pred[MAX_N], curr[MAX_N];
    int N,L;
    FILE *fin;
    int n;
    int offsetAddress;
    int t;
    int i, j, k;
    int size, rank;
    double  transTime = 0 ,tempCurrentTime, beginTime;

    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(rank == 0)
    {
        fin = fopen("dataIn.txt","r");
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
            fscanf(fin,"%d",&L);
            if ((L < 1)||(L > MAX_L))
            {
                puts("Input the length of array h is out of range!");
                exit(-1);
            }
            for(i = 0; i < L; i ++) fscanf(fin,"%d", h+i);
            for(i = 0; i < L; i ++) fscanf(fin,"%d", g+i);

            /* read array c from dataIn.txt */
            fscanf(fin,"%d",&N);
            if ((N < 1)||(N > MAX_N))
            {
                puts("Input the length of array c is out of range!");
                exit(-1);
            }
            for(i = 0; i < N; i ++) fscanf(fin,"%d", pred+i);
            /* if  2^(t-1) < N < 2^t , then attach zero to expend N to 2^t */

            puts("One dimensional wavelet transform\'s input data from file dataIn.txt");
            printf("h[] : ");
            print(h, 0, L);
            printf("g[] : ");
            print(g, 0, L);
            printf("c[] : ");
            print(pred, 0, N);
            printf("\n");
        }
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(h, L, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(g, L, MPI_INT, 0, MPI_COMM_WORLD);

    n = N;
    offsetAddress = 0;
    t = 0;
    while(n > 1)
    {
        MPI_Bcast(pred, n, MPI_INT, 0, MPI_COMM_WORLD);
        n /=2;

        int jLength =  n/size;
        int jStart  =  n%size+rank*jLength;

        if (rank == 0)
        {
            for(j = 0; j < jStart+jLength; j ++)
            {
                curr[j] = 0;
                d[offsetAddress+j] = 0;
                for(k = 0; k < L; k ++)
                {
                    curr[j] += h[k]*pred[(2*j+k)%(2*n)];
                    d[offsetAddress+j] += g[k]*pred[(2*j+k)%(2*n)];
                }
            }
            /* curr[0 .. jStart-1] ==> pred[0 .. jStart-1] */
            for(j = 0; j < jStart; j ++) pred[j] = curr[j];
        }
        else
        {
            for(j = jStart; j < jStart+jLength; j ++)
            {
                curr[j] = 0;
                d[offsetAddress+j] = 0;
                for(k = 0; k < L; k ++)
                {
                    curr[j] += h[k]*pred[(2*j+k)%(2*n)];
                    d[offsetAddress+j] += g[k]*pred[(2*j+k)%(2*n)];
                }
            }
        }

        /* gather into &(pred+jStart) */
        MPI_Gather( (curr+jStart), jLength, MPI_INT, (pred+n%size), jLength, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather( (d+offsetAddress+jStart),jLength, MPI_INT, (d+offsetAddress+jStart), jLength, MPI_INT, 0, MPI_COMM_WORLD);

        /* store in c[] */
        if (rank == 0)
        {
            for(j = 0, k = offsetAddress; j < n;  j++, k ++) c[k] = pred[j];
        }

        if (rank == 0)
        {
            t ++;
            printf("Running after loop  %d\n", t);
            printf("c[] : ");
            print(c, offsetAddress, n);
            printf("d[] : ");
            print(d, offsetAddress, n);
            printf("\n");
        }
        offsetAddress += n;

    }

    MPI_Finalize();
}
