// generate a random matrix and save to file
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <string.h>
#include <time.h>

int N;
void init_matrix(int** A) {
    // A and B are matrix of size nxn
    srand((unsigned int)time(0));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = rand() % 10;
        }
    }
}

int main(int argc, char** argv)
{
    N = atoi(argv[1]); // read N from command line
    printf("%d\n", N);
    int** A = (int**)malloc(sizeof(int*) * N);
    for (int i = 0; i < N;i++) {
        A[i] = (int*)malloc(sizeof(int) * N);
        memset(A[i], 0, sizeof(int) * N);
    }
    // generate matrix
    init_matrix(A);
    printf("%d %d\n", A[0][0], N);
    // write matrix to file
    FILE* fp;
    if ((fp = fopen("dataIn.txt", "wb")) != NULL) {
        fwrite(&N, sizeof(int), 1, fp);
        fwrite(&N, sizeof(int), 1, fp);
        // Write each line of the matrix
        // for (int i = 0; i < N; i++) {
        //     for (int j = 0; j < N; j++) {
        //         fprintf(fp, "%d ", A[i][j]);
        //     }
        //     fprintf(fp, "\n");
        // }
        // Write the entire matrix at once
        fwrite(A, sizeof(int), N * N, fp);
        fclose(fp);
    }
    else {
        printf("file open error\n");
        exit(1);
    }

    return 0;
}