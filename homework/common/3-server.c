#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define MAX_RANDOM 10

void parameterServer(int rank, int nproc, int P) {
    int size = nproc - P;
    int* received_data = (int*)malloc(size * sizeof(int));
    memset(received_data, 0, sizeof(int) * size);
    int count = 0;
    // receive data from my workers
    for (int i = 0; i < nproc; i++) {
        if (i >= P && i % P == rank) {
            // this worker is my worker to serve
            int tmp = 0;
            MPI_Recv(&tmp, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            received_data[count] = tmp;
            count++;
        }
    }
    // calculate mean
    int sum = 0;
    for (int i = 0; i < count; i++) {
        sum += received_data[i];
    }
    double average = (double)sum / count;
    printf("Parameter %d calculated average: %lf\n", rank, average);
    // send mean to my workers
    for (int i = 0; i < nproc; i++) {
        if (i >= P && i % P == rank) {
            // this worker is my worker to serve
            MPI_Send(&average, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
    }
    free(received_data);
}

void worker(int rank, int P) {
    // generate data
    int data = rand() % MAX_RANDOM;
    // send data to server
    int server = rank % P; // my parameter server's rank id
    MPI_Send(&data, 1, MPI_INT, server, 0, MPI_COMM_WORLD);
    // receive mean value from server
    double average;
    MPI_Recv(&average, 1, MPI_DOUBLE, server, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Worker %d received average: %lf\n", rank, average);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    srand(rank);

    int P = 3; // Number of parameter server processes
    int Q = nproc - P; // Number of worker processes

    if (rank < P) {
        parameterServer(rank, nproc, P);
    }
    else {
        worker(rank, P);
    }

    MPI_Finalize();

    return 0;
}
