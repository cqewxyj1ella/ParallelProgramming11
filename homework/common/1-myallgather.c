#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int MPI_Allgather_my(
    const int* sendbuf, int sendcount, MPI_Datatype sendtype,
    int* recvbuf, int recvcount, MPI_Datatype recvtype,
    MPI_Comm comm) {
    // by default, regard send datatype as int
    int rank, nproc;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);
    MPI_Status status;
    // MPI_Gather, gather all data to rank 0
    if (rank == 0) {
        // copy send data from rank 0 to rank 0's buffer
        if (sendtype == MPI_INT) {
            recvbuf[0] = sendbuf[0];
        }
        // recv send data from other ranks to rank 0's buffer
        for (int i = 1; i < nproc; i++) {
            // It's very important that recv tag should align to send tag!!!!!!
            MPI_Recv(&(recvbuf[i]), recvcount, recvtype, i, i, comm, &status);
        }
    }
    else {
        // send original data to rank 0, tagged with own rank
        MPI_Send(sendbuf, sendcount, sendtype, 0, rank, comm);
    }

    // MPI_Scatter, scatter all data from rank 0 to others
    if (rank == 0) {
        for (int i = 1; i < nproc; i++) {
            MPI_Send(recvbuf, recvcount * nproc, recvtype, i, i, comm);
        }
    }
    else {
        // receive from rank 0, tagged with own rank
        MPI_Recv(recvbuf, recvcount * nproc, recvtype, 0, rank, comm, &status);
    }
}

int main(int argc, char* argv[]) {
    int rank, nproc;
    int isend, irecv[32];
    double start_time, end_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // test my allgather function
    isend = rank + 1;
    start_time = MPI_Wtime();
    MPI_Allgather_my(&isend, 1, MPI_INT, irecv, 1, MPI_INT,
        MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("total time of my_allgather is: %f\n", end_time - start_time);
    }

    // test original allgather function
    isend = rank + 1;
    start_time = MPI_Wtime();
    MPI_Allgather(&isend, 1, MPI_INT, irecv, 1, MPI_INT,
        MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("total time of allgather is: %f\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}
