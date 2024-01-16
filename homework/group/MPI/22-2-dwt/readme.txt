1. compile:
mpicc dwt.c -o dwt

2. run:
mpirun -np 4 dwt

3. result:
One dimensional wavelet transform's input data from file dataIn.txt
h[] : -1  0  1  2
g[] : 1  1  0  3
c[] : 1  4  2  5  4  6  7  8

Running after loop  1
c[] : 11  14  19  2
d[] : 20  25  34  27

Running after loop  2
c[] : 12  20
d[] : 31  63

Running after loop  3
c[] : 40
d[] : 92

