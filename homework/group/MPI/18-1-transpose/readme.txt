1.compile
mpicc transpose.c -o transpose -lm
gcc -o gen_matrix gen_matrix.c

2.run
mpirun -np 4 transpose
./gen_matrix 3

3.result
Input of file "dataIn.txt"
3       3
1.000000        2.000000        3.000000
3.000000        4.000000        5.000000
5.000000        6.000000        7.000000

Output of Matrix AT
1.000000        3.000000        5.000000
2.000000        4.000000        6.000000
3.000000        5.000000        7.000000

Whole running time    = 0.012299 seconds
Distribute data time  = 0.011884 seconds
Parallel compute time = 0.000415 seconds
