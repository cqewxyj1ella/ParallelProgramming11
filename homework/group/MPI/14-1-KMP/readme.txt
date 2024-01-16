1.
gen_ped.c: 生成模式串的程序。
通过命令行参数传递文本串长度(Strlen)、模式串长度(Pedlen)、随机数种子(Seed)和模式文件名(Pattern_File)。
在给定的长度和种子下，生成一个模式串，然后将其重复拼接，最后存储到文件中。

gcc gen_ped.c -o gen_ped
./gen_ped 3 2 1 pattern.dat
3: total length of pattern string
2: prefix length of pattern string
1: random seed

2.
kmp.c: 使用MPI并行计算框架实现的KMP算法程序。
通过命令行参数传递文本串长度(m)和模式串长度(n)。
程序将文本串划分成若干个节点，每个节点执行KMP算法，并通过MPI通信协作。在最后输出每个节点的匹配结果。

mpicc kmp.c -o kmp -lm
mpirun -np 3 ./kmp 18 2 3
3: nproc
18: total length of text string
2: prefix length of pattern string
3: total length of pattern string


3. 输出结果: 根据你提供的示例，gen_ped生成了模式串"qmq"，而kmp程序通过MPI框架在3个节点上执行KMP算法。
输出结果显示了每个节点的文本串以及匹配结果。"+"表示匹配成功的位置，"-"表示匹配失败的位置。
The Text on node 0 is asasas .
The Text on node 1 is qmqmqm .
The Text on node 2 is ypypyp .
This is the match result on node 0
(0)  -
(1)  -
(2)  -
(3)  -
This is the match result on node 1
(4)  -
(5)  -
(6)  +
(7)  -
(8)  +
(9)  -
This is the match result on node 2
(10)  -
(11)  -
(12)  -
(13)  -
(14)  -
(15)  -