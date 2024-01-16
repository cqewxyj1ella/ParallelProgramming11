#include <malloc.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#define  MAX(m,n)    (m>n?m:n)

typedef struct {
	int pedlen;
	int psuffixlen;
	int pednum;
}pntype;


/* modified next map function */
void Next(char* W, int patlen, int* nextval, pntype* pped)
{
	int i, j, plen;
	int* next;

	if ((next = (int*)malloc((patlen + 1) * sizeof(int))) == NULL) {
		printf("no enough memory\n");
		exit(1);
	}

	/*计算next和nextval*/
	// 初始化第一个元素
	next[0] = nextval[0] = -1;
	j = 1;
	while (j <= patlen) {
		i = next[j - 1];
		// 在不匹配的情况下，回溯 i 直到找到匹配的前缀
		while (i != (-1) && W[i] != W[j - 1]) i = next[i];
		next[j] = i + 1;
		// 计算 nextval[j]，这是为了提高匹配效率
		if (j != patlen) {
			if (W[j] != W[i + 1])
				nextval[j] = i + 1;
			else
				nextval[j] = nextval[i + 1];
		}
		j++;
	}

	pped->pedlen = patlen - next[patlen];
	pped->pednum = (int)(patlen / pped->pedlen);
	pped->psuffixlen = patlen % pped->pedlen;

	free(next);
}

/* KMP算法实现 */
void kmp(char* T, char* W, int textlen, int patlen, int* nextval, pntype* pped, int prefix_flag, int matched_num, int* match, int* prefixlen)
{
	int i, j;

	i = matched_num;
	j = matched_num;
	// 主循环，处理文本串T和模式串W的匹配过程
	while (i < textlen)
	{
		// 如果启用了前缀优化，并且剩余匹配长度小于当前位置到文本末尾的长度，提前结束
		if ((prefix_flag == 1) && ((patlen - j) > (textlen - i))) {
			break;
		}
		// 在不匹配的情况下，回溯j直到找到匹配的前缀
		while (j != (-1) && W[j] != T[i])
			j = nextval[j];
		// 判断是否找到完整的匹配
		if (j == (patlen - 1)) {
			// 标记匹配位置
			match[i - (patlen - 1)] = 1;
			// 处理完整匹配后，如果有周期性的重复模式，跳过一定长度
			if (pped->pednum + pped->psuffixlen == 1)
				j = -1;
			else
				j = patlen - 1 - pped->pedlen;
		}
		j++;
		i++;
	}
	(*prefixlen) = j;
}

/* 重新构建模式串信息及其对应的next数组 */
void Rebuild_info(int patlen, pntype* pped, int* nextval, char* W)
{
	int i;
	// 如果周期模式的重复次数为1，则直接将模式串的后缀部分复制到模式串的末尾
	if (pped->pednum == 1)
		memcpy(W + pped->pedlen, W, pped->psuffixlen);
	else {
		// 否则，进行周期性的重复构建
		memcpy(W + pped->pedlen, W, pped->pedlen);
		for (i = 3; i <= pped->pednum; i++) {
			memcpy(W + (i - 1) * pped->pedlen, W, pped->pedlen);
			memcpy(nextval + (i - 1) * pped->pedlen, nextval + pped->pedlen, pped->pedlen * sizeof(int));
		}
		// 处理最后可能不完整的周期
		if (pped->psuffixlen != 0) {
			memcpy(W + (i - 1) * pped->pedlen, W, pped->psuffixlen);
			memcpy(nextval + (i - 1) * pped->pedlen, nextval + pped->pedlen, pped->psuffixlen * sizeof(int));
		}
	}
}

/* 生成随机字符串 */
void gen_string(int strlen, int pedlen, char* string, int seed)
{
	int suffixlen, num, i, j;

	srand(seed);
	// 生成模式串的前缀部分
	for (i = 0;i < pedlen;i++) {
		num = rand() % 26;
		string[i] = 'a' + num;
	}
	// 复制前缀以构建完整字符串
	for (j = 1;j < (int)(strlen / pedlen);j++)
		strncpy(string + j * pedlen, string, pedlen);
	// 处理可能不完整的后缀
	if ((suffixlen = strlen % pedlen) != 0)
		strncpy(string + j * pedlen, string, suffixlen);
}

/* 读取文件内容并获取文件大小 */
char* GetFile(char* filename, int number)
{
	FILE* fp;
	struct stat statbuf;
	if ((fp = fopen(filename, "r")) == NULL) {
		printf("Error open file %s\n", filename);
		exit(0);
	}
	fstat(fileno(fp), &statbuf);
	if (statbuf.st_size != number) {
		printf("pattern string length isn't right\n");
		exit(0);
	}

	char* place = NULL;
	if (((place) = (char*)malloc(sizeof(char) * (number + 1))) == NULL) {
		printf("Error alloc memory\n");
		exit(1);
	}
	size_t len = fread((place), sizeof(char), number, fp);
	if (len != number) {
		printf("Error in reading num\n");
		exit(0);
	}
	place[number] = '\0';
	printf("place: %s\n", place);
	fclose(fp);
	// (*number) = statbuf.st_size;
	return place;
}

/* 将节点的文本信息写入文件 */
void PrintFile_info(char* filename, char* T, int id)
{
	FILE* fp;
	int i;

	if ((fp = fopen(filename, "a")) == NULL) {
		printf("Error open file %s\n", filename);
		exit(0);
	}

	fprintf(fp, "The Text on node %d is %s .\n", id, T);

	fclose(fp);
}


/* 打印匹配结果到文件 */
void PrintFile_res(char* filename, int* t, int len, int init, int id)
{
	FILE* fp;
	int i;

	if ((fp = fopen(filename, "a")) == NULL) {
		printf("Error open file %s\n", filename);
		exit(0);
	}

	fprintf(fp, "This is the match result on node %d \n", id);
	for (i = 0; i <= len - 1; i++)
		if (t[i] == 1)
			fprintf(fp, "(%d)  +\n", i + init);
		else
			fprintf(fp, "(%d)  -\n", i + init);
	fclose(fp);
}

void main(int argc, char* argv[])
{
	char* T, * W;
	int* nextval, * match;
	int	textlen, patlen, pedlen, nextlen_send;
	pntype pped;
	int	i, myid, numprocs, prefixlen, ready;
	MPI_Status  status;
	double mpist, mpied;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	mpist = MPI_Wtime();
	nextlen_send = 0;
	ready = 1;
	/* 读取命令行参数并生成文本串 */
	textlen = atoi(argv[1]);
	textlen = textlen / numprocs;
	pedlen = atoi(argv[2]);
	patlen = atoi(argv[3]);
	if ((T = (char*)malloc(textlen * sizeof(char))) == NULL) {
		printf("no enough memory\n");
		exit(1);
	}
	memset(T, 0, textlen * sizeof(char));
	gen_string(textlen, pedlen, T, myid);
	// 主节点打印文本信息，并通知下一个节点
	if (myid == 0) {
		PrintFile_info("match_result", T, myid);
		if (numprocs > 1)
			MPI_Send(&ready, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	}
	else {
		MPI_Recv(&ready, 1, MPI_INT, myid - 1, myid - 1, MPI_COMM_WORLD, &status);
		PrintFile_info("match_result", T, myid);
		if (myid != numprocs - 1)
			MPI_Send(&ready, 1, MPI_INT, myid + 1, myid, MPI_COMM_WORLD);
	}

	printf("\n");

	if ((match = (int*)malloc(textlen * sizeof(int))) == NULL) {
		printf("no enough memory\n");
		exit(1);
	}

	/* 主节点读取模式串并初始化nextval数组 */
	if (myid == 0) {
		printf("processor num = %d \n", numprocs);
		printf("textlen = %d\n", textlen * numprocs);
		char* get_string = GetFile("pattern.dat", patlen);
		W = (char*)malloc(patlen * sizeof(char));
		strcpy(W, get_string);
		// W = "nwnedfkih";
		// W = "nwn";
		// patlen = strlen(W);
		printf("pattern : = %s, length is %d\n", W, patlen);
		// printf("patlen= %d\n", patlen);

		if ((nextval = (int*)malloc(patlen * sizeof(int))) == NULL) {
			printf("no enough memory\n");
			exit(1);
		}
		/* 计算模式串的next数组 */
		Next(W, patlen, nextval, &pped);
		if (numprocs > 1) {
			if (pped.pednum == 1)
				nextlen_send = patlen;
			else
				nextlen_send = pped.pedlen * 2;
		}
	}

	/* 广播模式串相关信息至其他节点 */
	if (numprocs > 1) {
		MPI_Bcast(&patlen, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (myid != 0)
			if (((nextval = (int*)malloc(patlen * sizeof(int))) == NULL)
				|| ((W = (char*)malloc(patlen * sizeof(char))) == NULL)) {
				printf("no enough memory\n");
				exit(1);
			}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&pped, 3, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&nextlen_send, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(nextval, nextlen_send, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(W, pped.pedlen, MPI_CHAR, 0, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/* 在单节点情况下直接执行KMP算法 */
	if (numprocs == 1) {
		kmp(T, W, textlen, patlen, nextval, &pped, 1, 0, match + patlen - 1, &prefixlen);
	}
	else {
		/* 在多节点情况下，首先进行模式串信息的重构 */
		if (myid != 0)
			Rebuild_info(patlen, &pped, nextval, W);
		/* 在中间节点和最后一个节点执行KMP算法 */
		if (myid != numprocs - 1)
			kmp(T, W, textlen, patlen, nextval, &pped, 0, 0, match + patlen - 1, &prefixlen);
		else
			kmp(T, W, textlen, patlen, nextval, &pped, 1, 0, match + patlen - 1, &prefixlen);

		MPI_Barrier(MPI_COMM_WORLD);

		/* 在多节点情况下，处理重构后的信息，确保匹配正确 */
		if (myid < numprocs - 1)
			MPI_Send(&prefixlen, 1, MPI_INT, myid + 1, 99, MPI_COMM_WORLD);

		if (myid > 0)
			MPI_Recv(&prefixlen, 1, MPI_INT, myid - 1, 99, MPI_COMM_WORLD, &status);

		MPI_Barrier(MPI_COMM_WORLD);
		/* 在多节点情况下，对匹配结果进行修正 */
		if ((myid > 0) && (prefixlen != 0))
			kmp(T - prefixlen, W, prefixlen + patlen - 1, patlen, nextval, &pped, 1, prefixlen, match + patlen - 1 - prefixlen, &prefixlen);

		MPI_Barrier(MPI_COMM_WORLD);
	}

	/* 打印匹配结果 */
	mpied = MPI_Wtime();

	if (myid == 0)
	{
		printf("TIME COST EVALUATED BY MPI: %f secs\n", mpied - mpist);
	}

	if (myid == 0) {
		PrintFile_res("match_result", match + patlen - 1, textlen - patlen + 1, 0, myid);
		if (numprocs > 1)
			MPI_Send(&ready, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	}
	else {
		MPI_Recv(&ready, 1, MPI_INT, myid - 1, myid - 1, MPI_COMM_WORLD, &status);
		PrintFile_res("match_result", match, textlen, myid * textlen - patlen + 1, myid);
		if (myid != numprocs - 1)
			MPI_Send(&ready, 1, MPI_INT, myid + 1, myid, MPI_COMM_WORLD);
	}

	free(T);
	free(W);
	free(nextval);
	MPI_Finalize();
}
