#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    if (argc < 4) {
        printf("Please provide three arguments.\n");
        return 1;
    }

    int size1 = atoi(argv[1]);
    int size2 = atoi(argv[2]);
    int seed = atoi(argv[3]);
    srand(seed);

    int* array1 = (int*)malloc(size1 * sizeof(int));
    int* array2 = (int*)malloc(size1 * sizeof(int));
    int* array3 = (int*)malloc(size2 * sizeof(int));

    // Generate arrays
    for (int i = 0; i < size1; i++) {
        array1[i] = rand() % 5;
        array2[i] = rand() % 8;
    }

    for (int i = 0; i < size2; i++) {
        array3[i] = rand() % 10;
    }

    // Store sizes and data into file
    FILE* file = fopen("dataIn.txt", "wb");
    if (file == NULL) {
        printf("Failed to open file.\n");
        return 1;
    }

    fwrite(&size1, sizeof(int), 1, file);
    fwrite(array1, sizeof(int), size1, file);
    fwrite(array2, sizeof(int), size1, file);
    fwrite(&size2, sizeof(int), 1, file);
    fwrite(array3, sizeof(int), size2, file);

    fclose(file);

    // Clean up
    free(array1);
    free(array2);
    free(array3);

    return 0;
}
