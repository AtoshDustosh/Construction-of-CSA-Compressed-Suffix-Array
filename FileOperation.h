#ifndef FILEOPERATION_H_INCLUDED
#define FILEOPERATION_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

int fnaDataSize(char* filePath);
int loadFnaData(char* filePath, int dataLength, char* T);


/**
 * Load ?.fna file and read data into an array T[] whose length is designated.
 *
 * NOTE: this method will append a character '$' to the end of the array.
 * NOTE.plus: this method is not responsible for requesting memory.
 *
 * @param filePath file path
 * @param arrayLength length of the DNA sequence plus a '$'
 * @param T[] array used for storing data (ending with '$')
 */
int loadFnaData(char* filePath, int arrayLength, char* T) {
    printf("Loading data ... \n");
    FILE* fp = fopen(filePath, "r");
    int ifData = 0;
    int fnaDataPointer = 0;

    if(fp != NULL && fnaDataPointer < arrayLength) {
        char ch = fgetc(fp);
        while(ch != EOF) {
            if(ch == '\n' && !ifData) {
                // according to format of *.fna file
                // data part is after the first line
                ifData = 1;
            }
            if(ifData && ch != '\n') {
                // if not '\n' and is already in the data part
                // store char into the array
                T[fnaDataPointer++] = lowerCase(ch);
            }
            ch = fgetc(fp);
        }
    } else {
        printf("failed to open file %s", filePath);
    }
    // add '$' to the end of the DNA seq
    T[fnaDataPointer++] = '$';
    free(fp);
    return fnaDataPointer;
}

/**
 * Get the size of data (DNA sequence) of a ?.fna file.
 *
 * @param filePath file path
 */
int fnaDataSize(char* filePath) {
    printf("Calculating size of the fna data file ... \n");
    FILE* fp = fopen(filePath, "r");
    int ifData = 0;   // if the file pointer is now in the data part
    int dataLength = 0;

    if(fp != NULL) {
        char ch = fgetc(fp);
        while(ch != EOF) {
            if(ch == '\n' && !ifData) {
                // according to format of *.fna file
                // data part is after the first line
                ifData = 1;
            }
            if(ifData && ch != '\n') {
                // if not '\n' and is already in the data part
                dataLength++;
            }
            ch = fgetc(fp);
        }
    } else {
        printf("failed to open file %s", filePath);
    }
    free(fp);
    return dataLength;
}

#endif // FILEOPERATION_H_INCLUDED
