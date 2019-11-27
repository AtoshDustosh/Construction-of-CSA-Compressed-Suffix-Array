#ifndef FILEOPERATION_H_INCLUDED
#define FILEOPERATION_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

long fnaDataSize(char* filePath);
long loadFnaData(char* filePath, long dataLength, char* T);


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
long loadFnaData(char* filePath, long arrayLength, char* T) {
    FILE* fp = fopen(filePath, "r");
    int ifData = 0;
    long fnaDataPointer = 0;

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
long fnaDataSize(char* filePath) {
    FILE* fp = fopen(filePath, "r");
    int ifData = 0;   // if the file pointer is now in the data part
    long dataLength = 0;

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
