#ifndef FILEOPERATION_H_INCLUDED
#define FILEOPERATION_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

long fnaDataSize(char* filePath);

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
    return dataLength;
}

#endif // FILEOPERATION_H_INCLUDED
