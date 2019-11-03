#include <stdio.h>
#include <stdlib.h>
#include "headfiletest.h"
#include "SABuildFunc.h"
#include "HelperFunction.h"
#include "FileOperation.h"

#define lineLength 70

char* FILEPATH = "NC_008253.fna";   // file path
long DATALENGTH = 0; // length of data (header not included)

// TODO rename these global variables ... T[] and SA[] is too~ simple and confusing
char DNASeq[];   // DNA sequence (?.fna file dat)

void testSet();
void readAndPrint();
// TODO add the lost claim of functions

void main() {

//    readAndPrint();
    long length = fnaDataSize(FILEPATH);
    printf("another method - data length: %ld\n", length);
}

// TODO add new functions




/////////////////////////////////// Test Functions Below ///////////////////////////////////

/**
 * Read all data in ?.fna file and print to console. (.fna file with specific data line length)
 */
void readAndPrint() {
    FILE* fp = fopen(FILEPATH, "r");
    int dataPart = 0;
    DATALENGTH = 0;
    if(fp != NULL) {
        char ch = fgetc(fp);
        char buffer[lineLength];
        int bufCount = 0;
        while(ch != EOF) {
            if(ch == '\n' && !dataPart) {
                dataPart = 1;
                printf(" #(... get into data part)");
            }
            if(dataPart && ch != '\n') {
//                printf("store character: %c\n", ch);
                buffer[bufCount++] = ch;
                DATALENGTH++;
            } else {
                printf("%c", ch);
            }
            if(bufCount != 0 && (bufCount % lineLength) == 0) {
                // print buffer content when buffer is full and reset the bufPointer
                printf("%s", buffer);
                bufCount = 0;
            }
            ch = fgetc(fp);
        }
    } else {
        printf("failed to open file %s", FILEPATH);
    }
    printf("dataLength: %ld\n", DATALENGTH);
}

/**
 * This is simply a debug function, mainly used when developing this program.
 * Or if you just want to see the effect or progress of some function, you can disable the comment mark,
 * and execute this function.
 */
void testSet() {
//    _CLanguageReview();
//    _quickSortTest();
//    _myStrLengthTest();
//    _suffixArrayQuickSortTest();
//    _compareSuffixTest();
//    _inverseSATest();
//    _psiArrayBuildTest();
}


