#include <stdio.h>
#include <stdlib.h>
#include "headfiletest.h"
#include "SABuildFunc.h"
#include "HelperFunction.h"
#include "FileOperation.h"

#define LINELENGTH 70

/*
 * Global variables.
 */
char* FILEPATH = "testdata_7000.txt";   // file path
long ARRAYLENGTH = 0; // length of DNA sequence array (ending with '$')
long PARTLENGTH = 0; // part length of DNA sequence array ( PARTLENGTH = log2(ARRAYLENGTH))
char* T = NULL; // DNA sequence of (A,C,G,T)

/*
 * Functions.
 */
void testProcess();
void testSet();
void readAndPrint();

// TODO add the lost claim of functions

int main() {
//    readAndPrint();
    testProcess();
    return 0;
}

// TODO add new functions

/**
 * Test steps of construction of CSA.
 *
 * <note> bug detected - array directly defined cannot be too big.
 * To solve this problem, use (long*)malloc(sizeof(long)*ARRAYLENGTH).
 */
void testProcess(){
    ARRAYLENGTH = fnaDataSize(FILEPATH);
    printf("data length: %ld\n", ARRAYLENGTH);

    // build T[] - DNA sequence array
    ARRAYLENGTH++; // get ready to add character '$' to the end of the DNA sequence
    char* temp = (char*)malloc(sizeof(char)*ARRAYLENGTH);
    T = temp;
    loadFnaData(FILEPATH, ARRAYLENGTH, temp);

    printf("DNA sequence - T[]: \n");
    for(long i = 0; i < ARRAYLENGTH; i++){
        //printf("%c", T[i]);
    }
    printf("\n");

    // build SA[] - suffix array
    long* SA = (int*)malloc(sizeof(long)*ARRAYLENGTH);
    for(long i = 0; i < ARRAYLENGTH; i++){
        SA[i] = i;
    }
    suffixArrayQuickSort(SA, T, 0, ARRAYLENGTH-1);
    printf("Suffix array - SA[]: \n");
    for(long i = 0; i < ARRAYLENGTH; i++){
        //printf("%ld\t%ld\t%c\n", i, SA[i], T[SA[i]]);
    }
    printf("\n");

    // build Psi[] - ... Psi array (I don't know how to describe it)
    long* SA_inverse = (int*)malloc(sizeof(long)*ARRAYLENGTH);
    long* Psi = (int*)malloc(sizeof(long)*ARRAYLENGTH);
    printf("Inverse suffix array - SA_inverse[]\n");
    inverseSA(SA, SA_inverse, ARRAYLENGTH);
    printf("Psi array - Psi[]: \n");
    psiArrayBuild(SA, SA_inverse, Psi, ARRAYLENGTH);
    for(long i = 0; i < ARRAYLENGTH; i++){
        printf("%ld\t%ld\t%c\n", i, Psi[i], T[Psi[i]]);
    }
    printf("\n");

    free(T);
    free(SA);
    free(SA_inverse);
    free(Psi);
}


/////////////////////////////////// Test Functions Below ///////////////////////////////////

/**
 * Read all data in ?.fna file and print to console. (.fna file with specific data line length)
 */
void readAndPrint() {
    FILE* fp = fopen(FILEPATH, "r");
    int dataPart = 0;
    ARRAYLENGTH = 0;
    if(fp != NULL) {
        char ch = fgetc(fp);
        char buffer[LINELENGTH];
        int bufCount = 0;
        while(ch != EOF) {
            if(ch == '\n' && !dataPart) {
                dataPart = 1;
                printf(" #(... get into data part)");
            }
            if(dataPart && ch != '\n') {
//                printf("store character: %c\n", ch);
                buffer[bufCount++] = ch;
                ARRAYLENGTH++;
            } else {
                printf("%c", ch);
            }
            if(bufCount != 0 && (bufCount % LINELENGTH) == 0) {
                // print buffer content when buffer is full and reset the bufPointer
                printf("%s", buffer);
                bufCount = 0;
            }
            ch = fgetc(fp);
        }
    } else {
        printf("failed to open file %s", FILEPATH);
    }
    printf("dataLength: %ld\n", ARRAYLENGTH);
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


