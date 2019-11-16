#ifndef MERGESTEP_H_INCLUDED
#define MERGESTEP_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SABuildFunc.h"
#include "HelperFunction.h"

void mergeStepA(char* T, long* SA, long arrayLength, long partLength, long partNum, long partIndex);
void mergeStepB(char* T, long* SA, long partLength, long partNum, long partIndex, long* order);
void mergeStepC();




/**
 * Sort suffixes (suf1, suf2, ..., sufl) to find the lex-order among themselves.
 *
 *
 * @param T DNA sequence (plus a '$')
 * @param SA use SA to store the startIndexes of suffixes sorted by lex-order
 * @param arrayLength length of T
 * @param partLength length of a part (n/log2(n))
 * @param partNum number of parts
 * @param partIndex index of this part
 */
void mergeStepA(char* T, long* SA, long arrayLength, long partLength, long partNum, long partIndex) {
    long i = 0;
    long startIndex_i = (partIndex - 1) * partLength;
    long startIndex_apostrophe = partIndex * partLength;

    char* T_i = (char*)malloc(sizeof(char) * partLength);
    long* localSA = (long*)malloc(sizeof(long) * partLength);

//    printf("index\tT_i[]\tlocalSA[]\n");
    for(i = startIndex_i; i < startIndex_apostrophe; i++) {
        long locali = i - startIndex_i;
        T_i[locali] = T[i];
        localSA[locali] = locali;
//        printf("%ld\t%c\t%ld\n", locali, T_i[locali], localSA[locali]);
    }

    suffixArrayQuickSort(localSA, T_i, 0, partLength - 1);

    // copy the sorted suffixes to SA
    printf("sorted suffixes(%ld): \n", partIndex);
    printf("i\tch\tSA[]\tch_SA\n");
    for(i = startIndex_i; i < startIndex_apostrophe; i++){
        long locali = i - startIndex_i;
        SA[i] = localSA[locali];
        printf("%ld\t%c\t%ld\t%c\n", locali, T_i[locali], localSA[locali], T_i[localSA[locali]]);
    }
}



void mergeStepB(char* T, long* SA, long partLength, long partNum, long partIndex, long* order) {
    long i = 0;
    long startIndex_i = (partIndex - 1) * partLength;
    long startIndex_apostrophe = partIndex * partLength;

}



void mergeStepC() {

}


#endif // MERGESTEP_H_INCLUDED
