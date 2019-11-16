#ifndef BASESTEP_H_INCLUDED
#define BASESTEP_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SABuildFunc.h"
#include "FileOperation.h"
#include "HelperFunction.h"

void baseStep(char* filePath, char* T, long* SA, long* SA_inverse, long* Psi,
              long arrayLength, long partLength, long partNum);



void baseStep(char* filePath, char* T, long* SA, long* SA_inverse, long* Psi,
              long arrayLength, long partLength, long partNum) {
    long bi = (partNum - 1) * partLength;
    long localLength = arrayLength - bi;
    long i = 0;
    char* localT = NULL;
    long* localSA = NULL;
    long* localSA_inverse = NULL;
    long* localPsi = NULL;

    // load the DNA sequence
    loadFnaData(filePath, arrayLength, T);
    printf("The whole DNA sequence plus a \'$\':\n");
    for(i = 0; i < arrayLength; i++) {
        printf("%c", T[i]);
    }

    printf("\n\nBeginning at index %ld, the last part T_(n/l) is:\n", bi);
    for(i = bi; i < arrayLength; i++) {
        printf("%c", T[i]);
    }

    printf("\n\n");

    localT = (char*)malloc(sizeof(char) * localLength);
    localSA = (long*)malloc(sizeof(long) * localLength);
    localSA_inverse = (long*)malloc(sizeof(long) * localLength);
    localPsi = (long*)malloc(sizeof(long) * localLength);

    for(i = bi; i < arrayLength; i++) {
        localT[i - bi] = T[i];
        localSA[i - bi] = i - bi;
    }

    suffixArrayQuickSort(localSA, localT, 0, localLength - 1);
//    for(i = 0; i < localLength; i++){
//        printf("%c\t%ld\n", localT[i], localSA[i]);
//    }
    inverseSAWhole(localSA, localSA_inverse, localLength);
    psiArrayBuildWhole(localSA, localSA_inverse, localPsi, localLength);

    // copy the results to the original arrays
    for(i = bi; i < arrayLength; i++) {
        long locali = i - bi;
        SA[i] = localSA[locali] + bi;
        SA_inverse[i] = localSA_inverse[locali] + bi;
        Psi[i] = localPsi[locali] + bi;
    }

    // print the sorted first increment
    printf("index\tch\tSA\tch_SA\tSA_inv\tPsi\tch_Psi\n");
//    for(i = 0; i < localLength; i++) {  // simple version
//        printf("%ld\t", i);
//        printf("%c\t", localT[i]);
//        printf("%ld\t", localSA[i]);
//        printf("%c\t", localT[localSA[i]]);
//        printf("%ld\t", localSA_inverse[i]);
//        printf("%ld\t", localPsi[i]);
//        printf("%c\t", localT[localPsi[i]]);
//        printf("\n");
//    }
    for(i = bi; i < arrayLength; i++) {  // complicated version
        printf("%ld\t", i);
        printf("%c\t", T[i]);
        printf("%ld\t", SA[i]);
        printf("%c\t", T[SA[i]]);
        printf("%ld\t", SA_inverse[i]);
        printf("%ld\t", Psi[i]);
        printf("%c\t", T[Psi[i]]);
        printf("\n");
    }

}

#endif // BASESTEP_H_INCLUDED
