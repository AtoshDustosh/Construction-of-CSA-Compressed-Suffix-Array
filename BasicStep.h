#ifndef BASESTEP_H_INCLUDED
#define BASESTEP_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SABuildFunc.h"
#include "FileOperation.h"
#include "HelperFunction.h"

void baseStep(char* filePath, char* T, int* SA, int* SA_inverse, int* Psi,
              int arrayLength, int partLength, int partNum);

void baseStep(char* filePath, char* T, int* SA, int* SA_inverse, int* Psi,
              int arrayLength, int partLength, int partNum) {
    int bi = (partNum - 1) * partLength;
    int localLength = arrayLength - bi;
    int i = 0;
    char* localT = NULL;
    int* localSA = NULL;
    int* localSA_inverse = NULL;
    int* localPsi = NULL;

    // load the DNA sequence
    loadFnaData(filePath, arrayLength, T);
    printf("The whole DNA sequence plus a \'$\':\n");
//    for(i = 0; i < arrayLength; i++) {
//        printf("%c", T[i]);
//    }

    localT = (char*)malloc(sizeof(char) * localLength);
    localSA = (int*)malloc(sizeof(int) * localLength);
    localSA_inverse = (int*)malloc(sizeof(int) * localLength);
    localPsi = (int*)malloc(sizeof(int) * localLength);

    printf("\n\nBeginning at index %d, the last part T_(n/l) is:\n", bi);
    for(i = bi; i < arrayLength; i++) {
//        printf("%c", T[i]);
        localT[i - bi] = T[i];
        localSA[i - bi] = i - bi;
    }
    printf("\n");
    printf("Length of the last part: %d\n", localLength);
    printf("\n");


    suffixArrayQuickSort(localSA, localT, 0, localLength - 1);
//    for(i = 0; i < localLength; i++){
//        printf("%c\t%d\n", localT[i], localSA[i]);
//    }
    inverseSAWhole(localSA, localSA_inverse, localLength);
    psiArrayBuildWhole(localSA, localSA_inverse, localPsi, localLength);

    // copy the results to the original arrays
    for(i = bi; i < arrayLength; i++) {
        int locali = i - bi;
        SA[i] = localSA[locali] + bi;
        SA_inverse[i] = localSA_inverse[locali] + bi;
        Psi[i] = localPsi[locali] + bi;
    }

    // print the sorted first increment
//    printf("index\tch\tSA\tch_SA\tSA_inv\tPsi\n");
//    for(i = bi; i < arrayLength; i++) {  // complicated version
//        printf("%d\t", i - bi);
//        printf("%c\t", T[i]);
//        printf("%d\t", SA[i] - bi);
//        printf("%c\t", T[SA[i]]);
//        printf("%d\t", SA_inverse[i] - bi);
//        printf("%d\t", Psi[i] - bi);
//        printf("\n");
//    }

    free(localSA);
    free(localSA_inverse);
    free(localT);
    free(localPsi);
}

#endif // BASESTEP_H_INCLUDED
