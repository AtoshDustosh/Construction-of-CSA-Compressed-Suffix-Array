#ifndef MERGESTEP_H_INCLUDED
#define MERGESTEP_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SABuildFunc.h"
#include "HelperFunction.h"

void mergeStepA(char* T, long* SA, long* SA_inverse, long arrayLength, long partLength,
                long partNum, long partIndex);
void mergeStepB(char* T, long* SA, long* Psi, long arrayLength, long partLength,
                long partNum, long partIndex, long* order);
void mergeStepC();




/**
 * Sort suffixes (suf1, suf2, ..., sufl) to find the lex-order among themselves.
 *
 * @param T DNA sequence (plus a '$')
 * @param SA use SA to store the startIndexes of suffixes sorted by lex-order
 * @param SA_inverse used to store inverse of SA_part of suffixes
 * @param arrayLength length of T
 * @param partLength length of a part (n/log2(n))
 * @param partNum number of parts
 * @param partIndex index of this part
 */
void mergeStepA(char* T, long* SA, long* SA_inverse, long arrayLength, long partLength,
                long partNum, long partIndex) {
    printf("Merge Step (a)\n");
    long i = 0;
    long bi_i = (partIndex - 1) * partLength;   // beginIndex_i
    long bi_apostrophe = partIndex * partLength;    // beginIndex_apostrophe

    char* T_i = (char*)malloc(sizeof(char) * partLength);
    long* localSA = (long*)malloc(sizeof(long) * partLength);
    long* localSA_inverse = (long*)malloc(sizeof(long) * partLength);

//    printf("index\tT_i[]\tlocalSA[]\n");
    for(i = bi_i; i < bi_apostrophe; i++) {
        long locali = i - bi_i;
        T_i[locali] = T[i];
        localSA[locali] = locali;
//        printf("%ld\t%c\t%ld\n", locali, T_i[locali], localSA[locali]);
    }

    suffixArrayQuickSort(localSA, T_i, 0, partLength - 1);
    inverseSAWhole(localSA, localSA_inverse, partLength);

    // copy sorted suffixes' indexes and their inverse to SA and SA_inverse
    printf("sorted suffixes(%ld): \n", partIndex);
    printf("begin_index: %ld\n", bi_i);
    printf("i\tch\tSA[]\tch_SA\tSA_inverse\n");
    for(i = bi_i; i < bi_apostrophe; i++) {
        long locali = i - bi_i;
        SA[i] = localSA[locali] + bi_i;
        SA_inverse[i] = localSA_inverse[locali] + bi_i;

        printf("%ld\t", i);
        printf("%c\t", T[i]);
        printf("%ld\t", SA[i]);
        printf("%c\t", T[SA[i]]);
        printf("%ld\t", SA_inverse[i]);
        printf("\n");
    }
}

/**
 * Calculate all order(suf_k, T') in an order k = l, l-1, ..., 1.
 * Note that suf[] is already sorted by lex-order.
 *
 * @param T DNA sequence (plus a '$')
 * @param SA use SA to store the startIndexes of suffixes sorted by lex-order
 * @param Psi Psi array of T
 * @param arrayLength length of T
 * @param partLength length of a part (n/log2(n))
 * @param partNum number of parts
 * @param partIndex index of this part
 * @param order order[] that stores result of func order(suf_k, T'), where suf_k is the k-th
 *      longest suffix of T_i and T' is combination of T_(i+1)...T_([n/l|+]).
 */
void mergeStepB(char* T, long* SA, long* Psi, long arrayLength, long partLength,
                long partNum, long partIndex, long* order) {
    printf("Merge Step (b)\n");
    long i = 0;
    long bi_i = (partIndex - 1) * partLength;   // beginIndex_i
    long bi_apostrophe = partIndex * partLength;    // beginIndex_apostrophe

    // T_apostrophe - start from bi_apostrophe and ends at arrayLength
    // note that lc and rc are not actual bounds of the character c
    long prevOrderValue = SA[bi_apostrophe];
    for(i = bi_apostrophe - 1; i >= bi_i; i--) {
        char c = T[i];  // the character in the formula
        long lc = bi_apostrophe;
        long rc = arrayLength - 1;
        long orderValue = 0;
        CSABinaryBoundSearch(T, SA, c, &lc, &rc);
        // T[SA[lc]] ~ T[SA[rc]] represents the field of c
        // implement of condition ¦×[b] -> ¦×[SA[b]], lc <= b <= rc
        printf("(%c) -> lc: %ld, rc: %ld", c, lc, rc);
        printf("\n");

        if(lc > rc) {
            // TODO I have doubt about this. I think it should be "orderValue = lc - bi_apostrophe;"
            orderValue = lc - 1 - bi_apostrophe;
        } else {
            // find the max b that satisfies condition that for any order(cX, T'), Psi_T'[b] <= order(X, T')
            long max_b = 0;
            CSABinarySearchOrderValue(SA, Psi, lc, rc, prevOrderValue, &max_b);
            orderValue = max_b - bi_apostrophe;
        }

//        printf("%ld, orderValue: %ld\n", i - bi_i, orderValue);

        order[i - bi_i] = orderValue;
        prevOrderValue = orderValue;
    }


    for(i = bi_apostrophe; i < arrayLength; i++) {
        printf("%c", T[SA[i]]);
    }

    printf("\n");

}



void mergeStepC() {

}


#endif // MERGESTEP_H_INCLUDED
