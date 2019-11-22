#ifndef MERGESTEP_H_INCLUDED
#define MERGESTEP_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SABuildFunc.h"
#include "HelperFunction.h"

void mergeStepA(char* T, long* SA, long arrayLength, long partLength, long partNum, long partIndex);
void mergeStepB(char* T, long* SA, long* Psi, long arrayLength, long partLength,
                long partNum, long partIndex, long* order);
void mergeStepC(char* T, long* SA, long* Psi, long arrayLength, long partLength,
                long partNum, long partIndex, long* order);

void _fgpsiFuncTest();


/**
 * Sort suffixes (suf1, suf2, ..., sufl) to find the lex-order among themselves.
 *
 * @param T DNA sequence (plus a '$')
 * @param SA use SA to store the startIndexes of suffixes sorted by lex-order
 * @param arrayLength length of T
 * @param partLength length of a part (n/log2(n))
 * @param partNum number of parts
 * @param partIndex index of this part
 */
void mergeStepA(char* T, long* SA, long arrayLength, long partLength, long partNum,
                long partIndex) {
    printf("Merge Step (a)\n");
    long i = 0;
    long bi_i = (partIndex - 1) * partLength;   // beginIndex_i
    long bi_apostrophe = partIndex * partLength;    // beginIndex_apostrophe

    char* T_i = (char*)malloc(sizeof(char) * partLength);
    long* localSA = (long*)malloc(sizeof(long) * partLength);

//    printf("index\tT_i[]\tlocalSA[]\n");
    for(i = bi_i; i < bi_apostrophe; i++) {
        long locali = i - bi_i;
        T_i[locali] = T[i];
        localSA[locali] = locali;
//        printf("%ld\t%c\t%ld\n", locali, T_i[locali], localSA[locali]);
    }

    suffixArrayQuickSort(localSA, T_i, 0, partLength - 1);

    // copy sorted suffixes' indexes and their inverse to SA and SA_inverse
    printf("sorted suffixes(%ld): \n", partIndex);
    printf("begin_index: %ld\n", bi_i);
    printf("i\tch\tSA[]\tch_SA\n");
    for(i = bi_i; i < bi_apostrophe; i++) {
        long locali = i - bi_i;
        SA[i] = localSA[locali] + bi_i;

        printf("%ld\t", i - bi_i);
        printf("%c\t", T[i]);
        printf("%ld\t", SA[i] - bi_i);
        printf("%c\t", T[SA[i]]);
        printf("\n");
    }
}

/**
 * Calculate all order(suf_k, T') in an order k = l, l-1, ..., 1.
 *
 * \note suf[] is already sorted by lex-order.
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
//        printf("(%c) -> lc: %ld, rc: %ld\t", c, lc, rc);

        if(lc > rc) {
            /**
             * \caution Still doubting about this, I think it should be "orderValue = lc - 0 -bi_apostrophe;"
             */
            orderValue = lc - 0 - bi_apostrophe;   // this is modified ... on my own will
        } else {
            // find the max b that satisfies condition that for any order(cX, T'), Psi_T'[b] <= order(X, T')
            long max_b = 0;
            CSABinarySearchOrderValue(SA, Psi, lc, rc, prevOrderValue, &max_b);
            orderValue = max_b - bi_apostrophe;
        }

//        printf("orderValue(%ld): %ld\t", i - bi_i, orderValue);

        order[i - bi_i] = orderValue;
        prevOrderValue = orderValue;

//        printf("\n");
    }

    printf("T\': ");
    for(i = bi_apostrophe; i < arrayLength; i++) {
        printf("%c", T[SA[i]]);
    }
    printf("\n");

    printf("For sorted suf_k[] ...\n");
    for(i = 0; i < partLength; i++) {
        printf("suf[%ld](%ld) ~ %c...\torder(suf[%ld], T\'): %ld\n", i, SA[i + bi_i], T[SA[i + bi_i]], i,
               order[SA[i + bi_i] - bi_i]);
    }



    printf("\n");

}


/**
 * Calculate the Psi function for T_iT'.
 *
 * \note you should build func f and g of the part first, and then calculate Psi based on them.
 *
 * @param T DNA sequence (plus a '$')
 * @param SA use SA to store the startIndexes of suffixes sorted by lex-order
 * @param Psi Psi array of T -- the results is stored in this array
 * @param arrayLength length of T
 * @param partLength length of a part (n/log2(n))
 * @param partNum number of parts
 * @param partIndex index of this part
 * @param order order[] that stores result of func order(suf_k, T'), where suf_k is the k-th
 *      longest suffix of T_i and T' is combination of T_(i+1)...T_([n/l|+]).
 */
void mergeStepC(char* T, long* SA, long* Psi, long arrayLength, long partLength,
                long partNum, long partIndex, long* order) {
    printf("Merge Step (c)\n");
    long i = 0;
    long j = 0;
    long bi_i = (partIndex - 1) * partLength;   // beginIndex_i
    long bi_apostrophe = partIndex * partLength;    // beginIndex_apostrophe

    long* fFunc = (long*)malloc(sizeof(long) * (arrayLength - bi_apostrophe));
    long* gFunc = (long*)malloc(sizeof(long) * partLength);
    long* psiFunc = (long*)malloc(sizeof(long) * (arrayLength - bi_i));

    long num = 0;
    long maxIndex = 0;

    // construction of func f
    printf("Calculating func f ...\n");
    num = 0;
    maxIndex = 0;
    for(i = 0; i < arrayLength - bi_apostrophe; i++) {
        printf("//****\n");
        for(j = maxIndex; j < partLength; j++) {
            long orderValue = order[SA[j + bi_i] - bi_i];

            printf("suf[%ld](%c)\t", j, T[SA[j + bi_i]]);
            printf("order(suf[%ld], T\'): %ld\n", j, orderValue);

            if(orderValue > i) {
                break;
            } else {
                num++;
                maxIndex++;
            }
        }

        printf("f[%ld](%c): %ld\n", SA[i + bi_apostrophe] - bi_apostrophe, T[SA[i + bi_apostrophe]],
               i + num);
        fFunc[SA[i + bi_apostrophe] - bi_apostrophe] = i + num;
    }

    printf("///////\n");
    printf("i\tch\tfunc f[]\n");
    for(i = 0; i < arrayLength - bi_apostrophe; i++) {
        printf("%ld\t%c\t%ld\n", i, T[i + bi_apostrophe], fFunc[i]);
    }

    // construction of func g
    printf("Calculating func g ...\n");
    num = 0;
    maxIndex = 0;
    for(i = 0; i < partLength; i++) {
        long orderValue_i = order[SA[i + bi_i] - bi_i];
        for(j = maxIndex; j < partLength; j++) {
            long orderValue_temp = order[SA[j + bi_i] - bi_i];
            /**
             * \caution still doubting about this. The original is "<=".
             */
            if(orderValue_temp < orderValue_i) {
                maxIndex++;
                num++;
            }
        }
//        printf("g[%ld](%c): %ld\n", SA[i + bi_i] - bi_i, T[SA[i + bi_i]], orderValue_i + num);
        gFunc[SA[i + bi_i] - bi_i] = orderValue_i + num;
    }

    printf("///////\n");
    printf("i\tch\tfunc g[]\n");
    for(i = 0; i < partLength; i++) {
        printf("%ld\t%c\t%ld\n", i, T[i + bi_i], gFunc[i]);
    }


    // construction of func Psi
    printf("Calculating func psi ...\n");
    long t = 0;
    for(t = 0, i = 0, j = 0; t < arrayLength - bi_i; t++) {
        if(t == gFunc[partLength - 1]) {
            psiFunc[t] = fFunc[Psi[0 + bi_apostrophe]];
        } else if(t == fFunc[i]) {
            psiFunc[t] = fFunc[Psi[i + bi_apostrophe]];
            i++;
        } else {
            psiFunc[t] = gFunc[j + 1];
            j++;
        }
    }

    printf("///////\n");
    printf("i\tch\tfunc psi[]\n");
    for(i = 0; i < arrayLength - bi_i; i++) {
        printf("%ld\t%c\t%ld\n", i, T[i + bi_i], psiFunc[i]);
    }


}


void _fgpsiFuncTest(){
    long i = 0;
    long j = 0;

    char* T_apostrophe = "cagac$";
    char* T_i = "gca";

    long fFuncLength = strlen(T_apostrophe);
    long gFuncLength = strlen(T_i);
    long* fFunc = (long*)malloc(sizeof(long) * fFuncLength);
    long* gFunc = (long*)malloc(sizeof(long) * gFuncLength);

    long num = 0;
    long maxIndex = 0;

    long order[3] = {5, 3, 1};

    // construction of func f
    printf("Calculating func f ...\n");
    num = 0;
    maxIndex = 0;
    for(i = 0; i < fFuncLength; i++) {
        printf("//****\n");
        for(j = maxIndex; j < gFuncLength; j++) {
            long orderValue = order[j];

            printf("suf[%ld](%c)\t", j, T_i[j]);
            printf("order(suf[%ld], T\'): %ld\n", j, orderValue);

            if(orderValue > i) {
                break;
            } else {
                num++;
                maxIndex++;
            }
        }

        printf("f[%ld](%c): %ld\n", i, T_apostrophe[i], i + num);
        fFunc[i] = i + num;
    }

    printf("///////\n");
    printf("i\tch\tfunc f[]\n");
    for(i = 0; i < fFuncLength; i++) {
        printf("%ld\t%c\t%ld\n", i, T_apostrophe[i], fFunc[i]);
    }

    // construction of func g
    printf("Calculating func g ...\n");
    num = 0;
    maxIndex = 0;
    for(i = 0; i < gFuncLength; i++) {
        long orderValue_i = order[i];
        for(j = maxIndex; j < gFuncLength; j++) {
            long orderValue_temp = order[j];
            /**
             * \caution still doubting about this. The original is "<=".
             */
            if(orderValue_temp <= orderValue_i) {
                maxIndex++;
                num++;
            }
        }
//        printf("g[%ld](%c): %ld\n", SA[i + bi_i] - bi_i, T[SA[i + bi_i]], orderValue_i + num);
        gFunc[i] = orderValue_i + num;
    }

    printf("///////\n");
    printf("i\tch\tfunc g[]\n");
    for(i = 0; i < gFuncLength; i++) {
        printf("%ld\t%c\t%ld\n", i, T_i[i], gFunc[i]);
    }


}


#endif // MERGESTEP_H_INCLUDED
