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

    /**
     *  \note this operation of localization results in the strange bug
     *     of order(cX, T') calculating.
     */
    char* T_i = (char*)malloc(sizeof(char) * (arrayLength - bi_i));
    long* localSA = (long*)malloc(sizeof(long) * (arrayLength - bi_i));

    // initialize array T_i and localSA
    for(i = bi_i; i < arrayLength; i++) {
        long locali = i - bi_i;
        T_i[locali] = T[i];
        localSA[locali] = locali;
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
 * \note suf[] is already sorted by lex-order
 * \note values of order[] correspond to the local index of T', not the global index of T
 * \note indexes of order[] correspond to the sorted suffixes' global indexes of T
 *
 * @param T DNA sequence (plus a '$')
 * @param SA use SA to store the startIndexes of suffixes sorted by lex-order
 * @param Psi Psi array of T
 * @param arrayLength length of T
 * @param partLength length of a part (n/log2(n))
 * @param partNum number of parts
 * @param partIndex index of this part
 * @param order order[] that stores result of func order(suf_k, T'), where suf_k is the
 *      k-th longest suffix of T_i and T' is combination of T_(i+1)...T_([n/l|+]).
 */
void mergeStepB(char* T, long* SA, long* Psi, long arrayLength, long partLength,
                long partNum, long partIndex, long* order) {
    /**
     * \TODO perhaps there exists bug in this block
     */
    printf("Merge Step (b)\n");
    long i = 0;
    long bi_i = (partIndex - 1) * partLength;   // beginIndex_i
    long bi_apostrophe = partIndex * partLength;    // beginIndex_apostrophe

    // T_apostrophe - start from bi_apostrophe and ends at arrayLength
    // note that lc and rc are not actual bounds of the character c
    long prevOrderValue = Psi[bi_apostrophe];    // "global" index value
    for(i = bi_apostrophe - 1; i >= bi_i; i--) {
        char c = T[i];  // the character in the formula
        long lc = bi_apostrophe;
        long rc = arrayLength - 1;
        long orderValue = 0;    // "local" index value
        CSABinaryBoundSearch(T, SA, c, &lc, &rc);
        // T[SA[lc]] ~ T[SA[rc]] represents the field of c
        // implement of condition ×[b] -> ×[SA[b]], lc <= b <= rc
        printf("(%c) -> lc: %ld, rc: %ld\t", c, lc - bi_apostrophe, rc - bi_apostrophe);

        if(lc > rc) {
            orderValue = lc - 1 - bi_apostrophe;
        } else {
            // find the max b that satisfies condition that for any order(cX, T'),
            //      Psi_T'[b] <= order(X, T')
            long max_b = 0;
            CSABinarySearchOrderValue(SA, Psi, lc, rc, prevOrderValue, &max_b);
            if(max_b == -1) {
                orderValue = lc - 1 - bi_apostrophe;
            } else {
                orderValue = max_b - bi_apostrophe;
            }
        }

        printf("prevOrderValue: %ld\torderValue(%ld): %ld\t", prevOrderValue - bi_apostrophe, i - bi_i,
               orderValue);

        order[i - bi_i] = orderValue;
        prevOrderValue = orderValue + bi_apostrophe;

        printf("\n");
    }

    printf("sorted T\': ");
    for(i = bi_apostrophe; i < arrayLength; i++) {
        printf("%c", T[SA[i]]);
    }
    printf("\n");

    printf("unsorted T\': ");
    for(i = bi_apostrophe; i < arrayLength; i++) {
        printf("%c", T[i]);
    }
    printf("\n");

    printf("T_i: ");
    for(i = bi_i; i < bi_apostrophe; i++) {
        printf("%c", T[i]);
    }
    printf("\n");

    printf("For sorted suf_k[] ...\n");
    for(i = 0; i < partLength; i++) {
        printf("suf[%ld](%ld) ~ %c...\torder(suf[%ld], T\'): %ld\n", i, SA[i + bi_i] - bi_i,
               T[SA[i + bi_i]], i, order[SA[i + bi_i] - bi_i]);
    }



    printf("\n");

}


/**
 * Calculate the Psi function for T_iT'.
 *
 * \note func f, g, psi corresponds to their local indexes of T', T_i and T_iT'
 *
 * @param T DNA sequence (plus a '$')
 * @param SA use SA to store the startIndexes of suffixes sorted by lex-order
 * @param Psi Psi array of T -- the results is stored in this array
 * @param arrayLength length of T
 * @param partLength length of a part (n/log2(n))
 * @param partNum number of parts
 * @param partIndex index of this part
 * @param order order[] that stores result of func order(suf_k, T'), where suf_k is the
 *      k-th longest suffix of T_i and T' is combination of T_(i+1)...T_([n/l|+]).
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
    // calculate by lex-order for the convenience of calculating #(order(suf_k, T') <= j)
    for(i = 0; i < arrayLength - bi_apostrophe; i++) {
        printf("//****\n");
        // make use of the increasing values of order with increasing lex-order
        for(j = maxIndex; j < partLength; j++) {
            long orderValue = order[SA[j + bi_i] - bi_i];

            printf("suf[%ld](%c)\t", j, T[SA[j + bi_i]]);
            printf("order(suf[%ld], T\'): %ld\n", j, orderValue);
            /**
             * \caution the original is ">"
             */
            if(orderValue >= i) {
                break;
            } else {
                num++;
                maxIndex++;
            }
        }

        printf("f[%ld](%c): %ld\n", SA[i + bi_apostrophe] - bi_apostrophe,
               T[SA[i + bi_apostrophe]], i + num);
        fFunc[SA[i + bi_apostrophe] - bi_apostrophe] = i + num;
    }

    printf("///////\n");
    printf("i\tch\tfunc f[]\n");
    for(i = 0; i < arrayLength - bi_apostrophe; i++) {
        printf("%ld\t%c\t%ld\n", i, T[SA[i + bi_apostrophe]],
               fFunc[SA[i + bi_apostrophe] - bi_apostrophe]);
    }

    // construction of func g
    /**
     * \TODO fix the bug - in the 2nd iteration, there exists duplicates in g func.
     * Maybe it has something to do with the bug in main func concerning "void testSet()"
     */
    printf("Calculating func g ...\n");
    num = 0;
    maxIndex = 0;
    // calculate by lex-order for the convenience of calculating #(suf_k <= suf_i)
    for(i = 0; i < partLength; i++) {
        long orderValue_i = order[SA[i + bi_i] - bi_i];
        long suffix_i = SA[i + bi_i];
        for(j = maxIndex; j < partLength; j++) {
            long suffix_k = SA[j + bi_i];
            if(compareSuffix(suffix_k, suffix_i, T) <= 0) {
                maxIndex++;
                num++;
            }
        }
        printf("g[%ld](%c): %ld\n", SA[i + bi_i] - bi_i, T[SA[i + bi_i]],
               orderValue_i + num);
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
            psiFunc[t] = fFunc[SA[Psi[0 + bi_apostrophe]] - bi_apostrophe];
        } else if(t == fFunc[SA[i + bi_apostrophe] - bi_apostrophe]) {
            psiFunc[t] = fFunc[SA[Psi[i + bi_apostrophe]] - bi_apostrophe];
            i++;
        } else {
            psiFunc[gFunc[j]] = gFunc[j + 1];
            j++;
        }
    }
    psiFunc[0] = gFunc[0];  // this is necessary for fixing the bug

    // update SA[] of T[] to get ready for next iteration
    printf("Update SA[] of T[] to get ready for next iteration. \n");
    for(i = bi_i; i < arrayLength; i++) {
        SA[i] = i;
    }
    suffixArrayQuickSort(SA, T, bi_i, arrayLength - 1);

    printf("///////\n");
    printf("i\tch\t");
    printf("SA[]\tch_SA\t");
    printf("Psi[]\t");
    printf("\n");
    for(i = bi_i; i < arrayLength; i++) {
//        printf("%ld\t%c\t%ld\n", i, T[i], psiFunc[i - bi_i]);
        // update the original Psi func
        Psi[i] = psiFunc[i - bi_i] + bi_i;
        printf("%ld\t%c\t", i - bi_i, T[i]);
        printf("%ld\t%c\t", SA[i] - bi_i, T[SA[i]]);
        printf("%ld\t", Psi[i] - bi_i);
        printf("\n");
    }


}


void _fgpsiFuncTest() {
    printf("************** _fgpsiFuncTest **************\n");
    long i = 0;
    long j = 0;

    char* T_apostrophe = "cagac$";
    char* T_i = "gca";

    long fFuncLength = strlen(T_apostrophe);
    long gFuncLength = strlen(T_i);
    long* fFunc = (long*)malloc(sizeof(long) * fFuncLength);
    long* gFunc = (long*)malloc(sizeof(long) * gFuncLength);
    long* psiFunc = (long*)malloc(sizeof(long) * (fFuncLength + gFuncLength));

    long num = 0;
    long maxIndex = 0;

    long order[3] = {5, 3, 1}; // order(suf_k, T') by lex-order of suffix[]
    long SA_i[3] = {2, 1, 0}; // the corresponding lex-order of suffix[] from T_i
    long SA[6] = {5, 3, 1, 4, 0, 2};
    long Psi[6] = {4, 3, 5, 0, 2, 1};

    // construction of func f
    printf("Calculating func f ...\n");
    num = 0;
    maxIndex = 0;
    // calculate by lex-order for the convenience of calculating #(order(suf_k, T') <= j)
    for(i = 0; i < fFuncLength; i++) {
        printf("//****\n");
        // make use of the increasing values of order with increasing lex-order
        for(j = maxIndex; j < gFuncLength; j++) {
            long orderValue = order[SA_i[j]];

            printf("suf[%ld](%c)\t", j, T_i[SA_i[j]]);
            printf("order(suf[%ld], T\'): %ld\n", j, orderValue);
            /**
             * \caution the original is ">"
             */
            if(orderValue >= i) {
                break;
            } else {
                num++;
                maxIndex++;
            }
        }

        printf("f[%ld](%c): %ld\n", SA[i], T_apostrophe[SA[i]], i + num);
        fFunc[SA[i]] = i + num;
    }

    printf("///////\n");
    printf("i\tch\tfunc f[]\n");
    for(i = 0; i < fFuncLength; i++) {
        printf("%ld\t%c\t%ld\n", i, T_apostrophe[SA[i]], fFunc[SA[i]]);
    }

    // construction of func g
    printf("Calculating func g ...\n");
    num = 0;
    maxIndex = 0;
    // calculate by lex-order for the convenience of calculating #(suf_k <= suf_j)
    for(i = 0; i < gFuncLength; i++) {
        long orderValue_i = order[SA_i[i]];
        long suffix_i = SA_i[i];
        for(j = maxIndex; j < gFuncLength; j++) {
            long suffix_k = SA_i[j];
            if(compareSuffix(suffix_k, suffix_i, T_i) <= 0) {
                maxIndex++;
                num++;
            }
        }
        printf("g[%ld](%c): %ld\n", SA_i[i], T_i[SA_i[i]], orderValue_i + num);
        gFunc[SA_i[i]] = orderValue_i + num;
    }

    printf("///////\n");
    printf("i\tch\tfunc g[]\n");
    for(i = 0; i < gFuncLength; i++) {
        printf("%ld\t%c\t%ld\n", i, T_i[i], gFunc[i]);
    }

    // construction of func Psi
    printf("Calculating func psi ...\n");
    long t = 0;
    for(t = 0, i = 0, j = 0; t < fFuncLength + gFuncLength; t++) {
        if(t == gFunc[gFuncLength - 1]) {
            psiFunc[t] = fFunc[SA[Psi[0]]];
        } else if(t == fFunc[SA[i]]) {
            psiFunc[t] = fFunc[SA[Psi[i]]];
            i++;
        } else {
            psiFunc[gFunc[j]] = gFunc[j + 1];
            j++;
        }
    }
    psiFunc[0] = gFunc[0];

    printf("///////\n");
    printf("i\tch\tfunc psi[]\n");
    for(i = 0; i < fFuncLength + gFuncLength; i++) {
        if(i < gFuncLength) {
            printf("%ld\t%c\t", i, T_i[i]);
        } else {
            printf("%ld\t%c\t", i, T_apostrophe[i - gFuncLength]);
        }
        printf("%ld\n", psiFunc[i]);
    }

    // update the SA[] of T[] to get ready for next iteration
    printf("Update SA[] of T[] to get ready for next iteration. \n");
    long SA_new[9] = {0};
    char* T = "gcacagac$";
    for(i = 0; i < 9; i++) {
        SA_new[i] = i;
    }
    suffixArrayQuickSort(SA_new, T, 0, 9 - 1);
    for(i = 0; i < 9; i++) {
        printf("%ld\t%c\t%ld\t%c\n", i, T[i], SA_new[i], T[SA_new[i]]);
    }

}


#endif // MERGESTEP_H_INCLUDED
