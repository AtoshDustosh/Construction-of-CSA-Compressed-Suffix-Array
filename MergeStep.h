#ifndef MERGESTEP_H_INCLUDED
#define MERGESTEP_H_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "SABuildFunc.h"
#include "HelperFunction.h"

void mergeStepA(char* T, int* SA, int* SA_inverse, int arrayLength, int partLength, int partIndex);
void mergeStepB(char* T, int* SA, int* Psi, int arrayLength, int partLength,
                int partIndex, int* order);
void mergeStepC(char* T, int* SA, int* SA_inverse, int* Psi, int arrayLength, int partLength,
                int partIndex, int* order);
void processFuncF(char* T, int* SA, int* order, int arrayLength, int partLength, int partIndex,
                  int* fFunc);
void processFuncG(char* T, int* SA, int* SA_inverse, int* order, int arrayLength, int partLength,
                  int partIndex, int* gFunc);
void processFuncPsi(char* T, int* SA, int* Psi, int arrayLength, int partLength, int partIndex,
                    int* psiFunc, int* fFunc, int* gFunc);
void _fgpsiFuncTest();


/**
 * Sort suffixes (suf1, suf2, ..., sufl) to find the lex-order among themselves.
 *
 * @param T DNA sequence (plus a '$')
 * @param SA use SA to store the startIndexes of suffixes sorted by lex-order
 * @param SA_inverse inverse of SA
 * @param arrayLength length of T
 * @param partLength length of a part (n/log2(n))
 * @param partIndex index of this part
 */
void mergeStepA(char* T, int* SA, int* SA_inverse, int arrayLength, int partLength, int partIndex) {
    printf("Merge Step (a)\n");
    int i = 0;
    int bi_i = (partIndex - 1) * partLength;   // beginIndex_i
    int bi_apostrophe = partIndex * partLength;    // beginIndex_apostrophe

    /**
     *  \note this operation of localization resulted in the strange bug
     *     of order(cX, T') calculating -- ignored all characters after
     *     T[bi_apostrophe - 1] when comparing suffixes
     */
    char* T_i = (char*)malloc(sizeof(char) * (arrayLength - bi_i));
    int* localSA = (int*)malloc(sizeof(int) * (arrayLength - bi_i));

    if(T_i == NULL || localSA == NULL) {
        printf("System memory not enough. \n");
        exit(-1);
    }

    // initialize array T_i and localSA
    for(i = bi_i; i < arrayLength; i++) {
        int locali = i - bi_i;
        T_i[locali] = T[i];
        localSA[locali] = locali;
    }

    suffixArrayQuickSort(localSA, T_i, 0, partLength - 1);

    // copy sorted suffixes' indexes and their inverse to SA and SA_inverse
//    printf("sorted suffixes(%d): \n", partIndex);
//    printf("begin_index: %d\n", bi_i);
//    printf("i\tch\tSA[]\tch_SA\n");
    for(i = bi_i; i < bi_apostrophe; i++) {
        int locali = i - bi_i;
        SA[i] = localSA[locali] + bi_i;
        SA_inverse[SA[i]] = i;

//        printf("%d\t", i - bi_i);
//        printf("%c\t", T[i]);
//        printf("%d\t", SA[i] - bi_i);
//        printf("%c\t", T[SA[i]]);
//        printf("\n");
    }

    free(T_i);
    free(localSA);
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
 * @param partIndex index of this part
 * @param order order[] that stores result of func order(suf_k, T'), where suf_k is the
 *      k-th intest suffix of T_i and T' is combination of T_(i+1)...T_([n/l|+]).
 */
void mergeStepB(char* T, int* SA, int* Psi, int arrayLength, int partLength,
                int partIndex, int* order) {
    printf("Merge Step (b)\n");
    int i = 0;
//    long startTime = 0;
//    long endTime = 0;
    int bi_i = (partIndex - 1) * partLength;   // beginIndex_i
    int bi_apostrophe = partIndex * partLength;    // beginIndex_apostrophe

    // T_apostrophe - start from bi_apostrophe and ends at arrayLength
    // note that lc and rc are not actual bounds of the character c
    int prevOrderValue = Psi[bi_apostrophe];    // "global" index value
    for(i = bi_apostrophe - 1; i >= bi_i; i--) {
        char c = T[i];  // the character in the formula
        int lc = bi_apostrophe;
        int rc = arrayLength - 1;
        int orderValue = 0;    // "local" index value
//        startTime = clock();
        CSABinaryBoundSearch(T, SA, c, &lc, &rc);
//        endTime = clock();
        // T[SA[lc]] ~ T[SA[rc]] represents the field of c
        // implement of condition ×[b] -> ×[SA[b]], lc <= b <= rc
//        printf("(%c) -> lc: %d, rc: %d\t", c, lc - bi_apostrophe, rc - bi_apostrophe);
//        printf("%d -> locate lc and rc time: %ld\t", (bi_apostrophe - 1 - i), endTime - startTime);

//        startTime = clock();
        if(lc > rc) {
            orderValue = lc - 1 - bi_apostrophe;
        } else {
            // find the max b that satisfies condition that for any order(cX, T'),
            //      Psi_T'[b] <= order(X, T')
            int max_b = 0;
            CSABinarySearchOrderValue(SA, Psi, lc, rc, prevOrderValue, &max_b);
            if(max_b == -1) {
                orderValue = lc - 1 - bi_apostrophe;
            } else {
                orderValue = max_b - bi_apostrophe;
            }
        }
//        endTime = clock();
//        printf("locate order value time: %ld\t", endTime - startTime);

//        printf("prevOrderValue: %d\torderValue(%d): %d\t", prevOrderValue - bi_apostrophe, i - bi_i,
//               orderValue);

        order[i - bi_i] = orderValue;
        prevOrderValue = orderValue + bi_apostrophe;

//        printf("\n");
    }

//    printf("sorted T\': ");
//    for(i = bi_apostrophe; i < arrayLength; i++) {
//        printf("%c", T[SA[i]]);
//    }
//    printf("\n");
//
//    printf("unsorted T\': ");
//    for(i = bi_apostrophe; i < arrayLength; i++) {
//        printf("%c", T[i]);
//    }
//    printf("\n");
//
//    printf("T_i: ");
//    for(i = bi_i; i < bi_apostrophe; i++) {
//        printf("%c", T[i]);
//    }
//    printf("\n");
//
//    printf("For sorted suf_k[] ...\n");
//    for(i = 0; i < partLength; i++) {
//        printf("suf[%d](%d) ~ %c...\torder(suf[%d], T\'): %d\n", i, SA[i + bi_i] - bi_i,
//               T[SA[i + bi_i]], i, order[SA[i + bi_i] - bi_i]);
//    }

//    printf("\n");

}


/**
 * Calculate the Psi function for T_iT'.
 *
 * \note func f, g, psi corresponds to their local indexes of T', T_i and T_iT'
 *
 * @param T DNA sequence (plus a '$')
 * @param SA use SA to store the startIndexes of suffixes sorted by lex-order
 * @param SA_inverse inverse of SA
 * @param Psi Psi array of T -- the results is stored in this array
 * @param arrayLength length of T
 * @param partLength length of a part (n/log2(n))
 * @param partIndex index of this part
 * @param order order[] that stores result of func order(suf_k, T'), where suf_k is the
 *      k-th longest suffix of T_i and T' is combination of T_(i+1)...T_([n/l|+]).
 */
void mergeStepC(char* T, int* SA, int* SA_inverse, int* Psi, int arrayLength, int partLength,
                int partIndex, int* order) {
    long startTime = 0;
    long endTime = 0;
    printf("Merge Step (c)\n");
    int i = 0;
    int bi_i = (partIndex - 1) * partLength;   // beginIndex_i
    int bi_apostrophe = partIndex * partLength;    // beginIndex_apostrophe

    int* fFunc = (int*)malloc(sizeof(int) * (arrayLength - bi_apostrophe));
    int* gFunc = (int*)malloc(sizeof(int) * partLength);
    int* psiFunc = (int*)malloc(sizeof(int) * (arrayLength - bi_i));

    if(fFunc == NULL || gFunc == NULL || psiFunc == NULL) {
        printf("System memory not enough. \n");
        exit(-1);
    }

    startTime = clock();
    // construction of func f
    processFuncF(T, SA, order, arrayLength, partLength, partIndex, fFunc);
    endTime = clock();
    printf("calculating func f takes time: %ld\n", endTime - startTime);

    startTime = clock();
    // construction of func g
    processFuncG(T, SA, SA_inverse, order, arrayLength, partLength, partIndex, gFunc);
    endTime = clock();
    printf("calculating func g takes time: %ld\n", endTime - startTime);

    startTime = clock();
    // construction of func Psi
    processFuncPsi(T, SA, Psi, arrayLength, partLength, partIndex, psiFunc, fFunc, gFunc);
    endTime = clock();
    printf("calculating func psi takes time: %ld\n", endTime - startTime);


    // update SA[] of T[] to get ready for next iteration
    // reconstruct the new SA[] from new Psi[] - according to mathematical theories
    printf("Update SA[] of T[] to get ready for next iteration. \n");
    /**
     * \note the original reconstruction method is to quick sort a newly-initialized SA[].
     *  But that appears very time-consuming. Thus I choose to make use of Psi[] for
     *  the reconstruction of SA[].
     */
    startTime = clock();
    int x = bi_i;
    SA[Psi[x]] = bi_i;
    x = Psi[bi_i];
    for(i = bi_i; i < arrayLength; i++) {
//        printf("%d -> x: %d, Psi[x]: %d, SA[Psi[x]]: %d\n",
//               i - bi_i, x - bi_i, Psi[x] - bi_i, SA[Psi[x]] - bi_i);
        if(x == bi_i) {
            continue;   // this is necessary - cannot be deleted
        }
        SA[Psi[x]] = SA[x] + 1;
        x = Psi[x];
    }
    endTime = clock();
    printf("update SA[] takes time: %ld\n", endTime - startTime);

    free(fFunc);
    free(gFunc);
    free(psiFunc);
}

/**
 * A function used for processing func f.
 */
void processFuncF(char* T, int* SA, int* order, int arrayLength, int partLength, int partIndex,
                  int* fFunc) {
    int i = 0;
    int j = 0;
    int bi_i = (partIndex - 1) * partLength;   // beginIndex_i
    int bi_apostrophe = partIndex * partLength;    // beginIndex_apostrophe

    // construction of func f
    printf("Calculating func f (length: %d) ...\n", arrayLength - bi_apostrophe);
    int num = 0;
    int maxIndex = 0;
    // calculate by lex-order for the convenience of calculating #(order(suf_k, T') <= j)
    for(i = 0; i < arrayLength - bi_apostrophe; i++) {
//        printf("//****\n");
        // make use of the increasing values of order with increasing lex-order
        for(j = maxIndex; j < partLength; j++) {
            int orderValue = order[SA[j + bi_i] - bi_i];

//            printf("suf[%d](%c)\t", j, T[SA[j + bi_i]]);
//            printf("order(suf[%d], T\'): %d\n", j, orderValue);
            /**
             * \caution the original is ">" but as index of this implementation starts from 0
             *      rather than 1, the program should be adjusted.
             */
            if(orderValue >= i) {
                break;
            } else {
                num++;
                maxIndex++;
            }
        }

//        printf("locate maxIndex: %d\tnum: %d\n", maxIndex, num);
//        printf("f[%d](%c): %d\n", SA[i + bi_apostrophe] - bi_apostrophe,
//               T[SA[i + bi_apostrophe]], i + num);
        fFunc[SA[i + bi_apostrophe] - bi_apostrophe] = i + num;
    }

    printf("///////\n");
    printf("i\tch\tfunc f[]\n");
    for(i = 0; i < arrayLength - bi_apostrophe; i++) {
        printf("%d\t%c\t%d\n", i, T[SA[i + bi_apostrophe]],
               fFunc[SA[i + bi_apostrophe] - bi_apostrophe]);
        i = i + (arrayLength - bi_apostrophe) / 10;
    }
}

/**
 * A function used for processing func g.
 */
void processFuncG(char* T, int* SA, int* SA_inverse, int* order, int arrayLength, int partLength,
                  int partIndex, int* gFunc) {
    int i = 0;
    int bi_i = (partIndex - 1) * partLength;   // beginIndex_i

    printf("Calculating func g (length: %d) ...\n", partLength);
    // calculate by lex-order for the convenience of calculating #(suf_k <= suf_i)
    for(i = 0; i < partLength; i++) {
        int orderValue_i = order[SA[i + bi_i] - bi_i];
        int num = SA_inverse[SA[i + bi_i]] - bi_i + 1;
//        printf("locate maxIndex: %d\tnum: %d\n", maxIndex, num);
//        printf("g[%d](%c): %d\n", SA[i + bi_i] - bi_i, T[SA[i + bi_i]],
//               orderValue_i + num);
        gFunc[SA[i + bi_i] - bi_i] = orderValue_i + num;
    }


    printf("///////\n");
    printf("i\tch\tSA\tch_SA\tg_SA[]\n");
    for(i = 0; i < partLength; i++) {
        printf("%d\t%c\t", i, T[i + bi_i]);
        printf("%d\t%c\t", SA[i + bi_i] - bi_i, T[SA[i + bi_i]]);
        printf("%d\t", gFunc[SA[i + bi_i] - bi_i]);
        printf("\n");
        i = i + partLength / 10;
    }


}


/**
 * A function used for processing func psi.
 */
void processFuncPsi(char* T, int* SA, int* Psi, int arrayLength, int partLength, int partIndex,
                    int* psiFunc, int* fFunc, int* gFunc) {
    int i = 0;
    int j = 0;
    int bi_i = (partIndex - 1) * partLength;   // beginIndex_i
    int bi_apostrophe = partIndex * partLength;    // beginIndex_apostrophe

    // construction of func Psi
    printf("Calculating func psi ...\n");
    int t = 0;
    for(t = 0; t < arrayLength - bi_i; t++) {
//        printf("t: %d, i: %d, j: %d\n", t, i, j);
        if(t == gFunc[partLength - 1]) {
//            int index1 = SA[Psi[0 + bi_apostrophe]] - bi_apostrophe;
//            if(index1 >= arrayLength - bi_apostrophe) {
//                printf("f func index error. \n");
//                exit(-2);
//            }

            psiFunc[t] = fFunc[SA[Psi[0 + bi_apostrophe]] - bi_apostrophe];
//            printf("... processed\n");
        } else if(i < arrayLength - bi_apostrophe
                  && t == fFunc[SA[i + bi_apostrophe] - bi_apostrophe]) {
//            int index1 = SA[Psi[0 + bi_apostrophe]] - bi_apostrophe;
//            int index2 = SA[Psi[i + bi_apostrophe]] - bi_apostrophe;
//            if(index1 >= arrayLength - bi_apostrophe || index2 >= arrayLength - bi_apostrophe) {
//                printf("f func index error. \n");
//                exit(-2);
//            }

            psiFunc[t] = fFunc[SA[Psi[i + bi_apostrophe]] - bi_apostrophe];
            i++;
//            printf("... processed\n");
        } else {
//            int index3 = gFunc[j] - bi_i;
//            if(index3 >= arrayLength - bi_i || (j + 1) > partLength) {
//                printf("g func index error. \n");
//                exit(-2);
//            }

            psiFunc[gFunc[j]] = gFunc[j + 1];
            j++;
//            printf("... processed\n");
        }
    }
    psiFunc[0] = gFunc[0];  // this is necessary for fixing the bug

    printf("///////\n");
    printf("i\tch\t");
    printf("SA[]\tch_SA\t");
    printf("Psi[]\t");
    printf("\n");
    j = bi_i;
    for(i = bi_i; i < arrayLength; i++) {
//         update the original Psi func
        Psi[i] = psiFunc[i - bi_i] + bi_i;
        if(i == j) {
            printf("%d\t%c\t", i - bi_i, T[i]);
            printf("%d\t%c\t", SA[i] - bi_i, T[SA[i]]);
            printf("%d\t", Psi[i] - bi_i);
            printf("\n");
            j = j + (arrayLength - bi_i) / 10;
        }
    }
}

void _fgpsiFuncTest() {
    printf("************** _fgpsiFuncTest **************\n");
    int i = 0;
    int j = 0;

    char* T_apostrophe = "gac$";
    char* T_i = "aca";

    int fFuncLength = strlen(T_apostrophe);
    int gFuncLength = strlen(T_i);
    int* fFunc = (int*)malloc(sizeof(int) * fFuncLength);
    int* gFunc = (int*)malloc(sizeof(int) * gFuncLength);
    int* psiFunc = (int*)malloc(sizeof(int) * (fFuncLength + gFuncLength));

    int num = 0;
    int maxIndex = 0;

    int order[3] = {1, 2, 1}; // order(suf_k, T') by lex-order of suffix[]
    int SA_i[3] = {2, 0, 1}; // the corresponding lex-order of suffix[] from T_i
    int SA[6] = {3, 1, 2, 0};
    int Psi[6] = {3, 2, 0, 1};

    // construction of func f
    printf("Calculating func f ...\n");
    num = 0;
    maxIndex = 0;
    // calculate by lex-order for the convenience of calculating #(order(suf_k, T') <= j)
    for(i = 0; i < fFuncLength; i++) {
        printf("//****\n");
        // make use of the increasing values of order with increasing lex-order
        for(j = maxIndex; j < gFuncLength; j++) {
            int orderValue = order[SA_i[j]];

            printf("suf[%d](%c)\t", j, T_i[SA_i[j]]);
            printf("order(suf[%d], T\'): %d\n", j, orderValue);
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

        printf("f[%d](%c): %d\n", SA[i], T_apostrophe[SA[i]], i + num);
        fFunc[SA[i]] = i + num;
    }

    printf("///////\n");
    printf("i\tch\tfunc f[]\n");
    for(i = 0; i < fFuncLength; i++) {
        printf("%d\t%c\t%d\n", i, T_apostrophe[SA[i]], fFunc[SA[i]]);
    }

    // construction of func g
    printf("Calculating func g ...\n");
    num = 0;
    maxIndex = 0;
    // calculate by lex-order for the convenience of calculating #(suf_k <= suf_j)
    for(i = 0; i < gFuncLength; i++) {
        int orderValue_i = order[SA_i[i]];
        int suffix_i = SA_i[i];
        for(j = maxIndex; j < gFuncLength; j++) {
            int suffix_k = SA_i[j];
            if(compareSuffix(suffix_k, suffix_i, T_i) <= 0) {
                maxIndex++;
                num++;
            }
        }
        printf("g[%d](%c): %d\n", SA_i[i], T_i[SA_i[i]], orderValue_i + num);
        gFunc[SA_i[i]] = orderValue_i + num;
    }

    printf("///////\n");
    printf("i\tch\tfunc g[]\n");
    for(i = 0; i < gFuncLength; i++) {
        printf("%d\t%c\t%d\n", i, T_i[i], gFunc[i]);
    }

    // construction of func Psi
    printf("Calculating func psi ...\n");
    int t = 0;
    for(t = 0, i = 0, j = 0; t < fFuncLength + gFuncLength; t++) {
        if(t == gFunc[gFuncLength - 1]) {
            psiFunc[t] = fFunc[SA[Psi[0]]];
        } else if(SA[i] < fFuncLength && t == fFunc[SA[i]]) {
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
            printf("%d\t%c\t", i, T_i[i]);
        } else {
            printf("%d\t%c\t", i, T_apostrophe[i - gFuncLength]);
        }
        printf("%d\n", psiFunc[i]);
    }

    // update the SA[] of T[] to get ready for next iteration
//    printf("Update SA[] of T[] to get ready for next iteration. \n");
//    int SA_new[9] = {0};
//    char* T = "gcacagac$";
//    for(i = 0; i < 9; i++) {
//        SA_new[i] = i;
//    }
//    suffixArrayQuickSort(SA_new, T, 0, 9 - 1);
//    for(i = 0; i < 9; i++) {
//        printf("%d\t%c\t%d\t%c\n", i, T[i], SA_new[i], T[SA_new[i]]);
//    }

}


#endif // MERGESTEP_H_INCLUDED
